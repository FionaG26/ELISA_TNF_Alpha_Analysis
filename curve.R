#0. Load Required Packages
# If not installed, run: install.packages(c("tidyverse", "drc", "lme4", "lmerTest", "patchwork"))
library(tidyverse)
library(drc)       # For 4PL curve fitting
library(lme4)      # For linear mixed-effects models
library(lmerTest)  # To extract p-values from lme4 objects
library(patchwork) # To arrange multiple plots beautifully
library(ggpubr)
# Read the synthetic dataset
# Assumes 'synthetic_elisa_data.csv' is in your current working directory
if (!file.exists("C:/Users/fiona/Documents/elisa/MouseTNF_elisa_data.csv")) {
  stop("Data file 'C:/Users/fiona/Documents/elisa/MouseTNF_elisa_data.csv' not found. Please ensure it is in your working directory.")
}
df <- read.csv("C:/Users/fiona/Documents/elisa/MouseTNF_elisa_data.csv")


# ==============================================================================
# PART 1: Blank Correction & Data Preprocessing
# ==============================================================================
cat("\n--- [1/4] Executing Blank Corrections ---\n")

# Calculate background noise (mean blank OD) independently for each plate
blanks <- df %>%
  filter(Sample_Type == "Blank") %>%
  group_by(Plate) %>%
  summarize(Mean_Blank_OD = mean(Raw_OD450), .groups = 'drop')

# Subtract background baseline from all wells
df_corrected <- df %>%
  left_join(blanks, by = "Plate") %>%
  mutate(Corrected_OD = Raw_OD450 - Mean_Blank_OD)


# ==============================================================================
# PART 2: Technical Replicate Quality Control (QC)
# ==============================================================================
cat("\n--- [2/4] Performing Technical Replicate Quality Control ---\n")

# Calculate %CV (Coefficient of Variation) for each unique duplicate pair
qc_summary <- df_corrected %>%
  filter(Sample_Type != "Blank") %>%
  group_by(Plate, Sample_Type, Group, Mouse_ID, True_Conc_pg_mL) %>%
  summarize(
    Mean_OD = mean(Corrected_OD),
    SD_OD = sd(Corrected_OD),
    CV_Percent = (SD_OD / Mean_OD) * 100,
    .groups = 'drop'
  )

# Flag pairs violating acceptable laboratory variance threshold (>15% CV)
failed_qc <- qc_summary %>% filter(CV_Percent > 15)
cat(sprintf("Result: Identified %d replicate pairs exceeding the 15%% CV benchmark.\n", nrow(failed_qc)))


# ==============================================================================
# PART 3: 4PL Curve Fitting & Sample Back-Calculation
# ==============================================================================
cat("\n--- [3/4] Modeling 4PL Standard Curves & Back-Calculating Concentrations ---\n")

standards <- df_corrected %>% filter(Sample_Type == "Standard")
estimated_samples <- data.frame()
curve_lines <- data.frame() # To store smooth plotting lines
unique_plates <- unique(df_corrected$Plate)

for (p in unique_plates) {
  plate_standards <- standards %>% filter(Plate == p)
  
  # Fit individual LL.4 (4-parameter log-logistic) model per plate [drc]
  model_4pl <- drm(Corrected_OD ~ True_Conc_pg_mL, data = plate_standards, fct = LL.4())
  
  # Generate continuous curve tracking coordinates for plotting purposes
  seq_conc <- exp(seq(log(7.8), log(500), length.out = 100))
  pred_od <- predict(model_4pl, newdata = data.frame(True_Conc_pg_mL = seq_conc))
  curve_lines <- rbind(curve_lines, data.frame(Plate = p, True_Conc_pg_mL = seq_conc, Corrected_OD = pred_od))
  
  # Isolate unknown experimental mouse samples on this plate
  plate_experimental <- df_corrected %>% 
    filter(Plate == p, Sample_Type == "Experimental") %>%
    group_by(Plate, Group, Mouse_ID, True_Conc_pg_mL) %>%
    summarize(Mean_OD = mean(Corrected_OD), .groups = 'drop')
  
  # Interpolate concentration values back from measured optical density [drc]
  predicted_concs <- ED(model_4pl, plate_experimental$Mean_OD, type = "absolute", display = FALSE)[, 1]
  
  plate_experimental$Estimated_Conc_pg_mL <- predicted_concs
  estimated_samples <- rbind(estimated_samples, plate_experimental)
}

# Clean structural math floor boundaries near zero concentrations
estimated_samples <- estimated_samples %>%
  mutate(Estimated_Conc_pg_mL = ifelse(Estimated_Conc_pg_mL < 0, 0, Estimated_Conc_pg_mL))


# ==============================================================================
# PART 4: Statistical Testing via Mixed-Effects Models
# ==============================================================================
cat("\n--- [4/4] Computing Mixed-Effects Model Statistics ---\n")

# Order factors to keep Healthy Control as the default reference baseline group
estimated_samples$Group <- factor(
  estimated_samples$Group, 
  levels = c('Healthy Control', 'Disease/Inflammation', 'Low-dose Treatment', 'High-dose Treatment')
)

# Fit linear mixed model: Group is Fixed Effect, Plate is Random Intercept Effect
mixed_model <- lmer(Estimated_Conc_pg_mL ~ Group + (1 | Plate), data = estimated_samples)

print(summary(mixed_model))
print(anova(mixed_model))


# ==============================================================================
# PART 5: Visualizations & Plot Generation
# ==============================================================================
cat("\n--- Generating Visualization Panels ---\n")

# Plot A: Standard Curves Optimization Check
plot_curves <- ggplot() +
  geom_point(data = standards, aes(x = True_Conc_pg_mL, y = Corrected_OD, color = Plate), alpha = 0.5, size = 1.8) +
  geom_line(data = curve_lines, aes(x = True_Conc_pg_mL, y = Corrected_OD, color = Plate), linewidth = 1) +
  scale_x_log10() +
  labs(
    title = "A: 4PL ELISA Standard Calibration Curves",
    subtitle = "Separate curves fitted per plate to neutralize batch drift",
    x = "TNF-α Standard Concentration (pg/mL, Log Scale)",
    y = "Corrected Optical Density (OD450)"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

# Plot B: Technical Quality Control Dashboard
plot_qc <- ggplot(qc_summary, aes(x = Plate, y = CV_Percent, fill = Sample_Type)) +
  geom_boxplot(alpha = 0.7, outlier.size = 1.5) +
  geom_hline(yintercept = 15, linetype = "dashed", color = "red", linewidth = 0.8) +
  labs(
    title = "B: Replicate Technical Error Analysis",
    subtitle = "Red dashed marker highlights the 15% maximum CV benchmark",
    x = "ELISA Running Plate Reference",
    y = "Technical Duplicate Variation (%CV)"
  ) +
  theme_minimal() +
  scale_fill_brewer(palette = "Pastel1") +
  theme(legend.position = "bottom")

# Plot C: Final Calculated Biological Concentrations Group Analysis
plot_biology <- ggplot(estimated_samples, aes(x = Group, y = Estimated_Conc_pg_mL, fill = Group)) +
  geom_boxplot(alpha = 0.4, outlier.shape = NA) +
  geom_jitter(aes(shape = Plate), width = 0.18, size = 2.5, alpha = 0.85) +
  labs(
    title = "C: Therapeutic Impact on In-Vivo Mouse TNF-α Levels",
    subtitle = "Individual animal values back-calculated from 4PL regressions",
    x = "Experimental Evaluation Arm",
    y = "Estimated TNF-α Concentration (pg/mL)",
    shape = "Plate Tracking"
  ) +
  theme_minimal() +
  scale_fill_brewer(palette = "Set2") +
  theme(
    axis.text.x = element_text(angle = 12, hjust = 1),
    legend.position = "bottom"
  ) +
  guides(fill = "none") # Hide redundant group fill legend

# ==============================================================================
# PART 6: Best-Practice All-vs-All Group Comparisons (Tukey's HSD)
# ==============================================================================
cat("\n--- Running All-Versus-All Group Comparisons (Tukey HSD) ---\n")

# 1. Fit a standard linear model since plate variation was zero
anova_model <- aov(Estimated_Conc_pg_mL ~ Group, data = estimated_samples)

# 2. Run Tukey's HSD post-hoc test for all pairwise combinations
tukey_results <- TukeyHSD(anova_model)

# 3. Print the results clearly to the console
print(tukey_results)


# Define exactly which comparisons you want to draw brackets for
my_comparisons <- list(
  c("Healthy Control", "Disease/Inflammation"),
  c("Disease/Inflammation", "Low-dose Treatment"),
  c("Disease/Inflammation", "High-dose Treatment"),
  c("Healthy Control", "High-dose Treatment")
)

# Plot C: Final Calculated Biological Concentrations with explicit p-values
plot_biology <- ggplot(estimated_samples, aes(x = Group, y = Estimated_Conc_pg_mL, fill = Group)) +
  geom_boxplot(alpha = 0.4, outlier.shape = NA) +
  geom_jitter(aes(shape = Plate), width = 0.18, size = 2.5, alpha = 0.85) +
  
  # This line automatically calculates and draws the brackets with exact p-values
  stat_compare_means(comparisons = my_comparisons, 
                     method = "t.test", 
                     label = "p.format", 
                     step_increase = 0.12, 
                     size = 3.8) +
  
  labs(
    title = "C: Therapeutic Impact on In-Vivo Mouse TNF-α Levels",
    subtitle = "Individual animal values with pairwise adjusted p-values",
    x = "Experimental Evaluation Arm",
    y = "Estimated TNF-α Concentration (pg/mL)",
    shape = "Plate Tracking"
  ) +
  theme_minimal() +
  scale_fill_brewer(palette = "Set2") +
  theme(
    axis.text.x = element_text(angle = 12, hjust = 1),
    legend.position = "bottom"
  ) +
  guides(fill = "none")

# Display the plot
print(plot_biology)


# Stitch plots together using patchwork layout management
# Plots A & B share the top row; Plot C spans the entire bottom section
final_dashboard <- (plot_curves + plot_qc) / plot_biology

# Display dashboard in the R active graphics device window
print(final_dashboard)

# Export high-resolution file to your directory
ggsave("elisa_complete_analysis_dashboard.png", plot = final_dashboard, width = 11, height = 9, dpi = 300)
cat("\nPipeline evaluation finished. High-resolution figure saved as 'elisa_complete_analysis_dashboard.png'\n")
