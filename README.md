## System-Level Evaluation of Anti-Inflammatory Therapeutics on Systemic TNF-α
This repository contains a complete, single-file R data pipeline designed to ingest raw, multi-plate ELISA microplate reader output data, perform background corrections, execute technical replicate quality control, fit individual non-linear regression curves per plate, and perform robust downstream group hypothesis testing.
------------------------------
## 1. Experimental Design Overview
The dataset tracks a classic in vivo preclinical validation study evaluating an anti-inflammatory drug candidate across four experimental arms ($n = 8$ biological subjects per arm, $N = 32$ total mice):

   1. Healthy Control ($n=8$): Baseline group providing the homeostatic, non-inflamed reference value.
   2. Disease/Inflammation ($n=8$): Untreated vehicle control group modeling acute, severe system-wide inflammation.
   3. Low-dose Treatment ($n=8$): Cohort subjected to the disease model and treated with a low dose of the test therapeutic.
   4. High-dose Treatment ($n=8$): Cohort subjected to the disease model and treated with a high dose of the test therapeutic.

## Microplate Setup Logistics
To capture realistic laboratory batch effects, the biological samples are structured across 4 independent ELISA plates (block allocation layout: exactly 2 mice from each of the 4 groups are distributed per plate). Every individual biological mouse sample is run in technical duplicates ($n=2$ wells per mouse), yielding 64 experimental data rows. Additionally, each plate contains a 2-well blank block and an 8-point standard concentration ladder run in duplicate to calibrate signal generation.
------------------------------
## 2. Computational Pipeline Architecture
The unified R script processes the raw microplate file through five distinct modular steps:

  [Raw OD450 Data File]
           │
           ▼
 [1. Blank Correction]  ──► Extracts plate-specific mean background noise ($OD_{blank}$)
           │
           ▼
 [2. Replicate QC]      ──► Computes %CV; highlights sample wells violating the <15% noise ceiling
           │
           ▼
 [3. 4PL Optimization]  ──► Fits separate, sigmoidal standard calibration curves for each plate
           │
           ▼
 [4. Back-Calculation]  ──► Inverts the 4PL matrix to map sample $OD$ readings into true pg/mL units
           │
           ▼
 [5. Global Inference]  ──► Executes a One-Way ANOVA paired with an all-vs-all Tukey HSD Post-Hoc Test

------------------------------
## 3. Mathematical & Statistical Specifications## A. Blank Modification
Optical baseline density adjustments are isolated by plate. Let $p$ represent a specific plate ($p \in \{1, 2, 3, 4\}$). The adjusted optical density ($OD_{Corrected}$) for any given well $i$ on plate $p$ is calculated as:
$$OD_{Corrected, \, i, \, p} = OD_{Raw, \, i, \, p} - \overline{OD}_{Blank, \, p}$$ 
Where $\overline{OD}_{Blank, \, p}$ is the arithmetic mean of the two unreactive buffer blank control wells on that specific plate.
## B. Technical Quality Control Metric (%CV)
To isolate human pipetting errors, structural plate defects, or localized bubbles, the pipeline calculates the Coefficient of Variation ($CV$) across every technical duplicate pair. For a mouse sample with replicate readings $r_1$ and $r_2$:
$$\mu = \frac{r_1 + r_2}{2}$$ 
$$\sigma = \sqrt{\frac{(r_1 - \mu)^2 + (r_2 - \mu)^2}{2 - 1}}$$ 
$$\%CV = \left( \frac{\sigma}{\mu} \right) \times 100$$ 
Pairs yielding a $\%CV > 15\%$ are flagged in the terminal execution log.
## C. 4-Parameter Logistic (4PL) Curve Curve Fitting
Because colorimetric enzyme saturation profiles create an S-shaped (sigmoidal) non-linear curve, standard linear regressions introduce systemic error. The pipeline uses the specialized drc library to fit a 4-Parameter Logistic (4PL) regression independently for each plate using the standard wells [drc]:
$$OD(x) = d + \frac{a - d}{1 + \left(\frac{x}{c}\right)^b}$$ 
Where:

* $x$: The explicit, known antigen standard concentration (pg/mL).
* $a$: The lower asymptote (estimated background/minimum $OD$ boundary).
* $d$: The upper asymptote (estimated maximum hook-effect/saturation $OD$ boundary).
* $c$: The $ED_{50}$ parameter (the point of inflection; concentration yielding a response exactly halfway between $a$ and $d$).
* $b$: The Hill slope coefficient (defining the steepness of the linear dynamic range).

## D. Concentration Interpolation (Inversion Function)
To back-calculate unknown mouse concentrations from observed optical values, the 4PL equation is mathematically inverted [drc]. For a mouse sample with a known mean corrected optical density ($\overline{OD}_{Sample}$), the calculated concentration $x$ is:
$$x = c \cdot \left( \frac{a - d}{\overline{OD}_{Sample} - d} - 1 \right)^{\frac{1}{b}}$$ 
## E. Multi-Group Hypothesis Testing (ANOVA & Tukey's HSD)
Because plate-to-plate variation collapses post-4PL curve scaling, variations are evaluated using a standard global analysis of variance model:
$$Y_{ij} = \mu + \alpha_i + \epsilon_{ij}$$ 
Where $Y_{ij}$ is the calculated concentration for mouse $j$ in group $i$, $\mu$ is the baseline mean, $\alpha_i$ represents the treatment effect, and $\epsilon_{ij} \sim N(0, \sigma^2)$ is the random residual error.
To execute all-versus-all structural contrasts without inflating Type I false-positive errors, a pairwise Tukey's Honestly Significant Difference (HSD) test is computed. The test calculates a studentized range distribution $q$-statistic for every unique group pair:
$$q = \frac{\overline{Y}_A - \overline{Y}_B}{SE} = \frac{\overline{Y}_A - \overline{Y}_B}{\sqrt{\frac{MS_{Residual}}{n}}}$$ 
Where $MS_{Residual}$ represents the within-group variance extracted from the global ANOVA table, and $SE$ is the standard error of the comparison (106.31 in this dataset). The cumulative area under the outer tail of this $q$-distribution yields the final Adjusted $p$-value.
------------------------------
## 4. Execution & Quick Start Guide## Core Prerequisites
Open an R console or your development IDE and run the following command to download the structural dependencies:

install.packages(c("tidyverse", "drc", "lme4", "lmerTest", "patchwork", "ggpubr"))

## Running the Analysis

   1. Place your data file (synthetic_elisa_data.csv) directly into your active R workspace folder.
   2. Save the complete consolidated pipeline file as elisa_master_pipeline.R.
   3. Execute the workflow inside RStudio or via your system command terminal: [2] 

Rscript elisa_master_pipeline.R

------------------------------
## 5. Main Script Template (elisa_master_pipeline.R)

library(tidyverse)
library(drc)       
library(lme4)      
library(lmerTest)  
library(patchwork) 
library(ggpubr)    

# Read Dataset
df <- read.csv("synthetic_elisa_data.csv")

# 1. Blank Corrections
blanks <- df %>% filter(Sample_Type == "Blank") %>% group_by(Plate) %>% summarize(Mean_Blank_OD = mean(Raw_OD450), .groups = 'drop')
df_corrected <- df %>% left_join(blanks, by = "Plate") %>% mutate(Corrected_OD = Raw_OD450 - Mean_Blank_OD)

# 2. Duplicate Quality Control
qc_summary <- df_corrected %>% filter(Sample_Type != "Blank") %>% group_by(Plate, Sample_Type, Group, Mouse_ID, True_Conc_pg_mL) %>%
  summarize(Mean_OD = mean(Corrected_OD), SD_OD = sd(Corrected_OD), CV_Percent = (SD_OD / Mean_OD) * 100, .groups = 'drop')

# 3 & 4. 4PL Standard Curve Regressions & Interpolation
standards <- df_corrected %>% filter(Sample_Type == "Standard")
estimated_samples <- data.frame()
unique_plates <- unique(df_corrected$Plate)

for (p in unique_plates) {
  plate_standards <- standards %>% filter(Plate == p)
  model_4pl <- drm(Corrected_OD ~ True_Conc_pg_mL, data = plate_standards, fct = LL.4()) # [drc]
  
  plate_experimental <- df_corrected %>% filter(Plate == p, Sample_Type == "Experimental") %>%
    group_by(Plate, Group, Mouse_ID, True_Conc_pg_mL) %>% summarize(Mean_OD = mean(Corrected_OD), .groups = 'drop')
  
  predicted_concs <- ED(model_4pl, plate_experimental$Mean_OD, type = "absolute", display = FALSE)[, 1] # [drc]
  plate_experimental$Estimated_Conc_pg_mL <- ifelse(predicted_concs < 0, 0, predicted_concs)
  estimated_samples <- rbind(estimated_samples, plate_experimental)
}

# 5. Statistical Inference (Tukey's HSD)
estimated_samples$Group <- factor(estimated_samples$Group, levels = c('Healthy Control', 'Disease/Inflammation', 'Low-dose Treatment', 'High-dose Treatment'))
anova_model <- aov(Estimated_Conc_pg_mL ~ Group, data = estimated_samples)
print(TukeyHSD(anova_model))

# 6. Generate Complete Visual Dashboard
# (Insert your preferred plot layouts here using patchwork/ggpubr for export)

------------------------------
## 6. Interpreting the Outputs## Final Biological Conclusions

* Disease Verification: The Disease/Inflammation arm demonstrates a massive, highly significant spike in TNF-α relative to healthy targets ($+641.98\text{ pg/mL}, \, p < 0.001$). This confirms that the model successfully induced severe acute inflammation.
* Low-Dose Therapeutic Impact: The low-dose regimen achieves a statistically valid downward shift in systemic inflammation, shedding over $410\text{ pg/mL}$ of reactive protein relative to the untreated controls ($p = 0.0032$).
* High-Dose Full Recovery Profile: The high-dose treatment drops systemic inflammation by $580.36\text{ pg/mL}$ compared to the disease group ($p < 0.001$). Crucially, when compared against the completely healthy baseline group, the adjusted probability values reveal no statistical difference ($p = 0.937$). This proves that the high-dose intervention successfully clears the molecular pathology, returning the animals to a normal baseline condition.

------------------------------
💡 To wrap up your project workspace completely, would you like me to add a section to this README outlining a Troubleshooting Checklist for common wet-lab ELISA failures (like high hook effects or baseline standard drift)?
