# 🧬 System-Level Evaluation of Anti-Inflammatory Therapeutics on Systemic TNF-α

An end-to-end R-based analysis pipeline for processing ELISA optical density (OD450) data, performing plate-specific quality control, fitting nonlinear 4-parameter logistic (4PL) standard curves, estimating TNF-α concentrations, and evaluating treatment effects statistically.

> **Research question:** Does anti-inflammatory treatment reduce systemic TNF-α levels in a mouse model of inflammation?

---

## 📋 Table of Contents

* [1. Experimental Design](#1-experimental-design)
* [2. ELISA Microplate Design](#2-elisa-microplate-design)
* [3. Computational Analysis Pipeline](#3-computational-analysis-pipeline)
* [4. Mathematical and Statistical Methods](#4-mathematical-and-statistical-methods)

  * [4.1 Blank Correction](#41-blank-correction)
  * [4.2 Technical Replicate QC](#42-technical-replicate-quality-control)
  * [4.3 4-Parameter Logistic Regression](#43-4-parameter-logistic-4pl-regression)
  * [4.4 Concentration Back-Calculation](#44-concentration-back-calculation)
  * [4.5 Statistical Inference](#45-statistical-inference)
* [5. Installation and Quick Start](#5-installation-and-quick-start)
* [6. Main Analysis Workflow](#6-main-analysis-workflow)
* [7. Results and Interpretation](#7-results-and-interpretation)
* [8. Troubleshooting Guide](#8-troubleshooting-guide)
* [9. R Packages](#9-r-packages)

---

# 1. Experimental Design

The dataset simulates a preclinical in vivo study evaluating the efficacy of an anti-inflammatory treatment across four experimental groups.

| Experimental Group       | Description                         | Sample Size |
| ------------------------ | ----------------------------------- | ----------: |
| **Healthy Control**      | Non-inflamed baseline reference     |     $n = 8$ |
| **Disease/Inflammation** | Untreated disease model             |     $n = 8$ |
| **Low-dose Treatment**   | Disease model + low-dose treatment  |     $n = 8$ |
| **High-dose Treatment**  | Disease model + high-dose treatment |     $n = 8$ |

**Total biological subjects:**

$$
N = 32
$$

The primary outcome is **TNF-α concentration (pg/mL)** estimated from ELISA optical density measurements using plate-specific 4PL standard curves.

---

# 2. ELISA Microplate Design

To simulate realistic laboratory batch effects, biological samples are distributed across **four independent ELISA plates**.

Each plate contains:

* 2 mice from each experimental group
* 8 biological samples total
* 2 technical replicates per biological sample
* 2 blank wells
* 8-point TNF-α standard curve
* Duplicate standard measurements

### Plate Allocation

| Plate     | Healthy | Disease | Low Dose | High Dose | Total Mice |
| --------- | ------: | ------: | -------: | --------: | ---------: |
| Plate 1   |       2 |       2 |        2 |         2 |          8 |
| Plate 2   |       2 |       2 |        2 |         2 |          8 |
| Plate 3   |       2 |       2 |        2 |         2 |          8 |
| Plate 4   |       2 |       2 |        2 |         2 |          8 |
| **Total** |   **8** |   **8** |    **8** |     **8** |     **32** |

Each mouse is measured in technical duplicate wells:

$$
32 \text{ mice} \times 2 \text{ technical replicates} = 64 \text{ experimental wells}
$$

The experimental data are accompanied by blank wells and duplicate standard curves on each plate.

---

# 3. Computational Analysis Pipeline

The complete workflow transforms raw ELISA optical density readings into estimated TNF-α concentrations and statistical conclusions.

```text
                 Raw ELISA OD450 Data
                         │
                         ▼
              ┌──────────────────────┐
              │ 1. Blank Correction  │
              │ Plate-specific       │
              │ background removal   │
              └──────────┬───────────┘
                         │
                         ▼
              ┌──────────────────────┐
              │ 2. Replicate QC      │
              │ Calculate technical  │
              │ replicate %CV        │
              └──────────┬───────────┘
                         │
                         ▼
              ┌──────────────────────┐
              │ 3. 4PL Calibration   │
              │ Fit plate-specific   │
              │ nonlinear curves     │
              └──────────┬───────────┘
                         │
                         ▼
              ┌──────────────────────┐
              │ 4. Back-Calculation  │
              │ Convert OD450 into   │
              │ TNF-α concentration  │
              └──────────┬───────────┘
                         │
                         ▼
              ┌──────────────────────┐
              │ 5. Statistical       │
              │ Inference            │
              │ ANOVA + Tukey HSD    │
              └──────────┬───────────┘
                         │
                         ▼
                  Biological Insight
```

---

# 4. Mathematical and Statistical Methods

## 4.1 Blank Correction

Background optical density is estimated independently for each ELISA plate.

Let $p$ represent a specific plate:

$$
p \in {1,2,3,4}
$$

The corrected optical density for well $i$ on plate $p$ is:

$$
OD_{\text{Corrected},i,p}
=========================

## OD_{\text{Raw},i,p}

\overline{OD}_{\text{Blank},p}
$$

where:

* $OD_{\text{Raw},i,p}$ = raw OD450 measurement
* $\overline{OD}_{\text{Blank},p}$ = mean OD450 of the blank wells on plate $p$

This removes plate-specific background signal before standard curve fitting.

---

## 4.2 Technical Replicate Quality Control

Each biological sample is measured in duplicate.

For replicate measurements $r_1$ and $r_2$, the mean is:

$$
\mu =
\frac{r_1+r_2}{2}
$$

The sample standard deviation is:

$$
s =
\sqrt{
\frac{
(r_1-\mu)^2+(r_2-\mu)^2
}{
2-1
}
}
$$

The coefficient of variation is:

$$
%CV =
\left(
\frac{s}{\mu}
\right)
\times 100
$$

Technical replicate pairs exceeding the predefined threshold are flagged:

$$
%CV > 15%
$$

The flagged samples are retained for investigation rather than automatically removed from the dataset.

---

## 4.3 4-Parameter Logistic (4PL) Regression

ELISA assays typically produce a nonlinear, sigmoidal relationship between analyte concentration and optical density.

Therefore, a 4-parameter logistic model is fitted separately to each ELISA plate using the `drc` package.

The model is:

<div align="center">

$$
OD(x)
=====

d+
\frac{a-d}
{1+\left(\frac{x}{c}\right)^b}
$$

</div>

where:

* $x$ = known TNF-α concentration (pg/mL)
* $a$ = lower asymptote
* $d$ = upper asymptote
* $c$ = inflection point / $ED_{50}$
* $b$ = slope parameter

Separate curves are fitted for each plate to account for potential plate-specific differences in assay response.

---

## 4.4 Concentration Back-Calculation

Once the 4PL curve has been fitted, the measured OD of each unknown sample can be converted back into an estimated TNF-α concentration.

Given a sample mean corrected OD, $\overline{OD}_{\text{Sample}}$, the inverse relationship is:

$$
x
=

c
\left[
\frac{a-d}
{\overline{OD}_{\text{Sample}}-d}
-1
\right]^{1/b}
$$

where $x$ represents the estimated TNF-α concentration in pg/mL.

The `drc::ED()` function is used to perform this back-calculation programmatically.

Negative numerical estimates are constrained to zero:

```r
Estimated_Conc_pg_mL =
  ifelse(Estimated_Conc_pg_mL < 0, 0, Estimated_Conc_pg_mL)
```

> **Note:** In real ELISA analysis, samples below the assay's lower limit of quantification should generally be handled according to the assay validation protocol rather than automatically interpreted as true zero concentrations.

---

## 4.5 Statistical Inference

### One-Way ANOVA

After concentration estimation, group differences are initially evaluated using a one-way analysis of variance (ANOVA).

The model is:

$$
Y_{ij}
======

\mu
+
\alpha_i
+
\epsilon_{ij}
$$

where:

* $Y_{ij}$ = estimated TNF-α concentration for mouse $j$ in group $i$
* $\mu$ = overall mean
* $\alpha_i$ = effect of experimental group $i$
* $\epsilon_{ij}$ = residual error

The null hypothesis is:

$$
H_0:
\mu_1=\mu_2=\mu_3=\mu_4
$$

The alternative hypothesis is:

$$
H_A:
\text{At least one group mean differs}
$$

---

### Tukey's HSD Post-Hoc Test

When the global ANOVA indicates a significant group effect, Tukey's Honestly Significant Difference (HSD) test is used to perform all pairwise comparisons while controlling the family-wise error rate.

The analysis compares:

* Healthy Control vs Disease/Inflammation
* Healthy Control vs Low-dose Treatment
* Healthy Control vs High-dose Treatment
* Disease/Inflammation vs Low-dose Treatment
* Disease/Inflammation vs High-dose Treatment
* Low-dose Treatment vs High-dose Treatment

The Tukey procedure uses the studentized range distribution to calculate multiplicity-adjusted $p$-values.

---

# 5. Installation and Quick Start

## Prerequisites

Install R and RStudio, then install the required packages:

```r
install.packages(c(
  "tidyverse",
  "drc",
  "lme4",
  "lmerTest",
  "patchwork",
  "ggpubr"
))
```

## Project Structure

```text
ELISA_TNF_Alpha_Analysis/
│
├── data/
│   └── MouseTNF_elisa_data.csv
│
├── scripts/
│   └── elisa_master_pipeline.R
│
├── results/
│   └── elisa_complete_analysis_dashboard.png
│
├── README.md
│
└── .gitignore
```

## Running the Analysis

Clone the repository:

```bash
git clone https://github.com/YOUR_USERNAME/ELISA_TNF_Alpha_Analysis.git
```

Navigate to the project:

```bash
cd ELISA_TNF_Alpha_Analysis
```

Run the analysis:

```bash
Rscript scripts/elisa_master_pipeline.R
```

Alternatively, open `elisa_master_pipeline.R` in RStudio and run the script interactively.

---

# 6. Main Analysis Workflow

The core analysis is implemented in:

```text
scripts/elisa_master_pipeline.R
```

### Load Packages

```r
library(tidyverse)
library(drc)
library(lme4)
library(lmerTest)
library(patchwork)
library(ggpubr)
```

### Read Data

```r
df <- read.csv("data/MouseTNF_elisa_data.csv")
```

### Blank Correction

```r
blanks <- df %>%
  filter(Sample_Type == "Blank") %>%
  group_by(Plate) %>%
  summarize(
    Mean_Blank_OD = mean(Raw_OD450),
    .groups = "drop"
  )

df_corrected <- df %>%
  left_join(blanks, by = "Plate") %>%
  mutate(
    Corrected_OD = Raw_OD450 - Mean_Blank_OD
  )
```

### Technical Replicate QC

```r
qc_summary <- df_corrected %>%
  filter(Sample_Type != "Blank") %>%
  group_by(
    Plate,
    Sample_Type,
    Group,
    Mouse_ID,
    True_Conc_pg_mL
  ) %>%
  summarize(
    Mean_OD = mean(Corrected_OD),
    SD_OD = sd(Corrected_OD),
    CV_Percent = (SD_OD / Mean_OD) * 100,
    .groups = "drop"
  )
```

### 4PL Standard Curve Fitting

```r
model_4pl <- drm(
  Corrected_OD ~ True_Conc_pg_mL,
  data = plate_standards,
  fct = LL.4()
)
```

### Back-Calculate Unknown Concentrations

```r
predicted_concs <- ED(
  model_4pl,
  plate_experimental$Mean_OD,
  type = "absolute",
  display = FALSE
)[, 1]
```

### Statistical Testing

```r
anova_model <- aov(
  Estimated_Conc_pg_mL ~ Group,
  data = estimated_samples
)

TukeyHSD(anova_model)
```

---

# 7. Results and Interpretation

## Overall Treatment Effect

The analysis identified a significant difference in TNF-α concentration across the four experimental groups:

$$
F(3,28)=14.8,
\qquad
p<0.001
$$

This indicates that at least one experimental group differed significantly from the others.

### Pairwise Comparisons

| Comparison            | Mean Difference (pg/mL) | Adjusted $p$-value | Interpretation        |
| --------------------- | ----------------------: | -----------------: | --------------------- |
| Disease vs Healthy    |                 +641.98 |            < 0.001 | Significant increase  |
| Low Dose vs Healthy   |                 +231.39 |              0.155 | Not significant       |
| High Dose vs Healthy  |                  +61.62 |              0.937 | Not significant       |
| Low Dose vs Disease   |                 −410.59 |              0.003 | Significant reduction |
| High Dose vs Disease  |                 −580.36 |            < 0.001 | Significant reduction |
| High Dose vs Low Dose |                 −169.77 |              0.397 | Not significant       |

### Biological Interpretation

The results support the following interpretation:

1. **Disease induction:** The disease/inflammation group showed substantially higher TNF-α concentrations than healthy controls, indicating successful induction of an inflammatory phenotype in the simulated dataset.

2. **Low-dose treatment:** Low-dose treatment significantly reduced TNF-α concentrations relative to untreated disease.

3. **High-dose treatment:** High-dose treatment produced the largest numerical reduction in TNF-α relative to untreated disease.

4. **Recovery toward baseline:** TNF-α concentrations in the high-dose group were statistically indistinguishable from healthy controls ($p=0.937$).

> **Important:** Statistical similarity to healthy controls does not by itself prove complete biological recovery or clearance of molecular pathology. It indicates that, for the TNF-α outcome measured in this experiment, the two groups could not be statistically distinguished.

---

# 8. Troubleshooting Guide

## A. Poor 4PL Fit or Saturation

### What happens?

The highest standard concentrations flatten or deviate substantially from the expected sigmoidal response.

### Potential causes

* Signal saturation
* Excessively concentrated standards
* Incorrect dilution series
* Substrate incubation time
* Assay-specific hook effects

### Pipeline symptoms

The `drm()` function may produce:

* Convergence warnings
* Poor parameter estimates
* `Inf` values
* `NaN` values

### Possible solutions

**Laboratory:**

* Increase dilution of high-concentration standards.
* Reduce substrate incubation time.
* Verify standard preparation.

**Computational:**

Consider whether a 5-parameter logistic model is more appropriate:

```r
drm(
  Corrected_OD ~ True_Conc_pg_mL,
  data = plate_standards,
  fct = LL.5()
)
```

The 5PL model adds an asymmetry parameter and may better represent asymmetric calibration curves.

---

## B. High Technical Replicate CV

### What happens?

A large number of duplicate measurements exceed the predefined 15% CV threshold.

### Potential causes

* Pipetting variability
* Edge effects
* Inconsistent washing
* Air bubbles
* Inadequate mixing

### Possible solutions

**Laboratory:**

* Improve pipetting consistency.
* Use calibrated pipettes.
* Minimize edge effects.
* Standardize washing procedures.

**Computational:**

Flag high-CV samples for investigation rather than automatically deleting observations.

> **Best practice:** Avoid automatically removing individual wells solely because they have a high CV. Investigate the raw measurements and assay quality criteria first.

---

## C. High Plate-to-Plate Variation

### What happens?

Equivalent standards produce substantially different calibration curves across plates.

### Potential causes

* Temperature differences
* Incubation time differences
* Reagent lot variation
* Timing differences
* Plate handling effects

### Computational approach

If plate-level variation remains substantial after plate-specific calibration, consider a mixed-effects model:

```r
lmer(
  Estimated_Conc_pg_mL ~ Group + (1 | Plate),
  data = estimated_samples
)
```

This models **Group** as a fixed effect while treating **Plate** as a random intercept.

---

## D. Zero or Negative Back-Calculated Concentrations

### What happens?

Some samples produce estimates near or below zero.

### Potential cause

The sample OD may fall below the lower asymptote of the calibration curve and therefore outside the validated quantitative range of the assay.

### Recommended approach

In a real experimental dataset, investigate whether the sample is:

* Below the lower limit of detection (LLOD)
* Below the lower limit of quantification (LLOQ)
* Outside the validated standard curve range

Avoid automatically interpreting a negative back-calculated concentration as a true biological zero.

For this educational dataset, negative numerical estimates are constrained to zero:

```r
Estimated_Conc_pg_mL <-
  ifelse(
    Estimated_Conc_pg_mL < 0,
    0,
    Estimated_Conc_pg_mL
  )
```

---

# 9. R Packages

This project uses the following R packages:

| Package     | Purpose                                             |
| ----------- | --------------------------------------------------- |
| `tidyverse` | Data manipulation and visualization                 |
| `drc`       | Nonlinear dose-response and 4PL curve fitting       |
| `lme4`      | Linear mixed-effects models                         |
| `lmerTest`  | Statistical inference for mixed-effects models      |
| `patchwork` | Combining multiple `ggplot2` figures                |
| `ggpubr`    | Statistical annotations and publication-ready plots |

---

## 📊 Key Takeaway

This project demonstrates a complete workflow for transforming raw ELISA optical density measurements into interpretable biological results:

> **Raw OD450 → Quality Control → Blank Correction → 4PL Calibration → Concentration Estimation → Statistical Testing → Biological Interpretation**

The analysis highlights how careful statistical methodology can connect **laboratory measurements** with **quantitative biological conclusions**.
