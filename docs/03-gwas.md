---
editor_options: 
  markdown: 
    wrap: 72
---

# GWAS {#sec-gw}

1.  **Conduct GWAS Analysis and Collect All Results**

2.  **Preview the GWAS Results**

## **GWAS Analysis**

We used the [**GAPIT**](https://zzlab.net/GAPIT/gapit_help_document.pdf)
R package, which offers multiple GWAS models, to conduct GWAS analyses.

Four separate GWAS results were stored in files GWAS1 through GWAS4. In
each analysis, we applied five models simultaneously: GLM, MLM, MLMM,
FarmCPU, and BLINK.

<div class="rmdnote">
<p>Genotypic dataset (Geno.rds): 327 samples and 1.9 million SNPs.
Phenotypic dataset (Pheno.rds): 327 samples and 3 traits.</p>
</div>

<details>

<summary>➡️ **Working on R (RStudio)**</summary>

``` r
##### R CODES #####

# Set the working directory for initial setup
setwd("...change...")

# Load genotype and phenotype data
Geno = readRDS("Geno.rds")
Pheno = readRDS("Pheno.rds") 
# 4 traits: Control_AG_14DAS; Control_AG_7DAS; AG_14DAS; AG_7DAS

# Load the GAPIT library for GWAS analysis
library(GAPIT)

# Set the working directory for GWAS_PC2 analysis
setwd("...change.../GWAS_PC2")

# Perform GWAS analysis using different models and 2 principal components
GAPIT = GAPIT(
  Y = Pheno,
  G = Geno,
  model = c("GLM", "MLM", "MLMM", "FarmCPU", "Blink"),
  PCA.total = 2,
  Multiple_analysis = TRUE)

# Set the working directory for GWAS_PC3 analysis
setwd("...change.../GWAS_PC3")

# Perform GWAS analysis using different models and 3 principal components
GAPIT = GAPIT(
  Y = Pheno,
  G = Geno,
  model = c("GLM", "MLM", "MLMM", "FarmCPU", "Blink"),
  PCA.total = 3,
  Multiple_analysis = TRUE)

# Set the working directory for GWAS_PC4 analysis
setwd("...change.../GWAS_PC4")

# Perform GWAS analysis using different models and 4 principal components
GAPIT = GAPIT(
  Y = Pheno,
  G = Geno,
  model = c("GLM", "MLM", "MLMM", "FarmCPU", "Blink"),
  PCA.total = 4,
  Multiple_analysis = TRUE)

# Set the working directory for GWAS_PC5 analysis
setwd("...change.../GWAS_PC5")

# Perform GWAS analysis using different models and 5 principal components
GAPIT = GAPIT(
  Y = Pheno,
  G = Geno,
  model = c("GLM", "MLM", "MLMM", "FarmCPU", "Blink"),
  PCA.total = 5,
  Multiple_analysis = TRUE)
```

</details>

## **Preview the GWAS Results**

-   **GAPIT.Association.Filter_GWAS_results.csv**: All of the
    significant SNP markers from the selected models and phenotypic
    traits.

-   **GAPIT.Association.Manhattan_Geno.Model.Trait.pdf**: Genome-wide
    Manhattan plots for each model and specific phenotypic trait.

Here are the results: We only included significant SNP markers from
**FarmCPU** and **BLINK** (the newest models) that were associated with
**Logit_Control_AG_14DAS** (log-transformed data of the difference in
seed germination rates between control and treatment groups at 14 days).

| File name     | *GWAS_PC1* | *GWAS_PC2* | *GWAS_PC3* | *GWAS_PC4* |
|---------------|------------|------------|------------|------------|
| Number of PCs | 2          | 3          | 4          | 5          |
| FarmCPU       |            |            |            |            |
| BLINK         |            |            |            |            |
