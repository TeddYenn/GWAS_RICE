---
editor_options: 
  markdown: 
    wrap: 72
---

# GWAS {#sec-gw}

<div class="rmdtip">
<p>In this chapter, we conduct the genome-wide association study (GWAS)
to link genotypic variation with the phenotypic traits of interest. By
systematically scanning the genome for marker-trait associations, we
identify loci that may influence the traits in question. The chapter
concludes with initial interpretations of the significant associations
discovered.</p>
</div>

1.  **Conduct GWAS Analysis and Collect All Results**

2.  **Preview the GWAS Results**

## **GWAS Analysis**

We used the [**GAPIT**](https://zzlab.net/GAPIT/gapit_help_document.pdf)
R package, which offers multiple GWAS models, to conduct GWAS analyses.

Four separate GWAS results (`PCA.total = 2, 3, 4, 5`) were stored in
files `GWAS_PC2` through `GWAS_PC5`. Create the GWAS_PC2, GWAS_PC3,
GWAS_PC4, GWAS_PC5 folder in your directory. In each analysis, we
applied five models simultaneously: GLM, MLM, MLMM, FarmCPU, and BLINK.

<div class="rmdnote">
<p>Genotypic dataset (Geno.rds): 327 samples and 1.9 million SNPs.
Phenotypic dataset (Pheno.rds): 327 samples and 3 traits.</p>
</div>

<details>

<summary>**R** ➡️ **GWAS analysis with GAPIT**</summary>

``` r
##### R CODES #####

# Set the working directory for initial setup
setwd("...change...")

# Load genotype and phenotype data
Geno = readRDS("Geno.rds")
Pheno = readRDS("Pheno.rds") 
# 3 traits: AG_14DAS; Control_AG_14DAS; Logit_Control_AG_14DAS

# Load the GAPIT library for GWAS analysis 
install.packages("devtools")
devtools::install_github("jiabowang/GAPIT", force = TRUE)
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
**Control_AG_14DAS** (the difference in seed germination rates between
control and treatment groups at 14 days).

| Filename      | *GWAS_PC2* | *GWAS_PC3* | *GWAS_PC4* | *GWAS_PC5* |
|---------------|------------|------------|------------|------------|
| Number of PCs | 2          | 3          | 4          | 5          |
| FarmCPU       | 0          | 5          | 4          | 3          |
| BLINK         | 4          | 3          | 3          | 3          |

1.  ***GWAS_PC2***

    -   BLINK: Chr02: 24380499; Chr03: 17005606; Chr03: 32521809; Chr10:
        15165249

2.  ***GWAS_PC3***

    -   FarmCPU: Chr03: 17005606; Chr07: 18136942; Chr07: 22535556;
        Chr10: 19883513; Chr11: 18075487

    -   BLINK: Chr03: 17005606; Chr03: 32521809; Chr10: 15162418

3.  ***GWAS_PC4***

    -   FarmCPU: Chr03: 17000165; Chr06: 12865843; Chr06: 27384104;
        Chr10: 19883513

    -   BLINK: Chr03: 17005606; Chr03: 32521809; Chr9: 12290639

4.  ***GWAS_PC5***

    -   FarmCPU: Chr01: 37357384; Chr03: 17005606; Chr10: 19883513

    -   BLINK: Chr03: 17005606; Chr03: 32521809; Chr10: 15165249

+---------+--------------+--------------+--------------+--------------+
| SNP Pos | Chr03:       | Chr03:       | Chr10:       | Chr10:       |
|         | 17005606     | 32521809     | 15162418     | 19883513     |
+=========+==============+==============+==============+==============+
| SNP ID  | 96208338     | 111729982    | 308650519    | 313368783    |
+---------+--------------+--------------+--------------+--------------+
| FarmCPU | 3 [PC3\~5]   | 0            | 0            | 3 [PC3, 4,   |
|         |              |              |              | 5]           |
+---------+--------------+--------------+--------------+--------------+
| BLINK   | 3 [PC3\~5]\  | 4 [PC2\~5]   | 3 [PC2, 3,   | 0            |
|         | 1\*(Chr03:   |              | 5]           |              |
|         | 17005606)    |              |              |              |
|         | [PC2]        |              |              |              |
+---------+--------------+--------------+--------------+--------------+
