---
editor_options: 
  markdown: 
    wrap: 72
---

# Findings {#sec-findings}

-   **We used a subset of 311 samples for further analysis.** This
    subset was obtained by excluding 30 samples that had phenotypic data
    but lacked genotypic data, as well as 2,713 samples that had
    genotypic data but no corresponding phenotypic data.

-   **The final dataset comprises 1,972,824 SNPs (1,972K; 1.9M) and 311
    samples.** Criteria were applied to filter out certain SNPs, while
    efforts were made to retain all samples if their missing rates were
    within acceptable limits ( \< \~10%). Although some samples
    exhibited lower inbreeding coefficients (indicating higher
    heterozygosity), all 311 samples were retained in the final dataset
    for GWAS analysis.

-   **The deviations of measurement** **under 14-day (3.09 and 9.10) are
    smaller** than under 7-day (2.93 and 8.91)→ use measurements under
    14-day are better.

| Groups | Control (7-day) | Treatment (7-day) | Control (14-day) | Treatment (14-day) |
|--------|-----------------|-------------------|------------------|--------------------|
| Mean   | 96.49           | 11.88             | 97.02            | 16.98              |
| SD     | 3.09            | 9.10              | 2.93             | 8.91               |

<details>

<summary>**How to download the R packages**</summary>

### Step 1: Understand R Packages {.unnumbered}

R packages are collections of functions, data, and documentation that
extend R's capabilities. These packages are hosted on CRAN,
Bioconductor, or GitHub.

### Step 2: Install CRAN Packages {.unnumbered}

1.  **Basic Installation**: Use the `install.packages()` function to
    install packages from CRAN. For example: This command downloads and
    installs the `ggplot2` package.

    ``` r
    install.packages("ggplot2") 
    ```

    **Load the Package**: After installation, load the package using:
    Once loaded, you can use its functions.

    ``` r
    library(ggplot2) 
    ```

### Step 3: Install Bioconductor Packages {.unnumbered}

Bioconductor hosts many bioinformatics packages. To install, follow
these steps:

1.  **Install BiocManager**:

    ``` r
    install.packages("BiocManager") 
    ```

2.  **Install Bioinformatics Packages:** Use `BiocManager::install()` to
    install Bioconductor packages. Example: This installs the `ggtree`
    package.

    ``` r
    BiocManager::install("ggtree") 
    ```

3.  **Load the Package**:

    ``` r
    library(ggtree)
    ```

### Step 4: Install GitHub Packages {.unnumbered}

Some specialized packages are hosted on GitHub. To install these:

1.  **Install `remotes` or `devtools`**:

    ``` r
    install.packages("remotes") 
    ```

2.  **Install a Package from GitHub**: Use `remotes::install_github()`.
    For example: This installs the `circlize` package for circular
    visualization.

    ``` r
    remotes::install_github("TeddYenn/ShiNyP_Test")  
    ```

3.  **Load the Package**:

    ``` r
    library(ShiNyP_Test) 
    ```

### Step 5: Update R Packages {.unnumbered}

1.  To update all installed packages:

    ``` r
    update.packages(ask = FALSE) 
    ```

<!-- -->

2.  For Bioconductor packages:

    ``` r
    BiocManager::install(update = TRUE) 
    ```

### Step 6: Uninstall a Package (If Needed) {.unnumbered}

If you need to remove a package:

``` r
remove.packages("ggplot2") 
```

### Summary of Installation Methods {.unnumbered}

| Source       | Command Example                                   |
|--------------|---------------------------------------------------|
| CRAN         | `install.packages("ggplot2")`                     |
| Bioconductor | `BiocManager::install("ggtree")`                  |
| GitHub       | `remotes::install_github("TeddYenn/ShiNyP_Test")` |

</details>

------------------------------------------------------------------------

<details>

<summary>**How to download the R & RStudio**</summary>

HOW TO DOWNLOAD & USE PLINK on TAMU HPRC
