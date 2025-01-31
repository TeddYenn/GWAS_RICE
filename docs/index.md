---
title: "GWAS Tutorial"
subtitle: "For Anaerobic Germination in Rice"
author: "Yen-Hsiang (Teddy) Huang"
date: "2025-01-28"
knit: "bookdown::render_book"
site: bookdown::bookdown_site
output: bookdown::bs4_book
documentclass: book
bibliography: [book.bib]
biblio-style: apalike
link-citations: true
links-as-notes: true
colorlinks: true
github-repo: TeddYenn/GWAS_RICE
cover-image: images/cover.png
url: https://github.com/TeddYenn/GWAS_RICE
description: "A guide to GWAS analysis in rice, using the R and PLINK programs"
editor_options: 
  markdown: 
    wrap: 72
---



# Welcome to GWAS Tutorial {#sec-welcome-to-gwas-tutorial .unnumbered}

This is the [website](https://teddyenn.github.io/GWAS_RICE/) for **GWAS
analysis with R and PLINK programs**.

This work by **Yen-Hsiang (Teddy) Huang** is licensed under a
<a rel="license" href="https://www.gnu.org/licenses/gpl-3.0.html.en">GNU
General Public License</a>.

<a href="https://www.gnu.org/licenses/gpl-3.0.html.en">
![](images/index/clipboard-3960043298.png)</a>

```{=html}
<script>
  (function(i,s,o,g,r,a,m){i['GoogleAnalyticsObject']=r;i[r]=i[r]||function(){
  (i[r].q=i[r].q||[]).push(arguments)},i[r].l=1*new Date();a=s.createElement(o),
  m=s.getElementsByTagName(o)[0];a.async=1;a.src=g;m.parentNode.insertBefore(a,m)
  })(window,document,'script','https://www.google-analytics.com/analytics.js','ga');

  ga('create', 'UA-68765210-2', 'auto');
  ga('send', 'pageview');

</script>
```

# Preface {.unnumbered}

This GitHub page provides a hands-on tutorial on genome-wide association
studies (GWAS) as part of Teddy's and Lumi's research program at Texas
A&M University under the supervision of Dr. Endang Septiningsih and Dr.
Michael Thomson.

The tutorial focuses on standard workflow, codes, and resources needed
to perform a GWAS in rice and is divided into seven sections. For
absolute beginners, we’ve also included an introductory section on
[R](https://www.r-project.org/),
[RStudio](https://posit.co/download/rstudio-desktop/), and
[PLINK](https://www.cog-genomics.org/plink/). **-\> [Click
Here!](https://teddyenn.github.io/GWAS_RICE/preface.html#sec-tutorial-for-r-rstudio--plink)**

This book was written in [RStudio](http://www.rstudio.com/ide/) using
[bookdown](http://bookdown.org/). The
[website](https://teddyenn.github.io/GWAS_RICE/) is hosted via GitHub
under [TeddYenn's repository.](https://github.com/TeddYenn)

If you have any questions or suggestions, feel free to reach out via
email at [teddyenn2\@gmail.com](mailto:teddyenn2@gmail.com){.email}
(Teddy Huang).

## Outline {.unnumbered}

We begin by processing genotypic and phenotypic datasets, followed by
conducting GWAS and post-GWAS analyses, and conclude with findings and
interpretations of the results.

-   **Chapter** \@ref(sec-genotype-data) downloads the SNP data from
    SNP-Seek (3K RG), align the IDs, perform SNP QC, and generate the
    final GWAS data.
-   **Chapter** \@ref(sec-phenotype-data) explores the phenotypic data
    and generate the final GWAS data.
-   **Chapter** \@ref(sec-gw) conduct the genome-wide association study
    (GWAS) to link genotypic variation with the phenotypic traits of
    interest.
-   **Chapter** \@ref(sec-post-gwas) delve into post-GWAS analyses that
    further dissect and validate the significant associations identified
    in the GWAS.
-   **Chapter** \@ref(sec-extra) present additional exploratory analyses
    that extend beyond the primary GWAS objectives.

We provide the detailed findings and make some narratives about the
results.

-   **Chapter** \@ref(sec-findings) synthesize the key findings from the
    genotypic, phenotypic, GWAS, and post-GWAS analyses.
-   **Chapter** \@ref(sec-summary) summarize the overall outcomes of the
    study, emphasizing major insights gained and acknowledging any
    limitations encountered.

## Tutorial for R, RStudio & PLINK {#sec-tutorial-for-r-rstudio--plink .unnumbered}

<details>

<summary>[**How to download the R &
RStudio**]{style="color: white; background-color: grey; padding: 2px 4px; border-radius: 3px;"}.</summary>

<div class="rmdtip">
<p>For TAMU HPRC users, you don’t need to download RStudio. There is the
built-in RStudio can be launch online on HPRC.</p>
</div>

**Step 1: Download and Install R (Prerequisite)**

1.  Visit the official R Project website: <https://www.r-project.org/>.

    -   **Windows**: Click **Download R for Windows**, select "base,"
        and download the latest version.

    -   **MacOS**: Click **Download R for macOS** and choose the correct
        version for your system.

2.  Install R by double-clicking the downloaded installer and following
    the on-screen instructions. Use the default options unless specific
    needs arise.

**Step 2: Download RStudio**

1.  Go to the RStudio official download page:
    <https://posit.co/download/rstudio-desktop/>.

2.  Click **Download** under "RStudio Desktop - Open Source License."

3.  Select the version suitable for your operating system:

    -   **Windows**: `.exe` file

    -   **MacOS**: `.dmg` file

**Step 3: Install RStudio**

1.  Locate the downloaded file and double-click it to start the
    installation.

    -   **Windows**: Run the `.exe` installer and follow the wizard
        steps.

    -   **MacOS**: Drag the RStudio icon into the Applications folder.

2.  After installation, launch RStudio.

**Step 4: Launch RStudio**

1.  Open RStudio by clicking the shortcut created during installation or
    searching for "RStudio" in your system's application launcher.

2.  RStudio will automatically detect your R installation and link to
    it.

**Step 5: Familiarize Yourself with RStudio Interface**

RStudio has four main panels:

1.  **Console (Bottom-left)**: Where you run R commands.

2.  **Source (Top-left)**: For writing and editing scripts.

3.  **Environment/History (Top-right)**: Displays objects, variables,
    and command history.

4.  **Plots/Files/Help (Bottom-right)**: Displays plots, files, and R
    documentation.

------------------------------------------------------------------------

</details>

<details>

<summary>[**How to download the R
packages**]{style="color: white; background-color: grey; padding: 2px 4px; border-radius: 3px;"}.</summary>

**Step 1: Understand R Packages**

R packages are collections of functions, data, and documentation that
extend R's capabilities. These packages are hosted on CRAN,
Bioconductor, or GitHub.

**Step 2: Install CRAN Packages**

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

**Step 3: Install Bioconductor Packages**

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

**Step 4: Install GitHub Packages**

Some specialized packages are hosted on GitHub. To install these:

1.  **Install `remotes` or `devtools`**:

    ``` r
    install.packages("remotes") 
    ```

2.  **Install a Package from GitHub**: Use `remotes::install_github()`.
    For example: This installs the `GAPIT3` package for GWAS analysis.

    ``` r
    remotes::install_github("jiabowang/GAPIT3")  
    ```

3.  **Load the Package**:

    ``` r
    library(GAPIT3) 
    ```

**Step 5: Update R Packages**

1.  To update all installed packages:

    ``` r
    update.packages(ask = FALSE) 
    ```

<!-- -->

2.  For Bioconductor packages:

    ``` r
    BiocManager::install(update = TRUE) 
    ```

**Step 6: Uninstall a Package (If Needed)**

If you need to remove a package:

``` r
remove.packages("ggplot2") 
```

**Summary of Installation Methods**

| Source       | Command Example                               |
|--------------|-----------------------------------------------|
| CRAN         | `install.packages("ggplot2")`                 |
| Bioconductor | `BiocManager::install("ggtree")`              |
| GitHub       | `remotes::install_github("jiabowang/GAPIT3")` |

------------------------------------------------------------------------

</details>

<details>

<summary>[**How to download and use PLINK** **on TAMU
HPRC**]{style="color: white; background-color: grey; padding: 2px 4px; border-radius: 3px;"}.</summary>

**Step 1: Download the PLINK File**

1.  Visit the [PLINK website](https://www.cog-genomics.org/plink/1.9/).

2.  Download the **Linux 64-bit** version, which is compatible with the
    [FASTER](https://hprc.tamu.edu/kb/User-Guides/FASTER/Hardware/) and
    [Grace](https://hprc.tamu.edu/kb/User-Guides/Grace/Hardware/)
    system.

**Step 2: Extract and Upload the PLINK File to TAMU HPRC**

1.  Extract the downloaded PLINK file (e.g., `plink_linux_x86_64`).

2.  Upload the extracted PLINK folder to TAMU HPRC via the file transfer
    method you prefer. For example:

    -   Create an `Upload` folder in your home directory
        (`/home/<username>/Upload`) or scratch
        (`/scratch/user/<username>/Upload`).

    -   Place the PLINK folder (`plink_linux_x86_64`) inside this
        directory.

**Step 3: Open the Terminal on TAMU HPRC**

1.  Log in to the TAMU HPRC portal.

2.  Navigate to your working directory where the PLINK file was uploaded
    (e.g. `/home/<username>/Upload/plink_linux_x86_64`).

3.  Click **`>_ Open in Terminal`** to access the terminal interface.

**Step 4: Verify Your Current Directory**

1.  In the terminal, run the command after \$

    ```         
    pwd 
    ```

2.  Confirm that you are in the directory where the PLINK file is
    located. For example:

    ```         
     /home/<username>/Upload/plink_linux_x86_64 
    ```

**Step 5: Make the PLINK File Executable and Test It**

-   In the terminal, navigate to the PLINK folder if not already there:

    ```         
    cd /home/<username>/Upload/plink_linux_x86_6 
    ```

-   Run the following commands:

    ```         
    chmod +x plink
    ```

    ```         
    ./plink --help
    ```

    -   **`chmod +x plink`**: Grants the PLINK file executable
        permissions.
    -   **`./plink --help`**: Runs the PLINK executable and displays
        help information to confirm proper setup.

**Step 6: Upload Your Data and Use PLINK**

1.  Upload your data files (e.g., `.bed`, `.bim`, `.fam` files) into the
    working directory (`/home/<username>/Upload/plink_linux_x86_64`).

2.  Use PLINK commands to process your data.

**Example: Converting `.bfile` to `.vcf`**

1.  Assuming your `.bed`, `.bim`, and `.fam` files are named `data.bed`,
    `data.bim`, and `data.fam`, use the following command:

    ```         
    ./plink --bfile data --recode vcf --out data 
    ```

<!-- -->

2.  This will generate a VCF file named `data.vcf` in your current
    directory.

------------------------------------------------------------------------

</details>

## Acknowledgements {.unnumbered}

We are deeply grateful for the contributions, support, and perspectives
of individuals and organizations that have helped move this project
forward.

We would like to extend our sincere thanks to Dr. Endang Septiningsih
and Dr. Michael Thomson for hosting us at Texas A&M University. We also
greatly appreciate the support provided by the International Agriculture
Center, National Chung Hsing University, and the Higher Education Sprout
Project, Ministry of Education, Taiwan, for offering opportunities and
partial financial support for this research program.
