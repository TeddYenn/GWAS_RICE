---
editor_options: 
  markdown: 
    wrap: 72
---

# Summary {#sec-summary}

**Under Construction...**

1.  **Align Accession IDs:** Match the accession IDs between the
    phenotypic and SNP datasets.

2.  **SNP Data Processing and Quality Control:** Perform data generation
    and quality control on the SNP data for the subset of 311 samples.

## **Download Rice SNP data from SNP-Seek (3K RG)**

Data repository for the 3K RG, hosted by IRRI:
[https://snpseekv3.irri-e-extension.com/v2/download.zul](https://snpseekv3.irri-e-extension.com/v2/_download.zul){.uri}

-   **3K RG 18mio Base SNP Dataset:** A base SNP set consisting of
    approximately 18 million SNPs was derived from the initial \~29
    million biallelic SNPs by excluding those with an excessive number
    of heterozygous calls.

-   **3K RG 404k CoreSNP Dataset:** The Core SNP set was derived from
    the filtered SNP set through a two-step LD pruning procedure.

Please download the PLINK `.bed`, `.bim`, and `.fam` files for the two
datasets. After downloading, extract the files and rename them as
`IRRI_18M_3024` and `IRRI_404K_3024`.

<div class="rmdtip">
<p>The <code>%/%</code> operator does integer division
(<code>x %/% y</code> is equivalent to <code>floor(x/y)</code>) so the
index keeps track of which 80-line section of text we are counting up
negative and positive sentiment in.</p>
</div>

<div class="rmdnote">
<p>The reference genome for these datasets is ‘Nipponbare
MSU7/IRGSP1.0’.</p>
</div>

**Step 1: Pre-install Required Package**

``` r
install.packages("BiocManager") BiocManager::install(version = "3.19") BiocManager::install("qvalue")
```

``` r
##### R CODES #####   vv
```

## **Download Rice SNP data from SNP-Seek (3K RG)**
