# Extra Analysis and Visualization {#sec-extra}

## **Geographic Information**

Create a world map showcasing the percentage distribution of the origins of 327 rice samples, represented using a gradient color scale.

![](images/clipboard-1499951598.png){width="600"}

<details>

<summary>➡️ **Working on R (RStudio)**</summary>

``` r
##### R CODES #####

# Set the working directory for initial setup
setwd("...change...")

# Load the phenotype data from a CSV file
Pheno_data = read.csv("IRRI_311_Pheno_AG_Clean_ID_Sorted.csv")

# Extract the origin information and create a frequency table
origin = Pheno_data$Origin
table(origin)

# Create a data frame with region names, counts, and percentages
data = data.frame(
  "region" = names(table(origin)),  # Region names
  "number" = as.numeric(table(origin)),  # Number of occurrences per region
  "perc" = as.numeric(table(origin)) / length(Pheno_data$Origin)  # Percentage per region
)

# Load ggplot2 for visualization and map_data for world map data
library(ggplot2)
options(scipen = 999)  # Disable scientific notation for readability
world = map_data("world")  # Load world map data
head(world)  # Preview the world map data
write.csv(world, "world_map.csv")  # Save the world map data to a CSV file for reference

# Rename specific regions in the data frame for consistency with map data
data[4, 1] = "Brunei"
data[5, 1] = "Burkina Faso"
data[11, 1] = "Ivory Coast"
data[20, 1] = "Laos"
data[31, 1] = "Sierra Leone"
data[32, 1] = "Solomon Islands"
data[33, 1] = "Sri Lanka"
data[36, 1] = "Vietnam"
data[37, 1] = "Niger"

# Save the updated region information to a CSV file
write.csv(data, "origin_map.csv")

# Read the processed origin map data
origin_map = read.csv("origin_map.csv")

# Load additional libraries for data manipulation
library(dplyr)
library(stringr)
library(ggplot2)

# Merge the world map data with the origin map data by region
worldSubset = inner_join(world, origin_map, by = "region")
table(worldSubset$region)  # Check the frequency of regions in the subset

# Define a minimalist theme for the plot
plain = theme(
  axis.text = element_blank(),  # Remove axis text
  axis.line = element_blank(),  # Remove axis lines
  axis.ticks = element_blank(),  # Remove axis ticks
  panel.border = element_blank(),  # Remove panel border
  panel.grid = element_blank(),  # Remove grid lines
  axis.title = element_blank(),  # Remove axis titles
  panel.background = element_rect(fill = "white"),  # Set panel background to white
  plot.title = element_text(hjust = 0.5)  # Center the plot title
)

# Create a world map visualization with percentage-based fill
ggplot(data = worldSubset, mapping = aes(x = long, y = lat, group = group)) +
  # Add the base world map with light gray fill
  geom_polygon(data = world, aes(x = long, y = lat, group = group),
               fill = "#f7f6f6", color = "#999999", linewidth = 0.01) +
  # Fill regions based on percentage
  geom_polygon(aes(fill = perc)) +
  # Define a gradient color scale for the fill
  scale_fill_gradient(
    low = "#f2d489",  # Low percentage color
    high = "#e8542d",  # High percentage color
    name = "Percentage (%)",  # Legend title
    breaks = c(4, 8, 12),  # Breaks for the legend
    labels = c("4", "8", "12"),  # Labels for the legend
    guide = guide_colorbar(
      barwidth = 10,  # Width of the color bar
      barheight = 1.2  # Height of the color bar
    )
  ) +
  plain +  # Apply the minimalist theme
  # Customize the legend appearance and position
  theme(
    legend.title = element_text(size = 14),  # Legend title font size
    legend.text = element_text(size = 12),  # Legend text font size
    legend.position = "bottom"  # Position the legend at the bottom
  )

# Final plot size recommendation: 7' x 4'
```

</details>

------------------------------------------------------------------------

## Linkage Disequilibrium (LD) Decay Across Genome

Computation and visualization of LD decay. First, utilizes PLINK to thin SNPs (10%) and calculate pairwise LD (r²) within specified distance windows (1000 Kb). Then, Analyzes PLINK output to compute average r² across distance intervals (`window.size`), summarizing LD decay patterns. Finally, generates a smooth LD decay plot, highlighting the relationship between physical distance (kb) and LD (r²).

![](images/clipboard-2982375851.png){width="600"}

<details>

<summary>➡️ **Working on PLINK in the terminal**</summary>

``` r
##### PLINK COMMANDS #####

# Change your path into the directory of the 'IRRI_18M_311'
# cd ...your file path...

./plink --bfile IRRI_1.9M_327 --thin 0.1 --set-hh-missing --r2 --ld-window-r2 0 --ld-window 99999 --ld-window-kb 1000 --out IRRI_190K_327_LD
```

</details>

<details>

<summary>➡️ **Working on R (RStudio)**</summary>

``` r
##### R CODES #####

# Set the working directory for initial setup
setwd("...change...")

# Load the required library
library(data.table)

# Read the input data file containing LD information
LD = fread("IRRI_190K_311_LD.ld")

# Function to calculate LD decay
LD.decay = function(data, window.size){ 
  # Order the data by distance
  ld1 = data[order(data[, 5] - data[, 2]), ]
  distance = ld1[, 5] - ld1[, 2]  # Calculate distance between loci
  ld_distance = data.frame(distance, R2 = ld1$R2)  # Create a data frame with distance and R² values
  
  # Define blocks for the distance intervals
  block = seq(1, max(ld_distance[, 1]), by = window.size)
  ld2 = rep(0, length(block))  # Initialize an array to store mean R² values
  decay = data.frame(block, ld2)  # Create a data frame to store LD decay information
  
  # Initialize a progress bar
  pb = txtProgressBar(min = 0, max = length(block), style = 3)
  
  # Loop through each block to calculate mean R²
  for (i in seq_along(block)) {
    # Subset data within the current block range
    intv = subset(ld_distance, (ld_distance[, 1] >= block[i] & ld_distance[, 1] < (block[i] + window.size)))
    # Compute the mean R² for the current block
    decay[i, 2] = mean(intv[, 2], na.rm = TRUE)
    # Update the progress bar
    setTxtProgressBar(pb, i)
  }
  
  # Close the progress bar
  close(pb)
  
  return(decay)  # Return the LD decay data frame
}

# Run the LD decay function with a window size of 1000
LD.decay = LD.decay(data = LD, window.size = 1000)

# Save the LD decay data as an RDS file
saveRDS(LD.decay, "IRRI_190K_327_LD.rds")

# Load the saved LD decay data
LD.decay = readRDS("IRRI_190K_327_LD.rds")

# Load ggplot2 for visualization
library(ggplot2)

# Plot LD decay
ggplot(LD.decay) +
  labs(x = "Distance (kb)", y = expression(r^{2})) +  # Axis labels with proper formatting
  geom_line(aes(x = block, y = ld2), color = "#d3a02f", size = 2) +  # Add a line plot for LD decay
  scale_x_continuous(
    limits = c(0.0, 4.5*10^5),  # Set x-axis limits
    breaks = c(0, 1*10^5, 2*10^5, 3*10^5, 4*10^5),  # Define breaks
    labels = c("0", "100", "200", "300", "400"),  # Label breaks in kb
    expand = c(0.01, 0)  # Adjust axis expansion
  ) +
  scale_y_continuous(
    limits = c(0.0, 0.65),  # Set y-axis limits
    breaks = c(0, 0.2, 0.4, 0.6),  # Define breaks
    labels = c("0.0", "0.2", "0.4", "0.6"),  # Label breaks
    expand = c(0.03, 0)  # Adjust axis expansion
  ) +
  theme_classic() +  # Apply a classic theme
  theme(
    axis.title.x = element_text(size = 18),  # Customize x-axis title font size
    axis.title.y = element_text(size = 18, vjust = 2),  # Customize y-axis title font size and position
    axis.text.x = element_text(size = 16),  # Customize x-axis text font size
    axis.text.y = element_text(size = 16)   # Customize y-axis text font size
  )
  
# Final plot size recommendation: 6' x 4'
```

</details>

------------------------------------------------------------------------

## Genetic Pattern

![](images/clipboard-1482032998.png){width="600"}

![](images/clipboard-2271077002.png){width="600"}

<details>

<summary>➡️ **Working on PLINK in the terminal**</summary>

``` r
##### PLINK COMMANDS ##### 

# Change your path into the directory of the 'IRRI_18M_327'
# cd ...your file path...

system("plink --bfile IRRI_404K_3024 --pca 3024")
```

</details>

<details>

<summary>➡️ **Working on R (RStudio)**</summary>

``` r
##### R CODES #####

# Set the working directory for initial setup
setwd("...change...")

eigenval = read.table("plink.eigenval")
eigenval = as.numeric(eigenval[,1])
eigenval = eigenval[which(eigenval>0)]

##### Eigenvalue of PCA from plink #####
variance_proportion = eigenval / sum(eigenval) # PC1: 8.11%; PC2: 4.17%; PC3: 2.52%
variance_proportion[1:50]*100

library(ggplot2)
library(ggthemes)

plot_pca_variance = function(variance_proportion, num_pcs = 15, num_bins = 4) {
  
  if (!is.numeric(variance_proportion)) {
    stop("variance_proportion must be a numeric vector")
  }
  
  variance_data = data.frame(
    PCs = 1:num_pcs,
    PV = variance_proportion[1:num_pcs] * 100,  # Proportion of Variance
    CV = cumsum(variance_proportion[1:num_pcs] * 100)  # Cumulative Proportion
  )
  col = c(rep("#765635", num_bins), rep("#e5d5c5", num_pcs - num_bins))
  
  ylim.PV = c(0, 10)
  ylim.CP = c(0, 100)
  b = diff(ylim.PV) / diff(ylim.CP)
  a = ylim.PV[1] - b * ylim.CP[1]
  
  plot = ggplot(variance_data, aes(x = PCs, y = PV)) +
    geom_bar(stat = "identity", show.legend = FALSE, fill = col) +
    geom_line(aes(y = a + CV * b), color = "grey30", lwd = 1.5) +
    geom_point(aes(y = a + CV * b), color = "grey10", size = 4) +
    xlab("Principal components axis") +
    scale_x_continuous(breaks = c(1, 5, 10, 15, num_pcs)) +
    scale_y_continuous("Proportion of variance (%)", limits = c(0, 10), 
                       sec.axis = sec_axis(~(. - a) / b, name = "Cumulative proportion (%)")) +
    theme_base() +
    theme(
      axis.title.x = element_text(size = 16),
      axis.title.y.left = element_text(size = 16, color = "#765635"),
      axis.title.y.right = element_text(size = 16, color = "grey10"),
      axis.text.x = element_text(size = 14),
      axis.text.y.left = element_text(size = 13),
      axis.text.y.right = element_text(size = 13, color = "grey10")
    ) +
    theme(legend.position = "none")
  return(plot) # SIZE: A6
}
plot_pca_variance(variance_proportion, num_pcs = 20, num_bins = 3)

##### Eigenvector of PCA from plink #####
library(data.table)
eigenvec = fread("plink.eigenvec")
eigenvec = as.data.frame(eigenvec[, 2:6])
colnames(eigenvec) = c("ID", "PC1", "PC2", "PC3", "PC4")

setwd("C:/Users/teddy/Desktop/TAMU/P1/404K/327")
ID = read.csv("327in3024.csv")
data = cbind(ID[,c(2,4,5,6)], eigenvec)

# PCA Scatter Plot: PC1 vs PC2
library(ggplot2)
PC1v2 = ggplot(data, aes(x = PC1, y = PC2, color = GWAS)) +
  geom_point(size = 3, alpha = 0.5) +
  labs(x = "PC 1 (8.11%)", y = "PC 2 (4.17%)") +
  theme_classic() +
  theme(legend.position = "none") +
  theme(
    axis.title.x =   element_text(size = 16),
    axis.title.y =   element_text(size = 16),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10)) +
  scale_color_manual(values = c("grey","#bd9233"))
PC1v2 +
  geom_point(data = subset(data, GWAS == "TRUE"), aes(PC1, PC2), color = "#bd9233", size = 3, alpha = 0.5) +
  geom_point(data = subset(data, GWAS == "TRUE"), aes(PC1, PC2), color = "#bd9233", size = 3, shape = 20)
# 5*7

# PCA Scatter Plot: PC1 vs PC3
PC1v3 = ggplot(data, aes(x = PC1, y = PC3, color = GWAS)) +
  geom_point(size = 3, alpha = 0.5) +
  labs(x = "PC 1 (8.11%)", y = "PC 3 (2.52%)") +
  theme_classic() +
  theme(legend.position = "none") +
  theme(
    axis.title.x =   element_text(size = 16),
    axis.title.y =   element_text(size = 16),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10)) +
  scale_color_manual(values = c("grey","#bd9233"))
PC1v3 +
  geom_point(data = subset(data, GWAS == "TRUE"), aes(PC1, PC3), color = "#bd9233", size = 3, alpha = 0.5) +
  geom_point(data = subset(data, GWAS == "TRUE"), aes(PC1, PC3), color = "#bd9233", size = 3, shape = 20)

# PCA Scatter Plot: PC1 vs PC2 vs PC3
library(plotly)
cols = colnames(data)
text_content = paste0('"</br> ◉ ', cols, ': ", data$', cols, collapse = ', ')
text = as.formula(paste("~paste(", text_content, ")"))

plot_ly(data, x = ~PC1, y = ~PC2, z = ~PC3,
        type = "scatter3d", mode = "markers", 
        color = ~GWAS, colors = c("grey","#bd9233"),
        marker = list(size = c(2, 4)),
        text = text) %>%
  layout(title = NA,
         scene = list(
           xaxis = list(title = "PC1 (8.11%)", showline = TRUE, showgrid = TRUE, zeroline = TRUE, zerolinecolor = "grey"),
           yaxis = list(title = "PC2 (4.17%)", showline = TRUE, showgrid = TRUE, zeroline = TRUE, zerolinecolor = "grey"),
           zaxis = list(title = "PC3 (2.52%)", showline = TRUE, showgrid = TRUE, zeroline = TRUE, zerolinecolor = "grey")))
```

</details>

------------------------------------------------------------------------

## **Backup**

#### **Align Accession IDs**

Match the accession IDs between the phenotypic data and SNP datasets. We identified three issues related to inconsistencies in the spelling, formatting, and ranking of ID names between the phenotype data and genotype datasets:

1.  The ID names in the phenotype data must correspond precisely with those in the genotype datasets, using the Mapping Table as a reference.

2.  The format of ID names in the Mapping Table is inconsistent with the format used in the phenotype data. For example, consider a sample named `EX WUKARI (WILD)::IRGC 63214-3`:

    | ID in Phenotype | ID in [Mapping Table](https://3kricegenome.s3.us-east-1.amazonaws.com/kaust_irri_3k_16refs/3K_list_sra_ids.txt) | ID in Genotype |
    |------------------------|------------------------|------------------------|
    | EX WUKARI (WILD)::IRGC 63214-3 | EX_WUKARI\_(WILD)::IRGC_63214-1 | IRIS_313-11657 |
    |  | IRIS_313-11657 |  |

3.  In the genotype datasets, the rank order of ID names does not follow the correct sequence. For example, the correct order should be: `IRIS_313-8890`, `IRIS_313-9101`, and `IRIS_313-11657`, not `IRIS_313-11657`, `IRIS_313-8890`, and `IRIS_313-9101`.

<details>

<summary>➡️ **Working on R (RStudio)**</summary>

``` r
##### R CODES #####  # Set the working directory setwd("...your file path...")  # Load the datasets IRRI_3K = read.csv("IRRI_3K_ID.csv") # Mapping table  IRRI_341 = read.csv("IRRI_341_Pheno_AG_Original.csv") # Phenotypic Data  # Process designations in IRRI_3K dataset IRRI_3K[1:5,2] IRRI_3K_ID = gsub("_", " ", IRRI_3K[,2]) IRRI_341[1:5,1]  # Match designations between the two datasets loc = c() for (i in 1:dim(IRRI_341)[1]) {   match_index = which(IRRI_3K_ID == IRRI_341$Designation[i])   if (length(match_index) == 0) {     loc[i] = NA    } else {     loc[i] = match_index   } }  # Check for unmatched entries and summarize loc sum(is.na(loc))   # Create a new annotated data frame combining information from both datasets IRRI_341_ID = data.frame("Designation" = IRRI_341$Designation) for (i in 1:length(loc)) {   if (is.na(loc[i])) {     IRRI_341_ID$ID = NA   } else {     IRRI_341_ID$ID = IRRI_3K$X3K_DNA_IRIS_UNIQUE_ID[loc]     IRRI_341_ID$Origin = IRRI_3K$Country_Origin_updated[loc]     IRRI_341_ID$SRA = IRRI_3K$SRA[loc]     IRRI_341_ID$Control_5DAS = IRRI_341$PERCENT_SEEDGER_5DAS     IRRI_341_ID$Control_7DAS = IRRI_341$PERCENT_SEEDGER_7DAS     IRRI_341_ID$Control_14DAS = IRRI_341$PERCENT_SEEDGER_14DAS     IRRI_341_ID$AG_7DAS = IRRI_341$PERCENT_SEEDGER_7DAS_AG     IRRI_341_ID$AG_14DAS = IRRI_341$PERCENT_SEEDGER_14DAS_AG   } }  # Preview the annotated data frame head(IRRI_341_ID) tail(IRRI_341_ID)  # Save the annotated data frame to a CSV file write.csv(IRRI_341_ID, "IRRI_341_Pheno_AG_Original_ID.csv", row.names = F)  # Sort the annotated data by ID and identify specific entries IRRI_341_ID_sorted = IRRI_341_ID[order(IRRI_341_ID$ID), ] which(IRRI_341_ID_sorted$ID == "IRIS_313-7620") # 245 which(IRRI_341_ID_sorted$ID == "IRIS_313-9989") # 311 which(IRRI_341_ID_sorted$ID == "IRIS_313-10020") # 5 which(IRRI_341_ID_sorted$ID == "IRIS_313-15901") # 5  # Rearrange the data frame to exclude specific entries IRRI_311_ID_sorted = IRRI_341_ID_sorted[c(1:4, 245:311, 5:244), ]  # Save the sorted and rearranged data frame to a CSV file write.csv(IRRI_311_ID_sorted, "IRRI_311_Pheno_AG_Original_ID_Sorted.csv", row.names = F)  # Export a text file with the IDs for further analysis write.table(cbind(IRRI_311_ID_sorted$ID, IRRI_311_ID_sorted$ID), "IRRI_311_ID.txt", row.names = F, col.names = F, quote = F)
```

</details>

<details>

<summary>**Outputs**</summary>

-   **IRRI_341_Pheno_AG_Original_ID.csv**\
    Annotated phenotypic data with aligned and matched ID names.

-   **IRRI_311_Pheno_AG_Original_ID_Sorted.csv**\
    Phenotypic data that has been sorted and reorganized.

-   **IRRI_311_ID.txt**\
    A text file containing the list of IDs.

![](images/output.png){width="400"}

</details>
