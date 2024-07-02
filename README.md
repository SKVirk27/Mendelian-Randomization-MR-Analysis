# Mendelian-Randomization-MR-Analysis
Sure! Here's an updated version of the README file that includes references to other repositories for GWAS and GTEx data preparation, along with a note about the directory paths.

### README File Template

```markdown
# GWAS_to_MR_Pipeline

This repository contains scripts and data for formatting GWAS (Genome-Wide Association Study) summary statistics and performing Mendelian Randomization (MR) analysis using eQTL data.

## Overview

The goal of this project is to format GWAS summary statistics, calculate necessary metrics, lift over genomic coordinates to a consistent build (GRCh38), and perform MR analysis to identify significant associations between genetic variants and traits.

## Data Source

The input data for this project includes:
- GWAS summary statistics file (`HF_HRC_GWAS_UKBB_EUR_refGenome38.txt`)
- eQTL data file (`Open_filtered_df.txt`)

These files should be placed in the `data/` directory.

## Prerequisites

The following R packages are required to run the scripts in this repository:

- `dplyr`
- `stringr`
- `readr`
- `vroom`
- `tidyr`
- `tibble`
- `TwoSampleMR`
- `MVMR`
- `ieugwasr`
- `GwasDataImport`
- `ggplot2`

You can install these packages using the following command in R:

```r
install.packages(c("dplyr", "stringr", "readr", "vroom", "tidyr", "tibble", "TwoSampleMR", "MVMR", "ieugwasr", "GwasDataImport", "ggplot2"))
```

## Directory Structure

- `scripts/`: Contains the main analysis script.
- `data/`: Directory to store input data files.
- `output/`: Directory to store output files.

## Usage

1. **Download the Data**: Ensure your GWAS and eQTL data files are ready and place them in the `data/` directory.
2. **Run the Script**: Execute the R script located in the `scripts/` directory to perform the analysis.

## Step-by-Step Guide

1. **Load Necessary Libraries**: The script begins by loading the required R packages.
2. **Set Up Working Directory**: The working directory is set to the location of the script.
3. **Define Relative Paths**: The paths for input and output files are defined relative to the working directory.
4. **Create Output Directory**: The output directory is created if it doesn't exist.
5. **Load and Process eQTL Data**: The eQTL data is loaded and processed for MR analysis.
6. **Loop Through Genes for MR Analysis**: The script loops through each gene to perform MR analysis.
7. **Load GWAS Data**: The GWAS data is loaded and harmonized with the eQTL data.
8. **Perform MR Analysis**: MR analysis is conducted to identify significant associations.
9. **Save Results**: The results are saved to an RDS file.
10. **Combine and Visualize Significant Results**: Significant results from different tissues are combined and visualized using a bar plot.

## Example

```r
# Load necessary libraries
library(dplyr)
library(stringr)
library(readr)
library(vroom)
library(tidyr)
library(tibble)
library(TwoSampleMR)
library(MVMR)
library(ieugwasr)
library(GwasDataImport)
library(ggplot2)

# Set up the working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Define relative paths for input and output files
eqtl_file <- "data/Open_filtered_df.txt"
gwas_file <- "data/HF_HRC_GWAS_UKBB_EUR_refGenome38.txt"
mr_output_file <- "output/MR_results_Left_Ventricle_l.rds"
significant_output_file <- "output/significant_HR_MMP_l.rds"

# Create the output directory if it doesn't exist
output_dir <- "output"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Load eQTL data
Open_filtered_df <- read.table(eqtl_file, header = TRUE, sep = "\t")

# Create an empty list to store the results
results_list1 <- list()

# Get the unique genes
unique_genes <- unique(Open_filtered_df$hgnc_symbol)

# Loop through the genes and perform MR analysis
for (gene in unique_genes) {
  
  # Subset the data for the current gene
  eQTL_MR_data <- subset(Open_filtered_df, hgnc_symbol == gene, 
                         select = c(snp_col, beta_col, se_col, eaf_col, pval_col, effect_allele_col, other_allele_col))
  
  # Format the data for the current gene
  MR_eqtl <- format_data(eQTL_MR_data, type = "exposure",
                         snp_col = "snp_col",
                         beta_col = "beta_col",
                         se_col = "se_col",
                         effect_allele_col = "effect_allele_col",
                         other_allele_col = "other_allele_col",
                         eaf_col = 'eaf_col',
                         gene_col = 'hgnc_symbol',
                         pval_col = "pval_col") %>%
    mutate(exposure = paste(gene, "Left_Ventricle_metallopeptidase"))
  
  # Read in the outcome data
  outcome_dat <- read_outcome_data(
    snps = MR_eqtl$snp_col,
    filename = gwas_file,
    snp_col = "SNP",
    beta_col = "BETA",
    se_col = "SE",
    effect_allele_col = "A1",
    other_allele_col = "A2",
    eaf_col = "FRQ",
    pval_col = "P",
    units_col = "Unit",
    samplesize_col = "N"
  )
  
  # Harmonize the exposure and outcome data
  dat <- harmonise_data(exposure_dat = MR_eqtl, outcome_dat = outcome_dat)
  
  # Perform MR analysis
  res <- mr(dat)
  
  # Store the results in the list
  results_list1[[gene]] <- res
}

# Combine the results into one dataframe
results <- do.call(rbind, results_list1)

# Save the results to an RDS file
saveRDS(results, file = mr_output_file)

# Load MR results from different tissues
Artery_aorta_MR <- readRDS("output/MR_results_Artery_aorta_l.rds")
Artery_Coronary_MR <- readRDS("output/MR_results_Artery_Coronary_l.rds")
Artery_Tibial_MR <- readRDS("output/MR_results_Artery_tibial_l.rds")
Artery_Appendage_MR <- readRDS("output/MR_results_Artery_appendages_l.rds")
Left_ventricle_MR <- readRDS(mr_output_file)

# Combine significant results
significant_results <- do.call(rbind, lapply(list(Artery_aorta_MR, Artery_Coronary_MR, Artery_Tibial_MR, Artery_Appendage_MR, Left_ventricle_MR), function(x) x[x$pval < 0.05, ]))

# Save significant results
saveRDS(significant_results, significant_output_file)

# Create a data frame from the significant results data frame
df <- data.frame(exposure = significant_results$exposure,
                 b = significant_results$b,
                 se = significant_results$se,
                 pval = significant_results$pval)

# Create a bar plot
ggplot(significant_results, aes(x = exposure, y = b)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = b - se, ymax = b + se), width = 0.2) +
  scale_y_continuous(limits = c(min(significant_results$b - significant_results$se), max(significant_results$b + significant_results$se))) +
  theme_classic() +
  labs(x = "Exposure", y = "Effect size (b)") +
  theme(plot.title = element_text(hjust = 0.5))
```

## Additional Repositories

This project also relies on data and scripts from the following repositories:

- **GTEx eQTL Data Processing**: [GTEx_eQTL_Pipeline](https://github.com/SKVirk27/GTEx_eQTL_Pipeline) - This repository contains scripts for processing GTEx eQTL data and extracting significant variant-gene pairs.
- **GWAS Data Formatting**: [GWAS_Formatting_Pipeline](https://github.com/SKVirk27/MMP_GWAS_Formatting_38) - This repository includes scripts for formatting GWAS summary statistics and lifting over genomic coordinates to build 38.

Please refer to these repositories for detailed instructions on preparing the necessary input files for this analysis.

## Contact

If you have any questions or need further assistance, please contact:

- Simranjit Kaur Kang
- Simivk1991@gmail.com
```

