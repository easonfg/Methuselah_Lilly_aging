# Protein-Age Regression Analysis

This folder contains scripts and results for analyzing the relationship between protein expression (NPQ) and age in the Methuselah dataset.

## Overview

The analysis performs linear regression of age on each protein's NPQ expression, adjusts p-values for multiple testing using FDR (False Discovery Rate), and creates visualizations including a volcano plot.

## Script: `01.protein_age_regression_analysis.R`

### What it does:

1. **Data Loading**: Loads Methuselah inflammation and CNS data, along with metadata
2. **Data Processing**: Combines data and creates wide format for analysis
3. **Regression Analysis**: For each protein, performs `age ~ protein + sex` regression
4. **Fold Change Calculation**: Calculates log2 fold change between young (<65) and old (≥65) age groups
5. **Multiple Testing Correction**: Applies FDR (Benjamini-Hochberg) correction to p-values
6. **Visualization**: Creates volcano plot and distribution plots
7. **Results Export**: Saves all results to CSV files

### Model Details:

- **Regression Model**: `age ~ protein + sex`
- **Age Effect**: Coefficient for protein term (how much age changes per unit change in protein)
- **Fold Change**: Old vs. Young median expression ratio (log2 scale)
- **Significance**: FDR-adjusted p-values with thresholds at 0.05, 0.01, 0.001

### Output Files:

1. **`protein_age_regression_results.csv`** - Complete results for all proteins
2. **`volcano_plot.pdf/png`** - Main volcano plot showing age effect vs. fold change
3. **`distribution_plots.pdf/png`** - P-value, adjusted p-value, and fold change distributions
4. **`summary_statistics.csv`** - Summary statistics across all proteins
5. **`top_significant_proteins.csv`** - Top 20 most significant proteins

## Running the Analysis

```bash
cd protein_age_regression
Rscript 01.protein_age_regression_analysis.R
```

## Results Interpretation

### Volcano Plot:
- **X-axis**: Log2 fold change (Old/Young)
- **Y-axis**: -Log10(adjusted p-value)
- **Point size**: Absolute age effect magnitude
- **Colors**: Significance categories
- **Reference lines**: 
  - Red dashed: p < 0.05 threshold
  - Blue dashed: |log2FC| > 0.5 threshold

### Significance Categories:
- **Significant & High FC**: Both p < 0.05 and |log2FC| > 0.5
- **Significant**: Only p < 0.05
- **High FC**: Only |log2FC| > 0.5
- **Not Significant**: Neither threshold met

### Key Metrics:
- **Age Effect**: Regression coefficient (positive = older age associated with higher protein)
- **Fold Change**: Expression ratio between age groups
- **R-squared**: Model fit quality
- **Adjusted P-value**: FDR-corrected significance

## Data Requirements

The script expects:
- Methuselah inflammation data: `data/methuselah.alpha/inflammation/...`
- Methuselah CNS data: `data/methuselah.alpha/CNS/...`
- Metadata: `data/methuselah.alpha/sample_metadata.csv`

## Dependencies

Required R packages:
- tidyverse
- ggplot2
- gridExtra
- patchwork

## Notes

- Age threshold for fold change calculation is set to 65 years
- Minimum sample size requirement is 10 for regression analysis
- Sex is included as a covariate in all regression models
- Missing values are handled by complete case analysis
