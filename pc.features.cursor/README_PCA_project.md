# Methuselah Data Analysis: Original vs. PCA-Based Models

This project contains machine learning scripts for analyzing Methuselah data using both original protein measurements and Principal Component Analysis (PCA) dimensionality reduction.

## Project Structure

```
├── scripts/                    # Original machine learning scripts
│   ├── 01.01.elastic.net.methuselah.iter.R
│   ├── 01.02.svm.methuselah.iter.R
│   ├── 01.03.lightbgm.methuselah.iter.R
│   └── ... (other original scripts)
│
├── PCA_scripts/               # PCA-based machine learning scripts
│   ├── 01.01.elastic.net.methuselah.iter.PCA.R
│   ├── 01.02.svm.methuselah.iter.PCA.R
│   ├── 01.03.lightbgm.methuselah.iter.PCA.R
│   ├── 01.00.PCA.vs.original.comparison.R
│   └── README_PCA_scripts.md
│
├── data/                      # Data files (methuselah.alpha, dataverse_files)
├── figures/                   # Results from original scripts
└── PCA_results/               # Results from PCA-based scripts (created when running)
```

## Quick Start

### 1. Run Original Models
```bash
cd scripts
Rscript 01.01.elastic.net.methuselah.iter.R
Rscript 01.02.svm.methuselah.iter.R
Rscript 01.03.lightbgm.methuselah.iter.R
cd ..
```

### 2. Run PCA-Based Models
```bash
cd PCA_scripts
Rscript 01.01.elastic.net.methuselah.iter.PCA.R
Rscript 01.02.svm.methuselah.iter.PCA.R
Rscript 01.03.lightbgm.methuselah.iter.PCA.R
cd ..
```

### 3. Compare Results
```bash
cd PCA_scripts
Rscript 01.00.PCA.vs.original.comparison.R
cd ..
```

## What This Project Does

- **Original Scripts**: Train Elastic Net, SVM, and LightGBM models on raw protein measurements
- **PCA Scripts**: Perform PCA on protein data, then train the same models on PC components
- **Comparison**: Automatically compare performance between original and PCA-based approaches

## Key Benefits of PCA Approach

1. **Dimensionality Reduction**: From hundreds of proteins to ~10-50 PCs
2. **Noise Reduction**: Filters out measurement noise
3. **Computational Efficiency**: Faster training and prediction
4. **Interpretability**: PCs represent major patterns of variation

## Output Locations

- **Original Results**: `figures/` folder
- **PCA Results**: `PCA_results/` folder
- **Comparison Results**: `PCA_results/` folder

## Requirements

- R with packages: tidyverse, caret, glmnet, e1071, lightgbm, pROC, openxlsx
- Data files in expected locations (see individual README files)

## Documentation

- **Original Scripts**: See individual script headers and `scripts/` folder
- **PCA Scripts**: See `PCA_scripts/README_PCA_scripts.md` for detailed documentation
- **Comparison**: See `PCA_scripts/01.00.PCA.vs.original.comparison.R` for analysis details
