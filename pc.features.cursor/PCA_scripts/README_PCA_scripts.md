# PCA-Based Machine Learning Scripts for Methuselah Data

This folder contains modified versions of the original machine learning scripts that use Principal Component Analysis (PCA) instead of raw protein data.

**Note**: These scripts are now located in the `PCA_scripts/` folder, separate from the original scripts in the `scripts/` folder.

## Overview

The original scripts trained Elastic Net, SVM, and LightGBM models on the raw Methuselah dataset. These new scripts perform the same analysis but use PCA components as features instead of the original protein measurements.

## Files

### Main Scripts
- **`01.01.elastic.net.methuselah.iter.PCA.R`** - Elastic Net regression using PCA components
- **`01.02.svm.methuselah.iter.PCA.R`** - Support Vector Machine using PCA components  
- **`01.03.lightbgm.methuselah.iter.PCA.R`** - LightGBM using PCA components

### Analysis Scripts
- **`01.00.PCA.vs.original.comparison.R`** - Compares performance of PCA vs. original models

## Key Changes from Original Scripts

### 1. PCA Transformation
- Extracts protein data (excluding metadata columns)
- Performs PCA with scaling and centering
- Automatically determines number of PCs to retain (95% variance explained)
- Creates new dataset with PC scores instead of raw proteins

### 2. Data Structure
- Original: `SampleName`, `Sex`, `Age2`, `Protein1`, `Protein2`, ..., `ProteinN`
- PCA: `SampleName`, `Sex`, `Age2`, `PC1`, `PC2`, ..., `PCN`

### 3. Model Training
- Models now train on PC components instead of protein measurements
- Same cross-validation structure and hyperparameter tuning
- Same performance metrics (RMSE, correlation, AUC)

## Usage

### Step 1: Run Original Scripts
First, run the original scripts from the `scripts/` folder to establish baseline performance:
```bash
cd scripts
Rscript 01.01.elastic.net.methuselah.iter.R
Rscript 01.02.svm.methuselah.iter.R  
Rscript 01.03.lightbgm.methuselah.iter.R
cd ..
```

### Step 2: Run PCA-Based Scripts
Then run the PCA-based versions from the `PCA_scripts/` folder:
```bash
cd PCA_scripts
Rscript 01.01.elastic.net.methuselah.iter.PCA.R
Rscript 01.02.svm.methuselah.iter.PCA.R
Rscript 01.03.lightbgm.methuselah.iter.PCA.R
cd ..
```

### Step 3: Compare Results
Finally, run the comparison script from the `PCA_scripts/` folder:
```bash
cd PCA_scripts
Rscript 01.00.PCA.vs.original.comparison.R
cd ..
```

## Output Files

### PCA Results
- `PCA_results/PCA_results_[model]_[output.var].RData` - PCA objects and variance explained
- `PCA_results/PCA_loadings_[model]_top10_[output.var].csv` - Top 10 protein loadings for each PC

### Model Results
- `PCA_results/model.iter.12/[model].PCA/both/[output.var].model.[fold].csv` - Cross-validation results

### Comparison Results
- `PCA_results/PCA_vs_Original_Comparison.png` - Performance comparison plots
- `PCA_results/PCA_vs_Original_Summary.csv` - Summary statistics

## Benefits of PCA Approach

1. **Dimensionality Reduction**: Reduces from hundreds of proteins to ~10-50 PCs
2. **Noise Reduction**: Filters out measurement noise and technical variation
3. **Computational Efficiency**: Faster training and prediction
4. **Interpretability**: PCs represent major patterns of variation in the data
5. **Overfitting Prevention**: Fewer features reduce risk of overfitting

## PCA Parameters

- **Variance Threshold**: 95% (automatically determined)
- **Scaling**: Standard scaling (z-score) applied to proteins before PCA
- **Centering**: Mean centering applied before PCA
- **Missing Values**: Handled by imputation with column means

## Model-Specific Modifications

### Elastic Net
- Same alpha/lambda tuning grid
- Now trains on PC components instead of proteins
- Model type labeled as 'elastic.net.PCA'

### SVM
- RBF kernel with C and gamma tuning
- Scale parameter set to FALSE (data already scaled)
- Model type labeled as 'svm.PCA'

### LightGBM
- Same hyperparameter tuning grid
- Uses lgb.Dataset format for PC components
- Model type labeled as 'lightgbm.PCA'

## Troubleshooting

### Common Issues
1. **Missing packages**: Install required packages (e1071, lightgbm, pROC)
2. **Memory issues**: PCA reduces memory usage, but large datasets may still require sufficient RAM
3. **File paths**: Ensure data files are in the expected locations

### Performance Notes
- PCA computation time scales with number of proteins and samples
- Model training should be faster due to reduced dimensionality
- Cross-validation may take longer due to hyperparameter tuning

## Expected Results

- **RMSE**: May improve due to noise reduction, may worsen due to information loss
- **Correlation**: Should remain similar if PCA captures relevant variation
- **AUC**: Mortality prediction performance should be comparable
- **Interpretability**: PC loadings show which proteins contribute most to each component

## Next Steps

After running these scripts, consider:
1. Analyzing PC loadings to understand biological patterns
2. Testing different variance thresholds (90%, 99%)
3. Comparing with other dimensionality reduction methods (t-SNE, UMAP)
4. Investigating PC stability across cross-validation folds
