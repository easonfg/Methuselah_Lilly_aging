rm(list = ls())

library(tidyverse)
library(TAPGeneralization)
library(doParallel); library(foreach)
library(caret)
library(openxlsx)
library(pROC)

seed.i = 12
set.seed(seed.i)
# output.var = 'inflam'
# sex.i = 'M'
sex.i = 'both'

# for (output.var in c('inflam', 'cns', 'both')){
# for (output.var in c('both', 'cns', 'inflam')){

lightgbm.roc = list()
lightgbm.stats.res = list()
svm.stats.res = list()
svm.roc = list()
for (output.var in c('both')){
  # for (sex.i in c('F')){
  # for (sex.i in c('M', 'F')){
  if (output.var == 'inflam'){
    data <- read_csv("data/methuselah.alpha/inflammation/20240520_P-000412_Methuselah-Foundation_NULISAseq_TAP_Counts_Report_FULL.csv") %>%
      filter(!grepl("SC", SampleType))

  } else if (output.var == 'cns'){
    data <- read_csv("data/methuselah.alpha/CNS/20240520_P-000412_Methuselah-Foundation_NULISAseq_TAP_Counts_Report_FULL.csv") %>%
      filter(!grepl("SC", SampleType))
  } else if (output.var == 'both'){
    data.inflam <- read_csv("data/methuselah.alpha/inflammation/20240520_P-000412_Methuselah-Foundation_NULISAseq_TAP_Counts_Report_FULL.csv") %>%
      filter(!grepl("SC", SampleType))
    data.cns <- read_csv("data/methuselah.alpha/CNS/20240520_P-000412_Methuselah-Foundation_NULISAseq_TAP_Counts_Report_FULL.csv") %>%
      filter(!grepl("SC", SampleType))
    data = rbind(data.inflam, data.cns) %>% group_by(SampleName, Target, SampleType) %>%
      summarise(NPQ = mean(NPQ)) %>% ungroup()
  }

  metadata = read_csv('data/methuselah.alpha/sample_metadata.csv')
  metadata %>% group_by(Sex, Race) %>% summarise(n_distinct(SampleName))
  metadata

  full.data = data %>% left_join(metadata, by = c('SampleName', 'SampleType'))
  full.data.m = full.data %>% filter(Sex == 'M')
  full.data.f = full.data %>% filter(Sex == 'F')
  full.data
  # full.data = full.data %>% filter(Sex == sex.i)

  full.data.wide = full.data %>% select(SampleName, Sex, Age2, Target, NPQ) %>%
    pivot_wider(names_from = 'Target', values_from = 'NPQ')
  full.data.wide.m = full.data.m %>% select(SampleName, Sex, Age2, Target, NPQ) %>%
    pivot_wider(names_from = 'Target', values_from = 'NPQ')
  full.data.wide.f = full.data.f %>% select(SampleName, Sex, Age2, Target, NPQ) %>%
    pivot_wider(names_from = 'Target', values_from = 'NPQ')

  full.data.wide

  ################################################################################
  ### PCA TRANSFORMATION ###
  
  # Extract protein data for PCA (excluding metadata columns)
  protein_data <- full.data.wide %>% 
    select(-c(SampleName, Sex, Age2)) %>%
    as.matrix()
  
  # Check for missing values and handle them
  if(any(is.na(protein_data))) {
    # Replace NAs with column means
    protein_data <- apply(protein_data, 2, function(x) {
      x[is.na(x)] <- mean(x, na.rm = TRUE)
      return(x)
    })
  }
  
  # Perform PCA
  pca_result <- prcomp(protein_data, scale. = TRUE, center = TRUE)
  
  # Determine number of PCs to retain (explain 95% of variance)
  variance_explained <- cumsum(pca_result$sdev^2 / sum(pca_result$sdev^2))
  n_pcs <- which(variance_explained >= 0.95)[1]
  if(is.na(n_pcs)) n_pcs <- min(50, ncol(protein_data)) # fallback
  
  cat("Number of PCs retained:", n_pcs, "\n")
  cat("Variance explained by", n_pcs, "PCs:", round(variance_explained[n_pcs] * 100, 2), "%\n")
  
  # Extract PC scores
  pc_scores <- pca_result$x[, 1:n_pcs]
  colnames(pc_scores) <- paste0("PC", 1:n_pcs)
  
  # Create new dataset with PCs instead of raw proteins
  full.data.wide.pca <- data.frame(
    SampleName = full.data.wide$SampleName,
    Sex = full.data.wide$Sex,
    Age2 = full.data.wide$Age2,
    pc_scores
  )
  
  # Save PCA results for later analysis
  dir.create('PCA_results', recursive = TRUE, showWarnings = FALSE)
  save(pca_result, n_pcs, variance_explained, file = paste0('PCA_results/PCA_results_', output.var, '.RData'))
  
  # Use PCA data for modeling
  full.data.wide <- full.data.wide.pca
  
  ################################################################################
  ### scaling data (now scaling PCs instead of proteins)
  full.data.wide <- full.data.wide %>%
    mutate(across(
      # Select columns to scale: everything EXCEPT Name, Sex, Age
      .cols = -c(SampleName, Sex, Age2),
      # Apply the scale function to each selected column
      .fns = scale
    ))

  # View the result
  head(full.data.wide)

  ################################################################################
  ### classification

  full.data.wide$Sex = ifelse(full.data.wide$Sex == 'M', 1, 0)
  full.data.wide
  adjusted.wide = full.data.wide
  x.variables = adjusted.wide %>% select(-c(SampleName, Age2, ))
  logical.complete.cases = complete.cases(x.variables)
  x.variables = x.variables[logical.complete.cases,]
  y.response = adjusted.wide$Age2
  y.response = y.response[logical.complete.cases]
  meta.data = adjusted.wide %>% select(SampleName, Age2, Sex)

  adjusted.wide
  md3cv <- glmnet::cv.glmnet(as.matrix(x.variables),
                             y.response,
                             nfold = 10, type.measure = "mse", parallel = TRUE, alpha = 1)

  md3cv$cvm[md3cv$lambda == md3cv$lambda.1se]
  md3cv$cvm[md3cv$lambda == md3cv$lambda.min]

  glmnet_pred = predict(md3cv, as.matrix(x.variables))
  glmnet_pred
  mydata = cbind(glmnet_pred, y.response)
  plot(glmnet_pred ~ y.response)
  abline(coef = c(0,1))

  sqrt(mean((glmnet_pred - y.response)^2))

  return_features = function(coeff){
    top_features = data.frame(name = coeff@Dimnames[[1]][coeff@i + 1], coefficient = coeff@x)
    top_features = top_features[order(top_features$coefficient),]
    # print(top_features)
    return(top_features)
  }

  features = return_features(coef(md3cv, s = 'lambda.1se'))
  features
  features %>% dim()
  # write.csv(features %>% arrange(desc(abs(coefficient))), paste0('figures/classifications/', output.var, '/one.model.csv'))
  features %>% arrange(desc(abs(coefficient)))
  features %>% head()

  ### cross validation
  features.list = list()
  train.rmse.list = list()
  test.rmse.list = list()
  alpha.list = list()
  roc.train.list = list()
  roc.test.list = list()
  model.list = list()
  counter = 1

  test.fold <- caret::createFolds(y.response, k = 5, list = TRUE, returnTrain = FALSE)
  for (test.fold.i in 1:5){
    test.set.x = x.variables[test.fold[[test.fold.i]],]
    test.set.y = y.response[test.fold[[test.fold.i]]]
    remains.x.variables = x.variables[-test.fold[[test.fold.i]],]
    remains.y.response = y.response[-test.fold[[test.fold.i]]]
    test.set.meta = meta.data[test.fold[[test.fold.i]],]

    model.type = 'elastic.net.PCA'
    dir.create(paste0('PCA_results/model.iter.', seed.i, '/',  model.type, '/', sex.i, '/'), recursive = T)

    test = test.set.x %>% as.matrix()
    y_test = test.set.y

    flds <- caret::createFolds(remains.y.response, k = 5, list = TRUE, returnTrain = FALSE)
    fold.i = flds[[1]]
    vali.set.x = remains.x.variables[fold.i,]
    vali.set.y = remains.y.response[fold.i]
    train.set.x = remains.x.variables[-fold.i,]
    train.set.y = remains.y.response[-fold.i]

    train <- train.set.x %>% as.matrix()
    train
    vali <- vali.set.x %>% as.matrix()

    y_train <- train.set.y
    y_vali <- vali.set.y

    ### iter light gbm
    mortality.data = read.xlsx('data/dataverse_files/BoAC_plasma_metadata_TTD.xlsx') %>% as_tibble()
    mortality.data
    data.score.mortality = test.set.meta %>%     mutate(Plasma.ID = str_extract(SampleName, "[^_]+$")) %>%
      left_join(mortality.data %>% select(Plasma.ID, Vital), by = 'Plasma.ID')

    ################################################################################
    # Define the tuning grid
    tune_grid <- expand.grid(
      alpha = seq(0, 1, by = 0.1),  # Test alpha from 0 (Ridge) to 1 (Lasso)
      lambda = 10^seq(-3, 2, length = 100)  # Test lambda from 0.001 to 10
    )
    tune_grid

    manual.iter = mclapply(1:nrow(tune_grid),
                           mc.cores = parallel::detectCores() - 1,
                           function(j){
                             glm.manual = glmnet(
                               x = remains.x.variables,
                               y = remains.y.response,
                               alpha = tune_grid$alpha[j],
                               lambda = tune_grid$lambda[j]
                             )
                             new_pred <- predict(glm.manual, test.set.x %>% as.matrix())
                             new_pred
                             pred.train <- predict(glm.manual, remains.x.variables %>% as.matrix())
                             test.rmse = RMSE(new_pred, test.set.y)
                             test.rmse
                             error = RMSE(pred.train, remains.y.response)
                             test.cor = cor(new_pred, test.set.y)

                             data.score.mortality$score = new_pred
                             model <- lm(score ~ Age2 + Sex, data = data.score.mortality)
                             adjusted_scores <- residuals(model)
                             data.score.mortality$adjusted_score = adjusted_scores
                             score.roc_obj <- roc(data.score.mortality$Vital, data.score.mortality$adjusted_score)

                             auc.j = auc(score.roc_obj)
                             auc.j

                             data.frame(a = j, error = error, rmse = test.rmse, cor = test.cor, auc = auc.j)
                           })
    manual.iter.list = manual.iter %>% data.table::rbindlist()
    manual.iter.list = cbind(manual.iter.list, tune_grid)
    manual.iter.list %>% filter(rmse == min(rmse))
    manual.iter.list %>% arrange(auc) %>%
      mutate(a = fct_inorder(as.character(a))) %>%
      ggplot() +
      geom_point(aes(a, rmse), color = 'green') +
      geom_point(aes(a, error))

    write.csv(manual.iter.list, paste0('PCA_results/model.iter.', seed.i, '/',  model.type, '/', sex.i, '/', output.var,
                                 '.model.', test.fold.i, '.csv'))
    ################################################################################
    ################################################################################

  }

  # }
}

# Save PCA loadings for interpretation
pca_loadings <- data.frame(
  Protein = rownames(pca_result$rotation),
  pca_result$rotation[, 1:min(10, n_pcs)]
)
write.csv(pca_loadings, paste0('PCA_results/PCA_loadings_top10_', output.var, '.csv'), row.names = FALSE)

cat("PCA-based Elastic Net analysis completed!\n")
cat("Number of PCs used:", n_pcs, "\n")
cat("Variance explained:", round(variance_explained[n_pcs] * 100, 2), "%\n")
