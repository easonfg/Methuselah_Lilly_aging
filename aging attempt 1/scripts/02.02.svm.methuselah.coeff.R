# svm1 = read_csv('figures/model.iter.12/SVM/both/both.model.1.csv')
# svm5 = read_csv('figures/model.iter.12/SVM/both/both.model.5.csv')
# svm1
# svm1 %>% filter(rmse == min(rmse))
# svm5 %>% filter(rmse == min(rmse))
# svm1 %>% filter(auc == max(auc))
# svm5 %>% filter(auc == max(auc))
#
# cbind(svm1, svm1 %>% scale() %>% data.frame() %>% mutate(sum = -rmse + (auc))) %>% data.frame() %>%
# # cbind(svm1, svm1 %>% scale() %>% data.frame() %>% mutate(sum = rmse + (-auc))) %>%
#   filter(sum == max(sum))
#
#
# svm1 %>% arrange(rmse)
# svm5 %>% arrange(rmse)
#
#

library(tidyverse)
library(TAPGeneralization)
library(doParallel); library(foreach)
library(caret)
library(openxlsx)

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

  # lm_res <- lmNULISAseq(
  #   data = as.matrix(full.data.wide %>% select(-c(Sex, Age2)) %>% column_to_rownames('SampleName') %>% t()),
  #   sampleInfo = metadata,
  #   # sampleName_var = "SampleID",
  #   # modelFormula = "Group",
  #   sampleName_var = 'SampleName',
  #   modelFormula = 'Age2 + Sex',
  #   reduced_modelFormula = 'Sex',
  #   # modelFormula = 'Age2 + Sex + Age2*Sex',
  #   # reduced_modelFormula = 'Age2 + Sex',
  #   return_model_fits = T
  # )
  # lm_res$Fstats %>% filter(Ftest_pval_FDR<0.05) %>% dim()
  # lm_res$modelStats %>% head()
  # lm_res$modelStats %>% filter(Age2_pval_FDR < 0.05) %>% dim()
  # # lm_res$modelStats %>% filter(Age2.SexM_pval_FDR < 0.05)
  # lm_res$Fstats %>% arrange(Ftest_pval_FDR) %>% head()
  # lm_res$modelStats %>% arrange(Age2_pval_FDR) %>% head()


  ################################################################################
  ### classification



  # scaled.wide = full.data.wide  %>% column_to_rownames('SampleName') %>%
  #   mutate(Sex = ifelse(Sex == 'M', 1, 0)) %>%
  #   scale()
  # adjusted.wide = scaled.wide %>% data.frame()
  # adjusted.wide
  adjusted.wide = full.data.wide
  x.variables = adjusted.wide %>% select(-c(SampleName, Age2, Sex))
  logical.complete.cases = complete.cases(x.variables)
  x.variables = x.variables[logical.complete.cases,]
  y.response = adjusted.wide$Age2
  y.response = y.response[logical.complete.cases]
  meta.data = adjusted.wide %>% select(SampleName, Age2, Sex)


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

    # female
    # test.fold.f <- caret::createFolds(y.response.f, k = 6, list = TRUE, returnTrain = FALSE)
    # test.set.x.f = x.variables.f[test.fold.f[[1]],]
    # test.set.y.f = y.response.f[test.fold.f[[1]]]
    # remains.x.variables.f = x.variables.f[-test.fold.f[[1]],]
    # remains.y.response.f = y.response.f[-test.fold.f[[1]]]
    #
    # # male
    # test.fold.m <- caret::createFolds(y.response.m, k = 6, list = TRUE, returnTrain = FALSE)
    # test.set.x.m = x.variables.m[test.fold.m[[1]],]
    # test.set.y.m = y.response.m[test.fold.m[[1]]]
    # remains.x.variables.m = x.variables.m[-test.fold.m[[1]],]
    # remains.y.response.m = y.response.m[-test.fold.m[[1]]]
    #
    ###############################################################################
    ###############################################################################
    ### lightGBM


    test = test.set.x %>% as.matrix()
    y_test = test.set.y
    #
    # flds <- caret::createFolds(remains.y.response, k = 5, list = TRUE, returnTrain = FALSE)
    # flds

    # for (fold.i in flds){
    #   print(counter)

    flds <- caret::createFolds(remains.y.response, k = 5, list = TRUE, returnTrain = FALSE)
    fold.i = flds[[1]]
    vali.set.x = remains.x.variables[fold.i,]
    vali.set.y = remains.y.response[fold.i]
    train.set.x = remains.x.variables[-fold.i,]
    train.set.y = remains.y.response[-fold.i]

    train <- train.set.x %>% as.matrix()
    vali <- vali.set.x %>% as.matrix()

    y_train <- train.set.y
    y_vali <- vali.set.y

    ################################################################################
    #### SVM
    ################################################################################
    library(e1071)
    library(caret)

    model.type = 'SVM'
    dir.create(paste0('figures/model.iter.', seed.i, '/',  model.type, '/', sex.i, '/'), recursive = T)



    train <- train.set.x %>% as.matrix()
    vali = vali.set.x %>% as.matrix()
    test <- test.set.x %>% as.matrix()

    y_train <- train.set.y
    y_vali = vali.set.y
    y_test <- test.set.y
    #
    train.comb = cbind(train, age = y_train)

    hyper_grid <- expand.grid(
      # cost = 2^seq(-5,5,1),
      # gamma= 2^seq(-5,5,1)
      # gamma = seq(0.0001, 0.002, 0.0001),
      gamma = seq(0.00001, 0.002, 0.00001),
      cost =  seq(0.5, 7, 0.5)
      # cost =  seq(7, 15, 0.5)
    )


    svm1 = read_csv(paste0('figures/model.iter.12/SVM/both/both.model.', test.fold.i, '.csv'))
    sorted.i = svm1 %>% arrange(rmse)
    sorted.i
    # sub.sorted = sorted.i %>%  filter(rmse < signif(min(rmse), digits = 2))
    # sub.sorted

    # for (i in nrow(sub.sorted)){
    # iter.coeff.res = lapply(1:20,
    iter.coeff.res = mclapply(1:20, mc.cores = parallel::detectCores() - 1,
                              function(i){
      # j =  sorted.i[i, 'j', drop = T]
      # j

      # browser()
      m_svm_untuned <- svm(
        formula = age ~ .,
        data    = train.comb,
        gamma = sorted.i$gamma[i],
        cost = sorted.i$cost[i]
      )

      ################################################################################
      # ### vip method
      library(vip)
      library(future)    # For parallel backend
      library(furrr)     # For parallel mapping (used internally by vip)

      # 1. Set a parallel plan *BEFORE* calling vi_permute
      # Choose one of the following plans:

      # For multicore on Linux/Mac (most efficient)
      plan(multisession, workers = parallel::detectCores() - 1) # Leave one core free

      perm_importance <- vip::vip(m_svm_untuned,
                             method = "permute",
                             train = train.comb,
                             target = train.comb[, 'age'],
                             metric = "rmse",
                             nsim = 10,
                             pred_wrapper = function(object, newdata) {
                               predict(object, newdata)
                             }, num_features = ncol(train.comb) - 1, keep = TRUE)
      perm_importance

      perm_importance
      # Get importance values
      importance_df <- as.data.frame(perm_importance$data)
      importance_df
      importance_df <- importance_df[order(-importance_df$Importance), ]
      importance_df
      # print(importance_df)
      coeff.folder = (paste0('figures/model.iter.', seed.i, '/',  model.type, '/', sex.i, '/model.coeff.', test.fold.i, '/'))
      dir.create(coeff.folder, recursive = T)
      write.csv(importance_df, paste0(coeff.folder, 'iter.', i, '.csv'))


      # ################################################################################
      # # 2. Define a custom metric function (e.g., RMSE for regression)
      # metric_rmse <- function(actual, predicted) {
      #   sqrt(mean((actual - predicted)^2))
      # }
      #
      # # 3. Calculate Permutation Importance
      # # Note: This can be computationally expensive as it re-predicts many times.
      # perm_importance <- vip::vi_permute(
      #   object = m_svm_untuned,          # Your trained model
      #   train = train.comb,          # Training data (without target)
      #   target = train.comb[,'age'],  # Name of the target variable
      #   metric = metric_rmse,        # Custom metric function
      #   type = "difference",         # How to compare: 'difference' or 'ratio'
      #   nsim = 10,                   # Number of permutations to average over
      #   keep = TRUE                  # Allows for progress bar
      # )
      #
      #
      # # Permutation importance function for regression
      # permutation_importance_regression <- function(model, X, y, metric = "mse", n_repeats = 10) {
      #
      #   browser()
      #   # Baseline performance
      #   baseline_pred <- predict(model, X)
      #   baseline_perf <- calculate_performance(baseline_pred, y, metric)
      #
      #   importance_scores <- matrix(0, nrow = n_repeats, ncol = ncol(X))
      #   colnames(importance_scores) <- colnames(X)
      #
      #   for (i in 1:ncol(X)) {
      #     cat("Processing feature:", colnames(X)[i], "\n")
      #
      #     for (j in 1:n_repeats) {
      #       # Shuffle the feature
      #       X_shuffled <- X
      #       X_shuffled[, i] <- sample(X_shuffled[, i])
      #
      #       # Calculate performance with shuffled feature
      #       shuffled_pred <- predict(model, X_shuffled)
      #       shuffled_perf <- calculate_performance(shuffled_pred, y, metric)
      #
      #       # For error metrics (MSE, RMSE), higher is worse, so we take the increase
      #       importance_scores[j, i] <- shuffled_perf - baseline_perf
      #     }
      #   }
      #
      #   return(importance_scores)
      # }
      #
      # # Performance metric calculation
      # calculate_performance <- function(predicted, actual, metric) {
      #   switch(metric,
      #          "mse" = mean((actual - predicted)^2),
      #          "rmse" = sqrt(mean((actual - predicted)^2)),
      #          "mae" = mean(abs(actual - predicted)),
      #          "r2" = cor(actual, predicted)^2,
      #          stop("Unknown metric")
      #   )
      # }
      #
      # # Get feature importance
      # importance_scores <- permutation_importance_regression(m_svm_untuned, test.set.x, test.set.y, metric = 'rmse')
      # importance_scores %>% view
      # mse_importance <- create_importance_df(importance_scores, "mse")
      #
      # importance_scores
      # feature_names <- colnames(X)
      #
      # importance_df <- data.frame(
      #   feature = feature_names,
      #   importance = importance_scores
      # ) %>% arrange(desc(importance))
      #
      # print(importance_df)

      print(i)
    })

  }
}

