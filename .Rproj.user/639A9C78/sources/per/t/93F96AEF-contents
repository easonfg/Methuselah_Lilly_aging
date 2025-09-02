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

  # ## female
  # adjusted.wide.f = full.data.wide.f
  # x.variables.f = adjusted.wide.f %>% select(-c(SampleName, Age2, Sex))
  # logical.complete.cases = complete.cases(x.variables.f)
  # x.variables.f = x.variables.f[logical.complete.cases,]
  # y.response.f = adjusted.wide.f$Age2
  # y.response.f = y.response.f[logical.complete.cases]
  #
  # ## male
  # adjusted.wide.m = full.data.wide.m
  # x.variables.m = adjusted.wide.m %>% select(-c(SampleName, Age2, Sex))
  # logical.complete.cases = complete.cases(x.variables.m)
  # x.variables.m = x.variables.m[logical.complete.cases,]
  # y.response.m = adjusted.wide.m$Age2
  # y.response.m = y.response.m[logical.complete.cases]

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

  # unscaled.glmnet.pred = glmnet_pred * attr(scaled.wide, 'scaled:scale')[['Age2']] + attr(scaled.wide, 'scaled:center')[['Age2']]
  # unscaled.org.y = y.response * attr(scaled.wide, 'scaled:scale')[['Age2']] + attr(scaled.wide, 'scaled:center')[['Age2']]
  # sqrt(mean((unscaled.glmnet.pred - unscaled.org.y)^2))




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

    model.type = 'lightGBM'
    dir.create(paste0('figures/model.iter.', seed.i, '/',  model.type, '/', sex.i, '/'), recursive = T)

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

    train_lgb <- lightgbm::lgb.Dataset(as.matrix(train),label=y_train)
    vali_lgb <- lightgbm::lgb.Dataset.create.valid(train_lgb,vali,label = y_vali)

    ### iter light gbm
    mortality.data = read.xlsx('data/dataverse_files/BoAC_plasma_metadata_TTD.xlsx') %>% as_tibble()
    mortality.data
    data.score.mortality = test.set.meta %>%     mutate(Plasma.ID = str_extract(SampleName, "[^_]+$")) %>%
      left_join(mortality.data %>% select(Plasma.ID, Vital), by = 'Plasma.ID')

    # # #################
    # #grid search
    # #create hyperparameter grid
    # num_leaves =seq(20,28,1)
    # max_depth = round(log(num_leaves) / log(2),0)
    # num_iterations = seq(200,400,50)
    # early_stopping_rounds = round(num_iterations * .1,0)
    #
    #
    # hyper_grid <- expand.grid(max_depth = max_depth,
    #                           num_leaves =num_leaves,
    #                           num_iterations = num_iterations,
    #                           early_stopping_rounds=early_stopping_rounds,
    #                           learning_rate = seq(.45, .50, .005))

    #
    # # Define search grid
    # hyper_grid = expand.grid(
    #   num_leaves = c(20, 30, 40),
    #   max_depth = c(-1, 5, 10),
    #   learning_rate = c(0.01, 0.05),
    #   feature_fraction = c(0.6, 0.8),
    #   lambda_l1 = c(0, 0.1),
    #   lambda_l2 = c(0, 0.1),
    #   min_data_in_leaf = c(20, 50)
    # )
    # hyper_grid
    # hyper_grid %>% dim()

    # Define parameter grid for grid search
    param_grid <- expand.grid(
      num_leaves = c(15, 31, 63, 127),          # Controls model complexity
      max_depth = c(3, 5, 7, 9, -1),           # -1 means no limit
      learning_rate = c(0.01, 0.05, 0.1, 0.2), # Step size shrinkage
      n_estimators = c(100, 200, 500),         # Number of boosting iterations
      min_data_in_leaf = c(20, 50, 100),       # Minimum samples in leaf
      feature_fraction = c(0.6, 0.8, 1.0),     # Fraction of features used
      bagging_fraction = c(0.6, 0.8, 1.0),     # Fraction of data used
      bagging_freq = c(0, 5, 10),              # Frequency for bagging
      lambda_l1 = c(0, 0.1, 1, 10),            # L1 regularization
      lambda_l2 = c(0, 0.1, 1, 10)             # L2 regularization
    )
    param_grid

    browser()
    # Bayesian optimization for LightGBM
    lightgbm_bayesian_optimization <- function(X, y, n_iter = 30) {

      # Define parameter bounds
      bounds <- list(
        num_leaves = c(15L, 255L),
        max_depth = c(3L, 12L),
        learning_rate = c(0.01, 0.3),
        n_estimators = c(100L, 1000L),
        min_data_in_leaf = c(5L, 100L),
        feature_fraction = c(0.5, 1.0),
        lambda_l1 = c(0, 10),
        lambda_l2 = c(0, 10)
      )

      # Objective function
      obj_func <- function(num_leaves, max_depth, learning_rate,
                           n_estimators, min_data_in_leaf,
                           feature_fraction, lambda_l1, lambda_l2) {

        params <- list(
          objective = "regression",
          metric = "rmse",
          num_leaves = as.integer(num_leaves),
          max_depth = as.integer(max_depth),
          learning_rate = learning_rate,
          n_estimators = as.integer(n_estimators),
          min_data_in_leaf = as.integer(min_data_in_leaf),
          feature_fraction = feature_fraction,
          lambda_l1 = lambda_l1,
          lambda_l2 = lambda_l2,
          verbose = -1
        )
        # params <- list(
        #   objective = "regression",
        #   metric = "rmse",
        #   num_leaves = as.integer(param_grid$num_leaves[1]),
        #   max_depth = as.integer(param_grid$max_depth[1]),
        #   learning_rate = param_grid$learning_rate[1],
        #   n_estimators = as.integer(param_grid$n_estimators[1]),
        #   min_data_in_leaf = as.integer(param_grid$min_data_in_leaf[1]),
        #   feature_fraction = param_grid$feature_fraction[1],
        #   lambda_l1 = param_grid$lambda_l1[1],
        #   lambda_l2 = param_grid$lambda_l2[1],
        #   verbose = -1
        # )

        cv <- lightgbm::lgb.cv(
          params = params,
          data = lightgbm::lgb.Dataset(data = as.matrix(X), label = y),
          nfold = 5,
          early_stopping_rounds = 20
        )
        # all_cv_results <- data.table::as.data.table(cv$record_evals$valid$rmse)
        # all_cv_results
        # cv$record_evals$valid$rmse
        # cv$best_iter
        # cv$best_score

          # browser()
        return(list(Score = -cv$best_score))  # Negative because we want to minimize loss
      }

      # Run optimization
      opt_results <- ParBayesianOptimization::bayesOpt(
        FUN = obj_func,
        bounds = bounds,
        initPoints = 10,
        iters.n = n_iter,
        acq = "ei"
      )
      browser()

      return(opt_results)
    }

    # Usage example
    # opt_results <- lightgbm_bayesian_optimization(train, y_train, n_iter = 30)
    opt_results <- lightgbm_bayesian_optimization(train, y_train, n_iter = 1)
    best_params <- ParBayesianOptimization::getBestPars(opt_results)
    best_params
    all_results <- ParBayesianOptimization::getScore(opt_results)
    opt_results$scoreSummary %>% dim()
    opt_results$scoreSummary

    best_params


    #
    # create_optimal_lgb_regression_grid <- function(n_combinations = 50) {
    #
    #   # # Adjust grid size based on data size if provided
    #   # if (!is.null(data_size)) {
    #   #   if (data_size < 1000) n_combinations <- 20
    #   #   else if (data_size < 10000) n_combinations <- 30
    #   #   else n_combinations <- 50
    #   # }
    #
    #   # Create intelligent parameter grid for regression
    #   param_grid <- expand.grid(
    #     learning_rate = exp(seq(log(0.005), log(0.2), length.out = 8)),
    #     num_leaves = round(2^seq(4, 8, length.out = 6)),
    #     min_data_in_leaf = round(exp(seq(log(5), log(50), length.out = 6))),
    #     feature_fraction = seq(0.6, 1.0, length.out = 5),
    #     bagging_fraction = seq(0.6, 1.0, length.out = 5),
    #     bagging_freq = c(0, 1, 3, 5),
    #     lambda_l1 = exp(seq(log(0.001), log(5), length.out = 5)),
    #     lambda_l2 = exp(seq(log(0.001), log(5), length.out = 5)),
    #     min_gain_to_split = exp(seq(log(0.001), log(0.5), length.out = 4)),
    #     max_depth = c(-1, 6, 8, 10, 12),
    #     stringsAsFactors = FALSE
    #   )
    #
    #   # Sample intelligently
    #   if (nrow(param_grid) > n_combinations) {
    #     # set.seed(123)
    #     sampled_indices <- sample(seq_len(nrow(param_grid)), n_combinations)
    #     param_grid <- param_grid[sampled_indices, ]
    #   }
    #
    #   return(param_grid)
    # }
    #
    # # Create regression grid
    # hyper_grid <- create_optimal_lgb_regression_grid(n_combinations = 6000)
    # hyper_grid %>% dim()

    hyper_grid %>% distinct()
    # iter.model.res = lapply(1:nrow(hyper_grid),
    iter.model.res = mclapply(1:nrow(hyper_grid),
                              mc.cores = parallel::detectCores() - 1,
                              function(j){
                                # for (j in 1:nrow(hyper_grid)) {
                                # print(j)
                                light_gbn_tuned <- lightgbm::lgb.train(
                                  params = list(
                                    objective = "regression",
                                    metric = "rmse",
                                    max_depth = hyper_grid$max_depth[j],
                                    num_leaves =hyper_grid$num_leaves[j],
                                    num_iterations = hyper_grid$num_iterations[j],
                                    early_stopping_rounds=hyper_grid$early_stopping_rounds[j],
                                    verbose = -1,
                                    learning_rate = hyper_grid$learning_rate[j]
                                    #feature_fraction = .9
                                  ),
                                  valids = list(vali = vali_lgb),
                                  data = train_lgb
                                )

                                yhat_fit_tuned <- predict(light_gbn_tuned,train)
                                yhat_predict_tuned <- predict(light_gbn_tuned,vali)

                                rmse_fit.j <- RMSE(y_train,yhat_fit_tuned)
                                rmse_predict.j <- RMSE(y_vali,yhat_predict_tuned)
                                cat(j, "\n")

                                yhat_predict_final <- predict(light_gbn_tuned,test)
                                rmse_predict_final <- RMSE(y_test,yhat_predict_final)
                                test.cor = cor(y_test,yhat_predict_final)



                                model = lm(score ~ age, data = data.frame(score = yhat_predict_final, age = test.set.y))
                                adjusted_scores <- residuals(model)

                                data.score.mortality$adj.score = adjusted_scores


                                library(pROC)

                                roc_obj <- roc(data.score.mortality$Vital, data.score.mortality$adj.score)
                                auc(roc_obj)
                                auc.j = auc(roc_obj) %>% as.numeric()
                                auc.j

                                data.frame(j = j, error = rmse_fit.j, rmse = rmse_predict_final, cor = test.cor, auc = auc.j)
                              })
    iter.model = iter.model.res %>% data.table::rbindlist()
    iter.model
    iter.model = cbind(iter.model, hyper_grid)

    write.csv(iter.model, paste0('figures/model.iter.', seed.i, '/',  model.type, '/', sex.i, '/', output.var,
                                 '.model.', test.fold.i, '.csv'))
    ### lightGBM
    ################################################################################
    ################################################################################

  }

  # }
}


#
