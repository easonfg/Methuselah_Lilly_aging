rm(list = ls())

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

  full.data.wide
  ### scaling data
  full.data.wide <- full.data.wide %>%
    mutate(across(
      # Select columns to scale: everything EXCEPT Name, Sex, Age
      .cols = -c(SampleName, Sex, Age2),
      # Apply the scale function to each selected column
      .fns = scale
    ))

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

  full.data.wide$Sex = ifelse(full.data.wide$Sex == 'M', 1, 0)

  adjusted.wide = full.data.wide
  x.variables = adjusted.wide %>% select(-c(SampleName, Age2, ))
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
    train
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
    train.comb = cbind(remains.x.variables, age = remains.y.response)
#
#     m_svm_untuned <- svm(formula = age ~ .,
#                          data    = train.comb)
#
#     pred_svm_untuned <- predict(m_svm_untuned, test)
#
#     yhat <- pred_svm_untuned
#     y <- y_test
#     svm_stats_untuned <- postResample(yhat, y)
#     svm_stats_untuned

    # #####################################
    # #grid search
    # #create hyperparameter grid
    # hyper_grid <- expand.grid(
    #   # cost = 2^seq(-5,5,1),
    #   # gamma= 2^seq(-5,5,1)
    #   # gamma = seq(0.0001, 0.002, 0.0001),
    #   gamma = seq(0.00001, 0.0015, 0.00001),
    #   # cost =  seq(0.5, 7, 0.5)
    #   cost =  seq(7, 15, 0.5)
    #
    # )
    # hyper_grid
    # e <- NULL
    #
    # print(dim(hyper_grid))
    # # for(j in 1:nrow(hyper_grid)){
    # # e.res = mclapply(1:10, function(j){
    # e.res = mclapply(1:nrow(hyper_grid),
    #                  mc.cores = parallel::detectCores() - 1,
    #                  function(j){
    #
    #   m_svm_untuned <- svm(
    #     formula = age ~ .,
    #     data    = train.comb,
    #     gamma = hyper_grid$gamma[j],
    #     cost = hyper_grid$cost[j]
    #   )
    #
    #   pred_svm_untuned <-predict(
    #     m_svm_untuned,
    #     newdata = vali,
    #   )
    #
    #   yhat <- pred_svm_untuned
    #   y <- y_vali
    #   postResample(yhat, y)
    #   e[j] <- postResample(yhat, y)[1]
    #   cat(j, "\n")
    #   postResample(yhat, y)[1]
    #
    # })
    # # }
    # e.res.list = e.res %>% unlist()
    # e.res.list
    #
    # # which.min(e)  #minimum MSE
    # # e[which.min(e)]
    # # hyper_grid[which.min(e),]
    # #
    # # j = which.min(e)  #minimum MSE
    #
    #
    # j = 1
    # m_svm_untuned <- svm(
    #   formula = age ~ .,
    #   data    = train.comb,
    #   gamma = hyper_grid$gamma[j],
    #   cost = hyper_grid$cost[j]
    # )
    # pred_svm_untuned = predict(m_svm_untuned, newdata = test)
    # yhat <- pred_svm_untuned
    #
    # pdf(paste0('figures/classifications.', seed.i, '/',  model.type, '/', sex.i, '/', output.var,
    #            '.', test.fold.i, '.scatterplot.pdf'), h = 5, w = 5)
    # test.rmse = RMSE(y_test,pred_svm_untuned)
    # test.cor = cor(y_test,pred_svm_untuned)
    # p = data.frame(manual.pred = pred_svm_untuned, chronological.age = test.set.y) %>%
    #   ggplot(aes(chronological.age, manual.pred)) +
    #   geom_point() +
    #   geom_abline(intercept = 0, slope = 1, color = "red")  +
    #   annotate("text", x = -Inf, y = Inf, label = paste0('RMSE: ', round(test.rmse, 3)),
    #            hjust = -0.1, vjust = 2, size = 5, color = "blue") +
    #   annotate("text", x = -Inf, y = Inf, label = paste0('Pearson Cor: ', round(test.cor, 3)),
    #            hjust = -0.1, vjust = 3.5, size = 5, color = "blue") +
    #   xlab('Chronological Age') + ylab('Predicted Age') +
    #   theme_bw(base_size =15)
    # print(p)
    # dev.off()

#
#     model = lm(score ~ age, data = data.frame(score = pred_svm_untuned, age = test.set.y))
#     plot(pred_svm_untuned, test.set.y)
#     adjusted_scores <- residuals(model)
#     # adjusted_scores = adjusted_scores + data.score.mortality$Age2

    mortality.data = read.xlsx('data/dataverse_files/BoAC_plasma_metadata_TTD.xlsx') %>% as_tibble()
    mortality.data
    data.score.mortality = test.set.meta %>%     mutate(Plasma.ID = str_extract(SampleName, "[^_]+$")) %>%
      left_join(mortality.data %>% select(Plasma.ID, Vital), by = 'Plasma.ID')

    # data.score.mortality$adj.score = adjusted_scores

    library(pROC)

    # roc_obj <- roc(data.score.mortality$Vital, data.score.mortality$adj.score)
    # auc(roc_obj)

    # iter.model.res = lapply(1:nrow(hyper_grid), function(j){
    hyper_grid <- expand.grid(
      # cost = 2^seq(-5,5,1),
      # gamma= 2^seq(-5,5,1)
      # gamma = seq(0.0001, 0.002, 0.0001),
      gamma = seq(0.00001, 0.002, 0.00001),
      cost =  seq(0.5, 7, 0.5)
      # cost =  seq(7, 15, 0.5)

    )
    hyper_grid
    # browser()
    j = 1

    iter.model.res = mclapply(1:nrow(hyper_grid),
    # iter.model.res = mclapply(1,
                              mc.cores = parallel::detectCores() - 1,
                              function(j){
      print(j)
      m_svm_untuned <- svm(
        formula = age ~ .,
        data    = train.comb,
        gamma = hyper_grid$gamma[j],
        cost = hyper_grid$cost[j]
      )


      pred_svm_untuned <-predict(
        m_svm_untuned,
        # newdata = vali,
        newdata = train.comb,
      )

      yhat <- pred_svm_untuned
      y <- y_vali
      y <- train.comb$age
      e =       postResample(yhat, y)[1]
      e

      pred_svm_untuned = predict(m_svm_untuned, newdata = test.set.x)
      pred_svm_untuned
      test.rmse = RMSE(y_test,pred_svm_untuned)
      test.cor = cor(y_test,pred_svm_untuned)

      model = lm(score ~ age , data = data.frame(score = pred_svm_untuned, age = test.set.y))
      adjusted_scores <- residuals(model)
      data.score.mortality$adj.score = adjusted_scores
      roc_obj <- roc(data.score.mortality$Vital, data.score.mortality$adj.score)
      auc.j = auc(roc_obj) %>% as.numeric()

      data.frame(j = j, error = e, rmse = test.rmse, cor = test.cor, auc = auc.j)
    })
    iter.model = iter.model.res %>% data.table::rbindlist()
    iter.model = cbind(iter.model, hyper_grid)
    iter.model
    iter.model %>% filter(rmse == min(rmse))
    iter.model %>% arrange(auc) %>%
      mutate(a = j) %>%
      mutate(a = fct_inorder(as.character(a))) %>%
      ggplot() +
      geom_point(aes(a, rmse), color = 'green') +
      geom_point(aes(a, error), color = 'red')

    write.csv(iter.model, paste0('figures/model.iter.', seed.i, '/',  model.type, '/', sex.i, '/', output.var,
               '.model.', test.fold.i, '.csv'))
    # iter.model
    # iter.model$auc = iter.model$auc %>% as.numeric()
    # iter.model %>% filter(auc == max(auc))
    # iter.model %>% filter(error == min(error))
    # iter.model %>% filter(rmse == min(rmse))
    # iter.model = cbind(iter.model, hyper_grid)
    # iter.model
    #
    # iter.model %>%
    #   mutate(scale.auc = scale(auc), scale.error = scale(error),
    #          scale.rmse = scale(rmse)) %>%
    #   mutate(comb.score = (1 - scale.auc) + scale.error) %>%
    #   group_by(cost) %>%
    #   slice_min(comb.score, n = 1) %>% select(-c(scale.auc, scale.error, scale.rmse))
    #
    # iter.model %>%
    #   mutate(scale.auc = scale(auc), scale.error = scale(error),
    #          scale.rmse = scale(rmse)) %>%
    #   mutate(comb.score = (1 - auc) + scale.error) %>%
    #   slice_min(comb.score, n = 1) %>% view
    #
    #
    # iter.model %>%
    #   # pivot_longer(cols = c(rmse, auc), names_to = 'metric') %>%
    #   ggplot()+
    #   # geom_line(aes(gamma, auc))
    #   # geom_point(aes(gamma, auc, color = cost))
    #   geom_point(aes(gamma, error, color = cost))
    #
    # iter.model %>% head(110)
    # iter.model.long = iter.model %>%
    #   pivot_longer(cols = c(rmse, auc, error), names_to = 'metric')
    #
    # best.metrics.row = rbind(iter.model.long %>%
    #   filter(metric == 'auc') %>%
    #   group_by(cost, metric) %>%
    #   slice_max(value) %>% ungroup(),
    # iter.model.long %>%
    #   filter(metric != 'auc') %>%
    #   group_by(cost, metric) %>%
    #   slice_min(value) %>% ungroup()
    # )
    # best.metrics.row
    #
    #  iter.model %>%
    #   # .[1:1000,] %>%
    #   mutate(auc = auc * 100) %>%
    #   # filter(cost == 6) %>%
    #   # filter(cost == 15) %>%
    #   pivot_longer(cols = c(rmse, auc, error), names_to = 'metric') %>%
    #   # mutate(value = scale(value)) %>%
    #   ggplot()+
    #   # geom_line(aes(gamma, auc))
    #   # geom_point(aes(gamma, rmse, color = cost))
    #   # geom_point(aes(j, rmse, color = cost))
    #   # geom_point(aes(j, rmse, color = gamma))
    #   # geom_point(aes(j, auc, color = gamma))
    #   geom_point(aes(j, value, color = metric)) +
    #    geom_vline(data = best.metrics.row
    #               # %>% filter(cost == 15)
    #               , aes(xintercept = j, color = metric))
    #
    #  iter.model %>%
    #    # .[1:1000,] %>%
    #    mutate(auc = auc * 100) %>%
    #    # filter(cost == 6) %>%
    #    # filter(cost == 3) %>%
    #    pivot_longer(cols = c(rmse, auc, error), names_to = 'metric') %>%
    #    # mutate(value = scale(value)) %>%
    #    ggplot()+
    #    # geom_line(aes(gamma, auc))
    #    # geom_point(aes(gamma, rmse, color = cost))
    #    # geom_point(aes(j, rmse, color = cost))
    #    # geom_point(aes(j, rmse, color = gamma))
    #    # geom_point(aes(j, auc, color = gamma))
    #    geom_point(aes(j, value, color = metric)) +
    #    geom_vline(data = best.metrics.row , aes(xintercept = j, color = metric))
    #
    # iter.model
    # iter.model %>%
    #   ggplot()+
    #   # geom_line(aes(error, auc))
    #   # geom_line(aes(error, rmse))
    #   # geom_point(aes(j, error))
    #   geom_line(aes(j, auc))
    #   # geom_line(aes(j, rmse))
    #   # geom_line(aes(j, error))
    #
    #
    #
    # svm.roc[[test.fold.i]] = roc_obj
    #
    # data.score.mortality = data.score.mortality %>% mutate(Vital = ifelse(Vital == 'Deceased', 1, 0))
    #
    # age.model <- glm(Vital ~ Age2,
    #                  data = data.score.mortality,
    #                  family = binomial(link = "logit"))
    #
    # age.score.model <- glm(Vital ~ adj.score + Age2,
    #                        data = data.score.mortality,
    #                        family = binomial(link = "logit"))
    # score.model <- glm(Vital ~ adj.score,
    #                    data = data.score.mortality,
    #                    family = binomial(link = "logit"))
    #
    # data.score.mortality$age.predicted_prob <- predict(age.model, type = "response")
    # data.score.mortality$age.score.predicted_prob <- predict(age.score.model, type = "response")
    # data.score.mortality$score.predicted_prob <- predict(score.model, type = "response")
    #
    # library(pROC)
    # age.roc_obj <- roc(data.score.mortality$Vital, data.score.mortality$age.predicted_prob)
    # # age.ci <- ci.se(age.roc_obj, specificities = seq(0, 1, 0.01))  # CI at 100 points
    # age.score.roc_obj <- roc(data.score.mortality$Vital, data.score.mortality$age.score.predicted_prob)
    # # age.score.ci <- ci.se(age.score.roc_obj, specificities = seq(0, 1, 0.01))  # CI at 100 points
    #
    # pdf(paste0('figures/classifications.', seed.i, '/',  model.type, '/', sex.i, '/', output.var,
    #            '.', test.fold.i, '.auc.pdf'), h = 5, w = 5)
    # p = ggroc(list(age = age.roc_obj, age.score = age.score.roc_obj ), legacy.axes = TRUE) +
    #   annotate("text", x = 0.6, y = 0.3,
    #            label = paste("Age AUC =", round(auc(age.roc_obj), 3)), color = "red") +
    #   annotate("text", x = 0.6, y = 0.2,
    #            label = paste("Age + score AUC =", round(auc(age.score.roc_obj), 3)), color = "blue")+
    #   scale_fill_manual(values = c("age" = "red", "age.score" = "blue")) +
    #   labs(
    #     # title = "Comparison of ROC Curves with 95% Confidence Intervals",
    #     x = "1 - Specificity (False Positive Rate)",
    #     y = "Sensitivity (True Positive Rate)"
    #   ) +
    #   theme_minimal() +
    #   theme(legend.position = "bottom")
    # print(p)
    # dev.off()
    #
    # roc.summary = roc.test(age.roc_obj, age.score.roc_obj)
    # roc.summary$p.value
    #
    # lrt.test = anova(age.model, age.score.model, test = "LRT")  # p < 0.05 → significant improvement
    # lrt.test
    # # lrt.test = anova(score.model, age.score.model, test = "LRT")  # p < 0.05 → significant improvement
    # svm.stats.res[[test.fold.i]] = data.frame(auc.age = roc.summary$roc1$auc,
    #                                                auc.age.score = roc.summary$roc2$auc,
    #                                                delong.pval = roc.summary$p.value,
    #                                                lrt.pval = lrt.test$`Pr(>Chi)`[2])


    }

  # }
}


#
#
#
# plot.auc = function(all_rocs, model.name){
#   # First calculate AUCs for each fold
#   auc_values <- sapply(all_rocs, auc)
#
#   # Create named vector for legend labels
#   legend_labels <- paste0("Fold ", seq_along(all_rocs), " (AUC = ", round(auc_values, 3), ")")
#   names(all_rocs) <- legend_labels
#
#   # Calculate mean ROC curve
#   mean_roc <- roc(
#     response = unlist(lapply(all_rocs, function(x) x$response)),
#     predictor = unlist(lapply(all_rocs, function(x) x$predictor))
#   )
#   mean_roc
#
#   all_rocs
#   all_rocs.mean = all_rocs
#   all_rocs.mean$mean = mean_roc
#   names(all_rocs.mean)[length(all_rocs.mean)] = paste0("Mean ", " (AUC = ", round(mean_roc$auc, 3), ")")
#
#   # Create combined plot
#   pdf(paste0('figures/', model.name, '.roc.plots.pdf'))
#   p = ggroc(all_rocs.mean, aes = c("color", "linetype")) +
#     geom_line(size = 1) +  # Make lines thicker
#     geom_abline(slope = 1, intercept = 1, linetype = "dashed", color = "gray") +
#     scale_color_manual(
#       name = "CV Fold",
#       values = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", 'red')  # Colorblind-friendly
#     ) +
#     scale_linetype_manual(
#       name = "CV Fold",
#       values = rep("solid", length(all_rocs.mean))  # All solid lines (or vary if preferred)
#     ) +
#     labs(
#       title = "ROC Curves from Cross-Validation",
#       x = "False Positive Rate (1 - Specificity)",
#       y = "True Positive Rate (Sensitivity)"
#     ) +
#     theme_minimal() +
#     theme(
#       legend.position = c(0.8, 0.2),
#       plot.title = element_text(hjust = 0.5, face = "bold")
#     ) +
#     coord_equal()  # Ensure 1:1 aspect ratio for proper ROC visualization
#   print(p)
#   dev.off()
#
#
#
# }
#
#
# plot.auc(lightgbm.roc, 'lightGBM')
# plot.auc(svm.roc, 'svm')
#
# light.stats.df = lightgbm.stats.res %>% data.table::rbindlist()
# light.stats.df %>% write.csv(.,
#                        paste0('figures/classifications.', seed.i, '/',  'lightgbm', '/', sex.i, '/', output.var,
#                               '.stats.csv'))
#
#
# svm.stats.df = svm.stats.res %>% data.table::rbindlist()
# svm.stats.df %>% write.csv(.,
#                              paste0('figures/classifications.', seed.i, '/',  'svm', '/', sex.i, '/', output.var,
#                                     '.stats.csv'))
