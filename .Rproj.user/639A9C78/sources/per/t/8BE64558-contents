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
library(glmnet)

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

    model.type = 'elastic.net'
    dir.create(paste0('figures/model.iter.', seed.i, '/',  model.type, '/', sex.i, '/'), recursive = T)



    train <- train.set.x %>% as.matrix()
    vali = vali.set.x %>% as.matrix()
    test <- test.set.x %>% as.matrix()

    y_train <- train.set.y
    y_vali = vali.set.y
    y_test <- test.set.y
    #
    train.comb = cbind(train, age = y_train)

    tune_grid <- expand.grid(
      alpha = seq(0, 1, by = 0.1),  # Test alpha from 0 (Ridge) to 1 (Lasso)
      lambda = 10^seq(-3, 2, length = 100)  # Test lambda from 0.001 to 10
    )


    svm1 = read_csv(paste0('figures/model.iter.12/elastic.net/both/both.model.', test.fold.i, '.csv'))
    sorted.i = svm1 %>% arrange(rmse)
    sorted.i
    # sub.sorted = sorted.i %>%  filter(rmse < signif(min(rmse), digits = 2))
    # sub.sorted

    # for (i in nrow(sub.sorted)){
    iter.coeff.res = lapply(1:20,
    # iter.coeff.res = mclapply(1:20, mc.cores = parallel::detectCores() - 1,
                              function(i){
                                # j =  sorted.i[i, 'j', drop = T]
                                # j

                                # browser()
                                glm.manual = glmnet(
                                  x = remains.x.variables,
                                  y = remains.y.response,
                                  alpha = sorted.i$alpha[i],
                                  lambda = sorted.i$lambda[i]
                                )

                                return_features = function(coeff){
                                  top_features = data.frame(name = coeff@Dimnames[[1]][coeff@i + 1], coefficient = coeff@x)
                                  top_features = top_features[order(top_features$coefficient),]
                                  # print(top_features)
                                  return(top_features)
                                }

                                features = return_features(coef(glm.manual, s = 'lambda.1se'))
                                features

                                coeff.folder = (paste0('figures/model.iter.', seed.i, '/',  model.type, '/', sex.i, '/model.coeff.', test.fold.i, '/'))
                                dir.create(coeff.folder, recursive = T)
                                write.csv(features, paste0(coeff.folder, 'iter.', i, '.csv'))

                              })

  }
}

