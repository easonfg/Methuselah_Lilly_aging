

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


####where to select rmse, lowest rmse, random, or highest rmse
# for (region.i in c('top')){
for (region.i in c( 'bottom', 'random', 'top')){
# for (region.i in c('top', 'bottom', 'random')){

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

      # for (model.type in c('elastic.net', 'SVM', 'lightGBM'))
      for (model.type in c("SVM", 'elastic.net', 'lightGBM')) {
        train <- train.set.x %>% as.matrix()
        vali = vali.set.x %>% as.matrix()
        test <- test.set.x %>% as.matrix()

        y_train <- train.set.y
        y_vali = vali.set.y
        y_test <- test.set.y
        #
        train.comb = cbind(train, age = y_train)


        svm1 = read_csv(paste0('figures/model.iter.12/', model.type, '/both/both.model.', test.fold.i, '.csv'))

        if (region.i == 'top'){
          sorted.i = svm1 %>% arrange(rmse)
        } else if (region.i == 'bottom'){


          sorted.i = svm1 %>% arrange(desc(rmse))
        } else if (region.i == 'random') {

          sorted.i = svm1[sample(1:nrow(svm1), 20), ]
        }
        sorted.i

        # iter.coeff.res = lapply(1:20,
        iter.coeff.res = mclapply(1:20, mc.cores = parallel::detectCores() - 1,
          function(i){
            # j =  sorted.i[i, 'j', drop = T]
            # j

            # browser()

            ############################## elastic net##############################
            if (model.type == 'elastic.net') {
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
              features %>% filter(!str_detect(name, 'Intercept'))
              features

            } else if (model.type == 'SVM') {
              # For multicore on Linux/Mac (most efficient)
              # plan(multisession, workers = parallel::detectCores() - 1) # Leave one core free

              m_svm_untuned <- svm(
                formula = age ~ .,
                data    = train.comb,
                gamma = sorted.i$gamma[i],
                cost = sorted.i$cost[i]
              )

              perm_importance <- vip::vip(m_svm_untuned,
                                          method = "permute",
                                          train = train.comb,
                                          target = train.comb[, 'age'],
                                          metric = "rmse",
                                          nsim = 50,
                                          pred_wrapper = function(object, newdata) {
                                            predict(object, newdata)
                                          }, num_features = ncol(train.comb) - 1, keep = TRUE)

              importance_df <- as.data.frame(perm_importance$data)
              features <- importance_df[order(-importance_df$Importance), ]

            } else if (model.type == 'lightGBM') {
              train_lgb <- lightgbm::lgb.Dataset(as.matrix(train),label=y_train)
              vali_lgb <- lightgbm::lgb.Dataset.create.valid(train_lgb,vali,label = y_vali)

              light_gbn_tuned <- lightgbm::lgb.train(
                params = list(
                  objective = "regression",
                  metric = "rmse",
                  max_depth = sorted.i$max_depth[i],
                  num_leaves =sorted.i$num_leaves[i],
                  num_iterations = sorted.i$num_iterations[i],
                  early_stopping_rounds=sorted.i$early_stopping_rounds[i],
                  verbose = -1,
                  learning_rate = sorted.i$learning_rate[i]
                  #feature_fraction = .9
                ),
                valids = list(vali = vali_lgb),
                data = train_lgb
              )

              yhat_fit_tuned <- predict(light_gbn_tuned,train)
              yhat_predict_tuned <- predict(light_gbn_tuned,vali)
              rmse_fit.j <- RMSE(y_train,yhat_fit_tuned)
              rmse_fit.j

              features <- lightgbm::lgb.importance(light_gbn_tuned)
            }

            coeff.folder = (paste0('figures/model.iter.', region.i, '.rmse.', seed.i, '/',  model.type, '/', sex.i, '/model.coeff.', test.fold.i, '/'))
            dir.create(coeff.folder, recursive = T)
            write.csv(features, paste0(coeff.folder, 'iter.', i, '.csv'))

          })

      }



    }
  }

}

