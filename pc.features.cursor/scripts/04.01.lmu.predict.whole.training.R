rm(list = ls())

library(glmnet)
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

  full.data.wide$Sex = ifelse(full.data.wide$Sex == 'M', 1, 0)


  adjusted.wide = full.data.wide
  x.variables = adjusted.wide %>% select(-c(SampleName, Age2, ))
  logical.complete.cases = complete.cases(x.variables)
  x.variables = x.variables[logical.complete.cases,]
  y.response = adjusted.wide$Age2
  y.response = y.response[logical.complete.cases]
  meta.data = adjusted.wide %>% select(SampleName, Age2, Sex)


  test.fold <- caret::createFolds(y.response, k = 5, list = TRUE, returnTrain = FALSE)
  res.table = list()
  counter = 1
  for (test.fold.i in 1:6){
    # test.set.x = x.variables[test.fold[[test.fold.i]],]
    # test.set.y = y.response[test.fold[[test.fold.i]]]
    # remains.x.variables = x.variables[-test.fold[[test.fold.i]],]
    # remains.y.response = y.response[-test.fold[[test.fold.i]]]
    # test.set.meta = meta.data[test.fold[[test.fold.i]],]
    #
    #
    #
    # test = test.set.x %>% as.matrix()
    # y_test = test.set.y
    #
    # flds <- caret::createFolds(remains.y.response, k = 5, list = TRUE, returnTrain = FALSE)
    # fold.i = flds[[1]]
    # vali.set.x = remains.x.variables[fold.i,]
    # vali.set.y = remains.y.response[fold.i]
    # train.set.x = remains.x.variables[-fold.i,]
    # train.set.y = remains.y.response[-fold.i]
    #
    # # # train <- train.set.x %>% as.matrix()
    # # train <- x.variables %>% as.matrix()
    # # vali <- vali.set.x %>% as.matrix()
    # #
    # # # y_train <- train.set.y
    # # y_train <- y.response
    # # y_vali <- vali.set.y

    ################################################################################
    #### SVM
    ################################################################################
    library(e1071)
    library(caret)

    # model.type = 'SVM'
    # dir.create(paste0('figures/model.iter.', seed.i, '/',  model.type, '/', sex.i, '/'), recursive = T)

    # train <- train.set.x %>% as.matrix()
    train <- x.variables %>% as.matrix()
    # vali = vali.set.x %>% as.matrix()
    # test <- test.set.x %>% as.matrix()

    # y_train <- train.set.y
    y_train <- y.response
    # y_vali = vali.set.y
    # y_test <- test.set.y
    #
    train.comb = cbind(train, age = y_train)




    ################################################################################
    other1.cns = read.csv('data/other.data/LMU/20240906_P-000518-CNS_LMU_NULISAseq_TAP_Counts_Report_FULL.csv')
    other1.ifn = read.csv('data/other.data/LMU/20240906_P-000580-IF_LMU_NULISAseq_TAP_Counts_Report_FULL.csv')
    other1.cns %>% head()
    other1.ifn %>% head()

    matrix.type = 'Serum'
    # matrix.type = 'CSF'
    otherdata = rbind(
      other1.cns %>% select(SampleName, Target, NPQ) %>% filter(str_detect(SampleName, matrix.type)),
      other1.ifn %>% select(SampleName, Target, NPQ) %>% filter(str_detect(SampleName, matrix.type))
    )
    otherdata$SampleName
    str_extract(otherdata$SampleName, "LMU-NU-\\d+(?=_)")
    otherdata = otherdata %>%
      mutate(patient.name = str_extract(SampleName, "LMU-NU-\\d+(?=_)"))
    otherdata$SampleName
    otherdata$patient.name

    sample.annotate = openxlsx::read.xlsx('data/other.data/LMU/Sample_Annotations_Subtypes.xlsx')
    sample.annotate %>% head()
    sample.annotate$Age.at.Collection %>% range()

    otherdata
    otherdata.wide = otherdata %>% group_by(patient.name, Target) %>%
      summarise(NPQ = mean(NPQ)) %>% ungroup() %>%
      pivot_wider(id_cols = patient.name, values_from = NPQ, names_from = Target) %>%
      column_to_rownames('patient.name')
    otherdata.wide

    colnames(otherdata.wide) <- gsub("[^[:alnum:]]", "", colnames(otherdata.wide))


    sex.annotate = left_join(data.frame(Sample.ID = rownames(otherdata.wide)), sample.annotate, by = 'Sample.ID') %>%
      select(Sample.ID, Sex) %>%
      mutate(Sex = ifelse(str_detect(Sex, 'm'), 1, 0) )
    sex.annotate

    otherdata.wide

    otherdata.wide$Sex = sex.annotate$Sex
    otherdata.wide = otherdata.wide %>% filter(!is.na(Sex)) %>%  as.matrix()
    ################################################################################

    library(pROC)

    # browser()
    train.comb.reformat = train.comb
    colnames(train.comb.reformat) <- gsub("[^[:alnum:]]", "", colnames(train.comb.reformat))
    # colnames(vali) <- gsub("[^[:alnum:]]", "", colnames(vali))

    # exclude.proteins = setdiff(union(colnames(train.comb.reformat), colnames(otherdata.wide)),
    #                            intersect(colnames(train.comb.reformat), colnames(otherdata.wide)))
    # exclude.proteins = grep('age', exclude.proteins, invert = T, value = T)
    # exclude.proteins
    common.proteins = intersect(colnames(train.comb.reformat), colnames(otherdata.wide))
    common.proteins

    train.scaled = train.comb.reformat %>% data.frame() %>% select(all_of(common.proteins)) %>%
      select(-Sex) %>% scale()
    train.comb.reformat = train.scaled %>% data.frame() %>%
      mutate(age = train.comb.reformat[, 'age'], Sex = train.comb.reformat[, 'Sex']) %>%
      as.matrix()

    scale.means <- attr(train.scaled, "scaled:center")
    scale.sds <- attr(train.scaled, "scaled:scale")

    otherdata.wide.tmp <- otherdata.wide %>% data.frame() %>% select(all_of(common.proteins)) %>%
      select(-Sex) %>%
      scale(., center = scale.means, scale = scale.sds)
    # otherdata.wide.tmp$Sex = otherdata.wide[, 'Sex', drop = T]
    otherdata.wide.tmp = otherdata.wide.tmp %>% data.frame() %>% mutate(Sex = otherdata.wide[, 'Sex', drop = T])
    otherdata.wide = otherdata.wide.tmp %>% as.matrix()
    otherdata.wide



    model.list = c('lm',  'lightGBM','SVM', 'elastic.net')
    # model.list = c('lm')
    for (model.type in model.list){
      print(model.type)
      ##################################################################################
      if (model.type == 'lm'){
        if (test.fold.i > 5){
          svm1 = read_csv(paste0('figures/lm.res/lm.feature.importance.csv'))
          svm1$protein[1:20]
          lm.res = lm(as.formula(paste0('age ~ ',
                                        # paste(svm1$protein[1:20], collapse = '+'))),
                                        paste(
                                          # svm1 %>% filter(adj.pval < 0.01) %>% pull(protein),
                                          c(svm1$protein[1:20], 'Sex'),
                                          # svm1$protein[1],
                                          collapse = '+'))),
                      data = train.comb.reformat %>% data.frame())

          predict(lm.res, otherdata.wide %>% data.frame())


          model.training <- lm(Age ~ pred.age,
                               data = data.frame(pred.age = predict(lm.res, train.comb.reformat[, !colnames(train.comb.reformat) %in% "age"] %>% data.frame()),
                                                 Age = train.comb.reformat[, 'age']))
          model.training
          pred.res <- predict(lm.res, otherdata.wide %>% data.frame())
          pred.res
          # pred_svm_untuned = pred_svm_untuned - 5
          # pred.res = pred.res[,1]
          new.pred = coefficients(model.training)['pred.age'] * pred.res + coefficients(model.training)[1]
          new.pred

        } else {next}

      }

      ##################################################################################
      if (model.type == 'elastic.net'){
        if (test.fold.i < 6) {
          svm1 = read_csv(paste0('figures/model.iter.12/', model.type, '/both/both.model.', test.fold.i, '.csv'))
          sorted.res = svm1 %>% arrange(rmse)
          j = 1
        } else {
          svm1 = read_csv(paste0('figures/model.iter.12/', model.type, '/both/aggregated.performance.rank.csv'))
          sorted.res = svm1
          j = 1
        }
        glm.manual = glmnet(
          x =  train.comb.reformat[, !colnames(train.comb.reformat) %in% "age"],
          y = train.comb.reformat[, 'age'],
          alpha = sorted.res$alpha[j],
          lambda = sorted.res$lambda[j]
        )
        model.training <- lm(Age ~ pred.age,
                             data = data.frame(pred.age = predict(glm.manual, train.comb.reformat[, !colnames(train.comb.reformat) %in% "age"])[,1],
                                               Age = train.comb.reformat[, 'age']))
        model.training
        pred.res <-predict(
          glm.manual,
          otherdata.wide,
        )
        # pred_svm_untuned = pred_svm_untuned - 5
        pred.res = pred.res[,1]
        new.pred = coefficients(model.training)['pred.age'] * pred.res + coefficients(model.training)[1]
        new.pred

      }

      ##################################################################################

      if (model.type == 'SVM'){
        if (test.fold.i < 6) {
          svm1 = read_csv(paste0('figures/model.iter.12/', model.type, '/both/both.model.', test.fold.i, '.csv'))
          sorted.res = svm1 %>% arrange(rmse)
          j = 1
        } else {
          svm1 = read_csv(paste0('figures/model.iter.12/', model.type, '/both/aggregated.performance.rank.csv'))
          sorted.res = svm1
          j = 1
        }
        m_svm_untuned <- svm(
          formula = age ~ .,
          data    = train.comb.reformat,
          gamma = sorted.res$gamma[j],
          cost = sorted.res$cost[j]
        )
        model.training <- lm(Age ~ pred.age,
                             data = data.frame(pred.age = predict(m_svm_untuned, train.comb.reformat),
                                               Age = train.comb.reformat[, 'age']))
        model.training

        # browser()
        pred.res <-predict(
          m_svm_untuned,
          newdata = otherdata.wide,
        )
        # pred_svm_untuned = pred_svm_untuned - 5
        pred.res
        new.pred = coefficients(model.training)['pred.age'] * pred.res + coefficients(model.training)[1]
        new.pred
      }

      ##################################################################################
      if (model.type == 'lightGBM'){
        if (test.fold.i < 6) {
          svm1 = read_csv(paste0('figures/model.iter.12/', model.type, '/both/both.model.', test.fold.i, '.csv'))
          sorted.res = svm1 %>% arrange(rmse)
          j = 1
        } else {
          svm1 = read_csv(paste0('figures/model.iter.12/', model.type, '/both/aggregated.performance.rank.csv'))
          sorted.res = svm1
          j = 1
        }
        # browser()
        train_lgb <- lightgbm::lgb.Dataset(train.comb.reformat[, !colnames(train.comb.reformat) %in% "age"],label=train.comb.reformat[, 'age'])
        vali_lgb <- lightgbm::lgb.Dataset.create.valid(train_lgb, train.comb.reformat[, !colnames(train.comb.reformat) %in% "age"],label=train.comb.reformat[, 'age'])
        # vali_lgb <- lightgbm::lgb.Dataset.create.valid(train_lgb,vali,label = y_vali)
        light_gbn_tuned <- lightgbm::lgb.train(
          params = list(
            objective = "regression",
            metric = "rmse",
            max_depth = sorted.res$max_depth[j],
            num_leaves =sorted.res$num_leaves[j],
            num_iterations = sorted.res$num_iterations[j],
            early_stopping_rounds=sorted.res$early_stopping_rounds[j],
            verbose = -1,
            learning_rate = sorted.res$learning_rate[j]
            #feature_fraction = .9
          ),
          valids = list(vali = vali_lgb),
          data = train_lgb
        )
        yhat_fit_tuned <- predict(light_gbn_tuned,train.comb.reformat[, !colnames(train.comb.reformat) %in% "age"])

        model.training <- lm(Age ~ pred.age,
                             data = data.frame(pred.age = yhat_fit_tuned,
                                               Age = train.comb.reformat[, 'age']))
        model.training
        # browser()
        pred.res <-predict(
          light_gbn_tuned,
          otherdata.wide ,
        )
        # pred_svm_untuned = pred_svm_untuned - 5
        pred.res
        new.pred = coefficients(model.training)['pred.age'] * pred.res + coefficients(model.training)[1]
        new.pred
      }

      ################################################################################

      others.age = data.frame(SampleName = rownames(otherdata.wide)) %>%
        left_join(sample.annotate %>% rename(age = Age.at.Collection) %>%
                    select(Sample.ID, age) %>% distinct(), by = c('SampleName' = 'Sample.ID'))
      others.age %>% head()
      others.age$age


      pred.rmse = RMSE(pred.res,  others.age$age, na.rm = T)
      adj.pred.rmse = RMSE(new.pred, others.age$age, na.rm = T)
      model <- lm(pred.res ~ others.age$age)
      model.new <- lm(new.pred ~ others.age$age)
      r2 = summary(model)$r.squared
      adj.r2 = summary(model.new)$r.squared
      pred.cor = cor(pred.res,  others.age$age, use = 'c')
      adj.pred.cor = cor(new.pred, others.age$age, use = 'c')

      # plot(pred.res, others.age$age)
      # abline(coef = c(0,1))

      plot.data = rbind(data.frame(pred.age = pred.res, chronological.age = others.age$age, type = 'org'),
                        data.frame(pred.age = new.pred, chronological.age = others.age$age, type = 'adjusted')
      )
      p = plot.data %>%
        ggplot(aes(chronological.age, pred.age, color = type)) +
        geom_point() +
        geom_smooth(method = 'lm') +
        geom_abline(slope = 1, intercept = 0) +
        annotate("text", x = -Inf, y = Inf, label = paste0('RMSE: ', round(pred.rmse, 3)),
                 hjust = -0.1, vjust = 2, size = 5, color = "blue") +
        annotate("text", x = -Inf, y = Inf, label = paste0('Adj. RMSE: ', round(adj.pred.rmse, 3)),
                 hjust = -0.1, vjust = 3.5, size = 5, color = "blue") +
        annotate("text", x = -Inf, y = Inf, label = paste0('R^2: ', round(r2, 3)),
                 hjust = -0.1, vjust = 5, size = 5, color = "orange") +
        annotate("text", x = -Inf, y = Inf, label = paste0('Adj. R^2: ', round(adj.r2, 3)),
                 hjust = -0.1, vjust = 6.5, size = 5, color = "orange") +
        annotate("text", x = -Inf, y = Inf, label = paste0('Corr: ', round(pred.cor, 3)),
                 hjust = -0.1, vjust = 8, size = 5, color = "green") +
        annotate("text", x = -Inf, y = Inf, label = paste0('Adj. Corr: ', round(adj.pred.cor, 3)),
                 hjust = -0.1, vjust = 9.5, size = 5, color = "green") +
        theme_bw(base_size = 15)



      # browser()
      lmu.pred.folder = paste0('figures/prediction/lmu.whole.training/', matrix.type, '/')
      dir.create(lmu.pred.folder, recursive = T)
      pdf(paste0(lmu.pred.folder, 'rmse.',  model.type, '.', test.fold.i, '.pdf'))
      print(p)
      dev.off()


      res.table[[counter]] = data.frame(model = model.type, test.fold = test.fold.i,
                                        rmse = pred.rmse, adj.rmse = adj.pred.rmse,
                                        r2 = r2, adj.r2 = adj.r2,
                                        corr = pred.cor, adj.corr = adj.pred.cor
      )
      counter = counter + 1
    }
  }
  res.table.final = res.table %>% data.table::rbindlist()
  res.table.final = res.table.final %>% mutate(test.fold = ifelse(test.fold == 6, 'best.param', test.fold))
  lmu.pred.folder = paste0('figures/prediction/lmu.whole.training/', matrix.type, '/')
  write.csv(res.table.final, paste0(lmu.pred.folder, 'res.table.csv'))

}




