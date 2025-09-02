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

for (capped.i in c('include.capped', 'remove.capped')){
  alpha.data.inflam <- read_csv("data/methuselah.alpha/inflammation/20240520_P-000412_Methuselah-Foundation_NULISAseq_TAP_Counts_Report_FULL.csv") %>%
    filter(!grepl("SC", SampleType)) %>%
    mutate(subject.id = str_extract(SampleName, "(?<=_)[^_]*$"))
  alpha.data.inflam$SampleName
  alpha.data.cns <- read_csv("data/methuselah.alpha/CNS/20240520_P-000412_Methuselah-Foundation_NULISAseq_TAP_Counts_Report_FULL.csv") %>%
    filter(!grepl("SC", SampleType)) %>%
    mutate(subject.id = str_extract(SampleName, "(?<=_)[^_]*$"))
  org.alpha.data = rbind(alpha.data.inflam, alpha.data.cns) %>% group_by(SampleName, subject.id, Target, SampleType) %>%
    summarise(NPQ = mean(NPQ)) %>% ungroup()


  metadata = read_csv('data/methuselah.alpha/sample_metadata.csv')
  metadata %>% group_by(Sex, Race) %>% summarise(n_distinct(SampleName))
  metadata

  full.data = org.alpha.data %>% left_join(metadata, by = c('SampleName', 'SampleType'))
  full.data
  # full.data = full.data %>% filter(Sex == sex.i)


  org.full.data.wide = full.data %>% select(SampleName, Sex, Age2, Target, NPQ) %>%
    pivot_wider(names_from = 'Target', values_from = 'NPQ')
  org.full.data.wide = org.full.data.wide %>%   mutate(subject.id = str_extract(SampleName, "(?<=_)[^_]*$")) %>% mutate(SampleName = subject.id) %>%
    select(-subject.id)

  # alpha.data

  argo.data.inflam <- openxlsx::read.xlsx("data/methuselah.argo/20250819_P-000412_Methuselah-Foundation_NULISAseq_TAP_Counts_Report_FULL.xlsx", sheet = 1) %>%
    filter(!grepl("SC", SampleType)) %>%
    mutate(SampleName = str_replace(SampleName, "_[^_]*$", "")) %>%
    mutate(subject.id = str_extract(SampleName, "(?<=_)[^_]*$"))
  argo.data.cns <- openxlsx::read.xlsx("data/methuselah.argo/20250819_P-000412_Methuselah-Foundation_NULISAseq_TAP_Counts_Report_FULL[27].xlsx", sheet = 1) %>%
    filter(!grepl("SC", SampleType)) %>%
    mutate(SampleName = str_replace(SampleName, "_[^_]*$", "")) %>%
    mutate(subject.id = str_extract(SampleName, "(?<=_)[^_]*$"))
  argo.data = rbind(argo.data.inflam, argo.data.cns) %>% group_by(SampleName, subject.id, Target, SampleType) %>%
    summarise(NPQ = mean(NPQ)) %>% ungroup()


  metadata = read_csv('data/methuselah.alpha/sample_metadata.csv')

  argo.full.data = argo.data %>% left_join(metadata %>%   mutate(subject.id = str_extract(SampleName, "(?<=_)[^_]*$")),
                                             by = c('subject.id'))

  argo.full.data.wide = argo.full.data %>% select(subject.id, Sex, Age2, Target, NPQ) %>%
    pivot_wider(names_from = 'Target', values_from = 'NPQ')


  # alpha.data = alpha.data %>%
  #   mutate(subject.id = str_extract(SampleName, "(?<=_)[^_]*$"))
  # argo.data = argo.data %>%
  #   mutate(subject.id = str_extract(SampleName, "(?<=_)[^_]*$"))
  # argo.data
  # alpha.data

  argo.data.inflam %>% head()
  subject.intersect = intersect(argo.data.inflam$subject.id, alpha.data.inflam$subject.id)
  # intersect(argo.data.cns$subject.id, alpha.data.cns$subject.id)
  inflam.target.intersect = intersect(argo.data.inflam$Target, alpha.data.inflam$Target)
  cns.target.intersect = intersect(argo.data.cns$Target, alpha.data.cns$Target)
  inflam.target.intersect
  cns.target.intersect
  # setdiff(argo.data.inflam$Target, alpha.data.inflam$Target)
  # setdiff( alpha.data.inflam$Target,argo.data.inflam$Target)
  # target.intersect

  alpha.inflam.wide = alpha.data.inflam %>% select(subject.id, Target, NPQ) %>% pivot_wider(names_from = 'subject.id', values_from = 'NPQ') %>%
    column_to_rownames('Target') %>% .[inflam.target.intersect, ]
  argo.inflam.wide = argo.data.inflam %>% select(subject.id, Target, NPQ) %>% pivot_wider(names_from = 'subject.id', values_from = 'NPQ') %>%
    column_to_rownames('Target') %>% .[inflam.target.intersect, ]
  alpha.inflam.wide %>% dim()
  argo.inflam.wide %>% dim()
  alpha.inflam.wide
  argo.inflam.wide
  alpha.cns.wide = alpha.data.cns %>% select(subject.id, Target, NPQ) %>% pivot_wider(names_from = 'subject.id', values_from = 'NPQ') %>%
    column_to_rownames('Target') %>% .[cns.target.intersect, ]
  argo.cns.wide = argo.data.cns %>% select(subject.id, Target, NPQ) %>% pivot_wider(names_from = 'subject.id', values_from = 'NPQ') %>%
    column_to_rownames('Target') %>% .[cns.target.intersect, ]

  #################################################################################
  #### bridging
  inflam.bridge_terms <- apply(argo.inflam.wide[, subject.intersect] - alpha.inflam.wide[, subject.intersect], 1, median)
  inflam.bridge_terms
  alpha.inflam.wide
  alpha.inflam.wide.bridged <- alpha.inflam.wide + inflam.bridge_terms
  alpha.inflam.bridged.long = alpha.inflam.wide.bridged %>%
    rownames_to_column('Target') %>%
    pivot_longer(cols = !Target, names_to = 'subject.id', values_to = 'NPQ')
  alpha.inflam.bridged.long

  cns.bridge_terms <- apply(argo.cns.wide[, subject.intersect] - alpha.cns.wide[, subject.intersect], 1, median)
  cns.bridge_terms
  alpha.cns.wide
  alpha.cns.wide.bridged <- alpha.cns.wide + cns.bridge_terms
  alpha.cns.bridged.long = alpha.cns.wide.bridged %>%
    rownames_to_column('Target') %>%
    pivot_longer(cols = !Target, names_to = 'subject.id', values_to = 'NPQ')

  bridged.alpha.data = rbind(alpha.inflam.bridged.long, alpha.cns.bridged.long) %>%
    group_by(subject.id, Target) %>%
    summarise(NPQ = mean(NPQ)) %>% ungroup()

  org.alpha.data
  bridged.full.data = bridged.alpha.data %>% left_join(metadata %>%
                                                         mutate(subject.id = str_extract(SampleName, "(?<=_)[^_]*$")),
                                                       by = c('subject.id'))

  bridged.alpha.full.data.wide = bridged.full.data %>% select(subject.id, Sex, Age2, Target, NPQ) %>%
    pivot_wider(names_from = 'Target', values_from = 'NPQ')
  bridged.alpha.full.data.wide

  ###################################################################################
  #####QC
  ###################################################################################
  browser()
  # both.methu.data = inner_join(alpha.data, argo.data, by = c('subject.id', 'Target'), suffix = c('.alpha', '.argo'))
  both.methu.data = inner_join(alpha.data.inflam, argo.data.inflam, by = c('subject.id', 'Target'), suffix = c('.alpha', '.argo'))
  both.methu.data = inner_join(alpha.data.cns, argo.data.cns, by = c('subject.id', 'Target'), suffix = c('.alpha', '.argo'))


  both.methu.data = inner_join(
    alpha.inflam.wide.bridged %>% rownames_to_column('Target') %>% pivot_longer(cols = !Target, names_to = 'subject.id', values_to = 'NPQ'),
    argo.inflam.wide %>% rownames_to_column('Target') %>% pivot_longer(cols = !Target, names_to = 'subject.id', values_to = 'NPQ'),
    by = c('subject.id', 'Target'), suffix = c('.alpha', '.argo'))
  both.methu.data = inner_join(
    alpha.cns.wide.bridged %>% rownames_to_column('Target') %>% pivot_longer(cols = !Target, names_to = 'subject.id', values_to = 'NPQ'),
    argo.cns.wide %>% rownames_to_column('Target') %>% pivot_longer(cols = !Target, names_to = 'subject.id', values_to = 'NPQ'),
    by = c('subject.id', 'Target'), suffix = c('.alpha', '.argo'))

  both.methu.data
  # argo.inflam.wide
  # both.methu.data = inner_join(bridged.full.data, argo.data.inflam, by = c('subject.id', 'Target'), suffix = c('.alpha', '.argo'))
  both.methu.data %>%
    ggplot(aes(NPQ.alpha, NPQ.argo, color = Target)) +
    geom_point() +
    geom_abline(slope = 1) +
    theme(legend.position = 'None')
  #

  combined.wide = rbind(
  alpha.inflam.wide %>% t() %>% data.frame(),
  argo.inflam.wide %>% t() %>% data.frame()
  )
  res.pca = prcomp(combined.wide, scale = T)
  inflam.groups = as.factor(c(rep('alpha',ncol(alpha.inflam.wide)), rep('argo', ncol(argo.inflam.wide))))
  groups <-
  fviz_pca_ind(res.pca,
               col.ind = inflam.groups, # color by groups
               palette = c("#00AFBB",  "#FC4E07"),
               addEllipses = TRUE, # Concentration ellipses
               ellipse.type = "confidence",
               legend.title = "Groups",
               label = 'none',
               repel = TRUE
  ) %>% print()

  combined.wide = rbind(
    alpha.inflam.wide.bridged %>% t() %>% data.frame(),
    argo.inflam.wide %>% t() %>% data.frame()
  )
  res.pca = prcomp(combined.wide, scale = T)
  inflam.groups = as.factor(c(rep('alpha',ncol(alpha.inflam.wide)), rep('argo', ncol(argo.inflam.wide))))
  groups <-
    fviz_pca_ind(res.pca,
                 col.ind = inflam.groups, # color by groups
                 palette = c("#00AFBB",  "#FC4E07"),
                 addEllipses = TRUE, # Concentration ellipses
                 ellipse.type = "confidence",
                 legend.title = "Groups",
                 label = 'none',
                 repel = TRUE
    ) %>% print()


  ################################################################################
  ################################################################################
  ################################################################################
  ################################################################################

  for (harmonized.i in c('bridged', 'not.bridged')){
  # for (harmonized.i in c('not.bridged')){

    browser()
    if (harmonized.i == 'bridged'){
      adjusted.wide = bridged.alpha.full.data.wide %>% rename(SampleName = subject.id)
    } else if (harmonized.i == 'not.bridged'){
      adjusted.wide = org.full.data.wide
    }
    colnames(adjusted.wide) <- gsub("[^[:alnum:]]", "", colnames(adjusted.wide))
    ### remove bridging samples from training
    adjusted.wide = adjusted.wide %>% anti_join(argo.full.data.wide %>% select(subject.id), by = c('SampleName' = 'subject.id'))
    ### remove ones that are 90 (capped)
    if (capped.i == 'remove.capped'){
      adjusted.wide = adjusted.wide %>% filter(Age2 != 90)
    }

    # adjusted.wide = full.data.wide
    x.variables = adjusted.wide %>% select(-c(SampleName, Age2, Sex))
    logical.complete.cases = complete.cases(x.variables)
    x.variables = x.variables[logical.complete.cases,]
    y.response = adjusted.wide$Age2
    y.response = y.response[logical.complete.cases]
    meta.data = adjusted.wide %>% select(SampleName, Age2, Sex)


    test.fold <- caret::createFolds(y.response, k = 5, list = TRUE, returnTrain = FALSE)
    for (test.fold.i in 1:5){
      test.set.x = x.variables[test.fold[[test.fold.i]],]
      test.set.y = y.response[test.fold[[test.fold.i]]]
      remains.x.variables = x.variables[-test.fold[[test.fold.i]],]
      remains.y.response = y.response[-test.fold[[test.fold.i]]]
      test.set.meta = meta.data[test.fold[[test.fold.i]],]



      test = test.set.x %>% as.matrix()
      y_test = test.set.y

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

      # model.type = 'SVM'
      # dir.create(paste0('figures/model.iter.', seed.i, '/',  model.type, '/', sex.i, '/'), recursive = T)

      # train <- train.set.x %>% as.matrix()
      train <- x.variables %>% as.matrix()

      vali = vali.set.x %>% as.matrix()
      test <- test.set.x %>% as.matrix()

      # y_train <- train.set.y
      y_train <- y.response

      y_vali = vali.set.y
      y_test <- test.set.y
      #
      train.comb = cbind(train, age = y_train)




      ################################################################################
      # other1.cns = read.csv('data/other.data/Basel/20240724_P-000507_University-of-Basel_NULISAseq_TAP_Counts_Report_FULL.csv')
      # other1.cns$SampleName %>% unique()
      # other1.ifn = read.csv('data/other.data/Basel/20240725_P-000507_University-of-Basel_NULISAseq_TAP_Counts_Report_FULL.csv')
      # other1.cns %>% head()
      # other1.ifn %>% head()
      # # other1.cns %>% filter(str_detect(SampleName, 'D_03_4031-6531_2017-07-04')) %>% select(SampleName) %>% distinct()
      # # other1.cns %>% filter(str_detect(SampleName, 'D_09')) %>% select(SampleName) %>% distinct()
      # other1.ifn %>% head()
      # otherdata = rbind(
      #   other1.cns %>% select(SampleName, SampleType, Target, NPQ),
      #   other1.ifn %>% select(SampleName, SampleType, Target, NPQ)
      # ) %>% group_by(SampleName, Target, SampleType) %>%
      #   summarise(NPQ = mean(NPQ)) %>% ungroup()
      # otherdata$SampleName
      # # str_extract(otherdata$SampleName, "LMU-NU-\\d+(?=_)")
      # # otherdata = otherdata %>%
      # #   mutate(patient.name = str_extract(SampleName, "LMU-NU-\\d+(?=_)"))
      #
      # otherdata = otherdata %>% select(SampleName) %>% distinct() %>% left_join(normalized.long) %>%
      #   mutate(SampleType = 'Sample')
      #
      #
      # otherdata.wide = otherdata %>% pivot_wider(id_cols = SampleName, values_from = 'NPQ', names_from = 'Target') %>%
      #   column_to_rownames('SampleName')
      # browser()

      ### remove capped at 90 samples
      if (capped.i == 'remove.capped'){
        argo.full.data.wide.alt =  argo.full.data.wide %>% filter(Age2 != 90)
      } else if (capped.i == 'include.capped'){
        argo.full.data.wide.alt =  argo.full.data.wide
      }
      sample.annotate = argo.full.data.wide.alt %>% select(subject.id, Age2) %>%
        rename(SampleName = subject.id , age = Age2)
      otherdata.wide = argo.full.data.wide.alt
      colnames(otherdata.wide) <- gsub("[^[:alnum:]]", "", colnames(otherdata.wide))
      otherdata.wide = otherdata.wide %>% select(-c(Sex, Age2)) %>%
        column_to_rownames('subjectid')

      otherdata.wide

      # sample.annotate = read.csv('data/other.data/Basel/Pira_manifest.csv')
      ################################################################################

      library(pROC)

      # browser()
      train.comb.reformat = train.comb
      colnames(train.comb.reformat) <- gsub("[^[:alnum:]]", "", colnames(train.comb.reformat))

      protein.list = intersect(colnames(otherdata.wide), grep('age', colnames(train.comb.reformat), invert = T, value = T))
      # otherdata.wide = otherdata.wide %>% select(grep('age', colnames(train.comb.reformat), invert = T, value = T))
      otherdata.wide = otherdata.wide %>% select(protein.list)
      train.comb.reformat = train.comb.reformat[, c(protein.list, 'age')]
      vali = vali[,protein.list]



      model.list = c(  'lightGBM','SVM', 'elastic.net')
      for (model.type in model.list){
        print(model.type)

        ##################################################################################
        if (model.type == 'elastic.net'){
          svm1 = read_csv(paste0('figures/model.iter.12/elastic.net/both/both.model.', test.fold.i, '.csv'))
          sorted.res = svm1 %>% arrange(rmse)
          j = 1
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
            # otherdata.wide %>% select(-DDC) %>% as.matrix(),
            otherdata.wide %>% as.matrix()
          )
          # pred_svm_untuned = pred_svm_untuned - 5
          pred.res = pred.res[,1]
          new.pred = coefficients(model.training)['pred.age'] * pred.res + coefficients(model.training)[1]
          new.pred

        }

        ##################################################################################

        if (model.type == 'SVM'){
          svm1 = read_csv(paste0('figures/model.iter.12/SVM/both/both.model.', test.fold.i, '.csv'))
          sorted.res = svm1 %>% arrange(rmse)
          j = 1
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
          svm1 = read_csv(paste0('figures/model.iter.12/lightGBM/both/both.model.', test.fold.i, '.csv'))
          sorted.res = svm1 %>% arrange(rmse)
          j = 1
          # browser()
          train_lgb <- lightgbm::lgb.Dataset(train.comb.reformat[, !colnames(train.comb.reformat) %in% "age"],label=train.comb.reformat[, 'age'])
          vali_lgb <- lightgbm::lgb.Dataset.create.valid(train_lgb,vali,label = y_vali)
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
          pred.res <-predict(
            light_gbn_tuned,
            # otherdata.wide %>% select(-DDC) %>% as.matrix(),
            otherdata.wide %>% as.matrix()
          )
          # pred_svm_untuned = pred_svm_untuned - 5
          pred.res
          new.pred = coefficients(model.training)['pred.age'] * pred.res + coefficients(model.training)[1]
          new.pred
        }

        ################################################################################

        others.age = data.frame(SampleName = rownames(otherdata.wide)) %>%
          left_join(sample.annotate %>% select(SampleName, age) %>% distinct(), by = 'SampleName')
        others.age %>% head()
        others.age$age

        pred.rmse = RMSE(pred.res,  others.age$age, na.rm = T)
        adj.pred.rmse = RMSE(new.pred, others.age$age, na.rm = T)
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
          theme_bw(base_size = 15)


        basel.pred.folder = paste0('figures/prediction/methuselah.argo/', harmonized.i, '/', capped.i, '/')
        dir.create(basel.pred.folder, recursive = T)
        pdf(paste0(basel.pred.folder, 'rmse.',  model.type, '.', test.fold.i, '.pdf'))
        print(p)
        dev.off()
      }
    }

  }

}
