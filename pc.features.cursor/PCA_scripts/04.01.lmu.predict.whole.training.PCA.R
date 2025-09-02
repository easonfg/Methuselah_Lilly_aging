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
  save(pca_result, n_pcs, variance_explained, file = paste0('PCA_results/PCA_results_LMU_prediction_', output.var, '.RData'))
  
  # Use PCA data for modeling
  full.data.wide <- full.data.wide.pca

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
    # y_vali <- vali.set.y
    # y_test <- test.set.y
    #
    train.comb = cbind(train, age = y_train)

    ################################################################################
    other1.cns = read.csv('data/other.data/LMU/20240906_P-000518-CNS_LMU_NULISAseq_TAP_Counts_Report_FULL.csv')
    other1.ifn = read.csv('data/other.data/LMU/20240906_P-000580-IF_LMU_NULISAseq_TAP_Counts_Report_FULL.csv')
    other1.cns %>% head()
    other1.ifn %>% head()

    # matrix.type = 'Serum'
    matrix.type = 'CSF'
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
    ### APPLY PCA TRANSFORMATION TO LMU DATA ###
    
    # Load the PCA results from Methuselah data to get the transformation matrix
    load(paste0('PCA_results/PCA_results_LMU_prediction_', output.var, '.RData'))
    
    # Get the protein names that were used for PCA (from the rotation matrix)
    methuselah_proteins <- rownames(pca_result$rotation)
    cat("Number of proteins used in Methuselah PCA:", length(methuselah_proteins), "\n")
    
    # Get the protein columns from LMU data (excluding Sex)
    protein_cols_lmu <- colnames(otherdata.wide)[!colnames(otherdata.wide) %in% c("Sex")]
    cat("Number of proteins in LMU data:", length(protein_cols_lmu), "\n")
    
    # Clean protein names by removing symbols and special characters for better matching
    methuselah_proteins_clean <- gsub("[^[:alnum:]]", "", methuselah_proteins)
    protein_cols_lmu_clean <- gsub("[^[:alnum:]]", "", protein_cols_lmu)
    
    # Create mapping between clean and original names
    methuselah_clean_to_original <- setNames(methuselah_proteins, methuselah_proteins_clean)
    lmu_clean_to_original <- setNames(protein_cols_lmu, protein_cols_lmu_clean)
    
    # Find common proteins using cleaned names
    common_proteins_clean <- intersect(methuselah_proteins_clean, protein_cols_lmu_clean)
    cat("Number of common proteins (after cleaning names):", length(common_proteins_clean), "\n")
    
    # Map back to original protein names
    common_proteins_methuselah <- methuselah_clean_to_original[common_proteins_clean]
    common_proteins_lmu <- lmu_clean_to_original[common_proteins_clean]
    cat("Number of common proteins (original names):", length(common_proteins_methuselah), "\n")
    
    if(length(common_proteins_methuselah) > 0) {
      cat("Common proteins found. Applying proper PCA transformation.\n")
      
      # Extract the common proteins from LMU data using LMU protein names
      lmu_protein_data <- otherdata.wide[, common_proteins_lmu, drop = FALSE]
      
      # Get the PCA loadings for common proteins
      pca_loadings_common <- pca_result$rotation[common_proteins_methuselah, 1:n_pcs, drop = FALSE]
      
      # IMPORTANT: We need to center and scale the LMU data using the SAME parameters as Methuselah
      # The PCA object stores the centering and scaling parameters
      # For new data, we need to apply: (X_new - center) / scale
      
      # Get the centering and scaling parameters from the original PCA
      # Note: prcomp() with scale.=TRUE and center=TRUE stores these in the PCA object
      # But we need to extract them from the original protein data that was used for PCA
      
      # For now, let's use the predict() function which handles this automatically
      # This is the proper way to apply PCA transformation to new data
      
      # First, we need to ensure the LMU data has the same structure as the original
      # Create a temporary dataset with the same proteins as Methuselah (fill missing with 0)
      lmu_proteins_complete <- matrix(0, nrow = nrow(lmu_protein_data), ncol = length(methuselah_proteins))
      colnames(lmu_proteins_complete) <- methuselah_proteins
      
      # Fill in the common proteins
      lmu_proteins_complete[, common_proteins_methuselah] <- as.matrix(lmu_protein_data)
      
      # Now apply the PCA transformation using predict()
      lmu_pc_scores <- predict(pca_result, newdata = lmu_proteins_complete)
      
      # Take only the first n_pcs components
      lmu_pc_scores <- lmu_pc_scores[, 1:n_pcs, drop = FALSE]
      colnames(lmu_pc_scores) <- paste0("PC", 1:n_pcs)
      
      # Create new LMU dataset with PCs
      lmu_pc_data <- data.frame(
        SampleName = rownames(otherdata.wide),
        lmu_pc_scores,
        Sex = otherdata.wide[, "Sex"]
      )
      
      # Convert to matrix format for consistency
      otherdata.wide <- as.matrix(lmu_pc_data[, -1]) # Remove SampleName, keep PCs and Sex
      rownames(otherdata.wide) <- lmu_pc_data$SampleName
      
      cat("Successfully transformed LMU data to", n_pcs, "PCs using proper PCA transformation\n")
      
    } else {
      cat("No common proteins found. Cannot apply PCA transformation.\n")
      next
    }

    ################################################################################

    library(pROC)

    # browser()
    train.comb.reformat = train.comb
    colnames(train.comb.reformat) <- gsub("[^[:alnum:]]", "", colnames(train.comb.reformat))
    # colnames(vali) <- gsub("[^[:alnum:]]", "", colnames(train.comb.reformat))
    
    # IMPORTANT: Since we're now using PCA-transformed data, we don't need additional scaling
    # The PCA transformation already handles the scaling from the original protein data
    # We just need to ensure the column names match and handle the Sex column
    
    # Get common columns between training data (PCs) and LMU data (PCs)
    common.proteins = intersect(colnames(train.comb.reformat), colnames(otherdata.wide))
    cat("Common columns between training and LMU data:", length(common.proteins), "\n")
    
    # Select only the common columns (PCs) from training data, but keep the age column
    train.comb.reformat = train.comb.reformat %>% data.frame() %>% 
      select(all_of(common.proteins), age) %>%
      as.matrix()
    
    # Select only the common columns (PCs) from LMU data
    otherdata.wide = otherdata.wide %>% data.frame() %>% 
      select(all_of(common.proteins)) %>%
      as.matrix()
    
    cat("Training data shape:", dim(train.comb.reformat), "\n")
    cat("LMU data shape:", dim(otherdata.wide), "\n")

    model.list = c('lm',  'lightGBM','SVM', 'elastic.net')
    # model.list = c('lm')
    for (model.type in model.list){
      print(model.type)
      ##################################################################################
      if (model.type == 'lm'){
        if (test.fold.i > 5){
          # For PCA-based models, we'll use the top PC features instead of protein features
          # Since we have 165 PCs, let's use the first 20 PCs plus Sex
          pc_features <- paste0("PC", 1:20)
          lm.res = lm(as.formula(paste0('age ~ ',
                                        paste(c(pc_features, 'Sex'),
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
          svm1 = read_csv(paste0('PCA_results/model.iter.12/elastic.net.PCA/both/both.model.', test.fold.i, '.csv'))
          sorted.res = svm1 %>% arrange(rmse)
          j = 1
        } else {
          svm1 = read_csv(paste0('PCA_results/model.iter.12/elastic.net.PCA/both/aggregated.performance.rank.csv'))
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
          svm1 = read_csv(paste0('PCA_results/model.iter.12/SVM.PCA/both/both.model.', test.fold.i, '.csv'))
          sorted.res = svm1 %>% arrange(rmse)
          j = 1
        } else {
          svm1 = read_csv(paste0('PCA_results/model.iter.12/SVM.PCA/both/aggregated.performance.rank.csv'))
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
          svm1 = read_csv(paste0('PCA_results/model.iter.12/lightgbm.PCA/both/both.model.', test.fold.i, '.csv'))
          sorted.res = svm1 %>% arrange(rmse)
          j = 1
        } else {
          svm1 = read_csv(paste0('PCA_results/model.iter.12/lightgbm.PCA/both/aggregated.performance.rank.csv'))
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
      lmu.pred.folder = paste0('PCA_results/prediction/lmu.whole.training.PCA/', matrix.type, '/')
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
  lmu.pred.folder = paste0('PCA_results/prediction/lmu.whole.training.PCA/', matrix.type, '/')
  write.csv(res.table.final, paste0(lmu.pred.folder, 'res.table.csv'))

}

cat("PCA-based LMU prediction completed!\n")
cat("Number of PCs used:", n_pcs, "\n")
cat("Variance explained:", round(variance_explained[n_pcs] * 100, 2), "%\n")
