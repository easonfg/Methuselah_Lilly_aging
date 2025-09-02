library(tidyverse)
library(TAPGeneralization)
library(doParallel); library(foreach)
library(caret)
library(openxlsx)
library(glmnet)
library(factoextra)

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
  org.full.data.wide$Sex = ifelse(org.full.data.wide$Sex == 'M', 1, 0)

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
  alpha.cns.wide.bridged.long = alpha.cns.wide.bridged %>%
    rownames_to_column('Target') %>%
    pivot_longer(cols = !Target, names_to = 'subject.id', values_to = 'NPQ')

  bridged.alpha.data = rbind(alpha.inflam.bridged.long, alpha.cns.wide.bridged.long) %>%
    group_by(subject.id, Target) %>%
    summarise(NPQ = mean(NPQ)) %>% ungroup()

  org.alpha.data
  bridged.full.data = bridged.alpha.data %>% left_join(metadata %>%
                                                         mutate(subject.id = str_extract(SampleName, "(?<=_)[^_]*$")),
                                                       by = c('subject.id'))

  bridged.alpha.full.data.wide = bridged.full.data %>% select(subject.id, Sex, Age2, Target, NPQ) %>%
    pivot_wider(names_from = 'Target', values_from = 'NPQ')
  bridged.alpha.full.data.wide
  bridged.alpha.full.data.wide$Sex = ifelse(bridged.alpha.full.data.wide$Sex == 'M', 1, 0)

  ###################################################################################
  #####QC
  ###################################################################################
  # browser()
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

  # browser()
  combined.wide = rbind(
  # alpha.inflam.wide %>% t() %>% data.frame(),
  # argo.inflam.wide %>% t() %>% data.frame()
  alpha.cns.wide %>% t() %>% data.frame(),
  argo.cns.wide %>% t() %>% data.frame()
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
    # alpha.inflam.wide.bridged %>% t() %>% data.frame(),
    # argo.inflam.wide %>% t() %>% data.frame()
    alpha.cns.wide.bridged %>% t() %>% data.frame(),
    argo.cns.wide %>% t() %>% data.frame()
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

    # browser()
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

    ################################################################################
    ### PCA TRANSFORMATION ###
    
    # Extract protein data for PCA (excluding metadata columns)
    protein_data <- adjusted.wide %>% 
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
    adjusted.wide.pca <- data.frame(
      SampleName = adjusted.wide$SampleName,
      Sex = adjusted.wide$Sex,
      Age2 = adjusted.wide$Age2,
      pc_scores
    )
    
    # Save PCA results for later analysis
    dir.create('PCA_results', recursive = TRUE, showWarnings = FALSE)
    save(pca_result, n_pcs, variance_explained, file = paste0('PCA_results/PCA_results_Argo_prediction_', harmonized.i, '_', capped.i, '.RData'))
    
    # Use PCA data for modeling
    adjusted.wide <- adjusted.wide.pca

    # adjusted.wide = full.data.wide
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
      # remains.y.response = y.response[-test.fold[[test.fold.i]],]
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
      # train <- train.set.x %>% as.matrix()
      # vali <- vali.set.x %>% as.matrix()
      #
      # y_train <- train.set.y
      # y_vali <- vali.set.y

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
      otherdata.wide = otherdata.wide %>% select(-c( Age2)) %>%
        column_to_rownames('subjectid')

      otherdata.wide
      otherdata.wide$Sex = ifelse(otherdata.wide$Sex == 'M', 1, 0)

      ################################################################################
      ### APPLY PCA TRANSFORMATION TO ARGO DATA ###
      
      # Load the PCA results from Methuselah data to get the transformation matrix
      load(paste0('PCA_results/PCA_results_Argo_prediction_', harmonized.i, '_', capped.i, '.RData'))
      
          # Get the protein names that were used for PCA (from the rotation matrix)
    methuselah_proteins <- rownames(pca_result$rotation)
    cat("Number of proteins used in Methuselah PCA:", length(methuselah_proteins), "\n")
    
    # Get the protein columns from Argo data (excluding Sex)
    protein_cols_argo <- colnames(otherdata.wide)[!colnames(otherdata.wide) %in% c("Sex")]
    cat("Number of proteins in Argo data:", length(protein_cols_argo), "\n")
    
    # Clean protein names by removing symbols and special characters for better matching
    methuselah_proteins_clean <- gsub("[^[:alnum:]]", "", methuselah_proteins)
    protein_cols_argo_clean <- gsub("[^[:alnum:]]", "", protein_cols_argo)
    
    # Create mapping between clean and original names
    methuselah_clean_to_original <- setNames(methuselah_proteins, methuselah_proteins_clean)
    argo_clean_to_original <- setNames(protein_cols_argo, protein_cols_argo_clean)
    
    # Find common proteins using cleaned names
    common_proteins_clean <- intersect(methuselah_proteins_clean, protein_cols_argo_clean)
    cat("Number of common proteins (after cleaning names):", length(common_proteins_clean), "\n")
    
        # Map back to original protein names
    common_proteins_methuselah <- methuselah_clean_to_original[common_proteins_clean]
    common_proteins_argo <- argo_clean_to_original[common_proteins_clean]
    cat("Number of common proteins (original names):", length(common_proteins_methuselah), "\n")
    
    if(length(common_proteins_methuselah) > 0) {
      cat("Common proteins found. Applying proper PCA transformation.\n")
      
      # Extract the common proteins from Argo data using Argo protein names
      argo_protein_data <- otherdata.wide[, common_proteins_argo, drop = FALSE]
        
              # Get the PCA loadings for common proteins
      pca_loadings_common <- pca_result$rotation[common_proteins_methuselah, 1:n_pcs, drop = FALSE]
      
      # IMPORTANT: We need to center and scale the Argo data using the SAME parameters as Methuselah
      # The PCA object stores the centering and scaling parameters
      # For new data, we need to apply: (X_new - center) / scale
      # For now, let's use the predict() function which handles this automatically
      # This is the proper way to apply PCA transformation to new data
      
      # First, we need to ensure the Argo data has the same structure as the original
      # Create a temporary dataset with the same proteins as Methuselah (fill missing with 0)
      argo_proteins_complete <- matrix(0, nrow = nrow(argo_protein_data), ncol = length(methuselah_proteins))
      colnames(argo_proteins_complete) <- methuselah_proteins
      
      # Fill in the common proteins
      argo_proteins_complete[, common_proteins_methuselah] <- as.matrix(argo_protein_data)
        
        # Now apply the PCA transformation using predict()
        argo_pc_scores <- predict(pca_result, newdata = argo_proteins_complete)
        
        # Take only the first n_pcs components
        argo_pc_scores <- argo_pc_scores[, 1:n_pcs, drop = FALSE]
        colnames(argo_pc_scores) <- paste0("PC", 1:n_pcs)
        
        # Create new Argo dataset with PCs
        argo_pc_data <- data.frame(
          SampleName = rownames(otherdata.wide),
          argo_pc_scores,
          Sex = otherdata.wide[, "Sex"]
        )
        
        # Convert to matrix format for consistency
        otherdata.wide <- as.matrix(argo_pc_data[, -1]) # Remove SampleName, keep PCs and Sex
        rownames(otherdata.wide) <- argo_pc_data$SampleName
        
        cat("Successfully transformed Argo data to", n_pcs, "PCs using proper PCA transformation\n")
        
      } else {
        cat("No common proteins found. Cannot apply PCA transformation.\n")
        next
      }

      # sample.annotate = read.csv('data/other.data/Basel/Pira_manifest.csv')
      ################################################################################

      library(pROC)

      # # browser()
      # train.comb.reformat = train.comb
      # colnames(train.comb.reformat) <- gsub("[^[:alnum:]]", "", colnames(train.comb.reformat))
      #
      # protein.list = intersect(colnames(otherdata.wide), grep('age', colnames(train.comb.reformat), invert = T, value = T))
      # # otherdata.wide = otherdata.wide %>% select(grep('age', colnames(train.comb.reformat), invert = T, value = T))
      # otherdata.wide = otherdata.wide %>% select(protein.list)
      # train.comb.reformat = train.comb.reformat[, c(protein.list, 'age')]
      # # vali = vali[,protein.list]

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

      # IMPORTANT: Since we're now using PCA-transformed data, we don't need additional scaling
      # The PCA transformation already handles the scaling from the original protein data
      # We just need to ensure the column names match and handle the Sex column
      
      # Get common columns between training data (PCs) and Argo data (PCs)
      common.proteins = intersect(colnames(train.comb.reformat), colnames(otherdata.wide))
      cat("Common columns between training and Argo data:", length(common.proteins), "\n")
      
      # Select only the common columns (PCs) from training data, but keep the age column
      train.comb.reformat = train.comb.reformat %>% data.frame() %>% 
        select(all_of(common.proteins), age) %>%
        as.matrix()
      
      # Select only the common columns (PCs) from Argo data
      otherdata.wide = otherdata.wide %>% data.frame() %>% 
        select(all_of(common.proteins)) %>%
        as.matrix()
      
      cat("Training data shape:", dim(train.comb.reformat), "\n")
      cat("Argo data shape:", dim(otherdata.wide), "\n")

      model.list = c('lm',  'lightGBM','SVM', 'elastic.net')
      for (model.type in model.list){
        print(model.type)
        ##################################################################################
        if (model.type == 'lm'){
          if (test.fold.i > 5){
            # For PCA-based models, we'll use the top PC features instead of protein features
            # Since we have n_pcs PCs, let's use the first 20 PCs plus Sex
            pc_features <- paste0("PC", 1:min(20, n_pcs))
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
          # vali_lgb <- lightgbm::lgb.Dataset.create.valid(train_lgb,vali,label = y_vali)
          vali_lgb <- lightgbm::lgb.Dataset.create.valid(train_lgb, train.comb.reformat[, !colnames(train.comb.reformat) %in% "age"],label=train.comb.reformat[, 'age'])
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
        argo.pred.folder = paste0('PCA_results/prediction/methuselah.argo.PCA/', harmonized.i, '/', capped.i, '/')
        dir.create(argo.pred.folder, recursive = T)
        pdf(paste0(argo.pred.folder, 'rmse.',  model.type, '.', test.fold.i, '.pdf'))
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
    argo.pred.folder = paste0('PCA_results/prediction/methuselah.argo.PCA/', harmonized.i, '/', capped.i, '/')
    write.csv(res.table.final, paste0(argo.pred.folder, 'res.table.csv'))
  }

}

cat("PCA-based Methuselah Argo prediction completed!\n")
cat("Number of PCs used:", n_pcs, "\n")
cat("Variance explained:", round(variance_explained[n_pcs] * 100, 2), "%\n")
