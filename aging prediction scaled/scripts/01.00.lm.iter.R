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

  full.data.wide
  ### scaling data
  full.data.wide <- full.data.wide %>%
    mutate(across(
      # Select columns to scale: everything EXCEPT Name, Sex, Age
      .cols = -c(SampleName, Sex, Age2),
      # Apply the scale function to each selected column
      .fns = scale
    ))

  # View the result
  head(full.data.wide)
  colnames(full.data.wide)
  colnames(full.data.wide) <- gsub("[^[:alnum:]]", "", colnames(full.data.wide))


  plot.full.data.wide = full.data.wide %>% select(-c(SampleName, Sex)) %>% scale()
  all.lm.res = lapply(full.data.wide %>% select(-c(SampleName, Sex, Age2)) %>% colnames(), function(prot.i){
    # browser()
    # lm(full.data.wide[, 'Age2'], full.data.wide[, prot.i])
    lm.res = lm(as.formula(paste0('Age2 ~ ', prot.i)), data = full.data.wide)
    summary(lm.res)$coefficients[2,] %>% data.frame() %>% t() %>% data.frame() %>%
      mutate(protein = prot.i)
    # lm(Age2 ~ ACHE, data = full.data.wide)

    # # plot.full.data.wide %>% data.frame()
    # plot.full.data.wide %>%
    #   ggplot() +
    #     geom_point(aes(Age2, !!sym(prot.i))) +
    #     geom_smooth(method = 'lm', aes(Age2, !!sym(prot.i)))

  })
  all.lm.res = all.lm.res %>% data.table::rbindlist()
  all.lm.res = all.lm.res %>% mutate(adj.pval = p.adjust(Pr...t..))
  all.lm.res %>% filter(protein == 'ACHE')
  all.lm.res %>% arrange(adj.pval) %>% view
  dir.create('figures/lm.res/')
  all.lm.res %>% arrange(adj.pval) %>%
    write.csv(., 'figures/lm.res/lm.feature.importance.csv')
}

# lm(Age2 ~ ACHE + CCL11, data = full.data.wide)
