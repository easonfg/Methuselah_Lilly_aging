 library(tidyverse)
 library(caret)
dein.info <- openxlsx::read.xlsx("~/Downloads/OneDrive_1_8-25-2025/20250620-1127_Bay1_CNS Disease Panel 120 _RBD_Donor Information for Alamar.xlsx", sheet = 2)
dein.info
dein.info = dein.info %>% select(Deidentified, Age) %>% distinct() %>% filter(!str_detect(Deidentified, 'SC'))
dein.info %>% head()


list.files('figures/prediction/lilly/bridged/include.capped/', full.names = T)
all.res = lapply(list.files('figures/prediction/lilly/bridged/include.capped/', full.names = T), function(model.i){
# all.res = lapply(list.files('figures/prediction/lilly/not.bridged//include.capped/', full.names = T), function(model.i){
# all.res = lapply(list.files('figures/prediction/lilly/bridged/remove.capped//', full.names = T), function(model.i){
  print(model.i %>% basename())
  model1.res = read.csv(model.i)
  model1.res %>% head()
  comb.res = dein.info %>% inner_join(model1.res, by = c('Deidentified' = 'X'))
  comb.res %>% head()
  data.frame(model = basename(model.i), RMSE(comb.res$pred.age, comb.res$Age))

  # comb.res %>%
  #   ggplot(aes(Age, pred.age)) +
  #     geom_smooth(method = 'lm') +
  #     geom_abline(slope = 1, intercept = 0) +
  #   geom_point()
})
all.res %>% data.table::rbindlist() %>% print()
all.res %>% unlist() %>% mean()


en1 = read.csv('figures/model.iter.12/elastic.net/both/both.model.1.csv')
svm1 = read.csv('figures/model.iter.12/SVM/both/both.model.1.csv')
lgbm1 = read.csv('figures/model.iter.12/lightGBM/both/both.model.1.csv')
en1 %>% arrange(desc(auc)) %>% head()
svm1 %>% arrange(desc(auc)) %>% head()
lgbm1 %>% arrange(desc(auc)) %>% head()

en1 %>% arrange(rmse) %>% head()
svm1 %>% arrange(rmse) %>% head()
lgbm1 %>% arrange(rmse) %>% head()

en1 %>% arrange(error) %>% head()
svm1 %>% arrange(error) %>% head()
lgbm1 %>% arrange(error) %>% head()
