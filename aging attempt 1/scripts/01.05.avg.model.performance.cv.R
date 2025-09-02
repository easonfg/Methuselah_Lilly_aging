### average the model performance across the 5 fold cv to pick the best hyperparameter


# detach(package:plyr)
# library(dplyr)
library(tidyverse)

# model.i = 'elastic.net'
# model.i = 'lightGBM'
for (model.i in c('elastic.net', 'SVM', 'lightGBM')){
  model.i

  if (model.i == 'elastic.net'){
    iter.i = 'a'
  } else if (model.i == 'SVM') {
    iter.i = 'j'
  } else if (model.i == 'lightGBM') {
    iter.i = 'j'
  }

  model.performance.folder = paste0('figures/model.iter.12/', model.i, '/both/')

  all.res = lapply(list.files(model.performance.folder,
                              pattern = 'model.*.csv', full.names = T), function(model.i){
    # all.res = lapply(list.files('figures/prediction/lilly/not.bridged//include.capped/', full.names = T), function(model.i){
    # all.res = lapply(list.files('figures/prediction/lilly/bridged/remove.capped//', full.names = T), function(model.i){
    print(model.i %>% basename())
    model1.res = read.csv(model.i)
    basename(model.i)
    iter.i = stringr::str_extract_all(basename(model.i), "\\b\\d+\\b")
    model1.res %>% mutate(iter =  iter.i)
  })
  all.res = all.res %>% data.table::rbindlist()
  all.res
  all.res$lambda %>% unique()

  # all.res %>% group_by(alpha, lambda) %>%
  #   summarise(mean.rmse = mean(rmse)) %>% arrange(mean.rmse)

  all.res %>%
    group_by(iter) %>%
    mutate(rank = rank(rmse)) %>% ungroup() %>%
    # group_by(alpha, lambda) %>%
    group_by(across(-c(X, !!sym(iter.i), error, rmse, cor, auc, iter, rank))) %>%
    summarise(mean.rmse = mean(rmse), mean.rank = mean(rank), list.rank = paste(rank, collapse = ', ')) %>%
    ungroup() %>%
    mutate(rmse.rank = rank(mean.rmse), ranks.rank = rank(mean.rank)) %>%
    mutate(mean.final.rank = (rmse.rank + ranks.rank)/2) %>%
    arrange(mean.final.rank) %>%
    write.csv(., paste0(model.performance.folder, 'aggregated.performance.rank.csv'))


}
