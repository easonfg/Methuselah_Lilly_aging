model.type = 'elastic.net'
# model.type = 'lightGBM'
model.type = 'SVM'
test.fold.i = 1
iter.model = read_csv(paste0('figures/model.iter.12/', model.type, '/both/both.model.', test.fold.i, '.csv'))
iter.model %>% filter(rmse == min(rmse))
# iter.model %>% select(alpha, lambda) %>% distinct()

iter.model.scale = cbind(iter.model, iter.model %>% scale()) %>% data.frame() %>% tibble()
iter.model.scale = iter.model.scale %>% mutate(sum.metric = -rmse.1+auc.1 )
iter.model.scale %>% filter(sum.metric == max(sum.metric))
# iter.model.scale %>% select(!contains('.1')) %>% arrange(desc(sum.metric)) %>% view

iter.model.long = iter.model %>%
  pivot_longer(cols = c(rmse, auc, error), names_to = 'metric')
iter.model.long

# best.metrics.row = rbind(iter.model.long %>%
#                            filter(metric == 'auc') %>%
#                            group_by(alpha, lambda) %>%
#                            slice_max(value) %>% ungroup(),
#                          iter.model.long %>%
#                            filter(metric != 'auc') %>%
#                            group_by(alpha, lambda) %>%
#                            slice_min(value) %>% ungroup()
# )
# best.metrics.row
#
# iter.model %>%
#   # .[1:1000,] %>%
#   mutate(auc = auc * 100) %>%
#   pivot_longer(cols = c(rmse, auc, error), names_to = 'metric') %>%
#   ggplot()+
#   geom_point(aes(a, value, color = metric)) +
#   geom_vline(data = best.metrics.row
#              , aes(xintercept = j, color = metric))

iter.model.scale %>%
  arrange(desc(sum.metric)) %>%
  # .[1:100,] %>%
  # arrange(auc) %>%
  arrange(auc) %>%
  mutate(j = 1:n()) %>%
  mutate(auc = auc * 100) %>%
  pivot_longer(cols = c(rmse, auc, error), names_to = 'metric') %>%
  ggplot()+
  ylim(c(0, 80)) +
  geom_point(aes(j, value, color = metric))

# +
#   geom_vline(data = best.metrics.row
#              , aes(xintercept = j, color = metric))

