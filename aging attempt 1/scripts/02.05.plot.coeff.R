library(tidyverse)
library(parallel)
region.i = 'bottom'
region.i = 'random'
region.i = 'top'

dir.create('figures/coeff', recursive = T)

for (region.i in c('top', 'bottom', 'random')){
  # for (region.i in c('top')){
  # for (model.type in c('elastic.net', 'SVM', 'lightGBM')){
  # model.list = c('elastic.net',  'lightGBM')
  model.list = c(  'lightGBM','SVM', 'elastic.net')
  # model.list = c('SVM')
  # model.list = c(  'elastic.net')
  # for (model.type in model.list){
  all.ranks.list = lapply(model.list, function(model.type){
    print(region.i)
    print(model.type)

    if (model.type == 'elastic.net'){
      var.col = 'name'
      coeff.col = 'coefficient'
    } else if (model.type == 'SVM') {
      var.col = 'Variable'
      coeff.col = ''
    } else if (model.type == 'lightGBM') {
      var.col = 'Feature'
      coeff.col = 'Gain'
    }

    all.lgbm.iter.coeffs = mclapply(list.files(paste0('figures/model.iter.', region.i, '.rmse.12/', model.type, '/both/'),
                                               pattern = 'model.coeff', full.names = T),
                                    function(coeff.i){

                                      mclapply(list.files(coeff.i, full.names = T),
                                               function(iter.i){
                                                 iter.csv = read_csv(iter.i)
                                                 iter.csv %>% mutate(model.n = basename(coeff.i), iter.n = basename(iter.i))
                                               }) %>% data.table::rbindlist()
                                    })

    all.lgbm.iter.coeffs
    # all.lgbm.iter.coeffs[[1]] %>% filter(iter.n == 'iter.1.csv')
    all.lgbm.iter.coeffs = all.lgbm.iter.coeffs %>% data.table::rbindlist()
    all.lgbm.iter.coeffs = all.lgbm.iter.coeffs %>% group_by(model.n, iter.n) %>%
      arrange(desc(!!sym(coeff.col)), .by_group = TRUE) %>% # Sort descending within group
      mutate(sequential_rank = row_number()) %>% ungroup()

    # all.lgbm.iter.coeffs %>% filter(Feature == 'GFAP') %>% view
    # all.lgbm.iter.coeffs %>% filter(Feature == 'CCL11') %>% view

    all.lgbm.iter.coeffs

    # browser()
    lgbm.coeff.rank = all.lgbm.iter.coeffs %>% ungroup() %>% group_by( model.n, !!sym(var.col)) %>%
      summarise(rank = mean(sequential_rank), sd.rank = sd(sequential_rank), freq = n())
    # all.lgbm.iter.coeffs %>% filter(Feature == 'AGER')
    lgbm.coeff.rank

    # Step 1: Calculate group means
    group_means <- lgbm.coeff.rank %>%
      group_by(!!sym(var.col)) %>%
      summarise(mean_value = mean(rank, na.rm = TRUE))
    group_means

    lgbm.coeff.rank = lgbm.coeff.rank %>% left_join(group_means, by = var.col)

    interested.feature = lgbm.coeff.rank %>% arrange(mean_value) %>% pull(!!sym(var.col)) %>% unique() %>% .[1:30]
    p = lgbm.coeff.rank %>%
      filter(!!sym(var.col) %in% interested.feature) %>%
      arrange(mean_value) %>%
      mutate(!!sym(var.col) := fct_inorder(!!sym(var.col))) %>%
      # .[1:50, ] %>%
      ggplot(., aes(!!sym(var.col), rank, fill = model.n)) +
      # geom_col(width = 0.7, alpha = 0.8) +  # Bar plot
      geom_col(position = position_dodge(width = 0.8), width = 0.7) +
      geom_errorbar(aes(ymin = rank - sd.rank, ymax = rank + sd.rank),
                    position = position_dodge(0.9),  # Must match geom_col position!
                    width = 0.2,            # Width of error bar ends
                    size = 0.8,             # Thickness of error bars
                    color = "black") +
      labs(title = paste0(model.type, ', ', region.i),
           x = "Target",
           y = "Mean Rank") +
      coord_cartesian(ylim = c(0, 100)) +
      theme_minimal(base_size = 30) +
      theme(legend.position = "none") +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


    pdf(paste0('figures/coeff/', model.type, '.', region.i, '.pdf'), w = 15)
    print(p)
    dev.off()

    lgbm.coeff.rank = lgbm.coeff.rank %>% mutate(model.type = model.type)
    colnames(lgbm.coeff.rank)[2] = 'name'
    lgbm.coeff.rank
  })

  all.ranks = all.ranks.list %>% data.table::rbindlist()
  # all.ranks %>% filter(model.type == 'elastic.net', name == 'AGER') %>%
  #   group_by(model.type, name) %>%
  #   summarise(sd.rank = sd(rank), rank = mean(rank),  freq = sum(freq))
  all.ranks = all.ranks %>% ungroup() %>% group_by(model.type, name) %>%
    summarise(sd.rank = sd(rank), rank = mean(rank),  freq = sum(freq))

  all.ranks
  for (model.i in model.list){
    rank.p = all.ranks %>%
      mutate(model.type = factor(model.type)) %>%
      mutate(model.type = relevel(model.type, ref = model.i)) %>%
      group_by(model.type) %>% arrange(rank) %>%
      mutate(name = fct_inorder(name)) %>%
      # .[1:50, ] %>%
      slice_head(n = 30) %>% # Will return all rows if group has <50 rows
      ggplot(., aes(name, rank, fill = model.type)) +
      geom_col(position = position_dodge(width = 0.8), width = 0.7) +
      geom_errorbar(aes(ymin = rank - sd.rank, ymax = rank + sd.rank),
                    position = position_dodge(0.9),  # Must match geom_col position!
                    width = 0.2,                    # Width of error bar ends
                    size = 0.8,                     # Thickness of error bars
                    color = "black") +
      scale_fill_manual(values = c("SVM" = "blue",
                                   "elastic.net" = "red",
                                   "lightGBM" = "green")) +
      labs(title = "Bar Plot with Error Bars (±SD)",
           x = "Group",
           y = "Mean Value") +
      theme_minimal(base_size = 30) +
      # theme(legend.position = "none") +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

    pdf(paste0('figures/coeff/', region.i, '.ref.', model.i, '.pdf'), w = 15)
    print(rank.p)
    dev.off()

  }
}
