library(dplyr)
library(tidyr)
library(ggplot2)


plot_box <- function(df) {
  df <- df %>%
    group_by(expr, change_feature, change_list, within_group, type) %>% 
    summarise(mean = mean(mean), .groups = 'drop')
  plot <- df %>% ggplot(aes(x = change_list, y = mean, color = expr, group = expr)) +
    geom_point() +
    geom_line() + 
    theme_bw() + 
    theme(legend.position = "bottom", legend.title = element_blank(), axis.title.x = element_blank()) +
    xlab("Change values") +
    ylab("Running Time (ms)") +
    facet_grid(type + within_group ~ change_feature, 
               labeller = as_labeller(c('p' = "Number of features p", 'n' = "Sample size n",
                                        'TRUE' = "Within-group sparsity", 'FALSE' = "Without-group sparsity",
                                        'linear' = "Linear regression", 'logistic' = "Logistic regression")), switch = "x") +
    scale_colour_viridis_d() +
    scale_y_log10() 
  plot
}

