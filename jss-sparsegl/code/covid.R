library(sparsegl)
library(tidyverse)

data("trust_experts")

y <- trust_experts[,"y"]
x <- trust_experts[,-1]
gr <- c(rep(1:7, times = c(8, 51, 5, 4, 8, 10, 10)))
fit <- sparsegl(x, y, gr)
er <- estimate_risk(fit, x, y, approx_df = FALSE)
cc <- coef(fit, s = er$lambda[which.min(er$BIC)])
states <- tibble(state = rownames(cc)[10:60], coef = cc[10:60]) %>%
  mutate(state_name = tolower(covidcast::abbr_to_name(state, TRUE)))
states_map <- map_data("state")
g <- ggplot(states, aes(map_id = state_name)) +
  geom_map(aes(fill = coef), map = states_map) +
  expand_limits(x = states_map$long, y = states_map$lat) +
  scale_fill_gradient2(
    low = "darkorange", high = "darkblue",
    trans = scales::pseudo_log_trans(),
    breaks = c(-10, -5, -2, 0, 2, 5, 10, 20),
  ) +
  coord_map(projection = "albers", parameters = c(30,40)) +
  theme_void() +
  theme(legend.position = "bottom", legend.key.width = unit(2, "cm"),
        legend.title = element_blank())
