library(splines)
library(dplyr)
library(sparsegl)
library(tibble)
library(ggplot2)
df <- 10

data("trust_experts", package = "sparsegl")
trust_experts <- trust_experts %>%
  mutate(across(
    where(is.factor),
    ~ `attr<-`(.x, "contrasts", contr.sum(nlevels(.x), FALSE, TRUE))
  ))

x <- Matrix::sparse.model.matrix(
  ~ 0 + region + age + gender + raceethnicity + period +
    bs(cli, df = df) + bs(hh_cmnty_cli, df = df),
  data = trust_experts, drop.unused.levels = TRUE)

gr <- sapply(trust_experts, function(x) ifelse(is.factor(x), nlevels(x), NA))
gr <- rep(seq(ncol(trust_experts) - 1), times = c(gr[!is.na(gr)], df, df))
fit <- cv.sparsegl(x, trust_experts$trust_experts, gr)

cc <- coef(fit, s = "lambda.1se")
reg <- which(substr(rownames(cc), 1, nchar("region")) == "region")
states <- tibble(state = rownames(cc)[reg], coef = cc[reg]) |>
  mutate(state = substring(state, nchar("region") + 1),
         state_name = tolower(covidcast::abbr_to_name(state, TRUE)))
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
