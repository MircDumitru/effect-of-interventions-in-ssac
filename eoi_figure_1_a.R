rm(list = ls())

library(tidyverse)

source("generate_randomness.R")
source("plot_trajectories.R")
source("generate_threshold_treatment.R")


set.seed(40)
tb_randomness  <- generate_randomness(N = 30)
head(tb_randomness)


tb_never <- generate_threshold_treatment(tb_randomness, beta = Inf)
tb_threshold <- generate_threshold_treatment(tb_randomness, beta = 0)
tb_always <- generate_threshold_treatment(tb_randomness, beta = -Inf)


head(tb_never)
head(tb_threshold)
head(tb_always)


# Plot trajectories for all three strategies
plot_trajectories(tb_never, size = 12)
ggsave(
  "figures/never.pdf",
  plot = plot_trajectories(tb_never, size = 22),
  width = 6,
  height = 4
)

plot_trajectories(tb_threshold, size = 12)
ggsave(
  "figures/threshold.pdf",
  plot = plot_trajectories(tb_threshold, size = 22),
  width = 6,
  height = 4
)

plot_trajectories(tb_always, size = 12)
ggsave(
  "figures/always.pdf",
  plot = plot_trajectories(tb_always, size = 22),
  width = 6,
  height = 4
)