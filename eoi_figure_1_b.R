rm(list = ls())

library(tidyverse)

source("plot_trajectories_observational.R")


#--------------------------------------------------------
out_tb5_1 <- readRDS("results/result_tau_060_sigma_eps_0150.rds")
out_tb5_1$tb_cohort
mean(out_tb5_1$tb_cohort$A)
set.seed(1)
random_sample <- sample(
  unique(out_tb5_1$tb_cohort$id), 
  30, 
  replace = FALSE
)
tb_logisitc_low_noise <- out_tb5_1$tb_cohort |> 
  filter(id %in% random_sample)
# Plot trajectories for observational data
plot_trajectories_observational(tb_logisitc_low_noise, size = 12)
ggsave(
  "figures/logistic_persistent_low_noise.pdf",
  plot = plot_trajectories_observational(tb_logisitc_low_noise, size = 22),
  width = 6,
  height = 4
)
#--------------------------------------------------------


#--------------------------------------------------------
out_tb5_2 <- readRDS("results/result_tau_100_sigma_eps_0250.rds")
out_tb5_2$tb_cohort
mean(out_tb5_2$tb_cohort$A)
set.seed(1)
random_sample <- sample(
  unique(out_tb5_2$tb_cohort$id), 
  30, 
  replace = FALSE
)
tb_logisitc_moderate_noise <- out_tb5_2$tb_cohort |> 
  filter(id %in% random_sample)
# Plot trajectories for observational data
plot_trajectories_observational(tb_logisitc_moderate_noise, size = 12)
ggsave(
  "figures/logistic_persistent_moderate_noise.pdf",
  plot = plot_trajectories_observational(tb_logisitc_moderate_noise, size = 22),
  width = 6,
  height = 4
)
#--------------------------------------------------------


#--------------------------------------------------------
out_tb5_3 <- readRDS("results/result_tau_200_sigma_eps_0500.rds")
out_tb5_3$tb_cohort
mean(out_tb5_3$tb_cohort$A)
set.seed(1)
random_sample <- sample(
  unique(out_tb5_3$tb_cohort$id), 
  30, 
  replace = FALSE
)
tb_logisitc_high_noise <- out_tb5_3$tb_cohort |> 
  filter(id %in% random_sample)

# Plot trajectories for observational data
plot_trajectories_observational(tb_logisitc_high_noise, size = 12)
ggsave(
  "figures/logistic_persistent_high_noise.pdf",
  plot = plot_trajectories_observational(tb_logisitc_high_noise, size = 22),
  width = 6,
  height = 4
)
#--------------------------------------------------------
