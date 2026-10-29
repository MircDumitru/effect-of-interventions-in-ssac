rm(list = ls())

library(tidyverse)
library(knitr)

source("generate_randomness_fixed_covariates.R")
source("beta_optimization.R")

omega_grid <- c(seq(0, 1, 0.1), 1.5, 3)


#--------------------------------------------------------
data <- generate_randomness_fixed_covariates(
  C1 = 0.5, 
  C2 = 0
)
table <- purrr::map_dfr(
  omega_grid,
  ~ beta_optimization(
    data = data,
    gamma_C1 = 0.3, gamma_C2 = 1,
    gamma_A = -4,  
    omega = .x
  )
)
table
table_results <- table |>
  mutate(
    across(
      where(is.numeric),
      ~ round(.x, 3)
    )
  )
cat(
  knitr::kable(table_results, format = "latex", booktabs = TRUE),
  file = "tables/table_optim_C1_05_C2_0_unk_kable.tex"
)
cat("Table C1 = 0.5, C2 = 0: done")
#--------------------------------------------------------


#--------------------------------------------------------
data <- generate_randomness_fixed_covariates(
  C1 = 0.5, 
  C2 = 1
)
table <- purrr::map_dfr(
  omega_grid,
  ~ beta_optimization(
    data = data,
    gamma_C1 = 0.3, gamma_C2 = 1,
    gamma_A = -4,      
    omega = .x
  )
)
table
table_results <- table |>
  mutate(
    across(
      where(is.numeric),
      ~ round(.x, 3)
    )
  )
cat(
  knitr::kable(table_results, format = "latex", booktabs = TRUE),
  file = "tables/table_optim_C1_05_C2_1_unk_kable.tex"
)
cat("Table C1 = 0.5, C2 = 1: done")
#--------------------------------------------------------


#--------------------------------------------------------
data <- generate_randomness_fixed_covariates(
  C1 = -0.5, 
  C2 = 0
)
table <- purrr::map_dfr(
  omega_grid,
  ~ beta_optimization(
    data = data,
    gamma_C1 = 0.3, gamma_C2 = 1,
    gamma_A = -4,      
    omega = .x
  )
)
table
table_results <- table |>
  mutate(
    across(
      where(is.numeric),
      ~ round(.x, 3)
    )
  )
cat(
  knitr::kable(table_results, format = "latex", booktabs = TRUE),
  file = "tables/table_optim_C1_m05_C2_0_unk_kable.tex"
)
cat("Table C1 = -0.5, C2 = 0: done")
#--------------------------------------------------------


#--------------------------------------------------------
data <- generate_randomness_fixed_covariates(
  C1 = -0.5, 
  C2 = 1
)
table <- purrr::map_dfr(
  omega_grid,
  ~ beta_optimization(
    data = data,
    gamma_C1 = 0.3, gamma_C2 = 1,
    gamma_A = -4,      
    omega = .x
  )
)
table
table_results <- table |>
  mutate(
    across(
      where(is.numeric),
      ~ round(.x, 3)
    )
  )
cat(
  knitr::kable(table_results, format = "latex", booktabs = TRUE),
  file = "tables/table_optim_C1_m05_C2_1_unk_kable.tex"
)
cat("Table C1 = -0.5, C2 = 1: done")
#--------------------------------------------------------