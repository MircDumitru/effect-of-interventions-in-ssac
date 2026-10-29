rm(list = ls())

library(tidyverse)
library(knitr)

source("one_BMMM_estimation_bootstrap.R")

set.seed(1)

n_sim <- 1000 # 1000
results <- map_dfr(
  1:n_sim,
  ~ one_BMMM_estimation_bootstrap(
      N = 1000, # 1000 
      B = 100 # 100
    )
)

# True parameters
true_params <- c(
  mu0 = -2, mu1 = 1,
  gamma_A = -4,
  gamma_C1 = 0.3, gamma_C2 = 1,
  sigma_mu0 = 1, sigma_mu1 = 0.5,
  tau = 2, sigma_eps = 0.5
)

parameters <- c(
  "mu0", "mu1",
  "gamma_A",
  "gamma_C1", "gamma_C2",
  "sigma_mu0", "sigma_mu1",
  "tau", "sigma_eps"
)

# Table
table <- purrr::map_dfr(
  parameters,
  function(par){
    est <- results[[par]]
    out <- tibble(
      Parameter = par,
      Bias = mean(est - true_params[[par]], na.rm = TRUE),
      `Empirical SD` = sd(est, na.rm = TRUE)
    )    
    # fixed effects
    sd_name <- paste0("sd_", par)
    if(sd_name %in% names(results)){
      est_sd <- results[[sd_name]]
      out <- out |>
        mutate(
          `Estimated SD` = mean(est_sd, na.rm = TRUE),
          Coverage = 
            mean(
              est - 1.96 * est_sd <= true_params[[par]] &
              est + 1.96 * est_sd >= true_params[[par]],
              na.rm = TRUE
            )
        )
    } else {
      out <- out |>
        mutate(
          `Estimated SD` = NA_real_,
          Coverage = NA_real_
        )
    }    
    # bootstrap sd
    boot_name <- paste0("boot_sd_", par)
    if(boot_name %in% names(results)){
      boot_sd <- results[[boot_name]]
      out <- out |>
        mutate(
          `Bootstrap SD` =
            mean(
              boot_sd,
              na.rm = TRUE
            ),
          `Bootstrap Coverage` =
            mean(
              est - 1.96 * boot_sd <= true_params[[par]] &
              est + 1.96 * boot_sd >= true_params[[par]],
              na.rm = TRUE
            )
        )
    } else {
      out <- out |>
        mutate(
          `Bootstrap SD` = NA_real_,
          `Bootstrap Coverage` = NA_real_
        )
    }
    out
  }
)

table_final <- table |>
  mutate(
    across(
      where(is.numeric),
      ~ round(.x, 4)
    )
  )

table_final

table_results <- table_final |>
  mutate(
    across(
      where(is.numeric),
      ~ round(.x, 3)
    )
  )

table_results

write_csv(
  table_results, 
  "results/table_1.csv"
)

saveRDS(
  table_results,
  "results/table_1.rds"
)


cat(
  knitr::kable(
    table_results,
    format = "latex",
    booktabs = TRUE
  ),
  file = "tables/table_1_kable.tex"
)
