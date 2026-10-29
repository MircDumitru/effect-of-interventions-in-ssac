rm(list = ls())

library(tidyverse)
library(knitr)

source("results_from_observational_cohort.R")

#--------------------------------------------------------
make_noise_filename <- function(
  tau, 
  sigma_eps, 
  prefix = "result"
){
  tau_txt <- gsub("\\.", "", sprintf("%.2f", tau))
  eps_txt <- gsub("\\.", "", sprintf("%.3f", sigma_eps))
  paste0(prefix, "_tau_", tau_txt, "_sigma_eps_", eps_txt, ".csv")
}
results_dir <- "results"
#--------------------------------------------------------

# Table 5 for tau = 0.6, sigma_eps = 0.15
#--------------------------------------------------------
tau <- 0.6
sigma_eps <- 0.15
table_5_1 <- results_from_observational_cohort(
  seed = 2, N_cohort = 1000,
  tau = tau, sigma_eps = sigma_eps,
  use_double_counting = FALSE
)
table_5_1$tb_cohort$A  |> mean()
table_5_1$table_final
table_5_1$parameter_check
table_5_1$noise
file_name <- make_noise_filename(tau, sigma_eps)
write_csv(
  table_5_1$table_final,
  file.path(results_dir, file_name)
)
saveRDS(
  table_5_1,
  file.path(
    results_dir,
    sub("\\.csv$", ".rds", file_name)
  )
)
#--------------------------------------------------------

# Table 5 for tau = 1, sigma_eps = 0.25
#--------------------------------------------------------
tau <- 1
sigma_eps <- 0.25
table_5_2 <- results_from_observational_cohort(
  seed = 2, N_cohort = 1000,
  tau = tau, sigma_eps = sigma_eps,
  use_double_counting = FALSE
)
table_5_2$tb_cohort$A  |> mean()
table_5_2$table_final
table_5_2$parameter_check
table_5_2$noise
file_name <- make_noise_filename(tau, sigma_eps)
write_csv(
  table_5_2$table_final,
  file.path(results_dir, file_name)
)
saveRDS(
  table_5_2,
  file.path(
    results_dir,
    sub("\\.csv$", ".rds", file_name)
  )
)
#--------------------------------------------------------

# Table 5 for tau = 2, sigma_eps = 0.5
#--------------------------------------------------------
tau <- 2
sigma_eps <- 0.5
table_5_3 <- results_from_observational_cohort(
  seed = 5, N_cohort = 1000,
  tau = tau, sigma_eps = sigma_eps,
  use_double_counting = FALSE
)
table_5_3$tb_cohort$A  |> mean()
table_5_3$table_final
table_5_3$parameter_check
table_5_3$noise
file_name <- make_noise_filename(tau, sigma_eps)
write_csv(
  table_5_3$table_final,
  file.path(results_dir, file_name)
)
saveRDS(
  table_5_3,
  file.path(
    results_dir,
    sub("\\.csv$", ".rds", file_name)
  )
)
#--------------------------------------------------------

# Final table 5 (biding all three tables/results)
#--------------------------------------------------------
table_5_final <- bind_rows(
  table_5_1$table_final,
  table_5_2$table_final,
  table_5_3$table_final
) |>  
  mutate(
    quantity = recode(      
      quantity,
      beta = "$\\beta$",
      risk = "Risk",
      mu0 = "$\\mu_{0i}$",
      mu1 = "$\\mu_{1i}$",
      sigma_mu0 = "$\\sigma_{\\mu_{0i}}$",
      sigma_mu1 = "$\\sigma_{\\mu_{1i}}$"
    )
  ) |>
  rename(Parameter = quantity) |> 
  mutate(
    across(where(is.numeric), round, 3)
  )
table_5_final
#--------------------------------------------------------

# Save results
#--------------------------------------------------------
cat(
  kable(
    table_5_final, 
    format = "latex", 
    booktabs = TRUE, 
    escape = FALSE
  ),
  file = "tables/table_5_kable.tex"
)
#--------------------------------------------------------
