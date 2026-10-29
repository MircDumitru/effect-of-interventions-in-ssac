rm(list = ls())

library(tidyverse)
library(knitr)

source("generate_randomness_fixed_covariates.R")
source("generate_threshold_treatment.R")
source("compute_estimates_analytical.R")

set.seed(1)

C1 <- 0.5 
C2 <- 1
N <- 10000
J <- 100
mu0 <- -2.5
mu1 <- 1
sigma_mu0 <- 1
sigma_mu1 <- 1
tau = 2
sigma_eps = 0.5

prior_mean_correct <- c(mu0, mu1)
prior_cov_correct <- diag(c(sigma_mu0^2, sigma_mu1^2))

prior_mean_misspec <- c(-2, 2)
prior_cov_misspec <- diag(c(0.5, 2))


tb_randomness <- generate_randomness_fixed_covariates(
  C1 = C1, C2 = C2, 
  N = N, J = J, 
  mu0 = mu0, mu1 = mu1, 
  sigma_mu0 = sigma_mu0, sigma_mu1 = sigma_mu1, 
  tau = tau, sigma_eps = sigma_eps
)
# tb_randomness


gamma_C1 = 0.3
gamma_C2 = 1
gamma_A = -4
# gamma_A = 0 
beta = 0
# beta = Inf
tb  <- generate_threshold_treatment(
  tb_randomness,
  gamma_C1 = gamma_C1, gamma_C2 = gamma_C2,
  gamma_A = gamma_A, beta = beta
)
# tb 
tb <- tb |>
  left_join(
    tb_randomness |>
      distinct(id, mu0_i, mu1_i),
    by = "id"
  )

designs <- list(
  hor_3_samp_1 = 1:3,
  hor_5_samp_1 = 1:5,
  hor_10_samp_2 = seq(1,10,by=2),  
  hor_10_samp_1 = 1:10,  
  hor_100_samp_5 = seq(1,100,by=5),
  hor_100_samp_2 = seq(1,100,by=2),  
  hor_100_samp_1 = 1:100  
)


nu0 <- prior_mean_correct
Omega0 <- prior_cov_correct
results_correct <- bind_rows(  
  lapply(
    names(designs),
    function(nm){
      compute_estimates_analytical(        
        tb = tb,
        times_keep = designs[[nm]],
        design_name = nm,
        Omega0 = Omega0, nu0 = nu0,
        gamma_C1 = gamma_C1, gamma_C2 = gamma_C2,
        gamma_A = gamma_A,
        C1 = C1, C2 = C2,
        tau = tau, sigma_eps = sigma_eps
      )
    }
  )
)
results_correct

nu0 <- prior_mean_misspec
Omega0 <- prior_cov_misspec
results_misspec <- bind_rows(  
  lapply(
    names(designs),
    function(nm){
      compute_estimates_analytical(
        tb = tb,
        times_keep = designs[[nm]],
        design_name = nm,
        Omega0 = Omega0, nu0 = nu0,
        gamma_C1 = gamma_C1, gamma_C2 = gamma_C2,
        gamma_A = gamma_A,
        C1 = C1, C2 = C2,
        tau = tau, sigma_eps = sigma_eps
      )
    }
  )
)
results_misspec


table <- bind_rows(results_correct, results_misspec)
table_final <- table |>
  mutate(
    across(
      where(is.numeric),
      ~ round(.x, 3)
    )
  )
table_final


cat(
  knitr::kable(
    table_final,
    format = "latex",
    booktabs = TRUE
  ),
  file = "tables/table_4_kable.tex"
)