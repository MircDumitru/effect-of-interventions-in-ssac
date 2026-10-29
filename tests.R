rm(list = ls())

library(tidyverse)

source("results_from_observational_cohort.R")
# a <- results_from_observational_cohort()
# a


# source("code_r/generate_randomness.R")
# source("code_r/generate_logistic_treatment_persistent.R")
# source("code_r/estimate_BMMM.R")
# source("code_r/generate_randomness_fixed_covariates.R")
# source("code_r/beta_optimization_2.R")
# source("code_r/update_beta_vector.R")
# source("code_r/beta_optimization_vector_beta_candidate.R")
# source("code_r/generate_threshold_treatment_vector_beta.R")
# source("code_r/risk_computation_vector_beta_candidate.R")
# source("code_r/posterior_update_2.R")
# source("code_r/plot_trajectories.R")
# source("code_r/results_from_observational_cohort.R")

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
  file = "tables/table_5_kable_2.tex"
)
#--------------------------------------------------------





















# source("generate_randomness.R")
# source("generate_logistic_treatment_persistent.R")
# source("estimate_BMMM.R")
# source("generate_randomness_fixed_covariates.R")
# source("beta_optimization.R")
# source("generate_threshold_treatment_vector_beta.R")
# source("beta_optimization_vector_candidate.R")
# source("posterior_update.R")

# verbose = TRUE
# seed = 40

# # cohort
# N_cohort = 1000
# J_cohort = 10

# # true cohort parameters
# mu0 = -2
# mu1 = 1
# sigma_mu0 = 1
# sigma_mu1 = 0.5
# sigma_c1 = 1
# p_c2 = 0.6
# tau = 0.5
# sigma_eps = 0.25

# gamma_A = -4
# gamma_C1 = 0.3
# gamma_C2 = 1

# # adaptation
# J_adapt = 50
# omega = 0.4
# C1_new = 0.5
# C2_new = 1
# mu0_true = -2.5
# mu1_true = 0.5

# update_times = c(1, 2, 3, 5, 10, 15, 20)

# use_double_counting = FALSE


# #--------------------------------------------------------
# tb_noise <- generate_randomness(
#   N = N_cohort,
#   J = J_cohort,
#   mu0 = mu0,
#   mu1 = mu1,
#   sigma_mu0 = sigma_mu0,
#   sigma_mu1 = sigma_mu1,
#   sigma_c1 = sigma_c1,
#   p_c2 = p_c2,
#   tau = tau,
#   sigma_eps = sigma_eps
# )

# tb <- generate_logistic_treatment_persistent(
#   tb_noise,
#   gamma_C1 = gamma_C1,
#   gamma_C2 = gamma_C2,
#   gamma_A = gamma_A
# )

# fit_estimates <- estimate_BMMM(tb)

# mu0_estim <- fit_estimates$mu0
# mu1_estim <- fit_estimates$mu1
# gamma_A_estim <- fit_estimates$gamma_A
# gamma_C1_estim <- fit_estimates$gamma_C1
# gamma_C2_estim <- fit_estimates$gamma_C2
# sigma_mu0_estim <- fit_estimates$sigma_mu0
# sigma_mu1_estim <- fit_estimates$sigma_mu1
# sigma_eps_estim <- fit_estimates$sigma_eps
# tau_estim <- fit_estimates$tau

# parameter_check <- tibble(
#   parameter = c(
#     "mu0", "mu1", 
#     "gamma_A", 
#     "gamma_C1", "gamma_C2",
#     "sigma_mu0", "sigma_mu1", 
#     "sigma_eps", "tau"
#   ),
  
#   true = c(
#     mu0, mu1, 
#     gamma_A, 
#     gamma_C1, gamma_C2,
#     sigma_mu0, sigma_mu1,
#     sigma_eps, tau
#   ),
  
#   estimated = c(
#     mu0_estim, mu1_estim, 
#     gamma_A_estim, 
#     gamma_C1_estim, gamma_C2_estim,
#     sigma_mu0_estim, sigma_mu1_estim,
#     sigma_eps_estim, tau_estim
#   ),

#   bias = true - estimated
# )
# parameter_check
# #--------------------------------------------------------


# # CRA-KP
# #--------------------------------------------------------
# if(verbose){
#   message(
#     "CRA-KP computation: tau = ", 
#     tau, ", sigma_eps = ", sigma_eps
#     )
# }

# N_crakp = 5000

# tb_crakp <- generate_randomness_fixed_covariates(
#   C1 = C1_new, C2 = C2_new,
#   N = N_crakp, J = J_adapt,
#   mu0 = mu0_true, mu1 = mu1_true,
#   sigma_mu0 = 0, sigma_mu1 = 0,
#   tau = tau_estim, sigma_eps = sigma_eps_estim
# )
# head(tb_crakp)

# res_crakp <- beta_optimization(
#   data = tb_crakp,
#   gamma_C1 = gamma_C1_estim,
#   gamma_C2 = gamma_C2_estim,
#   gamma_A = gamma_A_estim,
#   omega = omega
# )
# res_crakp

# table_crakp <- tibble(
#   quantity = c(
#     "beta", "risk", 
#     "mu0", "mu1", 
#     "sigma_mu0", "sigma_mu1"
#   ),
#   CRAKP = c(
#     res_crakp$beta, res_crakp$risk,
#     mu0_true, mu1_true,
#     0, 0
#   )
# )
# table_crakp
# #--------------------------------------------------------


# # CRA-UP
# #--------------------------------------------------------
# if(verbose){
#   message(
#     "CRA-UP computation: tau = ", 
#     tau, ", sigma_eps = ", sigma_eps
#     )
# }

# N_craup = 5000

# tb_craup <- generate_randomness_fixed_covariates(
#   C1 = C1_new, C2 = C2_new,
#   N = N_craup, J = J_adapt,
#   mu0 = mu0_estim, mu1 = mu1_estim,
#   sigma_mu0 = sigma_mu0_estim, sigma_mu1 = sigma_mu1_estim,
#   tau = tau_estim, sigma_eps = sigma_eps_estim
# )
# head(tb_craup)

# res_craup <- beta_optimization(
#   data = tb_craup,
#   gamma_C1 = gamma_C1_estim,
#   gamma_C2 = gamma_C2_estim,
#   gamma_A = gamma_A_estim,
#   omega = omega
# )
# res_craup

# table_craup <- tibble(
#   quantity = c(
#     "beta", "risk",
#     "mu0", "mu1",
#     "sigma_mu0", "sigma_mu1"
#   ),
#   CRAUP = c(
#     res_craup$beta, res_craup$risk,
#     mu0_estim, mu1_estim,
#     sigma_mu0_estim, sigma_mu1_estim
#   )
# )
# table_craup


# # One true patient
# #--------------------------------------------------------
# tb_true <- generate_randomness_fixed_covariates(
#   C1 = C1_new, C2 = C2_new,
#   N = 1, J = J_adapt,
#   mu0 = mu0_true, mu1 = mu1_true,
#   sigma_mu0 = 0, sigma_mu1 = 0,
#   tau = tau_estim, sigma_eps = sigma_eps_estim
# )
# head(tb_true)

# beta_vec <- rep(res_craup$beta, J_adapt)
# beta_vec

# tb_current <- generate_threshold_treatment_vector_beta(
#   tb = tb_true,
#   beta_vec = beta_vec,
#   gamma_C1 = gamma_C1_estim,
#   gamma_C2 = gamma_C2_estim,
#   gamma_A = gamma_A_estim
# )
# tb_current


# # DTDR
# #--------------------------------------------------------
# N_dtdr = 5000

# prior_mean <- c(mu0_estim, mu1_estim)

# prior_cov <- matrix(
#   c(sigma_mu0_estim^2, 0,
#     0, sigma_mu1_estim^2),
#   nrow = 2,
#   byrow = TRUE
# )

# prior_mean_current <- prior_mean
# prior_cov_current <- prior_cov

# results <- vector("list", length(update_times))

# for(k in seq_along(update_times)){
#   t_update <- update_times[k]  
#   if(verbose){    
#     message(      
#       "DTDR computation: tau = ", tau,
#       ", sigma_eps = ", sigma_eps,
#       ", observed up to time = ", t_update
#     )
#   }      

#   tb_observed <- tb_current |>
#     filter(time <= t_update)

#   if(use_double_counting){
#     pm <- prior_mean_current
#     pc <- prior_cov_current
#     } else{
#         pm <- prior_mean
#         pc <- prior_cov
#   }

#   post <- posterior_update(    
#     data = tb_observed,
#     prior_mean = pm, prior_cov = pc,
#     gamma_C1 = gamma_C1_estim, gamma_C2 = gamma_C2_estim,
#     gamma_A = gamma_A_estim,
#     C1 = C1_new, C2 = C2_new,
#     tau = tau_estim, sigma_eps = sigma_eps_estim
#   )

#   tb_post <- generate_randomness_fixed_covariates(
#     C1 = C1_new,
#     C2 = C2_new,
#     N = N_dtdr,
#     J = J_adapt,
#     mu0 = post$posterior_mean[1],
#     mu1 = post$posterior_mean[2],
#     sigma_mu0 = sqrt(post$posterior_cov[1, 1]),
#     sigma_mu1 = sqrt(post$posterior_cov[2, 2]),
#     tau = tau_estim,
#     sigma_eps = sigma_eps_estim
#     )

#   res <- beta_optimization_vector_candidate(
#     data = tb_post,
#     beta_vec = beta_vec,
#     update_time = t_update,
#     gamma_C1 = gamma_C1_estim,
#     gamma_C2 = gamma_C2_estim,
#     gamma_A = gamma_A_estim,
#     omega = omega
#     )

#   results[[k]] <- tibble(
#     update_time = t_update,
#     beta = res$beta,
#     risk = res$risk,
#     mu0 = post$posterior_mean[1],
#     mu1 = post$posterior_mean[2],
#     sigma_mu0 = sqrt(post$posterior_cov[1, 1]),
#     sigma_mu1 = sqrt(post$posterior_cov[2, 2])
#   )

#   tb_before <- tb_current

#   beta_vec <- update_beta_vector(
#     beta_vec = beta_vec,
#     beta_new = res$beta,
#     update_time = t_update + 1
#   )

#   tb_current <- generate_threshold_treatment_vector_beta(
#     tb = tb_true,
#     beta_vec = beta_vec,
#     gamma_C1 = gamma_C1_estim,
#     gamma_C2 = gamma_C2_estim,
#     gamma_A = gamma_A_estim
#   )

#     stopifnot(
#       isTRUE(
#         all.equal(
#           tb_before |>
#             filter(time < t_update) |>
#             arrange(time),
#           tb_current |>
#             filter(time < t_update) |>
#             arrange(time)
#         )
#       )
#     )

#   prior_mean_current <- post$posterior_mean
#   prior_cov_current <- post$posterior_cov
# }

# results_dtdr <- bind_rows(results)
# results_dtdr

# table_dtdr <- results_dtdr |>
#   mutate(update_time = paste0("DTDR_", update_time)) |>
#   pivot_longer(
#     cols = -update_time,
#     names_to = "quantity",
#     values_to = "value"
#   ) |>
#   pivot_wider(
#     names_from = update_time,
#     values_from = value
#   )

# table_dtdr

# table_final <- table_crakp |>
#   left_join(table_craup, by = "quantity") |>
#   left_join(table_dtdr, by = "quantity")
# table_final


# # list(
# #   table_final = table_final,
# #   results_dtdr = results_dtdr,
# #   parameter_check = parameter_check,
# #   fit_estimates = fit_estimates,
# #   res_crakp = res_crakp,
# #   res_craup = res_craup,
# #   beta_vec = beta_vec,
# #   tb_cohort = tb,
# #   tb_true = tb_true,
# #   tb_final = tb_current,
# #   tb_crakp = tb_crakp,
# #   tb_craup = tb_craup,
# #   tb_post = tb_post,
# #   noise = tibble(
# #     tau_true = tau,
# #     sigma_eps_true = sigma_eps,
# #     tau_estimated = tau_estim,
# #     sigma_eps_estimated = sigma_eps_estim
# #   )
# # )