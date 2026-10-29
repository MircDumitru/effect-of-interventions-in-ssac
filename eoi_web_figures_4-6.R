rm(list = ls())

library(tidyverse)

source("generate_randomness.R")
source("generate_logistic_treatment_persistent.R")
source("estimate_BMMM.R")
source("generate_randomness_fixed_covariates.R")
source("beta_optimization.R")
source("update_beta_vector.R")
source("beta_optimization_vector_candidate.R")
source("generate_threshold_treatment_vector_beta.R")
source("risk_computation_vector_candidate.R")
source("posterior_update.R")



# Global settings
# -------------------------------------------------------------
snr_values <- c(0, 5, 10, 15, 20)

## All updates are retained internally because the DTDR strategy at
## time 20 depends on the preceding sequential threshold updates.
update_times_all <- c(1, 2, 3, 5, 10, 15, 20)

## Only these update times are reported in the principal figures.
update_times_keep <- c(1, 5, 20)

n_patient_rep <- 100

N_cohort <- 1000
J_cohort <- 10
N_predictive <- 500
N_evaluation <- 5000
J_adapt <- 50

omega <- 0.4
eta <- 1.7
eps_near <- 0.005

C1_new <- 0.5
C2_new <- 1

mu0_true <- -2.5
mu1_true <- 0.5

gamma_C1_true <- 0.3
gamma_C2_true <- 1
gamma_A_true <- -4


# Convert target SNR to (tau, sigma_eps)
# -------------------------------------------------------------
make_noise_grid <- function(
  snr_db,
  p_treated = 0.53,
  tau_eps_ratio = 4,
  mu1_i = 0.5,
  C1 = 0.5,
  C2 = 1,
  gamma_C1 = 0.3,
  gamma_C2 = 1,
  gamma_A = -4
) {
  d0 <- mu1_i + gamma_C1 * C1 + gamma_C2 * C2
  d1 <- d0 + gamma_A

  signal_power <- (1 - p_treated) * d0^2 + p_treated * d1^2

  tibble(snr_db = snr_db) |>
    mutate(
      sigma_eps = sqrt(signal_power / ((tau_eps_ratio^2 + 2) * 10^(snr_db / 10))),
      tau = tau_eps_ratio * sigma_eps
    )
}

noise_grid <- make_noise_grid(snr_values)

print(noise_grid)


# Evaluate a constant threshold on a common true-patient sample
# -------------------------------------------------------------
evaluate_beta <- function(
  beta,
  tb_eval,
  fit,
  omega = 0.4,
  eta = 1.7
) {
  beta_vec <- rep(beta, max(tb_eval$time))

  dat <- generate_threshold_treatment_vector_beta(
    tb = tb_eval,
    beta_vec = beta_vec,
    gamma_C1 = fit$gamma_C1,
    gamma_C2 = fit$gamma_C2,
    gamma_A = fit$gamma_A
  )

  risk_Y <- mean(dat$Y > eta)
  risk_trt <- omega * mean(dat$A)

  risk_Y + risk_trt
}


# Fit one cohort model for each SNR
# -------------------------------------------------------------
prepare_snr_scenario <- function(
  snr_db,
  tau,
  sigma_eps,
  cohort_seed,
  evaluation_seed
) {    
  # Learning cohort
  # -------------------------------------------------------------
  set.seed(cohort_seed)
  tb_noise <- generate_randomness(
    N = N_cohort,
    J = J_cohort,
    mu0 = -2,
    mu1 = 1,
    sigma_mu0 = 1,
    sigma_mu1 = 0.5,
    sigma_c1 = 1,
    p_c2 = 0.6,
    tau = tau,
    sigma_eps = sigma_eps
  )
  tb_cohort <- generate_logistic_treatment_persistent(
    tb_noise,
    gamma_C1 = gamma_C1_true,
    gamma_C2 = gamma_C2_true,
    gamma_A = gamma_A_true
  )
  fit <- estimate_BMMM(tb_cohort)

  # CROUP optimization under the estimated population model
  # -------------------------------------------------------------
  set.seed(evaluation_seed)
  tb_croup <- generate_randomness_fixed_covariates(
    C1 = C1_new,
    C2 = C2_new,
    N = N_evaluation,
    J = J_adapt,
    mu0 = fit$mu0,
    mu1 = fit$mu1,
    sigma_mu0 = fit$sigma_mu0,
    sigma_mu1 = fit$sigma_mu1,
    tau = fit$tau,
    sigma_eps = fit$sigma_eps
  )

  res_croup <- beta_optimization(
    data = tb_croup,
    gamma_C1 = fit$gamma_C1,
    gamma_C2 = fit$gamma_C2,
    gamma_A = fit$gamma_A,
    omega = omega,
    eta = eta
  )
  
  # Common true-patient evaluation sample
  # -------------------------------------------------------------
  set.seed(evaluation_seed + 1L)
  tb_eval <- generate_randomness_fixed_covariates(
    C1 = C1_new,
    C2 = C2_new,
    N = N_evaluation,
    J = J_adapt,
    mu0 = mu0_true,
    mu1 = mu1_true,
    sigma_mu0 = 0,
    sigma_mu1 = 0,
    tau = fit$tau,
    sigma_eps = fit$sigma_eps
  )

  # CROKP threshold and minimum risk on the common sample
  res_oracle <- beta_optimization(
    data = tb_eval,
    gamma_C1 = fit$gamma_C1,
    gamma_C2 = fit$gamma_C2,
    gamma_A = fit$gamma_A,
    omega = omega,
    eta = eta
  )

  # Evaluate the CROUP-selected threshold on the same true-patient
  # sample used for CROKP and DTDR performance comparisons.
  risk_croup_true_patient <- evaluate_beta(
    beta = res_croup$beta,
    tb_eval = tb_eval,
    fit = fit,
    omega = omega,
    eta = eta
  )

  list(
    snr_db = snr_db,
    tau_true = tau,
    sigma_eps_true = sigma_eps,

    fit = fit,

    beta_croup = res_croup$beta,

    # Criterion minimized internally by CROUP under the estimated
    # population predictive distribution.
    risk_croup_predictive = res_croup$risk,

    # Risk incurred by the CROUP threshold for the true new patient.
    risk_croup_true_patient = risk_croup_true_patient,

    beta_oracle = res_oracle$beta,
    risk_oracle = res_oracle$risk,

    tb_eval = tb_eval
  )
}

scenario_list <- pmap(
  noise_grid,
  function(snr_db, tau, sigma_eps) {
    prepare_snr_scenario(
      snr_db = snr_db,
      tau = tau,
      sigma_eps = sigma_eps,
      cohort_seed = 1000L + as.integer(snr_db),
      evaluation_seed = 2000L + as.integer(snr_db)
    )
  }
)

names(scenario_list) <- paste0("snr_", snr_values)



## One new-patient DTDR realization
# -------------------------------------------------------------
run_one_patient <- function(
  scenario,
  replicate_id,
  patient_seed
) {
  fit <- scenario$fit

  set.seed(patient_seed)

  # One true new-patient path
  # -------------------------------------------------------------
  tb_true <- generate_randomness_fixed_covariates(
    C1 = C1_new,
    C2 = C2_new,
    N = 1,
    J = J_adapt,
    mu0 = mu0_true,
    mu1 = mu1_true,
    sigma_mu0 = 0,
    sigma_mu1 = 0,
    tau = fit$tau,
    sigma_eps = fit$sigma_eps
  )

  # Before individual observations are available, treatment starts
  # from the CROUP threshold.
  beta_vec <- rep(scenario$beta_croup, J_adapt)

  tb_current <- generate_threshold_treatment_vector_beta(
    tb = tb_true,
    beta_vec = beta_vec,
    gamma_C1 = fit$gamma_C1,
    gamma_C2 = fit$gamma_C2,
    gamma_A = fit$gamma_A
  )

  prior_mean <- c(fit$mu0, fit$mu1)

  prior_cov <- diag(c(fit$sigma_mu0^2, fit$sigma_mu1^2))

  output <- vector(
    mode = "list",
    length = length(update_times_all)
  )

  # CROUP-to-CROKP excess risk for the true new patient
  croup_oracle_gap <- scenario$risk_croup_true_patient - scenario$risk_oracle

  for (k in seq_along(update_times_all)) {
    t_update <- update_times_all[k]
      
    # Posterior update
    # -------------------------------------------------------------

    tb_observed <- tb_current |>
      filter(time <= t_update)

    ## The original cohort-derived prior is combined with all
    ## observations available at the current update time.
    ## This avoids counting earlier observations more than once.
    post <- posterior_update(
      data = tb_observed,
      prior_mean = prior_mean,
      prior_cov = prior_cov,
      gamma_C1 = fit$gamma_C1,
      gamma_C2 = fit$gamma_C2,
      gamma_A = fit$gamma_A,
      C1 = C1_new,
      C2 = C2_new,
      tau = fit$tau,
      sigma_eps = fit$sigma_eps
    )

    # Posterior-predictive optimization
    # -------------------------------------------------------------

    set.seed(patient_seed + 10000L + k)

    tb_post <- generate_randomness_fixed_covariates(
      C1 = C1_new,
      C2 = C2_new,
      N = N_predictive,
      J = J_adapt,
      mu0 = post$posterior_mean[1],
      mu1 = post$posterior_mean[2],
      sigma_mu0 = sqrt(post$posterior_cov[1, 1]),
      sigma_mu1 = sqrt(post$posterior_cov[2, 2]),
      tau = fit$tau,
      sigma_eps = fit$sigma_eps
    )

    res <- beta_optimization_vector_candidate(
      data = tb_post,
      beta_vec = beta_vec,
      update_time = t_update,
      gamma_C1 = fit$gamma_C1,
      gamma_C2 = fit$gamma_C2,
      gamma_A = fit$gamma_A,
      omega = omega,
      eta = eta
    )

    # True-patient evaluation of selected threshold
    # -------------------------------------------------------------

    risk_true <- evaluate_beta(
      beta = res$beta,
      tb_eval = scenario$tb_eval,
      fit = fit,
      omega = omega,
      eta = eta
    )

    excess_risk <- risk_true - scenario$risk_oracle

    relative_excess_risk <- excess_risk / scenario$risk_oracle

    gap_closed <- if (is.na(croup_oracle_gap) || abs(croup_oracle_gap) < 1e-8) {
      NA_real_
    } else {
        1 - excess_risk / croup_oracle_gap
    }

    output[[k]] <- tibble(
      snr_db = scenario$snr_db,
      replicate = replicate_id,
      update_time = t_update,

      beta_hat = res$beta,
      beta_oracle = scenario$beta_oracle,
      beta_croup = scenario$beta_croup,

      beta_abs_error = abs(res$beta - scenario$beta_oracle),

      risk_true = risk_true,
      risk_oracle = scenario$risk_oracle,
      risk_croup_true_patient = scenario$risk_croup_true_patient,
      risk_croup_predictive = scenario$risk_croup_predictive,

      excess_risk = excess_risk,
      relative_excess_risk = relative_excess_risk,
      gap_closed = gap_closed,

      near_optimal = excess_risk <= eps_near,

      mu0_hat = post$posterior_mean[1],
      mu1_hat = post$posterior_mean[2],

      sd_mu0 = sqrt(post$posterior_cov[1, 1]),
      sd_mu1 = sqrt(post$posterior_cov[2, 2]),

      error_mu0 = post$posterior_mean[1] - mu0_true,

      error_mu1 = post$posterior_mean[2] - mu1_true,

      cover_mu0 = abs(post$posterior_mean[1] - mu0_true ) <= 1.96 * sqrt(post$posterior_cov[1, 1]),

      cover_mu1 = abs(post$posterior_mean[2] - mu1_true) <= 1.96 * sqrt(post$posterior_cov[2, 2])
    )

    # Update only future thresholds
    # -------------------------------------------------------------

    beta_vec <- update_beta_vector(
      beta_vec = beta_vec,
      beta_new = res$beta,
      update_time = t_update + 1
    )

    tb_before <- tb_current

    tb_current <- generate_threshold_treatment_vector_beta(
      tb = tb_true,
      beta_vec = beta_vec,
      gamma_C1 = fit$gamma_C1,
      gamma_C2 = fit$gamma_C2,
      gamma_A = fit$gamma_A
    )

    # The observed history, including the current update visit,
    # must remain unchanged after changing future thresholds.
    stopifnot(
      isTRUE(
        all.equal(
          tb_before |>
            filter(time <= t_update) |>
            arrange(time),

          tb_current |>
            filter(time <= t_update) |>
            arrange(time)
        )
      )
    )
 }

 bind_rows(output) |>
  filter(update_time %in% update_times_keep)
}



# Run Monte Carlo study
# -------------------------------------------------------------

mc_results <- map_dfr(
  seq_along(scenario_list),
  function(s) {    
    scenario <- scenario_list[[s]]
    map_dfr(
      seq_len(n_patient_rep),
        function(r) {
          message("SNR = ", scenario$snr_db, " dB; replicate = ", r, "/", n_patient_rep)
          run_one_patient(
            scenario = scenario,
            replicate_id = r,
            patient_seed = 100000L + 1000L * s + r
          )
        }
    )
  }
)

saveRDS(
  mc_results,
  "results/dtdr_mc_by_snr_corrected.rds"
)



# Summaries
# -------------------------------------------------------------
risk_summary <- mc_results |>
  group_by(
    snr_db,
    update_time
  ) |>
    summarise(
      mean_excess_risk = mean(excess_risk),
      median_excess_risk = median(excess_risk),
      q025_excess_risk = quantile(excess_risk, 0.025),
      q975_excess_risk = quantile(excess_risk, 0.975),
      mean_relative_excess = mean(relative_excess_risk),
      mean_gap_closed = mean(gap_closed, na.rm = TRUE),
      median_gap_closed = median(gap_closed, na.rm = TRUE),
      probability_near_optimal = mean(excess_risk <= eps_near),
      .groups = "drop"
    )

parameter_summary <- mc_results |>
  group_by(
    snr_db,
    update_time
  ) |>
  summarise(
      bias_mu0 = mean(error_mu0),
      rmse_mu0 =  sqrt(mean(error_mu0^2)),
      coverage_mu0 =   mean(cover_mu0),
      bias_mu1 = mean(error_mu1),
      rmse_mu1 = sqrt(mean(error_mu1^2)),
      coverage_mu1 = mean(cover_mu1),
      .groups = "drop"
    )


# Figure: excess risk
# -------------------------------------------------------------
risk_plot_data <- risk_summary |>
  mutate(
    update_time = factor(
      update_time,
      levels = c(1, 5, 20),
      labels = c("DTDR 1 obs.", "DTDR 5 obs.", "DTDR 20 obs.")
    )
  )

p_excess_risk <- ggplot(
  risk_plot_data,
  aes(
    x = snr_db,
    y = mean_excess_risk,
    colour = update_time,
    fill = update_time,
    group = update_time
  )
) +
  geom_ribbon(
  aes(
    ymin = q025_excess_risk,
    ymax = q975_excess_risk
  ),
  alpha = 0.12,
  colour = NA
  ) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.6) +
  scale_x_continuous(breaks = snr_values) +
  scale_y_continuous(
    labels = scales::label_number(accuracy = 0.0001)
  ) +
  scale_colour_brewer(palette = "Blues", direction = -1) +
  scale_fill_brewer(palette = "Blues", direction = -1) +
  labs(
      x = "Signal-to-noise ratio (dB)",
      y = "Mean excess risk",
      colour = NULL,
      fill = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")
p_excess_risk

ggsave(
    "figures/web_fig_4_dtdr_excess_risk_by_snr_corrected.pdf",
    p_excess_risk,
    width = 6,
    height = 4
)


# Parameter-estimation diagnostic
# -------------------------------------------------------------

parameter_plot_data <- parameter_summary |>
  select(
    snr_db,
    update_time,
    rmse_mu0,
    rmse_mu1
  ) |>
  pivot_longer(
    cols = starts_with("rmse_"),
    names_to = "parameter",
    values_to = "rmse"
  ) |>
  mutate(
    parameter = recode(
      parameter,
      rmse_mu0 = "mu[0*i]",
      rmse_mu1 = "mu[1*i]"
    ),
    update_time = factor(
      update_time,
      levels = c(1, 5, 20),
      labels = c("1 observation", "5 observations", "20 observations")
    )
  )


p_parameter_rmse <- ggplot(
  parameter_plot_data,
  aes(
    x = snr_db,
    y = rmse,
    colour = update_time,
    group = update_time
  )
) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  facet_wrap(
    ~ parameter,
    scales = "fixed",
    labeller = label_parsed
  ) +
  scale_x_continuous(breaks = snr_values) +
  scale_colour_brewer(
    palette = "Blues",
    direction = -1
  ) +
  labs(
    x = "Signal-to-noise ratio (dB)",
    y = "Root mean squared error",
    colour = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

p_parameter_rmse

ggsave(
  "figures/web_fig_5_dtdr_parameter_rmse_by_snr.pdf",
  p_parameter_rmse,
  width = 6,
  height = 4
)



# Figure: full tolerance sensitivity
# -------------------------------------------------------------
ecdf_plot_data <- mc_results |>
  mutate(
    update_time = factor(
      update_time,
      levels = c(1, 5, 20),
      labels = c("1 observation", "5 observations", "20 observations")
    ),
    snr_label = factor(
      snr_db,
      levels = snr_values,
      labels = paste0(snr_values," dB")
    )
  )

p_excess_ecdf <- ggplot(
  ecdf_plot_data,
  aes(
    x = excess_risk,
    colour = update_time,
    group = update_time
  )
) +
  stat_ecdf(
    geom = "step",
    linewidth = 1.1
  ) +
  geom_vline(
    xintercept = eps_near,
    linetype = "dashed",
    linewidth = 0.7
  ) +
  facet_wrap(
    ~ snr_label,
    ncol = 3
  ) +
  scale_x_continuous(
    breaks = c(0, 0.0025, 0.005),
    labels = c("0", "0.0025", "0.005")
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.25),
    labels = scales::percent_format()
  ) +
  scale_colour_brewer(
    palette = "Blues",
    direction = -1
  ) +
  labs(
    x = expression(zeta ~ "(excess-risk tolerance)"),
    y = expression(Pr(R(hat(beta)[j]) - R(beta^star) <= zeta)),
    colour = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

p_excess_ecdf

ggsave(
  "figures/web_fig_6_dtdr_excess_risk_ecdf.pdf",
  p_excess_ecdf,
  width = 6,
  height = 4
)
