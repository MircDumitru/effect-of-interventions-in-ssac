rm(list = ls())

library(tidyverse)
library(knitr)

source("risk_computation.R")
source("update_beta_vector.R")
source("risk_computation_vector_candidate.R")
source("generate_threshold_treatment_vector_beta.R")

out_tb5_1 <- readRDS("results/result_tau_060_sigma_eps_0150.rds")
out_tb5_2 <- readRDS("results/result_tau_100_sigma_eps_0250.rds")
out_tb5_3 <- readRDS("results/result_tau_200_sigma_eps_0500.rds")

make_risk_curve_one_noise <- function(
  out,
  noise_label,
  beta_margin = 3,
  beta_step = 0.025,
  eps_tol = 0.005,
  dtdr_time = 20
){
  fit <- out$fit_estimates

  gamma_C1_use <- fit$gamma_C1
  gamma_C2_use <- fit$gamma_C2
  gamma_A_use  <- fit$gamma_A

  dtdr_row <- out$results_dtdr |>
    filter(update_time == dtdr_time)

  stopifnot(nrow(dtdr_row) == 1)

  J_adapt <- max(out$tb_crakp$time)

  ## Reference optima from Table 5 / saved RDS
  ref_tbl <- tibble(
    method = c("CRA-KP", "CRA-UP", paste0("DTDR ", dtdr_time)),
    beta_ref = c(out$res_crakp$beta, out$res_craup$beta, dtdr_row$beta),
    risk_ref = c(out$res_crakp$risk, out$res_craup$risk, dtdr_row$risk)
  )

  beta_min <- min(ref_tbl$beta_ref) - beta_margin
  beta_max <- max(ref_tbl$beta_ref) + beta_margin

  beta_grid <- sort(unique(c(seq(beta_min, beta_max, by = beta_step), ref_tbl$beta_ref)))

  tb_skp_curve  <- out$tb_crakp
  tb_spdp_curve <- out$tb_craup
  tb_dtdr_curve <- out$tb_post

  ## Beta vector before the DTDR update at dtdr_time
  beta_vec_before <- rep(out$res_craup$beta, J_adapt)

  previous_updates <- out$results_dtdr |>
    filter(update_time < dtdr_time) |>
    arrange(update_time)

  if(nrow(previous_updates) > 0){
    for(ii in seq_len(nrow(previous_updates))){
      beta_vec_before <- update_beta_vector(
        beta_vec = beta_vec_before,
        beta_new = previous_updates$beta[ii],
        update_time = previous_updates$update_time[ii] + 1
      )
    }
  }

  
  compute_curve_constant <- function(data, method){
    map_dfr(beta_grid, function(b){
      r <- risk_computation(
        data = data,
        gamma_C1 = gamma_C1_use,
        gamma_C2 = gamma_C2_use,
        gamma_A = gamma_A_use,
        beta = b,
        omega = 0.4
      )

      tibble(
        noise_setting = noise_label,
        method = method,
        beta = b,
        risk = as.numeric(r),
        risk_Y = attr(r, "risk_Y"),
        risk_trt = attr(r, "risk_trt")
      )
    }
    )
  }

  compute_curve_dtdr <- function(data, method){
    map_dfr(beta_grid, function(b){
      r <- risk_computation_vector_candidate(
        data = data,
        beta_vec = beta_vec_before,
        beta_candidate = b,
        update_time = dtdr_time,
        gamma_C1 = gamma_C1_use,
        gamma_C2 = gamma_C2_use,
        gamma_A = gamma_A_use,
        omega = 0.4
      )

      tibble(
        noise_setting = noise_label,
        method = method,
        beta = b,
        risk = as.numeric(r),
        risk_Y = attr(r, "risk_Y"),
        risk_trt = attr(r, "risk_trt")
      )
    }
    )
  }

  curves <- bind_rows(
    compute_curve_constant(tb_skp_curve, "CRA-KP"),
    compute_curve_constant(tb_spdp_curve, "CRA-UP"),
    compute_curve_dtdr(tb_dtdr_curve, paste0("DTDR ", dtdr_time))
  ) |>
    left_join(ref_tbl, by = "method") |>
    mutate(
      excess_risk = risk - risk_ref,
      near_optimal = excess_risk <= eps_tol,
      method = factor(
        method,
        levels = c("CRA-KP", "CRA-UP", paste0("DTDR ", dtdr_time))
      )
    )
  curves
}

plot_risk_curve_one_noise <- function(risk_curve){
  refs <- risk_curve |>
    distinct(method, beta_ref, risk_ref)

  ggplot(risk_curve, aes(x = beta, y = risk, colour = method)) +    
    geom_line(aes(group = method), linewidth = 1.75, alpha = 0.45) +
    geom_line(
      data = risk_curve |> filter(near_optimal),
      aes(group = method),
      linewidth = 2,
      alpha = 1
    ) +
    geom_point(
      data = refs,
      aes(x = beta_ref, y = risk_ref, colour = method),
      size = 3,
      inherit.aes = FALSE
    ) +
    scale_color_brewer(palette = "Dark2") +
    labs(
      x = expression(beta),
      y = "Risk",
      colour = NULL
    ) +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 23),
      axis.text = element_text(size = 23),
      legend.text = element_text(size = 23),
      legend.title = element_text(size = 23),
      legend.position = "bottom",
      plot.title = element_blank(),
      plot.subtitle = element_blank()
    )
}

risk_curve_low <- make_risk_curve_one_noise(
  out = out_tb5_1,
  noise_label = "Low-noise setting",
  beta_margin = 3,
  eps_tol = 0.005,
  dtdr_time = 20
)

risk_curve_moderate <- make_risk_curve_one_noise(
  out = out_tb5_2,
  noise_label = "Moderate-noise setting",
  beta_margin = 3,
  eps_tol = 0.005,
  dtdr_time = 20
)

risk_curve_high <- make_risk_curve_one_noise(
  out = out_tb5_3,
  noise_label = "High-noise setting",
  beta_margin = 3,
  eps_tol = 0.005,
  dtdr_time = 20
)

p_risk_low <- plot_risk_curve_one_noise(risk_curve_low)
p_risk_moderate  <- plot_risk_curve_one_noise(risk_curve_moderate)
p_risk_high  <- plot_risk_curve_one_noise(risk_curve_high)

p_risk_low
p_risk_moderate
p_risk_high

ggsave(
  "figures/web_fig_3_risk_curve_low_noise.pdf",
  p_risk_low,
  width = 6,
  height = 4
)

ggsave(
  "figures/web_fig_3_risk_curve_moderate_noise.pdf",
  p_risk_moderate,
  width = 6,
  height = 4
)

ggsave(
  "figures/web_fig_3_risk_curve_high_noise.pdf",
  p_risk_high,
  width = 6,
  height = 4
)

eps_tol <- 0.005

plateau_summary <- bind_rows(
  risk_curve_low, 
  risk_curve_moderate, 
  risk_curve_high
) |>
  filter(near_optimal) |>
  group_by(noise_setting, method) |>
  summarise(
    `Optimal threshold` = first(beta_ref),
    `Minimum risk` = first(risk_ref),
    `Near-optimal left` = min(beta),
    `Near-optimal right` = max(beta),
    `Interval width` = max(beta) - min(beta),
    .groups = "drop"
    ) |>
    rename(
      `Noise setting` = noise_setting,
      Method = method
    ) |>
    mutate(
      across(where(is.numeric), ~ round(.x, 3))
    )

plateau_summary

eps_txt <- gsub("\\.", "", sprintf("%.4f", eps_tol))

cat(
  kable(
    plateau_summary,
    format = "latex",
    booktabs = TRUE
  ),
  file = file.path(
    "tables",
    paste0("web_table_1_eps_tol_", eps_txt, "_kable.tex")
  )
)
