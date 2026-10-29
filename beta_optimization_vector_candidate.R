source("risk_computation_vector_candidate.R")

beta_optimization_vector_candidate <- function(
  data,
  beta_vec, update_time,
  gamma_C1, gamma_C2,
  gamma_A,
  omega, eta = 1.7
){
  obj <- function(beta){
    risk_computation_vector_candidate(
      data = data,
      beta_vec = beta_vec, beta_candidate = beta,
      update_time = update_time,
      gamma_C1 = gamma_C1, gamma_C2 = gamma_C2,
      gamma_A = gamma_A,
      omega = omega, eta = eta
    )
  }

  opt1 <- optimize(obj, interval = c(-50, -5))
  opt2 <- optimize(obj, interval = c(-5, 1))
  opt3 <- optimize(obj, interval = c(1, 50))

  vals <- c(opt1$objective, opt2$objective, opt3$objective)

  best <- list(opt1, opt2, opt3)[[which.min(vals)]]

  beta_hat <- best$minimum

  final_risk <- risk_computation_vector_candidate(
    data = data,
    beta_vec = beta_vec, beta_candidate = beta_hat,
    update_time = update_time,
    gamma_C1 = gamma_C1, gamma_C2 = gamma_C2,
    gamma_A = gamma_A,
    omega = omega, eta = eta
    )

  tibble(
    omega = omega,
    beta = beta_hat,
    cost_Y = attr(final_risk, "risk_Y"),
    cost_trt = attr(final_risk, "risk_trt"),
    risk = final_risk,
    pct_treated = attr(final_risk, "pct_trt")
    )
}