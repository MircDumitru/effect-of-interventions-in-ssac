source("update_beta_vector.R")

risk_computation_vector_candidate <- function(
  data,
  beta_vec, beta_candidate,
  update_time,
  gamma_C1, gamma_C2,
  gamma_A,
  omega, eta = 1.7
){
  beta_vec_candidate <- update_beta_vector(
    beta_vec = beta_vec,
    beta_new = beta_candidate,
    update_time = update_time + 1
  )

  dat <- generate_threshold_treatment_vector_beta(
    tb = data,
    beta_vec = beta_vec_candidate,
    gamma_C1 = gamma_C1,
    gamma_C2 = gamma_C2,
    gamma_A = gamma_A
  )

  risk_Y <- mean(dat$Y > eta)
  risk_trt <- omega * mean(dat$A)

  risk_total <- risk_Y + risk_trt

  attr(risk_total, "risk_Y") <- risk_Y
  attr(risk_total, "risk_trt") <- risk_trt
  attr(risk_total, "pct_trt") <- mean(dat$A)

  risk_total
}