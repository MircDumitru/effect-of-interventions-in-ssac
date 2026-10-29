source("generate_threshold_treatment.R")

risk_computation <- function(
  data,
  gamma_C1, gamma_C2,
  gamma_A,
  beta, omega,
  eta = 1.7
){  
  dat <- generate_threshold_treatment(
    data, 
    gamma_C1 = gamma_C1,
    gamma_C2 = gamma_C2,
    gamma_A = gamma_A,
    beta = beta
  )

  risk_Y <- mean(dat$Y > eta)
  risk_trt <- omega * mean(dat$A)
  risk_total <- risk_Y + risk_trt

  attr(risk_total, "risk_Y") <- risk_Y
  attr(risk_total, "risk_trt") <- risk_trt
  attr(risk_total, "pct_trt") <- mean(dat$A)

  risk_total
}