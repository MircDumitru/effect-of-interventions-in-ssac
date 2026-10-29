source("generate_randomness.R")
source("generate_logistic_treatment_persistent.R")
source("estimate_BMMM.R")
source("bootstrap_BMMM.R")

one_BMMM_estimation_bootstrap <- function(
  N = 1000, 
  B = 100
){
  tb_noise <- generate_randomness(N = N)
  tb <- generate_logistic_treatment_persistent(tb_noise)
  fit_estimates <- estimate_BMMM(tb)
  boot_estimates <- bootstrap_BMMM(tb, B = B)
  bind_cols(fit_estimates, boot_estimates)
}