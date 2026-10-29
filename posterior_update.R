posterior_update <- function(
  data,
  prior_mean, prior_cov,
  gamma_C1 = 0.331, gamma_C2 = 1.055, 
  gamma_A = -4.010,
  C1 = 0.5, C2 = 1,
  tau = 2.05, sigma_eps = 0.410
){
  times <- data$time - 1
  X <- cbind(1, times)
  
  # cumulative treatment effect
  trt_effect <- c(0, cumsum(data$A[-1]))
  cvec <- (gamma_C1*C1 + gamma_C2*C2)*times + 
    gamma_A * trt_effect
  
  Sigma <- outer(times, times, pmin) * tau^2
  diag(Sigma) <- diag(Sigma) + sigma_eps^2
  Sigma_inv <- solve(Sigma)

  Omega0_inv <- solve(prior_cov)
  Omega_post <- solve(t(X) %*% Sigma_inv %*% X + Omega0_inv)
  nu_post <- as.vector(
    Omega_post %*% 
    (
      t(X) %*% Sigma_inv %*% (data$Z - cvec) + 
      Omega0_inv %*% prior_mean
    )
  )
  
  list(
    posterior_mean = nu_post,
    posterior_cov  = Omega_post
  )
}