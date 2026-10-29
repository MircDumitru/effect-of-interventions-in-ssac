generate_randomness_fixed_covariates <- function(
  C1 = 0.5, C2 = 0,
  N = 1000, J = 10,
  mu0 = -2, mu1 = 1,
  sigma_mu0 = 1, sigma_mu1 = 0.5,
  sigma_c1 = 1, p_c2 = 0.6,
  tau = 2, sigma_eps = 0.5    
){
  n_total <- N * J
  
  # Time(s) and id(s)
  time <- rep(seq_len(J), N)
  id   <- rep(seq_len(N), each = J)

  # Random effects (one draw per subject)
  mu0_i <- rnorm(N, mu0, sigma_mu0)
  mu1_i <- rnorm(N, mu1, sigma_mu1)

  # Delta Brownians (i.e. Normal) and measurement error (Normal)
  deltaB  <- rnorm(n_total, 0, tau)
  eps <- rnorm(n_total, 0, sigma_eps)

  # Uniform draws (for logisitc approach)
  unif <- runif(n_total)

  # Baseline covariates
  # C2 <- rep(rbinom(N, 1, p_c2), each = J)
  # C1 <- rep(rnorm(N, 0, sigma_c1), each = J)

  tibble(
      id = id,
      time = time,
      mu0_i = rep(mu0_i, each = J),
      mu1_i = rep(mu1_i, each = J),
      C1 = C1, 
      C2 = C2, 
      deltaB = deltaB,
      eps = eps,
      unif = unif
  )
}