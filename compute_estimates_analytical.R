compute_estimates_analytical <- function(
  tb, times_keep,
  design_name,
  Omega0, nu0,
  gamma_C1 = 0.3, gamma_C2 = 1,
  gamma_A = 0,
  C1 = 0.5, C2 = 1,
  tau = 2, sigma_eps = 0.5
){
  tb_full <- tb |>
    arrange(id, time) |>
    group_by(id) |>
    mutate(
      time0 = time - 1,
      trt_effect = c(0, cumsum(A[-1]))
    ) |>
    ungroup()

    tb_sub <- tb_full |>
      filter(time %in% times_keep)

    ids <- sort(unique(tb_sub$id))

    mu0_hat <- numeric(length(ids))
    mu1_hat <- numeric(length(ids))
    sd_mu0_hat <- numeric(length(ids))
    sd_mu1_hat <- numeric(length(ids))
    true_mu0 <- numeric(length(ids))
    true_mu1 <- numeric(length(ids))

    for(ii in seq_along(ids)){
      data_i <- tb_sub |>
        filter(id == ids[ii]) |>
        arrange(time)

        times <- data_i$time0

        X <- cbind(1, times)

        cvec <- (gamma_C1 * C1 + gamma_C2 * C2) * times +
          gamma_A * data_i$trt_effect

        Sigma <- outer(times, times, pmin) * tau^2
        diag(Sigma) <- diag(Sigma) + sigma_eps^2

        Sigma_inv <- solve(Sigma)
        Omega0_inv <- solve(Omega0)

        Omega_post <- solve(
          t(X) %*% Sigma_inv %*% X +
          Omega0_inv
        )

        nu_post <- as.vector(
          Omega_post %*%
            (
              t(X) %*% Sigma_inv %*% (data_i$Z - cvec) +
              Omega0_inv %*% nu0
            )
        )

        mu0_hat[ii] <- nu_post[1]
        mu1_hat[ii] <- nu_post[2]

        sd_mu0_hat[ii] <- sqrt(Omega_post[1, 1])
        sd_mu1_hat[ii] <- sqrt(Omega_post[2, 2])

        true_mu0[ii] <- data_i$mu0_i[1]
        true_mu1[ii] <- data_i$mu1_i[1]
    }

    tibble(
        design = design_name,

        bias_mu0 = mean(mu0_hat - true_mu0),
        post_sd_mu0 = mean(sd_mu0_hat),
        empir_sd_bias_mu0 = sd(mu0_hat - true_mu0),

        bias_mu1 = mean(mu1_hat - true_mu1),
        post_sd_mu1 = mean(sd_mu1_hat),
        empir_sd_bias_mu1 = sd(mu1_hat - true_mu1)
    )
}