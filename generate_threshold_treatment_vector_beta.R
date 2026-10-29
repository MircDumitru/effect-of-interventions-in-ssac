generate_threshold_treatment_vector_beta <- function(
  tb,
  beta_vec,
  gamma_C1 = 0.331,
  gamma_C2 = 1.055,
  gamma_A  = -4.010
){
  ids <- sort(unique(tb$id))
  J <- max(tb$time)

  stopifnot(length(beta_vec) >= J)

  tb <- tb |>
    arrange(id, time)

  Y <- numeric(nrow(tb))
  Z <- numeric(nrow(tb))
  A <- numeric(nrow(tb))

  for(id_now in ids){
    rows <- which(tb$id == id_now)

    ## time 1
    r1 <- rows[1]

    Y[r1] <- tb$mu0_i[r1]
    Z[r1] <- Y[r1] + tb$eps[r1]
    A[r1] <- 0

    ## time 2
    if(J >= 2){
      A[rows[2]] <- 0
    }

    if(J >= 2){
      for(t in 2:J){
        r_now  <- rows[t]
        r_prev <- rows[t - 1]

        Y[r_now] <- Y[r_prev] + 
          (
            tb$mu1_i[r_now] +
            gamma_C1 * tb$C1[r_now] +
            gamma_C2 * tb$C2[r_now] +
            gamma_A  * A[r_now]
          ) +
          tb$deltaB[r_now]

        Z[r_now] <- Y[r_now] + tb$eps[r_now]

        if(t < J){
          A[rows[t + 1]] <- as.numeric(Z[r_now] > beta_vec[t])
          }
      }
    }
  }

  tibble(
    id = tb$id,
    time = tb$time,
    C1 = tb$C1,
    C2 = tb$C2,
    Y = Y,
    Z = Z,
    A = A
  )
}