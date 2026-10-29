generate_logistic_treatment_persistent <- function(
  tb,    
  gamma_C1 = 0.3, gamma_C2 = 1,
  gamma_A = -4
){
  N <- length(unique(tb$id))
  J <- nrow(tb) / length(unique(tb$id))

  Y <- numeric(nrow(tb))
  Z <- numeric(nrow(tb))
  A <- numeric(nrow(tb))

  k <- 1

  for(i in seq_len(N)){
    Y[k] <- tb$mu0_i[k]
    Z[k] <- Y[k] + tb$eps[k]
    A[k] <- 0

    k <- k + 1
    A[k] <- 0

    for(t in 2:J){
      Y[k] <- Y[k-1] +
        (
          tb$mu1_i[k] +
          gamma_C1 * tb$C1[k] +
          gamma_C2 * tb$C2[k] +
          gamma_A * A[k]
        ) +
        tb$deltaB[k]

      Z[k] <- Y[k] + tb$eps[k]

      # treatment assignment
      if (t < J){
        p <- plogis(-3 + 2 * Z[k] + 0.3 * tb$C1[k] + 0.5 * tb$C2[k])
        A[k+1] <- ifelse(A[k] == 1, 1, tb$unif[k] <= p)
      }

      k <- k + 1
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
