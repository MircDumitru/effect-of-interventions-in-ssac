generate_threshold_treatment <- function(
  tb,        
  gamma_C1 = 0.3, gamma_C2 = 1,
  gamma_A = -4, 
  beta = 0
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

      if(t < J){
        A[k + 1] <- as.numeric(Z[k] > beta)
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
