library(covBM)

estimate_BMMM <- function(dat){
  dat <- dat |>
  group_by(id) |>
  mutate(
    cumA = cumsum(A),
    tnew = time - 1
  ) |>
  ungroup()

  dat <- as.data.frame(dat)

  assign("dat", dat, envir = .GlobalEnv)

  BM_model <- try(
    lmeBM(
      fixed = Z ~ 1 + tnew + tnew*C1 + tnew*C2 + cumA - C1 - C2,
      random = ~ 1 + tnew | id,
      data = dat,
      covariance = covBM(form = ~ tnew | id)
    ),
    silent = TRUE
  )

  if(inherits(BM_model, "try-error")){
    return(
      tibble(
        mu0 = NA_real_,
        mu1 = NA_real_,
        gamma_A = NA_real_,
        gamma_C1 = NA_real_,
        gamma_C2 = NA_real_,

        sd_mu0 = NA_real_,
        sd_mu1 = NA_real_,
        sd_gamma_A = NA_real_,
        sd_gamma_C1 = NA_real_,
        sd_gamma_C2 = NA_real_,

        sigma_mu0 = NA_real_,
        sigma_mu1 = NA_real_,
        sigma_eps = NA_real_,
        tau = NA_real_
      )
    )
  }


  fixed_eff <- fixef(BM_model)
  sd_fixed_eff <- sqrt(diag(BM_model$varFix))
  var_corr <- VarCorr(BM_model)

  sigma_mu0 <- as.numeric(var_corr[1, "StdDev"])
  sigma_mu1 <- as.numeric(var_corr[2, "StdDev"])
  sigma_eps <- as.numeric(var_corr[3, "StdDev"])

  kappa <- sqrt(exp(coef(BM_model$modelStruct$corStruct)))
  tau <- sigma(BM_model) * kappa

  tibble(
    mu0 = unname(fixed_eff["(Intercept)"]),
    mu1 = unname(fixed_eff["tnew"]),
    gamma_A = unname(fixed_eff["cumA"]),
    gamma_C1 = unname(fixed_eff["tnew:C1"]),
    gamma_C2 = unname(fixed_eff["tnew:C2"]),

    sd_mu0 = unname(sd_fixed_eff["(Intercept)"]),
    sd_mu1 = unname(sd_fixed_eff["tnew"]),
    sd_gamma_A = unname(sd_fixed_eff["cumA"]),
    sd_gamma_C1 = unname(sd_fixed_eff["tnew:C1"]),
    sd_gamma_C2 = unname(sd_fixed_eff["tnew:C2"]),

    sigma_mu0 = sigma_mu0,
    sigma_mu1 = sigma_mu1,
    sigma_eps = sigma_eps,
    tau = tau
  )
}