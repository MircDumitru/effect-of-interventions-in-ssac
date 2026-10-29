library(covBM)

source("bootstrap_subjects.R")

bootstrap_BMMM <- function(
  dat, 
  B = 100
){
  dat <- dat |>
    group_by(id) |>
    mutate(
      cumA = cumsum(A),
      tnew = time - 1
    ) |>
    ungroup()

  dat <- as.data.frame(dat)
  boot_estimates <- vector("list", B)

  for(b in seq_len(B)){
    cat("Bootstrap", b, "/", B, "\n")
    boot_dat <- bootstrap_subjects(dat)
    assign("boot_dat", boot_dat, envir = .GlobalEnv)

    fit <- try(
      lmeBM(
        fixed = Z ~ 1 + tnew + tnew*C1 + tnew*C2 + cumA - C1 - C2,
        random = ~ 1 + tnew | id,
        data = boot_dat,
        covariance = covBM(form = ~ tnew | id)
      ),
      silent = TRUE
    )

    if(inherits(fit, "try-error")){

      boot_estimates[[b]] <- tibble(
        mu0 = NA_real_,
        mu1 = NA_real_,
        gamma_A = NA_real_,
        gamma_C1 = NA_real_,
        gamma_C2 = NA_real_,
        sigma_mu0 = NA_real_,
        sigma_mu1 = NA_real_,
        sigma_eps = NA_real_,
        tau = NA_real_
      )
      next
    }

    fixed_eff <- fixef(fit)
    var_corr <- VarCorr(fit)

    sigma_mu0 <- as.numeric(var_corr[1,"StdDev"])
    sigma_mu1 <- as.numeric(var_corr[2,"StdDev"])
    sigma_eps <- as.numeric(var_corr[3,"StdDev"])

    kappa <- sqrt(exp(coef(fit$modelStruct$corStruct)))
    tau <- sigma(fit) * kappa

    boot_estimates[[b]] <- tibble(

      mu0 = unname(fixed_eff["(Intercept)"]),
      mu1 = unname(fixed_eff["tnew"]),
      gamma_A = unname(fixed_eff["cumA"]),
      gamma_C1 = unname(fixed_eff["tnew:C1"]),
      gamma_C2 = unname(fixed_eff["tnew:C2"]),

      sigma_mu0 = sigma_mu0,
      sigma_mu1 = sigma_mu1,
      sigma_eps = sigma_eps,
      tau = tau
    )
  }

  boot_estimates <- bind_rows(
    boot_estimates
  )

  tibble(
    boot_sd_mu0 = sd(boot_estimates$mu0, na.rm = TRUE),
    boot_sd_mu1 = sd(boot_estimates$mu1, na.rm = TRUE),
    boot_sd_gamma_A = sd(boot_estimates$gamma_A, na.rm = TRUE),
    boot_sd_gamma_C1 = sd(boot_estimates$gamma_C1, na.rm = TRUE),
    boot_sd_gamma_C2 = sd(boot_estimates$gamma_C2, na.rm = TRUE),
    
    boot_sd_sigma_mu0 = sd(boot_estimates$sigma_mu0, na.rm = TRUE),
    boot_sd_sigma_mu1 = sd(boot_estimates$sigma_mu1, na.rm = TRUE),
    boot_sd_sigma_eps = sd(boot_estimates$sigma_eps, na.rm = TRUE),
    boot_sd_tau = sd(boot_estimates$tau, na.rm = TRUE)
  )
}
