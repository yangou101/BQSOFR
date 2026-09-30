# Back-transform existing posterior draws; this function never fits a model.
# FUI/Raw: supply fpca_reparam to recover the full FPCA curve intercept.
# Without fpca_reparam, auto keeps the fitted (centered-FPCA) predictor.
# Joint: auto preserves the original full-X intercept correction.
# Rebuild FPCA using the SAME input matrix, row/time order, package versions
# and settings as the fit. Prefer a retained original FPCA object if available.
# The saved fitted beta basis, not a newly oriented beta basis, is used with
# posterior betaf/betar. The integration rule is the fitted sum(...)/T.
backtransform_scalar_coefs_fui_raw <- function(
    fit,
    Y_use,
    Z_use,
    W_use,
    numeric_Z = c("AgeYR", "HEI"),
    factor_Z  = c("Gender", "Race", "HealthCondt2"),
    reg_formula,
    tau0 = NULL,
    intercept_type = c("auto", "fitted", "fullX"),
    fpca_reparam = NULL
) {
  
  library(rstan)

  intercept_type <- match.arg(intercept_type)
  first_attr <- function(keys) {
    for (key in keys) {
      value <- attr(fit, key, exact = TRUE)
      if (!is.null(value)) return(value)
    }
    NULL
  }
  valid_tau <- function(x) {
    is.numeric(x) && length(x) == 1L && is.finite(x) && x > 0 && x < 1
  }
  saved_tau <- Filter(Negate(is.null), lapply(
    c("tau0", "loo_tau", "tau_compare"),
    function(key) attr(fit, key, exact = TRUE)
  ))
  if (length(saved_tau) && !all(vapply(saved_tau, valid_tau, logical(1)))) {
    stop("Saved quantile metadata must contain single numeric values in (0, 1).")
  }
  if (is.null(tau0) && length(saved_tau)) tau0 <- saved_tau[[1L]]
  if (!valid_tau(tau0)) {
    stop("Supply tau0 as a single number in (0, 1), matching the fitted quantile.")
  }
  if (length(saved_tau) && any(abs(unlist(saved_tau) - tau0) > 1e-10)) {
    stop("tau0 disagrees with the fitted quantile metadata.")
  }

  is_fui <- identical(attr(fit, "loo_method", exact = TRUE), "FUI") ||
    !is.null(attr(fit, "B_beta_compare", exact = TRUE))
  if (intercept_type == "auto") {
    intercept_type <- if (is_fui && is.null(fpca_reparam)) "fitted" else "fullX"
  }
  if (is_fui && intercept_type == "fullX" && is.null(fpca_reparam)) {
    stop("Supply fpca_reparam = build_fpca_reparam(W_used_for_this_fit) for the FUI/Raw fullX intercept.")
  }
  x_cov_saved <- attr(fit, "X_cov", exact = TRUE)
  if (intercept_type == "fullX" && !is.null(x_cov_saved) &&
      ncol(x_cov_saved) > 0L &&
      !isTRUE(attr(fit, "x_covariate_in_outcome", exact = TRUE))) {
    stop("This Joint fit uses residual PA only. Use intercept_type = 'fitted'; a fullX conversion would also change scalar coefficients.")
  }
  
  
  ## ============================================================
  ## 1. Reproduce the exact analysis sample
  ## ============================================================
  
  Z_df <- as.data.frame(Z_use)
  
  for (nm in factor_Z) {
    Z_df[[nm]] <- as.factor(Z_df[[nm]])
  }
  
  idx <- complete.cases(Y_use, Z_df) &
    apply(
      W_use,
      1,
      function(x) all(is.finite(x))
    )
  
  Y0 <- as.numeric(Y_use[idx])
  Z0 <- Z_df[idx, , drop = FALSE]

  saved_y <- attr(fit, "loo_y", exact = TRUE)
  if (!is.null(saved_y) &&
      (length(saved_y) != length(Y0) ||
       !isTRUE(all.equal(as.numeric(saved_y), Y0, check.attributes = FALSE)))) {
    stop("The supplied data do not reproduce the fitted Y sample and row order.")
  }
  
  
  ## ============================================================
  ## 2. Means / SDs used in standardization
  ## ============================================================
  
  y_mean <- first_attr(c("Y_mean", "Y_mean_compare", "loo_y_mean"))
  y_sd   <- first_attr(c("Y_sd", "Y_sd_compare", "loo_y_sd"))
  
  if (is.null(y_mean)) {
    y_mean <- mean(Y0)
  }
  
  if (is.null(y_sd)) {
    y_sd <- sd(Y0)
  }
  if (length(y_mean) != 1L || !is.finite(y_mean) ||
      length(y_sd) != 1L || !is.finite(y_sd) || y_sd <= 0) {
    stop("Invalid fitted outcome mean or SD.")
  }
  
  
  ## continuous covariate means
  z_mean <- sapply(
    Z0[, numeric_Z, drop = FALSE],
    mean
  )
  
  ## continuous covariate SDs
  z_sd <- sapply(
    Z0[, numeric_Z, drop = FALSE],
    sd
  )
  
  
  ## ============================================================
  ## 3. Reconstruct exact scalar design matrix
  ## ============================================================
  
  Z_numeric_scaled <- as.data.frame(
    scale(
      Z0[, numeric_Z, drop = FALSE],
      center = TRUE,
      scale = TRUE
    )
  )
  
  Z_factor <- Z0[, factor_Z, drop = FALSE]
  
  Z_final_model <- data.frame(
    Z_numeric_scaled,
    Z_factor
  )
  
  df2 <- data.frame(
    Y = (Y0 - y_mean) / y_sd,
    Z_final_model
  )
  
  
  ## scalar covariates
  z_terms <- setdiff(
    colnames(df2),
    "Y"
  )
  
  
  ## Variables appearing inside interactions
  z_main_no_inter <- setdiff(
    z_terms,
    c(
      "Race",
      "Gender",
      "AgeYR",
      "HEI"
    )
  )
  
  
  ## EXACT same formula as fitted model
  #form_Z <- paste0(
  #  "Y ~ ", 
  #  paste(c(
  #      z_main_no_inter,
  #      "Race * Gender",
  #      "Race * AgeYR",
  #      "HEI * Race"
  #    ),
  #    collapse = " + "
  #  )
  #)
  
  #form_Z <- "Y ~ HealthCondt2 + Race * Gender + Race * AgeYR"
  form_Z <- reg_formula
  
  saved_formula <- attr(fit, "loo_formula", exact = TRUE)
  if (!is.null(saved_formula)) {
    current_terms <- stats::terms(stats::as.formula(form_Z))
    fitted_terms <- stats::terms(stats::as.formula(saved_formula))
    if (!identical(attr(current_terms, "term.labels"), attr(fitted_terms, "term.labels")) ||
        !identical(attr(current_terms, "intercept"), attr(fitted_terms, "intercept"))) {
      stop("reg_formula must match the fitted scalar formula, including term order.")
    }
  }

  ## Prefer the exact column names saved by the fitting function.
  b_names <- attr(fit, "scalar_names_compare", exact = TRUE)
  if (is.null(b_names)) {
    saved_X <- attr(fit, "loo_X", exact = TRUE)
    if (!is.null(colnames(saved_X))) b_names <- colnames(saved_X)[-1L]
  }
  if (is.null(b_names)) {
    Dt2 <- brms::make_standata(
      brms::bf(as.formula(form_Z), quantile = tau0),
      data = df2,
      family = brms::asym_laplace()
    )
    b_names <- colnames(Dt2$X)[-1L]
  }
  
  
  ## ============================================================
  ## 4. Extract posterior draws
  ##
  ## IMPORTANT:
  ## a_x       = coefficients of population mean alpha_0(t)
  ## b_func_gq = functional beta mapped into Bx coordinates
  ## ============================================================
  
  post <- rstan::extract(
    fit,
    pars = c("b", "b_Intercept",
             if (intercept_type == "fullX") {
               if (is_fui) c("betaf", "betar") else c("a_x", "b_func_gq")
             }),
    permuted = TRUE
  )
  
  
  B_std <- as.matrix(post$b)
  
  stopifnot(
    ncol(B_std) == length(b_names)
  )
  
  colnames(B_std) <- b_names
  
  
  ## FUI has no sampled a_x or b_func_gq. Joint full-X correction uses both.
  if (intercept_type == "fullX" && !is_fui) {
    a_x_draws <- as.matrix(post$a_x)
    b_func_draws <- as.matrix(post$b_func_gq)
    stopifnot(
      nrow(a_x_draws) == nrow(B_std),
      nrow(b_func_draws) == nrow(B_std),
      ncol(a_x_draws) == ncol(b_func_draws)
    )
  }
  
  
  ## ============================================================
  ## 5. Restore Y scale for scalar coefficients
  ## ============================================================
  
  B_orig <- B_std * y_sd
  
  
  ## ============================================================
  ## 6. Restore Age / HEI original units
  ##
  ## Age main + Race x Age   : divide by SD(Age)
  ## HEI main + Race x HEI   : divide by SD(HEI)
  ## ============================================================
  
  for (j in seq_along(b_names)) {
    
    pieces <- strsplit(
      b_names[j],
      ":",
      fixed = TRUE
    )[[1]]
    
    
    for (nm in numeric_Z) {
      
      if (nm %in% pieces) {
        
        B_orig[, j] <-
          B_orig[, j] /
          z_sd[nm]
      }
    }
  }
  
  
  ## ============================================================
  ## 7. Correct Race main effects
  ##
  ## Because:
  ##
  ## Age_std = (Age - mean_age) / sd_age
  ## HEI_std = (HEI - mean_hei) / sd_hei
  ##
  ## Race main effects must absorb the centering terms from
  ## Race x Age and Race x HEI.
  ## ============================================================
  
  race_main <- grep(
    "^Race[^:]+$",
    colnames(B_orig),
    value = TRUE
  )
  
  
  find_interaction <- function(
    race_name,
    numeric_name
  ) {
    
    candidates <- c(
      paste0(race_name, ":", numeric_name),
      paste0(numeric_name, ":", race_name)
    )
    
    hit <- intersect(
      candidates,
      colnames(B_orig)
    )
    
    if (length(hit) == 0) {
      return(NA_character_)
    }
    
    hit[1]
  }
  
  
  for (race_name in race_main) {
    
    ## --------------------------
    ## Race x Age
    ## --------------------------
    
    nm_age <- find_interaction(
      race_name,
      "AgeYR"
    )
    
    if (!is.na(nm_age)) {
      
      B_orig[, race_name] <-
        B_orig[, race_name] -
        B_orig[, nm_age] *
        z_mean["AgeYR"]
    }
    
    
    ## --------------------------
    ## Race x HEI
    ## --------------------------
    
    nm_hei <- find_interaction(
      race_name,
      "HEI"
    )
    
    if (!is.na(nm_hei)) {
      
      B_orig[, race_name] <-
        B_orig[, race_name] -
        B_orig[, nm_hei] *
        z_mean["HEI"]
    }
  }
  
  
  ## ============================================================
  ## 8. FUNCTIONAL POPULATION MEAN alpha_0(t)
  ##
  ## Joint fitted outcome model uses:
  ##
  ##   gamma_i' b_func
  ##
  ## but full latent curve is:
  ##
  ##   X_i = a_x + gamma_i
  ##
  ## therefore:
  ##
  ##   gamma_i' b_func
  ##      =
  ##   X_i' b_func - a_x' b_func
  ##
  ## Thus the intercept for the FULL-X representation is:
  ##
  ##   b_Intercept - a_x' b_func
  ##
  ## IMPORTANT: calculate draw-by-draw.
  ## FUI/Raw instead use Xi_fix * b_func, with
  ## Xhat_i(t) = fpca.fit$mu(t) + Phi(t) * Xi_fix[i,].
  ## Its mean correction is sum(mu(t) * beta_std(t)) / T.
  ## ============================================================
  
  intercept_fitted_std <- as.numeric(post$b_Intercept)
  alpha0_func_contrib_std <- NULL
  intercept_fullX_std <- NULL
  fpca_basis_linkage_error <- NULL
  functional_mean_curve <- NULL
  if (intercept_type == "fullX") {
    if (is_fui) {
      ## Use the exact saved basis corresponding to the existing draws.
      beta_basis <- attr(fit, "B_beta_compare", exact = TRUE)
      Kf <- attr(fit, "beta_Kf_compare", exact = TRUE)
      Kr <- attr(fit, "beta_Kr_compare", exact = TRUE)
      mu_fpca <- as.numeric(fpca_reparam$fpca.fit$mu)
      Phi <- fpca_reparam$Phi
      scores <- fpca_reparam$Xi_fix
      if (!is.matrix(W_use) || !is.matrix(beta_basis) ||
          !is.matrix(Phi) || !is.matrix(scores) ||
          length(Kf) != 1L || length(Kr) != 1L ||
          !is.finite(Kf) || !is.finite(Kr) || Kf < 1L || Kr < 1L) {
        stop("FUI/Raw fullX requires an N x T W matrix, saved beta basis/dimensions, and the build_fpca_reparam() result.")
      }
      T_num <- nrow(beta_basis)
      if (T_num != ncol(W_use) || length(mu_fpca) != T_num ||
          nrow(Phi) != T_num || ncol(beta_basis) != Kf + Kr ||
          nrow(scores) != length(Y0) || ncol(scores) != ncol(Phi) ||
          any(!is.finite(mu_fpca)) || any(!is.finite(beta_basis)) ||
          any(!is.finite(Phi)) || any(!is.finite(scores))) {
        stop("FPCA reconstruction dimensions/values do not match this fitted sample and beta basis.")
      }
      ## Check the user's reconstructed linkage against the saved fitted basis.
      G_rebuilt <- cbind(t(fpca_reparam$X_mat_f), t(fpca_reparam$X_mat_r))
      G_saved_basis <- crossprod(Phi, beta_basis) / T_num
      if (!identical(dim(G_rebuilt), dim(G_saved_basis)) ||
          any(!is.finite(G_rebuilt))) {
        stop("Reconstructed X_mat_f/X_mat_r do not match the fitted beta dimensions.")
      }
      fpca_basis_linkage_error <- max(abs(G_rebuilt - G_saved_basis)) /
        max(1, max(abs(G_saved_basis)))
      if (fpca_basis_linkage_error > 1e-6) {
        stop("Reconstructed FPCA/spline linkage differs from the saved fit. Check the input W, time order, and package versions/settings.")
      }
      saved_J <- attr(fit, "fpca_components_compare", exact = TRUE)
      saved_ev <- attr(fit, "fpca_evalues_compare", exact = TRUE)
      if (!is.null(saved_J) && !identical(as.integer(saved_J), as.integer(ncol(Phi)))) {
        stop("Reconstructed FPCA component count differs from the fitted object.")
      }
      if (!is.null(saved_ev) &&
          !isTRUE(all.equal(as.numeric(saved_ev),
                            as.numeric(fpca_reparam$fpca.fit$evalues),
                            tolerance = 1e-6, check.attributes = FALSE))) {
        stop("Reconstructed FPCA eigenvalues differ from the fit; use the exact fitted W and FPCA settings.")
      }
      beta_coef <- cbind(matrix(post$betaf, ncol = Kf), matrix(post$betar, ncol = Kr))
      if (nrow(beta_coef) != nrow(B_std)) stop("Posterior draws are not aligned.")
      ## beta_std(t) = beta_basis(t,) %*% c(betaf, betar).
      ## Compute the mean contribution without allocating draws x T curves.
      mean_basis_integral <- as.numeric(crossprod(beta_basis, mu_fpca)) / T_num
      alpha0_func_contrib_std <- as.numeric(beta_coef %*% mean_basis_integral)
      functional_mean_curve <- mu_fpca
    } else {
      alpha0_func_contrib_std <- rowSums(a_x_draws * b_func_draws)
    }
    intercept_fullX_std <- intercept_fitted_std - alpha0_func_contrib_std
    intercept_used_std <- intercept_fullX_std
  } else {
    ## Retain the actual functional predictor used during fitting.
    intercept_used_std <- intercept_fitted_std
  }
  
  
  ## ============================================================
  ## 9. Back-transform intercept to original BMI scale
  ## ============================================================
  
  intercept_orig <-
    y_mean +
    y_sd *
    intercept_used_std
  
  
  ## Undo Age centering
  if ("AgeYR" %in% colnames(B_orig)) {
    
    intercept_orig <-
      intercept_orig -
      B_orig[, "AgeYR"] *
      z_mean["AgeYR"]
  }
  
  
  ## Undo HEI centering
  if ("HEI" %in% colnames(B_orig)) {
    
    intercept_orig <-
      intercept_orig -
      B_orig[, "HEI"] *
      z_mean["HEI"]
  }
  
  
  ## ============================================================
  ## 10. Final original-scale posterior draws
  ## ============================================================
  
  coef_orig <- cbind(
    Intercept = intercept_orig,
    B_orig
  )
  
  
  ## ============================================================
  ## 11. Posterior summaries
  ## ============================================================
  
  summ <- function(x) {
    
    c(
      Mean = mean(x),
      SD = sd(x),
      
      `2.5%` = unname(
        quantile(
          x,
          0.025
        )
      ),
      
      `50%` = unname(
        quantile(
          x,
          0.50
        )
      ),
      
      `97.5%` = unname(
        quantile(
          x,
          0.975
        )
      )
    )
  }
  
  
  tab <- as.data.frame(
    t(
      apply(
        coef_orig,
        2,
        summ
      )
    )
  )
  
  
  ## ============================================================
  ## 12. b[1] ... b[p] mapping
  ## ============================================================
  
  mapping <- data.frame(
    
    parameter = paste0(
      "b[",
      seq_along(b_names),
      "]"
    ),
    
    term = b_names
  )
  
  
  ## ============================================================
  ## 13. Useful intercept diagnostics
  ## ============================================================
  
  intercept_gamma_std <-
    as.numeric(
      post$b_Intercept
    )
  
  
  intercept_gamma_Yscale <-
    y_mean +
    y_sd *
    intercept_gamma_std
  
  
  alpha0_func_contrib_Yscale <- if (is.null(alpha0_func_contrib_std)) NULL else
    y_sd * alpha0_func_contrib_std
  
  
  diagnostic_draws <- list(
    "b_Intercept: fitted functional predictor, standardized Y" = intercept_fitted_std,
    "selected intercept before Age/HEI correction, standardized Y" = intercept_used_std
  )
  if (intercept_type == "fullX") {
    diagnostic_draws[["alpha0 functional contribution, standardized Y"]] <- alpha0_func_contrib_std
    diagnostic_draws[["alpha0 functional contribution, original Y scale"]] <- alpha0_func_contrib_Yscale
  }
  intercept_diagnostics <- data.frame(
    quantity = names(diagnostic_draws),
    mean = vapply(diagnostic_draws, mean, numeric(1)),
    sd = vapply(diagnostic_draws, sd, numeric(1)),
    row.names = NULL
  )
  
  
  ## ============================================================
  ## 14. Return
  ## ============================================================
  
  list(
    tau0 = tau0,
    intercept_type = intercept_type,
    intercept_fitted_std = intercept_fitted_std,
    intercept_used_std = intercept_used_std,
    functional_mean_curve = functional_mean_curve,
    fpca_basis_linkage_error = fpca_basis_linkage_error,
    
    mapping = mapping,
    
    table = tab,
    
    draws = coef_orig,
    
    ## Standardization information
    Y_mean = y_mean,
    Y_sd = y_sd,
    Z_mean = z_mean,
    Z_sd = z_sd,
    
    ## Functional alpha_0 contribution
    alpha0_func_contrib_std =
      alpha0_func_contrib_std,
    
    alpha0_func_contrib_Yscale =
      alpha0_func_contrib_Yscale,
    
    ## Different intercept parameterizations
    intercept_gamma_std =
      intercept_gamma_std,
    
    intercept_fullX_std =
      intercept_fullX_std,
    
    intercept_original =
      intercept_orig,
    
    intercept_diagnostics =
      intercept_diagnostics
  )
}

# Same FPCA and mgcv reparameterization supplied by the user.
# Run separately for FUI curves and Raw daily-mean curves.
build_fpca_reparam <- function(W) {
  if (!is.matrix(W) || !is.numeric(W) || any(!is.finite(W))) {
    stop("W must be the finite numeric N x T matrix used in the fitted model.")
  }
  Xt <- W
  n_num <- nrow(Xt)
  nt <- ncol(Xt)
  tind <- seq(0, 1, length.out = nt)

  fpca.fit <- refund::fpca.sc(Xt)
  Xi_fix <- fpca.fit$scores
  Phi <- fpca.fit$efunctions
  J_num <- ncol(Phi)

  wmat <- I(Xt)
  lmat <- I(matrix(1 / nt, nrow = n_num, ncol = nt))
  tmat <- I(matrix(tind, nrow = n_num, ncol = nt, byrow = TRUE))
  dat_list <- list(tmat = tmat, lmat = lmat, wmat = wmat)
  knots <- NULL
  object <- mgcv::s(tmat, by = lmat * wmat, bs = "cc", k = 10)
  dk <- mgcv:::ExtractData(object, data = dat_list, knots = knots)
  splinecons <- mgcv:::smooth.construct.cc.smooth.spec(object, dk$data, dk$knots)

  Psi_mat <- splinecons$X
  S_mat <- splinecons$S[[1]]
  rank <- splinecons$rank
  M_num <- nrow(Psi_mat)
  K_num <- ncol(Psi_mat)
  if (M_num != nt || nrow(Phi) != nt) {
    stop("The FPCA/spline basis must have one row per fitted time point.")
  }

  X_mat_t <- matrix(0, nrow = J_num, ncol = K_num)
  for (j in 1:J_num) {
    for (k in 1:K_num) {
      X_mat_t[j,k] <- sum(Phi[,j] * Psi_mat[,k]) / M_num
    }
  }
  maXX <- norm(Psi_mat, type = "I")^2
  maS <- norm(S_mat, type = "I") / maXX
  S_mat <- S_mat / maS
  eig <- eigen(S_mat, symmetric = TRUE)
  E <- rep(1, ncol(X_mat_t))
  E[1:rank] <- sqrt(pmax(eig$values[1:rank], 1e-12))
  X_mat_t <- X_mat_t %*% eig$vectors
  col.norm <- colSums(X_mat_t^2)
  col.norm <- col.norm / (E^2)
  av.norm <- mean(col.norm[1:rank])
  if (rank < ncol(X_mat_t)) {
    for (i in (rank + 1):ncol(X_mat_t)) {
      E[i] <- sqrt(col.norm[i] / av.norm)
    }
  }
  if (any(!is.finite(E)) || any(E <= 0)) stop("Invalid beta-basis scaling E.")
  X_mat_t <- t(t(X_mat_t) / E)
  X_mat_r <- t(X_mat_t[, 1:rank, drop = FALSE])
  X_mat_f <- t(X_mat_t[, (rank + 1):ncol(X_mat_t), drop = FALSE])

  list(fpca.fit = fpca.fit, Xi_fix = Xi_fix, Phi = Phi,
       X_mat_f = X_mat_f, X_mat_r = X_mat_r, rank = rank)
}
