# Optional latent-PA covariates:
# Set BOTH x_numeric_Z = NULL and x_factor_Z = NULL for no measurement covariates.
# Then K_xcov = 0, alpha_x has no free elements, and x_score = a_x + gamma.
# x_covariate_in_outcome has no effect in that mode (covariate contribution = 0).
# Existing covariate fits and all priors, including U ~ Beta(4,4), are retained.
#
# 2026-09-23: latent-PA covariate extension.
# The measurement model now allows x_score_i = a_x + X_cov_i alpha_x + gamma_i.
# Use x_numeric_Z / x_factor_Z to define Age/Gender/etc. for latent PA.
# x_covariate_in_outcome = FALSE keeps beta(t) on residual PA gamma only;
# TRUE uses X_cov alpha_x + gamma in the scalar-on-function outcome term.
#
# 2026-09-21: equivalent likelihood acceleration; fast_likelihood = TRUE.
# Measurement likelihood uses subject means AND within-subject sums of squares.
# GAL shared terms are evaluated once per joint log-density evaluation.
# Set fast_likelihood = FALSE to evaluate the original likelihood expressions.
# The posterior target is the same up to parameter-independent constants;
# floating-point summation and realized MCMC trajectories may differ.
# Numerical density/gradient checks were done outside Stan; no R/Stan benchmark
# was available here. Assess performance using effective samples per second.
#
# scalar_formula interface added for FUI/Joint design alignment.
# Entrypoint: fit_flfoqr_fui_bounded_gamma_betaf_xcov().
# Explicit scalar_formula supplies the scalar main effects/interactions.
# NULL retains the legacy hard-coded interaction formula for compatibility.
# The measurement basis, beta basis, priors, and pointwise LOO density
# are retained. The disabled examples use the projected tau=0.1/0.5 priors.
# These outcome-derived anchors are for prior sensitivity analysis.
# Full-pipeline cross-validation must rebuild outcome-derived anchors within
# training folds; saved log_lik_y alone does not account for that selection.
# source() defines functions only; it does not start MCMC.

# Y-only LOO extension of your supplied Joint model.
# The likelihood, priors, function arguments and existing summaries are retained.
# log_lik_y scores only Y, conditional on sampled latent gamma.
# The W likelihood remains in fitting; removing Y_i retains W_i.
# 03_loo_compare.R applies the original-Y-scale density Jacobian.

## BETAF-CENTER EXTENSION (2026-09-13)
## Added beta_fixed_prior_mean; default 0 reproduces the original prior.
## gamma bounds, a_x prior, measurement/outcome likelihoods, and the
## Half-normal sigma_beta prior are retained from the supplied source.
## New main function: fit_flfoqr_fui_bounded_gamma_betaf_xcov().
## source() defines functions only. The run example at the end is disabled.
## A nonzero betaf center is a prior sensitivity setting, not a guaranteed
## vertical displacement of the posterior beta curve.
## ============================================================

## FLFOQR: shared gamma variance and subject-specific omega variance.
## USER-SPECIFIED STRUCTURE (in the existing unit-normalized Bx coordinates):
##   gamma_ir | sigma_gamma ~ independent Normal(0, sigma_gamma^2), all i,r.
##   omega_ijr | sigma_omega[i] ~ independent Normal(0, sigma_omega[i]^2).
##   sigma_gamma is ONE scalar, shared across all subjects and directions.
##   sigma_omega is a vector of length N; each subject shares its SD over r,j.
##   gamma_ir values remain different across subjects/directions.
## omega is analytically integrated out in the measurement likelihood.
## FUI-based SENSITIVITY variant: sigma_gamma has a bounded normal prior.
## Defaults: Normal(location=26.4, SD=2.64), restricted to [21.12,31.68].
## 26.4 is a proxy from reconstructed FUI curves in the SAME Bx coordinates.
## The +/-20% bounds are chosen diagnostic bounds, NOT a FUI confidence interval.
## FUI estimates from these same data do not supply independent prior evidence.
## A bound on the population SD does not bound gamma_ir or X_i(t) themselves.
## Omega/e priors remain HN(15)/HN(10); the likelihood and sharing are retained.
## Sharing is imposed in the retained UNIT basis; basis/priors from a paper
## in another coordinate scaling are not automatically numerically equivalent.
## Source this file once, then call fit_flfoqr_fui_bounded_gamma_betaf_xcov(...).
## Source only defines functions; the example at the bottom is disabled.
## beta(t) retains the explicit mgcv penalty eigen-reparameterization.
## ============================================================
## Helper functions
## ============================================================

#Transform the integral into a weighted summation
trapezoid_weights <- function(grid) {    
  # 1. Ensure the input is treated as a numeric vector (handles integers/factors
  grid <- as.numeric(grid)
  # 2. Validate grid size: At least two points are needed to define an interval
  if (length(grid) < 2) stop("grid must have at least two points")
  # 3. Validate monotonicity: diff(grid) <= 0 catches decreasing or duplicate points
  if (any(diff(grid) <= 0)) stop("grid must be strictly increasing")
  # 4. Pre-allocate a numeric vector of zeros to store the weights
  w <- numeric(length(grid))
  
  #    w_1 = (x_2 - x_1) / 2
  w[1] <- (grid[2] - grid[1]) / 2
  
  #    w_n = (x_n - x_(n-1)) / 2
  w[length(grid)] <- (grid[length(grid)] - grid[length(grid) - 1]) / 2
  
  # 7. Interior nodes: If grid length > 2, compute weights for internal points
  if (length(grid) > 2) {
    w[2:(length(grid) - 1)] <-
      (grid[3:length(grid)] - grid[1:(length(grid) - 2)]) / 2
  }
  w
}

## Construct a second-order difference penalty.
##
## D2 applies second differences to the original B-spline coefficients.
## Therefore,
##
##   P = D2' D2
##
## penalizes curvature: constant and approximately linear coefficient
## patterns form the unpenalized null space, whereas increasingly
## oscillatory patterns receive larger penalties.
second_diff_penalty <- function(K, ridge = 0) {
  if (K <= 2) return(diag(K) * max(ridge, 1))
  ## Second-difference operator on K spline coefficients.
  D2 <- diff(diag(K), differences = 2)
  ## Quadratic penalty matrix. In the current basis construction,
  ## ridge = 0 here so that the original penalty structure is retained.
  crossprod(D2) + ridge * diag(K)
}

## Convert an arbitrary spline basis B0 into a weighted-orthogonal,
## penalty-ordered basis B.
##
## Given the original basis B0, penalty P, and weight matrix W,
## the function eigendecomposes
##
##   Cw = W^(1/2) B0 (P + ridge I)^(-1) B0' W^(1/2).
##
## The eigenvectors define functional directions ordered from smooth
## to rough. The returned basis satisfies
##
##   B' W B = diag(d).
##
## Hence the columns are orthogonal, although they are not normalized
## to have unit norm. The projection operator is
##
##   Proj = diag(1/d) B' W,
##
## so that Proj %*% f gives the coefficients of a discretized curve f
## in the B basis.

orthogonalize_basis_any <- function(
    B0, #basis matrix
    P = NULL, #penalty matrix
    grid = NULL, 
    grid_weights = NULL, #Integral weight
    
    ridge = 1e-8, 
    ## Added to P before inversion because a second-difference penalty
    ## is singular. Ridge also determines the scale of null-space directions.
    
    tol = 1e-10,  # The threshold for whether eigenvalue is zero.
    max_rank = NULL #number of basis directions that can be retained
) {
  B0 <- as.matrix(B0)
  if (!is.numeric(B0)) stop("B0 must be numeric")
  
  # check Dimensions
  Tgrid <- nrow(B0)
  K0 <- ncol(B0)
  
  if (is.null(grid_weights)) {
    grid_weights <-
      if (!is.null(grid)) trapezoid_weights(grid) else rep(1, Tgrid)
  }
  
  grid_weights <- as.numeric(grid_weights)
  
  if (length(grid_weights) != Tgrid) {
    stop("grid_weights must have length nrow(B0)")
  }
  
  if (any(grid_weights <= 0)) {
    stop("grid_weights must be positive")
  }
  
  if (is.null(max_rank)) max_rank <- K0
  max_rank <- min(max_rank, K0, Tgrid)
  
  #Weighted space transformation
  sqrt_w <- sqrt(grid_weights)
  inv_sqrt_w <- 1 / sqrt_w
  
  if (is.null(P)) {
    
    #no penalty: weighted QR branch
    qrobj <- qr(B0 * sqrt_w, tol = tol)
    rnk <- min(qrobj$rank, max_rank)
    
    Q <- qr.Q(
      qrobj,
      complete = FALSE
    )[, seq_len(rnk), drop = FALSE]
    
    B <- Q * inv_sqrt_w # weighted-orthonormal basis
    d <- rep(1, rnk)
    
  } else {
    
    ##penalty-based eigen branch
    P <- as.matrix(P)
    
    if (!all(dim(P) == c(K0, K0))) {
      stop("P must be K x K where K = ncol(B0)")
    }
    
    P <- 0.5 * (P + t(P))
    
    #Regularize the penalty and find its inverse.
    P_inv <- solve(
      P + ridge * diag(K0)
    )
    
    #Constructing a weighted kernel/operator
    ##TÃ—T symmetric positive semidefinite matrix
    ## 1. Original basis space: B0
    ## 2. smoothness penaltyï¼šPï¼›
    ## 3. functional inner productï¼šW
    Cw <-
      (B0 * sqrt_w) %*%
      P_inv %*%
      t(B0 * sqrt_w)
    
    Cw <- 0.5 * (Cw + t(Cw))
    
    #CW=UÎ›UâŠ¤
    ## Eigenvalues are returned from largest to smallest.
    ## Large d_r corresponds to a smoother, weakly penalized direction.
    ## Small d_r corresponds to a rougher, more strongly penalized direction.
    eig <- eigen(
      Cw,
      symmetric = TRUE
    )
    
    #How many directions to retain
    keep <- which(eig$values > tol)
    
    if (length(keep) == 0) {
      stop("No positive eigenvalues retained; check B0/P/tol")
    }
    
    keep <-
      keep[
        seq_len(
          min(
            length(keep),
            max_rank
          )
        )
      ]
    
    d <- as.numeric(eig$values[keep])
    
    U <-
      eig$vectors[
        ,
        keep,
        drop = FALSE
      ]
    
    ## Construct the final basis
    B <-
      (U * rep(
        inv_sqrt_w,
        times = length(keep)
      )) %*%
      diag(
        sqrt(d),
        nrow = length(d)
      )
  }
  
  G <- crossprod(
    B,
    B * grid_weights
  )
  
  if (
    max(
      abs(
        G -
        diag(
          diag(G),
          nrow = ncol(G)
        )
      )
    ) > 1e-6
  ) {
    warning(
      "Orthogonalized basis check failed: ",
      "off-diagonal inner products not near zero"
    )
  }
  
  d <- as.numeric(diag(G))
  
  ## Convert the curve to a basis coefficient.
  Proj <-
    diag(
      1 / d,
      nrow = length(d)
    ) %*%
    t(B * grid_weights)
  
  ## Final output object
  out <- list(
    B = B, 
    ## Final weighted-orthogonal, penalty-ordered basis.
    
    B_orth = B,
    
    d = d,
    ## Weighted squared basis norms: diag(B' W B).
    
    Proj = Proj, 
    #projection matrix of the curve onto the basis coefficient
    
    grid_weights = grid_weights,
    
    Gram = G 
    ## Weighted Gram matrix B' W B, which should be diagonal.
  )
  
  class(out) <- "orthogonal_basis"
  
  out
}

## ============================================================
## Construct beta(t) exactly in the mgcv penalty parameterization
## used by the FPCA scalar-on-function model.
##
## Raw cyclic/cubic spline representation:
##
##   beta(t) = Psi(t)' b.
##
## If S = U diag(lambda) U' is the mgcv smoothing penalty, define
##
##   B_beta(t) = Psi(t) U E^{-1}.
##
## The positive-penalty columns are random/shrunk directions and
## the penalty-null-space columns are fixed directions.  The returned
## B_beta and G_beta place fixed columns first, followed by random
## columns, so their coefficient vector is c(betaf, betar).
## ============================================================

build_mgcv_beta_basis <- function(
    time_grid,
    quad_w,
    Bx,
    k = 10,
    beta_bs = c("cc", "cr")
) {
  beta_bs <- match.arg(beta_bs)
  time_grid <- as.numeric(time_grid)
  quad_w <- as.numeric(quad_w)
  Bx <- as.matrix(Bx)
  k <- as.integer(k)
  
  if (length(time_grid) < 4L) {
    stop("time_grid must contain at least four points.")
  }
  
  if (length(quad_w) != length(time_grid)) {
    stop("quad_w must have the same length as time_grid.")
  }
  
  if (any(!is.finite(quad_w)) || any(quad_w <= 0)) {
    stop("quad_w must contain finite positive weights.")
  }
  
  if (nrow(Bx) != length(time_grid)) {
    stop("nrow(Bx) must equal length(time_grid).")
  }
  
  if (k < 4L) {
    stop("k must be at least 4 for an mgcv cubic spline basis.")
  }
  
  ## Construct the raw time-domain beta basis and its mgcv penalty.
  ## Supplying the full [0,1] boundary for a cyclic spline preserves
  ## the intended 24-hour period even when the last observed grid point
  ## is (T-1)/T rather than 1.
  t_beta <- time_grid
  beta_spec <- mgcv::s(
    t_beta,
    bs = beta_bs,
    k = k
  )
  
  beta_knots <-
    if (identical(beta_bs, "cc")) {
      list(t_beta = c(0, 1))
    } else {
      NULL
    }
  
  splinecons <- mgcv::smooth.construct(
    beta_spec,
    data = list(t_beta = t_beta),
    knots = beta_knots
  )
  
  Psi_mat <- as.matrix(splinecons$X)
  S_mat_raw <- as.matrix(splinecons$S[[1]])
  rank_beta <- as.integer(splinecons$rank)
  K_raw <- ncol(Psi_mat)
  
  if (rank_beta < 1L || rank_beta >= K_raw) {
    stop(
      "The beta penalty must have both a positive-penalty space ",
      "and a non-empty null space."
    )
  }
  
  ## Match the penalty scaling used by mgcv::smoothCon before the
  ## natural-parameter transformation.
  maXX <- norm(Psi_mat, type = "I")^2
  maS <- norm(S_mat_raw, type = "I") / maXX
  
  if (!is.finite(maS) || maS <= 0) {
    stop("The mgcv beta penalty has a non-finite or zero scale.")
  }
  
  S_mat <- S_mat_raw / maS
  eig_beta <- eigen(S_mat, symmetric = TRUE)
  U_beta <- eig_beta$vectors
  
  ## First rank_beta eigenvalues are the positive penalty eigenvalues;
  ## eigen() returns them in decreasing order for a symmetric matrix.
  ## Keep the same positive-eigenvalue floor as the FPCA reference code.
  eigen_floor <- 1e-12
  
  E_beta <- rep(1, K_raw)
  E_beta[seq_len(rank_beta)] <-
    sqrt(
      pmax(
        eig_beta$values[seq_len(rank_beta)],
        eigen_floor
      )
    )
  
  ## Raw X-beta linkage before penalty reparameterization:
  ##
  ##   G0 = integral Bx(t) Psi(t)' dt.
  G0_beta <- crossprod(
    Bx,
    sweep(
      Psi_mat,
      1,
      quad_w,
      "*"
    )
  )
  
  B_rot <- Psi_mat %*% U_beta
  G_rot <- G0_beta %*% U_beta
  
  ## This is the E construction used in the FPCA implementation:
  ## scale the penalty-null-space columns to the average design-column
  ## norm of the positive-penalty space.
  col_norm <- colSums(G_rot^2) / E_beta^2
  av_norm <- mean(col_norm[seq_len(rank_beta)])
  
  if (!is.finite(av_norm) || av_norm <= 0) {
    stop("The penalized beta directions have zero linkage to Bx.")
  }
  
  fixed_idx <- seq.int(rank_beta + 1L, K_raw)
  random_idx <- seq_len(rank_beta)
  
  fixed_norm <- col_norm[fixed_idx]
  
  if (any(!is.finite(fixed_norm)) || any(fixed_norm <= 0)) {
    stop("A beta null-space direction has zero linkage to Bx.")
  }
  
  E_beta[fixed_idx] <- sqrt(fixed_norm / av_norm)
  
  ## Apply the identical U E^{-1} map to the time-domain basis and
  ## to the X-beta linkage.  This identity is essential:
  ##
  ##   crossprod(Bx, Wq B_beta) = G_beta.
  B_penalty_order <- sweep(B_rot, 2, E_beta, "/")
  G_penalty_order <- sweep(G_rot, 2, E_beta, "/")
  
  B_beta_f <- B_penalty_order[, fixed_idx, drop = FALSE]
  B_beta_r <- B_penalty_order[, random_idx, drop = FALSE]
  G_beta_f <- G_penalty_order[, fixed_idx, drop = FALSE]
  G_beta_r <- G_penalty_order[, random_idx, drop = FALSE]
  
  ## Fixed directions are stored first to match c(betaf, betar).
  B_beta <- cbind(B_beta_f, B_beta_r)
  G_beta <- cbind(G_beta_f, G_beta_r)
  
  G_beta_direct <- crossprod(
    Bx,
    sweep(
      B_beta,
      1,
      quad_w,
      "*"
    )
  )
  
  linkage_error <- max(abs(G_beta_direct - G_beta))
  linkage_tol <- 1e-10 * max(1, max(abs(G_beta)))
  
  if (!is.finite(linkage_error) || linkage_error > linkage_tol) {
    stop("The transformed beta basis and G_beta linkage are inconsistent.")
  }
  
  Gram_beta <- crossprod(
    B_beta,
    sweep(
      B_beta,
      1,
      quad_w,
      "*"
    )
  )
  
  list(
    beta_bs = beta_bs,
    k_requested = k,
    splinecons = splinecons,
    Psi_mat = Psi_mat,
    S_mat_raw = S_mat_raw,
    S_mat = S_mat,
    penalty_scale = maS,
    rank_beta = rank_beta,
    U_beta = U_beta,
    E_beta = E_beta,
    fixed_idx_penalty_order = fixed_idx,
    random_idx_penalty_order = random_idx,
    B_penalty_order = B_penalty_order,
    G_penalty_order = G_penalty_order,
    B_beta_f = B_beta_f,
    B_beta_r = B_beta_r,
    B_beta = B_beta,
    G_beta_f = G_beta_f,
    G_beta_r = G_beta_r,
    G_beta = G_beta,
    K_beta_f = ncol(B_beta_f),
    K_beta_r = ncol(B_beta_r),
    K_beta = ncol(B_beta),
    Gram_beta = Gram_beta,
    d_beta = diag(Gram_beta),
    linkage_error = linkage_error
  )
}

## ============================================================
## Post-hoc reconstruction of beta(t), ORIGINAL Y scale
##
## V9 IMPORTANT:
## beta(t) is reconstructed with the transformed mgcv B_beta whose
## columns are ordered c(fixed, random), exactly like b_beta in Stan.
## ============================================================

flfoqr_beta_curve <- function(
    fit,
    probs = c(0.025, 0.5, 0.975)
) {
  
  B_beta <- attr(
    fit,
    "B_beta"
  )
  
  Y_sd <- attr(
    fit,
    "Y_sd"
  )
  
  tg <- attr(
    fit,
    "time_grid"
  )
  
  if (is.null(B_beta)) {
    stop(
      "B_beta is missing. ",
      "This beta reconstruction requires the explicit mgcv-beta model."
    )
  }
  
  if (is.null(Y_sd)) {
    stop("Y_sd is missing from fit attributes.")
  }
  
  draws <-
    rstan::extract(
      fit,
      pars = "b_beta",
      permuted = TRUE
    )$b_beta
  
  draws <- as.matrix(draws)
  
  if (
    ncol(draws) !=
    ncol(B_beta)
  ) {
    stop(
      "Dimension mismatch: ncol(b_beta draws) = ",
      ncol(draws),
      " but ncol(B_beta) = ",
      ncol(B_beta)
    )
  }
  
  ## posterior draws x time
  beta_draws <-
    draws %*%
    t(B_beta)
  
  ## Back-transform from standardized Y to original Y scale.
  beta_draws <-
    beta_draws *
    Y_sd
  
  qs <-
    apply(
      beta_draws,
      2,
      quantile,
      probs = probs,
      na.rm = TRUE
    )
  
  data.frame(
    t = tg,
    mean = colMeans(
      beta_draws,
      na.rm = TRUE
    ),
    lower = qs[1, ],
    median = qs[2, ],
    upper = qs[3, ]
  )
}

## ============================================================
## Post-hoc reconstruction of a subject's latent curve.
## gamma is stored as a DEVIATION, so a_x must be added back.
## ============================================================

## Extract posterior mean latent deviations without relying on printed order.
## The name is retained for compatibility with V6/V7 analysis scripts.
## rstan::extract() returns draws x N x K_func.
flfoqr_xi_mean <- function(
    fit
) {
  latent_name <-
    if (identical(attr(fit, "outcome_xi_mode"), "joint") &&
        grepl("v8|v9", attr(fit, "model_version"))) {
      "gamma"
    } else {
      "xi"
    }
  
  xi_draws <-
    rstan::extract(
      fit,
      pars = latent_name,
      permuted = TRUE
    )[[latent_name]]
  
  if (
    length(
      dim(
        xi_draws
      )
    ) != 3L
  ) {
    stop(
      "The extracted latent-score object must have dimensions draws x N x K_func."
    )
  }
  
  apply(
    xi_draws,
    c(2, 3),
    mean
  )
}

flfoqr_latent_curve <- function(
    fit,
    i
) {

  Bx <- attr(fit, "Bx")
  X_cov <- attr(fit, "X_cov")

  if (is.null(Bx)) stop("Bx is missing from fit attributes.")
  if (is.null(X_cov)) stop("X_cov is missing; use the xcov model fit.")
  if (i < 1L || i > nrow(X_cov)) stop("i is outside the retained-subject range.")

  a_x_hat <- rstan::summary(fit, pars = "a_x")$summary[, "mean"]
  gamma_hat <- flfoqr_xi_mean(fit)

  cov_score_i <- rep(0, ncol(Bx))
  if (ncol(X_cov) > 0L) {
    alpha_draws <- rstan::extract(fit, pars = "alpha_x", permuted = TRUE)$alpha_x
    alpha_hat <- apply(alpha_draws, c(2, 3), mean)
    cov_score_i <- as.numeric(X_cov[i, , drop = FALSE] %*% alpha_hat)
  }

  as.numeric(
    Bx %*% (a_x_hat + cov_score_i + gamma_hat[i, ])
  )
}

## Reconstruct alpha_q(t): the time-varying association between a scalar
## latent-PA covariate and the true PA trajectory.  Numeric x covariates were
## standardized, so their curves are per 1-SD increase.  Dummy-factor effects
## retain their one-unit group contrast despite centering of the design column.
flfoqr_xcov_effect_curve <- function(
    fit,
    covariate,
    probs = c(0.025, 0.5, 0.975)
) {
  Bx <- attr(fit, "Bx")
  tg <- attr(fit, "time_grid")
  nm <- attr(fit, "x_covariate_names")

  if (is.null(Bx) || is.null(nm)) stop("Latent-PA covariate attributes are missing.")
  if (length(nm) == 0L) stop("This fit has no measurement-model covariates.")

  q <- if (is.character(covariate)) match(covariate, nm) else as.integer(covariate)
  if (length(q) != 1L || is.na(q) || q < 1L || q > length(nm)) {
    stop("covariate must be one valid x_covariate_names entry or column index.")
  }

  alpha_draws <- rstan::extract(fit, pars = "alpha_x", permuted = TRUE)$alpha_x
  a_q <- alpha_draws[, q, , drop = FALSE]
  dim(a_q) <- c(dim(alpha_draws)[1], dim(alpha_draws)[3])
  curve_draws <- a_q %*% t(Bx)

  qs <- apply(curve_draws, 2, quantile, probs = probs, na.rm = TRUE)
  data.frame(
    covariate = nm[q],
    t = tg,
    mean = colMeans(curve_draws),
    lower = qs[1, ],
    median = qs[2, ],
    upper = qs[3, ]
  )
}

## ============================================================
## Main function
## ============================================================

fit_flfoqr_fui_bounded_gamma_betaf_xcov <- function(
    Y,
    W,
    Z,
    a = NULL,
    Nsim,
    tau0 = 0.5,
    numeric_Z = NULL,
    factor_Z = NULL,

    ## ---- NEW: scalar covariates for latent physical activity X_i(t) ----
    ## Example: x_numeric_Z = "AgeYR", x_factor_Z = "Gender".
    ## Set both to NULL (the default) for no measurement-model covariates.
    ## Continuous x covariates are standardized; the final design columns are
    ## centered so a_x remains the mean latent-PA score vector.
    x_numeric_Z = NULL,
    x_factor_Z = NULL,
    x_covariate_prior_sd = 10,

    ## FALSE: outcome functional term uses residual PA gamma only.
    ## TRUE:  outcome functional term uses X-covariate contribution + gamma.
    ## This switch lets you fit both sensitivity models discussed with Dr. Zoh.
    x_covariate_in_outcome = FALSE,
    
    ## ---- measurement basis ----
    K_func = 12,
    basis_df = NULL,
    basis_degree = 3,
    basis_ridge = 0.01,
    basis_tol = 1e-4,
    basis_weights = c(
      "unit",
      "trapezoid"
    ),
    
    ## ---- explicit mgcv beta(t) basis ----
    ## K_beta is the mgcv k argument.  For bs="cc", k=10 normally
    ## produces 9 coefficients: 1 fixed null-space direction and
    ## 8 penalized random directions.
    K_beta = 10,
    beta_bs = c("cc", "cr"),
    
    periodic_time = FALSE,
    
    ## ---- prior scales ----
    a_x_prior_sd = 100,
    sigma_gamma_prior_mean = 26.4,
    sigma_gamma_prior_sd = 2.64,
    sigma_gamma_lower = 21.12,
    sigma_gamma_upper = 31.68,
    sigma_omega_prior_scale = 15,
    sigma_e_prior_scale = 10,
    beta_fixed_prior_sd = 0.5,
    beta_random_prior_scale = 0.5,
    chains = 2,
    cores = 2,
    seed = 1123,
    warmup = NULL,
    adapt_delta = 0.95,
    max_treedepth = 12,
    beta_fixed_prior_mean = 0,
    scalar_formula = NULL,
    fast_likelihood = TRUE
) {
  if (!is.logical(fast_likelihood) || length(fast_likelihood) != 1L ||
      is.na(fast_likelihood)) {
    stop("fast_likelihood must be TRUE or FALSE.")
  }

  if (!is.logical(x_covariate_in_outcome) ||
      length(x_covariate_in_outcome) != 1L ||
      is.na(x_covariate_in_outcome)) {
    stop("x_covariate_in_outcome must be TRUE or FALSE.")
  }

  if (!is.numeric(x_covariate_prior_sd) ||
      length(x_covariate_prior_sd) != 1L ||
      !is.finite(x_covariate_prior_sd) ||
      x_covariate_prior_sd <= 0) {
    stop("x_covariate_prior_sd must be one finite positive number.")
  }
  
  ################################################################################
  #                                                                              #
  # Function input and parameter checking                                                      
  #                                                                              #
  ################################################################################
  
  cat(
    "\n#--- Bayesian Quantile FLFOQR: FUI-bounded shared gamma / subject omega ---\n"
  )
  
  cat(
    "#--- weighted-orthonormal X basis; ONE sigma_gamma shared over all i and r ---\n"
  )
  
  cat(
    "#--- beta uses mgcv spline/penalty eigen-reparameterization ---\n"
  )
  
  cat(
    "#--- omega = sigma_omega_i (shared over r); ",
    "sigma_e shared in-span and off-span ---\n\n"
  )
  
  require(rstan)
  require(brms)
  require(splines)
  require(mgcv)
  require(MASS)
  
  rstan_options(
    auto_write = TRUE
  )
  
  options(
    mc.cores = cores
  )
  
  basis_weights <-
    match.arg(
      basis_weights
    )
  
  beta_bs <-
    match.arg(
      beta_bs
    )
  
  if (is.null(warmup)) {
    warmup <-
      floor(
        Nsim / 2
      )
  }
  
  if (
    warmup < 1L ||
    warmup >= Nsim
  ) {
    stop(
      "warmup must be >=1 and < Nsim."
    )
  }
  
  if (
    !is.finite(adapt_delta) ||
    adapt_delta <= 0 ||
    adapt_delta >= 1
  ) {
    stop(
      "adapt_delta must lie strictly between 0 and 1."
    )
  }
  
  if (max_treedepth < 1L) {
    stop(
      "max_treedepth must be a positive integer."
    )
  }
  
  gamma_args <- list(
    sigma_gamma_prior_mean = sigma_gamma_prior_mean,
    sigma_gamma_prior_sd = sigma_gamma_prior_sd,
    sigma_gamma_lower = sigma_gamma_lower,
    sigma_gamma_upper = sigma_gamma_upper
  )
  if (!all(vapply(gamma_args, function(x) {
    is.numeric(x) && length(x) == 1L && is.finite(x)
  }, logical(1)))) {
    stop("Gamma prior arguments must each be one finite number.")
  }
  if (!(sigma_gamma_lower >= 0 &&
        sigma_gamma_lower < sigma_gamma_prior_mean &&
        sigma_gamma_prior_mean < sigma_gamma_upper &&
        sigma_gamma_prior_sd > 0)) {
    stop("Require 0 <= gamma lower < normal location < gamma upper, and prior SD > 0.")
  }
  cat("\nBounded gamma SD prior: Normal(location=", sigma_gamma_prior_mean,
      ", SD=", sigma_gamma_prior_sd, "), restricted to [", sigma_gamma_lower,
      ", ", sigma_gamma_upper, "].\n", sep = "")
  cat("Bounds apply to the shared SD, not to individual gamma coefficients.\n")
  
  if (!is.numeric(beta_fixed_prior_mean) ||
      !is.null(dim(beta_fixed_prior_mean)) ||
      length(beta_fixed_prior_mean) < 1L ||
      any(!is.finite(beta_fixed_prior_mean))) {
    stop("beta_fixed_prior_mean must be a finite numeric scalar or vector.")
  }
  
  prior_scales <-
    c(
      a_x_prior_sd,
      sigma_gamma_prior_sd,
      sigma_omega_prior_scale,
      sigma_e_prior_scale,
      beta_fixed_prior_sd,
      beta_random_prior_scale
    )
  
  if (
    any(!is.finite(prior_scales)) ||
    any(prior_scales <= 0)
  ) {
    stop(
      "All prior scale arguments must be finite and strictly positive."
    )
  }
  
  K_func <-
    as.integer(
      K_func
    )
  
  if (
    K_func < 4L
  ) {
    stop(
      "K_func must be >= 4 for a cubic B-spline basis."
    )
  }
  
  if (
    is.null(
      basis_df
    )
  ) {
    basis_df <-
      K_func
  }
  
  basis_df <-
    as.integer(
      basis_df
    )
  
  if (
    basis_df <
    K_func
  ) {
    stop(
      "basis_df must be at least K_func."
    )
  }
  
  K_beta <-
    as.integer(
      K_beta
    )
  
  if (
    K_beta < 4L
  ) {
    stop(
      "K_beta must be >= 4 for an mgcv cubic beta basis."
    )
  }
  
  ################################################################################
  #                                                                              #
  # GAL helper functions                                                         #
  #                                                                              #
  ################################################################################
  
  ## ============================================================
  ## 1. GAL helper functions
  ## ============================================================
  
  GamF <- function(
    gam,
    p0
  ) {
    (
      2 *
        pnorm(
          -abs(gam)
        ) *
        exp(
          0.5 *
            gam^2
        ) -
        p0
    )^2
  }
  
  GamBnd <- function(
    p0
  ) {
    
    Re1 <-
      optimize(
        GamF,
        interval = c(
          -30,
          30
        ),
        p0 = 1 - p0
      )
    
    Re2 <-
      optimize(
        GamF,
        interval = c(
          -30,
          30
        ),
        p0 = p0
      )
    
    c(
      -abs(
        Re1$minimum
      ),
      abs(
        Re2$minimum
      ),
      Re1$objective,
      Re2$objective
    )
  }
  
  Bd <-
    round(
      GamBnd(
        tau0
      )[1:2] *
        0.99,
      4
    )
  
  ## ============================================================
  ## 2. Input checks and complete-case filtering
  ## ============================================================
  
  Y <-
    as.numeric(
      Y
    )
  
  Z_df <-
    as.data.frame(
      Z
    )
  
  if (
    is.null(
      numeric_Z
    )
  ) {
    numeric_Z <-
      character(0)
  }
  
  if (
    is.null(
      factor_Z
    )
  ) {
    factor_Z <-
      character(0)
  }
  
  if (
    length(
      factor_Z
    ) > 0L
  ) {
    
    for (
      nm in factor_Z
    ) {
      
      Z_df[[nm]] <-
        as.factor(
          Z_df[[nm]]
        )
    }
  }

  ## NEW: covariates used to explain the latent physical-activity process.
  if (is.null(x_numeric_Z)) x_numeric_Z <- character(0)
  if (is.null(x_factor_Z))  x_factor_Z  <- character(0)

  x_covariate_vars <- unique(c(x_numeric_Z, x_factor_Z))
  has_x_covariates <- length(x_covariate_vars) > 0L

  missing_xcov <- setdiff(x_covariate_vars, colnames(Z_df))
  if (length(missing_xcov) > 0L) {
    stop("Latent-PA covariates not found in Z: ", paste(missing_xcov, collapse = ", "))
  }

  if (length(intersect(x_numeric_Z, x_factor_Z)) > 0L) {
    stop("A latent-PA covariate cannot be listed in both x_numeric_Z and x_factor_Z.")
  }

  if (length(x_factor_Z) > 0L) {
    for (nm in x_factor_Z) Z_df[[nm]] <- as.factor(Z_df[[nm]])
  }
  
  if (
    length(Y) !=
    nrow(Z_df)
  ) {
    stop(
      "length(Y) and nrow(Z) are not equal."
    )
  }
  
  if (
    length(
      dim(W)
    ) == 3L
  ) {
    
    if (
      dim(W)[1] !=
      length(Y)
    ) {
      stop(
        "dim(W)[1] and length(Y) are not equal."
      )
    }
    
  } else if (
    length(
      dim(W)
    ) == 2L
  ) {
    
    if (
      nrow(W) !=
      length(Y)
    ) {
      stop(
        "nrow(W) and length(Y) are not equal."
      )
    }
    
  } else {
    
    stop(
      "W must be either N x Tn or N x Tn x J."
    )
  }
  
  idx <-
    complete.cases(
      Y,
      Z_df
    ) &
    apply(
      W,
      1,
      function(x) {
        all(
          is.finite(x)
        )
      }
    )
  
  Y_use <-
    Y[
      idx
    ]
  
  Z_use <-
    Z_df[
      idx,
      ,
      drop = FALSE
    ]
  
  W_use <-
    if (
      length(
        dim(W)
      ) == 3L
    ) {
      
      W[
        idx,
        ,
        ,
        drop = FALSE
      ]
      
    } else {
      
      W[
        idx,
        ,
        drop = FALSE
      ]
    }
  
  cat(
    "After filtering:  N =",
    length(Y_use),
    " Z:",
    paste(
      dim(Z_use),
      collapse = " x "
    ),
    " W:",
    paste(
      dim(W_use),
      collapse = " x "
    ),
    "\n\n"
  )
  
  ## ============================================================
  ## 3. Standardize outcome and scalar covariates
  ## ============================================================
  
  Y_mean <-
    mean(
      Y_use,
      na.rm = TRUE
    )
  
  Y_sd <-
    sd(
      Y_use,
      na.rm = TRUE
    )
  
  if (
    is.finite(Y_sd) &&
    Y_sd > 1e-9
  ) {
    
    Y_std <-
      as.numeric(
        (
          Y_use -
            Y_mean
        ) /
          Y_sd
      )
    
  } else {
    
    Y_std <-
      as.numeric(
        Y_use -
          Y_mean
      )
    
    Y_sd <-
      1
  }
  
  Z_numeric <-
    Z_use[
      ,
      numeric_Z,
      drop = FALSE
    ]
  
  Z_factor <-
    Z_use[
      ,
      factor_Z,
      drop = FALSE
    ]
  
  if (
    length(
      factor_Z
    ) > 0L
  ) {
    
    for (
      nm in factor_Z
    ) {
      
      Z_factor[[nm]] <-
        as.factor(
          Z_factor[[nm]]
        )
    }
  }
  
  Z_final_model <-
    if (
      ncol(
        Z_numeric
      ) > 0L
    ) {
      
      data.frame(
        as.data.frame(
          scale(
            Z_numeric,
            center = TRUE,
            scale = TRUE
          )
        ),
        Z_factor
      )
      
    } else {
      
      Z_factor
    }

  ## ============================================================
  ## 3b. NEW: latent-PA covariate design matrix
  ##
  ## In the Bx score coordinates the latent process is
  ##   x_score_i = a_x + X_cov_i * alpha_x + gamma_i.
  ##
  ## Numeric covariates are standardized.  Then every design column is
  ## centered.  Centering keeps a_x interpretable as the sample-average
  ## latent PA score vector and reduces posterior dependence with alpha_x.
  ## ============================================================

  Xcov_df <- Z_use[, x_covariate_vars, drop = FALSE]

  x_numeric_center <- numeric(0)
  x_numeric_scale <- numeric(0)

  if (length(x_numeric_Z) > 0L) {
    x_numeric_center <- vapply(
      Xcov_df[, x_numeric_Z, drop = FALSE], mean, numeric(1), na.rm = TRUE
    )
    x_numeric_scale <- vapply(
      Xcov_df[, x_numeric_Z, drop = FALSE], sd, numeric(1), na.rm = TRUE
    )
    if (any(!is.finite(x_numeric_scale)) || any(x_numeric_scale <= 1e-9)) {
      stop("Every x_numeric_Z variable must have a finite, non-zero SD.")
    }
    for (nm in x_numeric_Z) {
      Xcov_df[[nm]] <-
        (as.numeric(Xcov_df[[nm]]) - x_numeric_center[[nm]]) / x_numeric_scale[[nm]]
    }
  }

  if (length(x_factor_Z) > 0L) {
    for (nm in x_factor_Z) Xcov_df[[nm]] <- as.factor(Xcov_df[[nm]])
  }

  if (has_x_covariates) {
    xcov_formula <- stats::reformulate(x_covariate_vars)
    X_cov_raw <- stats::model.matrix(xcov_formula, data = Xcov_df)
    X_cov_raw <- X_cov_raw[, colnames(X_cov_raw) != "(Intercept)", drop = FALSE]

    if (ncol(X_cov_raw) < 1L) {
      stop("The latent-PA covariate design matrix has no non-intercept columns.")
    }

    x_covariate_design_center <- colMeans(X_cov_raw)
    X_cov <- sweep(X_cov_raw, 2, x_covariate_design_center, "-")
    storage.mode(X_cov) <- "double"
    K_xcov <- as.integer(ncol(X_cov))
    x_covariate_names <- colnames(X_cov)
  } else {
    X_cov <- matrix(0, nrow = nrow(Z_use), ncol = 0L)
    K_xcov <- 0L
    x_covariate_names <- character(0)
    x_covariate_design_center <- numeric(0)
  }

  cat("\nLatent-PA covariate design columns (centered):\n")
  print(x_covariate_names)
  if (!has_x_covariates) {
    cat("No measurement-model covariates; x_score = a_x + gamma.\n")
  }
  cat("x_covariate_in_outcome =", x_covariate_in_outcome, "\n\n")
  
  ## ============================================================
  ## 4. Reshape repeated W to M x T
  ## ============================================================
  
  if (
    length(
      dim(W_use)
    ) == 3L
  ) {
    
    N <-
      dim(W_use)[1]
    
    Tn <-
      dim(W_use)[2]
    
    J <-
      dim(W_use)[3]
    
    M <-
      N *
      J
    
    W_obs <-
      matrix(
        NA_real_,
        nrow = M,
        ncol = Tn
      )
    
    ##Subject ID mapping
    id <-
      rep(
        seq_len(N),
        each = J
      )
    
    rr <-
      1L
    
    for (
      i in seq_len(N)
    ) {
      
      for (
        j in seq_len(J)
      ) {
        
        ## Write array as matrix
        W_obs[
          rr,
        ] <-
          W_use[
            i,
            ,
            j
          ]
        
        rr <-
          rr +
          1L
      }
    }
    
  } else {
    
    N <-
      nrow(
        W_use
      )
    
    Tn <-
      ncol(
        W_use
      )
    
    J <-
      1L
    
    M <-
      N
    
    W_obs <-
      as.matrix(
        W_use
      )
    
    id <-
      seq_len(N)
  }
  
  cat(
    "Repeated functional predictor:  N =",
    N,
    " Tn =",
    Tn,
    " J =",
    J,
    " M =",
    M,
    "\n\n"
  )
  
  if (
    J < 2L
  ) {
    stop(
      "At least two replicates per subject are required."
    )
  }
  
  J_i <-
    as.numeric(
      tabulate(
        id,
        nbins = N
      )
    )
  
  if (
    any(
      J_i < 2L
    )
  ) {
    stop(
      "Every retained subject must have at least two replicates."
    )
  }
  
  cat("Outcome xi mode: fully joint (gamma is sampled in both likelihoods)\n\n")
  
  ## ============================================================
  ## 5. Time grid and quadrature weights for outcome integral
  ## ============================================================
  
  if (
    isTRUE(
      periodic_time
    )
  ) {
    
    time_grid <-
      (
        0:
          (
            Tn -
              1L
          )
      ) /
      Tn
    
    quad_w <-
      rep(
        1 / Tn,
        Tn
      )
    
  } else {
    
    time_grid <-
      seq(
        0,
        1,
        length.out = Tn
      )
    
    quad_w <-
      rep(
        1,
        Tn
      )
    
    quad_w[1] <-
      0.5
    
    quad_w[Tn] <-
      0.5
    
    quad_w <-
      quad_w /
      (
        Tn -
          1
      )
  }
  
  W_mean_curve_emp <-
    colMeans(
      W_obs
    )
  
  cat(
    "Observed population mean curve range:",
    round(
      range(
        W_mean_curve_emp
      ),
      4
    ),
    "\n\n"
  )
  
  ## ============================================================
  ## 6. Construct the X basis
  ##
  ## The subject-specific latent deviation is represented as
  ##
  ##   gamma_i(t) = sum_r xi_ir Bx_r(t).
  ##
  ## Bx is constructed from a cubic B-spline basis and a
  ## second-difference smoothness penalty. Its directions are made
  ## orthogonal under the measurement inner product.
  ## ============================================================
  
  ## Evaluate Bx and B_beta on exactly the same physical grid.  This is
  ## essential when periodic_time=TRUE, because that grid ends at
  ## (Tn-1)/Tn rather than duplicating midnight at t=1.
  tau_x <-
    time_grid
  
  if (
    basis_weights ==
    "unit"
  ) {
    
    #Euclidean inner product
    ## Use the discrete Euclidean inner product for the measurement model:
    ##
    ##   <f,g>_measurement = sum_t f(t)g(t).
    ##
    ## This choice gives an exact Euclidean in-span/off-span
    ## decomposition of the observed curves.
    w_basis <-
      rep(
        1,
        Tn
      )
    
  } else {
    
    ## Alternative numerical integration weights.
    w_basis <-
      trapezoid_weights(
        tau_x
      )
    
    warning(
      "basis_weights='trapezoid': ",
      "off-span df is then only approximate."
    )
  }
  
  ################################################################################
  #                                                                              #
  # B-spline X basis (measurement spline space)                                                 
  #                                                                              #
  ################################################################################
  
  ## Create the original cubic B-spline basis.
  ## Rows = time points; columns = spline basis functions.
  ##   rows    = 1440 observed time points
  ##   columns = 20 original spline functions.
  
  B0_x <-
    as.matrix(
      splines::bs(
        tau_x,
        df = basis_df,
        degree = basis_degree,
        intercept = TRUE
      )
    )
  
  ## Create the second-difference penalty matrix.
  P_x <-
    second_diff_penalty(
      ncol(B0_x),
      ridge = 0
    )
  
  ##  ob_x: X basis Penalty orthogonalization
  ## Transform B0_x into penalty-ordered directions satisfying
  ##
  ##   Bx' W_basis Bx = diag(d_x).
  ##
  ## The ridge is introduced here when P_x is inverted.
  ob_x <-
    orthogonalize_basis_any(
      B0 = B0_x,
      P = P_x,
      grid = tau_x,
      grid_weights = w_basis,
      ridge = basis_ridge,
      tol = basis_tol,
      max_rank = K_func
    )
  
  ## Keep the original penalty-ordered directions for diagnostics, then
  ## normalize each direction.  The V9 coordinates satisfy Bx' W Bx = I.
  ## Sharing of sigma_gamma across subjects/directions is imposed in THESE
  ## unit coordinates; sigma_omega_i is shared over directions within subject.
  Bx_orthogonal <-
    ob_x$B
  
  d_x_orthogonal <-
    as.numeric(ob_x$d)
  
  Bx <- sweep(
    Bx_orthogonal,
    2,
    sqrt(d_x_orthogonal),
    "/"
  )
  
  K_func_eff <-
    ncol(
      Bx
    )
  
  if (
    K_func_eff <
    K_func
  ) {
    
    warning(
      "Only ",
      K_func_eff,
      " directions survived; K_func reduced from ",
      K_func,
      "."
    )
    
    K_func <-
      as.integer(
        K_func_eff
      )
  }
  
  ## Projection coefficients for a weighted-orthonormal basis are simply
  ## <W, Bx_r>.  The retained subspace is unchanged from V7.
  d_x <- rep(1, K_func)
  Proj_x <- t(Bx * w_basis)
  
  ## Confirm weighted orthonormality, not merely orthogonality.
  Gram_x <-
    crossprod(
      Bx,
      Bx * w_basis
    )
  
  max_unit_gram_error <-
    max(
      abs(
        Gram_x - diag(K_func)
      )
    )
  
  ## Largest off-diagonal inner product (reported for continuity).
  max_offdiag_x <-
    max(
      abs(
        Gram_x -
          diag(
            diag(Gram_x),
            nrow = K_func
          )
      )
    )
  
  ## Stop if the basis is not sufficiently orthogonal.
  if (
    max_unit_gram_error > 1e-8
  ) {
    stop(
      "Unit-normalized measurement basis failed Bx' W Bx = I."
    )
  }
  
  cat(
    "Measurement basis:  B0 =",
    paste(
      dim(B0_x),
      collapse = " x "
    ),
    "  Bx =",
    paste(
      dim(Bx),
      collapse = " x "
    ),
    "  weights =",
    basis_weights,
    "\n"
  )
  
  cat(
    "max error of Bx'WBx from I =",
    format(
      max_unit_gram_error,
      scientific = TRUE
    ),
    "\n"
  )
  
  cat(
    "original orthogonal d_x (before unit normalization) =",
    paste(
      signif(
        d_x_orthogonal,
        4
      ),
      collapse = ", "
    ),
    "\n\n"
  )
  
  ## ============================================================
  ## 7. Projection and off-span decomposition
  ## ============================================================
  
  ##RSS_off=RSS_totalâˆ’RSS_in-span.
  
  ## Project each observed daily curve into the Bx coordinates:
  ##
  ##   Wstar_mr = <W_m, Bx_r>, because Bx' W Bx = I.
  ##
  ## Wstar has M = N*J rows and K_func columns.
  Wstar <-
    W_obs %*%
    t(
      Proj_x
    )
  
  ## Algebraic regression check against the old non-unit coordinates:
  ## Wstar_unit = Wstar_orthogonal * sqrt(d_old).
  Wstar_orthogonal <-
    W_obs %*% t(ob_x$Proj)
  
  unit_score_transform_error <-
    max(
      abs(
        Wstar -
          sweep(Wstar_orthogonal, 2, sqrt(d_x_orthogonal), "*")
      )
    )
  
  if (unit_score_transform_error > 1e-8 * max(1, max(abs(Wstar)))) {
    stop("Unit-score coordinate transformation check failed.")
  }
  
  #Off-span decomposition
  rss_total <-
    sum(
      sweep(
        W_obs^2,
        2,
        w_basis,
        "*"
      )
    )
  
  rss_in_span <-
    sum(
      sweep(
        Wstar^2,
        2,
        d_x,
        "*"
      )
    )
  
  rss_off_raw <-
    rss_total -
    rss_in_span
  
  if (
    rss_off_raw <
    -1e-8 *
    max(
      1,
      rss_total
    )
  ) {
    stop(
      "Off-span RSS negative beyond tolerance."
    )
  }
  
  rss_off <-
    max(
      rss_off_raw,
      0
    )
  
  n_off <-
    as.integer(
      M *
        (
          Tn -
            K_func
        )
    )
  
  if (
    n_off < 1L
  ) {
    stop(
      "Tn must be strictly larger than K_func."
    )
  }
  
  offspan_fraction <-
    if (
      rss_total > 0
    ) {
      
      rss_off /
        rss_total
      
    } else {
      
      0
    }
  
  sigma_e_offspan_moment <-
    sqrt(
      rss_off /
        n_off
    )
  
  cat(
    "Wstar:",
    paste(
      dim(Wstar),
      collapse = " x "
    ),
    "\n"
  )
  
  cat(
    "per-direction SD =",
    paste(
      signif(
        apply(
          Wstar,
          2,
          sd
        ),
        3
      ),
      collapse = ", "
    ),
    "\n"
  )
  
  cat(
    "off-span fraction =",
    round(
      offspan_fraction,
      6
    ),
    "  off-span moment sigma_e =",
    round(
      sigma_e_offspan_moment,
      6
    ),
    "\n\n"
  )
  
  ## ============================================================
  ## 8. Construct the mgcv penalty-reparameterized beta(t) basis
  ##
  ## X keeps the existing measurement basis Bx.  Only beta changes:
  ##
  ##   beta(t) = B_beta_f(t) betaf + B_beta_r(t) betar,
  ##   betar   = sigma_beta z_betar.
  ##
  ## Both B_beta and its linkage to Bx receive the same U E^{-1}
  ## transformation from the FPCA scalar-on-function implementation.
  ## ============================================================
  
  beta_basis <-
    build_mgcv_beta_basis(
      time_grid = time_grid,
      quad_w = quad_w,
      Bx = Bx,
      k = K_beta,
      beta_bs = beta_bs
    )
  
  ## Raw mgcv time basis and penalty objects.
  Psi_mat <- beta_basis$Psi_mat
  S_mat_beta_raw <- beta_basis$S_mat_raw
  S_mat_beta <- beta_basis$S_mat
  U_beta <- beta_basis$U_beta
  E_beta <- beta_basis$E_beta
  rank_beta <- beta_basis$rank_beta
  
  ## Reparameterized time-domain beta bases.  Fixed/null-space columns
  ## are stored first; positive-penalty random columns follow.
  B_beta_f <- beta_basis$B_beta_f
  B_beta_r <- beta_basis$B_beta_r
  B_beta <- beta_basis$B_beta
  
  ## Matching linkage matrices:
  ##
  ##   G_beta_f = integral Bx(t) B_beta_f(t)' dt,
  ##   G_beta_r = integral Bx(t) B_beta_r(t)' dt.
  G_beta_f <- beta_basis$G_beta_f
  G_beta_r <- beta_basis$G_beta_r
  G_beta <- beta_basis$G_beta
  
  ## Direct old-coordinate linkage used only to verify that the coordinate
  ## change leaves every functional linear predictor invariant.
  G_beta_orthogonal_direct <-
    crossprod(
      Bx_orthogonal,
      sweep(B_beta, 1, quad_w, "*")
    )
  
  G_beta_orthogonal_from_unit <-
    sweep(G_beta, 1, sqrt(d_x_orthogonal), "*")
  
  unit_linkage_transform_error <-
    max(abs(G_beta_orthogonal_direct - G_beta_orthogonal_from_unit))
  
  if (
    unit_linkage_transform_error >
    1e-8 * max(1, max(abs(G_beta_orthogonal_direct)))
  ) {
    stop("Unit-coordinate G_beta transformation check failed.")
  }
  
  K_beta_f <- as.integer(beta_basis$K_beta_f)
  if (!(length(beta_fixed_prior_mean) %in% c(1L, K_beta_f))) {
    stop("beta_fixed_prior_mean must have length 1 or K_beta_f.")
  }
  beta_fixed_prior_mean <- rep(beta_fixed_prior_mean, length.out = K_beta_f)
  K_beta_r <- as.integer(beta_basis$K_beta_r)
  K_beta <- as.integer(beta_basis$K_beta)
  
  Gram_beta <- beta_basis$Gram_beta
  d_beta <- as.numeric(beta_basis$d_beta)
  beta_linkage_error <- beta_basis$linkage_error
  
  ## Coupling diagnostics in the transformed beta coordinates.
  G_row_norm <- sqrt(rowSums(G_beta^2))
  G_col_norm <- sqrt(colSums(G_beta^2))
  G_singular <- svd(G_beta, nu = 0, nv = 0)$d
  
  G_row_rel <-
    if (max(G_row_norm) > 0) {
      G_row_norm / max(G_row_norm)
    } else {
      rep(0, length(G_row_norm))
    }
  
  G_col_rel <-
    if (max(G_col_norm) > 0) {
      G_col_norm / max(G_col_norm)
    } else {
      rep(0, length(G_col_norm))
    }
  
  cat(
    "Explicit mgcv beta basis:",
    beta_bs,
    "with requested k =",
    beta_basis$k_requested,
    "\n"
  )
  
  cat(
    "Psi_mat =",
    paste(dim(Psi_mat), collapse = " x "),
    "  B_beta =",
    paste(dim(B_beta), collapse = " x "),
    "\n"
  )
  
  cat(
    "penalty rank =",
    rank_beta,
    "  K_beta_f =",
    K_beta_f,
    "  K_beta_r =",
    K_beta_r,
    "\n"
  )
  
  cat(
    "G_beta dimension =",
    paste(dim(G_beta), collapse = " x "),
    "  linkage error =",
    format(beta_linkage_error, scientific = TRUE),
    "\n"
  )
  
  cat(
    "G row norms relative to max:\n",
    paste(signif(G_row_rel, 4), collapse = ", "),
    "\n"
  )
  
  cat(
    "G column norms relative to max (fixed first, then random):\n",
    paste(signif(G_col_rel, 4), collapse = ", "),
    "\n\n"
  )
  
  ## ============================================================
  ## 9. Scalar outcome standata
  ## ============================================================
  
  df2 <- data.frame(
    Y = Y_std,
    Z_final_model
  )
  
  ## ============================================================
  ## Add interactions through the brms/model-matrix interface
  ## ============================================================
  
  ## All scalar variables currently available
  z_terms <- setdiff(
    colnames(df2),
    "Y"
  )
  
  ## These variables will be introduced through * below,
  ## so remove them from the ordinary main-effect list
  z_main_no_inter <- setdiff(
    z_terms,
    c(
      "Race",
      "Gender",
      "AgeYR",
      "HEI"
    )
  )
  
  ## Final model:
  ## main effects +
  ## Race x Gender +
  ## Race x AgeYR +
  ## Race x HEI
  form_Z <- paste0(
    "Y ~ ",
    paste(
      c(
        z_main_no_inter,
        "Race * Gender",
        "Race * AgeYR",
        "HEI * Race"
      ),
      collapse = " + "
    )
  )
  
  # An explicit scalar_formula overrides the legacy interaction formula above.
  # Numeric predictors have already been standardized in Z_final_model.
  # Accept Y ~ ... or a one-sided ~ ... without changing its RHS terms.
  # R formula reference:
  # https://stat.ethz.ch/R-manual/R-devel/library/stats/html/formula.html
  if (!is.null(scalar_formula)) {
    if (inherits(scalar_formula, "formula")) {
      form_Z <- scalar_formula
    } else if (is.character(scalar_formula) &&
               length(scalar_formula) == 1L && !is.na(scalar_formula)) {
      form_Z <- stats::as.formula(scalar_formula, env = parent.frame())
    } else {
      stop("scalar_formula must be a formula or one formula string.")
    }
    if (length(form_Z) == 2L) {
      form_Z <- stats::as.formula(
        call("~", as.name("Y"), form_Z[[2L]]),
        env = environment(form_Z)
      )
    }
    if (length(form_Z) != 3L || !identical(form_Z[[2L]], as.name("Y"))) {
      stop("Use response Y in scalar_formula, or supply a one-sided formula.")
    }
    scalar_terms <- stats::terms(form_Z, data = df2)
    if (attr(scalar_terms, "intercept") != 1L) {
      stop("This Stan model requires an intercept in scalar_formula.")
    }
    if (length(attr(scalar_terms, "offset")) > 0L) {
      stop("Offsets are not implemented in this Stan outcome predictor.")
    }
  }

  cat("\nFinal scalar outcome formula:\n")
  print(form_Z)
  
  Dt2 <- brms::make_standata(
    brms::bf(
      as.formula(form_Z),
      quantile = tau0
    ),
    data = df2,
    family = brms::asym_laplace()
  )
  
  cat("\nScalar design-matrix columns:\n")
  print(colnames(Dt2$X))
  ## ============================================================
  ## 10. Stan data
  ## ============================================================
  
  Dt2$Bd <-
    Bd
  
  Dt2$tau0 <-
    tau0
  
  Dt2$M <-
    M
  
  Dt2$K_func <-
    as.integer(
      K_func
    )
  
  Dt2$K_beta <-
    as.integer(
      K_beta
    )
  
  Dt2$K_beta_f <-
    as.integer(
      K_beta_f
    )
  
  Dt2$K_beta_r <-
    as.integer(
      K_beta_r
    )
  
  Dt2$G_beta_f <-
    G_beta_f
  
  Dt2$G_beta_r <-
    G_beta_r
  
  Dt2$id <-
    as.integer(
      id
    )
  
  Dt2$Wstar <-
    Wstar
  
  Dt2$rss_off <-
    as.numeric(
      rss_off
    )
  
  Dt2$n_off <-
    as.integer(
      n_off
    )
  
  Dt2$a_x_prior_sd <-
    as.numeric(
      a_x_prior_sd
    )
  
  Dt2$sigma_gamma_prior_mean <- as.numeric(sigma_gamma_prior_mean)
  Dt2$sigma_gamma_prior_sd <- as.numeric(sigma_gamma_prior_sd)
  Dt2$sigma_gamma_lower <- as.numeric(sigma_gamma_lower)
  Dt2$sigma_gamma_upper <- as.numeric(sigma_gamma_upper)
  
  Dt2$sigma_omega_prior_scale <-
    as.numeric(
      sigma_omega_prior_scale
    )
  
  Dt2$sigma_e_prior_scale <-
    as.numeric(
      sigma_e_prior_scale
    )
  
  Dt2$beta_fixed_prior_mean <- array(
    as.numeric(beta_fixed_prior_mean), dim = c(K_beta_f)
  )
  
  Dt2$beta_fixed_prior_sd <-
    as.numeric(
      beta_fixed_prior_sd
    )
  
  Dt2$beta_random_prior_scale <-
    as.numeric(
      beta_random_prior_scale
    )

  ## NEW: latent-PA covariate model data
  Dt2$K_xcov <- as.integer(K_xcov)
  Dt2$X_cov <- X_cov
  Dt2$x_covariate_prior_sd <- as.numeric(x_covariate_prior_sd)
  Dt2$x_covariate_in_outcome <- as.integer(x_covariate_in_outcome)
  
  ## ============================================================
  ## 11. Stan model
  ##
  ## Unit-coordinate shared-gamma hierarchy:
  ## gamma_ir ~ N(0, sigma_gamma^2), omega_ijr ~ N(0, sigma_omega_i^2).
  ## beta has its own independent mgcv basis and K_beta dimension.
  ## ============================================================
  
  Dt2$fast_likelihood <- as.integer(fast_likelihood)

  StanCode_v9_fui_bounded_gamma <- "
functions {

  real calculate_r_gamma(real gamma) {
    real log_r = log(2) + normal_lcdf(-abs(gamma) | 0, 1) + square(gamma) / 2.0;
    return exp(log_r);
  }

  real calculate_p(real gamma, real p0) {
    real r_gam = calculate_r_gamma(gamma);
    if (gamma < 0) return 1.0 + ((p0 - 1.0) / r_gam);
    else return p0 / r_gam;
  }

  real asymmetric_laplace_lpdf(real eps, real mu, real sigma, real p) {
    real ystar = (eps - mu) / sigma;
    real rho = ystar * (p - (ystar < 0 ? 1.0 : 0.0));
    return log(p) + log(1.0 - p) - log(sigma) - rho;
  }

  real gal_p0_lpdf(real eps, real mu, real sigma, real gamma, real p0) {

    if (abs(gamma) < 1e-9) {
      return asymmetric_laplace_lpdf(eps | mu, sigma, p0);
    }

    real estar = (eps - mu) / sigma;
    real p = calculate_p(gamma, p0);
    real pg_plus  = p - (gamma > 0 ? 1.0 : 0.0);
    real pg_minus = p - (gamma < 0 ? 1.0 : 0.0);
    real abs_g = abs(gamma);
    real ratio_mp = pg_minus / pg_plus;

    real arg1 = (-pg_plus * (estar / abs_g)) + (ratio_mp * abs_g);
    real arg2 = ratio_mp * abs_g;

    real log_cdf_arg1 = normal_lcdf(arg1 | 0, 1);
    real log_cdf_arg2 = normal_lcdf(arg2 | 0, 1);

    real log_abs_part1;
    real sgn_part1;

    if ((log_cdf_arg1 - log_cdf_arg2) < 0) {
      log_abs_part1 = log_cdf_arg2
                      + log1m_exp(log_cdf_arg1 - log_cdf_arg2)
                      + ((-pg_minus * estar)
                         + ((square(gamma) / 2.0) * square(ratio_mp)))
                      + log((estar / gamma) > 0 ? 1.0 : 0.0);
      sgn_part1 = -1.0;
    } else {
      log_abs_part1 = log_cdf_arg1
                      + log1m_exp(log_cdf_arg2 - log_cdf_arg1)
                      + ((-pg_minus * estar)
                         + ((square(gamma) / 2.0) * square(ratio_mp)))
                      + log((estar / gamma) > 0 ? 1.0 : 0.0);
      sgn_part1 = 1.0;
    }

    real ind_term = ((estar / gamma) > 0) ? 1.0 : 0.0;
    real arg3 = -abs_g + ((pg_plus * estar * ind_term) / abs_g);
    real log_part2 = normal_lcdf(arg3 | 0, 1)
                     + ((-pg_plus * estar) + (square(gamma) / 2.0));

    real val;

    if (sgn_part1 == 1.0) {
      val = log(2.0) + log(p) + log(1.0 - p) - log(sigma)
            + log_sum_exp(log_abs_part1, log_part2);
    } else {
      if (log_part2 > log_abs_part1) {
        val = log(2.0) + log(p) + log(1.0 - p) - log(sigma)
              + log_part2
              + log1m_exp(log_abs_part1 - log_part2);
      } else {
        val = negative_infinity();
      }
    }

    return val;
  }

  // Same normalized GAL density as the scalar function, with shared
  // expressions evaluated once. Keep scalar gal_p0_lpdf for pointwise LOO.
  real gal_p0_sum_lpdf(vector y, vector mu, real sigma, real gamma, real p0) {
    int Ny = num_elements(y);
    vector[Ny] ll;
    real abs_g = abs(gamma);

    if (abs_g < 1e-9) {
      real log_norm = log(p0) + log(1.0 - p0) - log(sigma);
      for (n in 1:Ny) {
        real ystar = (y[n] - mu[n]) / sigma;
        real rho = ystar * (p0 - (ystar < 0 ? 1.0 : 0.0));
        ll[n] = log_norm - rho;
      }
    } else {
      real p = calculate_p(gamma, p0);
      real pg_plus = p - (gamma > 0 ? 1.0 : 0.0);
      real pg_minus = p - (gamma < 0 ? 1.0 : 0.0);
      real gamma_sq_half = square(gamma) / 2.0;
      real ratio_mp = pg_minus / pg_plus;
      real arg2 = ratio_mp * abs_g;
      real log_cdf_arg2 = normal_lcdf(arg2 | 0, 1);
      real shift1 = gamma_sq_half * square(ratio_mp);
      real log_norm = log(2.0) + log(p) + log(1.0 - p) - log(sigma);

      for (n in 1:Ny) {
        real estar = (y[n] - mu[n]) / sigma;
        real arg1 = (-pg_plus * (estar / abs_g)) + arg2;
        real log_cdf_arg1 = normal_lcdf(arg1 | 0, 1);
        real log_abs_part1;
        real sgn_part1;
        real ind_term = ((estar / gamma) > 0) ? 1.0 : 0.0;
        real arg3 = -abs_g + ((pg_plus * estar * ind_term) / abs_g);
        real log_part2 = normal_lcdf(arg3 | 0, 1)
                        + ((-pg_plus * estar) + gamma_sq_half);

        if ((log_cdf_arg1 - log_cdf_arg2) < 0) {
          log_abs_part1 = log_cdf_arg2
                          + log1m_exp(log_cdf_arg1 - log_cdf_arg2)
                          + ((-pg_minus * estar) + shift1)
                          + log((estar / gamma) > 0 ? 1.0 : 0.0);
          sgn_part1 = -1.0;
        } else {
          log_abs_part1 = log_cdf_arg1
                          + log1m_exp(log_cdf_arg2 - log_cdf_arg1)
                          + ((-pg_minus * estar) + shift1)
                          + log((estar / gamma) > 0 ? 1.0 : 0.0);
          sgn_part1 = 1.0;
        }

        if (sgn_part1 == 1.0) {
          ll[n] = log_norm + log_sum_exp(log_abs_part1, log_part2);
        } else {
          if (log_part2 > log_abs_part1) {
            ll[n] = log_norm + log_part2
                    + log1m_exp(log_abs_part1 - log_part2);
          } else {
            ll[n] = negative_infinity();
          }
        }
      }
    }
    return sum(ll);
  }
}

data {
  int<lower=1> N;
  vector[N] Y;

  int<lower=1> K;
  matrix[N, K] X;
  int prior_only;
  int<lower=0, upper=1> fast_likelihood;

  vector[2] Bd;
  real<lower=0, upper=1> tau0;

  int<lower=1> M; //Total number of day-level curves
  int<lower=1> K_func; //number of measurement directions

  // NEW: subject-level covariates for the latent PA process.
  int<lower=0> K_xcov;
  matrix[N, K_xcov] X_cov;
  real<lower=0> x_covariate_prior_sd;
  int<lower=0, upper=1> x_covariate_in_outcome;

  int<lower=1> K_beta;
  int<lower=1> K_beta_f;
  int<lower=1> K_beta_r;

  array[M] int<lower=1, upper=N> id; // The subject of the m-th curve
  matrix[M, K_func] Wstar; //projected observed score

  real<lower=0> rss_off; //residual sum of squares
  int<lower=1> n_off; //residual sum of squares

  // Link the mgcv fixed and positive-penalty beta directions to Bx.
  matrix[K_func, K_beta_f] G_beta_f;
  matrix[K_func, K_beta_r] G_beta_r;

  real<lower=0> a_x_prior_sd;         
  real<lower=0> sigma_gamma_lower;
  real<lower=sigma_gamma_lower> sigma_gamma_upper;
  real<lower=sigma_gamma_lower, upper=sigma_gamma_upper> sigma_gamma_prior_mean;
  real<lower=0> sigma_gamma_prior_sd;
  real<lower=0> sigma_omega_prior_scale;
  real<lower=0> sigma_e_prior_scale;
  vector[K_beta_f] beta_fixed_prior_mean;
  real<lower=0> beta_fixed_prior_sd;
  real<lower=0> beta_random_prior_scale;
}

transformed data {
  int Kc = K - 1;
  matrix[N, Kc] Xc;
  vector[Kc] means_X;
  vector[N] measurement_count = rep_vector(0.0, N);
  matrix[N, K_func] measurement_mean = rep_matrix(0.0, N, K_func);
  vector[N] measurement_within_ss = rep_vector(0.0, N);

  for (i in 2:K) {
    means_X[i - 1] = mean(X[, i]);
    Xc[, i - 1] = X[, i] - means_X[i - 1];
  }

  // Data-only work: once per chain, not once per leapfrog gradient.
  // Use centered residuals to compute within-subject variation stably.
  if (fast_likelihood == 1) {
    for (m in 1:M) {
      measurement_count[id[m]] += 1.0;
      measurement_mean[id[m]] += Wstar[m];
    }
    for (i in 1:N) {
      if (measurement_count[i] > 0) {
        measurement_mean[i] /= measurement_count[i];
      }
    }
    for (m in 1:M) {
      measurement_within_ss[id[m]] +=
        dot_self(Wstar[m] - measurement_mean[id[m]]);
    }
  }
}

parameters {
  vector[Kc] b;
  real Intercept_centered;
  real log_sigma;
  real<lower=0, upper=1> U;

  vector[K_func] a_x;          // Population mean coefficients

  // NEW: each row is one scalar-covariate effect across Bx directions.
  matrix[K_xcov, K_func] alpha_x;

  matrix[N, K_func] z_gamma;   // Non-centered residual subject effects
  real<lower=sigma_gamma_lower, upper=sigma_gamma_upper> sigma_gamma; // ONE shared SD
  vector<lower=0>[N] sigma_omega; // One day-level SD per subject

  // Shared white-noise SD: same sigma_e in- and off-span.
  real<lower=0> sigma_e;

  // mgcv penalty-null-space coefficients.
  vector[K_beta_f] betaf;

  // Non-centered coefficients for the positive-penalty directions.
  vector[K_beta_r] z_betar;
  real<lower=0> sigma_beta;
}

transformed parameters {
  real Intercept;
  real gam;
  real<lower=0> sigma;

  vector[K_beta] b_beta;
  vector[K_beta_r] betar;

  vector[K_func] b_func;

  matrix[N, K_func] gamma;
  matrix[N, K_func] x_cov_effect;
  matrix[N, K_func] x_score;
  matrix[N, K_func] x_outcome_dev;

  // Residual subject deviations after accounting for latent-PA covariates.
  gamma = sigma_gamma * z_gamma;

  // Covariate-dependent latent-PA mean in Bx score coordinates.
  if (K_xcov > 0) {
    x_cov_effect = X_cov * alpha_x;
  } else {
    x_cov_effect = rep_matrix(0.0, N, K_func);
  }
  x_score = rep_matrix(a_x', N) + x_cov_effect + gamma;

  // Sensitivity switch:
  // 0 = beta(t) multiplies residual PA gamma only;
  // 1 = beta(t) multiplies covariate-dependent deviation + gamma.
  if (x_covariate_in_outcome == 1) {
    x_outcome_dev = x_cov_effect + gamma;
  } else {
    x_outcome_dev = gamma;
  }

  Intercept = Intercept_centered;

  gam = (Bd[2] - Bd[1]) * U + Bd[1];

  sigma = exp(log_sigma);

  // Positive-penalty beta coefficients use a non-centered hierarchy.
  betar = sigma_beta * z_betar;

  // B_beta was stored fixed first and random second in R.
  b_beta[1:K_beta_f] = betaf;
  b_beta[(K_beta_f + 1):K_beta] = betar;

  // Functional coefficient in the latent Bx-score coordinates.
  b_func = G_beta_f * betaf + G_beta_r * betar;
}

model {

  // ---------------- scalar outcome priors ----------------

  target +=
    student_t_lpdf(
      Intercept_centered |
      3,
      0,
      2.5
    );

  target +=
    normal_lpdf(
      log_sigma |
      log(0.5),
      1.0
    );

  target +=
    normal_lpdf(
      b |
      0,
      2.5
    );

  target +=
    beta_lpdf(
      U |
      4,
      4
    );

  // ---------------- functional measurement priors ---------------

   // Population mean prior
  target +=
    normal_lpdf(
      a_x |
      0,
      a_x_prior_sd
    );

  // NEW: functional covariate effects on the latent PA score vector.
  if (K_xcov > 0) {
    target += normal_lpdf(
      to_vector(alpha_x) |
      0,
      x_covariate_prior_sd
    );
  }

  // gamma_ir | sigma_gamma ~ N(0, sigma_gamma^2), non-centered, all i,r.
  target += std_normal_lpdf(to_vector(z_gamma));
  
  // Variance priors
  
  // Truncated normal: parameter bounds supply the truncation.
  // Bounds/location/SD are fixed data, so its normalization is constant.
  target += normal_lpdf(sigma_gamma | sigma_gamma_prior_mean, sigma_gamma_prior_sd);

  target +=
    normal_lpdf(
      sigma_omega |
      0,
      sigma_omega_prior_scale
    );

  target +=
    normal_lpdf(
      sigma_e |
      0,
      sigma_e_prior_scale
    );

  // ---------------- beta(t) prior ----------------
  //
  // betaf spans the mgcv penalty null space.  The positive-penalty
  // directions are betar = sigma_beta * z_betar.
  target +=
    normal_lpdf(
      betaf |
      beta_fixed_prior_mean,
      beta_fixed_prior_sd
    );

  target +=
    std_normal_lpdf(
      z_betar
    );

  target +=
    normal_lpdf(
      sigma_beta |
      0,
      beta_random_prior_scale
    );

  // ---------------- functional measurement likelihood ----------------
  //
  // Wstar_ij,r | x_score_i,r
  //   ~ N(
  //       x_score_i,r,
  //       sigma_omega_i^2 + sigma_e^2
  //     )

  {
    if (fast_likelihood == 1) {
      vector[N] measurement_var = square(sigma_omega) + square(sigma_e);
      vector[N] measurement_ss = measurement_within_ss
        + measurement_count .* rows_dot_self(measurement_mean - x_score);
      // Exact Gaussian log likelihood up to its fixed log(2*pi) constant.
      // Includes day-to-day variation: this is not a mean-only likelihood.
      target += -0.5 * (
        K_func * dot_product(measurement_count, log(measurement_var))
        + sum(measurement_ss ./ measurement_var)
      );
    } else {
      vector[M] measurement_sd = sqrt(
        square(sigma_omega[id]) + square(sigma_e)
      );
      for (r in 1:K_func) {
        Wstar[, r] ~ normal(x_score[id, r], measurement_sd);
      }
    }

    // Shared sigma_e is also identified by the off-span residual.
    
    //Off-span likelihood
    target +=
      -n_off *
        log(
          sigma_e
        )
      -
      0.5 *
        rss_off /
        square(
          sigma_e
        );
  }

  // ---------------- scalar-on-function quantile outcome ----------------
  //
  // gamma is the DEVIATION score matrix in the unit X basis.
  // The mgcv beta coefficients are mapped into the X-score coordinates by
  // G_beta_f * betaf + G_beta_r * betar.

  if (prior_only == 0) {

    vector[N] mu =
      Intercept_centered
      +
      Xc *
        b;

    mu += x_outcome_dev * b_func;

    if (fast_likelihood == 1) {
      target += gal_p0_sum_lpdf(Y | mu, sigma, gam, tau0);
    } else {
      for (n in 1:N) {
        target += gal_p0_lpdf(Y[n] | mu[n], sigma, gam, tau0);
      }
    }
  }
}

generated quantities {

  vector[N] log_lik_y;
  real b_Intercept;

  vector[K_func] b_func_gq =
    b_func;

  vector[K_beta] b_beta_gq =
    b_beta;

  b_Intercept =
    Intercept_centered
    -
    dot_product(
      means_X,
      b
    );

  {
    // Same predictor as the outcome likelihood.
    // Do not add a_x, measurement densities, priors, or lp__ here.
    vector[N] mu_y = Intercept_centered + Xc * b + x_outcome_dev * b_func;
    for (n in 1:N) {
      log_lik_y[n] = gal_p0_lpdf(Y[n] | mu_y[n], sigma, gam, tau0);
    }
  }
}
"

## ============================================================
## 12. Initial values from variance spectrum
##
## Moment initialization is computed directly in the unit-score coordinates.
## ============================================================

subj_mean_Wstar <-
  rowsum(
    Wstar,
    group = id
  ) /
  J_i

a_x_init <-
  colMeans(
    subj_mean_Wstar
  )

## NEW: initialize alpha_x by ridge-stabilized multivariate least squares
## on subject-level mean projected PA scores.  X_cov is centered, so a_x_init
## remains the empirical mean score vector.
if (K_xcov > 0L) {
  Xcov_cross <- crossprod(X_cov)
  Xcov_ridge <- 1e-6 * max(1, mean(diag(Xcov_cross)))
  alpha_x_init <- solve(
    Xcov_cross + Xcov_ridge * diag(K_xcov),
    crossprod(
      X_cov,
      sweep(subj_mean_Wstar, 2, a_x_init, "-")
    )
  )
} else {
  alpha_x_init <- matrix(0, nrow = 0L, ncol = K_func)
}

## Residual subject-specific latent-score initialization.
xi_init <- sweep(subj_mean_Wstar, 2, a_x_init, "-") - X_cov %*% alpha_x_init
xi_init <- sweep(xi_init, 2, colMeans(xi_init), "-")

sigma_e_init <-
  sigma_e_offspan_moment

if (
  !is.finite(
    sigma_e_init
  ) ||
  sigma_e_init <= 1e-6
) {
  
  sigma_e_init <-
    max(
      0.10 *
        sd(
          as.vector(
            W_obs
          )
        ),
      0.1
    )
}

## ------------------------------------------------------------
## Direction-specific empirical variance decomposition
## ------------------------------------------------------------

ss_within <-
  numeric(
    K_func
  )

for (
  i in seq_len(N)
) {
  
  rows_i <-
    which(
      id ==
        i
    )
  
  resid_i <-
    sweep(
      Wstar[
        rows_i,
        ,
        drop = FALSE
      ],
      2,
      subj_mean_Wstar[
        i,
      ],
      "-"
    )
  
  ss_within <-
    ss_within +
    colSums(
      resid_i^2
    )
}

within_r <-
  ss_within /
  (
    M -
      N
  )

## Raw between-subject variance before explaining PA by scalar covariates.
between_r_raw <-
  apply(
    subj_mean_Wstar,
    2,
    var
  )

## NEW: gamma is the residual subject effect after X_cov * alpha_x, so its
## empirical between-subject diagnostic must use the covariate-adjusted scores.
subj_mean_Wstar_resid <-
  sweep(subj_mean_Wstar, 2, a_x_init, "-") - X_cov %*% alpha_x_init

between_r <-
  apply(
    subj_mean_Wstar_resid,
    2,
    var
  )

# These empirical direction-specific moments remain diagnostics only.
# They do NOT become separate sigma_gamma parameters in this shared model.
empirical_signal_variance_unclipped <-
  between_r - within_r * mean(1 / J_i)
xi_var_r <- pmax(empirical_signal_variance_unclipped, 1e-8)

# Retain the pooled empirical moment as a diagnostic; it may exceed the bound.
sigma_gamma_moment_init <- sqrt(max(mean(empirical_signal_variance_unclipped), 1e-8))
# Start at the prior location, guaranteed to be inside the permitted interval.
# sigma_gamma is still sampled, rather than fixed at this starting value.
sigma_gamma_init <- sigma_gamma_prior_mean

## In the unit basis: within_r = E(sigma_omega_i^2) + sigma_e^2.
omega_var_r <-
  pmax(
    within_r -
      sigma_e_init^2,
    1e-8
  )

sigma_omega_base <-
  sqrt(
    mean(omega_var_r)
  )

sigma_xi_r_init <-
  sqrt(
    xi_var_r
  )

## ------------------------------------------------------------
## Subject-specific omega scale initialization
## ------------------------------------------------------------

subj_ratio <-
  numeric(
    N
  )

for (
  i in seq_len(N)
) {
  
  rows_i <-
    which(
      id ==
        i
    )
  
  resid_i <-
    sweep(
      Wstar[
        rows_i,
        ,
        drop = FALSE
      ],
      2,
      subj_mean_Wstar[
        i,
      ],
      "-"
    )
  
  v_i <-
    colMeans(
      resid_i^2
    ) *
    J_i[i] /
    (
      J_i[i] -
        1
    )
  
  subj_ratio[i] <-
    sqrt(
      max(
        mean(
          pmax(
            v_i -
              sigma_e_init^2,
            1e-8
          )
        ),
        1e-4
      )
    )
}

sigma_omega_i_init <-
  subj_ratio

## ------------------------------------------------------------
## Reliability initialization
## ------------------------------------------------------------

# In the unit basis, measurement variance is the same over r within subject.
# A shared gamma SD therefore gives the same measurement reliability over r.
noise_in_mean_i <- (sigma_omega_i_init^2 + sigma_e_init^2) / J_i
reliability_i_init <- sigma_gamma_init^2 / (sigma_gamma_init^2 + noise_in_mean_i)
reliability_init <- matrix(reliability_i_init, nrow = N, ncol = K_func)

reliability_init <-
  pmin(
    pmax(
      reliability_init,
      0
    ),
    1
  )

if (
  any(
    !is.finite(
      reliability_init
    )
  )
) {
  stop(
    "Non-finite reliability_init."
  )
}

## ------------------------------------------------------------
## mgcv beta-basis overlap / leverage diagnostics
##
## G_beta is rectangular in general, so diagonal-based v4 diagnostics
## are no longer appropriate.
##
## For X direction r:
##   sigma_gamma_init * ||G_beta[r, ]||
## summarizes how strongly that latent direction can communicate with
## the beta space.
##
## For beta direction s:
##   sigma_gamma_init * sqrt(sum_r G_beta[r,s]^2)
## summarizes its outcome leverage through the latent-X distribution.
## ------------------------------------------------------------

x_to_beta_leverage <- sigma_gamma_init * G_row_norm
G_weighted <- sigma_gamma_init * G_beta

beta_coef_eta_scale <-
  sqrt(
    colSums(
      G_weighted^2
    )
  )

if (
  any(!is.finite(beta_coef_eta_scale)) ||
  max(beta_coef_eta_scale) <= 0
) {
  stop(
    "All effective beta-direction scales are zero/non-finite; beta cannot link to the outcome."
  )
}

beta_leverage_table <-
  data.frame(
    beta_direction = seq_len(K_beta),
    beta_block = c(
      rep("fixed", K_beta_f),
      rep("random", K_beta_r)
    ),
    eta_scale = beta_coef_eta_scale,
    G_column_norm = G_col_norm,
    d_beta = d_beta
  )

prior_sd_b <-
  c(
    rep(
      beta_fixed_prior_sd,
      K_beta_f
    ),
    rep(
      beta_random_prior_scale,
      K_beta_r
    )
  )

eta_prior_sd <-
  sqrt(
    sum(
      (
        beta_coef_eta_scale *
          prior_sd_b
      )^2
    ) + sigma_gamma_init^2 *
      sum(as.vector(G_beta_f %*% beta_fixed_prior_mean)^2)
  )

cat(
  "Initialization from the variance spectrum:\n"
)

cat(
  "sigma_e_init       =",
  round(
    sigma_e_init,
    6
  ),
  "\n"
)

cat(
  "sigma_omega_base   =",
  round(
    sigma_omega_base,
    6
  ),
  "\n"
)

cat("shared sigma_gamma_init =", signif(sigma_gamma_init, 4), "\n")
cat("median subject reliability at initialization =",
    signif(median(reliability_i_init), 3), "\n\n")

cat(
  "Shared-gamma X-direction overlap leverage (sigma_gamma * ||G_row_r||):\n"
)

cat(
  paste(
    signif(
      x_to_beta_leverage,
      3
    ),
    collapse = ", "
  ),
  "\n"
)

cat(
  "Shared-gamma mgcv beta coefficient eta scales (fixed first, then random):\n"
)

cat(
  paste(
    signif(
      beta_coef_eta_scale,
      3
    ),
    collapse = ", "
  ),
  "\n"
)

cat(
  "implied approximate prior SD of the functional term in eta =",
  round(
    eta_prior_sd,
    3
  ),
  " (Y is standardized)\n\n"
)

if (
  min(G_col_rel) < 1e-4
) {
  warning(
    "At least one beta direction has extremely weak overlap with the X space. ",
    "Inspect G_col_norm / singular values before interpreting beta."
  )
}

if (
  eta_prior_sd > 20
) {
  warning(
    "The beta prior implies a very large functional effect. ",
    "Consider reducing beta_fixed_prior_sd / beta_random_prior_scale."
  )
}

xi_dev_init <-
  reliability_init *
  xi_init

outcome_intercept_init <-
  as.numeric(
    quantile(
      Y_std,
      probs = tau0,
      na.rm = TRUE
    )
  )

if (
  abs(
    tau0 -
    0.1
  ) < 1e-8
) {
  
  U_init_center <-
    0.27
  
  outcome_sigma_init <-
    0.23
  
} else {
  
  U_init_center <-
    0.50
  
  outcome_sigma_init <-
    0.20
}

stan_vector_init <- function(
    x
) {
  
  x <-
    as.numeric(
      x
    )
  
  array(
    x,
    dim = length(x)
  )
}

init_fun <- function(
    chain_id = 1
) {
  
  set.seed(
    seed +
      1000L *
      chain_id
  )
  
  Kc_now <-
    as.integer(
      Dt2$K -
        1L
    )
  
  a_x_chain <-
    a_x_init +
    rnorm(
      K_func,
      0,
      0.01 *
        max(
          sd(
            a_x_init
          ),
          1e-3
        )
    )
  
  # Scalar initialization is required by the scalar Stan parameter.
  gamma_interval_width <- sigma_gamma_upper - sigma_gamma_lower
  gamma_init_fraction <- (sigma_gamma_init - sigma_gamma_lower) / gamma_interval_width
  gamma_init_fraction <- pmin(pmax(gamma_init_fraction, 0.05), 0.95)
  # Jitter on an unconstrained logit scale, then return strictly inside the bounds.
  sigma_gamma_chain <- as.numeric(sigma_gamma_lower + gamma_interval_width *
                                    plogis(qlogis(gamma_init_fraction) + rnorm(1, 0, 0.15)))
  z_gamma_chain <- xi_dev_init / sigma_gamma_chain +
    matrix(rnorm(N * K_func, 0, 0.02), N, K_func)
  
  sigma_beta_chain <-
    0.05 *
    exp(
      rnorm(
        1,
        0,
        0.05
      )
    )
  
  out <-
    list(
      U =
        plogis(
          qlogis(
            U_init_center
          ) +
            rnorm(
              1,
              0,
              0.03
            )
        ),
      
      log_sigma =
        log(
          outcome_sigma_init
        ) +
        rnorm(
          1,
          0,
          0.03
        ),
      
      b =
        stan_vector_init(
          rnorm(
            Kc_now,
            0,
            0.02
          )
        ),
      
      Intercept_centered =
        outcome_intercept_init +
        rnorm(
          1,
          0,
          0.02
        ),
      
      a_x =
        stan_vector_init(
          a_x_chain
        ),

      alpha_x =
        alpha_x_init +
        matrix(
          rnorm(K_xcov * K_func, 0, 0.01),
          nrow = K_xcov,
          ncol = K_func
        ),
      
      z_gamma =
        z_gamma_chain,
      
      sigma_gamma = as.numeric(sigma_gamma_chain),
      
      sigma_omega =
        stan_vector_init(
          sigma_omega_i_init *
            exp(
              rnorm(
                N,
                0,
                0.03
              )
            )
        ),
      
      sigma_e =
        sigma_e_init *
        exp(
          rnorm(
            1,
            0,
            0.03
          )
        ),
      
      ## mgcv penalty-null-space coefficients.
      betaf =
        stan_vector_init(
          rnorm(
            K_beta_f,
            beta_fixed_prior_mean,
            min(0.02, beta_fixed_prior_sd)
          )
        ),
      
      ## Standard-normal coordinates for positive-penalty directions.
      z_betar =
        stan_vector_init(
          rnorm(
            K_beta_r,
            0,
            0.05
          )
        ),
      
      sigma_beta =
        sigma_beta_chain
    )
  
  out
}

## ============================================================
## 13. Fit Stan
## ============================================================

pars_keep <-
  c(
    "b",
    "Intercept_centered",
    "Intercept",
    "b_Intercept",
    "log_sigma",
    "sigma",
    "U",
    "gam",
    "a_x",
    if (K_xcov > 0L) "alpha_x",
    "gamma",
    "sigma_gamma",
    "sigma_omega",
    "sigma_e",
    "betaf",
    "z_betar",
    "betar",
    "sigma_beta",
    "b_beta",
    "b_beta_gq",
    "b_func_gq",
    "log_lik_y"
  )

fit_tmp <-
  rstan::stan(
    model_code = StanCode_v9_fui_bounded_gamma,
    data = Dt2,
    iter = Nsim,
    warmup = warmup,
    chains = chains,
    refresh = 50,
    init = init_fun,
    pars = pars_keep,
    include = TRUE,
    control = list(
      adapt_delta = adapt_delta,
      max_treedepth =
        as.integer(
          max_treedepth
        )
    ),
    cores = cores,
    seed = seed
  )

## ============================================================
## 14. Sampler diagnostics
## ============================================================

sampler_params <-
  rstan::get_sampler_params(
    fit_tmp,
    inc_warmup = FALSE
  )

sampler_check <-
  do.call(
    rbind,
    lapply(
      seq_along(
        sampler_params
      ),
      function(
    ch
      ) {
        
        x <-
          sampler_params[[ch]]
        
        data.frame(
          chain = ch,
          mean_leapfrog = mean(x[, "n_leapfrog__"]),
          median_leapfrog = median(x[, "n_leapfrog__"]),
          
          divergent =
            sum(
              x[
                ,
                "divergent__"
              ]
            ),
          
          hit_max_treedepth =
            sum(
              x[
                ,
                "treedepth__"
              ] >=
                max_treedepth
            ),
          
          mean_treedepth =
            mean(
              x[
                ,
                "treedepth__"
              ]
            ),
          
          mean_accept =
            mean(
              x[
                ,
                "accept_stat__"
              ]
            ),
          
          mean_stepsize =
            mean(
              x[
                ,
                "stepsize__"
              ]
            )
        )
      }
    )
  )

cat(
  "\nPer-chain sampler diagnostics:\n"
)

print(
  sampler_check,
  row.names = FALSE
)

## ============================================================
## 15. Shared gamma variance and subject-specific omega diagnostics
## ============================================================

# Exactly ONE sigma_gamma row. Mean/sd below describe the posterior of an SD.
sigma_gamma_summary_raw <- rstan::summary(
  fit_tmp, pars = "sigma_gamma", probs = c(0.025, 0.5, 0.975)
)$summary
sigma_gamma_summary <- data.frame(
  parameter = "sigma_gamma",
  sigma_gamma_summary_raw["sigma_gamma", , drop = FALSE],
  row.names = NULL, check.names = FALSE
)

variance_draws_shared <- rstan::extract(
  fit_tmp, pars = c("sigma_gamma", "sigma_omega", "sigma_e"), permuted = TRUE
)
sigma_gamma_draw <- as.numeric(variance_draws_shared$sigma_gamma)
n_post_draws <- length(sigma_gamma_draw)
stopifnot(nrow(sigma_gamma_summary_raw) == 1L, n_post_draws >= 2L)
gamma_var_draw <- sigma_gamma_draw^2
gamma_bound_width <- sigma_gamma_upper - sigma_gamma_lower
sigma_gamma_bound_check <- data.frame(
  lower = sigma_gamma_lower,
  upper = sigma_gamma_upper,
  posterior_min = min(sigma_gamma_draw),
  posterior_max = max(sigma_gamma_draw),
  fraction_in_bottom_5pct_of_interval = mean(
    sigma_gamma_draw <= sigma_gamma_lower + 0.05 * gamma_bound_width),
  fraction_in_top_5pct_of_interval = mean(
    sigma_gamma_draw >= sigma_gamma_upper - 0.05 * gamma_bound_width)
)
scalar_draw_summary <- function(x) {
  data.frame(mean = mean(x), sd = sd(x),
             lower = unname(quantile(x, 0.025)),
             median = unname(quantile(x, 0.5)),
             upper = unname(quantile(x, 0.975)))
}
sigma_gamma_variance_summary <- data.frame(
  parameter = "sigma_gamma^2", scalar_draw_summary(gamma_var_draw)
)

# D x N draws: each subject's omega SD is shared over directions and days.
sigma_omega_draw <- matrix(variance_draws_shared$sigma_omega,
                           nrow = n_post_draws, ncol = N)
omega_summary_raw <- rstan::summary(
  fit_tmp, pars = "sigma_omega", probs = c(0.025, 0.5, 0.975)
)$summary
omega_names <- paste0("sigma_omega[", seq_len(N), "]")
sigma_omega_subject_summary <- data.frame(
  subject_index = seq_len(N),
  omega_summary_raw[omega_names, , drop = FALSE],
  variance_mean = colMeans(sigma_omega_draw^2),
  row.names = NULL, check.names = FALSE
)

noise_mean_variance_post <- sweep(
  sweep(sigma_omega_draw^2, 1, variance_draws_shared$sigma_e^2, "+"),
  2, J_i, "/"
)
signal_variance_post <- matrix(gamma_var_draw, nrow = n_post_draws, ncol = N)
reliability_subject_draw <- signal_variance_post /
  (signal_variance_post + noise_mean_variance_post)
reliability_subject_summary <- data.frame(
  subject_index = seq_len(N),
  posterior_mean = colMeans(reliability_subject_draw),
  lower = apply(reliability_subject_draw, 2, quantile, probs = 0.025),
  median = apply(reliability_subject_draw, 2, quantile, probs = 0.5),
  upper = apply(reliability_subject_draw, 2, quantile, probs = 0.975)
)
average_reliability_draw <- rowMeans(reliability_subject_draw)
reliability_summary <- data.frame(
  quantity = "Subject-averaged residual-gamma reliability (shared over directions)",
  scalar_draw_summary(average_reliability_draw)
)

# The empirical MOM spectrum may differ over r. The fitted shared model has
# ONE reliability per subject; its subject average is identical over r.
# Retain this compatibility table, explicitly labeling the repeated posterior.
empirical_mom_reliability <- pmin(pmax(xi_var_r / pmax(between_r, 1e-8), 0), 1)
empirical_variance_direction_summary <- data.frame(
  direction = seq_len(K_func), within_variance = within_r,
  between_mean_variance_raw = between_r_raw,
  between_mean_variance_adjusted = between_r,
  empirical_residual_signal_variance = xi_var_r,
  empirical_signal_SD = sigma_xi_r_init,
  empirical_MOM = empirical_mom_reliability
)
reliability_direction_summary <- data.frame(
  direction = seq_len(K_func), empirical_MOM = empirical_mom_reliability,
  posterior_mean = rep(reliability_summary$mean, K_func),
  lower = rep(reliability_summary$lower, K_func),
  median = rep(reliability_summary$median, K_func),
  upper = rep(reliability_summary$upper, K_func),
  posterior_shared_over_directions = TRUE
)

# Population between-subject variance of gamma_i(t), ORIGINAL PA units squared.
# This is different from the posterior uncertainty of a particular gamma_i(t).
basis_variance_weight <- rowSums(Bx^2)
gamma_variance_curve_summary <- data.frame(
  t = tau_x,
  mean = basis_variance_weight * mean(gamma_var_draw),
  lower = basis_variance_weight * unname(quantile(gamma_var_draw, 0.025)),
  median = basis_variance_weight * unname(quantile(gamma_var_draw, 0.5)),
  upper = basis_variance_weight * unname(quantile(gamma_var_draw, 0.975))
)
variance_parameter_counts <- data.frame(
  parameter = c("sigma_gamma", "sigma_omega", "sigma_e"),
  count = c(1L, N, 1L),
  sharing = c("all subjects and all directions", "one per subject; shared over directions/days",
              "all observations")
)
cat("\nVariance parameter structure:\n")
print(variance_parameter_counts, row.names = FALSE)
cat("\nONE shared gamma SD posterior:\n")
print(sigma_gamma_summary, row.names = FALSE, digits = 4)
cat("\nONE shared gamma variance posterior:\n")
print(sigma_gamma_variance_summary, row.names = FALSE, digits = 4)
cat("\nPosterior position within the hard SD bounds:\n")
print(sigma_gamma_bound_check, row.names = FALSE, digits = 4)
cat("\nSubject-averaged residual-gamma reliability (same over directions):\n")
print(reliability_summary, row.names = FALSE, digits = 4)

## ============================================================
## 16. Attach useful objects
## ============================================================

attr(
  fit_tmp,
  "model_version"
) <-
  paste0(
    "M2_PENALTY_ORTHO_v9_shared_gamma_joint_fui_bounded_xcov"
  )

attr(
  fit_tmp,
  "functional_math"
) <-
  paste0(
    "Bx' W Bx = I; Wstar_ij,r ~ N(a_r + Xcov_i alpha_r + gamma_ir, ",
    "sigma_omega_i^2 + sigma_e^2); ",
    "gamma_ir ~ N(0, sigma_gamma^2) after latent-PA covariate adjustment, ",
    "ONE sigma_gamma shared over all i,r; ",
    "X(t) uses the unit-normalized penalty-ordered basis; ",
    "beta(t) uses the explicit mgcv penalty basis Psi U E^{-1}; ",
    "betar = sigma_beta z_betar; ",
    "outcome functional term uses gamma only when x_covariate_in_outcome=FALSE, ",
    "and X_cov alpha_x + gamma when TRUE"
  )

attr(
  fit_tmp,
  "x_basis"
) <-
  "bspline_second_diff_penalty_ordered_unit_normalized"

attr(
  fit_tmp,
  "beta_basis"
) <-
  "mgcv_penalty_eigen_reparameterized"

attr(
  fit_tmp,
  "K_func"
) <-
  K_func

attr(
  fit_tmp,
  "basis_df"
) <-
  basis_df

attr(
  fit_tmp,
  "basis_degree"
) <-
  basis_degree

attr(
  fit_tmp,
  "basis_ridge"
) <-
  basis_ridge

attr(
  fit_tmp,
  "basis_tol"
) <-
  basis_tol

attr(
  fit_tmp,
  "basis_weights"
) <-
  basis_weights

attr(
  fit_tmp,
  "B0_x"
) <-
  B0_x

attr(
  fit_tmp,
  "P_x"
) <-
  P_x

attr(
  fit_tmp,
  "ob_x"
) <-
  ob_x

attr(
  fit_tmp,
  "Bx"
) <-
  Bx

attr(
  fit_tmp,
  "Bx_orthogonal"
) <-
  Bx_orthogonal

attr(
  fit_tmp,
  "d_x_orthogonal"
) <-
  d_x_orthogonal

attr(
  fit_tmp,
  "B"
) <-
  Bx

attr(
  fit_tmp,
  "d_x"
) <-
  d_x

attr(
  fit_tmp,
  "d"
) <-
  d_x

attr(
  fit_tmp,
  "Proj_x"
) <-
  Proj_x

attr(fit_tmp, "unit_score_transform_error") <-
  unit_score_transform_error

attr(fit_tmp, "unit_linkage_transform_error") <-
  unit_linkage_transform_error

attr(
  fit_tmp,
  "w_basis"
) <-
  w_basis

## Explicit mgcv beta-basis objects
attr(
  fit_tmp,
  "K_beta"
) <-
  K_beta

attr(
  fit_tmp,
  "K_beta_requested"
) <-
  beta_basis$k_requested

attr(
  fit_tmp,
  "K_beta_f"
) <-
  K_beta_f

attr(
  fit_tmp,
  "K_beta_r"
) <-
  K_beta_r

attr(
  fit_tmp,
  "beta_bs"
) <-
  beta_bs

attr(
  fit_tmp,
  "beta_basis_object"
) <-
  beta_basis

attr(
  fit_tmp,
  "Psi_mat"
) <-
  Psi_mat

attr(
  fit_tmp,
  "S_mat_beta_raw"
) <-
  S_mat_beta_raw

attr(
  fit_tmp,
  "S_mat_beta"
) <-
  S_mat_beta

attr(
  fit_tmp,
  "U_beta"
) <-
  U_beta

attr(
  fit_tmp,
  "E_beta"
) <-
  E_beta

attr(
  fit_tmp,
  "rank_beta"
) <-
  rank_beta

attr(
  fit_tmp,
  "B_beta_f"
) <-
  B_beta_f

attr(
  fit_tmp,
  "B_beta_r"
) <-
  B_beta_r

attr(
  fit_tmp,
  "B_beta"
) <-
  B_beta

attr(
  fit_tmp,
  "d_beta"
) <-
  d_beta

attr(
  fit_tmp,
  "Gram_beta"
) <-
  Gram_beta

attr(
  fit_tmp,
  "G_beta_f"
) <-
  G_beta_f

attr(
  fit_tmp,
  "G_beta_r"
) <-
  G_beta_r

attr(
  fit_tmp,
  "G_beta"
) <-
  G_beta

attr(
  fit_tmp,
  "beta_linkage_error"
) <-
  beta_linkage_error

attr(
  fit_tmp,
  "beta_leverage_table"
) <-
  beta_leverage_table

attr(
  fit_tmp,
  "G_row_norm"
) <-
  G_row_norm

attr(
  fit_tmp,
  "G_row_rel"
) <-
  G_row_rel

attr(
  fit_tmp,
  "G_col_norm"
) <-
  G_col_norm

attr(
  fit_tmp,
  "G_col_rel"
) <-
  G_col_rel

attr(
  fit_tmp,
  "G_singular"
) <-
  G_singular

attr(
  fit_tmp,
  "x_to_beta_leverage"
) <-
  x_to_beta_leverage

attr(
  fit_tmp,
  "beta_coef_eta_scale"
) <-
  beta_coef_eta_scale

attr(
  fit_tmp,
  "eta_prior_sd"
) <-
  eta_prior_sd

## measurement / variance objects
attr(
  fit_tmp,
  "rss_off"
) <-
  rss_off

attr(
  fit_tmp,
  "n_off"
) <-
  n_off

attr(
  fit_tmp,
  "offspan_fraction"
) <-
  offspan_fraction

attr(
  fit_tmp,
  "sigma_e_offspan_moment"
) <-
  sigma_e_offspan_moment

attr(fit_tmp, "sigma_gamma_init") <- sigma_gamma_init
attr(fit_tmp, "sigma_gamma_moment_init") <- sigma_gamma_moment_init
attr(fit_tmp, "sigma_gamma_bound_check") <- sigma_gamma_bound_check
attr(fit_tmp, "gamma_prior") <- list(
  family = "truncated_normal",
  normal_location = sigma_gamma_prior_mean,
  normal_sd = sigma_gamma_prior_sd,
  lower = sigma_gamma_lower,
  upper = sigma_gamma_upper,
  gamma_variance_bounds = c(lower = sigma_gamma_lower^2, upper = sigma_gamma_upper^2),
  interpretation = "FUI-calibrated sensitivity analysis; bounds are chosen, not a FUI CI"
)
attr(fit_tmp, "sigma_gamma_shared") <- TRUE
attr(fit_tmp, "sigma_gamma_summary") <- sigma_gamma_summary
attr(fit_tmp, "sigma_gamma_variance_summary") <- sigma_gamma_variance_summary
attr(fit_tmp, "sigma_omega_subject_summary") <- sigma_omega_subject_summary
attr(fit_tmp, "reliability_subject_summary") <- reliability_subject_summary
attr(fit_tmp, "reliability_summary") <- reliability_summary
attr(fit_tmp, "empirical_variance_direction_summary") <- empirical_variance_direction_summary
attr(fit_tmp, "gamma_variance_curve_summary") <- gamma_variance_curve_summary
attr(fit_tmp, "variance_parameter_counts") <- variance_parameter_counts
attr(fit_tmp, "variance_structure") <- list(
  gamma = "shared_across_subjects_and_directions",
  omega = "subject_specific_shared_across_directions_and_days",
  epsilon = "shared_global"
)
attr(fit_tmp, "tau0") <- tau0

attr(
  fit_tmp,
  "reliability_direction_summary"
) <-
  reliability_direction_summary

attr(
  fit_tmp,
  "empirical_gamma_sd_by_direction"
) <-
  sigma_xi_r_init

attr(
  fit_tmp,
  "within_r"
) <-
  within_r

attr(
  fit_tmp,
  "between_r"
) <-
  between_r

attr(fit_tmp, "between_r_raw") <- between_r_raw

attr(
  fit_tmp,
  "reliability_init"
) <-
  reliability_init

attr(
  fit_tmp,
  "Wstar"
) <-
  Wstar

attr(
  fit_tmp,
  "id"
) <-
  id

attr(
  fit_tmp,
  "J_i"
) <-
  J_i

attr(
  fit_tmp,
  "W_mean_curve_emp"
) <-
  W_mean_curve_emp

attr(
  fit_tmp,
  "quad_w"
) <-
  quad_w

attr(
  fit_tmp,
  "time_grid"
) <-
  time_grid

## NEW: latent-PA covariate bookkeeping.
attr(fit_tmp, "X_cov") <- X_cov
attr(fit_tmp, "has_x_covariates") <- has_x_covariates
attr(fit_tmp, "x_covariate_names") <- x_covariate_names
attr(fit_tmp, "x_numeric_Z") <- x_numeric_Z
attr(fit_tmp, "x_factor_Z") <- x_factor_Z
attr(fit_tmp, "x_numeric_center") <- x_numeric_center
attr(fit_tmp, "x_numeric_scale") <- x_numeric_scale
attr(fit_tmp, "x_covariate_design_center") <- x_covariate_design_center
attr(fit_tmp, "x_covariate_prior_sd") <- x_covariate_prior_sd
attr(fit_tmp, "x_covariate_in_outcome") <- x_covariate_in_outcome

attr(
  fit_tmp,
  "periodic_time"
) <-
  periodic_time

attr(fit_tmp, "outcome_xi_mode") <- "joint"

attr(
  fit_tmp,
  "Y_mean"
) <-
  Y_mean

attr(
  fit_tmp,
  "Y_sd"
) <-
  Y_sd

## Record the exact beta priors used by this fit.
attr(fit_tmp, "beta_prior") <- list(
  family = "Half-normal on sigma_beta",
  betaf_mean = beta_fixed_prior_mean,
  betaf_sd = beta_fixed_prior_sd,
  beta_random_prior_scale = beta_random_prior_scale,
  betar_conditional = "Independent Normal(0, sigma_beta^2)",
  direction_multipliers = "None"
)
attr(fit_tmp, "beta_fixed_prior_summary") <- data.frame(
  t = time_grid,
  mean = as.vector(Y_sd * (B_beta_f %*% beta_fixed_prior_mean)),
  sd = as.vector(Y_sd * beta_fixed_prior_sd * sqrt(rowSums(B_beta_f^2)))
)

attr(
  fit_tmp,
  "prior_scales"
) <-
  list(
    a_x_prior_sd =
      a_x_prior_sd,
    
    sigma_gamma_prior_mean = sigma_gamma_prior_mean,
    sigma_gamma_prior_sd = sigma_gamma_prior_sd,
    sigma_gamma_lower = sigma_gamma_lower,
    sigma_gamma_upper = sigma_gamma_upper,
    
    sigma_omega_prior_scale =
      sigma_omega_prior_scale,
    
    sigma_e_prior_scale =
      sigma_e_prior_scale,
    
    beta_fixed_prior_mean = beta_fixed_prior_mean,
    beta_fixed_prior_sd =
      beta_fixed_prior_sd,
    
    beta_random_prior_scale =
      beta_random_prior_scale
  )

attr(
  fit_tmp,
  "warmup"
) <-
  warmup

attr(
  fit_tmp,
  "adapt_delta"
) <-
  adapt_delta

attr(
  fit_tmp,
  "max_treedepth"
) <-
  max_treedepth

attr(
  fit_tmp,
  "sampler_check"
) <-
  sampler_check

# LOO bookkeeping; no change to the fitted posterior.
attr(fit_tmp, "loo_y") <- as.numeric(Y_use)
attr(fit_tmp, "loo_input_row") <- which(idx)
attr(fit_tmp, "loo_y_mean") <- Y_mean
attr(fit_tmp, "loo_y_sd") <- Y_sd
attr(fit_tmp, "loo_tau") <- tau0
attr(fit_tmp, "loo_X") <- Dt2$X
attr(fit_tmp, "loo_formula") <-
  if (inherits(form_Z, "formula")) paste(deparse(form_Z), collapse = " ") else form_Z
attr(fit_tmp, "loo_log_lik_scale") <- "standardized_Y"
attr(fit_tmp, "loo_method") <- "Joint"
attr(fit_tmp, "fast_likelihood") <- fast_likelihood
attr(fit_tmp, "elapsed_time") <- rstan::get_elapsed_time(fit_tmp)

return(
  fit_tmp
)
}

## ============================================================
## Example: use the projected tau=0.1/0.5 betaf priors.
## This block is disabled: source() never reruns MCMC.
## fast_likelihood = TRUE is the default; FALSE evaluates the original form.
if (FALSE) {
  # source() only defines functions. Copy these calls to your R session.
  # Use the same formula01 that was supplied to the corresponding FUI fits.
  # Both calls use the gamma settings in the user's current examples.
  fit_joint01_betaf_loo <- fit_flfoqr_fui_bounded_gamma_betaf_xcov(
    Y = Y_use, W = W_use, Z = Z_use,
    Nsim = 1000, warmup = 500, tau0 = 0.1,
    numeric_Z = c("AgeYR", "HEI"),
    factor_Z = c("Gender", "Race", "HealthCondt2"),
    scalar_formula = formula01,
    x_numeric_Z = "AgeYR", x_factor_Z = "Gender",
    x_covariate_prior_sd = 10, x_covariate_in_outcome = FALSE,
    K_func = 12, basis_weights = "unit",
    K_beta = 10, beta_bs = "cc", periodic_time = FALSE,
    sigma_gamma_prior_mean = 20,
    sigma_gamma_prior_sd = 0.5,
    sigma_gamma_lower = 18,
    sigma_gamma_upper = 21,
    a_x_prior_sd = 100,
    sigma_omega_prior_scale = 15,
    sigma_e_prior_scale = 10,
    beta_fixed_prior_mean = 0.08067051,
    beta_fixed_prior_sd = 0.01349938,
    beta_random_prior_scale = 1,
    chains = 2, cores = 2, seed = 1123
  )

  fit_joint05_betaf_loo <- fit_flfoqr_fui_bounded_gamma_betaf_xcov(
    Y = Y_use, W = W_use, Z = Z_use,
    Nsim = 1000, warmup = 500, tau0 = 0.5,
    numeric_Z = c("AgeYR", "HEI"),
    factor_Z = c("Gender", "Race", "HealthCondt2"),
    scalar_formula = formula01,
    x_numeric_Z = "AgeYR", x_factor_Z = "Gender",
    x_covariate_prior_sd = 10, x_covariate_in_outcome = FALSE,
    K_func = 12, basis_weights = "unit",
    K_beta = 10, beta_bs = "cc", periodic_time = FALSE,
    sigma_gamma_prior_mean = 20,
    sigma_gamma_prior_sd = 0.5,
    sigma_gamma_lower = 18,
    sigma_gamma_upper = 21,
    a_x_prior_sd = 100,
    sigma_omega_prior_scale = 15,
    sigma_e_prior_scale = 10,
    beta_fixed_prior_mean = 0.07642642,
    beta_fixed_prior_sd = 0.01460296,
    beta_random_prior_scale = 1,
    chains = 2, cores = 2, seed = 1123
  )

  attr(fit_joint01_betaf_loo, "beta_prior")
  attr(fit_joint05_betaf_loo, "beta_prior")
  attr(fit_joint01_betaf_loo, "loo_formula")
  attr(fit_joint05_betaf_loo, "loo_formula")
}
