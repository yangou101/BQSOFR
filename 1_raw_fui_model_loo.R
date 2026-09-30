# FUI / Raw: second-stage Bayesian functional quantile regression.
#
# Retained functions:
#   BQBayes.MCMCSampleSFPCA_raw_fui(): fit one supplied N x T matrix.
#   .rfc_beta(): reconstruct beta(t) on the original Y scale, with
#                posterior means and pointwise 95% credible intervals.
#
# Inputs must already use the same complete subjects and time-column order:
#   FUI: W = w_fui_use (already estimated and aligned FUI curves).
#   Raw: W = W_raw     (daily means prepared from the N x T x J array).
# This file does not perform first-stage FUI estimation or data alignment.
# Pass scalar_formula = formula01 for the desired main effects/interactions.
# scalar_formula = NULL uses main effects only.



BQBayes.MCMCSampleSFPCA_raw_fui <- function(
  Y, W, Z, a = NULL, Nsim,
  tau0 = 0.9, numeric_Z, factor_Z,
  seed = 1123, chains = 2, cores = 2,
  warmup = floor(Nsim / 2), adapt_delta = 0.99, max_treedepth = 10,
  scalar_formula = NULL
){

  .rfc_require(c("rstan", "brms", "refund", "mgcv"))
  if (!is.matrix(W) || !is.numeric(W)) stop("W must be a numeric N x T matrix.")
  if (nrow(W) != length(Y) || nrow(W) != nrow(Z)) stop("Y, Z and W have different subject counts.")
  if (any(!is.finite(Y)) || any(!is.finite(W)) || any(!complete.cases(Z))) {
    stop("Use the same complete cohort for both fits before calling this function.")
  }
  if (!is.finite(sd(Y)) || sd(Y) <= 1e-9) stop("Y must have positive variance.")
  if (Nsim <= warmup || warmup < 1 || Nsim != as.integer(Nsim)) stop("Invalid Nsim/warmup.")
  if (tau0 <= 0 || tau0 >= 1) stop("tau0 must be strictly between 0 and 1.")
  cat("\n#--- Initial values...\n")

  GamF <- function(gam, p0){
    (2*pnorm(-abs(gam))*exp(.5*gam^2) - p0)^2
  }
  GamBnd <- function(p0){
    Re1 <- optimize(GamF, interval = c(-30, 30), p0 = 1-p0)
    Re2 <- optimize(GamF, interval = c(-30, 30), p0 = p0)
    c(-abs(Re1$minimum), abs(Re2$minimum), Re1$objective, Re2$objective)
  }

  Y <- as.numeric(Y)
  Y_original_loo <- Y
  val_mean <- mean(Y, na.rm = TRUE)
  val_sd   <- sd(Y,   na.rm = TRUE)
  if (is.finite(val_sd) && val_sd > 1e-9) Y <- (Y - val_mean)/val_sd else Y <- Y - val_mean

  Z_df <- as.data.frame(Z)
  if(length(factor_Z) > 0){
    for(nm in factor_Z) Z_df[[nm]] <- as.factor(Z_df[[nm]])
  }

  for (nm in numeric_Z) {
    if (!is.numeric(Z_df[[nm]]) || any(!is.finite(Z_df[[nm]])) ||
        !is.finite(sd(Z_df[[nm]])) || sd(Z_df[[nm]]) <= 0) {
      stop("Numeric covariate must be finite and have positive SD: ", nm)
    }
  }
  Z_numeric <- Z_df[, numeric_Z, drop = FALSE]
  Z_factor  <- Z_df[, factor_Z,  drop = FALSE]
  Z_final   <- Z_factor
  if (ncol(Z_numeric) > 0) {
    Z_numeric_scaled <- as.data.frame(scale(Z_numeric, center = TRUE, scale = TRUE))
    Z_final <- data.frame(Z_numeric_scaled, Z_factor)
  }

  Xt<-W

  n_num <- nrow(Xt)
  nt    <- ncol(Xt)

  tind <- seq(0, 1, length.out = nt)

  idx <- complete.cases(Y, Z_final) & apply(is.finite(Xt), 1, all)
  if (!all(idx)) {
    Y      <- Y[idx]
    Z_final <- Z_final[idx, , drop = FALSE]
    Xt     <- Xt[idx, , drop = FALSE]
    n_num  <- nrow(Xt)
  }

  cat("\n#--- FPCA...\n")
  
  fpca.fit <- refund::fpca.sc(Xt)

  J_num<-dim(fpca.fit$efunctions)[2]

  Xi_fix <- fpca.fit$scores

  wmat <- I(Xt) # predict  matrix Xi(t)
  lmat <- I(matrix(1/nt, ncol = nt, nrow = n_num))#The Integration Weights
  tmat <- I(matrix(tind, ncol = nt, nrow = n_num, byrow = TRUE))

  dat_list <- list(tmat = tmat, lmat = lmat, wmat = wmat)

  knots <- NULL

  object <- mgcv::s(tmat, by = lmat * wmat, bs = "cc", k = 10)

  dk <- mgcv:::ExtractData(object, data = dat_list, knots = knots)
  splinecons <- mgcv:::smooth.construct.cc.smooth.spec(object, dk$data, dk$knots)

  Psi_mat <- splinecons$X # Time-domain beta basis before penalty transformation.
  if (nrow(Psi_mat) != nt || nrow(fpca.fit$efunctions) != nt) {
    stop("The original mgcv/FPCA construction did not return one basis row per time point.")
  }
  S_mat   <- splinecons$S[[1]] #The Penalty Matrix
  rank    <- splinecons$rank #degrees of freedom subject to the smoothing penalty

  M_num <- dim(splinecons$X)[1]
  J_num = dim(fpca.fit$efunctions)[2]
  K_num <- dim(splinecons$X)[2]
  
  X_mat_t = matrix(nrow=J_num, ncol= K_num)
for(j in 1:J_num){
  for(k in 1:K_num) {
    X_mat_t[j,k] = sum(fpca.fit$efunctions[,j] * Psi_mat[,k]) / M_num
  }
}

  maXX <- norm(Psi_mat, type="I")^2
  maS  <- norm(S_mat,  type="I") / maXX
  S_mat <- S_mat / maS

  eig <- eigen(S_mat, symmetric = TRUE)

  E <- rep(1, ncol(X_mat_t))
  E[1:rank] <- sqrt(pmax(eig$values[1:rank], 1e-12))

  X_mat_t <- X_mat_t %*% eig$vectors

  col.norm <- colSums(X_mat_t^2)
  col.norm <- col.norm / (E^2)
  av.norm  <- mean(col.norm[1:rank])

  if (rank < ncol(X_mat_t)) {
    for (i in (rank + 1):ncol(X_mat_t)) {
      E[i] <- sqrt(col.norm[i] / av.norm)
    }
  }

  if (any(!is.finite(E)) || any(E <= 0)) {
    stop("A beta basis direction has zero or invalid linkage to the FPCA space.")
  }
  X_mat_t <- t(t(X_mat_t) / E)

  B_penalty <- sweep(Psi_mat %*% eig$vectors, 2, E, "/")
  beta_order <- c(seq.int(rank + 1L, K_num), seq_len(rank))
  B_beta_compare <- B_penalty[, beta_order, drop = FALSE]
  G_direct <- crossprod(fpca.fit$efunctions, B_beta_compare) / M_num
  G_expected <- X_mat_t[, beta_order, drop = FALSE]
  if (max(abs(G_direct - G_expected)) > 1e-8 * max(1, max(abs(G_expected)))) {
    stop("Beta reconstruction is inconsistent with the fitted functional term.")
  }

  X_mat_r <- t(X_mat_t[, 1:rank, drop = FALSE])                 # Kr x J. Random Effects Part
  X_mat_f <- t(X_mat_t[, (rank + 1):ncol(X_mat_t), drop = FALSE]) # Kf x J. Fixed Effects Part

  Bd <- round(GamBnd(tau0)[1:2] * 0.99, 4)

  StandCodeVGalLin <- "
functions {

  real calculate_r_gamma(real gamma) {
    real log_r = log(2) + normal_lcdf(-abs(gamma) | 0, 1) + (square(gamma) / 2.0);
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
      log_abs_part1 = log_cdf_arg2 + log1m_exp(log_cdf_arg1 - log_cdf_arg2)
                      + ((-pg_minus * estar) + ((square(gamma) / 2.0) * square(ratio_mp)))
                      + log((estar / gamma) > 0 ? 1.0 : 0.0);
      sgn_part1 = -1.0;
    } else {
      log_abs_part1 = log_cdf_arg1 + log1m_exp(log_cdf_arg2 - log_cdf_arg1)
                      + ((-pg_minus * estar) + ((square(gamma) / 2.0) * square(ratio_mp)))
                      + log((estar / gamma) > 0 ? 1.0 : 0.0);
      sgn_part1 = 1.0;
    }

    real ind_term = ((estar / gamma) > 0) ? 1.0 : 0.0;
    real arg3 = -abs_g + ((pg_plus * estar * ind_term) / abs_g);
    real log_part2 = normal_lcdf(arg3 | 0, 1) + ((-pg_plus * estar) + (square(gamma) / 2.0));

    real val;

    if (sgn_part1 == 1.0) {
      val = log(2.0) + log(p) + log(1.0 - p) - log(sigma)
            + log_sum_exp(log_abs_part1, log_part2);
    } else {
      if (log_part2 > log_abs_part1) {
        val = log(2.0) + log(p) + log(1.0 - p) - log(sigma)
              + log_part2 + log1m_exp(log_abs_part1 - log_part2);
      } else {
        val = negative_infinity();
      }
    }
    return val;
  }

}

data {
  int<lower=1> N;
  vector[N] Y;
  int<lower=1> K;
  matrix[N, K] X;
  int prior_only;

  vector[2] Bd;
  real<lower=0,upper=1> tau0;

  int<lower=1> J_num;
  int<lower=1> Kf;
  int<lower=1> Kr;

  matrix[N, J_num] Xi_fix;
  matrix[Kf, J_num] X_mat_f;
  matrix[Kr, J_num] X_mat_r;
}

transformed data {
  int Kc = K - 1;
  matrix[N, Kc] Xc;
  vector[Kc] means_X;

  for (i in 2:K) {
    means_X[i - 1] = mean(X[, i]);
    Xc[, i - 1] = X[, i] - means_X[i - 1];
  }
}

parameters {
  vector[Kc] b;
  real Intercept;
  //real<lower=0> sigma;
  real log_sigma; // added this line u=11 可以不要

  vector[Kf] betaf;
  real<lower=0> sigma_betar2;
  vector[Kr] z_betar;

  real<lower=0, upper=1> U;
}

transformed parameters {

real<lower=0> sigma = exp(log_sigma); // added this line u=11 可以不要

  real<lower=0> sigma_betar = sqrt(sigma_betar2);
  real gam;
  vector[Kr] betar;

  gam = (Bd[2] - Bd[1]) * U + Bd[1];
  betar = sigma_betar * z_betar;
}

model {
  real lprior = 0;
  vector[J_num] b_func;
  vector[N] mu;

  lprior += student_t_lpdf(Intercept | 3, 0, 2.5);
 // lprior += student_t_lpdf(sigma | 3, 0, 2.5)
   //         - student_t_lccdf(0 | 3, 0, 2.5);
  lprior += normal_lpdf(log_sigma | log(0.5), 1.0);         

  lprior += normal_lpdf(b | 0, 2.5);
  //lprior += beta_lpdf(U | 1, 1);
  lprior += beta_lpdf(U | 4, 4);
  lprior += normal_lpdf(betaf | 0, 1);
  lprior += std_normal_lpdf(z_betar);
  lprior += inv_gamma_lpdf(sigma_betar2 | 0.001, 0.001);

  target += lprior;

  if (!prior_only) {
    b_func = (X_mat_f' * betaf) + (X_mat_r' * betar);
    mu = Intercept + Xc * b + Xi_fix * b_func;

    for (n in 1:N) {
      target += gal_p0_lpdf(Y[n] | mu[n], sigma, gam, tau0);
    }
  }
}

generated quantities {
  real b_Intercept = Intercept - dot_product(means_X, b);
  vector[N] log_lik_y;
  {
    // Exactly the outcome predictor and normalized GAL density in model {}.
    vector[J_num] b_func_y = (X_mat_f' * betaf) + (X_mat_r' * betar);
    vector[N] mu_y = Intercept + Xc * b + Xi_fix * b_func_y;
    for (n in 1:N) {
      log_lik_y[n] = gal_p0_lpdf(Y[n] | mu_y[n], sigma, gam, tau0);
    }
  }
}
"

  df2 <- data.frame(Y = Y, Z_final)
  form_Z <- if (is.null(scalar_formula)) {
    stats::reformulate(colnames(Z_final), response = "Y")
  } else {
    stats::as.formula(scalar_formula)
  }
  if (!identical(form_Z[[2L]], quote(Y)) ||
      !all(all.vars(form_Z) %in% names(df2))) {
    stop("scalar_formula must have the untransformed response Y and only selected Z variables.")
  }
  formula_terms <- stats::terms(form_Z)
  rhs_calls <- setdiff(all.names(form_Z[[3L]], functions = TRUE, unique = TRUE),
                       all.vars(form_Z[[3L]]))
  if (attr(formula_terms, "intercept") != 1L ||
      length(attr(formula_terms, "offset")) ||
      !all(rhs_calls %in% c("+", "-", "*", ":", "^", "/", "("))) {
    stop("scalar_formula supports an intercept, ordinary scalar main effects and interactions; offsets and special terms are not supported.")
  }
  cat("\nScalar outcome formula:\n")
  print(form_Z)

  Dt2 <- brms::make_standata(
    brms::bf(as.formula(form_Z), quantile = tau0),
    data   = df2,
    family = brms::asym_laplace()
  )
  if (Dt2$K < 2L || !all(Dt2$X[, 1L] == 1) ||
      nrow(Dt2$X) != length(Y) || any(!is.finite(Dt2$X))) {
    stop("The scalar design must retain all subjects, with an intercept in column 1 and at least one scalar predictor.")
  }

  Dt2$Bd    <- Bd
  Dt2$tau0  <- tau0

  Dt2$Kf      <- NROW(X_mat_f)
  Dt2$Kr      <- NROW(X_mat_r)
  Dt2$X_mat_f <- X_mat_f
  Dt2$X_mat_r <- X_mat_r
  
  Dt2$J_num   <- J_num
  Dt2$M_num   <- M_num
  
  Xt_c <- sweep(Xt, 2, fpca.fit$mu, FUN = "-")   # N x T

  Dt2$W_mat <- Xt_c
  Dt2$T_num <- ncol(Xt_c)

  Dt2$Xi_fix  <- Xi_fix

vec1 <- function(n, val = 0) array(val, dim = as.integer(n))
mat1 <- function(nr, nc, val = 0) array(val, dim = c(as.integer(nr), as.integer(nc)))

init_fun <- function() {
  Kc_now <- as.integer(Dt2$K - 1)

  list(
    U = 0.5,
    log_sigma = log(0.5),
    sigma_betar2 = 0.2^2,
    betaf   = vec1(Dt2$Kf, 0),
    b       = vec1(Kc_now, 0),
    z_betar = vec1(Dt2$Kr, 0),
    Intercept = 0
  )
}

fit_tmp <- rstan::stan(
  model_code = StandCodeVGalLin,
  data  = Dt2,
  iter  = Nsim,
  warmup = warmup,
  chains = chains,
  refresh = 50,
  init = init_fun,
  control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
  cores = cores,
  seed  = seed
)

  attr(fit_tmp, "B_beta_compare") <- B_beta_compare
  attr(fit_tmp, "time_grid_compare") <- tind
  attr(fit_tmp, "Y_mean_compare") <- val_mean
  attr(fit_tmp, "Y_sd_compare") <- val_sd
  attr(fit_tmp, "fpca_components_compare") <- J_num
  attr(fit_tmp, "fpca_evalues_compare") <- fpca.fit$evalues
  attr(fit_tmp, "beta_Kf_compare") <- Dt2$Kf
  attr(fit_tmp, "beta_Kr_compare") <- Dt2$Kr
  attr(fit_tmp, "scalar_names_compare") <- colnames(Dt2$X)[-1L]
  attr(fit_tmp, "tau_compare") <- tau0
  attr(fit_tmp, "stage2_priors_compare") <- c(
    betaf_SD = 1, sigma_betar2_IG_shape = 0.001,
    sigma_betar2_IG_scale = 0.001
  )
  # Metadata makes cohort, scale and scalar-design checks possible before LOO.
  attr(fit_tmp, "loo_y") <- Y_original_loo[idx]
  attr(fit_tmp, "loo_input_row") <- which(idx)
  attr(fit_tmp, "loo_y_mean") <- val_mean
  attr(fit_tmp, "loo_y_sd") <- val_sd
  attr(fit_tmp, "loo_tau") <- tau0
  attr(fit_tmp, "loo_X") <- Dt2$X
  attr(fit_tmp, "loo_formula") <- paste(deparse(form_Z), collapse = " ")
  attr(fit_tmp, "loo_log_lik_scale") <- "standardized_Y"
  attr(fit_tmp, "loo_method") <- "FUI"
  attr(fit_tmp, "max_treedepth") <- max_treedepth
  return(fit_tmp)
}

# Posterior beta on ORIGINAL Y scale, using each fitted model's exact basis.
.rfc_beta <- function(fit, method) {
  dr <- rstan::extract(fit, pars = c("betaf", "betar"), permuted = TRUE)
  Kf <- attr(fit, "beta_Kf_compare")
  Kr <- attr(fit, "beta_Kr_compare")
  cf <- matrix(dr$betaf, ncol = Kf)
  cr <- matrix(dr$betar, ncol = Kr)
  B <- attr(fit, "B_beta_compare")
  stopifnot(nrow(cf) == nrow(cr), ncol(B) == Kf + Kr)
  beta_draws <- (cbind(cf, cr) %*% t(B)) * attr(fit, "Y_sd_compare")
  qs <- apply(beta_draws, 2, quantile, probs = c(0.025, 0.5, 0.975))
  tg <- attr(fit, "time_grid_compare")
  data.frame(method = method, t = tg, hour = 24 * tg,
    mean = colMeans(beta_draws), lower = qs[1L, ],
    median = qs[2L, ], upper = qs[3L, ])
}
