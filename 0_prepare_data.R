# Prepare common Y/Z and Joint, FUI, and Raw functional inputs.
# Run in a fresh R session. Required packages: lme4, dplyr, mgcv.
# Warray must be ordered exactly as df_sample$SEQN in GAL.RData.
# This preserves the original FUI calibration population (all of Warray).
# No outcome model fitting, scaling, PA summaries, or BMI outlier exclusion is done here.


library(emmeans)
library(kableExtra)
library(papeR)
library(table1)
library(tidyr)
library(dplyr)
library(lubridate)
library(Hmisc)
library(consort)
library(readxl)
library(ggplot2)
library(gtsummary)
library(data.table)
library(lavaan) 
library(labelled)
library(lmtest)
library(car)
library(brms)
library(rstan)
library(splines)
library(readr)
library(plyr)
library(mombf)
library(LaplacesDemon)
library(data.table)



load("/Users/yangou/Documents/GitHub/BQSOFR/GAL.RData")

# MEM/FUI function, retained without algorithm changes.
FUI = function(model="gaussian",smooth=FALSE, data,silent = FALSE){
  
  ##' @param model is the model for apporixation 
  ##   model = gaussian, poisson 
  ##' @param data is the observed measurement of X(t) with repeated measures, which is in 3 dimension [n, t, m_w]
  ##' @param smooth whether conduct the smoothing step or not (Step II)
  ##' @param silent whether to show descriptions of each step
  
  ## The output of this function is the approximation of the true function covariate X(t),
  ##    where each column represent each time point and each row represent each subject
  
  
  
  library(lme4) ## mixed models
  #library(refund) ## fpca.face
  library(dplyr) ## organize lapply results
  #library(progress) ## display progress bar
  library(mgcv) ## smoothing in step 2
  #library(mvtnorm) ## joint CI
  #library(parallel) ## mcapply
  
  ### create a dataframe for analysis 
  m_w = dim(data)[3] ## number of repeated measures
  t = dim(data)[2] ## number of observations on the functional domain (time points)
  n = dim(data)[1] ## number of subjects 
  M = NULL
  i=1
  while(i <= m_w){
    M = rbind(M,data[,,i])
    i=i+1
  } ### convert the 3d array to 2d matrix 
  
  data.M = data.frame(M = I(M),seqn= rep(1:n, m_w) )
  
  ##########################################################################################
  ## Step 1 (Massive Univariate analysis)
  ##########################################################################################
  if(silent == FALSE) print("Step 1: Massive Univariate Mixed Models")
  
  if(model =="gaussian"){
    fit = data.M %>%
      dplyr::select(M) %>%
      apply( 2, function(s) {
        temp= data.frame(M=s, seqn=data.M$seqn)
        mod =  suppressMessages(lmer(M ~ (1 | seqn), data = temp, na.action = na.exclude))
        return(predict(mod))
      })
  }else{
    fit = data.M %>%
      dplyr::select(M) %>%
      apply( 2, function(s) {
        temp= data.frame(M=s, seqn=data.M$seqn)
        mod =  suppressMessages(glmer(M ~ (1 | seqn), data = temp, family = model))
        return(predict(mod))
      })
    
  }
  
  approx= fit[!duplicated(data.M$seqn),] ## each column represents one time point 
  
  ##########################################################################################
  ## Step 2 (Smoothing)
  ##########################################################################################
  
  if(smooth==TRUE){
    if(silent == FALSE) print("Step 2: Smoothing")
    nknots <- min(round(t/4), 35) ## number of knots for penalized splines smoothing
    argvals = seq(0,1, length.out = t) ## locations of observations on the functional domain
    approx= t(apply(approx, 1, function(x) gam(x ~ s(argvals, bs = "cr", k = (nknots + 1)), method = "REML")$fitted.values))
  }
  
  return(approx)
}

# 1. Validate source objects and subject IDs.
df_base <- as.data.frame(df_sample)
hei_base <- as.data.frame(HEI_2013_2014)

numeric_vars_Z <- c("AgeYR", "HEI")
factor_vars_Z <- c("Gender", "Race", "HealthCondt2")

stopifnot(
  all(c("SEQN", "BMI", "AgeYR", factor_vars_Z) %in% names(df_base)),
  "SEQN" %in% names(hei_base),
  length(dim(Warray)) == 3L,
  is.numeric(Warray),
  dim(Warray)[1] == nrow(df_base),
  dim(Warray)[1] >= 2L,
  dim(Warray)[2] >= 4L,
  dim(Warray)[3] >= 2L
)

id_W <- df_base$SEQN
if (anyNA(id_W) || anyDuplicated(id_W)) {
  stop("df_sample$SEQN must be nonmissing and unique.")
}
if (anyNA(hei_base$SEQN) || anyDuplicated(hei_base$SEQN)) {
  stop("HEI_2013_2014$SEQN must be nonmissing and unique.")
}

# 2. Add HEI 
hei_score_column <- 16L
stopifnot(ncol(hei_base) >= hei_score_column)
message("HEI source column: ", names(hei_base)[hei_score_column])

idx_HEI <- match(id_W, hei_base$SEQN)
has_HEI_record <- !is.na(idx_HEI)
df_sampleJt00 <- df_base[has_HEI_record, , drop = FALSE]
df_sampleJt00$HEI <- hei_base[[hei_score_column]][idx_HEI[has_HEI_record]]

# 3. Match functional rows and define the common analysis sample.
df_sampleWt0 <- df_sampleJt00
Y <- df_sampleWt0$BMI
Z_final <- df_sampleWt0[, c(numeric_vars_Z, factor_vars_Z), drop = FALSE]
id_final <- df_sampleWt0$SEQN
idx_W <- match(id_final, id_W)
stopifnot(!anyNA(idx_W))

if (!is.numeric(Y) || !all(vapply(Z_final[numeric_vars_Z], is.numeric, logical(1)))) {
  stop("BMI, AgeYR and HEI must be numeric; check their original coding.")
}

idx <- complete.cases(Y, Z_final)
if (!any(idx)) stop("No subjects remain after complete-case selection.")

# Any separately justified exclusion rule must be added to idx BEFORE this step.
Y_use <- Y[idx]
Z_use <- Z_final[idx, , drop = FALSE]
id_use <- id_final[idx]
rows_use <- idx_W[idx]

# 4. Joint: retain ALL repeated days. Dimensions N x T x J.
W_use <- Warray[rows_use, , , drop = FALSE]

if (!all(is.finite(Y_use)) ||
    !all(vapply(Z_use[numeric_vars_Z], function(x) all(is.finite(x)), logical(1)))) {
  stop("Non-finite BMI, AgeYR or HEI values remain; inspect them before fitting.")
}
if (!all(is.finite(W_use))) {
  stop("W_use contains NA/NaN/Inf. Define a common missing-data rule before proceeding.")
}

# 5. Raw: average over the THIRD dimension (days). Dimensions N x T.
W_raw <- apply(W_use, c(1, 2), mean)

# 6. MEM/FUI

w_fui <- as.matrix(FUI(
  model = "gaussian", smooth = T, data = Warray, silent = FALSE
))



stopifnot(identical(dim(w_fui), dim(Warray)[1:2]))
w_fui_use <- w_fui[rows_use, , drop = FALSE]
if (!all(is.finite(w_fui_use))) stop("FUI returned non-finite values.")

# 7. Attach IDs and check common dimensions.
# These labels record the alignment above; they cannot prove the source row order.
names(Y_use) <- as.character(id_use)
rownames(Z_use) <- as.character(id_use)
rownames(W_raw) <- as.character(id_use)
rownames(w_fui_use) <- as.character(id_use)
dimnames(W_use)[[1]] <- as.character(id_use)

stopifnot(
  length(Y_use) == nrow(Z_use),
  length(Y_use) == dim(W_use)[1],
  length(Y_use) == nrow(W_raw),
  length(Y_use) == nrow(w_fui_use),
  dim(W_use)[2] == ncol(W_raw),
  dim(W_use)[2] == ncol(w_fui_use),
  identical(as.character(id_use), as.character(id_W[rows_use]))
)

print(list(
  N_Y = length(Y_use),
  Z_dim = dim(Z_use),
  W_joint_dim = dim(W_use),
  W_raw_dim = dim(W_raw),
  W_fui_dim = dim(w_fui_use),
  N_source = nrow(df_base),
  N_no_HEI_record = sum(!has_HEI_record),
  N_incomplete_Y_Z = sum(!idx)
))

# Existing model calls keep their original arguments:
# Joint: Y = Y_use, Z = Z_use, W = W_use
# FUI:   Y = Y_use, Z = Z_use, W = w_fui_use
# Raw:   Y = Y_use, Z = Z_use, W = W_raw
