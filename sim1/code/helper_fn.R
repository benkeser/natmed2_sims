#############################################################################
#                               HELPER FUNCTIONS                            #
#############################################################################

#' Generate observed data for observational simulation to confirm 
#' theoretical properties of estimators.
#' 
#' There are two binary covariates (\code{W1}, \code{W2}), 
#' a binary treatment (\code{A}), a binomial mediator (\code{S}),
#' an indicator of having outcome measured (\code{C}), a 
#' binary outcome (\code{Y}), and a two-phase sampling indicator (\code{R}).
#' Individuals with missing \code{Y} or \code{S} have their 
#' missing values replace with \code{-999}.
#' 
#' @param n Sample size
#' @param p1 Probability W1 = 1
#' @param p2 Probability W2 = 1
#' @return A \code{data.frame} with the observed data.

sim_data = function(n, p1 = 0.5, p2 = 0.5){
  # generate W(covariates)
  W1 = rbinom(n, 1, p1)
  W2 = rbinom(n, 1, p2)
  # generate A(vaccine or not)
  A = rbinom(n, 1, plogis(W1 - W2))
  # generate S(immune response)
  S = rbinom(n, 2, plogis(-1 + W1 / 4  - W2 / 3 + A / 2))
  # generate Y(infect or not)
  Y = rbinom(n, 1, plogis(-2 + A / 2 + W1 / 2 - S / 2))
  
  # add censoring
  C = rbinom(n, 1, plogis(2 + W1 / 2 - W2 / 3))
  # arbitrary fill in
  Y[C == 0] = -999
  
  # case-cohort sampling
  R = rbinom(n, 1, 0.25)
  R[Y == 1] = 1
  
  # remove S from R == 0 people
  S[R == 0] = -999
  
  return(data.frame(W1 = W1, W2 = W2, A = A, S = S, Y = Y, C = C, R = R))
}

#' Compute the true value of E[Y(1, S(0))] based on the data
#' generating process encoded by \code{sim_data}.
#' 
#' @param n Sample size (should be very large to have accurate approximation)
#' @param p1 Probability W1 = 1
#' @param p2 Probability W2 = 1
#' @param return A numeric
get_truth = function(n = 1e6, p1 = 0.5, p2 = 0.5){
  # W(covariates)
  W1 = rbinom(n, 1, p1)
  W2 = rbinom(n, 1, p2)
  # S(immune response)
  # S1: under vaccine; S0: under placebo
  S1 = rbinom(n, 2, plogis(-1 + W1 / 4  - W2 / 3 + 1 / 2))
  S0 = rbinom(n, 2, plogis(-1 + W1 / 4  - W2 / 3 + 0 / 2))

  # generate Y(infect or not)
  # Y(1,S(1)), Y(1,S(0)), Y(0,S(1)), Y(0,S(0))
  Y1S1 = rbinom(n, 1, plogis(-2 + 1 / 2 + W1 / 2 - S1 / 2))
  Y1S0 = rbinom(n, 1, plogis(-2 + 1 / 2 + W1 / 2 - S0 / 2))
  Y0S1 = rbinom(n, 1, plogis(-2 + 0 / 2 + W1 / 2 - S1 / 2))
  Y0S0 = rbinom(n, 1, plogis(-2 + 0 / 2 + W1 / 2 - S0 / 2))
  return(c(
    mean(Y1S0)
  ))
}

#' Compute the variance of the efficient influence function of
#' E[Y(1, S(0))].
#' 
#' The function simulates a large data set and uses either known (in 
#' case of outcome regression and treatment/censoring probabilities)
#' nuisance parameter values or consistent estimates thereof (in the case
#' of the sequential regression quantities and the mean of the full
#' data efficient influence function).
#' 
#' @param n Sample size (should be very large to have accurate approximation)
#' @param p1 Probability W1 = 1
#' @param p2 Probability W2 = 1
#' @return A numeric

var_eif_ey1s0 <- function(n = 1e6, p1 = 0.5, p2 = 0.5){
  # generate W(covariates)
  W1 = rbinom(n, 1, p1)
  W2 = rbinom(n, 1, p2)
  # generate A(vaccine or not)
  A = rbinom(n, 1, plogis(W1 - W2))
  # generate S(immune response)
  S = rbinom(n, 2, plogis(-1 + W1 / 4  - W2 / 3 + A / 2))
  # generate Y(infect or not)
  Y = rbinom(n, 1, plogis(-2 + A / 2 + W1 / 2 - S / 2))
  
  # add censoring
  C = rbinom(n, 1, plogis(2 + W1 / 2 - W2 / 3))
  # arbitrary fill in
  Y[C == 0] = -999
  
  # case-cohort sampling
  R = rbinom(n, 1, 0.25)
  R[C * Y == 1] = 1
  
  # remove S from R == 0 people
  # S[R == 0] = -999

  Qbar0 <- plogis(-2 + 1 / 2 + W1 / 2 - S / 2)
  gAW0 <- plogis(W1 - W2)
  gC0 <- plogis(2 + W1 / 2 - W2 / 3)
  gR0 <- rep(0.25, n); gR0[C * Y == 1] <- 1
  gAWS0 <- rep(NA, n)
  Qbarbar0 <- rep(NA, n)
  w1_vec <- w2_vec <- s_vec <- gAWS_vec <- Qbarbar0_vec <- Qbar0_vec <- gAW0_vec <- gC0_vec <- NULL
  for(w1 in c(0,1)){
    for(w2 in c(0,1)){
      Qbarbar0[W1 == w1 & W2 == w2] <- mean(Qbar0[A == 0 & W1 == w1 & W2 == w2])
      for(s in 0:2){
        gAWS0[W1 == w1 & W2 == w2 & S == s] <- mean(A[W1 == w1 & W2 == w2 & S == s])
        w1_vec <- c(w1_vec, w1)
        w2_vec <- c(w2_vec, w2)
        s_vec <- c(s_vec, s)
        gAWS_vec <- c(gAWS_vec, mean(A[W1 == w1 & W2 == w2 & S == s]))
        Qbarbar0_vec <- c(Qbarbar0_vec, mean(Qbar0[A == 0 & W1 == w1 & W2 == w2]))
        Qbar0_vec <- c(Qbar0_vec, mean(Y[A == 1 & C == 1 & W1 == w1 & W2 == w2 & S == s]))
        gAW0_vec <- c(gAW0_vec, mean(A[W1 == w1 & W2 == w2]))
        gC0_vec <- c(gC0_vec, mean(C[W1 == w1 & W2 == w2]))
      }
    }
  }
  np_df <- data.frame(w1 = w1_vec, w2 = w2_vec, s = s, gAW = gAW0_vec, gAWS = gAWS_vec, gC = gC0_vec,
                      Qbar = Qbar0_vec, Qbarbar = Qbarbar0_vec)

  psi <- mean(Qbarbar0)

  full_data_eif <- as.numeric(A == 1 & C == 1) / (gAW0 * gC0) * gAW0 / (1 - gAW0) * (1 - gAWS0) / gAWS0 * (Y - Qbar0) + 
                    as.numeric(A == 0) / (1 - gAW0) * (Qbar0 - Qbarbar0) + 
                      Qbarbar0 - psi

  mean_full_data_eif <- rep(NA, n)
  blah <- w1_vec <- w2_vec <- a_vec <- c_vec <- y_vec <- NULL                  
  for(w1 in c(0,1)){
    for(w2 in c(0,1)){
      for(a in c(0,1)){
          a_vec <- c(a_vec, a)
          w1_vec <- c(w1_vec, w1)
          w2_vec <- c(w2_vec, w2)
          c_vec <- c(c_vec, 0)
          y_vec <- c(y_vec, -999)
          mean_full_data_eif[W1 == w1 & W2 == w2 & A == a & C == 0] <- mean(full_data_eif[R == 1 & W1 == w1 & W2 == w2 & A == a & C == 0])
          blah <- c(blah, mean(full_data_eif[R == 1 & W1 == w1 & W2 == w2 & A == a & C == 0]))
          a_vec <- c(a_vec, a)
          w1_vec <- c(w1_vec, w1)
          w2_vec <- c(w2_vec, w2)
          c_vec <- c(c_vec, 1)
          y_vec <- c(y_vec, 1)
          mean_full_data_eif[W1 == w1 & W2 == w2 & A == a & C == 1 & Y == 1] <- mean(full_data_eif[R == 1 & W1 == w1 & W2 == w2 & A == a & C == 1 & Y == 1])
          blah <- c(blah, mean(full_data_eif[R == 1 & W1 == w1 & W2 == w2 & A == a & C == 1 & Y == 1]))
          a_vec <- c(a_vec, a)
          w1_vec <- c(w1_vec, w1)
          w2_vec <- c(w2_vec, w2)
          c_vec <- c(c_vec, 1)
          y_vec <- c(y_vec, 0)
          mean_full_data_eif[W1 == w1 & W2 == w2 & A == a & C == 1 & Y == 0] <- mean(full_data_eif[R == 1 & W1 == w1 & W2 == w2 & A == a & C == 1 & Y == 0])
          blah <- c(blah, mean(full_data_eif[R == 1 & W1 == w1 & W2 == w2 & A == a & C == 1 & Y == 0]))
        }
      }
    }
    np_df2 <- data.frame(w1 = w1_vec, w2 = w2_vec, s = s_vec, a = a_vec, 
                         c = c_vec, y = y_vec, QD = blah)

  obs_data_eif <- R / gR0 * full_data_eif + (1 - R / gR0) * mean_full_data_eif
  return(var(obs_data_eif))
}


# simulate the data by setting different seeds
# X is a data.frame or list where 
#   $seed is the seed to set
#   $setting is the regression formula setting
#   $n is the sample size

#' Run one iteration of simulation for a specified set of parameters
#' 
#' Given parameter specification \code{X} (see below for description), 
#' draw a sample of data and compute estimators of E[Y(1, S(0))].
#' 
#' The parameter specification is a list with named elements \code{n} (sample size) 
#' and \code{setting}. The \code{setting} should itself be a list with 
#' name elements providing the regression formula for the various nuisance
#' quantities. Lists that hold the settings considered in the simulation can be created 
#' using \code{create_setting} function.
#' 
#' @param X A list with named elements \code{n} (sample size) and \code{setting} (described above).
#' @param p1 Probability W1 = 1
#' @param p2 Probability W2 = 1
#' @return A six length numeric. The first three entries are respectively, the point estimate
#' and 95% confidence interval for E[Y(1, S(0))]; the second three entries are the same, but the
#' lazy version.

one_sim = function(X, p1 = 0.5, p2 = 0.5){
  dat = sim_data(n = X$n, p1 = p1, p2 = p2)
  
  fit = natmed2(
    W = data.frame(W1 = dat$W1, W2 = dat$W2), 
    A = dat$A, R = dat$R, S = dat$S, C = dat$C, Y = dat$Y,
    glm_gR = X$settings$glm_gR,
    glm_gA = X$settings$glm_gA,
    glm_gC = X$settings$glm_gC,
    glm_gAS = X$settings$glm_gAS,
    glm_QY_WAS = X$settings$glm_QY_WAS,
    glm_QY_WACY = X$settings$glm_QY_WACY,
    glm_QY_W = X$settings$glm_QY_W,
    glm_QY_WA = X$settings$glm_QY_WA,
    glm_QD_WACY = X$settings$glm_QD_WACY,
    glm_QD_WACY_lazy = X$settings$glm_QD_WACY_lazy
  )
  return(as.numeric(c(fit$risk[3,1:3], fit$risk_lazy[3,1:3])))
}

#' Define a simulation setting where some nuisance quantities are inconsistently estimated
#' 
#' If a nuisance quantity is not included in \code{which_wrong}, the proper specification is 
#' given for the regression function. If a nuisance quantity is included in \code{which_wrong},
#' an intercept-only model is used instead. 
#' 
#' @param which_wrong A character vector indicating which nuisance regressions are inconsistently estimated.
#' @return A named list that can be passed to \code{one_sim} function to run a simulation 
#' under a certain pattern of nuisance parameter (mis)specification.
create_setting <- function(which_wrong){
  list(
    glm_gR = ifelse(!("gR" %in% which_wrong), "CY11", "1"),
    glm_gC = ifelse(!("gC" %in% which_wrong), "W1 + W2", "1"),
    glm_gAS = ifelse(!("gAS" %in% which_wrong), "W1*I(S == 1)*W2 + W1*I(S == 2)*W2","1"),
    glm_gA = ifelse(!("gA" %in% which_wrong), "W1 + W2", "1"),
    glm_QY_WAS = ifelse(!("QY_WAS" %in% which_wrong), "A + W1 + S","1"),
    glm_QY_WACY = ifelse(!("QY_WACY" %in% which_wrong), "W1*W2*CY11 + W1*W2*CY10", "1"),
    glm_QY_W = ifelse(!("QY_W" %in% which_wrong), "W1*W2", "1"),
    glm_QY_WA = ifelse(!("QY_WA" %in% which_wrong), "A + W1", "1"),
    glm_QD_WACY = ifelse(!("QD_WACY" %in% which_wrong), "W1*W2*CY11 + W1*W2*CY10", "1"),
    glm_QD_WACY_lazy = ifelse(!("QD_WACY_lazy" %in% which_wrong), "W1*W2*CY11*A + W1*W2*CY10*A", "1")
  )  
}

#' Compute summaries of simulation results for a given sample size and setting
#' 
#' Calculates the bias (times square root of sample size), the standard error 
#' (times square root of sample size), coverage of 95% confidence interval,
#' and ratio of actual standard error divided by theoretical standard error based
#' on the standard deviation of the influence function.
#' @param x A data.frame containing the formatted results. See construction of the 
#' object \code{df_final}. 
#' @param lazy Boolean indicating which estimators to summarize.

summarize <- function(x, lazy = FALSE){
  if(!lazy){
    bias_est <- mean(x$est - x$truth) * sqrt(x$n[1])
    se_est <- sd(x$est* sqrt(x$n[1]))
    ratio <- sd(x$est* sqrt(x$n[1])) / sqrt(x$var_eif[1])
    cover <- mean(x$cover)
  }else{
    bias_est <- mean(x$est_lazy - x$truth) * sqrt(x$n[1])
    se_est <- sd(x$est_lazy* sqrt(x$n[1]))
    ratio <- sd(sqrt(x$n) * x$est_lazy) / sqrt(x$var_eif[1])
    cover <- mean(x$cover_lazy)
  }
  
  return(c(bias_est, se_est, cover, ratio))
}