#############################################################################
#                               SIMULATION CODE                             #
#############################################################################

# specify the path
here::i_am("code/run_sim.R")
# load helper functions
source(here::here("code", "helper_fn.R"))

# install package if needed
if(!require("survtmle")){
  devtools::install_github("benkeser/survtmle@mediation")
}

# libraries
library(natmed2)
library(future.apply)

# SuperLearner library
SL_library <-  c("SL.glm", "SL.step.interaction", "SL.mean", "SL.earth")

# fitting function
fitting = function(X, cens_rate, covid_rate, set_t0, study_stop, version){
  # set seed
  set.seed(X)
  # simulate data
  dat <- make_ows_data_survival(cens_rate = cens_rate, covid_rate = covid_rate, study_stop = study_stop)
  
  # get sampling probabilities
  glm_fit <- glm(random_subcohort ~ age*race*risk*vax,
                 data = data)
  samp_prob <- glm_fit$fitted.values
  samp_prob[data$ftype == 1] <- 1
  
  # estimates using glm
  if(version == "glm"){
    # estimate total effect
    fit1 <- survtmle::hazard_tmle(
      ftime = data$ftime,
      ftype = data$ftype,
      trt = data$vax,
      adjustVars = data[ , c("age", "race", "risk")],
      t0 = set_t0,
      glm.ctime = "age + race + risk",
      glm.ftime = "trt + age + race + risk"
    )
    
    # estimate mediation parameter
    fit2 <- survtmle::hazard_tmle(
      ftime = data$ftime,
      ftype = data$ftype,
      trt = data$vax,
      adjustVars = data[ , c("age", "race", "risk")],
      mediator = data[ , c("ab"), drop = FALSE],
      mediatorTrtVal = 0,
      trtOfInterest = 1,
      mediatorSampProb = samp_prob,
      mediatorInCensMod = FALSE,
      t0 = set_t0,
      glm.ctime = "age + race + risk",
      glm.ftime = "trt + age + race + risk + ab",
      glm.mediator = "age*race*risk",
      glm.trtMediator = "ab + race + risk + age",
      # should be SL.eif to be more flexible?
      glm.eif = "age + race + risk + trt + ftype + ftime"
    )
  }
  
  # estimates using SuperLearner
  if(version == "SL"){
    # estimate total effect
    fit1 <- survtmle::hazard_tmle(
      ftime = data$ftime,
      ftype = data$ftype,
      trt = data$vax,
      adjustVars = data[ , c("age", "race", "risk")],
      t0 = set_t0,
      glm.ctime = "age + race + risk",
      SL.ftime = SL_library
    )
    
    # estimate mediation parameters
    fit2 <- survtmle::hazard_tmle(
      ftime = data$ftime,
      ftype = data$ftype,
      trt = data$vax,
      adjustVars = data[ , c("age", "race", "risk")],
      mediator = data[ , c("ab"), drop = FALSE],
      mediatorTrtVal = 0,
      trtOfInterest = 1,
      mediatorSampProb = samp_prob,
      mediatorInCensMod = FALSE,
      t0 = set_t0,
      glm.ctime = "age + race + risk",
      glm.ftime = "trt + age + race + risk + ab",
      glm.mediator = "age*race*risk",
      glm.trtMediator = "ab + race + risk + age",
      # should be SL.eif to be more flexible?
      glm.eif = "age + race + risk + trt + ftype + ftime"
    )
  }
  
  if(version == "SL"){
    fit = natmed2(
      W = data.frame(W1 = dat$age, W2 = dat$race, W3 = dat$risk), 
      A = dat$vax, R = dat$measure_ab, S = dat$ab, C = as.numeric(!is.na(dat$covid)), Y = dat$covid,
      gRn = gRn,
      glm_gA = ".",  
      glm_gC = "1", 
      SL_gC = NULL,
      glm_gAS = NULL, 
      SL_gAS = SL_library,
      glm_QY_WAS = NULL, 
      SL_QY_WAS = SL_library, 
      glm_QY_WACY = NULL,
      SL_QY_WACY = SL_library,
      glm_QY_W = NULL,
      SL_QY_W = SL_library, 
      glm_QY_WA = NULL,
      SL_QY_WA = SL_library,
      glm_QD_WACY = NULL,
      SL_QD_WACY = SL_library,
      glm_QD_WACY_lazy = NULL,
      SL_QD_WACY_lazy = SL_library,
      lazy = lazy)
  }
  
  if(lazy == FALSE){
    return(list(total = fit$eff[1 ,2:4], direct = fit$eff[2 ,2:4], 
                indirect = fit$eff[3 ,2:4], estimators = fit$risk[,1], 
                cil = fit$risk[,2], ciu = fit$risk[,3], cov = fit$cov))
  }else{
    return(list(total = fit$eff_lazy[1 ,2:4], direct = fit$eff_lazy[2 ,2:4], 
                indirect = fit$eff_lazy[3 ,2:4], estimators = fit$risk_lazy[,1], 
                cil = fit$risk_lazy[,2], ciu = fit$risk_lazy[,3] , cov = fit$cov_lazy))
  }
}

# bias, confidence interval
sim_fit = function(size, rate, lazy, version){
  
  # plan(multisession)
  # use parallel computing to get the result
  result = future_sapply(X = 1:size, FUN = fitting, rate = rate, lazy = lazy, version = version)
  # save the result
  setwd(paste0(path, "/output"))
  save(result, file = paste0("data_version=", version, "_rate=", rate,"_lazy=", lazy, ".RData"))
  
  #### part I: risk estimator and coverage ####
  
  # get truth
  truth = get_ows_truth1(covid_rate = rate)
  
  # combine point estimator results together [size(1000) simulations*(4+4)]
  pt_est_matrix = Reduce(cbind, result[4,])
  # average the 1000 simulations and get 4 point estimators + 4 point estimators
  pt_est_avg = rowMeans(pt_est_matrix)
  # calculate bias
  bias = pt_est_avg - truth
  
  # combine lower CI results together [size(1000) simulations*(4+4)]
  cil_matrix = Reduce(rbind, result[5,])
  # combine higher CI results together [size(1000) simulations*(4+4)]
  ciu_matrix = Reduce(rbind, result[6,])
  
  # for each row, calculate whether the point estimator is in [lower CI, higher CI]
  cil_less_than_truth = apply(cil_matrix, 1, function(x){x < truth})
  # best scenario: TRUE TRUE TRUE TRUE
  ciu_greater_than_truth = apply(ciu_matrix, 1, function(x){x > truth})
  # best scenario: TRUE TRUE TRUE TRUE
  
  # calculate the mean coverage rate for 4 point estimators
  # highest value 1(best); lowest value -1(worst)
  coverage = rowMeans(cil_less_than_truth + ciu_greater_than_truth - 1)
  
  #### part II: estimate indirect effect and prop. mediated and their coverage ####
  
  # get truth
  truth = get_ows_truth2(covid_rate = rate)
  
  # bias of indirect effect
  indirects = Reduce(rbind, result[3,])
  # save the estimates
  setwd(paste0(path, "/output"))
  save(indirects, file = paste0("indirects_version=", version, "_rate=", rate,"_lazy=", lazy, ".RData"))
  
  # calculate bias
  indirect_avg = mean(indirects[,1], na.rm = T)
  bias_indirect = indirect_avg - truth$indirect
  
  # confidence interval of indirect effect
  # lower bound
  cil = indirects[,2]
  # upper bound
  ciu = indirects[,3]
  # coverage
  coverage_indirect = mean((cil<truth$indirect) + (ciu>truth$indirect) - 1, na.rm = T)
  
  # bias of proportion mediated 1-log(dir)/log(total)
  totals = Reduce(rbind, result[1,])
  directs = Reduce(rbind, result[2,])
  # save the estimates
  setwd(paste0(path, "/output"))
  save(directs, file = paste0("directs_version=", version, "_rate=", rate,"_lazy=", lazy, ".RData"))
  save(totals, file = paste0("totals_version=", version, "_rate=", rate,"_lazy=", lazy, ".RData"))
  
  prop_med = 1 - log(directs[,1])/log(totals[,1])
  prop_med_avg = mean(prop_med, na.rm = T)
  
  # bias
  truth_propmed = 1 - log(truth$direct)/log(truth$total)
  bias_propmed = prop_med_avg - truth_propmed
  
  # confidence interval of proportion mediated
  estimators = Reduce(rbind, result[4,])
  g = matrix(c(log(directs[,1])/(log(totals[,1]))^2*1/estimators[,1], log(indirects[,1])/(log(totals[,1]))^2*1/estimators[,2],
               -1/(estimators[,3]*log(totals[,1])), rep(0, 1000)), nrow = 4, byrow = T)
  sd = c()
  for(i in 1:size){
    sd[i] = sqrt(t(g[,i]) %*% result[7,][i][[1]] %*% g[,i])
  }
  
  # calculate confidence interval
  cil = prop_med - 1.96*sd
  ciu = prop_med + 1.96*sd
  
  # coverage
  coverage_propmed = mean((cil<truth_propmed) + (ciu>truth_propmed) - 1, na.rm = T)
  
  return(list(bias = bias, coverage = coverage, comb = c(bias_indirect, coverage_indirect, bias_propmed, coverage_propmed)))
}

# simulation
# get environment variables
# if you want to see commands in output file
options(echo=TRUE) 
# import args
args=(commandArgs(TRUE))
# split and assign
arguments = matrix(unlist(strsplit(args, '[=,]')),ncol=2,byrow = T)
# iter
if(nrow(arguments) == 1){
  assign(arguments[1,1],as.numeric(arguments[1,2]))
}

# simulation parameters
parameter_grid = expand.grid(
  # covid_rate value to be considered
  rates = c(-5, -4.1, -3.6, -3.3, -3.1),
  # X
  size = 1000,
  # lazies
  lazy = c(TRUE, FALSE), 
  # version
  version = c("glm main", "glm inter", "SL")
)

# execute job ##################

# library
library(natmed2)
library(future)
library(future.apply)
library(SuperLearner)

# do the simulation for row iter of parameter_grid
result = sim_fit(size = parameter_grid$size[iter], rate = parameter_grid$rates[iter], lazy = parameter_grid$lazy[iter], version = parameter_grid$version[iter])

# save the result
save(result, file = paste0(path, "/scratch/result_version=", parameter_grid$version[iter],
                           "_rate=", parameter_grid$rates[iter], "_lazy=", parameter_grid$lazy[iter], ".RData"))