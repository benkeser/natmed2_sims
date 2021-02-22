#############################################################################
#                               SIMULATION CODE                             #
#############################################################################

# specify the path
path = "..."

source(paste0(path, "/helper_fn.R"))

# libraries
library(natmed2)
library(future.apply)

# SuperLearner library
SL_library =  c("SL.glm", "SL.step.interaction", "SL.mean", "SL.earth")

# fitting function
fitting = function(X, rate, lazy, version){
  # set seed
  set.seed(X)
  # simulate data
  dat = make_ows_data(covid_rate = rate)
  
  # calculate gRn
  gRn_fit = glm(random_subcohort ~ age*race*risk*vax, data = dat)
  gRn = predict(gRn_fit)
  gRn[dat$covid == 1] = 1
  
  # estimates using main-terms glm
  if(version == "glm main"){
    fit = natmed2(
      W = data.frame(W1 = dat$age, W2 = dat$race, W3 = dat$risk), 
      A = dat$vax, R = dat$measure_ab, S = dat$ab, C = as.numeric(!is.na(dat$covid)), Y = dat$covid,
      gRn = gRn,
      glm_gC = "W1 + W2 + W3",
      glm_gAS = "W1 + W2 + W3 + S",
      glm_QY_WAS = "A + W1 + W2 + W3 + S",
      glm_QY_WACY = "W1 + W2 + W3 + CY11 + CY10",
      glm_QY_W = "W1 + W2 + W3",
      glm_QY_WA = "A + W1 + W2 + W3",
      glm_QD_WACY = "W1 + W2 + W3 + CY11 + CY10",
      glm_QD_WACY_lazy = "A + W1 + W2 + W3 + CY11 + CY10", 
      lazy = lazy)
  }
  
  # estimates using interactions glm
  if(version == "glm inter"){
    fit = natmed2(
      W = data.frame(W1 = dat$age, W2 = dat$race, W3 = dat$risk), 
      A = dat$vax, R = dat$measure_ab, S = dat$ab, C = as.numeric(!is.na(dat$covid)), Y = dat$covid,
      gRn = gRn,
      glm_gC = "W1 + W2 + W3",
      glm_gAS = "W1*S + W2 + W3",
      glm_QY_WAS = "A + W1 + W2 + W3 + S",
      glm_QY_WACY = "W1*W2*W3*CY11 + W1*W2*W3*CY10",
      glm_QY_W = "W1*W2*W3",
      glm_QY_WA = "A*W1*W2*W3",
      glm_QD_WACY = "W1*W2*W3*CY11 + W1*W2*W3*CY10",
      glm_QD_WACY_lazy = "W1*W2*W3*CY11 + W1*W2*W3*CY10", 
      lazy = lazy)
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
  setwd(paste0(path, "/data"))
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
  setwd(paste0(path, "/data"))
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
  setwd(paste0(path, "/data"))
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



