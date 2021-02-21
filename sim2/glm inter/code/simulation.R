# library
library(natmed2)
library(future.apply)
# load R script
source("/home/jran2/vaccine/analysis/trial10_glm_inter/covid_med_sim.R")

# fitting function
fitting = function(X, rate){
  # set seed
  set.seed(X)
  # simulate data
  dat = make_ows_data(covid_rate = rate)
  
  # calculate gRn
  gRn_fit = glm(random_subcohort ~ age*race*risk*vax, data = dat)
  gRn = predict(gRn_fit)
  gRn[dat$covid == 1] = 1
  
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
    lazy = TRUE)
  
  return(rbind(fit$risk[,1:3], fit$risk_lazy[,1:3]))
}

# bias, confidence interval
sim_fit = function(size, rate){
  
  plan(multisession)
  # use parallel computing to get the result
  result = future_sapply(X = 1:size, FUN = fitting, rate = rate)
  # # save the result
  # save(result, file = paste0("/home/jran2/vaccine/analysis/trial10_glm_inter/scratch/initial_result_size=", parameter_grid$size[iter],
  #                            "_rate=", parameter_grid$rates[iter], ".RData"))

  # get truth
  truth = get_ows_truth(covid_rate = rate)
  
  # combine point estimator results together [size(1000) simulations*(4+4)]
  pt_est_matrix = Reduce(rbind, result[1,])
  # average the 1000 simulations and get 4 point estimators + 4 point estimators
  pt_est_avg = colMeans(pt_est_matrix)
  # calculate bias
  bias = pt_est_avg - truth
  
  # combine lower CI results together [size(1000) simulations*(4+4)]
  cil_matrix = Reduce(rbind, result[2,])
  # combine higher CI results together [size(1000) simulations*(4+4)]
  ciu_matrix = Reduce(rbind, result[3,])
  
  # for each row, calculate whether the point estimator is in [lower CI, higher CI]
  cil_less_than_truth = apply(cil_matrix, 1, function(x){x < truth})
  # best scenario: TRUE TRUE TRUE TRUE
  ciu_greater_than_truth = apply(ciu_matrix, 1, function(x){x > truth})
  # best scenario: TRUE TRUE TRUE TRUE
  
  # calculate the mean coverage rate for 4 point estimators
  # highest value 1(best); lowest value -1(worst)
  coverage = rowMeans(cil_less_than_truth + ciu_greater_than_truth - 1)
  
  return(list(bias = bias, coverage = coverage))
}

