# merge simulation results
setwd("/home/jran2/vaccine/JnJ/result")
# load helper functions
source("/home/jran2/vaccine/JnJ/rscript/helper_fn.R")

# simulation parameters
parameter_grid_merge = expand.grid(
  # cens_rate values to be considered
  # cens_rate = seq(-10, -5, 1),
  # covid_rate value to be considered
  # covid_rate = c(-10, -8, 0.1),
  # t0
  set_t0 = 65,
  # study_stop
  # study_stop = 67,
  # seed
  # X = 1:1000,
  # version
  version = c("glm", "SL")
)

# calculate bias and coverage of estimated total, in/direct effect and proportion mediated
for(set_t0 in parameter_grid_merge$set_t0){
  # get truth
  truth = get_ows_truth2_survival(t0 = set_t0)
  # total, indirect, direct effect
  true_effects = c(truth$total, truth$indirect, truth$direct)
  for(version in parameter_grid_merge$version){
    # combine estimated effect results together
    # total, indirect, direct effect
    pt_est_matrix = matrix(nrow = 0, ncol = 3)
    cil_matrix = matrix(nrow = 0, ncol = 3)
    ciu_matrix = matrix(nrow = 0, ncol = 3)
    # read in simulation results
    for(X in 1:1000){
      filename = paste0("X", X, "_t", set_t0, "_version_", version, ".RData")
      if(file.exists(filename)){
        load(filename)
        # store information
        pt_est_matrix = rbind(pt_est_matrix, result$est)
        cil_matrix = rbind(cil_matrix, result$cil)
        ciu_matrix = rbind(ciu_matrix, result$ciu)
      }
    }

    # average over simulations 
    pt_est_avg = colMeans(pt_est_matrix)
    # calculate bias
    bias = pt_est_avg - true_effects
    
    # for each column, calculate whether the estimated effects are in [lower CI, higher CI]
    cil_less_than_truth = t(apply(cil_matrix, 1, function(x){x < true_effects}))
    # best scenario: TRUE TRUE TRUE in each row
    ciu_greater_than_truth = t(apply(ciu_matrix, 1, function(x){x > true_effects}))
    
    # calculate the mean coverage rate for estimators
    # highest value 1(best); lowest value -1(worst)
    coverage = colMeans(cil_less_than_truth + ciu_greater_than_truth - 1)
  }
  
  # bias and coverage 
  # total, indirect, direct effect
  out = list(bias = bias,
             coverage = coverage)
  
  save(out, file = paste0("/home/jran2/vaccine/JnJ/final_result/result_t", set_t0, "_version_", version, ".RData"))
}



