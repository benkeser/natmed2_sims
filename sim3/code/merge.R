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
  # get truth of cumulative incidence parameters
  # ey00, ey11, ey10
  truth_est = get_ows_truth1_survival(t0 = set_t0)[c(2,1,3)]
  
  # get truth of estimated effects
  truth_eff = get_ows_truth2_survival(t0 = set_t0)
  # total, indirect, direct effect, proportion mediated
  true_effects = c(truth_eff$total, truth_eff$indirect, truth_eff$direct, 1 - log(truth_eff$direct)/log(truth_eff$total))
  
  for(version in parameter_grid_merge$version){
    # combine estimated cumulative incidence parameters together
    pt_est_matrix = matrix(nrow = 0, ncol = 3)
    cil_est_matrix = matrix(nrow = 0, ncol = 3)
    ciu_est_matrix = matrix(nrow = 0, ncol = 3)
    # combine estimated effect results together
    # total, indirect, direct effect, proportion mediated
    pt_eff_matrix = matrix(nrow = 0, ncol = 4)
    cil_eff_matrix = matrix(nrow = 0, ncol = 4)
    ciu_eff_matrix = matrix(nrow = 0, ncol = 4)
    # read in simulation results
    for(X in 1:1000){
      print(X)
      filename = paste0("X", X, "_t", set_t0, "_version_", version, ".RData")
      if(file.exists(filename)){
        load(filename)
        # store information
        # cumulative incidence parameters
        result_est = result$out_est
        pt_est_matrix = rbind(pt_est_matrix, result_est$est)
        cil_est_matrix = rbind(cil_est_matrix, result_est$cil)
        ciu_est_matrix = rbind(ciu_est_matrix, result_est$ciu)
        # estimated effects
        result_eff = result$out_eff
        pt_eff_matrix = rbind(pt_eff_matrix, result_eff$est)
        cil_eff_matrix = rbind(cil_eff_matrix, result_eff$cil)
        ciu_eff_matrix = rbind(ciu_eff_matrix, result_eff$ciu)
      }
    }
    
    # average over simulations 
    pt_est_avg = apply(pt_est_matrix, MARGIN = 2, mean)
    pt_eff_avg = apply(pt_eff_matrix, MARGIN = 2, mean)
    # use median
    # pt_est_avg = apply(pt_est_matrix, MARGIN = 2, median)
    # pt_eff_avg = apply(pt_eff_matrix, MARGIN = 2, median)
    
    # calculate bias
    bias_est = pt_est_avg - truth_est
    bias_eff = pt_eff_avg - true_effects
    
    # for each column, calculate whether the estimated effects are in [lower CI, higher CI]
    cil_est_less_than_truth = t(apply(cil_est_matrix, 1, function(x){x < truth_est}))
    cil_eff_less_than_truth = t(apply(cil_eff_matrix, 1, function(x){x < true_effects}))
    # best scenario: TRUE TRUE TRUE in each row
    ciu_est_greater_than_truth = t(apply(ciu_est_matrix, 1, function(x){x > truth_est}))
    ciu_eff_greater_than_truth = t(apply(ciu_eff_matrix, 1, function(x){x > true_effects}))
    
    # calculate the mean coverage rate for estimators
    # highest value 1(best); lowest value -1(worst)
    coverage_est = colMeans(cil_est_less_than_truth + ciu_est_greater_than_truth - 1)
    coverage_eff = colMeans(cil_eff_less_than_truth + ciu_eff_greater_than_truth - 1)
    
    # bias and coverage 
    # total, indirect, direct effect
    out = list(truth_est = truth_est,
               pt_est_avg = pt_est_avg,
               bias_est = bias_est,
               coverage_est = coverage_est,
               true_effects = true_effects,
               pt_eff_avg = pt_eff_avg,
               bias_eff = bias_eff,
               coverage_eff = coverage_eff)
    
    save(out, file = paste0("/home/jran2/vaccine/JnJ/final_result/result_t", set_t0, "_version_", version, ".RData"))
  }
}



