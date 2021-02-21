#############################################################################
#                               SIMULATION CODE                             #
#############################################################################

source(here::here("code", "helper_fn.R"))

# libraries
library(future.apply)

# number of cases for vaccine and placebo
covid_vax_control = function(X, rate){
  # set seed
  set.seed(X)
  # simulate data
  dat = make_ows_data(covid_rate = rate)
  
  # calculation
  return(c(sum(dat$vax == 1 & dat$covid == 1), sum(dat$vax == 0 & dat$covid == 1)))
}

# calculate true 1-VE, indirect, direct effect, prop. mediated, 
# and average number of cases for vaccine and placebo over 1000 simulations
# for specific covid_rate values
mat = matrix(0, nrow = 0, ncol = 9)
for(rate in c(-5, -4.1, -3.6, -3.3, -3.1)){
  # get the true estimates
  result1 = get_ows_truth1(covid_rate = rate)
  result2 = get_ows_truth2(covid_rate = rate)
  
  # 1-VE
  one_minus_VE = 1-result1[1]/result1[2]
  
  # indirect
  indirect = result2$indirect
  
  # direct
  direct = result2$direct
  
  # prop med
  total = result2$total
  prop_med = 1 - log(direct)/log(total)
  
  # average number of cases for vaccine and placebo over 1000 simulations
  num_vax_control = future_sapply(X = 1:1000, FUN = covid_vax_control, rate = rate, future.seed = NULL)
  # n1 in vax
  avg_n1 = mean(num_vax_control[1,])
  sd_n1 = sd(num_vax_control[1,])
  # n0 in placebo
  avg_n0 = mean(num_vax_control[2,])
  sd_n0 = sd(num_vax_control[2,])
  
  mat = rbind(mat, c(rate, one_minus_VE, indirect, direct, prop_med, avg_n1, sd_n1, avg_n0, sd_n0))
}

save(mat, file = here::here("output", "table2.RData"))
