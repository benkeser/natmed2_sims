# sim
# get environment variables
# if you want to see commands in output file
options(echo=TRUE) 
# import args
args=(commandArgs(TRUE))
print(args)
# split and assign
arguments = matrix(unlist(strsplit(args, '[=,]')),ncol=2,byrow = T)
# job_control
assign(arguments[1,1],arguments[1,2]) 
# iter
if(nrow(arguments) == 2){
  assign(arguments[2,1],as.numeric(arguments[2,2]))
}

print(job_control)
print(iter)

# simulation parameters
parameter_grid = expand.grid(
  # covid_rate value to be considered
  rates = c(-5, -4.4, -4.1, -3.8, -3.6, -3.5, -3.3, -3.2, -3.1, -3),
  # X
  size = 1000
)

# source in functions 
source("/home/jran2/vaccine/analysis/trial10_glm_inter/simulation.R")

# execute job ##################
if(job_control == "run"){
  # library
  library(natmed2)
  library(future)
  library(future.apply)
  library(SuperLearner)
  
  # do the simulation for row iter of parameter_grid
  # fit the model
  result = sim_fit(size = parameter_grid$size[iter], rate = parameter_grid$rates[iter])

  # save the result
  save(result, file = paste0("/home/jran2/vaccine/analysis/trial10_glm_inter/scratch/result_size=",parameter_grid$size[iter],
                               "_rate=",parameter_grid$rates[iter], ".RData"))
}
