#############################################################################
#                               SIMULATION CODE                             #
#############################################################################

# specify the path
here::i_am("code/run_sim.R")
# load helper functions
source(here::here("code", "helper_fn.R"))

# on cluster
# source("/home/jran2/vaccine/JnJ/rscript/helper_fn.R")

# libraries
library(survtmle)
library(SuperLearner)

# SuperLearner library
SL_library <-  c("SL.glm", "SL.step.interaction", "SL.mean", "SL.earth")

# fitting function
fitting = function(X, cens_rate, covid_rate, set_t0, study_stop, version){
  # set seed
  set.seed(X)
  # simulate data
  data <- make_ows_data_survival(cens_rate = cens_rate, covid_rate = covid_rate, study_stop = study_stop)
  # data <- make_ows_data_survival(n = 1000, cens_rate = cens_rate, covid_rate = covid_rate, study_stop = study_stop)
  
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
      SL.ftime = SL_library,
      cvControl = list()
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
      SL.ftime = SL_library,
      glm.mediator = "age*race*risk",
      SL.trtMediator = SL_library,
      SL.eif = SL_library,
      cvControl = list()
    )
  }
  
  # compute total, in/direct effects, ci's based on fitted results
  # three rows: total, indirect, direct
  # three columns: est, cil, ciu
  out = compute_mediation_params(normal_survtmle_fit = fit1, mediation_survtmle_fit = fit2)
  
  return(out)
}

# simulation parameters
parameter_grid = expand.grid(
  # cens_rate values to be considered
  # cens_rate = seq(-10, -5, 1),
  # covid_rate value to be considered
  # covid_rate = c(-10, -8, 0.1),
  # t0
  set_t0 = 65,
  # study_stop
  # study_stop = 67,
  # seed
  X = 1:1000,
  # version
  version = c("glm", "SL")
)

# get environment variables
args = commandArgs(trailingOnly = TRUE)
iter = as.numeric(args[1])
# parameters
set_t0 = parameter_grid[iter,]$set_t0
X = parameter_grid[iter,]$X
version = parameter_grid[iter,]$version
# fit and save
result = fitting(X = X, cens_rate = -8, covid_rate = -9.3, set_t0 = set_t0, study_stop = 67, version = version)
save(result, file = paste0("/home/jran2/vaccine/JnJ/result/X", X,
                           "_t", set_t0, "_version_", version, ".RData"))

