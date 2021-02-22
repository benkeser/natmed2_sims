#############################################################################
#                     MAKING TABLE AND PLOTS FOR LATEX                      #
#############################################################################

# library
library(here)

# set working directory
setwd(here::here("output"))

#---------------------------------------------------------------
# calculate bias&coverage of indirect effects and prop. mediated
#---------------------------------------------------------------

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









#----------------------------------------------
# calculate bias/coverage of the risk estimates
#----------------------------------------------

# simulation parameters
parameter_grid = expand.grid(
  # covid_rate value to be considered
  rates = c(-5, -4.4, -4.1, -3.8, -3.6, -3.5, -3.3, -3.2, -3.1, -3)
)

# get the files
result.glmmain = matrix(nrow = 10, ncol = 8)
for(i in 1:nrow(parameter_grid)){
  tmp = get(load(paste0("/Users/scarlett/Box\ Sync/vaccine/project/simulation/simulations/scenario\ 12\ analysis/glm\ inter/risk\ result/result_size=", 1000,
                        "_rate=",parameter_grid$rates[i],".RData")))
  result.glmmain[i,] = tmp$coverage
}
result.glmmain = cbind(parameter_grid, result.glmmain)
result.glmmain = result.glmmain[,-c(6,7)]
colnames(result.glmmain) = c("covid_rate", "ey11", "ey00", "ey10", "ey01", "ey10 lazy", "ey01 lazy")

# print
library(xtable)
print(xtable(result.glmmain, digits=c(0, 1,6,6,6,6,6,6)),include.rownames = F)
print(xtable(result.glmmain, digits=c(0, 1,4,4,4,4,4,4)),include.rownames = F)