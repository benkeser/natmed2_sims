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
  # version
  version = c("glm main", "glm inter", "SL")
)

# combine the results to generate table
mat = matrix(nrow = nrow(parameter_grid), ncol = 8)
for(i in 1:nrow(parameter_grid)){
  # get the files
  result_lazy = get(load(here::here("output",paste0("/result_version=", parameter_grid$version[i],
                                                    "_rate=", parameter_grid$rates[i], "_lazy=", TRUE, ".RData"))))
  result_nonlazy = get(load(here::here("output",paste0("/result_version=", parameter_grid$version[i],
                                                       "_rate=", parameter_grid$rates[i], "_lazy=", FALSE, ".RData"))))
  # extract the results
  mat[i,] = round(c(result_lazy$comb, result_nonlazy$comb),2)
}
mat = cbind(parameter_grid, mat)

# print the table
# column for multirow describing scenario
n_sample_sizes = 5
add_multirow_scen <- function(scenario){
  c(paste0("\\multirow{", n_sample_sizes, "}{6em}{", scenario, "}"), rep("", n_sample_sizes - 1))
}

versions = c("GLM main", "GLM inter", "SuperLearner")
version_rows = c(sapply(versions, add_multirow_scen))

# row names
addtorow = list()
addtorow$pos = list(-1, 0, 0)
addtorow$command = c("& \\multicolumn{4}{c}{$\\psi_{n,1}(1,0)$} & \\multicolumn{4}{c}{$\\psi_{n,2}(1,0)$} & \\\\ \n",
                     "& \\multicolumn{2}{c}{Indirect Effect} & \\multicolumn{2}{c}{Prop. Mediated} & \\multicolumn{2}{c}{Indirect Effect} & \\multicolumn{2}{c}{Prop. Mediated} & \\\\ \n",
                     "$\\alpha$ & Bias & Coverage & Bias & Coverage & Bias & Coverage & Bias & Coverage & Method \\\\ \n")

df = data.frame(mat[,c(1, 3:10)], Method = version_rows)

print(xtable(df), include.rownames = F, 
      include.colnames = F,
      sanitize.colnames.function = identity,
      sanitize.text.function = identity,
      add.to.row = addtorow)






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