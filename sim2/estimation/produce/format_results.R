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
addtorow$command = c("& \\multicolumn{4}{c}{$\\psi_{n,1}^{+}$} & \\multicolumn{4}{c}{$\\psi_{n,2}^{+}$} & \\\\ \n",
                     "& \\multicolumn{2}{c}{Indirect Effect} & \\multicolumn{2}{c}{Prop. Mediated} & \\multicolumn{2}{c}{Indirect Effect} & \\multicolumn{2}{c}{Prop. Mediated} & \\\\ \n",
                     "$\\alpha$ & Bias & Coverage & Bias & Coverage & Bias & Coverage & Bias & Coverage & Method \\\\ \n")

df = data.frame(mat[,c(1, 3:10)], Method = version_rows)

print(xtable(df, digits=c(0,1,3,3,3,3,3,3,3,3,0)), include.rownames = F, 
      hline.after = c(-1, 0, seq(n_sample_sizes, nrow(df), by = n_sample_sizes)),
      include.colnames = F,
      sanitize.colnames.function = identity,
      sanitize.text.function = identity,
      add.to.row = addtorow)

#---------------------------------------------------------
# calculate bias/coverage of the risk estimator E(Y(1), 0)
#---------------------------------------------------------

# simulation parameters
parameter_grid = expand.grid(
  # covid_rate value to be considered
  rates = c(-5, -4.1, -3.6, -3.3, -3.1),
  # version
  version = c("glm main", "glm inter", "SL")
)

# combine the results to generate table 3
mat = matrix(nrow = nrow(parameter_grid), ncol = 6)
for(i in 1:nrow(parameter_grid)){
  # get the files
  result_lazy = get(load(here::here("output",paste0("/result_version=", parameter_grid$version[i],
                                                    "_rate=", parameter_grid$rates[i], "_lazy=", TRUE, ".RData"))))
  result_nonlazy = get(load(here::here("output",paste0("/result_version=", parameter_grid$version[i],
                                                       "_rate=", parameter_grid$rates[i], "_lazy=", FALSE, ".RData"))))
  
  # calculate avg(sd)/true(sd) for lazy = TRUE
  data_lazy = get(load(here::here("output",paste0("/data_version=", parameter_grid$version[i],
                                                  "_rate=", parameter_grid$rates[i], "_lazy=", TRUE, ".RData"))))
  # combine point estimator results together 
  pt_est_matrix = Reduce(rbind, result[4,])
  # combine lower CI results together 
  cil_matrix = Reduce(rbind, result[5,])
  # combine higher CI results together 
  ciu_matrix = Reduce(rbind, result[6,])
  # calculate sd
  sd_matrix = (pt_est_matrix - cil_matrix)/1.96
  avg_sd = apply(sd_matrix, 2, mean, na.rm = T)
  # calculate true estimate
  sds = apply(pt_est_matrix, 2, sd, na.rm = T)
  # avg(sd)/true sd
  avg_true_sd_lazy = avg_sd/sds
  
  # calculate avg(sd)/true(sd) for lazy = FALSE
  data_nonlazy = get(load(here::here("output",paste0("/data_version=", parameter_grid$version[i],
                                                     "_rate=", parameter_grid$rates[i], "_lazy=", FALSE, ".RData"))))
  # combine point estimator results together
  pt_est_matrix = Reduce(rbind, result[4,])
  # combine lower CI results together
  cil_matrix = Reduce(rbind, result[5,])
  # combine higher CI results together
  ciu_matrix = Reduce(rbind, result[6,])
  # calculate sd
  sd_matrix = (pt_est_matrix - cil_matrix)/1.96
  avg_sd = apply(sd_matrix, 2, mean, na.rm = T)
  # calculate true estimate
  sds = apply(pt_est_matrix, 2, sd, na.rm = T)
  # avg(sd)/true sd
  avg_true_sd_nonlazy = avg_sd/sds
  
  # extract the results
  mat[i,] = round(as.numeric(rbind(result_lazy$bias, result_lazy$coverage, avg_true_sd_lazy, 
                                   result_nonlazy$bias, result_nonlazy$coverage, avg_true_sd_nonlazy)[,3]),4)
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
addtorow$pos = list(-1, 0)
addtorow$command = c("& \\multicolumn{3}{c}{$\\psi_{n,1}^{+}(1,0)$} & \\multicolumn{3}{c}{$\\psi_{n,2}^{+}(1,0)$} & \\\\ \n",
                     "$\\alpha$ & Bias & Coverage & $\\frac{\\text{Est. std.}}{\\text{True. std.}}$ & Bias & Coverage & $\\frac{\\text{Est. std.}}{\\text{True. std.}}$ & Method \\\\ \n")

df = data.frame(mat[,c(1, 3:8)], Method = version_rows)

print(xtable(df, digits=c(0,1,4,4,4,4,4,4,0)), include.rownames = F, 
      hline.after = c(-1, 0, seq(n_sample_sizes, nrow(df), by = n_sample_sizes)),
      include.colnames = F,
      sanitize.colnames.function = identity,
      sanitize.text.function = identity,
      add.to.row = addtorow)
