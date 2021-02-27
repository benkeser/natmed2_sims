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
  result_lazy = get(load(here::here("output",paste0("result_version=", parameter_grid$version[i],
                                                    "_rate=", parameter_grid$rates[i], "_lazy=", TRUE, ".RData"))))
  result_nonlazy = get(load(here::here("output",paste0("result_version=", parameter_grid$version[i],
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
  result_lazy = get(load(here::here("output",paste0("result_version=", parameter_grid$version[i],
                                                    "_rate=", parameter_grid$rates[i], "_lazy=", TRUE, ".RData"))))
  result_nonlazy = get(load(here::here("output",paste0("result_version=", parameter_grid$version[i],
                                                       "_rate=", parameter_grid$rates[i], "_lazy=", FALSE, ".RData"))))
  
  # calculate avg(sd)/true(sd) for lazy = TRUE
  data_lazy = get(load(here::here("output",paste0("data_version=", parameter_grid$version[i],
                                                  "_rate=", parameter_grid$rates[i], "_lazy=", TRUE, ".RData"))))
  # combine point estimator results together 
  pt_est_matrix = Reduce(rbind, result[4,])
  # combine lower CI results together 
  cil_matrix = Reduce(rbind, result[5,])
  # combine higher CI results together 
  ciu_matrix = Reduce(rbind, result[6,])
  # calculate sd
  sd_matrix = (pt_est_matrix - cil_matrix)/1.96
  median_sd = apply(sd_matrix, 2, median, na.rm = T)
  # calculate true estimate
  sds = apply(pt_est_matrix, 2, sd, na.rm = T)
  # median(sd)/true sd
  med_true_sd_lazy = median_sd/sds
  
  # calculate avg(sd)/true(sd) for lazy = FALSE
  data_nonlazy = get(load(here::here("output",paste0("data_version=", parameter_grid$version[i],
                                                     "_rate=", parameter_grid$rates[i], "_lazy=", FALSE, ".RData"))))
  # combine point estimator results together
  pt_est_matrix = Reduce(rbind, result[4,])
  # combine lower CI results together
  cil_matrix = Reduce(rbind, result[5,])
  # combine higher CI results together
  ciu_matrix = Reduce(rbind, result[6,])
  # calculate sd
  sd_matrix = (pt_est_matrix - cil_matrix)/1.96
  median_sd = apply(sd_matrix, 2, median, na.rm = T)
  # calculate true estimate
  sds = apply(pt_est_matrix, 2, sd, na.rm = T)
  # avg(sd)/true sd
  med_true_sd_nonlazy = median_sd/sds
  
  # extract the results
  mat[i,] = round(as.numeric(rbind(result_lazy$bias, result_lazy$coverage, med_true_sd_lazy, 
                                   result_nonlazy$bias, result_nonlazy$coverage, med_true_sd_nonlazy)[,3]),4)
}
# add a column: true values
true_values = sapply(unique(parameter_grid$rates), FUN = get_ows_truth1, n = 1e6)[3,]
# combine
mat = cbind(parameter_grid, true_values, mat)

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
addtorow$command = c("& &\\multicolumn{3}{c}{$\\psi_{n,1}^{+}(1,0)$} & \\multicolumn{3}{c}{$\\psi_{n,2}^{+}(1,0)$} & \\\\ \n",
                     "$\\alpha$ & True & Bias & Coverage & $\\frac{\\text{Est. std.}}{\\text{True. std.}}$ & Bias & Coverage & $\\frac{\\text{Est. std.}}{\\text{True. std.}}$ & Method \\\\ \n")

df = data.frame(mat[,c(1, 3:9)], Method = version_rows)

print(xtable(df, digits=c(0,1,4,4,4,4,4,4,4,0)), include.rownames = F, 
      hline.after = c(-1, 0, seq(n_sample_sizes, nrow(df), by = n_sample_sizes)),
      include.colnames = F,
      sanitize.colnames.function = identity,
      sanitize.text.function = identity,
      add.to.row = addtorow)

#--------------------------------------
# plot indiret effect and prop mediated
#--------------------------------------

# libraries
library(ggplot2)
library(ggthemes)
library(grid)

# simulation parameters
parameter_grid = expand.grid(
  # covid_rate value to be considered
  rates = c(-5, -4.1, -3.6, -3.3, -3.1),
  # version
  version = c("glm main", "glm inter", "SL")
)

##### plot indirect effect #####

# generate indirect dataframe for plot
gen_indirect = function(data){
  # indirect effects
  indirects = Reduce(rbind, data[3,])
  est = indirects[,1]
  # calculate bias
  indirect_avg = mean(est, na.rm = T)
  bias_indirect = indirect_avg - truth$indirect
  # confidence interval of indirect effect
  # lower bound
  cil = indirects[,2]
  # upper bound
  ciu = indirects[,3]
  # coverage
  coverage_indirect = mean((cil<truth$indirect) + (ciu>truth$indirect) - 1, na.rm = T)
  # generate plot
  df = data.frame(est = est, cil = cil, ciu = ciu)
  # df = df[!is.na(ciu),]
  index = order(df$est)
  df = df[index,]
  # log scale df
  df = log(df[,1:3])
  
  # color = #BB8FCE for a colored version, color = #575657 for a black-and-white version 
  p = ggplot(data=df, aes(x=1:nrow(df))) +
    geom_errorbar(aes(ymax = ciu, ymin = cil), color = "#575657", alpha = 1, width = 0.1, size = 0.1) +
    geom_point(aes(y=est), shape=20, size=0.1, color = "#575657") +
    geom_abline(intercept=log(truth$indirect), slope=0,color="#17202A",lwd=0.5, linetype = 2) +
    ylab("log scale range") +
    xlab("simulation") + theme_bw() + ggtitle(paste0("rate=", rate, " bias=", round(bias_indirect,4), " coverage=", round(coverage_indirect,2)))
  
  return(p)
}

# generate plots
for(version in unique(parameter_grid$version)){
  i = 0
  # create two lists for storing the pics
  p_lazy = list()
  p_nonlazy = list()
  for(rate in unique(parameter_grid$rates)){
    i = i + 1
    # get true value
    truth = get_ows_truth2(covid_rate = rate)
    # load the data
    data_lazy = get(load(here::here("output",paste0("/data_version=", version,
                                                    "_rate=", rate, "_lazy=", TRUE, ".RData"))))
    data_nonlazy = get(load(here::here("output",paste0("/data_version=", version,
                                                       "_rate=", rate, "_lazy=", FALSE, ".RData"))))
    # generate plot for indirect when lazy = TRUE
    p_lazy[[i]] = gen_indirect(data_lazy)
    # generate dataframe for plot indirect when lazy = FALSE
    p_nonlazy[[i]] = gen_indirect(data_nonlazy)
  }
  if(version == "glm main"){
    g = arrangeGrob(p_lazy[[1]]+ylim(c(-100, 100)), p_nonlazy[[1]]+ylim(c(-100, 100)), 
                    p_lazy[[2]]+ylim(c(-100, 100)), p_nonlazy[[2]]+ylim(c(-100, 100)),
                    p_lazy[[3]]+ylim(c(-100, 100)), p_nonlazy[[3]]+ylim(c(-100, 100)),
                    p_lazy[[4]]+ylim(c(-100, 100)), p_nonlazy[[4]]+ylim(c(-8, 8)),
                    p_lazy[[5]]+ylim(c(-100, 100)), p_nonlazy[[5]]+ylim(c(-8, 8)),
                    nrow = 5, ncol = 2)
  }else if(version == "glm inter"){
    g = arrangeGrob(p_lazy[[1]]+ylim(c(-100, 100)), p_nonlazy[[1]]+ylim(c(-100, 100)), 
                    p_lazy[[2]]+ylim(c(-100, 100)), p_nonlazy[[2]]+ylim(c(-100, 100)),
                    p_lazy[[3]]+ylim(c(-100, 100)), p_nonlazy[[3]]+ylim(c(-100, 100)),
                    p_lazy[[4]]+ylim(c(-50, 50)), p_nonlazy[[4]]+ylim(c(-100, 100)),
                    p_lazy[[5]]+ylim(c(-8, 8)), p_nonlazy[[5]]+ylim(c(-8, 8)),
                    nrow = 5, ncol = 2)
  }else{
    g = arrangeGrob(p_lazy[[1]]+ylim(c(-100, 100)), p_nonlazy[[1]]+ylim(c(-100, 100)), 
                    p_lazy[[2]]+ylim(c(-100, 100)), p_nonlazy[[2]]+ylim(c(-8, 8)),
                    p_lazy[[3]]+ylim(c(-100, 100)), p_nonlazy[[3]]+ylim(c(-8, 8)),
                    p_lazy[[4]]+ylim(c(-100, 100)), p_nonlazy[[4]]+ylim(c(-8, 8)),
                    p_lazy[[5]]+ylim(c(-100, 100)), p_nonlazy[[5]]+ylim(c(-8, 8)),
                    nrow = 5, ncol = 2)
  }
  ggsave(filename = here::here("output", paste0("indirect_", version,".png")), g, width = 10, height = 15)
}

## plot the inset-plot in the paper for indirect zoom-in plot(glm interaction)
# parameters
parameter_grid = expand.grid(
  # covid_rate value to be considered
  rates = c(-5, -3.1),
  # version
  version = c("glm inter"),
  # value
  value = c(TRUE, FALSE)
)

p_original = list()
p_zoom_in = list()

for(i in 1:nrow(parameter_grid)){
  # rate
  rate = parameter_grid$rates[i]
  # version
  version = parameter_grid$version[i]
  # value
  value = parameter_grid$value[i]
  
  # get true value
  truth = get_ows_truth2(covid_rate = rate)
  # load the data
  data = get(load(here::here("output",paste0("data_version=", version,
                                             "_rate=", rate, "_lazy=", value, ".RData"))))
  # orginal plot
  if(rate == -5){
    p_original[[i]] = gen_indirect(data) + ylim(c(-100, 100))
  }else{
    p_original[[i]] = gen_indirect(data) + ylim(c(-8, 8))
  }
  
  # generate plot for indirect effect
  p_zoom_in[[i]] = gen_indirect(data_lazy) + xlab("") + ylab("") + ggtitle("") + 
    theme(plot.margin = unit(c(-0.4,0,-0.4,-0.4), "cm")) + ylim(c(-8, 8))
}

# combine the orginal and zoom_in plots together
p1 = p_original[[1]] + annotation_custom(ggplotGrob(p_zoom_in[[1]]), xmin = 1, xmax = 750, 
                                         ymin = 10, ymax = 100)
p2 = p_original[[2]]
p3 = p_original[[3]] + annotation_custom(ggplotGrob(p_zoom_in[[3]]), xmin = 1, xmax = 750, 
                                         ymin = 10, ymax = 100)
p4 = p_original[[4]]

g = arrangeGrob(p1, p3,
                p2, p4,
                nrow = 2, ncol = 2)
# save the plot
ggsave(filename = here::here("output", paste0("inset_plot.png")), g, width = 10, height = 6)

#### plot prop mediated ####

# simulation parameters
parameter_grid = expand.grid(
  # covid_rate value to be considered
  rates = c(-5, -4.1, -3.6, -3.3, -3.1),
  # version
  version = c("glm main", "glm inter", "SL")
)

# generate indirect dataframe for plot
gen_prop_med = function(data){
  # calculate proportion mediated 1-log(dir)/log(total)
  totals = Reduce(rbind, data[1,])
  directs = Reduce(rbind, data[2,])
  indirects = Reduce(rbind, data[3,])
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
  for(i in 1:1000){
    sd[i] = sqrt(t(g[,i]) %*% result[7,][i][[1]] %*% g[,i])
  }
  # calculate confidence interval
  cil = prop_med - 1.96*sd
  ciu = prop_med + 1.96*sd
  
  # coverage
  coverage_propmed = mean((cil<truth_propmed) + (ciu>truth_propmed) - 1, na.rm = T)
  
  # generate plot
  df = data.frame(est = prop_med, cil = cil, ciu = ciu)
  # df = df[!is.na(ciu),]
  index = order(df$est)
  df = df[index,]
  
  # original scale
  # can not use log scale since prop med having negative values
  # color = #BB8FCE for a colored version, color = #575657 for a black-and-white version 
  p = ggplot(data=df, aes(x=1:nrow(df))) +
    geom_errorbar(aes(ymax = ciu, ymin = cil), color = "#575657", alpha = 1, width = 0.1, size = 0.1) +
    geom_point(aes(y=est), shape=20, size=0.1, color = "#575657") +
    geom_abline(intercept=truth_propmed, slope=0,color="#17202A",lwd=0.5, linetype = 2) +
    ylab("range") +
    xlab("simulation") + theme_bw() + ggtitle(paste0("rate=", rate, " bias=", round(bias_propmed,4), " coverage=", round(coverage_propmed,2)))
  
  return(p)
}

# generate plots
for(version in unique(parameter_grid$version)){
  i = 0
  # create two lists for storing the pics
  p_lazy = list()
  p_nonlazy = list()
  for(rate in unique(parameter_grid$rates)){
    i = i + 1
    # get true value
    truth = get_ows_truth2(covid_rate = rate)
    # load the data
    data_lazy = get(load(here::here("output",paste0("/data_version=", version,
                                                    "_rate=", rate, "_lazy=", TRUE, ".RData"))))
    data_nonlazy = get(load(here::here("output",paste0("/data_version=", version,
                                                       "_rate=", rate, "_lazy=", FALSE, ".RData"))))
    # generate plot for indirect when lazy = TRUE
    p_lazy[[i]] = gen_prop_med(data_lazy)
    # generate dataframe for plot indirect when lazy = FALSE
    p_nonlazy[[i]] = gen_prop_med(data_nonlazy)
  }
  if(version == "glm main"){
    g = arrangeGrob(p_lazy[[1]]+ylim(c(-20, 20)), p_nonlazy[[1]]+ylim(c(-20, 20)), 
                    p_lazy[[2]]+ylim(c(-20, 20)), p_nonlazy[[2]]+ylim(c(-20, 20)),
                    p_lazy[[3]]+ylim(c(-20, 20)), p_nonlazy[[3]]+ylim(c(-20, 20)),
                    p_lazy[[4]]+ylim(c(-2, 2)), p_nonlazy[[4]]+ylim(c(-2, 2)),
                    p_lazy[[5]]+ylim(c(-2, 2)), p_nonlazy[[5]]+ylim(c(-2, 2)),
                    nrow = 5, ncol = 2)
  }else if(version == "glm inter"){
    g = arrangeGrob(p_lazy[[1]]+ylim(c(-20, 20)), p_nonlazy[[1]]+ylim(c(-20, 20)), 
                    p_lazy[[2]]+ylim(c(-20, 20)), p_nonlazy[[2]]+ylim(c(-20, 20)),
                    p_lazy[[3]]+ylim(c(-20, 20)), p_nonlazy[[3]]+ylim(c(-20, 20)),
                    p_lazy[[4]]+ylim(c(-20, 20)), p_nonlazy[[4]]+ylim(c(-20, 20)),
                    p_lazy[[5]]+ylim(c(-2, 2)), p_nonlazy[[5]]+ylim(c(-2, 2)),
                    nrow = 5, ncol = 2)
  }else{
    g = arrangeGrob(p_lazy[[1]]+ylim(c(-20, 20)), p_nonlazy[[1]]+ylim(c(-20, 20)), 
                    p_lazy[[2]]+ylim(c(-2, 2)), p_nonlazy[[2]]+ylim(c(-2, 2)),
                    p_lazy[[3]]+ylim(c(-2, 2)), p_nonlazy[[3]]+ylim(c(-2, 2)),
                    p_lazy[[4]]+ylim(c(-2, 2)), p_nonlazy[[4]]+ylim(c(-2, 2)),
                    p_lazy[[5]]+ylim(c(-2, 2)), p_nonlazy[[5]]+ylim(c(-2, 2)),
                    nrow = 5, ncol = 2)
  }
  ggsave(filename = here::here("output", paste0("prop med_", version,".png")), g, width = 10, height = 15)
}

## plot the inset-plot in the paper for indirect zoom-in plot(glm interaction)
# parameters
parameter_grid = expand.grid(
  # covid_rate value to be considered
  rates = c(-5, -3.1),
  # version
  version = c("glm inter"),
  # value
  value = c(TRUE, FALSE)
)

p_original = list()
p_zoom_in = list()

for(i in 1:nrow(parameter_grid)){
  # rate
  rate = parameter_grid$rates[i]
  # version
  version = parameter_grid$version[i]
  # value
  value = parameter_grid$value[i]
  
  # get true value
  truth = get_ows_truth2(covid_rate = rate)
  # load the data
  data = get(load(here::here("output",paste0("data_version=", version,
                                             "_rate=", rate, "_lazy=", value, ".RData"))))
  # orginal plot
  if(rate == -5){
    p_original[[i]] = gen_prop_med(data) + ylim(c(-20, 20))
  }else{
    p_original[[i]] = gen_prop_med(data) + ylim(c(-2, 2))
  }
  
  # generate plot for indirect effect
  p_zoom_in[[i]] = gen_prop_med(data_lazy) + xlab("") + ylab("") + ggtitle("") + 
    theme(plot.margin = unit(c(-0.4,0,-0.4,-0.4), "cm")) + ylim(c(-2, 2))
}

# combine the orginal and zoom_in plots together
p1 = p_original[[1]] + annotation_custom(ggplotGrob(p_zoom_in[[1]]), xmin = 250, xmax = 1000, 
                                         ymin = 2.5, ymax = 20)
p2 = p_original[[2]]
p3 = p_original[[3]] + annotation_custom(ggplotGrob(p_zoom_in[[3]]), xmin = 250, xmax = 1000, 
                                         ymin = 2.5, ymax = 20)
p4 = p_original[[4]]

g = arrangeGrob(p1, p3,
                p2, p4,
                nrow = 2, ncol = 2)
# save the plot
ggsave(filename = here::here("output", paste0("inset_plot.png")), g, width = 10, height = 6)
