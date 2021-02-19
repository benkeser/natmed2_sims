#############################################################################
#                               SIMULATION CODE                             #
#############################################################################

source(here::here("code", "helper_fn.R"))
library(natmed2)
library(future)
library(future.apply)
plan(multicore)

# Generate all the scenarios for the simulation
# in all of the following scenarios, the sampling probabilities gR are correct
# because these are assumed to be known.

# all nuisance parameters are correct
all_correct <- create_setting(which_wrong = "none")
# outcome regressions are correct
gR_Qbars_correct <- create_setting(which_wrong = c("gC", "gAS", "gA", "QD_WACY",  "QD_WACY_lazy"))
# outcome regressions are correct and EIF is correct
QD_Qbars_correct <- create_setting(which_wrong = c("gC", "gAS", "gA"))
# propensity scores and mediator are correct
gR_gs_correct <- create_setting(which_wrong = c("QY_WAS", "QY_W", "QY_WA", "QD_WACY","QD_WACY_lazy"))
# propensity scores and mediator are correct and EIF is correct
QD_gs_correct <- create_setting(which_wrong = c("QY_WAS", "QY_W", "QY_WA"))
# outcome regression and propensity scores are correct
gR_QbargAC_correct <- create_setting(which_wrong = c("gAS", "QY_W", "QD_WACY","QD_WACY_lazy"))
# outcome regression and propensity scores are correct and EIF is correct
QD_QbargAC_correct <- create_setting(which_wrong = c("gAS", "QY_W"))

all_settings <- list(
  all_correct = all_correct, 
  gR_Qbars_correct = gR_Qbars_correct, 
  QD_Qbars_correct = QD_Qbars_correct, 
  gR_gs_correct = gR_gs_correct, 
  QD_gs_correct = QD_gs_correct, 
  gR_QbargAC_correct = gR_QbargAC_correct, 
  QD_QbargAC_correct = QD_QbargAC_correct
)

# define a grid of parameters for the simulation
# 1000 replicates at each sample size in each setting
parms <- expand.grid(seed = 1:1000,
                     n = c(500, 1000, 2000, 4000, 8000), 
                     settings = all_settings)

# parallelized for loop over parameter grid
rslt <- future_apply(X = parms, MARGIN = 1, FUN = one_sim, future.seed = 1)

# format the results
df <- data.frame(t(rslt))
colnames(df) <- c("est", "cil", "ciu", "est_lazy", "cil_lazy", "ciu_lazy")
df$truth <- get_truth(n = 1e7)
df$var_eif <- var_eif_ey1s0(n = 1e7)
df$cover <- df$truth < df$ciu & df$truth > df$cil
df$cover_lazy <- df$truth < df$ciu_lazy & df$truth > df$cil_lazy
df_parms <- data.frame(n = parms$n, setting = names(parms$settings))
df_final <- data.frame(df_parms, df, stringsAsFactors = FALSE)

saveRDS(df_final, file = here::here("output", "sim_result.rds"))