#############################################################################
#                               MAKING TABLE FOR LATEX                      #
#############################################################################
source(here::here("code", "helper_fn.R"))
library(xtable)

# load results
df_final <- readRDS(here::here("output", "sim_result.rds"))
# summarize results for both estimators across each setting by sample size
setting_names <- unique(df_final$setting)
summ_list <- vector(mode = "list", length = length(setting_names))
summ_list_lazy <- vector(mode = "list", length = length(setting_names))
ct <- 0
for(i in setting_names){
  ct <- ct + 1
  summ_list[[ct]] <- summarize_eff_rslt <- by(df_final[df_final$setting == i,], 
                         df_final$n[df_final$setting == i], 
                         summarize, lazy = FALSE)
  summ_list_lazy[[ct]] <- summarize_eff_rslt <- by(df_final[df_final$setting == i,], 
                       df_final$n[df_final$setting == i], 
                       summarize, lazy = TRUE)
}

# format into a table
summ_df <- lapply(summ_list, function(x){ 
  df <- Reduce(rbind, x)
})
summ_df <- data.frame(Reduce(rbind, summ_df))
colnames(summ_df) <- c("Bias", "SE", "Cov.", "Ratio")
summ_df_lazy <- lapply(summ_list_lazy, function(x){ 
  df <- Reduce(rbind, x)
})
summ_df_lazy <- Reduce(rbind, summ_df_lazy)
colnames(summ_df_lazy) <- c("Bias", "SE", "Cov.", "Ratio")

# column for multirow describing scenario
n_sample_sizes <- 5
add_multirow_scen <- function(scenario){
  c(paste0("\\multirow{", n_sample_sizes, "}{6em}{", scenario, "}"), rep("", n_sample_sizes - 1))
}


scenarios <- c("All", 
               "$\\bar{Q}_Y$, $\\bar{Q}_{\\bar{Q}_Y}$", 
               "$\\bar{Q}_Y$, $\\bar{Q}_{\\bar{Q}_Y}$, $\\bar{Q}_{D_{P_X}}$, $\\tilde{Q}_{D_{P_X}}$",
               "$g_{A \\mid W}$, $g_{A \\mid W, S}$, $g_C$",
               "$g_{A \\mid W}$, $g_{A \\mid W, S}$, $g_C$, $\\bar{Q}_{D_{P_X}}$, $\\tilde{Q}_{D_{P_X}}$",
               "$\\bar{Q}_Y$, $g_{A \\mid W}$, $g_C$", 
               "$\\bar{Q}_Y$, $g_{A \\mid W}$, $g_C$, $\\bar{Q}_{D_{P_X}}$, $\\tilde{Q}_{D_{P_X}}$" 
               )
scenario_rows <- c(sapply(scenarios, add_multirow_scen))

tab_df <- data.frame(Setting = scenario_rows, 
                     `$n$` = c(500, 1000, 2000, 4000, 8000), 
                     summ_df_lazy, 
                     summ_df, check.names = FALSE)

xtab_df <- xtable(tab_df, align = "c|r|rrrrr|rrrr",
                  digits = c(1, 1, 0, 2, 2, 2, 2, 2, 2, 2, 2))
print_xtab_df <- print(
  xtab_df, 
  include.rownames = FALSE, 
  hline.after = c(0, seq(n_sample_sizes, nrow(xtab_df), by = n_sample_sizes)),
  sanitize.colnames.function = identity,
  sanitize.text.function = identity,
  add.to.row = list(list(-1),  " & & \\multicolumn{4}{|c|}{$\\psi_{n,1}(1,0)$} & \\multicolumn{4}{|c|}{$\\psi_{n,2}(1,0)$}\\\\ \n \\hline "),
  file = here::here("output", "latex_table.tex")
)