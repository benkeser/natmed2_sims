#############################################################################
#                               MAKING TABLE FOR LATEX                      #
#############################################################################


library(xtable)

# load results
load(here::here("output", "table2.RData"))

# create a dataframe
mat = round(mat, 2)
df = data.frame(`$\\alpha$` = mat[,1],
                "1-VE" = mat[,2],
                `$\\frac{\\psi(1,1)}{\\psi(1,0)}$` = mat[,3],
                `$\\frac{\\psi(1,0)}{\\psi(0,0)}$` = mat[,4],
                "Prop. Mediated" = mat[,5],
                `$\\bar{n_1}(sd)$` = paste0(mat[,6],"(",mat[,7],")"),
                `$\\bar{n_2}(sd)$` = paste0(mat[,8],"(",mat[,9],")"),
                check.names = FALSE)

# generate the table
library(xtable)
print(xtable(df), include.rownames = F, 
      sanitize.colnames.function = identity,
      sanitize.text.function = identity)
