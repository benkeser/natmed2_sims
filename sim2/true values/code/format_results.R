#############################################################################
#                               MAKING TABLE FOR LATEX                      #
#############################################################################

library(here)
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
                `$\\bar{n}_1(sd)$` = paste0(round(mat[,6],1),"(",round(mat[,7],1),")"),
                `$\\bar{n}_2(sd)$` = paste0(round(mat[,8],1),"(",round(mat[,9],1),")"),
                check.names = FALSE)

# xtable
library(xtable)
print(xtable(df), include.rownames = F, 
      sanitize.colnames.function = identity,
      sanitize.text.function = identity)
