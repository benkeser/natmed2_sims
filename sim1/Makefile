## all            : sim_results, latex_table, and pdf_table
all: sim_results latex_table pdf_table

## latex_table    : format sim results into latex table saved in output/latex_table.tex
latex_table: sim_results format_code
	Rscript code/format_results.R
format_code: code/format_results.R code/helper_fn.R

## sim_results    : run simulation and produces output/sim_result.rds
sim_results: sim_code
	Rscript code/run_sim.R
sim_code: code/run_sim.R code/helper_fn.R

## pdf_table      : compile standalone pdf with simulation results
pdf_table: latex_table code/table.tex
	pdflatex -halt-on-error -output-directory output code/table.tex && rm output/table.aux output/table.log

## clean          : delete contents of output directory
clean:
	rm -f output/*

help: Makefile
	@sed -n 's/^##//p' $<

.PHONY: format_code sim_code latex_table results pdf_table help