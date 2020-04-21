.PHONY: report.pdf plots
plots:
	./plot-euler.py
	./mcmc1d.py
	./hmc1d.py

	
report.pdf: report.tex 
	latexmk -pdf -shell-escape report.tex && latexmk -c
