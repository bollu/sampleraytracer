report.pdf: report.tex
	./plot-euler.py
	./mcmc1d.py
	./hmc1d.py
	latexmk -pdf -shell-escape report.tex && latexmk -c
