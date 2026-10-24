# Minimal makefile to build onboarding guide

default: pdf

# Put it first so that "make" without argument is like "make help".
help:
	echo "run `make pdf` to build the onboarding guide."

clean:
	rm -f *.aux *.bbl *.log *.out *.gz *.toc *.blg *.lot *.lof *.run.xml *.snm *.nav *-blx.bib onboarding.pdf

pdf: source/onboarding.tex
	pdflatex source/onboarding.tex
	pdflatex source/onboarding.tex

.PHONY: help Makefile clean

