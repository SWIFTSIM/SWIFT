main_file = swift

default: all

all: ${main_file}.pdf 

Figures/Gravity/correction.png Figures/Gravity/error.png Figures/Gravity/force.png: Figures/Gravity/force.py
	python3 Figures/Gravity/force.py
	mv correction.png Figures/Gravity
	mv error.png Figures/Gravity
	mv force.png Figures/Gravity

${main_file}.pdf: ${main_file}.tex Figures/Gravity/correction.png Figures/Gravity/error.png Figures/Gravity/force.png
	# pdflatex ${main_file}.tex
	# bibtex ${main_file}
	pdflatex ${main_file}.tex
	pdflatex ${main_file}.tex
clean:
	rm -f ${main_file}.pdf *.aux *.log *.out *.dvi *.bbl *.blg *.toc *.gz *.lof *.lot Figures/Gravity/correction.png Figures/Gravity/error.png Figures/Gravity/force.png

.PHONY: clean

