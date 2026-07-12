#!/bin/bash
echo "Generating PDF..."

pdflatex -jobname=GEARRT GEARRT.tex
bibtex GEARRT
pdflatex -jobname=GEARRT GEARRT.tex
pdflatex -jobname=GEARRT GEARRT.tex
