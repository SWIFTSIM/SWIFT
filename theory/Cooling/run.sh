#!/bin/bash
echo "Generating PDF..."
pdflatex -jobname=eagle_cooling eagle_cooling.tex
bibtex eagle_cooling.aux
pdflatex -jobname=eagle_cooling eagle_cooling.tex
pdflatex -jobname=eagle_cooling eagle_cooling.tex
