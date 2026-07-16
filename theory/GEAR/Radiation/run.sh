#!/bin/bash
echo "Generating PDF..."
pdflatex -jobname=radiation radiation.tex
bibtex radiation.aux
pdflatex -jobname=radiation radiation.tex
pdflatex -jobname=radiation radiation.tex
