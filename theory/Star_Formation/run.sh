#!/bin/bash
echo "Generating PDF..."
pdflatex -jobname=starform starformation_standalone.tex
bibtex starform.aux
pdflatex -jobname=starform starformation_standalone.tex
pdflatex -jobname=starform starformation_standalone.tex
