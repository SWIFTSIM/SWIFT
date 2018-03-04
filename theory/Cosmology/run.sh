#!/bin/bash
echo "Generating PDF..."
pdflatex -jobname=cosmology cosmology_standalone.tex
bibtex cosmology.aux
pdflatex -jobname=cosmology cosmology_standalone.tex
pdflatex -jobname=cosmology cosmology_standalone.tex
