#!/bin/bash
echo "Generating figures..."
python plot_potential.py
echo "Generating PDF..."
pdflatex -jobname=fmm fmm_standalone.tex
bibtex fmm.aux
pdflatex -jobname=fmm fmm_standalone.tex
pdflatex -jobname=fmm fmm_standalone.tex
