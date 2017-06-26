#!/bin/bash
if [! -e potential.pdf ]
then
    echo "Generating figures..."
    python plot_potential.py
fi
echo "Generating PDF..."
pdflatex -jobname=fmm fmm_standalone.tex
bibtex fmm.aux
pdflatex -jobname=fmm fmm_standalone.tex
pdflatex -jobname=fmm fmm_standalone.tex
