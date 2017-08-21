#!/bin/bash
if [ ! -e potential.pdf ]
then
    echo "Generating 1st figure..."
    python plot_potential.py
fi
if [ ! -e potential_short.pdf ]
then
    echo "Generating 2nd figures..."
    python plot_mesh.py
fi
echo "Generating PDF..."
pdflatex -jobname=fmm fmm_standalone.tex
bibtex fmm.aux
pdflatex -jobname=fmm fmm_standalone.tex
pdflatex -jobname=fmm fmm_standalone.tex
