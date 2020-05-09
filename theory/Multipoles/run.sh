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
if [ ! -e alpha_powers.pdf ]
then
    echo "Generating derivative figures..."
    python plot_derivatives.py
fi
if [ ! -e mac_potential.pdf ]
then
    echo "Generating derivative figures..."
    python3 plot_mac_potential.py
fi
echo "Generating PDF..."
pdflatex -jobname=fmm fmm_standalone.tex
bibtex fmm.aux
pdflatex -jobname=fmm fmm_standalone.tex
pdflatex -jobname=fmm fmm_standalone.tex
