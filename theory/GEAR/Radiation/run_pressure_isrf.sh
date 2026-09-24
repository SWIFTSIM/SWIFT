#!/bin/bash
echo "Generating PDF..."
PDFLATEX="pdflatex -interaction=nonstopmode -halt-on-error -jobname=radiation_pressure_isrf"
$PDFLATEX radiation_pressure_isrf.tex || { echo "pdflatex failed, see radiation_pressure_isrf.log"; exit 1; }
bibtex radiation_pressure_isrf.aux
$PDFLATEX radiation_pressure_isrf.tex || { echo "pdflatex failed, see radiation_pressure_isrf.log"; exit 1; }
$PDFLATEX radiation_pressure_isrf.tex || { echo "pdflatex failed, see radiation_pressure_isrf.log"; exit 1; }
$PDFLATEX radiation_pressure_isrf.tex || { echo "pdflatex failed, see radiation_pressure_isrf.log"; exit 1; }
