#!/bin/bash
echo "Generating PDF..."
pdflatex -jobname=radiation_pressure_isrf radiation_pressure_isrf.tex
bibtex radiation_pressure_isrf.aux
pdflatex -jobname=radiation_pressure_isrf radiation_pressure_isrf.tex
pdflatex -jobname=radiation_pressure_isrf radiation_pressure_isrf.tex
