#!/bin/bash

# --figures: regenerate the schematic PDFs in figures/ (needs matplotlib)
# before building the document. Default build uses the committed PDFs.
if [ "$1" = "--figures" ]; then
  echo "Regenerating figures..."
  (cd figures && for f in fig*.py; do python3 "$f" || exit 1; done) || exit 1
fi

# -----------------------------------------------------------------------------
echo "Generating PDF for the HII photoionisation algorithm..."
PDFLATEX="pdflatex -interaction=nonstopmode -halt-on-error -jobname=radiation"
$PDFLATEX radiation.tex || { echo "pdflatex failed, see radiation.log"; exit 1; }
bibtex radiation.aux
$PDFLATEX radiation.tex || { echo "pdflatex failed, see radiation.log"; exit 1; }
$PDFLATEX radiation.tex || { echo "pdflatex failed, see radiation.log"; exit 1; }

# -----------------------------------------------------------------------------
echo "Generating PDF for the radiation pressure and the ISRF subgri models..."
PDFLATEX="pdflatex -interaction=nonstopmode -halt-on-error -jobname=radiation_pressure_isrf"
$PDFLATEX radiation_pressure_isrf.tex || { echo "pdflatex failed, see radiation_pressure_isrf.log"; exit 1; }
bibtex radiation_pressure_isrf.aux
$PDFLATEX radiation_pressure_isrf.tex || { echo "pdflatex failed, see radiation_pressure_isrf.log"; exit 1; }
$PDFLATEX radiation_pressure_isrf.tex || { echo "pdflatex failed, see radiation_pressure_isrf.log"; exit 1; }
$PDFLATEX radiation_pressure_isrf.tex || { echo "pdflatex failed, see radiation_pressure_isrf.log"; exit 1; }
