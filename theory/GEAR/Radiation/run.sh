#!/bin/bash
# --figures: regenerate the schematic PDFs in figures/ (needs matplotlib)
# before building the document. Default build uses the committed PDFs.
if [ "$1" = "--figures" ]; then
  echo "Regenerating figures..."
  (cd figures && for f in fig*.py; do python3 "$f" || exit 1; done) || exit 1
fi
echo "Generating PDF..."
pdflatex -jobname=radiation radiation.tex
bibtex radiation.aux
pdflatex -jobname=radiation radiation.tex
pdflatex -jobname=radiation radiation.tex
