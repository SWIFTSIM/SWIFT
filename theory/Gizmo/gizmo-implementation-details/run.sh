#!/bin/bash
echo "Generating PDF..."
pdflatex -jobname=gizmo-implementation-details gizmo-implementation-details.tex
bibtex gizmo-implementation-details
pdflatex -jobname=gizmo-implementation-details gizmo-implementation-details.tex
pdflatex -jobname=gizmo-implementation-details gizmo-implementation-details.tex
