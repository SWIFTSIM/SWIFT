#!/bin/bash
# Build the combined parabolic + hyperbolic diffusion theory document,
# running pdflatex/bibtex enough times to resolve cross-references and
# citations.
set -e

build() {
  local jobname=$1
  local texfile=$2
  pdflatex -interaction=nonstopmode -jobname="$jobname" "$texfile"
  bibtex "$jobname"
  pdflatex -interaction=nonstopmode -jobname="$jobname" "$texfile"
  pdflatex -interaction=nonstopmode -jobname="$jobname" "$texfile"
}

build GEAR_FVPM_Diffusion_Theory GEAR_FVPM_Diffusion_Theory.tex
