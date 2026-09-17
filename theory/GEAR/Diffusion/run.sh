#!/bin/bash
# Build the combined parabolic + hyperbolic diffusion theory document,
# running pdflatex/bibtex enough times to resolve cross-references and
# citations. Run with "clean" to remove generated build artifacts
# (aux/log/bbl/blg/toc/pdf and the auto/ directory) instead of building.
set -e

jobname=GEAR_FVPM_Diffusion_Theory
texfile=GEAR_FVPM_Diffusion_Theory.tex

clean() {
  rm -f "$jobname".aux "$jobname".log "$jobname".bbl "$jobname".blg \
        "$jobname".toc "$jobname".pdf
  rm -rf auto
}

build() {
  pdflatex -interaction=nonstopmode -jobname="$jobname" "$texfile"
  bibtex "$jobname"
  pdflatex -interaction=nonstopmode -jobname="$jobname" "$texfile"
  pdflatex -interaction=nonstopmode -jobname="$jobname" "$texfile"
}

if [ "$1" = "clean" ]; then
  clean
else
  build
fi
