#!/bin/bash
python potential.py
pdflatex -jobname=fmm fmm_standalone.tex
bibtex fmm.aux
pdflatex -jobname=fmm fmm_standalone.tex
pdflatex -jobname=fmm fmm_standalone.tex
