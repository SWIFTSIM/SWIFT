#!/bin/bash
python3 kernels.py
pdflatex -jobname=kernel_definitions kernel_definitions_standalone.tex
bibtex kernel_definitions.aux
pdflatex -jobname=kernel_definitions kernel_definitions_standalone.tex
pdflatex -jobname=kernel_definitions kernel_definitions_standalone.tex
