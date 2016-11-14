#!/bin/bash
python kernels.py
pdflatex -jobname=kernel_definitions kernel_definitions_standalone.tex
bibtex kernel_definitions.aux
pdflatex -jobname=kernel_definitions kernel_definitions_standalone.tex
pdflatex -jobname=kernel_definitions kernel_definitions_standalone.tex
