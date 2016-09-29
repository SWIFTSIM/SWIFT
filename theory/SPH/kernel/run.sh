#!/bin/bash
python kernels.py
pdflatex kernel_definitions.tex
bibtex kernel_definitions.aux
pdflatex kernel_definitions.tex
pdflatex kernel_definitions.tex
