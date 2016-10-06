#!/bin/bash
python kernels.py
pdflatex -jobname=eos eos_standalone.tex
bibtex eos.aux
pdflatex -jobname=eos eos_standalone.tex
pdflatex -jobname=eos eos_standalone.tex
