#!/bin/bash
pdflatex uldm.tex
bibtex uldm.aux 
pdflatex uldm.tex
pdflatex uldm.tex
