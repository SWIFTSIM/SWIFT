#!/bin/bash

pdflatex pasc_paper.tex
bibtex pasc_paper.aux
pdflatex pasc_paper.tex
pdflatex pasc_paper.tex
