#!/bin/bash
pdflatex sph_flavours.tex
bibtex sph_flavours.aux 
pdflatex sph_flavours.tex
pdflatex sph_flavours.tex
