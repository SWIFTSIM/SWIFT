#!/bin/bash
python3 plotSoundspeed.py

pdflatex -jobname=sph_flavours sph_flavours_standalone.tex
bibtex sph_flavours.aux 
pdflatex -jobname=sph_flavours sph_flavours_standalone.tex
pdflatex -jobname=sph_flavours sph_flavours_standalone.tex
