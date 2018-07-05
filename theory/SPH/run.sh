#!/bin/bash
cd Kernels
python kernels.py
cp kernels.pdf ..
cp kernel_derivatives.pdf ..
cd ..
cd Flavours
python plotSoundspeed.py
cp sedov_blast_soundspeed.pdf ..
cd ..
pdflatex swift_sph.tex
bibtex swift_sph.aux 
pdflatex swift_sph.tex
pdflatex swift_sph.tex
