# First plot the IMF image
python3 ./plot_imf_image.py

#Then, run latex to produce the output
pdflatex -jobname=gear_imf_sampling gear_imf_sampling.tex

