folder=.
snaps=$(shell ls ${folder}/sedov_*.hdf5)
imgs=$(patsubst ${folder}/sedov_%.hdf5,${folder}/SedovCooling_%.png,$(snaps))

all: ${folder}/profile.png ${folder}/SedovCooling.mp4 ${folder}/energy.png

${folder}/energy.png: ${folder}/statistics.txt
	python3 plot_energy.py ${folder}/statistics.txt ${folder}/energy.png

${folder}/profile.png: ${folder}/profile.txt
	python3 plot_time_profile.py ${folder}/profile.txt ${folder}/profile.png

${folder}/profile.txt: $(snaps)
	python3 get_time_profile.py $(snaps) ${folder}/profile.txt

${folder}/SedovCooling.mp4: $(imgs)
	ffmpeg -y -framerate 10 -pattern_type glob -i "${folder}/SedovCooling_*.png" -c:v libx264 ${folder}/SedovCooling.mp4

${folder}/SedovCooling_%.png: ${folder}/sedov_%.hdf5
	python3 plot_profile.py $< $@
