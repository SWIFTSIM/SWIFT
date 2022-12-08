folder=run
pwd=$(shell pwd)

${folder}/profile.png: ${folder}/sedov_0100.hdf5
	make -f make_visualisations.make -j 16 folder=${folder}

${folder}/sedov_0100.hdf5: ${folder}/swift ${folder}/sedov.hdf5 ${folder}/coolingtables/redshifts.dat ${folder}/sedov.yml
	cd ${folder}; ./swift --threads 8 --hydro --cooling --limiter sedov.yml

${folder}/sedov.hdf5: sedov.hdf5
	ln -s ${pwd}/sedov.hdf5 ${folder}/sedov.hdf5

${folder}/coolingtables/redshifts.dat: coolingtables/redshifts.dat
	ln -s ${pwd}/coolingtables ${folder}/coolingtables

${folder}/swift: ../../../swift
	mkdir -p ${folder}; ln -s ${pwd}/../../../swift ${folder}/swift

${folder}/sedov.yml: sedov.yml
	ln -s ${pwd}/sedov.yml ${folder}/sedov.yml

sedov.hdf5: glassCube_64.hdf5
	python3 makeIC.py

glassCube_64.hdf5:
	bash getGlass.sh

coolingtables/redshifts.dat:
	bash ../getEagleCoolingTable.sh
