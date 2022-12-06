
profile.png: sedov_0100.hdf5
	make -f make_visualisations.make -j 16

sedov_0100.hdf5: sedov.hdf5 coolingtables/redshifts.dat
	../../../swift --threads 8 --hydro --cooling --limiter sedov.yml

sedov.hdf5: glassCube_64.hdf5
	python3 makeIC.py

glassCube_64.hdf5:
	bash getGlass.sh

coolingtables/redshifts.dat:
	bash ../getEagleCoolingTable.sh
