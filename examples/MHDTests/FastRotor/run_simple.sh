#!/bin/bash

# Generate the initial conditions if they are not present.
if [ ! -e FastRotor.hdf5 ]
then
    echo "Fetching Glass Files..."
    ./getGlass.sh
    echo "Generating the ICs"
    python ./makeIC_schemes.py 
fi

SCHEME_ID=("VeP" "ODI" "FDI")
SCHEME_IDD=("vp" "odi" "fdi")
####################
case $# in
	1)
	   echo "Using Dirs as TEST"
	   DIRS="TEST"
	   WHAT=$1
        ;;
	2) 
	   DIRS=$2
	   WHAT=$1
	   echo "Using Dirs as $DIRS"
	;;
        *)
	   echo "Usage $0 [which] [DIRECTORY]"
	   echo "[what scheme]:"
	   echo "vep: vector potentials"
           echo "odi: Oresti's direct induction"
           echo "fdi: simple direct induction"
           echo "[FOLDER_TAIL]:"
           echo "trailer naming of folders"
	   echo ""
	   exit 
	;;
esac
if [ ! -e $DIRS ]
then 
   echo "Folder $DIRS does not exist"
   mkdir $DIRS
fi
cd $DIRS
if [ ! -e ../$WHAT ]
then 
   echo "$WHAT does not exist!"
   echo "CHECK!"
   exit
fi
cp ../$WHAT ./sw_run
cat <<-EOF > ./run.sh
	#!/bin/bash
	# Run SWIFT
	./sw_run --hydro --threads=12 ../FastRotor_simple.yml 2>&1 > out.log 
	
	# Plot the evolution
	python3 ../plot_schemes.py 0 16 2>&1 > plot.log
	EOF
chmod u+x ./run.sh
./run.sh &
cd ..

