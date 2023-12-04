#!/bin/bash

# Generate the initial conditions if they are not present.
if [ ! -e SWIFT_ICs_L0150N0064.hdf5 ]
then
     echo "Fetching initial conditions for the small cosmological volume example..."
     ./getIC.sh
fi

if [ ! -e SWIFT_MHD_ICs_L0150N0064.hdf5 ]
then
     echo "Generating the Bfield"
     python ../makeIC.py SWIFT_ICs_L0150N0064.hdf5 SWIFT_MHD_ICs_L0150N0064.hdf5
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
#####################
case $WHAT in
	vep)
	  SCHEME_Nr=0
	;;
	odi)
	  SCHEME_Nr=1
	;;
	fdi)
	  SCHEME_Nr=2
	;;
	all)
	  SCHEME_Nr=( 0 1 2 )
	;;
	*)
	  echo $WHAT" wrong scheme"
	  exit 2
	;;
esac
ADDOPT="--with-hydro=minimal"
SCHEME_DIRS=("VeP_$DIRS" "ODI_$DIRS" "FDI_$DIRS")

for J in ${SCHEME_Nr[@]}
do
   echo $J
   ID=${SCHEME_ID[$J]}
   IDD=${SCHEME_IDD[$J]}
   DIR=${SCHEME_DIRS[$J]}
   if [ ! -e $DIR ]
   then 
      echo "Folder $DIR does not exist"
      mkdir $DIR
   fi
   cd $DIR
   if [ ! -e sw_$ID ]
   then 
      cur_dir=`pwd`
      cd ../../../../../
      pwd
      ./TestAllMHD.sh $IDD $ADDOPT
      cd $cur_dir
      cp ../../../../../sw_$ID .
   fi
   cat <<-EOF > ./run.sh
	#!/bin/bash
	# Run SWIFT
	./sw_$ID  --cosmology --hydro --self-gravity --threads=16 ../sw.yml &> output.log &
	
	# Plot the evolution
	#python3 ../plot_all.py 0 41 2>&1 > plot.log
	EOF
   chmod u+x ./run.sh
   ./run.sh &
   cd ..
done

