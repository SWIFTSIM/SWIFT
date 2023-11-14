#!/bin/bash

# Generate the initial conditions if they are not present.
if [ ! -e MagneticBlastWave_LR.hdf5 ]
then
    echo "Fetching Glass Files..."
    ./getGlass.sh
    echo "Generating the ICs"
    python ./makeIC.py 
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
      cd ../../../../
      pwd
      ./TestAllMHD.sh $IDD
      cd $cur_dir
      cp ../../../../sw_$ID .
   fi
   #../../../../sw_VeP --hydro --threads=16 ../BW_schemes.yml 2>&1 > out.log 
   cat <<-EOF > ./run.sh
	#!/bin/bash
	# Run SWIFT
	./sw_$ID --hydro --threads=16 ../BW_schemes.yml 2>&1 > out.log 
	
	# Plot the temperature evolution
	python3 ../plot_all.py 0 41 > plot.log
	EOF
   chmod u+x ./run.sh
   ./run.sh &
   cd ..
done

