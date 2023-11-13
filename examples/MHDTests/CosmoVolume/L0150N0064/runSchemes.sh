#!/bin/bash

 # Generate the initial conditions if they are not present.
if [ ! -e SWIFT_ICs_L0150N0064.hdf5 ]
then
    echo "Fetching initial conditions for the small cosmological volume example..."
    ./getIC.sh
    echo "Generating the Bfield"
    python ../makeIC.py SWIFT_ICs_L0150N0064.hdf5 SWIFT_MHD_ICs_L0150N0064.hdf5
fi

SCHEME_ID=("VeP" "ODI" "FDI")
SCHEME_IDD=("vp" "odi" "fdi")
SCHEME_Nr=( 0 1 2 )
case $# in
	1)
	   echo "Using Dirs as TEST"
	   DIRS="TEST"
        ;;
	2) 
	   DIRS=$2
	   echo "Using Dirs as $DIRS"
	;;
        *)
	   echo "Usage $0 [what] [DIRECTORY]"
	   echo ""
	   exit 
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
      cd ../../../../../
      pwd
      ./TestAllMHD.sh $IDD
      cd $cur_dir
      cp ../../../../../sw_$ID .
   fi
	#sw_$ID --cosmology --hydro --self-gravity --threads=16 ../sw_param_L0150N0064.yml 2>&1 | tee output.log
   cat <<-EOF > ./run.sh
	#!/bin/bash
	sw_$ID --cosmology --hydro --self-gravity --threads=16 ../sw_param_L0150N0064.yml 2>&1 > output.log
	EOF
   chmod u+x ./run.sh
   ./run.sh &
   cd ..
done

