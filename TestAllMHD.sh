#!/bin/bash

echo "Testing MHD schemes"

LOG_FILE="MHD_Pass.log"
rm $LOG_FILE

SCHEME=("vector-potential" "direct-induction" "direct-induction-fede" "none" "direct-induction")
SCHEME_ID=("VeP" "ODI" "FDI" "MDI" "ODI2")
SCHEME_NAME=("Vector Potential" "Direct Induction (Orestis)" "Direct Induction (Fede)" "Direct Induction (Matthieu)" "Direct Induction 2")
# VP ~ Vector Potential
# DO ~ Direct Induction Orestis
# DF ~ Direct Induction Fede
# DM ~ Direct Induction Matthieu
# D2 ~ Direct Induction Orestis 2

#BASE_CONF="--with-kernel=wendland-C4 --disable-hand-vec"
BASE_CONF=""
GIT_BBASE=(` git status -b | head -n 1 `)
GIT_BBASE=${GIT_BBASE[2]}
case $1 in
   all)
	varN=( 0 1 2 3 4 )
	echo ${varN[@]}
	echo "testing all schemes: "${SCHEME_NAME[*]}
   ;;
   vp)
	varN=0
	echo "scheme chosen: "${SCHEME_NAME[$varN]}
   ;;
   odi)
	varN=1
	echo "scheme chosen: "${SCHEME_NAME[$varN]}
   ;;
   fdi)
	varN=2
	echo "scheme chosen: "${SCHEME_NAME[$varN]}
   ;;
   mdi)
	varN=3
	echo "scheme chosen: "${SCHEME_NAME[$varN]}
	BASE_CONF=$BASE_CONF+"--with-sph=minimal"
	git switch MHD_canvas_Matthieu
	git pull 
   ;;
   odi2)
	varN=4
	echo "scheme chosen: "${SCHEME_NAME[$varN]}
	git switch karapiperis/ODI2
	git pull
   ;;
   *)
	echo "Usage $0 [what] [configure extra parameters]"
	echo "[what]:"
	echo "all: compile all schemes"
	echo " vp: vector potenial scheme"
	echo "odi: Direct Induction Orestis scheme"
	echo "fdi: Direct Induction Federico scheme"
	echo "mdi: Direct Induction Matthieu scheme"
	echo "odi: Direct Induction Orestis 2 scheme"
	exit
   ;;
esac

git status -bs > $LOG_FILE
echo ${varN[*]}
for J in ${varN[@]}
do
   echo $J
   ID=${SCHEME_ID[$J]}
   rm sw_$ID*  MHD_$ID.log
   echo "Compiling "${SCHEME_NAME[$J]}" version" 

   ./configure --with-spmhd=${SCHEME[$J]} $BASE_CONF $2 > MHD_$ID.log
   make -j 32 >> MHD_$ID.log
   mv swift sw_$ID
   mv swift_mpi sw_$ID"_mpi"
	
   echo "MHD: "$ID" OK!" >> $LOG_FILE
   git switch $GIT_BBASE
done
exit
