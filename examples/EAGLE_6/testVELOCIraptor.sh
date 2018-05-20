#!/bin/bash

SCRIPT_PATH=/home/james/Work/VELOCIraptor-dev/VELOCIraptor-STF/stf/examples
TREEFROG_PATH=/home/james/Work/VELOCIraptor-STF/stf/bin
VELOCIRAPTOR_PATH=/home/james/Work/VELOCIraptor-STF/stf

# List each test that should be run
declare -a DM_TEST_LIST=(6dfof_dmonly_sub)
declare -a GAS_TEST_LIST=(6dfof_gas_sub)
INFO_FILE_TEXT="$TREEFROG_PATH/treefrog \n-s 2 -d 2 -F 0 -B 2 -m -1 -N 1 -T -1 \n2 # No. of outputs to compare"

# Run catalog comparison for DM only runs
for TEST in "${DM_TEST_LIST[@]}"
do

  # Remove old comparison files
  rm catcomp*

  # Create output directory
  OUTPUT=halo_1_mpi_1_$TEST
  
  if [ ! -d $OUTPUT ]; then mkdir $OUTPUT; fi

  # Run test using SWIFT + VELOCIraptor
  mpirun -n 1 ../swift_mpi -G -t 8 eagle_6.yml -x -n 5 -P StructureFinding:basename:./$OUTPUT/stf -P StructureFinding:config_file_name:./stf_input_$TEST.cfg

  # Run test using VELOCIraptor
  mpirun -np 1 $VELOCIRAPTOR_PATH/bin/stf-gasstarbh -I 2 -i eagle_0000 -C $VELOCIRAPTOR_PATH/vel_input_$TEST.cfg -o ./vel_outputs_new/vel_$TEST
 
  # Create info file for python comparison script
  echo -e $INFO_FILE_TEXT > infoFile_$TEST.txt
  echo "vel_$TEST vel_outputs_new/vel_$TEST" >> infoFile_$TEST.txt
  echo "stf_$TEST $OUTPUT/stf_0000.VELOCIraptor" >> infoFile_$TEST.txt

  # Run comparison script on VELOCIraptor output
  if python $SCRIPT_PATH/catalogcomparisontolerancecheck.py infoFile_$TEST.txt toIfile.txt
  then
    echo "Catalog comparison passed: "$TEST
  else
    echo "Catalog comparison failed: "$TEST
    exit 1
  fi
  
  echo "------------"

done

# Run catalog comparison for DM + GAS runs
for TEST in "${GAS_TEST_LIST[@]}"
do

  # Remove old comparison files
  rm catcomp*

  # Create output directory
  OUTPUT=halo_1_mpi_1_$TEST

  if [ ! -d $OUTPUT ]; then mkdir $OUTPUT; fi
  
  # Run test using SWIFT + VELOCIraptor
  mpirun -n 1 ../swift_mpi -s -G -t 8 eagle_6.yml -x -n 5 -P StructureFinding:basename:./$OUTPUT/stf -P StructureFinding:config_file_name:./stf_input_$TEST.cfg

  # Run test using VELOCIraptor
  mpirun -np 1 $VELOCIRAPTOR_PATH/bin/stf-gasstarbh -I 2 -i eagle_0000 -C $VELOCIRAPTOR_PATH/vel_input_$TEST.cfg -o ./vel_outputs_new/vel_$TEST
  
  # Create info file for python comparison script
  echo -e $INFO_FILE_TEXT > infoFile_$TEST.txt
  echo "vel_$TEST vel_outputs_new/vel_$TEST" >> infoFile_$TEST.txt
  echo "stf_$TEST $OUTPUT/stf_0000.VELOCIraptor" >> infoFile_$TEST.txt

  # Run comparison script on VELOCIraptor output
  if python $SCRIPT_PATH/catalogcomparisontolerancecheck.py infoFile_$TEST.txt toIfile.txt
  then
    echo "Catalog comparison passed: "$TEST
  else
    echo "Catalog comparison failed: "$TEST
    exit 1
  fi
  
  echo "------------"

done

exit $?
