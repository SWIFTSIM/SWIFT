#!/bin/bash

VELOCIRAPTOR_PATH=/cosma5/data/dp004/dc-will2/VELOCIraptor-STF/stf
SCRIPT_PATH=$VELOCIRAPTOR_PATH/examples
TREEFROG_PATH=$VELOCIRAPTOR_PATH/bin

# List each test that should be run
declare -a DM_TEST_LIST=(6dfof_dmonly_sub)
declare -a GAS_TEST_LIST=(6dfof_gas_sub)
INFO_FILE_TEXT="$TREEFROG_PATH/treefrog \n-s 2 -d 2 -F 0 -B 2 -m -1 -N 1 -T -1\n2 # No. of outputs to compare"

# Check for command line arguments
while getopts 'g:s:' opt ; do
  case $opt in
    g) RUN_DM=$OPTARG ;;
    s) RUN_GAS=$OPTARG ;;
  esac
done

# Run catalog comparison for DM only runs
if [ "$RUN_DM" = "1" ]; then
  for TEST in "${DM_TEST_LIST[@]}"
  do

    # Create output directory
    OUTPUT=halo_1_mpi_1_$TEST

    if [ ! -d $OUTPUT ]; then mkdir $OUTPUT; fi

    # Remove old comparison files
    rm catcomp*
    rm $OUTPUT/stf* 
    rm vel_outputs_new/vel_$TEST*

    # Run test using SWIFT + VELOCIraptor
    echo "Running: mpirun -np 1 ../swift_mpi -G -t 8 eagle_6.yml -x -n 5 -P StructureFinding:basename:./$OUTPUT/stf -P StructureFinding:config_file_name:./stf_input_$TEST.cfg -P Snapshots:basename:./eagle_dmonly"
    mpirun -np 1 ../swift_mpi -G -t 8 eagle_6.yml -x -n 5 -P StructureFinding:basename:./$OUTPUT/stf -P StructureFinding:config_file_name:./stf_input_$TEST.cfg -P Snapshots:basename:./eagle_dmonly

    # Run test using VELOCIraptor
    echo "Running: mpirun -np 1 $VELOCIRAPTOR_PATH/bin/stf-gas -I 2 -i eagle_dmonly_0000 -C $VELOCIRAPTOR_PATH/vel_input_$TEST.cfg -o ./vel_outputs_new/vel_$TEST"
    mpirun -np 1 $VELOCIRAPTOR_PATH/bin/stf-gas -I 2 -i eagle_dmonly_0000 -C $VELOCIRAPTOR_PATH/vel_input_$TEST.cfg -o ./vel_outputs_new/vel_$TEST

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
fi

# Run catalog comparison for DM + GAS runs
if [ "$RUN_GAS" = "1" ]; then
  for TEST in "${GAS_TEST_LIST[@]}"
  do

    # Create output directory
    OUTPUT=halo_1_mpi_1_$TEST

    if [ ! -d $OUTPUT ]; then mkdir $OUTPUT; fi

    # Remove old comparison files
    rm catcomp*
    rm $OUTPUT/stf* 
    rm vel_outputs_new/vel_$TEST*

    # Run test using SWIFT + VELOCIraptor
    echo "Running: mpirun -np 1 ../swift_mpi -s -G -t 8 eagle_6.yml -x -n 5 -P StructureFinding:basename:./$OUTPUT/stf -P StructureFinding:config_file_name:./stf_input_$TEST.cfg -P Snapshots:basename:./eagle_gas"
    mpirun -np 1 ../swift_mpi -s -G -t 8 eagle_6.yml -x -n 5 -P StructureFinding:basename:./$OUTPUT/stf -P StructureFinding:config_file_name:./stf_input_$TEST.cfg -P Snapshots:basename:./eagle_gas

    # Run test using VELOCIraptor
    echo "Running: mpirun -np 1 $VELOCIRAPTOR_PATH/bin/stf-gas -I 2 -i eagle_gas_0000 -C $VELOCIRAPTOR_PATH/vel_input_$TEST.cfg -o ./vel_outputs_new/vel_$TEST"
    mpirun -np 1 $VELOCIRAPTOR_PATH/bin/stf-gas -I 2 -i eagle_gas_0000 -C $VELOCIRAPTOR_PATH/vel_input_$TEST.cfg -o ./vel_outputs_new/vel_$TEST

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
fi

exit $?
