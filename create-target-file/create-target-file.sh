#!/bin/bash -l

# Assumes we run SLURM; provision for SGE will be added afterwards

# For use in several known cluster environments: [Durham cosma8] for now, 
#    perhaps [UCL's Myriad, Kathleen, Young] later on, 
#    [ARCHER2], etc.

# Minimal error checking throughout

# Workflow:
#  - this script is called from within a submitted jobscript.
#  - tries to determine architecture [cosma8 | Kathleen etc] 
#  - tries to determine jobscript directives relevant to how many
#    nodes, mpi ranks/tasks etc are requested, and some other runtime parameters
#  - in order to produce a [target.tgt] file, to be passed on to the SCOTCH module.


# Process the jobscript from which we are called, strip comments and gather 
#  switches / parameters relevant to SCOTCH, remove '=' from those lines

sed -e '/^\s\{0,\}#\ /d' $1 | \
  sed -n -r -e '/^#SBATCH/p;/\-\-with\-arch/{s/^.*\-\-with\-arch\=/\-\-with\-arch\=/;s/[ ].*$//;p}' | \
  sed 's/\=/\ /' | \
  sed 's/^#SBATCH\ //' > target-specs.txt

# Parse the lines in [target-specs.txt], determine SCOTCH-relevant parameters, 
#  consult cluster architecture, and if successful, produce the specs file [target.tgt];
#  if not, return some error, so we may return that error to the caller script.


