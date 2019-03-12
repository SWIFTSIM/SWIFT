.. VELOCIraptor stand alone 
   Folkert Nobels 12th October 2018

Stand alone VELOCIraptor configuration
======================================


.. toctree::    
   :maxdepth: 2    
   :hidden:    
   :caption: Contents: 
   
Besides running VELOCIraptor on the fly when using SWIFT, it is also possible
to run VELOCIraptor alone without using SWIFT. In this section we explain how 
VELOCIraptor can be run stand alone without using SWIFT.

Setting up VELOCIraptor
-----------------------

The first step is setting up VELOCIraptor, this requires us to download the 
git repository as::
  
  git clone https://github.com/pelahi/VELOCIraptor-STF

Similar to the SWIFT with VELOCIraptor configuration, we can use the 
master to analyse individual snapshots. We can use this branch
by doing::

  cd VELOCIraptor-STF
  git fetch

Again we need to configure VELOCIraptor::

  cmake . -DVR_USE_GAS=ON

In this case, we do not need the SWIFT interface, therefore we can drop
this option (disabled by default).

Compiling VELOCIraptor
----------------------

Compiling goes completely different as compared to the on the fly halo finder
configuration with SWIFT. In this case we can compile the code as::

  make 

After this an executable is created (``VELOCIraptor-stf/stf``).

Running VELOCIraptor on a Snapshot
----------------------------------

After the code is compile the next step is using VELOCIraptor on a single 
snapshot of a simulation. The code has several options which can be used, which
can be displayed by running a terminal command of an invalid letter like::

  ./stf -h

which gives the information about the usage of the command::

  VELOCIraptor/STF running with MPI. Number of mpi threads: 1
  VELOCIraptor/STF running with OpenMP. Number of openmp threads: 8
  USAGE:

  -C <configuration file (overrides other options)> 
  -I <input format [Gadget (Default) 1, HDF (if implemented)2, TIPSY 3, RAMSES 4, HDF 2, NCHILADA 5>
  -i <input file> 
  -s <number of files per output for gadget input 1 [default]>
  -Z <number of threads used in parallel read (1)>
  -o <output filename>
   ===== EXTRA OPTIONS FOR GADGET INPUT ====== 
  -g <number of extra sph/gas blocks for gadget>
  -s <number of extra star blocks for gadget>
  -b <number of extra bh blocks for gadget>
   ===== EXTRA OPTIONS REQUIRED FOR RAMSES INPUT ====== 
  -t <ramses snapnumber>

After this we can run VELOCIraptor on a snapshot as::
  
  ./stf -i input -o output -C configfile.txt
