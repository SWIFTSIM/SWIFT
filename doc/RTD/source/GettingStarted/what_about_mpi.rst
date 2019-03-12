.. What about MPI? Running SWIFT on more than one node
   Josh Borrow, 5th April 2018

What about MPI? Running SWIFT on more than one node
===================================================

After compilation, you will be left with two binaries. One is called ``swift``,
and the other ``swift_mpi``. Current wisdom is to run ``swift`` if you are only
using one node (i.e. without any interconnect), and one MPI rank per NUMA
region using ``swift_mpi`` for anything larger. You will need some GADGET-2
HDF5 initial conditions to run SWIFT, as well as a compatible yaml
parameter file.
