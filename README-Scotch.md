Information on how to run SWIFT with Scotch mapping and the test environment used on Cosma 8. Code has been tested with Scotch version > 7.0

Last update 18th August 2023.


Obtaining Scotch as not installed system wide on Cosma 8.
----------------

**Scotch** is publicly available under the CeCILL-C free software license, as described [here](https://gitlab.inria.fr/scotch/scotch/blob/master/LICENSE_en.txt). The license itself is available [here](https://gitlab.inria.fr/scotch/scotch/-/blob/master/doc/CeCILL-C_V1-en.txt).

To use the lastest version of **Scotch**, please clone the master branch:

    git clone git@gitlab.inria.fr:scotch/scotch.git

Tarballs of the **Scotch** releases are available [here](https://gitlab.inria.fr/scotch/scotch/-/releases).

Instructions for installing locally on Cosma 8, please ammend as appropriate. 
----------------
_Environment_
```
    module load cosma/2018 python/3.6.5 intel_comp/2022.1.2 compiler openmpi/4.1.1 fftw/3.3.9 parallel_hdf5/1.12.0 parmetis/4.0.3-64bit gsl/2.5
    module load cmake
    module load bison
```
Navigate to the Scotch directory and carry out the following commands

```
    mkdir build 
    cd build
    cmake -DCMAKE_INSTALL_PREFIX=/path-to-install-dir ..
    make -j5
    make install
```

Configure SWIFT with Scotch
----------------

Follow the usual installation [instructions](https://gitlab.cosma.dur.ac.uk/swift/swiftsim/-/blob/master/INSTALL.swift) but if Scotch installed locally the added `--with-scotch=\path-to-scotch` flag will need to be passed to `./configure`

Running with Scotch
----------------

Scotch decomposes the SWIFT spatial domain and maps it to the available compute - taking into consideration the communication cost between components of the architecture. In order for this to be carried out the user needs to generate an appropriate architecture file. This architecture file should mirror the set up of the cluster being used. Scotch provides optimised architecture files which capture most HPC set ups. As we will be targetting NUMA regions on Cosma 8 we have modelled the architecture as a `tleaf` structure. 

In the following examples it is assumed that one mpi rank is mapped to each Cosma 8 NUMA region. This enforces that `cpus-per-task=16` is defined in the SLURM submission script. The Cosma 8 nodes consist of 8 NUMA regions per node, with 4 NUMA regions per socket. Example `tleaf`files for various setups are given below, where the intrasocket communication cost between NUMA regions is set at _5_, intranode but across sockets is set at _10_ and the internode cost is set at _1000_. These weightings are estimated values but have been shown to give satisfactory results in the testcases explored.

| Number of nodes | Number of MPI ranks | tleaf                   |
| --------------- | ------------------- | ----------------------- |
| 1               | 2                   | tleaf 1 2 5             |
| 1               | 8                   | tleaf 2 2 10 4 5        |
| 2               | 16                  | tleaf 3 2 1000 2 10 4 5 |
| 4               | 32                  | tleaf 3 4 1000 2 10 4 5 |
| 8               | 64                  | tleaf 3 8 1000 2 10 4 5 |

The user needs to define this tleaf structure and save it as `target.tgt` in the directory they will run SWIFT from. Ongoing work focuses on automatically generating this target architecture upon run time. 


With OpenMPI the `mpirun` option `--map-by numa` has been found to be optimal. This is in contrast to previously suggested `--bind-to none` on the cosma8 [site](https://www.dur.ac.uk/icc/cosma/support/cosma8/).

Scotch details
----------------

Scotch carries out the mapping using various strategies which are outlined in the documentation. The Scotch strategy is passed to the Mapping functions. For the runs carried out here it was found that the global flag `SCOTCH_STRATBALANCE` and a imbalance ratio of `0.05` worked best. These values are passed to `SCOTCH_stratGraphMapBuild`. 

One issue with Scotch is that when the number of mpi ranks is comparable to the dimensionality for the SWIFT system the optimally mapping strategy doesn't map to all available NUMA regions. At present this isn't handled explicity in the code and the paritition reverts to a vectorised or previous partitioning.

The SWIFT edge and vertex weights are estimated in the code, however edge weights are not symmetric. This causes an issue with SWIFT therefore before building the SCOTCH Graph the edge weigths are updated to equal to the sum of the two associated edge weights.  

