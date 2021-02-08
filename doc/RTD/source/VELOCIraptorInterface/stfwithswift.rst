.. SWIFT with VELOCIraptor
   Folkert Nobels 12th October 2018


Configuring SWIFT with VELOCIraptor
===================================

.. toctree::    
   :maxdepth: 2    
   :hidden:    
   :caption: Contents:

In the following three paragraphs we will explain how to setup VELOCIraptor,
how to compile it and how to compile SWIFT with VELOCIraptor. 


Setting up VELOCIraptor
-----------------------

Before we can run SWIFT with VELOCIraptor we first need to download
VELOCIraptor. This can be done by cloning the repository on GitHub::

  git clone https://github.com/ICRAR/VELOCIraptor-STF.git

The SWIFT interface is in the master branch of VELOCIraptor so nothing is more
is needed besides fetching the latest version of the `NBodyLib` that the code
relies upon::
  
  cd VELOCIraptor-STF 
  git submodule update --init --recursive

To get VELOCIraptor working with SWIFT simply use::

  mkdir build
  cd build
  cmake ../ -DVR_USE_HYDRO=ON -DVR_USE_SWIFT_INTERFACE=ON -DCMAKE_CXX_FLAGS="-fPIC" -DCMAKE_BUILD_TYPE=Release 

The first parameter activates the processing of gas, stars and black holes. It
can be omitted for simulations evolving only dark matter.

If you wish to run swift without MPI, you will need to add ``-DVR_MPI=OFF``.

Compiling VELOCIraptor
----------------------

After we downloaded the files and made a configuration file we can compile
VELOCIraptor as follows::

  make -j

After the compilation of your code, you will find a static library ``libvelociraptor.a``,
that is required to run SWIFT with VELOCIraptor.
Note that VELOCIraptor needs a serial version of the
HDF5 library, not a parallel build.

Compiling SWIFT
---------------
The next part is compiling SWIFT with VELOCIraptor and assumes you already
downloaded SWIFT from the GitLab repository. This can be done by running

.. code:: bash
  
  ./autogen.sh 
  ./configure --with-velociraptor=/path/to/VELOCIraptor-STF/build/src
  make 

In which ``./autogen.sh`` only needs to be run once after the code is cloned
from the GitLab_, and ``/path/to/`` is the path to the ``VELOCIraptor-STF``
directory on your machine. Note that this path must be absolute (i.e. relative
to ``/``. In general ``./configure`` can be run with other options as desired.
After this we can run SWIFT with VELOCIraptor, but for this we first need to
add several lines to the yaml file of our simulation

    
.. code:: YAML

   StructureFinding:      
     config_file_name:     vrconfig_3dfof_subhalos_SO_hydro.cfg
     basename:             haloes
     scale_factor_first:   0.02
     delta_time:           1.02

In which we specify the ``.cfg`` file that is used by VELOCIraptor and the 
other parameters which SWIFT needs to use. In the case of 
the Small Cosmological Volume DMO example we can run a simulation with halo
finder as::

  cd examples/SmallCosmoVolume_DM 
  ../swift --cosmology --hydro --self-gravity --velociraptor --threads=8 small_cosmo_volume_dm.yml

Which activates the VELOCIraptor interface.


.. _GitLab: https://gitlab.cosma.dur.ac.uk/swift/swiftsim
