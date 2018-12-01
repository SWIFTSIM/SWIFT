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
VELOCIraptor. This can be done by cloning the repository on GitHub_::

  git clone https://github.com/pelahi/VELOCIraptor-STF

Currently the best version that works with SWIFT is the swift-interface branch
of VELOCIraptor, to get this branch use::

  cd VELOCIraptor-STF 
  git fetch 
  git checkout swift-interface

To get the default that works with SWIFT simply copy the SWIFT template file in
the ``Makefile.config``::

  cd stf 
  cp Makefile.config.SWIFT-template Makefile.config

Depending on your compiler you want to change the first 20 lines of your
``Makefile.config`` to work with your compiler and whether you want to use MPI
or not. 


Compiling VELOCIraptor
----------------------

After we downloaded the files and made a configuration file we can compile
VELOCIraptor as follows::

  make lib 
  make libstf

After the compilation of your code, there is an additional folder created in
the ``VELOCIraptor-stf/stf`` directory called ``lib`` this directory has the
library of VELOCIraptor and is required to run SWIFT with
VELOCIraptor. Note that VELOCIraptor needs a serial version of the
HDF5 library, not a parallel build.

Compiling SWIFT
---------------
The next part is compiling SWIFT with VELOCIraptor and assumes you already
downloaded SWIFT from the GitLab_, this can be done by running::

  ./autogen.sh 
  ./configure --with-velociraptor=/path/to/VELOCIraptor-STF/stf/lib 
  make 

In which ``./autogen.sh`` only needs to be run once after the code is cloned
from the GitLab_, and ``/path/to/`` is the path to the ``VELOCIraptor-STF``
directory on your machine. In general ``./configure`` can be run with other
options as desired. After this we can run SWIFT with VELOCIraptor, but for this
we first need to add several lines to the yaml file of our simulation::

    
  #structure finding options
  StructureFinding:
  config_file_name:     stf_input_6dfof_dmonly_sub.cfg
  basename:             ./stf
  output_time_format:   1
  scale_factor_first:   0.02
  delta_time:           1.02

In which we specify the ``.cfg`` file that is used by VELOCIraptor and the 
other parameters which SWIFT needs to use. In the case of 
the Small Cosmological Volume DMO example we can run a simulation with halo
finder as::

  cd examples/SmallCosmoVolume_DM 
  ../swift --cosmology --hydro --self-gravity --velociraptor --threads=8 small_cosmo_volume_dm.yml

Which activates the VELOCIraptor interface.


.. _GitHub: https://github.com/pelahi/VELOCIraptor-STF
.. _GitLab: https://gitlab.cosma.dur.ac.uk/swift/swiftsim
