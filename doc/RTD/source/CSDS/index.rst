Continuous Simulation Data Stream (CSDS)
========================================

The CSDS is a particle based output (e.g. snapshot) that takes into account the large difference of timescale.
If you have any questions, a slack channel is available for it in SWIFT's slack.

To run it, you will need to use the configuration option ``--enable-csds`` and the run time argument ``--csds``.
Currently the CSDS is implemented only for GEAR, Gadget2 and the default modules, but can be easily extended to the other schemes by adding the CSDS structure to the particles and implementing the IO functions (see ``src/hydro/Gadget2/hydro_part.h``, ``src/hydro/Gadget2/hydro_csds.c`` and ``src/hydro/Gadget2/hydro_csds.h``).
The main parameters of the CSDS are ``CSDS:delta_step`` and ``CSDS:index_mem_frac`` that define the time accuracy of the CSDS and the number of index files.
The first parameter defines the number of active steps that a particle is doing before writing and the second defines the total storage size of the index files as a fraction of the dump file.

For reading, the python wrapper is available through the configuration option ``--with-python``.
I recommend running the SedovBlast_3D with the CSDS and then using the example ``csds/examples/reader_example.py``.
This file is kept up to date with the most recent changes and includes a call to all the existing functions.
If you need some extra information, a doc string is provided for the class ``csds.Reader`` and all its methods.

If you wish to obtain a snapshot from the CSDS, a script is available in ``csds/examples/create_snapshot.py``.
