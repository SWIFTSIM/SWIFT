Logger Output
=============

The logger is a particle based output (e.g. snapshot) that takes into account the large difference of timescale.
If you have any questions, a slack channel is available for it in SWIFT's slack.

To run it, you will need to use the configuration option ``--enable-logger`` and the run time argument ``--logger``.
Currently the logger is implemented only for GEAR, Gadget2 and the default modules, but can be easily extended to the other schemes by adding the logger structure to the particles and implementing the IO functions (see ``src/hydro/Gadget2/hydro_part.h``, ``src/hydro/Gadget2/hydro_logger.c`` and ``src/hydro/Gadget2/hydro_logger.h``).
The main parameters of the logger are ``Logger:delta_step`` and ``Logger:index_mem_frac`` that define the time accuracy of the logger and the number of index files.
The first parameter defines the number of active steps that a particle is doing before writing and the second defines the total storage size of the index files as a fraction of the dump file.

For reading, the python wrapper is available through the configuration option ``--with-python``.
I recommend running the SedovBlast_3D with the logger and then using the example ``logger/examples/reader_example.py``.
This file is kept up to date with the most recent changes and includes a call to all the existing functions.
If you need some extra information, a doc string is provided for the class ``logger.Reader`` and all its methods.

If you wish to obtain a snapshot from the logger, a script is available in ``logger/examples/create_snapshot.py``.
