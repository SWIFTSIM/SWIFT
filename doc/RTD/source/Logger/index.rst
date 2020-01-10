Logger Output
=============

The logger is a particle based output (e.g. snapshot) that takes into account the large difference of timescale.
If you have any question, a slack channel is available for it in SWIFT's slack.

To run it, you will need to use the configuration option ``--enable-logger``.
Currently the logger is implemented only for Gadget2 and the default gravity / stars, but can be easily extended to the other schemes by adding the logger structure to the particles (see ``src/hydro/Gadget2/hydro_part.h``).
The main parameters of the logger are ``Logger:delta_step`` and ``Logger:index_mem_frac`` that define the time accuracy of the logger and the number of index files.
The first parameter defines the number of active steps that a particle is doing before writing and the second defines the total storage size of the index files as function of the dump file.

Unfortunately, the API is not really developed yet. Therefore if you wish to dump another field, you will need to trick the logger by replacing a field in the ``logger_log_part`` function.

For reading, the python wrapper is available through the configuration option ``--with-python``. Once compiled, you will be able to use the file ``logger/examples/reader_example.py``.
The first argument is the basename of the index file and the second one is the time requested.
During the first reading, the library is manipulating the dump file and therefore it should not be killed and may take a bit more time than usual.
