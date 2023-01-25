.. Runtime Options
   Josh Borrow, 5th April 2018

Runtime Options and Parameter Files
===================================

SWIFT requires a number of runtime options to run and get any sensible output.
For instance, just running the ``swift`` binary will not use any SPH or gravity;
the particles will just sit still!

A list of  command line options can be found by running the compiled binary with
the ``-h`` or ``--help`` flag:

.. code-block:: bash

   ./swift --help


You will also need to specify a number of runtime parameters that are dependent 
on your compile-time configuration in a parameter file. A list of all of these 
parameters can be found in ``examples/parameter_example.yml``, and you can check 
out examples in the ``examples/`` directory.
