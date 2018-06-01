.. Runtime Options
   Josh Borrow, 5th April 2018

Runtime Options
===============

SWIFT requires a number of runtime options to run and get any sensible output.
For instance, just running the ``swift`` binary will not use any SPH or gravity;
the particles will just sit still!

Below is a list of the runtime options and when they should be used. The same list
can be found by typing ``./swift -h``.

+ ``-a``: Pin runners using processor affinity.
+ ``-c``: Run with cosmological time integration.
+ ``-C``: Run with cooling.
+ ``-d``: Dry run. Read the parameter file, allocate memory but does not read
  the particles from ICs and exit before the start of time integration. Allows
  user to check validity of parameter and IC files as well as memory limits.
+ ``-D``: Always drift all particles even the ones far from active particles.
  This emulates Gadget-[23] and GIZMO's default behaviours.
+ ``-e``: Enable floating-point exceptions (debugging mode).
+ ``-f``: {int} Overwrite the CPU frequency (Hz) to be used for time measurements.
+ ``-g``: Run with an external gravitational potential.
+ ``-G``: Run with self-gravity.
+ ``-M``: Reconstruct the multipoles every time-step.
+ ``-n``: {int} Execute a fixed number of time steps. When unset use the
  time_end parameter to stop.
+ ``-o``: {str} Generate a default output parameter file.
+ ``-P``: {sec:par:val} Set parameter value and overwrites values read from the
  parameters file. Can be used more than once.
+ ``-s``: Run with hydrodynamics.
+ ``-S``: Run with stars.
+ ``-t``: {int} The number of threads to use on each MPI rank. Defaults to 1 if
  not specified.
+ ``-T``: Print timers every time-step.
+ ``-v``: [12] Increase the level of verbosity: 1, MPI-rank 0 writes, 2, All
  MPI-ranks write.
+ ``-y``: {int} Time-step frequency at which task graphs are dumped.
+ ``-Y``: {int} Time-step frequency at which threadpool tasks are dumped.
+ ``-h``: Print a help message and exit.
