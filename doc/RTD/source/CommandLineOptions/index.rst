.. Command line options
   Matthieu Schaller, 21st October 2018

.. _cmdline-options:

Command line options
====================

SWIFT requires a number of runtime options to run and get any sensible output.
For instance, just running the ``swift`` binary will not use any SPH or gravity;
the particles will just sit still!

Below is a list of the command line options and when they should be used. The same list
can be found by typing ``./swift -h``::

    -h, --help                        show this help message and exit

  Simulation options:

    -b, --feedback                    Run with stars feedback.
    -c, --cosmology                   Run with cosmological time integration.
    --temperature                     Run with temperature calculation. 
    -C, --cooling                     Run with cooling (also switches on --with-temperature).
    -D, --drift-all                   Always drift all particles even the ones
                                      far from active particles. This emulates
                                      Gadget-[23] and GIZMO's default behaviours.
    -F, --star-formation	      Run with star formation.
    -g, --external-gravity            Run with an external gravitational potential.
    -G, --self-gravity                Run with self-gravity.
    -M, --multipole-reconstruction    Reconstruct the multipoles every time-step.
    -s, --hydro                       Run with hydrodynamics.
    -S, --stars                       Run with stars.
    -x, --velociraptor                Run with structure finding.
    --limiter                         Run with time-step limiter.

  Control options:

    -a, --pin                         Pin runners using processor affinity.
    -d, --dry-run                     Dry run. Read the parameter file, allocates
                                      memory but does not read the particles
                                      from ICs. Exits before the start of time
                                      integration. Checks the validity of
                                      parameters and IC files as well as memory
                                      limits.
    -e, --fpe                         Enable floating-point exceptions (debugging
                                      mode).
    -f, --cpu-frequency=<str>         Overwrite the CPU frequency (Hz) to be
                                      used for time measurements.
    -n, --steps=<int>                 Execute a fixed number of time steps.
                                      When unset use the time_end parameter
                                      to stop.
    -o, --output-params=<str>         Generate a default output parameter
                                      file.
    -P, --param=<str>                 Set parameter value, overiding the value
                                      read from the parameter file. Can be used
                                      more than once {sec:par:value}.
    -r, --restart                     Continue using restart files.
    -t, --threads=<int>               The number of threads to use on each MPI
                                      rank. Defaults to 1 if not specified.
    -T, --timers=<int>                Print timers every time-step.
    -v, --verbose=<int>               Run in verbose mode, in MPI mode 2 outputs
                                      from all ranks.
    -y, --task-dumps=<int>            Time-step frequency at which task analysis
                                      files and/or tasks are dumped.
    -Y, --threadpool-dumps=<int>      Time-step frequency at which threadpool
                                      tasks are dumped.
