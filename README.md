<a name="logo"/>
<div align="center">
<a href="https://www.swiftsim.com/" target="_blank">
<img src="https://swift.strw.leidenuniv.nl/SWIFT_banner.jpg" alt="SWIFT banner" width="1016" height="242"></img>
</a>
</div>

SWIFT: SPH WIth Fine-grained inter-dependent Tasking
====================================================

[![Build Status](https://gitlab.cosma.dur.ac.uk/jenkins/job/GNU%20SWIFT%20build/badge/icon)](https://gitlab.cosma.dur.ac.uk/jenkins/job/GNU%20SWIFT%20build/)

SWIFT is a gravity and SPH solver designed to run cosmological simulations
on peta-scale machines, scaling well up to 10's of thousands of compute
node.

More general information about SWIFT is available on the project
[webpages](http://www.swiftsim.com).

For information on how to _run_ SWIFT, please consult the onboarding guide
available [here](https://swift.strw.leidenuniv.nl/onboarding.pdf). This includes
dependencies, and a few examples to get you going.

We suggest that you use the latest release branch of SWIFT, rather than the
current master branch as this will change rapidly. We do, however, like to
ensure that the master branch will build and run.

This GitHub repository is designed to be an issue tracker, and a space for
the public to submit patches through pull requests. It is synchronised with
the main development repository that is available on the
[ICC](http://icc.dur.ac.uk)'s GitLab server which is available
[here](https://gitlab.cosma.dur.ac.uk/swift/swiftsim).

Please feel free to submit issues to this repository, or even pull
requests. We will try to deal with them as soon as possible, but as the
core development team is quite small this could take some time.

Disclaimer
----------

We would like to emphasise that SWIFT comes without any warranty of accuracy,
correctness or efficiency. As mentioned in the license, the software comes
`as-is` and the onus is on the user to get meaningful results. Whilst the
authors will endeavour to answer questions related to using the code, we
recommend users build and maintain their own copies. This documentation contains
the most basic information to get started. Reading it and possibly also the
source code is the best way to start running simulations.

The users are responsible to understand what the code is doing and for the
results of their simulation runs.

Note also that the values of the parameters given in the examples are only
indicative. We recommend users experiment by themselves and a campaign of
experimentation with various values is highly encouraged. Each problem will
likely require different values and the sensitivity to the details of the
physical model is something left to the users to explore.

Acknowledgment & Citation
-------------------------

The SWIFT code was last described in this paper:
https://ui.adsabs.harvard.edu/abs/2023arXiv230513380S.  The core solver, the
numerical methods as well as many extensions where described there. We ask users
running SWIFT for their research to please cite this paper when they present
their results.

In order to keep track of usage and measure the impact of the software, we
kindly ask users publishing scientific results using SWIFT to add the following
sentence to the acknowledgment section of their papers:

"The research in this paper made use of the SWIFT open-source
simulation code (http://www.swiftsim.com, Schaller et al. 2018)
version X.Y.Z."

with the version number set to the version used for the simulations and the
reference pointing to the ASCL entry of the code: https://ascl.net/1805.020.



Contribution Guidelines
-----------------------

The SWIFT source code uses a variation of the 'Google' formatting style.
The script 'format.sh' in the root directory applies the clang-format-13
tool with our style choices to all the SWIFT C source file. Please apply
the formatting script to the files before submitting a pull request.

Please check that the test suite still runs with your changes applied before
submitting a pull request and add relevant unit tests probing the correctness
of new modules. An example of how to add a test to the suite can be found by
considering the tests/testGreeting case.

Any contributions that fail any of the automated tests will not be accepted.
Contributions that include tests of the proposed modules (or any current ones!)
are highly encouraged.

Runtime parameters
------------------

```
 Welcome to the cosmological hydrodynamical code
    ______       _________________
   / ___/ |     / /  _/ ___/_  __/
   \__ \| | /| / // // /_   / /
  ___/ /| |/ |/ // // __/  / /
 /____/ |__/|__/___/_/    /_/
 SPH With Inter-dependent Fine-grained Tasking

 Version : 1.0.0
 Website: www.swiftsim.com
 Twitter: @SwiftSimulation

See INSTALL.swift for install instructions.

Usage: swift [options] [[--] param-file]
   or: swift [options] param-file
   or: swift_mpi [options] [[--] param-file]
   or: swift_mpi [options] param-file

Parameters:

    -h, --help                        show this help message and exit

  Simulation options:

    -b, --feedback                    Run with stars feedback.
    -c, --cosmology                   Run with cosmological time integration.
    --temperature                     Run with temperature calculation.
    -C, --cooling                     Run with cooling (also switches on --temperature).
    -D, --drift-all                   Always drift all particles even the ones
                                      far from active particles. This emulates
                                      Gadget-[23] and GIZMO's default behaviours.
    -F, --star-formation              Run with star formation.
    -g, --external-gravity            Run with an external gravitational potential.
    -G, --self-gravity                Run with self-gravity.
    -M, --multipole-reconstruction    Reconstruct the multipoles every time-step.
    -s, --hydro                       Run with hydrodynamics.
    -S, --stars                       Run with stars.
    -B, --black-holes                 Run with black holes.
    -k, --sinks                       Run with sink particles.
    -u, --fof                         Run Friends-of-Friends algorithm to
                                      perform black hole seeding.
    --lightcone                       Generate lightcone outputs.
    -x, --velociraptor                Run with structure finding.
    --line-of-sight                   Run with line-of-sight outputs.
    --limiter                         Run with time-step limiter.
    --sync                            Run with time-step synchronization
                                      of particles hit by feedback events.
    --csds                            Run with the Continuous Simulation Data
                                      Stream (CSDS).
    -R, --radiation                   Run with radiative transfer.
    --power                           Run with power spectrum outputs.

  Simulation meta-options:

    --quick-lyman-alpha               Run with all the options needed for the
                                      quick Lyman-alpha model. This is equivalent
                                      to --hydro --self-gravity --stars --star-formation
                                      --cooling.
    --eagle                           Run with all the options needed for the
                                      EAGLE model. This is equivalent to --hydro
                                      --limiter --sync --self-gravity --stars
                                      --star-formation --cooling --feedback
                                      --black-holes --fof.
    --gear                            Run with all the options needed for the
                                      GEAR model. This is equivalent to --hydro
                                      --limiter --sync --self-gravity --stars
                                      --star-formation --cooling --feedback.
    --agora                           Run with all the options needed for the
                                      GEAR model. This is equivalent to --hydro
                                      --limiter --sync --self-gravity --stars
                                      --star-formation --cooling --feedback.
                                      
  Control options:

    -a, --pin                         Pin runners using processor affinity.
    --nointerleave                    Do not interleave memory allocations across
                                      NUMA regions.
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
    -o, --output-params=<str>         Generate a parameter file with the options
                                      for selecting the output fields.
    -P, --param=<str>                 Set parameter value, overiding the value
                                      read from the parameter file. Can be used
                                      more than once {sec:par:value}.
    -r, --restart                     Continue using restart files.
    -t, --threads=<int>               The number of task threads to use on each
                                      MPI rank. Defaults to 1 if not specified.
    --pool-threads=<int>              The number of threads to use on each MPI
                                      rank for the threadpool operations.
                                      Defaults to the numbers of task threads
                                      if not specified.
    -T, --timers=<int>                Print timers every time-step.
    -v, --verbose=<int>               Run in verbose mode, in MPI mode 2 outputs
                                      from all ranks.
    -y, --task-dumps=<int>            Time-step frequency at which task graphs
                                      are dumped.
    --cell-dumps=<int>                Time-step frequency at which cell graphs
                                      are dumped.
    -Y, --threadpool-dumps=<int>      Time-step frequency at which threadpool
                                      tasks are dumped.
    --dump-tasks-threshold=<flt>      Fraction of the total step's time spent
                                      in a task to trigger a dump of the task plot
                                      on this step

See the file examples/parameter_example.yml for an example of parameter file.
```
