SWIFT: SPH WIth Fine-grained inter-dependent Tasking
====================================================

[![Build Status](https://gitlab.cosma.dur.ac.uk/jenkins/job/GNU%20SWIFT%20build/badge/icon)](https://gitlab.cosma.dur.ac.uk/jenkins/job/GNU%20SWIFT%20build/)

SWIFT is a gravity and SPH solver designed to run cosmological simulations
on peta-scale machines, scaling well up to 10's of thousands of compute
node.

More general information about SWIFT is available on the project
[webpages](http://www.swiftsim.com).

For information on how to _run_ SWIFT, please consult the onboarding guide
available [here](http://www.swiftsim.com/onboarding.pdf). This includes
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

Contribution Guidelines
-----------------------

The SWIFT source code uses a variation of the 'Google' formatting style.
The script 'format.sh' in the root directory applies the clang-format-5.0
tool with our style choices to all the SWIFT C source file. Please apply
the formatting script to the files before submitting a pull request.

Please check that the test suite still runs with your changes applied before
submitting a pull request and add relevant unit tests probing the correctness
of new modules. An example of how to add a test to the suite can be found by
considering the tests/testGreeting case.

Any contributions that fail any of the automated tests will not be accepted.
Contributions that include tests of the proposed modules (or any current ones!)
are highly encouraged.

```
 Welcome to the cosmological hydrodynamical code
    ______       _________________
   / ___/ |     / /  _/ ___/_  __/
   \__ \| | /| / // // /_   / /   
  ___/ /| |/ |/ // // __/  / /    
 /____/ |__/|__/___/_/    /_/     
 SPH With Inter-dependent Fine-grained Tasking

 Website: www.swiftsim.com
 Twitter: @SwiftSimulation

See INSTALL.swift for install instructions.

Usage: swift [OPTION]... PARAMFILE
       swift_mpi [OPTION]... PARAMFILE

Valid options are:
  -a                Pin runners using processor affinity.
  -c                Run with cosmological time integration.
  -C                Run with cooling.
  -d                Dry run. Read the parameter file, allocate memory but does not read
                    the particles from ICs and exit before the start of time integration.
                    Allows user to check validity of parameter and IC files as well as memory limits.
  -D                Always drift all particles even the ones far from active particles. This emulates
                    Gadget-[23] and GIZMO's default behaviours.
  -e                Enable floating-point exceptions (debugging mode).
  -f          {int} Overwrite the CPU frequency (Hz) to be used for time measurements.
  -g                Run with an external gravitational potential.
  -G                Run with self-gravity.
  -M                Reconstruct the multipoles every time-step.
  -n          {int} Execute a fixed number of time steps. When unset use the time_end parameter to stop.
  -P  {sec:par:val} Set parameter value and overwrites values read from the parameters file. Can be used more than once.
  -s                Run with hydrodynamics.
  -S                Run with stars.
  -t          {int} The number of threads to use on each MPI rank. Defaults to 1 if not specified.
  -T                Print timers every time-step.
  -v           [12] Increase the level of verbosity:
                    1: MPI-rank 0 writes,
                    2: All MPI-ranks write.
  -y          {int} Time-step frequency at which task graphs are dumped.
  -Y          {int} Time-step frequency at which threadpool tasks are dumped.
  -h                Print this help message and exit.

See the file examples/parameter_example.yml for an example of parameter file.
```
