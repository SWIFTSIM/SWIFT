A demo of making initial conditions for a planetary impact simulation using
SEAGen and WoMa [1,2], and running standard "settling" simulations. See WoMa's
documentation and tutorial at https://github.com/srbonilla/WoMa for more info.

First, the target and impactor planet profiles and particles are generated using
WoMa. Second, a "settling" simulation is run for each body in isolation to allow
any final relaxation of the particles to occur. Third, the settled particles are
loaded and given initial positions and velocities to set up the impact scenario.

The resulting initial conditions are used by the DemoImpact example. See its
README.md for more info, and to run the impact simulation.

The set up for these planets and the collision scenario is described in [3 (see
S2 and Fig. A1, etc)]. These bodies use ANEOS equations of state that include
entropy information, so for settling simulations the entropy of each particle
can be kept fixed to ensure adiabatic settling without viscosity heating etc.

The resolution, set by the number of ~equal-mass particles, can be controlled by
the `N` variables at the start of the scripts, with examples included for 10^5,
10^6 (default), and 10^7 particles. In the .yml input files for different
resolutions, in addition to the different file names, the gravitational
softening set by max_physical_baryon_softening is set to approximately the
minimum inter-particle separation. See the planetary and other sections of the
SWIFT documentation for details on the other .yml input parameters and other
simulation options.

[1] Kegerreis et al. 2019, MNRAS 487:4
[2] Ruiz-Bonilla et al. 2021, MNRAS 500:3
[3] Kegerreis et al. 2022, ApJL 937:2 L40


This example requires the code to be configured to use the planetary
hydrodynamics scheme and planetary equations of state, and fixed entropy:
    --with-hydro=planetary --with-equation-of-state=planetary --enable-planetary-fixed-entropy
