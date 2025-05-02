A demo planetary simulation of a giant impact, using the initial conditions
created by the DemoImpactInitCond example with SEAGen and WoMa [1,2].

The scenario is a so-called canonical-like Moon-forming impact of a ~Mars-mass
body onto the proto-Earth, at the mutual escape speed, with an impact angle of
45 degrees, set up as described in [3].

The initial conditions can be created via the DemoImpactInitCond example, see
its README.md and the WoMa documentation for more details. Or the premade
initial conditions can be downloaded directly with `get_init_cond.sh`.

The default resolution is ~10^5 particles, set by the `N` variables at the start
of the scripts. Note that here the simulation is cut short with an end time of
just 15 h, for testing. Note also that this type of graze-and-merge collision is
notoriously chaotic, so minor changes to the initial conditions can have large
effects on the outcome, and a higher resolution than the default for this demo
is required for convergence [3].

If you have just run the initial settling simulations, then please ensure you
have configured and compiled SWIFT without the fixed-entropy setting before
running an impact simulation. You may find it convenient to rename your SWIFT
executables according to their configuration.

[1] Kegerreis et al. 2019, MNRAS 487:4
[2] Ruiz-Bonilla et al. 2021, MNRAS 500:3
[3] Kegerreis et al. 2022, ApJL 937:2 L40


This example requires the code to be configured to use the planetary
hydrodynamics scheme and planetary equations of state:
    --with-hydro=planetary --with-equation-of-state=planetary
