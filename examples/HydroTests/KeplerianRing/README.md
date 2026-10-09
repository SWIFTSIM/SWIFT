Keplerian disc test (Hopkins 2015)
==================================

A cold, pressure-less gas disc on circular orbits around a point mass. The
disc should stay in its initial configuration indefinitely; numerical
(artificial) shear viscosity transports angular momentum, spreads the disc and
eventually breaks it up through the viscous instability. The set-up follows
Section 4.2.4 of Hopkins (2015, MNRAS 450, 53, "A new class of accurate,
mesh-free hydrodynamic simulation methods"), which itself is a variant of the
ring tests of Maddison et al. (1996), Imaeda & Inutsuka (2002), Cartwright et
al. (2009) and Cullen & Dehnen (2010). For the Cullen & Dehnen (2010) version
(a Gaussian ring with their units and particle placement) see the example
`KeplerianRing_CullenDehnen`.

Definitions (Hopkins 2015, eq. 34 and text):

+ surface density Sigma = 1 for 0.5 <= r <= 2 (the paper's profile has a
  (r/0.5)^3 rise inside and a (1 + (r - 2)/0.1)^-3 fall outside plus a 0.01
  floor; for particle codes Hopkins removes the particles inside r < 0.5 and
  outside r > 2, which is what `makeIC.py` does by default; `--background 1`
  keeps the outer floor as a 100:1 lower-density background filling the box);
+ vanishing, uniform pressure P = 1e-6, i.e. u = 1.5e-6 for gamma = 5/3
  (Mach ~ 700 at r = 1);
+ softened Keplerian potential Phi = -GM / sqrt(r^2 + eps^2) with GM = 1; the
  paper does not quote eps, this example uses eps = 0.05 (the IC velocities
  are the circular velocities of the softened potential);
+ particles on a uniform grid of spacing 4/256 across the square -2 < x < 2
  ("effective resolution 256^2", about 48000 particles in the disc), equal
  masses; `-n` sets this effective resolution;
+ evolved to t = 120, about 20 orbits at r = 1 (orbital period 2 pi r^1.5).

Hopkins measures the surface-density map at t = 120 and, for convergence, the
L1 velocity error at r = 1 after 10 orbits and the times at which the L1
velocity error at r = 1 exceeds 0.01 and the L1 density error exceeds 0.25.

Code setup
----------

The disc can be run with the 2D code (particles in the z = 0 plane, the point
mass at z = 0) or as a razor-thin disc with the 3D code (all particles in the
z = 5 plane). In both cases SWIFT has to be configured with the softened
point-mass potential:

    ./configure --with-ext-potential=point-mass-softened              # 3D planar
    ./configure --with-hydro-dimension=2 --with-ext-potential=point-mass-softened

The parameter file caps the smoothing length (`SPH:h_max`) because particles
ejected from the disc can lose all their neighbours, and uses a small minimal
time step: particles that fall onto the point mass orbit at the softening
radius with very short time steps, and the run must not abort on them (the
number of such particles is one measure of the viscosity of a scheme).

Initial conditions
------------------

    python3 makeIC.py --help

`-n` effective resolution (default 256), `-i` internal energy (default
1.5e-6), `-e` softening (default 0.05, also the inner exclusion radius unless
`-c` is larger), `-c` inner exclusion radius (default 0.5), `--background`
(default 0), `-m` GM (default 1; change `PointMassPotential:mass` in the
parameter file accordingly), `-g` generation method (`grid` for this test;
`spiral` and `gaussian` place the particles with the method of Cartwright et
al. 2009 and are kept for experiments).
