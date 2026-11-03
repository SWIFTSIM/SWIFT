Keplerian ring test (Cullen & Dehnen 2010)
==========================================

The 2D Keplerian ring of Cullen & Dehnen (2010, MNRAS 408, 669,
http://arxiv.org/abs/1006.1524), Section 4.3 and Figure 8. A cold gas ring
with a Gaussian surface density orbits a point mass; pressure forces are
negligible (sound speed 1e-3 of the orbital velocity) and the ring should stay
in its initial configuration. Artificial viscosity acting on the Keplerian
shear transports angular momentum, spreads the ring and makes it break up
through the viscous instability (Lyubarskij et al. 1994) within a few orbits.
Cullen & Dehnen use it to show that their viscosity switch with its shear
limiter keeps the ring intact for more than five orbits where standard SPH,
the Morris & Monaghan switch and the Balsara switch all fail.

Definitions (Cullen & Dehnen 2010, Section 4.3):

+ GM = 1000 for the central point mass; the orbital period at the ring
  centre is 2 pi and the orbital velocity 10;
+ Gaussian surface density centred on r = 10 with width (standard deviation)
  2.5;
+ N = 9745 equal-mass particles placed with the method of Cartwright,
  Stamatellos & Whitworth (2009);
+ sound speed c = 0.01 everywhere (Mach 1000 at the ring centre);
+ evolved for more than five orbits; their Figure 8 shows the particle
  positions at t = 10, 15 and a later time for each scheme.

Deviations in this example: SWIFT's point-mass potential is softened, with a
softening of 0.1 (1% of the ring radius) and circular velocities computed for
the softened potential; the ring sits in a 50 x 50 periodic box; the kernel and
neighbour number are SWIFT's defaults rather than the paper's cubic spline
with 13 neighbours; and the Gaussian is truncated at 2 widths (r in [5, 15],
95% of the mass, option `--truncate`). The paper's untruncated tails place a
few dozen isolated particles where the surface density is a per cent of the
peak, many spacings apart: no kernel resolves them, they orbit much faster
than the ring, fall onto the point mass within an orbit and then act on the
ring's inner edge through kernels the size of the gap. Cullen & Dehnen report
the resulting central features as artefacts; with SWIFT's kernels they break
the ring for every scheme, so the resolved annulus is the form of the test
that probes the viscosity rather than the tails. The example `KeplerianRing` is the related disc test of
Hopkins (2015).

Code setup
----------

The test runs with the 2D code and the softened point-mass potential:

    ./configure --with-hydro-dimension=2 --with-ext-potential=point-mass-softened

The parameter file uses a small minimal time step: particles that fall onto
the point mass orbit at the softening radius with very short time steps, and
the run must not abort on them (how many particles do so, and how early, is a
measure of the viscosity of a scheme). The smoothing length is capped
(`SPH:h_max`) for particles ejected from the ring.

Initial conditions
------------------

    python3 makeIC.py --help

`-n` number of particles (default 9745), `-c` sound speed (default 0.01),
`-e` softening (default 0.1; keep `PointMassPotential:softening` consistent),
`--truncate` Gaussian truncation in widths (default 2; 0 for the paper's untruncated ring),
`-m` GM (default 1000; keep `PointMassPotential:mass` consistent, G = 6.6743e-8
in the cgs internal units of the parameter file).
