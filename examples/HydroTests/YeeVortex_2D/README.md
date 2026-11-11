Yee isentropic vortex (2D)
==========================

The stationary isentropic vortex of Yee, Vinokur & Djomehri (2000) in the form
used by Pakmor et al., "Improving the convergence properties of the
moving-mesh code AREPO", MNRAS 455, 1134 (2016), arXiv:1503.00562, appendix A:

    T(r)   = T_inf - (gamma - 1) beta^2 / (8 gamma pi^2) exp(1 - r^2)
    rho(r) = T^(1 / (gamma - 1)),   u = T / (gamma - 1)
    v      = beta / (2 pi) exp((1 - r^2) / 2) (-y, x)

with gamma = 7/5, beta = 5 and T_inf = 1 on a periodic box of side 10 (the
paper uses [-5, 5]^2). The density is 0.49 at the centre and 1 far out, the
azimuthal velocity peaks at 0.80 at r = 1 (Mach 0.67) and the solution is
stationary, so the error at any time is the departure from the initial state.
The run goes to t = 10, about 1.3 rotations at r = 1. The paper uses the test
to measure the order of convergence of a scheme: it quotes the L1 norm of the
density error, (1/V) sum_i |rho_i - rho(x_i)| V_i, at t = 10 for resolutions
of 40^2 to 2560^2 (its figure 4, digitised in `pakmor2016_fig4_*.txt`: the
original MUSCL-Hancock + Green-Gauss Arepo is first order, the Runge-Kutta +
least-squares version second order).

Configuration
-------------

    ./configure --with-hydro-dimension=2 --with-adiabatic-index=7/5 --with-hydro=magma2 --with-kernel=quintic-spline
    ./configure --with-hydro-dimension=2 --with-adiabatic-index=7/5 --with-hydro=sphenix [--with-kernel=quintic-spline]

The parameter file carries SWIFT's generic eta 1.2348 and CFL 0.1; MAGMA's
recommended values are eta 1.5 and CFL 0.2, which the scripts below pass to
swift when SWIFT_ARGS is set, e.g.

    SWIFT_ARGS="-P SPH:resolution_eta:1.5 -P SPH:CFL_condition:0.2" ./run_convergence.sh

Initial conditions
------------------

`makeIC.py [N]` puts N^2 particles on a 2D glass scaled to the box (the
glassPlane_128 glass, tiled for N > 128; glassPlane_32 and 64 for the lower
resolutions) with particle masses rho(x) L^2 / N^2, so the spacing is uniform
and the density is carried by the masses, and sets u and v from the analytic
fields. `getGlass.sh` fetches the three glass files.

Running
-------

`run.sh` runs 128^2 particles on 4 threads (about a minute) and plots the
density, azimuthal velocity, internal energy and pressure against r at t = 10
with the L1 errors (`plotSolution.py <snapshot>`). `run_convergence.sh [N ...]`
runs 32^2 to 256^2 in `run_N<N>/` (256^2 takes about ten minutes on 4
threads) and calls `plotConvergence.py`, which plots the L1 errors against N
with the paper's Arepo curves and the kinetic-energy history of every run.
Three error measures are printed: the L1 of the particle densities, the L1 of
the azimuthally averaged density profile (the particle noise of the SPH density
estimate on a glass is a floor of about 2e-3 in L1(rho) that does not shrink
with N, and the profile removes it), and the L1 of the velocity.

What to expect
--------------

SPH schemes converge at first order on this smooth problem, and the error is
viscous: the total energy is conserved to 1e-5 but the artificial viscosity
thermalises part of the vortex's kinetic energy at a rate that falls with
resolution. With MAGMA (quintic spline, eta = 1.5, CFL 0.2, i.e. the
SWIFT_ARGS line above) the slope of L1(rho) is -0.76 (-0.91 for the profile) with errors two to three times those
of the paper's first-order Arepo at equal N, and 1.2 % of the kinetic energy
is lost by t = 10 at 512^2. SPHENIX with the Wendland C2 kernel at the
parameter file's eta = 1.2348 does not converge in L1(rho): that is only 17 neighbours in 2D,
where the Wendland density estimate is biased high by 1.9 % everywhere (the
self-contribution bias of Dehnen & Aly 2012), a floor of 1.8e-2 in L1(rho);
its velocity error converges at slope -0.4 and 13 % of the kinetic energy is
lost by t = 10 at 512^2. With the quintic spline at eta = 1.5 SPHENIX
converges at first order with errors 1.5 to 2.3 times MAGMA's. For a
like-for-like comparison of 2D schemes use the same kernel and eta.
