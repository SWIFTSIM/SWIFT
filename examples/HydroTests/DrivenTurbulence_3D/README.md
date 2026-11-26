Driven subsonic turbulence (Bauer & Springel 2012)
==================================================

Reproduction of the Mach ~ 0.3 driven-turbulence test of

  A. Bauer & V. Springel, "Subsonic turbulence in smoothed particle hydrodynamics
  and moving-mesh simulations", MNRAS 423, 2558 (2012), arXiv:1109.4413,
  https://doi.org/10.1111/j.1365-2966.2012.21058.x

Set-up (section 3 and Table 4 of the paper): a unit periodic box of gas with
rho = 1 and isothermal sound speed c_s = 1, initially at rest on a Cartesian
lattice, stirred by a solenoidal Ornstein-Uhlenbeck force field with
sigma = 0.014, update interval 0.005, correlation time t_s = 1, driven modes
with k in [6.27, 12.57] (the two largest scales of the box) and a k^-5/3
amplitude spectrum. A quasi-stationary state with an rms Mach number of about
0.3 is reached after t ~ 10 and the run ends at t = 25.6. The paper measures
the velocity power spectrum E(k) (normalised so that its integral is the
velocity variance) by nearest-particle sampling of the velocity on a grid
twice finer than the particle resolution, averaged over t = 10 to 25.6, and
plots k E(k) against k (their figures 5 and 6). The dissipation range of SPH
at 64^3 to 512^3, and the small-scale velocity noise that makes the SPH
spectrum rise again towards the Nyquist wavenumber, are the features of
interest.

Configuration
-------------

The gas must be nearly isothermal. Configure SWIFT with

    ./configure --with-forcing=driven-turbulence --with-adiabatic-index=1.001 [--with-hydro=...]

(gamma = 1.001 gives P = rho c_s^2 to 0.1 % with u = 1000 for c_s = 1).
GSL is required by the forcing. The run script uses 64^3 particles; change the
argument of makeIC.py for other resolutions.

Driving normalisation
---------------------

The forcing module is a port of Arepo's TURB_DRIVING (also used by Pakmor et
al. 2026, arXiv:2607.06660). It applies the force density
2 * amplitude_factor * sqrt(3/2) * (sum over the driven modes) in 3D for
purely solenoidal driving. With amplitude_factor = 1 the box reaches a
stationary rms Mach number of about 0.35 for SPH, above the 0.33 (Arepo) and
0.30 (SPH with time-variable viscosity) of the paper. With
amplitude_factor = 1/sqrt(3/2) = 0.8165, i.e. without the solenoidal-weight
normalisation, which the 2012 code did not have, SPHENIX at 32^3 and 64^3
gives 0.30 and MAGMA 0.33, matching the paper for the respective class of
scheme. The parameter file therefore uses 0.8165. The Mach number is
insensitive to the resolution (as the paper also notes), so a 32^3 run of a
few minutes is enough to check the calibration of a new scheme.

The time-step is limited by the forcing to the update interval of the field.

Outputs
-------

`plotSolution.py` computes the power spectrum of every snapshot with t >= 10
(`power_spectrum.py`, cached in `spectra/`), averages them, and plots k E(k)
and the Mach-number history against the curves of the paper, which were
digitised from its figures and stored next to the scripts:

  - `bs12_fig6top_A1.txt`, `bs12_fig6bot_S1.txt`, `bs12_fig6bot_S1-tav.txt`:
    64^3 Arepo, SPH (alpha = 1) and SPH with time-variable viscosity (fig. 6),
  - `bs12_fig5_A3.txt`, `bs12_fig5_S3.txt`, `bs12_fig5_S3-tav.txt`: the same
    at 256^3 (fig. 5),
  - `bs12_fig2_*.txt`: the Mach-number histories of the 256^3 runs (fig. 2).

Columns are k and k E(k) (t and Mach number for fig. 2).

Cost
----

With SPHENIX on 4 threads, 64^3 to t = 25.6 takes about 2.5 hours (16 000
steps, the CFL condition at 0.1 sets dt = 1.6e-3); 32^3 takes 10 minutes.
