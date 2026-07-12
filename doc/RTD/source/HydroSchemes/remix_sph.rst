.. REMIX SPH
    Thomas Sandnes, 13th May 2025

.. _remix_sph:

REMIX SPH
==============================================

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Contents:

REMIX is an SPH scheme designed to alleviate effects that typically suppress
mixing and instability growth at density discontinuities in SPH simulations
(Sandnes et al. 2025). REMIX addresses this problem by directly targeting sources
of kernel smoothing error and discretisation error, resulting in a generalised,
material-independent formulation that improves the treatment both of
discontinuities within a single material, for example in an ideal gas, and of
interfaces between dissimilar materials. The scheme combines:

+ An evolved density estimate to avoid the kernel smoothing error in the
  standard SPH integral density estimate;
+ Thermodynamically consistent, conservative equations of motion, with
  free functions chosen to limit zeroth-order error;
+ Linear-order reproducing kernels with grad-h terms and a vacuum interface
  treatment;
+ A "kernel normalising term" to avoid potential accumulation of error in
  the evolved density estimate, such that densities are ensured to remain 
  representative of the distribution of particle masses in the simulation volume;
+ Advanced artificial viscosity and diffusion schemes with linear reconstruction
  of quantities to particle midpoints, and a set of novel improvements to
  effectively switch between treatments for shock-capturing under compression and
  noise-smoothing in shearing regions.

To configure with this scheme, use

.. code-block:: bash

  ./configure --with-hydro=remix --with-equation-of-state=planetary


This scheme allows multiple materials,
meaning that different SPH particles can be assigned different
`equations of state <equations_of_state.html>`_ (EoS).
Every SPH particle then requires and carries the additional ``MaterialID`` flag
from the initial conditions file. This flag indicates the particle's material
and which EoS it should use. Note that configuring with
``--with-equation-of-state=planetary`` is required for this scheme, although
for simulations that use a single, ideal gas EoS, setting all MaterialIDs to
``0`` and including

.. code-block:: yaml

  EoS:
      planetary_use_idg_def:    1

in the parameter file are the only EoS-related additions needed compared with
other non-Planetary hydro schemes. Note also that since densities are evolved in
time, initial particle densities are required in initial conditions.

We additionally recommend configuring with ``--with-kernel=wendland-C2`` and with

.. code-block:: yaml

  SPH:
    resolution_eta:        1.487

in the parameter file for improved hydrodynamic behaviour and since this is the
configuration used for the validation simulations of Sandnes et al. (2025).


The current implementation of the REMIX hydro scheme has been validated for
planetary applications and various hydrodynamic test cases, and does not include 
all necessary functionality for e.g. cosmological simulations.

Default parameters used in the artificial viscosity and diffusion schemes and the
normalising term (see Sandnes et al. 2025) are:

.. code-block:: c

  #define const_remix_visc_alpha 1.5f
  #define const_remix_visc_beta 3.f
  #define const_remix_visc_epsilon 0.1f
  #define const_remix_visc_a 2.0f / 3.0f
  #define const_remix_visc_b 1.0f / 3.0f
  #define const_remix_difn_a_u 0.05f
  #define const_remix_difn_b_u 0.95f
  #define const_remix_difn_a_rho 0.05f
  #define const_remix_difn_b_rho 0.95f
  #define const_remix_norm_alpha 1.0f
  #define const_remix_slope_limiter_exp_denom 0.04f

These can be changed in ``src/hydro/REMIX/hydro_parameters.h``.
