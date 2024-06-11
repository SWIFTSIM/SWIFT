.. Parameter Description
   Matthieu Schaller, 21st October 2018

.. _Parameters_basics:

File format and basic information
---------------------------------

The parameter file uses a format similar to the `YAML format
<https://en.wikipedia.org/wiki/YAML>`_ but reduced to only the
elements required for the SWIFT parameters. Options are given by a
name followed by a column and the value of the parameter:

.. code:: YAML

   ICs:        santa_barbara.hdf5
   dt_max:     1.5
   shift:      [2., 4., 5.]

Comments can be inserted anywhere and start with a hash:

.. code:: YAML

   # Description of the physics
   viscosity_alpha:     2.0
   dt_max:              1.5     # seconds

A typical SWIFT parameter file is split into multiple sections that
may or may not be present depending on the different configuration
options. The sections start with a label and can contain any number of
parameters:

.. code:: YAML

   Cosmology:    # Planck13
     Omega_cdm:      0.2587481
     Omega_lambda:   0.693
     Omega_b:        0.0482519
     h:              0.6777
     a_begin:        0.0078125     # z = 127

The options can be integer values, floating point numbers, characters
or strings. If SWIFT expects a number and string is given, an error
will be raised. The code can also read an array of values:

.. code:: YAML

   shift:  [2., 4., 5.]

Some options in the parameter file are optional and
when not provided, SWIFT will run with the default value. However, if
a compulsory parameter is missing an error will be raised at
start-up.

Finally, SWIFT outputs two YAML files at the start of a run. The first one
``used_parameters.yml`` contains all the parameters that were used for this
run, **including all the optional parameters left unspecified with their
default values**. This file can be used to start an exact copy of the run. The
second file, ``unused_parameters.yml`` contains all the values that were not
read from the parameter file. This can be used to simplify the parameter file
or check that nothing important was ignored (for instance because the code is
not configured to use some options). Note that on restart a new file
``used_parameters.yml.stepno`` is created and any changed parameters will be
written to it.

The rest of this page describes all the SWIFT parameters, split by
section. A list of all the possible parameters is kept in the file
``examples/parameter_examples.yml``.

.. _Parameters_meta_data:

Meta Data
---------

The ``MetaData`` section contains basic information about the simulation. It
currently only contains one parameter: ``run_name``. This is a string of
characters describing the simulation. It is written to the snapshots' headers.

.. _Parameters_units:

Internal Unit System
--------------------

The ``InternalUnitSystem`` section describes the units used internally by the
code. This is the system of units in which all the equations are solved. All
physical constants are converted to this system and if the ICs use a different
system (see the snapshots' :ref:`ICs_units_label` section of the documentation)
the particle quantities will be converted when read in.

The system of units is described using the value of the 5 basic units
of any system with respect to the CGS system. Instead of using a unit
of time we use a unit of velocity as this is more intuitive. Users
hence need to provide:

* a unit of length: ``UnitLength_in_cgs``,
* a unit of mass: ``UnitMass_in_cgs``,
* a unit of velocity ``UnitVelocity_in_cgs``,
* a unit of electric current ``UnitCurrent_in_cgs``,
* a unit of temperature ``UnitTemp_in_cgs``.

All these need to be expressed with respect to their cgs counter-part
(i.e. :math:`cm`, :math:`g`, :math:`cm/s`, :math:`A` and :math:`K`). Recall
that there are no h-factors in any of SWIFT's quantities; we, for instance,
use :math:`cm` and not :math:`cm/h`.

For instance to use the commonly adopted system of 10^10 Msun as a
unit for mass, mega-parsec as a unit of length and km/s as a unit of
speed, we would use:

.. code:: YAML

   # Common unit system for cosmo sims
   InternalUnitSystem:
     UnitMass_in_cgs:     1.98848e43    # 10^10 M_sun in grams
     UnitLength_in_cgs:   3.08567758e24 # 1 Mpc in centimeters
     UnitVelocity_in_cgs: 1e5           # 1 km/s in centimeters per second
     UnitCurrent_in_cgs:  1             # 1 Ampere
     UnitTemp_in_cgs:     1             # 1 Kelvin

Note that there are currently no variables in any of the SWIFT physics
schemes that make use of the unit of electric current. There is also
no incentive to use anything else than Kelvin but that makes the whole
system consistent with any possible unit system.

If one is interested in using the more humorous `FFF unit
system <https://en.wikipedia.org/wiki/FFF_system>`_ one would use

.. code:: YAML

   # FFF unit system
   InternalUnitSystem:
     UnitMass_in_cgs:     40823.3133  # 1 Firkin (fir) in grams
     UnitLength_in_cgs:   20116.8     # 1 Furlong (fur) in cm
     UnitVelocity_in_cgs: 0.01663095  # 1 Furlong (fur) per Fortnight (ftn) in cm/s
     UnitCurrent_in_cgs:  1           # 1 Ampere
     UnitTemp_in_cgs:     1           # 1 Kelvin

The value of the physical constants in this system is left as an
exercise for the reader [#f1]_.

.. _Parameters_cosmology:

Cosmology
---------

When running a cosmological simulation, the section ``Cosmology`` sets the values of the
cosmological model. The expanded :math:`\Lambda\rm{CDM}` parameters governing the
background evolution of the Universe need to be specified here. These are:

* The reduced Hubble constant: :math:`h`: ``h``,
* The cold dark matter density parameter :math:`\Omega_cdm`: ``Omega_cdm``,
* The cosmological constant density parameter :math:`\Omega_\Lambda`: ``Omega_lambda``,
* The baryon density parameter :math:`\Omega_b`: ``Omega_b``,
* The radiation density parameter :math:`\Omega_r`: ``Omega_r``.

The last parameter can be omitted and will default to :math:`\Omega_r = 0`. Note
that SWIFT will verify on start-up that the matter content of the initial conditions
matches the cosmology specified in this section.

This section also specifies the start and end of the simulation expressed in
terms of scale-factors. The two parameters are:

* Initial scale-factor: ``a_begin``,
* Final scale-factor: ``a_end``.

Two additional optional parameters can be used to change the equation of
state of dark energy :math:`w(a)`. We use the evolution law :math:`w(a) =
w_0 + w_a (1 - a)`. The two parameters in the YAML file are:

* The :math:`z=0` dark energy equation of state parameter :math:`w_0`: ``w_0``
* The dark energy equation of state evolution parameter :math:`w_a`: ``w_a``

If unspecified these parameters default to the default
:math:`\Lambda\rm{CDM}` values of :math:`w_0 = -1` and :math:`w_a = 0`.

The radiation density :math:`\Omega_r` can also be specified by setting
an alternative optional parameter:

* The number of ultra-relativistic degrees of freedom :math:`N_{\rm{ur}}`:
  ``N_ur``.

The radiation density :math:`\Omega_r` is then automatically inferred from
:math:`N_{\rm{ur}}` and the present-day CMB temperature
:math:`T_{\rm{CMB},0}=2.7255` Kelvin. This parametrization cannot
be used together with :math:`\Omega_r`. If neither parameter is used, SWIFT
defaults to :math:`\Omega_r = 0`. Note that :math:`N_{\rm{ur}}` differs from
:math:`N_{\rm{eff}}`, the latter of which also includes massive neutrinos.

Massive neutrinos can be included by specifying the optional parameters:

* The number of massive neutrino species :math:`N_{\nu}`: ``N_nu``,
* A comma-separated list of neutrino masses in eV: ``M_nu_eV``,
* A comma-separated list of neutrino degeneracies: ``deg_nu``,
* The present-day neutrino temperature :math:`T_{\nu,0}`: ``T_nu_0``.

When including massive neutrinos, only ``N_nu`` and ``M_nu_eV`` are necessary.
By default, SWIFT will assume non-degenerate species and
:math:`T_{\nu,0}=(4/11)^{1/3}T_{\rm{CMB},0}`. Neutrinos do not contribute to
:math:`\Omega_m = \Omega_{\rm{cdm}} + \Omega_b` in our conventions.

For a Planck+13 cosmological model (ignoring radiation density as is
commonly done) and running from :math:`z=127` to :math:`z=0`, one would hence
use the following parameters:

.. code:: YAML

   Cosmology:        # Planck13 (EAGLE flavour)
     a_begin:        0.0078125     # z = 127
     a_end:          1.0           # z = 0
     h:              0.6777
     Omega_cdm:      0.2587481
     Omega_lambda:   0.693
     Omega_b:        0.0482519
     Omega_r:        0.            # (Optional)
     w_0:            -1.0          # (Optional)
     w_a:            0.            # (Optional)

When running a non-cosmological simulation (i.e. without the ``--cosmology`` run-time
flag) this section of the YAML file is entirely ignored.

.. _Parameters_gravity:

Gravity
-------

The behaviour of the self-gravity solver can be modified by the parameters
provided in the ``Gravity`` section. The theory document puts these parameters into the
context of the equations being solved. We give a brief overview here.

* The Plummer-equivalent co-moving softening length used for all dark matter particles :math:`\epsilon_{\rm com,DM}`: ``comoving_DM_softening``,
* The Plummer-equivalent co-moving softening length used for all baryon particles (gas, stars, BHs) :math:`\epsilon_{\rm com,bar}`: ``comoving_baryon_softening``,
* The Plummer-equivalent maximal physical softening length used for all dark matter particles :math:`\epsilon_{\rm max,DM}`: ``max_physical_DM_softening``,
* The Plummer-equivalent maximal physical softening length used for all baryon particles (gas, stars, BHs) :math:`\epsilon_{\rm max,bar}`: ``max_physical_baryon_softening``,

At any redshift :math:`z`, the Plummer-equivalent softening length used by
the code will be :math:`\epsilon=\min(\epsilon_{max},
\frac{\epsilon_{com}}{z+1})`. The same calculation is performed
independently for the dark matter and baryon particles. All the softening
quantities are expressed in internal units. Calculations that only involve
DM or baryons can leave the unused quantities out of the parameter
file. For non-cosmological runs, only the physical softening lengths need
to be supplied.

In case of zoom simulations, the softening of the additional, more massive, background
particles is specified via the parameter
``softening_ratio_background``. Since these particles will typically have
different masses to degrade the resolution away from the zoom-in region, the
particles won't have a single softening value. Instead, we specify the
fraction of the mean inter-particle separation to use. The code will then
derive the softening length of each particle assuming the mean density of
the Universe. That is :math:`\epsilon_{\rm background} =
f\sqrt[3]{\frac{m}{\Omega_m\rho_{\rm crit}}}`, where :math:`f` is the
user-defined value (typically of order 0.05).

The accuracy of the gravity calculation is governed by the following four parameters:

* The multipole acceptance criterion: ``MAC``
* The fixed opening angle used in the geometric MAC :math:`\theta_{\rm cr}`: ``theta_cr``,
* The accuracy criterion used in the adaptive MAC:  :math:`\epsilon_{\rm fmm}`: ``epsilon_fmm``,
* The time-step size pre-factor :math:`\eta`: ``eta``,

The first three parameters govern the way the Fast-Multipole method tree-walk is
done (see the theory documents for full details).  The ``MAC`` parameter can
take three values: ``adaptive``, ``geometric``, or ``gadget``. In the first
case, the tree recursion decision is based on the estimated accelerations that a
given tree node will produce, trying to recurse to levels where the fractional
contribution of the accelerations to the cell is less than :math:`\epsilon_{\rm
fmm}`. In the second case, a fixed Barnes-Hut-like opening angle
:math:`\theta_{\rm cr}` is used. The final case corresponds to the choice made
in the Gadget-4 code. It is an implementation using eq. 36 of `Springel et
al. (2021) <https://adsabs.harvard.edu/abs/2021MNRAS.506.2871S>`_.

The time-step of a given particle is given by :math:`\Delta t =
\sqrt{2\eta\epsilon_i/|\overrightarrow{a}_i|}`, where
:math:`\overrightarrow{a}_i` is the particle's acceleration and
:math:`\epsilon_i` its (spline) softening length. `Power et al. (2003)
<http://adsabs.harvard.edu/abs/2003MNRAS.338...14P>`_ recommend using
:math:`\eta=0.025`.

Two further parameters determine when the gravity tree is reconstructed:

* The tree rebuild frequency: ``rebuild_frequency``.
* The fraction of active particles to trigger a rebuild:
  ``rebuild_active_fraction``.

The tree rebuild frequency is an optional parameter defaulting to
:math:`0.01`. It is used to trigger the re-construction of the tree every
time a fraction of the particles have been integrated (kicked) forward in
time. The second parameter is also optional and determines a separate rebuild
criterion, based on the fraction of particles that is active at the
beginning of a step. This can be seen as a forward-looking version of the
first criterion, which can be useful for runs with very fast particles.
The second criterion is not used for values :math:`>1`, which is the default
assumption.


The last tree-related parameters are:

* Whether or not to use the approximate gravity from the FMM tree below the
  softening scale: ``use_tree_below_softening`` (default: 0)
* Whether or not the truncated force estimator in the adaptive tree-walk
  considers the exponential mesh-related cut-off:
  ``allow_truncation_in_MAC`` (default: 0)

These parameters default to good all-around choices. See the
theory documentation about their exact effects.

Simulations using periodic boundary conditions use additional parameters for the
Particle-Mesh part of the calculation. The last five are optional:

* The number cells along each axis of the mesh :math:`N`: ``mesh_side_length``,
* Whether or not to use a distributed mesh when running over MPI: ``distributed_mesh`` (default: ``0``),
* Whether or not to use local patches instead of direct atomic operations to
  write to the mesh in the non-MPI case (this is a performance tuning
  parameter): ``mesh_uses_local_patches`` (default: ``1``),
* The mesh smoothing scale in units of the mesh cell-size :math:`a_{\rm
  smooth}`: ``a_smooth`` (default: ``1.25``),
* The scale above which the short-range forces are assumed to be 0 (in units of
  the mesh cell-size multiplied by :math:`a_{\rm smooth}`) :math:`r_{\rm
  cut,max}`: ``r_cut_max`` (default: ``4.5``),
* The scale below which the short-range forces are assumed to be exactly Newtonian (in units of
  the mesh cell-size multiplied by :math:`a_{\rm smooth}`) :math:`r_{\rm
  cut,min}`: ``r_cut_min`` (default: ``0.1``),

For most runs, the default values can be used. Only the number of cells along
each axis needs to be specified. The remaining three values are best described
in the context of the full set of equations in the theory documents.

By default, SWIFT will replicate the mesh on each MPI rank. This means that a
single MPI reduction is used to ensure all ranks have a full copy of the density
field. Each node then solves for the potential in Fourier space independently of
the others. This is a fast option for small meshes. This technique is limited to
mesh with sizes :math:`N<1291` due to the limitations of MPI. Larger meshes need
to use the distributed version of the algorithm. The code then also needs to be
compiled with ``--enable-mpi-mesh-gravity``. That algorithm is slower for small
meshes but has no limits on the size of the mesh and truly huge Fourier
transforms can be performed without any problems. The only limitation is the
amount of memory on each node. The algorithm will use ``N^3 * 8 * 2 / M`` bytes
on each of the ``M`` MPI ranks.

As a summary, here are the values used for the EAGLE :math:`100^3~{\rm Mpc}^3`
simulation:

.. code:: YAML

   # Parameters for the self-gravity scheme for the EAGLE-100 box
   Gravity:
     eta:                    0.025
     MAC:                    adaptive
     theta_cr:               0.6
     epsilon_fmm:            0.001
     mesh_side_length:       2048
     distributed_mesh:       0
     comoving_DM_softening:         0.0026994  # 0.7 proper kpc at z=2.8.
     max_physical_DM_softening:     0.0007     # 0.7 proper kpc
     comoving_baryon_softening:     0.0026994  # 0.7 proper kpc at z=2.8.
     max_physical_baryon_softening: 0.0007     # 0.7 proper kpc
     rebuild_frequency:      0.01   # Default optional value
     a_smooth:          1.25        # Default optional value
     r_cut_max:         4.5         # Default optional value
     r_cut_min:         0.1         # Default optional value
     use_tree_below_softening: 0    # Default optional value
     allow_truncation_in_MAC:  0    # Default optional value

.. _Parameters_SPH:

SPH
---

The ``SPH`` section is used to set parameters that describe the SPH
calculations. There are some scheme-specific values that are detailed in the
:ref:`hydro` section. The common parameters are detailed below.

In all cases, users have to specify two values:

* The smoothing length in terms of mean inter-particle separation:
  ``resolution_eta``
* The CFL condition that enters the time-step calculation: ``CFL_condition``

These quantities are dimensionless. The first, ``resolution_eta``, specifies
how smooth the simulation should be, and is used here instead of the number
of neighbours to smooth over as this also takes into account non-uniform
particle distributions. A value of 1.2348 gives approximately 48 neighbours
in 3D with the cubic spline kernel. More information on the choices behind
these parameters can be found in
`Dehnen & Aly 2012 <https://ui.adsabs.harvard.edu/abs/2012MNRAS.425.1068D/>`_.


The second quantity, the CFL condition, specifies how accurate the time
integration should be and enters as a pre-factor into the hydrodynamics
time-step calculation. This factor should be strictly bounded by 0 and 1, and
typically takes a value of 0.1 for SPH calculations.

The next set of parameters deal with the calculation of the smoothing lengths
directly and are all optional:

* Whether to use or not the mass-weighted definition of the SPH number of
  neighbours: ``use_mass_weighted_num_ngb`` (Default: 0)
* The (relative) tolerance to converge smoothing lengths within:
  ``h_tolerance`` (Default: 1e-4)
* The maximal smoothing length in internal units: ``h_max`` (Default: FLT_MAX)
* The minimal allowed smoothing length in terms of the gravitational
  softening: ``h_min_ratio`` (Default: 0.0, i.e. no minimum)
* The maximal (relative) allowed change in volume over one time-step:
  ``max_volume_change`` (Default: 1.4)
* The maximal number of iterations allowed to converge the smoothing
  lengths: ``max_ghost_iterations`` (Default: 30)

These parameters all set the accuracy of the smoothing lengths in various
ways. The first one specified what definition of the local number density
of particles to use. By default, we use

.. math::
   n_i = \sum_j W(\|\mathbf{r}_i - \mathbf{r}_j\|, h_i)

but switching on the ``use_mass_weighted_num_ngb`` flag changes the
defintion to:

.. math::
   n_i = \frac{\rho_i}{m_i}

where the density has been computed in the traditional SPH way
(i.e. :math:`\rho_i = \sum_j m_j W(\|\mathbf{r}_i - \mathbf{r}_j\|,
h_i)`). Note that in the case where all the particles in the simulation
have the same mass, the two definitions lead to the same number density
value.

**We dot not recommend using this alternative neighbour number definition
in production runs.** It is mainly provided for backward compatibility with
earlier simulations.

The second one, the relative tolerance for the smoothing length, specifies
the convergence criteria for the smoothing length when using the
Newton-Raphson scheme. This works with the maximal number of iterations,
``max_ghost_iterations`` (so called because the smoothing length calculation
occurs in the ghost task), to ensure that the values of the smoothing lengths
are consistent with the local number density. We solve:

.. math::
   (\eta \gamma)^{n_D} = n_i

with :math:`\gamma` the ratio of smoothing length to kernel support (this
is fixed for a given kernel shape), :math:`n_D` the number of spatial
dimensions, :math:`\eta` the value of ``resolution_eta``, and :math:`n_i`
the local number density. We adapt the value of the smoothing length,
:math:`h`, to be consistent with the number density.

The maximal smoothing length, by default, is set to ``FLT_MAX``, and if set
prevents the smoothing length from going beyond ``h_max`` (in internal units)
during the run, irrespective of the above equation. The minimal smoothing
length is set in terms of the gravitational softening, ``h_min_ratio``, to
prevent the smoothing length from going below this value in dense
environments. This will lead to smoothing over more particles than specified
by :math:`\eta`.

The optional parameter ``particle_splitting`` (Default: 0) activates the
splitting of overly massive particles into 2. By switching this on, the code
will loop over all the particles at every tree rebuild and split the particles
with a mass above a fixed threshold into two copies that are slightly shifted
(by a randomly orientated vector of norm :math:`0.2h`). Their masses and other
relevant particle-carried quantities are then halved. The mass threshold for
splitting is set by the parameter ``particle_splitting_mass_threshold`` which is
specified using the internal unit system. The IDs of the newly created particles
can be either drawn randomly by setting the parameter ``generate_random_ids``
(Default: 0) to :math:`1`. When this is activated, there is no check that the
newly generated IDs do not clash with any other pre-existing particle. If this
option is set to :math:`0` (the default setting) then the new IDs are created in
increasing order from the maximal pre-existing value in the simulation, hence
preventing any clash.

The final set of parameters in this section determine the initial and minimum
temperatures of the particles.

* The initial temperature of all particles: ``initial_temperature`` (Default:
  InternalEnergy from the initial conditions)
* The minimal temperature of any particle: ``minimal_temperature`` (Default: 0)
* The mass fraction of hydrogen used to set the initial temperature:
  ``H_mass_fraction`` (Default: 0.755)
* The ionization temperature (from neutral to ionized) for primordial gas,
  again used in this conversion: ``H_ionization_temperature`` (Default: 1e4)

These parameters, if not present, are set to the default values. The initial
temperature is used, along with the hydrogen mass fraction and ionization
temperature, to set the initial internal energy per unit mass (or entropy per
unit mass) of the particles.

Throughout the run, if the temperature of a particle drops below
``minimal_temperature``, the particle has energy added to it such that it
remains at that temperature. The run is not terminated prematurely. The
temperatures specified in this section are in internal units.

The full section to start a typical cosmological run would be:

.. code:: YAML

   SPH:
     resolution_eta:                     1.2
     CFL_condition:                      0.1
     h_tolerance:                        1e-4
     h_min_ratio:                        0.1
     h_max:                              1.    # U_L
     initial_temperature:                273   # U_T
     minimal_temperature:                100   # U_T
     H_mass_fraction:                    0.755
     H_ionization_temperature:           1e4   # U_T
     particle_splitting:                 1 
     particle_splitting_mass_threshold:  5e-3  # U_M

.. _Parameters_Stars:

Stars
-----

The ``Stars`` section is used to set parameters that describe the Stars
calculations when doing feedback or enrichment. Note that if stars only act
gravitationally (i.e. SWIFT is run *without* ``--feedback``) no parameters
in this section are used. 

The first four parameters are related to the neighbour search:

* The (relative) tolerance to converge smoothing lengths within:
  ``h_tolerance`` (Default: same as SPH scheme)
* The maximal smoothing length in internal units: ``h_max`` (Default: same
  as SPH scheme)
* The minimal allowed smoothing length in terms of the gravitational
  softening: ``h_min_ratio`` (Default: same as SPH scheme)
* The maximal (relative) allowed change in volume over one time-step:
  ``max_volume_change`` (Default: same as SPH scheme)

These four parameters are optional and will default to their SPH equivalent
if left unspecified. That is the value specified by the user in that
section or the default SPH value if left unspecified there as well.

The next four parameters govern the time-step size choices for star
particles. By default star particles get their time-step sizes set
solely by the condition based on gravity. Additional criteria can be
applied by setting some of the following parameters (the actual
time-step size is then the minimum of this criterion and of the gravity
criterion):

* Time-step size for young stars in Mega-years:
  ``max_timestep_young_Myr`` (Default: FLT_MAX)
* Time-step size for old stars in Mega-years: ``max_timestep_old_Myr``
  (Default: FLT_MAX)
* Age transition from young to old in Mega-years:
  ``timestep_age_threshold_Myr`` (Default: FLT_MAX)
* Age above which no time-step limit is applied in Mega-years:
  ``timestep_age_threshold_unlimited_Myr`` (Default: 0)

Star particles with ages above the unlimited threshold only use the
gravity condition. Star particles with ages below that limit use
either the young or old time-step sizes based on their ages. These
parameters effectively allow for three different age brackets with the
last age bracket imposing no time-step length.

The remaining parameters can be used to overwrite the birth time (or
scale-factor), birth density and birth temperatures (if these
quantities exist) of the stars that were read from the ICs. This can
be useful to start a simulation with stars already of a given age or
with specific (uniform and non-0) properties at birth. The parameters
are:

* Whether or not to overwrite *all* the birth times: ``overwrite_birth_time``
  (Default: 0)
* The value to use: ``birth_time``
* Whether or not to overwrite *all* the birth densities: ``overwrite_birth_density``
  (Default: 0)
* The value to use: ``birth_density``
* Whether or not to overwrite *all* the birth temperatures: ``overwrite_birth_temperature``
  (Default: 0)
* The value to use: ``birth_temperature``

Note that if the birth time is set to ``-1`` then the stars will never
enter any feedback or enrichment loop. When these values are not
specified, SWIFT will start and use the birth times specified in the
ICs. If no values are given in the ICs, the stars' birth times will be
zeroed, which can cause issues depending on the type of run performed.

.. _Parameters_Sinks:

Sinks
-----

Currently, there are two models for the sink particles, the Default model and the GEAR one. Their parameters are described below. To choose a model, configure the code with ``--with-sink=<model>``, where ``<model>`` can be ``none`` or ``GEAR``. To run with sink particles, add the option ``--sinks``.
Below you will find the description of the ``none`` which is the default model. For ``GEAR`` model, please refer to :ref:`sink_GEAR_model`.

By default, the code is configured with ``--with-sink=none``. Then, the ``DefaultSink`` section is used to set parameters that describe the sinks in this model. The unique parameter is the sink accretion radius (also called cut-off radius): ``cut_off_radius``.

Note that this model does not create sink particles or accrete gas. 

The full section is:

.. code:: YAML

  DefaultSink:
    cut_off_radius:        1e-3       # Cut off radius of the sink particles (in internal units). This parameter should be adapted with the resolution..


.. _Parameters_time_integration:

Time Integration
----------------

The ``TimeIntegration`` section is used to set some general parameters related to time
integration. In all cases, users have to provide a minimal and maximal time-step
size:

* Maximal time-step size: ``dt_max``
* Minimal time-step size: ``dt_min``

These quantities are expressed in internal units. All particles will have their
time-step limited by the maximal value on top of all the other criteria that may
apply to them (gravity acceleration, Courant condition, etc.). If a particle
demands a time-step size smaller than the minimum, SWIFT will abort with an
error message. This is a safe-guard against simulations that would never
complete due to the number of steps to run being too large. Note that in
cosmological runs, the meaning of these variables changes slightly. They do not
correspond to differences in time but in logarithm of the scale-factor. For
these runs, the simulation progresses in jumps of
:math:`\Delta\log(a)`. ``dt_max`` is then the maximally allowed change in
:math:`\Delta\log(a)` allowed for any particle in the simulation. This behaviour
mimics the variables of the smae name in the Gadget code.

When running a non-cosmological simulation, the user also has to provide the
time of the start and the time of the end of the simulation:

* Start time: ``time_begin``
* End time: ``time_end``

Both are expressed in internal units. The start time is typically set to ``0``
but SWIFT can handle any value here. For cosmological runs, these values are
ignored and the start- and end-points of the runs are specified by the start and
end scale-factors in the cosmology section of the parameter file.

Additionally, when running a cosmological volume, advanced users can specify the
value of the dimensionless pre-factor entering the time-step condition linked
with the motion of particles with respect to the background expansion and mesh
size. See the theory document for the exact equations. Note that we explicitly
ignore the ``Header/Time`` attribute in initial conditions files, and only read
the start and end times or scale factors from the parameter file.

* Dimensionless pre-factor of the maximal allowed displacement:
  ``max_dt_RMS_factor`` (default: ``0.25``)
* Whether or not only the gas particle masses should be considered for
  the baryon component of the calculation: ``dt_RMS_use_gas_only`` (default: ``0``)
  
These values rarely need altering. The second parameter is only
meaningful if a subgrid model produces star (or other) particles with
masses substantially smaller than the gas masses. See the theory
documents for the precise meanings.

A full time-step section for a non-cosmological run would be:

.. code:: YAML

  TimeIntegration:
    time_begin:   0    # Start time in internal units.
    time_end:     10.  # End time in internal units.
    dt_max:       1e-2
    dt_min:       1e-6

Whilst for a cosmological run, one would need:

.. code:: YAML

  TimeIntegration:
    dt_max:              1e-4
    dt_min:              1e-10
    max_dt_RMS_factor:   0.25     # Default optional value
    dt_RMS_use_gas_only: 0        # Default optional value
    
.. _Parameters_ICs:

Initial Conditions
------------------

The ``InitialConditions`` section of the parameter file contains all the options related to
the initial conditions. The main two parameters are

* The name of the initial conditions file: ``file_name``,
* Whether the problem uses periodic boundary conditions or not: ``periodic``.

The file path is relative to where the code is being executed. These
parameters can be complemented by some optional values to drive some
specific behaviour of the code.

* Whether to generate gas particles from the DM particles: ``generate_gas_in_ics`` (default: ``0``),
* Whether to activate an additional clean-up of the SPH smoothing lengths: ``cleanup_smoothing_lengths`` (default: ``0``)

The procedure used to generate gas particles from the DM ones is
outlined in the theory documents and is too long for a full
description here.  The cleaning of the smoothing lengths is an
expensive operation but can be necessary in the cases where the
initial conditions are of poor quality and the values of the smoothing
lengths are far from the values they should have.

When starting from initial conditions created for Gadget, some
additional flags can be used to convert the values from h-full to
h-free and remove the additional :math:`\sqrt{a}` in the velocities:

* Whether to re-scale all the fields to remove powers of h from the quantities: ``cleanup_h_factors`` (default: ``0``),
* Whether to re-scale the velocities to remove the :math:`\sqrt{a}` assumed by Gadget : ``cleanup_velocity_factors`` (default: ``0``).

The h-factors are self-consistently removed according to their units
and this is applied to all the quantities irrespective of particle
types. The correct power of ``h`` is always calculated for each
quantity.

Finally, SWIFT also offers these options:

* A factor to re-scale all the smoothing-lengths by a fixed amount: ``smoothing_length_scaling`` (default: ``1.``),
* A shift to apply to all the particles: ``shift`` (default: ``[0.0,0.0,0.0]``),
* Whether to replicate the box along each axis: ``replicate`` (default: ``1``).
* Whether to re-map the IDs to the range ``[0, N]`` and hence discard
  the original IDs from the IC file: ``remap_ids`` (default: ``0``).
  
The shift is expressed in internal units and will be written to the header of
the snapshots. The option to replicate the box is especially useful for
weak-scaling tests. When set to an integer >1, the box size is multiplied by
this integer along each axis and the particles are duplicated and shifted such
as to create exact copies of the simulation volume.

The remapping of IDs is especially useful in combination with the option to
generate increasing IDs when splitting gas particles as it allows for the
creation of a compact range of IDs beyond which the new IDs generated by
splitting can be safely drawn from. Note that, when ``remap_ids`` is
switched on the ICs do not need to contain a ``ParticleIDs`` field.

Both replication and remapping explicitly overwrite any particle IDs
provided in the initial conditions. This may cause problems for runs
with neutrino particles, as some models assume that that the particle
ID was used as a random seed for the Fermi-Dirac momentum. In this case,
the ``Neutrino:generate_ics`` option can be used to generate new initial
conditions based on the replicated or remapped IDs. See :ref:`Neutrinos`
for details.

* Name of a HDF5 group to copy from the ICs file(s): ``metadata_group_name`` (default: ``ICs_parameters``)

If the initial conditions generator writes a HDF5 group with the parameters
used to make the initial conditions, this group can be copied through to
the output snapshots by specifying its name.

The full section to start a DM+hydro run from Gadget DM-only ICs would
be:

.. code:: YAML

   InitialConditions:
     file_name:  my_ics.hdf5
     periodic:                    1
     cleanup_h_factors:           1
     cleanup_velocity_factors:    1
     generate_gas_in_ics:         1
     cleanup_smoothing_lengths:   1
     metadata_group_name:         ICs_parameters

.. _Parameters_constants:

Physical Constants
------------------

For some idealised test it can be useful to overwrite the value of some physical
constants; in particular the value of the gravitational constant and vacuum
permeability. SWIFT offers an optional parameter to overwrite the value of
:math:`G_N` and :math:`\mu_0`.

.. code:: YAML

   PhysicalConstants:
     G:    1
     mu_0: 1

Note that this set :math:`G` to the specified value in the internal system
of units. Setting a value of `1` when using the system of units (10^10 Msun,
Mpc, km/s) will mean that :math:`G_N=1` in these units [#f2]_ instead of the
normal value :math:`G_N=43.00927`. The same applies to :math:`\mu_0`.

This option is only used for specific tests and debugging. This entire
section of the YAML file can typically be left out. More constants may
be handled in the same way in future versions.

.. _Parameters_snapshots:

Snapshots
---------

The ``Snapshots`` section of the parameter file contains all the options related to
the dump of simulation outputs in the form of HDF5 :ref:`snapshots`. The main
parameter is the base name that will be used for all the outputs in the run:

* The base name of the HDF5 snapshots: ``basename``.

This name will then be appended by an under-score and 4 digits followed by
``.hdf5`` (e.g. ``base_name_1234.hdf5``). The 4 digits are used to label the
different outputs, starting at ``0000``. In the default setup the digits simply
increase by one for each snapshot. (See :ref:`Output_list_label` to change that
behaviour.)

The time of the first snapshot is controlled by the two following options:

* Time of the first snapshot (non-cosmological runs): ``time_first``,
* Scale-factor of the first snapshot (cosmological runs): ``scale_factor_first``.

One of those two parameters has to be provided depending on the type of run. In
the case of non-cosmological runs, the time of the first snapshot is expressed
in the internal units of time. Users also have to provide the difference in time
(or scale-factor) between consecutive outputs:

* Time difference between consecutive outputs: ``delta_time``.

In non-cosmological runs this is also expressed in internal units. For
cosmological runs, this value is *multiplied* to obtain the
scale-factor of the next snapshot. This implies that the outputs are
equally spaced in :math:`\log(a)` (See :ref:`Output_list_label` to have
snapshots not regularly spaced in time).

The location and naming of the snapshots is altered by the following options:

* Directory in which to write snapshots: ``subdir``.
  (default: empty string).

If this is set then the full path to the snapshot files will be generated
by taking this value and appending a slash and then the snapshot file name
described above - e.g. ``subdir/base_name_1234.hdf5``. The directory is
created if necessary. Note however, that the sub-directories are created
when writing the first snapshot of a given category; the onus is hence on
the user to ensure correct writing permissions ahead of that time. Any
VELOCIraptor output produced by the run is also written to this directory.

When running the code with structure finding activated, it is often
useful to have a structure catalog written at the same simulation time
as the snapshots. To activate this, the following parameter can be
switched on:

* Run VELOCIraptor every time a snapshot is dumped: ``invoke_stf``
  (default: ``0``).

This produces catalogs using the options specified for the stand-alone
VELOCIraptor outputs (see the section :ref:`Parameters_structure_finding`) but
with a base name and output number that matches the snapshot name
(e.g. ``stf_base_name_1234.hdf5``) irrespective of the name specified in the
section dedicated to VELOCIraptor. Note that the invocation of VELOCIraptor at
every dump is done additionally to the stand-alone dumps that can be specified
in the corresponding section of the YAML parameter file. When running with
_more_ calls to VELOCIraptor than snapshots, gaps between snapshot numbers will
be created to accommodate for the intervening VELOCIraptor-only catalogs.

It is also possible to run the FOF algorithm just before writing each snapshot.

* Run FOF every time a snapshot is dumped: ``invoke_fof``
  (default: ``0``).

See the section :ref:`Parameters_fof` for details of the FOF parameters.

It is also possible to run the power spectrum calculation just before writing
each snapshot.

* Run PS every time a snapshot is dumped: ``invoke_ps``
  (default: ``0``).

See the section :ref:`Parameters_ps` for details of the power spectrum parameters.

When running over MPI, users have the option to split the snapshot over more
than one file. This can be useful if the parallel-io on a given system is slow
but has the drawback of producing many files per time slice. This is activated
by setting the parameter:

* Distributed snapshots over MPI: ``distributed`` (default: ``0``).

This option has no effect when running the non-MPI version of the code. Note
also that unlike other codes, SWIFT does *not* let the users chose the number of
individual files over which a snapshot is distributed. This is set by the number
of MPI ranks used in a given run. The individual files of snapshot 1234 will
have the name ``base_name_1234.x.hdf5`` where when running on N MPI ranks, ``x``
runs from 0 to N-1. If HDF5 1.10.0 or a more recent version is available,
an additional meta-snapshot named ``base_name_1234.hdf5`` will be produced
that can be used as if it was a non-distributed snapshot. In this case, the
HDF5 library itself can figure out which file is needed when manipulating the
snapshot.

On Lustre filesystems [#f4]_ it is important to properly stripe files to achieve
a good writing speed. If the parameter ``lustre_OST_count`` is set to the number
of OSTs present on the system, then SWIFT will set the `stripe count` of each
distributed file to `1` and set each file's `stripe index` to the MPI rank
generating it modulo the OST count [#f5]_. If the parameter is not set then the
files will be created with the default system policy (or whatever was set for
the directory where the files are written). This parameter has no effect on
non-Lustre file systems and no effect if distributed snapshots are not used.

* The number of Lustre OSTs to distribute the single-striped distributed
  snapshot files over: ``lustre_OST_count`` (default: ``0``)


Users can optionally ask to randomly sub-sample the particles in the snapshots.
This is specified for each particle type individually:

* Whether to switch on sub-sampling: ``subsample``   
* Whether to switch on sub-sampling: ``subsample_fraction`` 

These are arrays of 7 elements defaulting to seven 0s if left unspecified. Each
entry corresponds to the particle type used in the initial conditions and
snapshots [#f3]_.  The ``subsample`` array is made of ``0`` and ``1`` to indicate which
particle types to subsample. The other array is a float between ``0`` and ``1``
indicating the fraction of particles to keep in the outputs.  Note that the
selection of particles is selected randomly for each individual
snapshot. Particles can hence not be traced back from output to output when this
is switched on.
  
Users can optionally specify the level of compression used by the HDF5 library
using the parameter:

* GZIP compression level of the HDF5 arrays: ``compression`` (default: ``0``).

The default level of ``0`` implies no compression and values have to be in the
range :math:`[0-9]`. This integer is passed to the i/o library and used for the
loss-less GZIP compression algorithm. The compression is applied to *all* the
fields in the snapshots. Higher values imply higher compression but also more
time spent deflating and inflating the data.  When compression is switched on
the SHUFFLE filter is also applied to get higher compression rates. Note that up
until HDF5 1.10.x this option is not available when using the MPI-parallel
version of the i/o routines.

When applying lossy compression (see :ref:`Compression_filters`), particles may
be be getting positions that are marginally beyond the edge of the simulation
volume. A small vector perpendicular to the edge can be added to the particles
to alleviate this issue. This can be switched on by setting the parameter
``use_delta_from_edge`` (default: ``0``) to ``1`` and the buffer size from the
edge ``delta_from_edge`` (default: ``0.``). An example would be when using
Mega-parsec as the base unit and using a filter rounding to the nearest 10
parsec (``DScale5``). Adopting a buffer of 10pc (``delta_from_edge:1e-5``) would
alleviate any possible issue of seeing particles beyond the simulation volume in
the snapshots. In all practical applications the shift would be << than the
softening.

Users can run a program after a snapshot is dumped to disk using the following
parameters:

* Use the extra command after snapshot creation: ``run_on_dump`` (default :``0``)
* Command to run after snapshot creation: ``dump_command`` (default: nothing)

These are particularly useful should you wish to submit a job for postprocessing
the snapshot after it has just been created. Your script will be invoked with
two parameters, the snapshot base-name, and the snapshot number that was just
output as a zero-padded integer. For example, if the base-name is "eagle" and
snapshot 7 was just dumped, with ``dump_command`` set to ``./postprocess.sh``,
then SWIFT will run ``./postprocess.sh eagle 0007``.

For some quantities, especially in the subgrid models, it can be advantageous to
start recording numbers at a fixed time before the dump of a snapshot. Classic
examples are an averaged star-formation rate or accretion rate onto BHs. For the
subgrid models that support it, the triggers can be specified by setting the
following parameters:

* for gas: ``recording_triggers_part`` (no default, array of size set by each subgrid model)
* for stars: ``recording_triggers_spart`` (no default, array of size set by each subgrid model)
* for BHs: ``recording_triggers_bpart`` (no default, array of size set by each subgrid model)

The time is specified in internal time units (See the :ref:`Parameters_units`
section) and a recording can be ignored by setting the parameter to ``-1``. Note
that the code will verify that the recording time is smaller than the gap in
between consecutive snapshot dumps and if the recording window is longer, it
will reduce it to the gap size between the snapshots.

Finally, it is possible to specify a different system of units for the snapshots
than the one that was used internally by SWIFT. The format is identical to the
one described above (See the :ref:`Parameters_units` section) and read:

* a unit of length: ``UnitLength_in_cgs`` (default: ``InternalUnitSystem:UnitLength_in_cgs``),
* a unit of mass: ``UnitMass_in_cgs`` (default: ``InternalUnitSystem:UnitMass_in_cgs``),
* a unit of velocity ``UnitVelocity_in_cgs`` (default: ``InternalUnitSystem:UnitVelocity_in_cgs``),
* a unit of electric current ``UnitCurrent_in_cgs`` (default: ``InternalUnitSystem:UnitCurrent_in_cgs``),
* a unit of temperature ``UnitTemp_in_cgs`` (default: ``InternalUnitSystem:UnitTemp_in_cgs``).

When unspecified, these all take the same value as assumed by the internal
system of units. These are rarely used but can offer a practical alternative to
converting data in the post-processing of the simulations.

For a standard cosmological run with structure finding activated, the
full section would be:

.. code:: YAML

   Snapshots:
     basename:            output
     scale_factor_first:  0.02    # z = 49
     delta_time:          1.02
     invoke_stf:          1

Showing all the parameters for a basic non-cosmological hydro test-case, one
would have:

.. code:: YAML

   Snapshots:
     basename:            sedov
     subdir:              snapshots
     time_first:          0.01
     delta_time:          0.005
     invoke_stf:          0
     invoke_fof:          1
     compression:         3
     distributed:         1
     lustre_OST_count:   48   # System has 48 Lustre OSTs to distribute the files over
     UnitLength_in_cgs:   1.  # Use cm in outputs
     UnitMass_in_cgs:     1.  # Use grams in outputs
     UnitVelocity_in_cgs: 1.  # Use cm/s in outputs
     UnitCurrent_in_cgs:  1.  # Use Ampere in outputs
     UnitTemp_in_cgs:     1.  # Use Kelvin in outputs
     subsample:           [0, 1, 0, 0, 0, 0, 1]   # Sub-sample the DM and neutrinos
     subsample_fraction:  [0, 0.01, 0, 0, 0, 0, 0.1]  # Write 1% of the DM parts and 10% of the neutrinos
     run_on_dump:         1
     dump_command:        ./submit_analysis.sh
     use_delta_from_edge: 1
     delta_from_edge:     1e-6  # Move particles away from the edge by 1-e6 of the length unit.
     recording_triggers_part: [1.0227e-4, 1.0227e-5]   # Recording starts 100M and 10M years before a snapshot
     recording_triggers_spart: [-1, -1]                # No recording
     recording_triggers_bpart: [1.0227e-4, 1.0227e-5]  # Recording starts 100M and 10M years before a snapshot


Some additional specific options for the snapshot outputs are described in the
following pages:

* :ref:`Output_list_label` (to have snapshots not evenly spaced in time or with
  non-regular labels),
* :ref:`Output_selection_label` (to select what particle fields to write).

.. _Parameters_line_of_sight:

Line-of-sight outputs
---------------------

The ``LineOfSight`` section of the parameter file contains all the options related to
the dump of simulation outputs in the form of HDF5 :ref:`line_of_sight` data to
be processed by the ``SpecWizard`` tool
(See `Theuns et al. 1998 <https://ui.adsabs.harvard.edu/abs/1998MNRAS.301..478T/>`_,
`Tepper-Garcia et al. 2011
<https://ui.adsabs.harvard.edu/abs/2011MNRAS.413..190T/>`_). The parameters are:

.. code:: YAML

   LineOfSight:
     basename:            los
     scale_factor_first:  0.02    # Only used when running in cosmological mode
     delta_time:          1.02
     time_first:          0.01    # Only used when running in non-cosmological mode
     output_list_on:      0       # Overwrite the regular output times with a list of output times
     num_along_x:         0
     num_along_y:         0
     num_along_z:         100
     allowed_los_range_x: [0, 100.]   # Range along the x-axis where LoS along Y or Z are allowed
     allowed_los_range_y: [0, 100.]   # Range along the y-axis where LoS along X or Z are allowed
     allowed_los_range_z: [0, 100.]   # Range along the z-axis where LoS along X or Y are allowed
     range_when_shooting_down_x: 100. # Range along the x-axis of LoS along x
     range_when_shooting_down_y: 100. # Range along the y-axis of LoS along y
     range_when_shooting_down_z: 100. # Range along the z-axis of LoS along z


.. _Parameters_light_cone:

Light Cone Outputs
---------------------

One or more light cone outputs can be configured by including ``LightconeX`` sections
in the parameter file, where X is in the range 0-7. It is also possible to include a
``LightconeCommon`` section for parameters which are the same for all lightcones. The
parameters for each light cone are:

* Switch to enable or disable a lightcone: ``enabled``

This should be set to 1 to enable the corresponding lightcone or 0 to disable it.
Has no effect if specified in the LightconeCommon section.

* Directory in which to write light cone output: ``subdir``

All light cone output files will be written in the specified directory.

* Base name for particle and HEALPix map outputs: ``basename``.

Particles will be written to files ``<basename>_XXXX.Y.hdf5``, where XXXX numbers the files
written by a single MPI rank and Y is the MPI rank index. HEALPix maps are written to files
with names ``<basename>.shell_X.hdf5``, where X is the index of the shell. The basename must
be unique for each light cone so it cannot be specified in the LightconeCommon section.

See :ref:`lightcone_adding_outputs_label` for information on adding new output quantities.

* Location of the observer in the simulation box, in internal units: ``observer_position``

* Size of in memory chunks used to store particles and map updates: ``buffer_chunk_size``

During each time step buffered particles and HEALPix map updates are stored in a linked
list of chunks of ``buffer_chunk_size`` elements. Additional chunks are allocated as needed.
The map update process is parallelized over chunks so the chunks should be small enough that
each MPI rank typically has more chunks than threads.

* Maximum amount of map updates (in MB) to send on each iteration: ``max_map_update_send_size_mb``

Flushing the map update buffer involves sending the updates to the MPI ranks with the affected
pixel data. Sending all updates at once can consume a large amount of memory so this parameter
allows updates to be applied over multiple iterations to reduce peak memory usage.

* Redshift range to output each particle type: ``z_range_for_<type>``

A two element array with the minimum and maximum redshift at which particles of type ``<type>``
will be output as they cross the lightcone. ``<type>`` can be Gas, DM, DMBackground, Stars, BH
or Neutrino. If this parameter is not present for a particular type then that type will not
be output.

* The number of buffered particles which triggers a write to disk: ``max_particles_buffered``

If an MPI rank has at least max_particles_buffered particles which have crossed the lightcone,
it will write them to disk at the end of the current time step.

* Size of chunks in the particle output file

This sets the HDF5 chunk size. Particle outputs must be chunked because the number of particles
which will be written out is not known when the file is created.

* Whether to use lossy compression in the particle outputs: ``particles_lossy_compression``

If this is 1 then the HDF5 lossy compression filter named in the definition of each particle
output field will be enabled. If this is 0 lossy compression is not applied.

* Whether to use lossless compression in the particle outputs: ``particles_gzip_level``

If this is non-zero the HDF5 deflate filter will be applied to lightcone particle output with
the compression level set to the specified value. 

* HEALPix map resolution: ``nside``

* Name of the file with shell radii: ``radius_file.txt``

This specifies the name of a file with the inner and outer radii of the shells used to make
HEALPix maps. It should be a text file with a one line header and then two comma separated columns
of numbers with the inner and outer radii. The units are determined by the header. The header must
be one of the following:

``# Minimum comoving distance, Maximum comoving distance``,
``# Minimum redshift, Maximum redshift``, or
``# Maximum expansion factor, Minimum expansion factor``. Comoving distances are in internal units.
The shells must be in ascending order of radius and must not overlap.

* Number of pending HEALPix map updates before the buffers are flushed: ``max_updates_buffered``

In MPI mode applying updates to the HEALPix maps requires communication and forces synchronisation
of all MPI ranks, so it is not done every time step. If any MPI rank has at least
``max_updates_buffered`` pending updates at the end of a time step, then all ranks will apply
their updates to the HEALPix maps.

* Which types of HEALPix maps to create: ``map_names_file``

This is the name of a file which specifies what quantities should be accumulated to HEALPix maps.
The possible map types are defined in the lightcone_map_types array in ``lightcone_map_types.h``.
See :ref:`lightcone_adding_outputs_label` if you'd like to add a new map type.

* Whether to distribute HEALPix maps over multiple files: ``distributed_maps``

If this is 0 then the code uses HDF5 collective writes to write each map to a single file. If this
is 1 then each MPI rank writes its part of the HEALPix map to a separate file.

The file contains two columns: the first column is the name of the map type and the second is the
name of the compression filter to apply to it. See io_compression.c for the list of compression
filter names. Set the filter name to ``on`` to disable compression.

* Whether to use lossless compression in the HEALPix map outputs: ``maps_gzip_level``

If this is non-zero the HDF5 deflate filter will be applied to the lightcone map output with
the compression level set to the specified value. 

The following shows a full set of light cone parameters for the case where we're making two
light cones which only differ in the location of the observer:

.. code:: YAML

  LightconeCommon:

    # Common parameters
    subdir:            lightcones
    buffer_chunk_size:      100000
    max_particles_buffered: 1000000
    hdf5_chunk_size:        10000
 
    # Redshift ranges for particle types
    z_range_for_Gas:           [0.0, 0.05]
    z_range_for_DM:            [0.0, 0.05]
    z_range_for_DMBackground:  [0.0, 0.05]
    z_range_for_Stars:         [0.0, 0.05]
    z_range_for_BH:            [0.0, 0.05]
    z_range_for_Neutrino:      [0.0, 0.05]
    
    # Healpix map parameters
    nside:                512
    radius_file:          ./shell_radii.txt
    max_updates_buffered: 100000
    map_names_file:       map_names.txt
    max_map_update_send_size_mb: 1.0
    distributed_maps:     0

    # Compression options
    particles_lossy_compression: 0
    particles_gzip_level:        6
    maps_gzip_level:             6

  Lightcone0:

    enabled:  1
    basename: lightcone0
    observer_position: [35.5, 78.12, 12.45]

  Lightcone0:

    enabled:  1
    basename: lightcone1
    observer_position: [74.2, 10.80, 53.59]
  

An example of the radius file::

  # Minimum comoving distance, Maximum comoving distance
  0.0,   50.0
  50.0,  100.0
  150.0, 200.0
  200.0, 400.0
  400.0, 1000.0

An example of the map names file::

  TotalMass         on
  SmoothedGasMass   on
  UnsmoothedGasMass on
  DarkMatterMass    on


.. _Parameters_eos:

Equation of State (EoS)
-----------------------

The ``EoS`` section contains options for the equations of state.
Multiple EoS can be used for :ref:`planetary`,
see :ref:`planetary_eos` for more information. 

To enable one or multiple EoS, the corresponding ``planetary_use_*:``
flag(s) must be set to ``1`` in the parameter file for a simulation,
along with the path to any table files, which are set by the 
``planetary_*_table_file:`` parameters.

For the (non-planetary) isothermal EoS, the ``isothermal_internal_energy:``
parameter sets the thermal energy per unit mass.

.. code:: YAML

   EoS:
     isothermal_internal_energy: 20.26784      # Thermal energy per unit mass for the case of isothermal equation of state (in internal units).
     barotropic_vacuum_sound_speed: 2e4        # Vacuum sound speed in the case of the barotropic equation of state (in internal units).
     barotropic_core_density:       1e-13      # Core density in the case of the barotropic equation of state (in internal units).
     # Select which planetary EoS material(s) to enable for use.
     planetary_use_idg_def:    0               # Default ideal gas, material ID 0
     planetary_use_Til_iron:       1           # Tillotson iron, material ID 100
     planetary_use_Til_granite:    1           # Tillotson granite, material ID 101
     planetary_use_Til_water:      0           # Tillotson water, material ID 102
     planetary_use_Til_basalt:     0           # Tillotson basalt, material ID 103
     planetary_use_HM80_HHe:   0               # Hubbard & MacFarlane (1980) hydrogen-helium atmosphere, material ID 200
     planetary_use_HM80_ice:   0               # Hubbard & MacFarlane (1980) H20-CH4-NH3 ice mix, material ID 201
     planetary_use_HM80_rock:  0               # Hubbard & MacFarlane (1980) SiO2-MgO-FeS-FeO rock mix, material ID 202
     planetary_use_SESAME_iron:    0           # SESAME iron 2140, material ID 300
     planetary_use_SESAME_basalt:  0           # SESAME basalt 7530, material ID 301
     planetary_use_SESAME_water:   0           # SESAME water 7154, material ID 302
     planetary_use_SS08_water:     0           # Senft & Stewart (2008) SESAME-like water, material ID 303
     planetary_use_ANEOS_forsterite:   0       # ANEOS forsterite (Stewart et al. 2019), material ID 400
     planetary_use_ANEOS_iron:         0       # ANEOS iron (Stewart 2020), material ID 401
     planetary_use_ANEOS_Fe85Si15:     0       # ANEOS Fe85Si15 (Stewart 2020), material ID 402
     # Tablulated EoS file paths.
     planetary_HM80_HHe_table_file:    ./EoSTables/HM80_HHe.txt
     planetary_HM80_ice_table_file:    ./EoSTables/HM80_ice.txt
     planetary_HM80_rock_table_file:   ./EoSTables/HM80_rock.txt
     planetary_SESAME_iron_table_file:     ./EoSTables/SESAME_iron_2140.txt
     planetary_SESAME_basalt_table_file:   ./EoSTables/SESAME_basalt_7530.txt
     planetary_SESAME_water_table_file:    ./EoSTables/SESAME_water_7154.txt
     planetary_SS08_water_table_file:      ./EoSTables/SS08_water.txt
     planetary_ANEOS_forsterite_table_file:    ./EoSTables/ANEOS_forsterite_S19.txt
     planetary_ANEOS_iron_table_file:          ./EoSTables/ANEOS_iron_S20.txt
     planetary_ANEOS_Fe85Si15_table_file:      ./EoSTables/ANEOS_Fe85Si15_S20.txt

.. _Parameters_ps:

Power Spectra Calculation
-------------------------


SWIFT can compute a variety of auto- and cross- power spectra at user-specified
intervals. The behaviour of this output type is governed by the ``PowerSpectrum``
section of the parameter file. The calculation is performed on a regular grid
(typically of size 256^3) and foldings are used to extend the range probed to
smaller scales.

The options are:

 * The size of the base grid to perform the PS calculation:
   ``grid_side_length``.
 * The number of grid foldings to use: ``num_folds``.
 * The factor by which to fold at each iteration: ``fold_factor`` (default: 4)
 * The order of the window function: ``window_order`` (default: 3)
 * Whether or not to correct the placement of the centre of the k-bins for small k values: ``shift_centre_small_k_bins`` (default: 1)

The window order sets the way the particle properties get assigned to the mesh.
Order 1 corresponds to the nearest-grid-point (NGP), order 2 to cloud-in-cell
(CIC), and order 3 to triangular-shaped-cloud (TSC). Higher-order schemes are not
implemented.

Finally, the quantities for which a PS should be computed are specified as a
list of pairs of values for the parameter ``requested_spectra``.  Auto-spectra
are specified by using the same type for both pair members. The available values
listed in the following table:

+---------------------+---------------------------------------------------+
| Name                | Description                                       |
+=====================+===================================================+
| ``matter``          | Mass density of all matter                        |
+---------------------+---------------------------------------------------+
| ``cdm``             | Mass density of all dark matter                   |
+---------------------+---------------------------------------------------+
| ``gas``             | Mass density of all gas                           |
+---------------------+---------------------------------------------------+
| ``starBH``          | Mass density of all stars and BHs                 |
+---------------------+---------------------------------------------------+
| ``neutrino``        | Mass density of all neutrinos                     |
+---------------------+---------------------------------------------------+
| ``neutrino1``       | Mass density of a random half of the neutrinos    |
+---------------------+---------------------------------------------------+
| ``neutrino2``       | Mass density of a the other half of the neutrinos |
+---------------------+---------------------------------------------------+
| ``pressure``        | Electron pressure                                 |
+---------------------+---------------------------------------------------+

A dark matter mass density auto-spectrum is specified as ``cdm-cdm`` and a gas
density - electron pressure cross-spectrum as ``gas-pressure``.

The ``neutrino1`` and ``neutrino2`` selections are based on the particle IDs and
are mutually exclusive. The particles selected in each half are different in
each output. Note that neutrino PS can only be computed when neutrinos are
simulated using particles.

SWIFT uses bins of integer :math:`k`, with bins :math:`[0.5,1.5]`, :math:`[1.5,2.5]` etc.  The
representative :math:`k` values used to be assigned to the bin centres (so k=1, 2, etc), which are
then transformed to physical :math:`k` by a factor :math:`kL/(2*pi)`. For the first few bins, only
few modes contribute to each bin. It is then advantageous to move the "centre" of the bin to the
actual location correponding to the mean of the contributing modes. The :math:`k` label of the bin
is thus shifted by a small amount. The way to calculate these shifts is to consider a 3D cube of
:math:`(kx,ky,kz)` cells and check which cells fall inside a spherical shell with boundaries
:math:`(i+0.5,i+1.5)`, then calculate the average :math:`k=\sqrt{kx^2+ky^2+kz^2}`. So for
:math:`i=0` there cells :math:`k=1` and 12 cells :math:`k=\sqrt(2)`, so the weighted k becomes
:math:`(6 * 1 + 12 * \sqrt(2)) / 18 = 1.2761424`. Note that only the first 7 (22) bins require a
correction larger than 1 (0.1) percent. We apply a correction to the first 128 terms. This
correction is activated when ``shift_centre_small_k_bins`` is switched on (that's the default
behaviour).

An example of a valid power-spectrum section of the parameter file looks like:

.. code:: YAML

  PowerSpectrum:
    grid_side_length:  256
    num_folds:         3
    requested_spectra: ["matter-matter", "cdm-cdm", "cdm-matter"] # Total-matter and CDM auto-spectra + CDM-total cross-spectrum

Some additional specific options for the power-spectra outputs are described in the
following pages:

* :ref:`Output_list_label` (to have PS not evenly spaced in time)


.. _Parameters_fof:

Friends-Of-Friends (FOF)
------------------------

The parameters are described separately on the page
:ref:`Fof_Parameter_Description_label` within the more general
:ref:`Friends_Of_Friends_label` description.

.. _Parameters_statistics:

Statistics
----------

Some additional specific options for the statistics outputs are described in the
following page:

* :ref:`Output_list_label` (to have statistics outputs not evenly spaced in time).

.. _Parameters_restarts:

Restarts
--------

SWIFT can write check-pointing files and restart from them. The behaviour of
this mechanism is driven by the options in the ``Restarts`` section of the YAML
parameter file. All the parameters are optional but default to values that
ensure a reasonable behaviour.

* Whether or not to enable the dump of restart files: ``enable`` (default:
  ``1``).

This parameter acts a master-switch for the check-pointing capabilities. All the
other options require the ``enable`` parameter to be set to ``1``.

* Whether or not to save a copy of the previous set of check-pointing files:
  ``save`` (default: ``1``),
* Whether or not to dump a set of restart file on regular exit: ``onexit``
  (default: ``0``),
* The wall-clock time in hours between two sets of restart files:
  ``delta_hours`` (default: ``5.0``).

Note that there is no buffer time added to the ``delta_hours`` value. If the
system's batch queue run time limit is set to 5 hours, the user must specify a
smaller value to allow for enough time to safely dump the check-point files.

* The sub-directory in which to store the restart files: ``subdir`` (default:
  ``restart``),
* The basename of the restart files: ``basename`` (default: ``swift``)

If the directory does not exist, SWIFT will create it.  When resuming a run,
SWIFT, will look for files with the name provided in the sub-directory specified
here. The files themselves are named ``basename_000001.rst`` where the basename
is replaced by the user-specified name and the 6-digits number corresponds to
the MPI-rank. SWIFT writes one file per MPI rank. If the ``save`` option has
been activated, the previous set of restart files will be named
``basename_000000.rst.prev``.

On Lustre filesystems [#f4]_ it is important to properly stripe files to achieve
a good writing and reading speed. If the parameter ``lustre_OST_count`` is set
to the number of OSTs present on the system, then SWIFT will set the `stripe
count` of each restart file to `1` and set each file's `stripe index` to the MPI
rank generating it modulo the OST count [#f5]_. If the parameter is not set then
the files will be created with the default system policy (or whatever was set
for the directory where the files are written). This parameter has no effect on
non-Lustre file systems.

* The number of Lustre OSTs to distribute the single-striped restart files over:
  ``lustre_OST_count`` (default: ``0``)

SWIFT can also be stopped by creating an empty file called ``stop`` in the
directory where the restart files are written (i.e. the directory speicified by
the parameter ``subdir``). This will make SWIFT dump a fresh set of restart file
(irrespective of the specified ``delta_time`` between dumps) and exit
cleanly. One parameter governs this behaviour:

* Number of steps between two checks for the presence of a ``stop`` file:
  ``stop_steps`` (default: ``100``).

The default value is chosen such that SWIFT does not need to poll the
file-system to often, which can take a significant amount of time on distributed
systems. For runs where the small time-steps take a much larger amount of time,
a smaller value is recommended to allow for a finer control over when the code
can be stopped.

Finally, SWIFT can automatically stop after a specified amount of wall-clock
time. The code can also run a command when exiting in this fashion, which can be
used, for instance, to interact with the batch queue system:

* Maximal wall-clock run time in hours: ``max_run_time`` (default: ``24.0``),
* Whether or not to run a command on exit: ``resubmit_on_exit`` (default:
  ``0``),
* The command to run on exit: ``resubmit_command`` (default: ``./resub.sh``).

Note that no check is performed on the validity of the command to run. SWIFT
simply calls ``system()`` with the user-specified command.

To run SWIFT, dumping check-pointing files every 6 hours and running for 24
hours after which a shell command will be run, one would use:

.. code:: YAML

  Restarts:
    enable:             1
    save:               1          # Keep copies
    onexit:             0
    subdir:             restart    # Sub-directory of the directory where SWIFT is run
    basename:           swift
    delta_hours:        5.0
    stop_steps:         100
    max_run_time:       24.0       # In hours
    lustre_OST_count:   48         # System has 48 Lustre OSTs to distribute the files over
    resubmit_on_exit:   1
    resubmit_command:   ./resub.sh

.. _Parameters_scheduler:

Scheduler
---------

The Scheduler section contains various parameters that control how the cell
tree is configured and defines some values for the related tasks.  In general
these should be considered as tuning parameters, both for speed and memory
use.

.. code:: YAML

   nr_queues: 0

Defines the number of task queues used. These are normally set to one per
thread and should be at least that number.

A number of parameters decide how the cell tree will be split into sub-cells,
according to the number of particles and their expected interaction count,
and the type of interaction. These are:

.. code:: YAML

  cell_max_size:             8000000
  cell_sub_size_pair_hydro:  256000000
  cell_sub_size_self_hydro:  32000
  cell_sub_size_pair_grav:   256000000
  cell_sub_size_self_grav:   32000
  cell_sub_size_pair_stars:  256000000
  cell_sub_size_self_stars:  32000
  cell_split_size:           400

when possible cells that exceed these constraints will be split into a further
level of sub-cells. So for instance a sub-cell should not contain more than
400 particles (this number defines the scale of most `N*N` interactions).

To control the number of self-gravity tasks we have the parameter:

.. code:: YAML

  cell_subdepth_diff_grav:   4

which stops these from being done at the scale of the leaf cells, of which
there can be a large number. In this case cells with gravity tasks must be at
least 4 levels above the leaf cells (when possible).

To control the depth at which the ghost tasks are placed, there are two
parameters (one for the gas, one for the stars). These specify the maximum
number of particles allowed in such a task before splitting into finer ones. A
similar parameter exists for the cooling tasks, which can be useful to tweak for
models in which the cooling operations are expensive. These three parameters
are:

.. code:: YAML

  engine_max_parts_per_ghost:    1000
  engine_max_sparts_per_ghost:   1000
  engine_max_parts_per_cooling: 10000


Extra space is required when particles are created in the system (to the time
of the next rebuild). These are controlled by:

.. code:: YAML

  cell_extra_parts:          0
  cell_extra_gparts:         0
  cell_extra_sparts:         400


The number of top-level cells is controlled by the parameter:

.. code:: YAML

  max_top_level_cells:       12

this is the number per dimension, we will have 12x12x12 cells. There must be
at least 3 top-level cells per dimension.

The number of top-level cells should be set so that the number of particles
per cell is not too large, this is particularly important when using MPI as
this defines the maximum size of cell exchange and also the size of non-local
cells (these are used for cell interactions with local cells), which can have
a large influence on memory use. Best advice for this is to at least scale for
additional nodes.

The memory used for holding the task and task-link lists needs to be
pre-allocated, but cannot be pre-calculated, so we have the two parameters:

.. code:: YAML

  tasks_per_cell:            0.0
  links_per_tasks:           10

which are guesses at the mean numbers of tasks per cell and number of links
per task. The tasks_per_cell value will be conservatively guessed when set to
0.0, but you will be able to save memory by setting a value. The way to get a
better estimate is to run SWIFT with verbose reporting on (```--verbose=1```)
and check for the lines that report the ```per cell``` or with MPI ``maximum
per cell``` values. This number can vary as the balance between MPI ranks
does, so it is probably best to leave some head room.

If these are exceeded you should get an obvious error message.

Finally the parameter:

.. code:: YAML

  mpi_message_limit:         4096

Defines the size (in bytes) below which MPI communication will be sent using
non-buffered calls. These should have lower latency, but how that works or
is honoured is an implementation question.


.. _Parameters_domain_decomposition:

Domain Decomposition:
---------------------

This section determines how the top-level cells are distributed between the
ranks of an MPI run. An ideal decomposition should result in each rank having
a similar amount of work to do, so that all the ranks complete at the same
time. Achieving a good balance requires that SWIFT is compiled with either the
ParMETIS or METIS libraries. ParMETIS is an MPI version of METIS, so is
preferred for performance reasons.

When we use ParMETIS/METIS the top-level cells of the volume are considered as
a graph, with a cell at each vertex and edges that connect the vertices to all
the neighbouring cells (so we have 26 edges connected to each vertex).
Decomposing such a graph into domains is known as partitioning, so in SWIFT we
refer to domain decomposition as partitioning.

This graph of cells can have weights associated with the vertices and the
edges. These weights are then used to guide the partitioning, seeking to
balance the total weight of the vertices and minimize the weights of the edges
that are cut by the domain boundaries (known as the edgecut). We can consider
the edge weights as a proxy for the exchange of data between cells, so
minimizing this reduces communication.

The Initial Partition:
^^^^^^^^^^^^^^^^^^^^^^

When SWIFT first starts it reads the initial conditions and then does an
initial distribution of the top-level cells. At this time the only information
available is the cell structure and, by geometry, the particles each cell
should contain. The type of partitioning attempted is controlled by the::

  DomainDecomposition:
    initial_type:

parameter. Which can have the values *memory*, *edgememory*, *region*, *grid* or
*vectorized*:

    * *edgememory*

    This is the default if METIS or ParMETIS is available. It performs a
    partition based on the memory use of all the particles in each cell.
    The total memory per cell is used to weight the cell vertex and all the
    associated edges. This attempts to equalize the memory used by all the
    ranks but with some consideration given to the need to not cut dense
    regions (by also minimizing the edge cut). How successful this
    attempt is depends on the granularity of cells and particles and the
    number of ranks, clearly if most of the particles are in one cell, or a
    small region of the volume, balance is impossible or difficult. Having
    more top-level cells makes it easier to calculate a good distribution
    (but this comes at the cost of greater overheads).

    * *memory*

    This is like *edgememory*, but doesn't include any edge weights, it should
    balance the particle memory use per rank more exactly (but note effects
    like the numbers of cells per rank will also have an effect, as that
    changes the need for foreign cells).

    * *region*

    The one other METIS/ParMETIS option is "region". This attempts to assign equal
    numbers of cells to each rank, with the surface area of the regions minimised.

If ParMETIS and METIS are not available two other options are possible, but
will give a poorer partition:

    * *grid*

    Split the cells into a number of axis aligned regions. The number of
    splits per axis is controlled by the::

       initial_grid

    parameter. It takes an array of three values. The product of these values
    must equal the number of MPI ranks. If not set a suitable default will be used.

    * *vectorized*

    Allocate the cells on the basis of proximity to a set of seed
    positions. The seed positions are picked every nranks along a vectorized
    cell list (1D representation). This is guaranteed to give an initial
    partition for all cases when the number of cells is greater equal to the
    number of MPI ranks, so can be used if the others fail. Don't use this.

If ParMETIS and METIS are not available then only an initial partition will be
performed. So the balance will be compromised by the quality of the initial
partition.

Repartitioning:
^^^^^^^^^^^^^^^

When ParMETIS or METIS is available we can consider adjusting the balance
during the run, so we can improve from the initial partition and also track
changes in the run that require a different balance. The initial partition is
usually not optimal as although it may have balanced the distribution of
particles it has not taken account of the fact that different particles types
require differing amounts of processing and we have not considered that we
also need to do work requiring communication between cells. This latter point
is important as we are running an MPI job, as inter-cell communication may be
very expensive.

There are a number of possible repartition strategies which are defined using
the::

  DomainDecomposition:
    repartition_type:

parameter. The possible values for this are *none*, *fullcosts*, *edgecosts*,
*memory*, *timecosts*.

    * *none*

    Rather obviously, don't repartition. You are happy to run with the
    initial partition.

    * *fullcosts*

    Use computation weights derived from the running tasks for the vertex and
    edge weights. This is the default.

    * *edgecosts*

    Only use computation weights derived from the running tasks for the edge
    weights.

    * *memory*

    Repeat the initial partition with the current particle positions
    re-balancing the memory use.

    * *timecosts*

    Only use computation weights derived from the running tasks for the vertex
    weights and the expected time the particles will interact in the cells as
    the edge weights. Using time as the edge weight has the effect of keeping
    very active cells on single MPI ranks, so can reduce MPI communication.

The computation weights are actually the measured times, in CPU ticks, that
tasks associated with a cell take. So these automatically reflect the relative
cost of the different task types (SPH, self-gravity etc.), and other factors
like how well they run on the current hardware and are optimized by the
compiler used, but this means that we have a constraint on how often we can
consider repartitioning, namely when all (or nearly all) the tasks of the
system have been invoked in a step. To control this we have the::

    minfrac:     0.9

parameter. Which defines the minimum fraction of all the particles in the
simulation that must have been actively updated in the last step, before
repartitioning is considered.

That then leaves the question of when a run is considered to be out of balance
and should benefit from a repartition. That is controlled by the::

    trigger:          0.05

parameter. This value is the CPU time difference between MPI ranks, as a
fraction, if less than this value a repartition will not be
done. Repartitioning can be expensive not just in CPU time, but also because
large numbers of particles can be exchanged between MPI ranks, so is best
avoided.

If you are using ParMETIS there additional ways that you can tune the
repartition process.

METIS only offers the ability to create a partition from a graph, which means
that each solution is independent of those that have already been made, that
can make the exchange of particles very large (although SWIFT attempts to
minimize this), however, using ParMETIS we can use the existing partition to
inform the new partition, this has two algorithms that are controlled using::

    adaptive:         1

which means use adaptive repartition, otherwise simple refinement. The
adaptive algorithm is further controlled by the::

    itr:              100

parameter, which defines the ratio of inter node communication time to data
redistribution time, in the range 0.00001 to 10000000.0. Lower values give
less data movement during redistributions. The best choice for these can only
be determined by experimentation (the gains are usually small, so not really
recommended).

Finally we have the parameter::

    usemetis:         0

Forces the use of the METIS API, probably only useful for developers.

**Fixed cost repartitioning:**

So far we have assumed that repartitioning will only happen after a step that
meets the `minfrac:` and `trigger:` criteria, but we may want to repartition
at some arbitrary steps, and indeed do better than the initial partition
earlier in the run. This can be done using *fixed cost* repartitioning.

Fixed costs are output during each repartition step into the file
`partition_fixed_costs.h`, this should be created by a test run of your
full simulation (with possibly with a smaller volume, but all the physics
enabled). This file can then be used to replace the same file found in the
`src/` directory and SWIFT should then be recompiled. Once you have that, you
can use the parameter::

    use_fixed_costs:  1

to control whether they are used or not. If enabled these will be used to
repartition after the second step, which will generally give as good a
repartition immediately as you get at the first unforced repartition.

Also once these have been enabled you can change the ``trigger`` value to
numbers greater than 2, and repartitioning will be forced every ``trigger``
steps. This latter option is probably only useful for developers, but tuning
the second step to use fixed costs can give some improvements.

.. _Parameters_structure_finding:

Structure finding (VELOCIraptor)
--------------------------------

This section describes the behaviour of the on-the-fly structure
finding using the VELOCIraptor library (see
:ref:`VELOCIraptor_interface`). The section is named
``StructureFinding`` and also governs the behaviour of the
structure finding code when invoked at snapshots dumping time via
the parameter ``Snapshots:invoke_stf``.

The main parameters are:

 * The VELOCIraptor parameter file to use for the run:
   ``config_file_name``,
 * The directory in which the structure catalogs will be written: ``basename``.

Both these parameters must always be specified when running SWIFT with
on-the-fly calls to the structure finding code. In particular, when
only running VELOCIraptor when snapshots are written, nothing more is
necessary and one would use:

.. code:: YAML

  Snapshots:
    invoke_stf:        1                              # We want VELOCIraptor to be called when snapshots are dumped.
    # ...
    # Rest of the snapshots properties
	  
  StructureFinding:
    config_file_name:  my_stf_configuration_file.cfg  # See the VELOCIraptor manual for the content of this file.
    basename:          ./haloes/                      # Write the catalogs in this sub-directory
     
If one additionally want to call VELOCIraptor at times not linked with
snapshots, the additional parameters need to be supplied.

The time of the first call is controlled by the two following options:

* Time of the first call to VELOCIraptor (non-cosmological runs): ``time_first``,
* Scale-factor of the first call to VELOCIraptor (cosmological runs): ``scale_factor_first``.

One of those two parameters has to be provided depending on the type of run. In
the case of non-cosmological runs, the time of the first call is expressed
in the internal units of time. Users also have to provide the difference in time
(or scale-factor) between consecutive outputs:

* Time difference between consecutive outputs: ``delta_time``.

In non-cosmological runs this is also expressed in internal units. For
cosmological runs, this value is *multiplied* to obtain the
scale-factor of the next call. This implies that the outputs are
equally spaced in :math:`\log(a)` (See :ref:`Output_list_label` to have
calls not regularly spaced in time).

Since VELOCIraptor produces many small output files when running with MPI,
it can be useful to make a separate directory for each output time:

* Base name of directory created for each VELOCIraptor output: ``subdir_per_output``
  (default: empty string).

If this is set then a new directory is created each time VELOCIraptor is run.
The directory name will be subdir_per_output followed by the same output number
used in the filenames. Note that this directory is relative to the ``subdir`` parameter
from the Snapshots section if that is set.

By default this is an empty string, which means that all VELOCIraptor outputs will
be written to a single directory.

Showing all the parameters for a basic cosmologica test-case, one would have:

.. code:: YAML

   StructureFinding:
    config_file_name:     my_stf_configuration_file.cfg  # See the VELOCIraptor manual for the content of this file.
    basename:             haloes                         # Base name for VELOCIraptor output files
    subdir_per_output:    stf                            # Make a stf_XXXX subdirectory for each output
    scale_factor_first:   0.1                            # Scale-factor of the first output
    delta_time:           1.1                            # Delta log-a between outputs


Gravity Force Checks
--------------------

By default, when the code is configured with ``--enable-gravity-force-checks``,
the "exact" forces of all active gparts are computed during each timestep.

To give a bit more control over this, you can select to only perform the exact
force computation during the timesteps that all gparts are active, and/or only
at the timesteps when a snapshot is being dumped, i.e.,

.. code:: YAML

  ForceChecks:
    only_when_all_active:   1    # Only compute exact forces during timesteps when all gparts are active.
    only_at_snapshots:      1    # Only compute exact forces during timesteps when a snapshot is being dumped.

If ``only_when_all_active:1`` and ``only_at_snapshots:1`` are enabled together,
and all the gparts are not active during the timestep of the snapshot dump, the
exact forces computation is performed on the first timestep at which all the
gparts are active after that snapshot output timestep.

Neutrinos
---------

The ``Neutrino`` section of the parameter file controls the behaviour of
neutrino particles (``PartType6``). This assumes that massive neutrinos have
been specified in the ``Cosmology`` section described above. Random
Fermi-Dirac momenta will be generated if ``generate_ics`` is used. The
:math:`\delta f` method for shot noise reduction can be activated with
``use_delta_f``. Finally, a random seed for the Fermi-Dirac momenta can
be set with ``neutrino_seed``.

For mode details on the neutrino implementation, refer to :ref:`Neutrinos`. 
A complete specification of the model looks like

.. code:: YAML

  Neutrino:
    generate_ics:  1    # Replace neutrino particle velocities with random Fermi-Dirac momenta at the start
    use_delta_f:   1    # Use the delta-f method for shot noise reduction
    neutrino_seed: 1234 # A random seed used for the Fermi-Dirac momenta


------------------------
    
.. [#f1] The thorough reader (or overly keen SWIFT tester) would find  that the speed of light is :math:`c=1.8026\times10^{12}\,\rm{fur}\,\rm{ftn}^{-1}`, Newton's constant becomes :math:`G_N=4.896735\times10^{-4}~\rm{fur}^3\,\rm{fir}^{-1}\,\rm{ftn}^{-2}` and Planck's constant turns into :math:`h=4.851453\times 10^{-34}~\rm{fur}^2\,\rm{fir}\,\rm{ftn}^{-1}`.


.. [#f2] which would translate into a constant :math:`G_N=1.5517771\times10^{-9}~cm^{3}\,g^{-1}\,s^{-2}` if expressed in the CGS system.

.. [#f3] The mapping is 0 --> gas, 1 --> dark matter, 2 --> background dark
	 matter, 3 --> sinks, 4 --> stars, 5 --> black holes, 6 --> neutrinos.

.. [#f4] https://wiki.lustre.org/Main_Page

.. [#f5] We add a per-output random integer to the OST value such that we don't
	 generate a bias towards low OSTs. This averages the load over all OSTs
	 over the course of a run even if the number of OSTs does not divide the
	 number of files and vice-versa.
