.. GEAR sub-grid model
   Loic Hausammann, 17th April 2020


GEAR model
===========

GEAR's model are mainly described in `Revaz \& Jablonka <https://ui.adsabs.harvard.edu/abs/2018A%26A...616A..96R/abstract>`_.
This model can be selected with the configuration option ``--with-subgrid=GEAR`` and run with the option ``--gear``. A few examples exist and can be found in ``examples/GEAR``. 

.. _gear_pressure_floor:

Pressure Floor
~~~~~~~~~~~~~~

In order to avoid the artificial collapse of unresolved clumps, a minimum in pressure is applied to the particles.
This additional pressure can be seen as the pressure due to unresolved hydrodynamics turbulence and is given by:

.. math::
    P_\textrm{Jeans} = \frac{\rho}{\gamma} \frac{4}{\pi} G h^2 \rho N_\textrm{Jeans}^{2/3}

where :math:`\rho` is the density, :math:`\gamma` the adiabatic index, :math:`G` is the gravitational constant,
:math:`h` the kernel support and :math:`N_\textrm{Jeans}` (``GEARPressureFloor:jeans_factor`` in the parameter file) is the number of particle required in order to resolve a clump.


This must be directly implemented into the hydro schemes, therefore only a subset of schemes (Gadget-2, SPHENIX and Pressure-Energy) have the floor available.
In order to implement it, you need equation 12 in `Hopkins 2013 <https://arxiv.org/abs/1206.5006>`_:

.. math::
   m_i \frac{\mathrm{d}v_i}{\mathrm{d}t} = - \sum_j x_i x_j \left[ \frac{P_i}{y_i^2} f_{ij} \nabla_i W_{ij}(h_i) + \frac{P_j}{y_j^2} f_{ji} \nabla_j W_{ji}(h_j) \right]

and simply replace the :math:`P_i, P_j` by the pressure with the floor (when the pressure is below the floor).
Here the :math:`x, y` are simple weights that should never have the pressure floor included even if they are related to the pressure (e.g. pressure-entropy).

.. code:: YAML

   GEARPressureFloor:
    jeans_factor: 10.       # Number of particles required to suppose a resolved clump and avoid the pressure floor.



.. _gear_grackle_cooling:

Cooling: Grackle
~~~~~~~~~~~~~~~~
   
Grackle is a chemistry and cooling library presented in `B. Smith et al. 2017 <https://ui.adsabs.harvard.edu/abs/2017MNRAS.466.2217S>`_ 
(do not forget to cite if used).  Four different modes are available:
equilibrium, 6 species network (H, H\\( ^+ \\), e\\( ^- \\), He, He\\( ^+ \\)
and He\\( ^{++} \\)), 9 species network (adds H\\(^-\\), H\\(_2\\) and
H\\(_2^+\\)) and 12 species (adds D, D\\(^+\\) and HD).  Following the same
order, the swift cooling options are ``grackle_0``, ``grackle_1``, ``grackle_2``
and ``grackle_3`` (the numbers correspond to the value of
``primordial_chemistry`` in Grackle).  It also includes some metal cooling (on/off with ``GrackleCooling:with_metal_cooling``), self-shielding
methods and UV background (on/off with ``GrackleCooling:with_UV_background``).  In order to use the Grackle cooling, you will need
to provide a HDF5 table computed by Cloudy (``GrackleCooling:cloudy_table``).


When starting a simulation without providing the different element fractions in the non equilibrium mode, the code supposes an equilibrium and computes them automatically.
The code uses an iterative method in order to find the correct initial composition and this method can be tuned with two parameters. ``GrackleCooling:max_steps`` defines the maximal number of steps to reach the convergence and ``GrackleCooling:convergence_limit`` defines the tolerance in the relative error.

In order to compile SWIFT with Grackle, you need to provide the options ``with-chemistry=GEAR`` and ``with-grackle=$GRACKLE_ROOT``
where ``$GRACKLE_ROOT`` is the root of the install directory (not the ``lib``). 

.. warning::
  The actual Grackle version fully supported by SWIFT is 3.2.1. It can be downloaded from 
  `the official Grackle git repository <https://github.com/grackle-project/grackle/archive/refs/tags/grackle-3.2.1.tar.gz>`_.
  However, this version still had a bug when using threadsafe functions. Alternately, it is possible to get a fixed version
  using `the following fork frozen for compatibility with SWIFT <https://github.com/mladenivkovic/grackle-swift>`_.


To compile it, run
the following commands from the root directory of Grackle:
``./configure; cd src/clib``.
Update the variables ``LOCAL_HDF5_INSTALL`` and ``MACH_INSTALL_PREFIX`` in
the file ``src/clib/Make.mach.linux-gnu``.
Finish with ``make machine-linux-gnu; make && make install``.
Note that we require the 64 bit float version of Grackle, which should be the default setting. 
(The precision can be set while compiling grackle with ``make precision-64``).
If you encounter any problem, you can look at the `Grackle documentation <https://grackle.readthedocs.io/en/latest/>`_

You can now provide the path given for ``MACH_INSTALL_PREFIX`` to ``with-grackle``.

In the parameters file, a few different parameters are available.
``GrackleCooling:redshift`` defines the redshift to use for the UV background (for cosmological simulation, it must be set to -1 in order to use the simulation's redshift) and ``GrackleCooling:provide_*_heating_rates`` can enable the computation of user provided heating rates (such as with the radiative transfer) in either volumetric or specific units.

For the feedback, it can be made more efficient by turning off the cooling during a few Myr (``GrackleCooling:thermal_time_myr``) for the particles touched by a supernovae.

The self shielding method is defined by ``GrackleCooling:self_shielding_method`` where 0 means no self shielding, > 0 means a method defined in Grackle (see Grackle documentation for more information) and -1 means GEAR's self shielding that simply turn off the UV background when reaching a given density (``GrackleCooling:self_shielding_threshold_atom_per_cm3``).

A maximal (physical) density must be set with the ``GrackleCooling:maximal_density_Hpcm3 parameter``. The density passed to Grackle is *the minimum of this density and the gas particle (physical) density*. A negative value (:math:`< 0`) deactivates the maximal density, i.e. there is no maximal density limit.
The purpose of this parameter is the following. The Cloudy tables provided by Grackle are limited in density (typically to  :math:`10^4 \; \mathrm{hydrogen \; atoms/cm}^3`). In high-resolution simulations, particles can have densities higher than :math:`10^4 \; \mathrm{hydrogen \; atoms/cm}^3`. This maximal density ensures that we pass a density within the interpolation ranges of the table, should the density exceed it.
It can be a solution to some of the following errors (with a translation of what the values mean):

.. code:: text

	  inside if statement solve rate cool:           0           0
	  MULTI_COOL iter >        10000  at j,k =           1           1
	  FATAL error (2) in MULTI_COOL
	  dt =  1.092E-08 ttmin =  7.493E-12
	  2.8E-19    // sub-cycling timestep
	  7.5E-12    // time elapsed (in the sub-cycle)
	  2.2E+25    // derivative of the internal energy
	  T

.. note::
   This problem is particularly relevant with metal cooling enabled. Another solution is to modify the tables. But one is not exempted from exceeding the table maximal density value, since Grackle does not check if the particle density is larger than the table maximal density.

Here is the complete section in the parameter file:

.. code:: YAML

  GrackleCooling:
    cloudy_table: CloudyData_UVB=HM2012.h5       # Name of the Cloudy Table (available on the grackle bitbucket repository)
    with_UV_background: 1                        # Enable or not the UV background
    redshift: 0                                  # Redshift to use (-1 means time based redshift)
    with_metal_cooling: 1                        # Enable or not the metal cooling
    provide_volumetric_heating_rates: 0          # (optional) User provide volumetric heating rates
    provide_specific_heating_rates: 0            # (optional) User provide specific heating rates
    max_steps: 10000                             # (optional) Max number of step when computing the initial composition
    convergence_limit: 1e-2                      # (optional) Convergence threshold (relative) for initial composition
    thermal_time_myr: 5                          # (optional) Time (in Myr) for adiabatic cooling after a feedback event.
    self_shielding_method: -1                    # (optional) Grackle (1->3 for Grackle's ones, 0 for none and -1 for GEAR)
    self_shielding_threshold_atom_per_cm3: 0.007 # Required only with GEAR's self shielding. Density threshold of the self shielding
    maximal_density_Hpcm3:         1e4                 # Maximal density (in atoms/cm^3) for cooling. Higher densities are floored to this value to ensure grackle works properly when interpolating beyond the cloudy_table maximal density. A value < 0 deactivates this parameter.
    HydrogenFractionByMass : 1.                  # Hydrogen fraction by mass (default is 0.76)
    use_radiative_transfer : 1                   # Arrays of ionization and heating rates are provided
    RT_heating_rate_cgs    : 0                   # heating         rate in units of / nHI_cgs
    RT_HI_ionization_rate_cgs  : 0               # HI ionization   rate in cgs [1/s]
    RT_HeI_ionization_rate_cgs : 0               # HeI ionization  rate in cgs [1/s]
    RT_HeII_ionization_rate_cgs: 0               # HeII ionization rate in cgs [1/s]
    RT_H2_dissociation_rate_cgs: 0               # H2 dissociation rate in cgs [1/s]
    volumetric_heating_rates_cgs: 0              # Volumetric heating rate in cgs  [erg/s/cm3]
    specific_heating_rates_cgs: 0                # Specific heating rate in cgs    [erg/s/g]
    H2_three_body_rate : 1                       # Specific the H2 formation three body rate (0->5,see Grackle documentation)
    H2_cie_cooling : 0                           # Enable/disable H2 collision-induced emission cooling from Ripamonti & Abel (2004)
    cmb_temperature_floor : 1                    # Enable/disable an effective CMB temperature floor
  
.. note::
   A simple example running SWIFT with Grackle can be find in ``examples/Cooling/CoolingBox``. A more advanced example combining heating and cooling (with heating and ionization sources) is given in ``examples/Cooling/CoolingHeatingBox``.


.. _gear_star_formation:

Star formation
~~~~~~~~~~~~~~

The star formation is done in two steps: first we check if a particle is in the star forming regime and then we use a stochastic approach to transform the gas particles into stars.

A particle is in the star forming regime if:
 - The velocity divergence is negative (:math:`\nabla\cdot v < 0`),
 - The temperature is lower than a threshold (:math:`T < T_t` where :math:`T_t` is defined with ``GEARStarFormation:maximal_temperature``),
 - The gas density is higher than a threshold (:math:`\rho > \rho_t` where :math:`\rho_t` is defined with ``GEARStarFormation:density_threshold``)
 - The particle reaches the pressure floor (:math:`\rho > \frac{\pi}{4 G N_\textrm{Jeans}^{2/3} h^2}\frac{\gamma k_B T}{\mu m_p}` where :math:`N_\textrm{Jeans}` is defined in the pressure floor).

If ``GEARStarFormation:star_formation_mode`` is set to ``agora``, the condition on the pressure floor is ignored. Its default value is ``default``.

A star will be able to form if a randomly drawn number is below :math:`\frac{m_g}{m_\star}\left(1 - \exp\left(-c_\star \Delta t / t_\textrm{ff}\right)\right)` where :math:`t_\textrm{ff}` is the free fall time, :math:`\Delta t` is the time step of the particle and :math:`c_\star` is the star formation coefficient (``GEARStarFormation:star_formation_efficiency``), :math:`m_g` the mass of the gas particle and :math:`m_\star` the mass of the possible future star. The mass of the star is computed from the average gas mass in the initial conditions divided by the number of possible stars formed per gas particle (``GEARStarFormation:n_stars_per_particle``). When we cannot have enough mass to form a second star (defined with the fraction of mass ``GEARStarFormation:min_mass_frac``), we fully convert the gas particle into a stellar particle. Once the star is formed, we move it a bit in a random direction and fraction of the smoothing length in order to avoid any division by 0.

Currently, only the following hydro schemes are compatible: SPHENIX and Gadget2.
Implementing the other hydro schemes is not complicated but requires some careful thinking about the cosmological terms in the definition of the velocity divergence (comoving vs non comoving coordinates and if the Hubble flow is included or not).

.. code:: YAML

  GEARStarFormation:
    star_formation_efficiency: 0.01   # star formation efficiency (c_*)
    maximal_temperature:  3e4         # Upper limit to the temperature of a star forming particle
    n_stars_per_particle: 4           # Number of stars that an hydro particle can generate
    min_mass_frac: 0.5                # Minimal mass for a stellar particle as a fraction of the average mass for the stellar particles.

Sink particles
~~~~~~~~~~~~~~

GEAR now implements sink particles for star formation. Instead of stochastically transforming gas particles into stars as is done in the star formation scheme above when some criteria are met, we transform a gas into a sink particle. The main property of the sink particle is its accretion radius. When gas particles within this accretion radius are eligible to be swallowed by the sink, we remove them and transfer their mass, momentum, angular momentum, chemistry properties, etc to the sink particle.

With the sink particles, the IMF splits into two parts: the continuous part and the discrete part. Those parts will correspond to two kinds of stars. Particles in the discrete part of the IMF represent individual stars. It means that discrete IMF-sampled stars have different masses. Particles in the continuous part represent a population of stars, all with the same mass.

The sink particle will randomly choose a target mass, accrete gas until it reaches this target mass and finally spawn a star. Then, the sink chooses a new target mass and repeats the same procedure. 

This was a brief overview of the model. More details can be found in :ref:`sink_GEAR_model` pages. Please refer to these for configuration, compilations and details of the model. 


Chemistry
~~~~~~~~~

In the chemistry, we are using the smoothed metallicity scheme that consists in using the SPH to smooth the metallicity of each particle over the neighbors. It is worth to point the fact that we are not exchanging any metals but only smoothing it. The parameter ``GEARChemistry:initial_metallicity`` set the (non smoothed) initial mass fraction of each element for all the particles and ``GEARChemistry:scale_initial_metallicity`` use the feedback table to scale the initial metallicity of each element according the Sun's composition.

.. code:: YAML

   GEARChemistry:
    initial_metallicity: 1         # Initial metallicity of the gas (mass fraction)
    scale_initial_metallicity: 1   # Should we scale the initial metallicity with the solar one?


.. _gear_feedback:
   
Feedback
~~~~~~~~

The feedback is composed of a few different models:
  - The initial mass function (IMF) defines the quantity of each type of stars,
  - The lifetime of a star defines when a star will explode (or simply die),
  - The supernovae of type II (SNII) defines the rates and yields,
  - The supernovae of type Ia (SNIa) defines the rates and yields,
  - The energy injection that defines how to inject the energy / metals into the particles.

Most of the parameters are defined inside a table (``GEARFeedback:yields_table``) but can be override with some parameters in the YAML file.
I will not describe theses parameters more than providing them at the end of this section.
Two different models exist for the supernovae (``GEARFeedback:discrete_yields``).
In the continuous mode, we integrate the quantities over the IMF and then explodes a floating point number of stars (can be below 1 in some cases).
In the discrete mode, we avoid the problem of floating points by rounding the number of supernovae (using a floor and randomly adding a supernovae depending on the fractional part) and then compute the properties for a single star at a time.

Initial mass function
^^^^^^^^^^^^^^^^^^^^^

GEAR is using the IMF model from `Kroupa (2001) <https://ui.adsabs.harvard.edu/abs/2001MNRAS.322..231K/abstract>`_.
We have a difference of 1 in the exponent due to the usage of IMF in mass and not in number.
We also restrict the mass of the stars to be inside :math:`[0.05, 50] M_\odot`.
Here is the default model used, but it can be easily adapted through the initial mass function parameters:

.. math::
  \xi(m) \propto m^{-\alpha_i}\, \textrm{where}\,
  \begin{cases}
   \alpha_0 = 0.3,\, & 0.01 \leq m / M_\odot < 0.08, \\
   \alpha_1 = 1.3,\, & 0.08 \leq m / M_\odot < 0.50, \\
   \alpha_2 = 2.3,\, & 0.50 \leq m / M_\odot < 1.00, \\
   \alpha_3 = 2.3,\, & 1.00 \leq m / M_\odot,
  \end{cases}


Lifetime
^^^^^^^^

The lifetime of a star in GEAR depends only on two parameters: first its mass and then its metallicity.

.. math::
   \log(\tau(m)) = a(Z) \log^2(m) + b(Z) \log(m) + c(Z) \\ \\
   a(Z) = -40.110 Z^2 + 5.509 Z + 0.7824 \\
   b(Z) = 141.929 Z^2 - 15.889 Z - 3.2557 \\
   c(Z) = -261.365 Z^2 + 17.073 Z + 9.8661

where :math:`\tau` is the lifetime in years, :math:`m` is the mass of the star (in solar mass) and Z the metallicity of the star.
The parameters previously given are the default ones, they can be modified in the parameters file.

Supernovae II
^^^^^^^^^^^^^

The supernovae rate is simply given by the number of stars massive enough that end their life at the required time.

.. math::
   \dot{N}_\textrm{SNII}(t) = \int_{M_l}^{M_u} \delta(t - \tau(m)) \frac{\phi(m)}{m} \mathrm{d}m

where :math:`M_l` and :math:`M_u` are the lower and upper mass limits for a star exploding in SNII, :math:`\delta` is the Dirac function and :math:`\phi` is the initial mass function (in mass).

The yields for SNII cannot be written in an analytical form, they depend on a few different tables that are based on the work of `Kobayashi et al. (2000) <https://ui.adsabs.harvard.edu/abs/2000ApJ...539...26K/abstract>`_ and `Tsujimoto et al. (1995) <https://ui.adsabs.harvard.edu/abs/1995MNRAS.277..945T/abstract>`_.

Supernovae Ia
^^^^^^^^^^^^^

The supernovae Ia are a bit more complicated as they involve two different stars.

.. math::
  \dot{N}_\textrm{SNIa}(t) = \left( \int_{M_{p,l}}^{M_{p,u}} \frac{\phi(m)}{m} \mathrm{d}m \right) \sum_i b_i \int_{M_{d,l,i}}^{M_{d,u,i}}
  \delta(t-\tau(m)) \frac{\phi_d(m)}{m}\mathrm{d}m

.. math::
   \phi_d(m) \propto m^{-0.35}

where :math:`M_{p,l}` and :math:`M_{p,u}` are the mass limits for a progenitor of a white dwarf, :math:`b_i` is the probability to have a companion and
:math:`M_{d,l,i}` and :math:`M_{d,u,i}` are the mass limits for each type of companion.
The first parenthesis represents the number of white dwarfs and the second one the probability to form a binary.

+------------------+--------------------+-------------------+------------------+
| Companion        |  :math:`M_{d,l,i}` | :math:`M_{d,u,i}` | :math:`b_i`      |
+==================+====================+===================+==================+
| Red giant        |   0.9              |    1.5            |    0.02          |
+------------------+--------------------+-------------------+------------------+
| Main sequence    |   1.8              |    2.5            |    0.05          |
+------------------+--------------------+-------------------+------------------+

The yields are based on the same papers than the SNII.

Energy injection
^^^^^^^^^^^^^^^^

All the supernovae (type II and Ia) inject the same amount of energy into the surrounding gas (``GEARFeedback:supernovae_energy_erg``) and distribute it according to the hydro kernel.
The same is done with the metals and the mass.


Generating a new table
^^^^^^^^^^^^^^^^^^^^^^

The feedback table is an HDF5 file with the following structure:

.. graphviz:: feedback_table.dot

where the solid (dashed) squares represent a group (a dataset) with the name of the object underlined and the attributes written below. Everything is in solar mass or without units (e.g. mass fraction or unitless constant).
In ``Data``, the attribute ``elts`` is an array of string with the element names (the last should be ``Metals``, it corresponds to the sum of all the elements), ``MeanWDMass`` is the mass of the white dwarfs
and ``SolarMassAbundances`` is an array of float containing the mass fraction of the different element in the sun.
In ``IMF``, ``n + 1`` is the number of part in the IMF, ``as`` are the exponent (``n+1`` elements), ``ms`` are the mass limits between each part (``n`` elements) and
``Mmin`` (``Mmax``) is the minimal (maximal) mass of a star.
In ``LifeTimes``, the coefficient are given in the form of a single table (``coeff_z`` with a 3x3 shape).
In ``SNIa``, ``a`` is the exponent of the distribution of binaries, ``bb1``  and ``bb2`` are the coefficient :math:`b_i` and the other attributes follow the same names than in the SNIa formulas.
The ``Metals`` group from the ``SNIa`` contains the name of each elements (``elts``) and the metal mass fraction ejected by each supernovae (``data``) in the same order. They must contain the same elements than in ``Data``.
Finally for the ``SNII``, the mass limits are given by ``Mmin`` and ``Mmax``. For the yields, the datasets required are ``Ej`` (mass fraction ejected [processed]), ``Ejnp`` (mass fraction ejected [non processed]) and one dataset for each element present in ``elts``. The datasets should all have the same size, be uniformly sampled in log and contains the attributes ``min`` (mass in log for the first element) and ``step`` (difference of mass in log between two elements).

.. code:: YAML

  GEARFeedback:
    supernovae_energy_erg: 0.1e51                            # Energy released by a single supernovae.
    yields_table: chemistry-AGB+OMgSFeZnSrYBaEu-16072013.h5  # Table containing the yields.
    discrete_yields: 0                                       # Should we use discrete yields or the IMF integrated one?
  GEARInitialMassFunction:
    number_function_part:  4                       # Number of different part in the IMF
    exponents:  [0.7, -0.8, -1.7, -1.3]            # Exponents of each part of the IMF
    mass_limits_msun:  [0.05, 0.08, 0.5, 1, 50]    # Limits in mass between each part of the IMF
  GEARLifetime:
   quadratic:  [-40.1107, 5.50992, 0.782432]  # Quadratic terms in the fit
   linear:  [141.93, -15.8895, -3.25578]      # Linear terms in the fit
   constant:  [-261.366, 17.0735, 9.86606]    # Constant terms in the fit
  GEARSupernovaeIa:
    exponent:  -0.35                      # Exponent for the distribution of companions
    min_mass_white_dwarf_progenitor:  3   # Minimal mass of a progenitor of white dwarf
    max_mass_white_dwarf_progenitor:  8   # Maximal mass of a progenitor of white dwarf
    max_mass_red_giant:  1.5              # Maximal mass for a red giant
    min_mass_red_giant:  0.9              # Minimal mass for a red giant
    coef_red_giant:  0.02                 # Coefficient for the distribution of red giants companions
    max_mass_main_sequence:  2.6          # Maximal mass for a main sequence star
    min_mass_main_sequence:  1.8          # Minimal mass for a main sequence star
    coef_main_sequence:  0.05             # Coefficient for the distribution of main sequence companions
    white_dwarf_mass:  1.38               # Mass of a white dwarf
  GEARSupernovaeII:
  interpolation_size:  200                # Number of elements for the interpolation of the data
