.. GEAR sub-grid model
   Loic Hausammann, 17th April 2020
   Darwin Roduit, 30th March 2025

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


.. _gear_star_formation:


Star formation
~~~~~~~~~~~~~~

The star formation is done in two steps: first we check if a particle is in the star forming regime and then we use a stochastic approach to transform the gas particles into stars.

A particle is in the star forming regime if:
 - The velocity divergence is negative (:math:`\nabla\cdot v < 0`),
 - The temperature is lower than a threshold (:math:`T < T_t` where :math:`T_t` is defined with ``GEARStarFormation:maximal_temperature_K``),
 - The gas density is higher than a threshold (:math:`\rho > \rho_t` where :math:`\rho_t` is defined with ``GEARStarFormation:density_threshold_Hpcm3``)
 - The particle reaches the pressure floor (:math:`\rho > \frac{\pi}{4 G N_\textrm{Jeans}^{2/3} h^2}\frac{\gamma k_B T}{\mu m_p}` where :math:`N_\textrm{Jeans}` is defined in the pressure floor).

If ``GEARStarFormation:star_formation_mode`` is set to ``agora``, the condition on the pressure floor is ignored. Its default value is ``default``.

A star will be able to form if a randomly drawn number is below :math:`\frac{m_g}{m_\star}\left(1 - \exp\left(-c_\star \Delta t / t_\textrm{ff}\right)\right)` where :math:`t_\textrm{ff}` is the free fall time, :math:`\Delta t` is the time step of the particle and :math:`c_\star` is the star formation coefficient (``GEARStarFormation:star_formation_efficiency``), :math:`m_g` the mass of the gas particle and :math:`m_\star` the mass of the possible future star. The mass of the star is computed from the average gas mass in the initial conditions divided by the number of possible stars formed per gas particle (``GEARStarFormation:n_stars_per_particle``). When we cannot have enough mass to form a second star (defined with the fraction of mass ``GEARStarFormation:min_mass_frac``), we fully convert the gas particle into a stellar particle. Once the star is formed, we move it a bit in a random direction and fraction of the smoothing length in order to avoid any division by 0.

Currently, only the following hydro schemes are compatible: SPHENIX, Gadget2, minimal SPH, Gasoline-2 and Pressure-Energy.
Implementing the other hydro schemes is not complicated but requires some careful thinking about the cosmological terms in the definition of the velocity divergence (comoving vs non comoving coordinates and if the Hubble flow is included or not).

.. code:: YAML

  GEARStarFormation:
    star_formation_efficiency: 0.01   # star formation efficiency (c_*)
    maximal_temperature_K:  3e4       # Upper limit to the temperature of a star forming particle
    density_threshold_Hpcm3:   10     # Density threshold (Hydrogen atoms/cm^3) for star formation
    n_stars_per_particle: 4           # Number of stars that an hydro particle can generate
    min_mass_frac: 0.5                # Minimal mass for a stellar particle as a fraction of the average mass for the stellar particles.

Initial Conditions
++++++++++++++++++

Note that if in the initial conditions, the time of formation of a stellar particle is given (``BirthTime``)
and set to a negative value, the stellar particle will provide no feedback.
A similar behavior will be obtained if the parameter ``Stars:overwrite_birth_time`` is set to 1 and
``Stars:birth_time`` to -1.


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

Configuring and compiling SWIFT with Grackle
++++++++++++++++++++++++++++++++++++++++++++

In order to compile SWIFT with Grackle, you need to provide the options ``with-chemistry=GEAR`` and ``with-grackle=$GRACKLE_ROOT``
where ``$GRACKLE_ROOT`` is the root of the install directory (not the ``lib``). 

.. warning::
    (State 2023) Grackle is experiencing current development, and the API is subject
    to changes in the future. For convenience, a frozen version is hosted as a fork
    on github here: https://github.com/mladenivkovic/grackle-swift .
    The version available there will be tried and tested and ensured to work with
    SWIFT.

    Additionally, that repository hosts files necessary to install that specific 
    version of grackle with spack.

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

Parameters
++++++++++

When starting a simulation without providing the different element fractions in the non equilibrium mode, the code supposes an equilibrium and computes them automatically.
The code uses an iterative method in order to find the correct initial composition and this method can be tuned with two parameters. ``GrackleCooling:max_steps`` defines the maximal number of steps to reach the convergence and ``GrackleCooling:convergence_limit`` defines the tolerance in the relative error.

In the parameters file, a few different parameters are available.

- ``GrackleCooling:redshift`` defines the redshift to use for the UV background (for cosmological simulation, it must be set to -1 in order to use the simulation's redshift).

- ``GrackleCooling:provide_*_heating_rates`` can enable the computation of user provided heating rates (such as with the radiative transfer) in either volumetric or specific units.

- Feedback can be made more efficient by turning off the cooling during a few Myr (``GrackleCooling:thermal_time_myr``) for the particles touched by a supernovae.

- The self shielding method is defined by ``GrackleCooling:self_shielding_method`` where 0 means no self shielding, > 0 means a method defined in Grackle (see Grackle documentation for more information) and -1 means GEAR's self shielding that simply turn off the UV background when reaching a given density (``GrackleCooling:self_shielding_threshold_atom_per_cm3``).

- The initial elemental abundances can be specified with ``initial_nX_to_nY_ratio``, e.g. if you want to specify an initial HII to H abundance. A negative value ignores the parameter. The complete list can be found below.

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


Complete parameter list
-----------------------

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
    maximal_density_Hpcm3:   1e4                 # Maximal density (in hydrogen atoms/cm^3) for cooling. Higher densities are floored to this value to ensure grackle works properly when interpolating beyond the cloudy_table maximal density. A value < 0 deactivates this parameter.
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
    H2_on_dust: 0                                # Flag to enable H2 formation on dust grains
    local_dust_to_gas_ratio : -1                 # The ratio of total dust mass to gas mass in the local Universe (-1 to use the Grackle default value). 
    cmb_temperature_floor : 1                    # Enable/disable an effective CMB temperature floor

    initial_nHII_to_nH_ratio:    -1              # initial nHII   to nH ratio (number density ratio). Value is ignored if set to -1.
    initial_nHeI_to_nH_ratio:    -1              # initial nHeI   to nH ratio (number density ratio). Value is ignored if set to -1.
    initial_nHeII_to_nH_ratio:   -1              # initial nHeII  to nH ratio (number density ratio). Value is ignored if set to -1.
    initial_nHeIII_to_nH_ratio:  -1              # initial nHeIII to nH ratio (number density ratio). Value is ignored if set to -1.
    initial_nDI_to_nH_ratio:     -1              # initial nDI    to nH ratio (number density ratio). Value is ignored if set to -1.
    initial_nDII_to_nH_ratio:    -1              # initial nDII   to nH ratio (number density ratio). Value is ignored if set to -1.
    initial_nHM_to_nH_ratio:     -1              # initial nHM    to nH ratio (number density ratio). Value is ignored if set to -1.
    initial_nH2I_to_nH_ratio:    -1              # initial nH2I   to nH ratio (number density ratio). Value is ignored if set to -1.
    initial_nH2II_to_nH_ratio:   -1              # initial nH2II  to nH ratio (number density ratio). Value is ignored if set to -1.
    initial_nHDI_to_nH_ratio:    -1              # initial nHDI   to nH ratio (number density ratio). Value is ignored if set to -1.

.. note::
   A simple example running SWIFT with Grackle can be find in ``examples/Cooling/CoolingBox``. A more advanced example combining heating and cooling (with heating and ionization sources) is given in ``examples/Cooling/CoolingHeatingBox``. ``examples/Cooling/CoolingWithPrimordialElements/`` runs a uniform cosmological box with imposed abundances and let them evolve down to redshift 0.
