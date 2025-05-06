.. GEAR sub-grid model feedback and stellar evolution
   Loic Hausammann, 17th April 2020
   Darwin Roduit, 30th March 2025

.. _gear_stellar_evolution_and_feedback:  

Stellar evolution and feedback
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

Supernova energy injection
^^^^^^^^^^^^^^^^^^^^^^^^^^

When a star goes into a supernova (type II and Ia), the stellar evolution determines how much energy (``GEARFeedback:supernovae_energy_erg`` TO BE CHECKED), mass and metals are released during the explosion. The energy can be ditributed as internal/thermal energy or as momentum. Thus, we need to distribute internal energy, momentum, mass and metals to the gas particles.  We will group all these in the “fluxes” term. 

We have two models for the distribution of these fluxes and the subgrid modelling of the supernovae : GEAR model and GEAR mechanical model. We describe the two schemes in the  :ref:`gear_sn_feedback_models` page.


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
