.. GEAR sub-grid model feedback and stellar evolution
   Loic Hausammann, 17th April 2020
   Darwin Roduit, 20th January 2026
   Rey Zachary, 18th February 2026

.. _gear_stellar_evolution_and_feedback:  

Stellar evolution and feedback
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The feedback prescription is composed of a few different models:
  - The initial mass function (IMF) defines the quantity of each type of star,
  - The lifetime of a star defines when a star will explode (or die),
  - The supernovae of type II (SNII) define the rates and yields,
  - The supernovae of type Ia (SNIa) define the rates and yields,
  - The energy injection that defines how to inject the energy / metals into the particles.
  - The Stellar Winds (SW) defines the energy and mass continuously ejected by stars until their death. 

Most of the parameters are defined inside a table (``GEARFeedback:yields_table``). To generate the table, we use `pychem <https://www.astro.unige.ch/~revazy/PyChem/>`_ python module. You can get a table by clicking on `this link <https://virgodb.cosma.dur.ac.uk/swift-webstorage/FeedbackTables/POPIIsw.h5>`_. Some examples in ``swiftsim/examples/`` use this table, e.g. ``swiftsim/examples/GEAR``, ``swiftsim/examples/IsolatedGalaxy/IsolatedGalaxy_multi_component/GEAR/``


Stellar evolution
-----------------

Initial mass function
^^^^^^^^^^^^^^^^^^^^^

GEAR is using the IMF model from `Kroupa (2001) <https://ui.adsabs.harvard.edu/abs/2001MNRAS.322..231K/abstract>`_.
We also restrict the mass of the stars to be inside :math:`[\mathrm{Mmin}, \mathrm{Mmax}] M_\odot`, where ``Mmin`` and ``Mmax`` are read from the table. In the default model, the IMF :math:`\xi(m)` that represents the *number of stars*  found in the mass range :math:`[m, m + \mathrm{d}m]` is (`Revaz (2012) <https://ui.adsabs.harvard.edu/abs/2012A%26A...538A..82R/abstract>`_):

.. math::
  \xi(m) \propto m^{\alpha},
  
where,
  
.. math::
  \alpha =  
  \begin{cases}
    -0.3,\, & 0.01 \leq m / M_\odot < 0.08, \\
    -1.8,\, & 0.08 \leq m / M_\odot < 0.50, \\
    -2.7,\, & 0.50 \leq m / M_\odot < 1.00, \\
    -2.3,\, & 1.00 \leq m / M_\odot.
  \end{cases}

Note that those values differs by -1 with respect to the ones provided in the GEAR tables, where the IMF :math:`\phi(m)` is defined as the *fraction of stellar mass* found in the mass range :math:`[m, m + \mathrm{d}m]`:

.. math::
  \phi(m) = m\,\xi(m) \propto m^{\alpha+1}.




Lifetime
^^^^^^^^

The lifetime of a star in GEAR depends only on two parameters: first its mass and then its metallicity. In GEAR, we use the following approximation `(Poirier, 2004) <https://theses.fr/2004STR13003>`_ and summarised in `Hausammann (2021) <https://infoscience.epfl.ch/entities/publication/3e6d2e54-a782-440a-86c3-05482e83794d>`_:

.. math::
   \log(\tau(m)) = a(Z) \log^2(m) + b(Z) \log(m) + c(Z)

where :math:`\tau` is the lifetime in years, :math:`m` is the mass of the star (in solar mass) and :math:`Z` is the metallicity of the star. The default table we provided uses the following values:

.. math::
   a(Z) = -40.110 Z^2 + 5.509 Z + 0.7824 \\
   b(Z) = 141.929 Z^2 - 15.889 Z - 3.2557 \\
   c(Z) = -261.365 Z^2 + 17.073 Z + 9.8661 \, .

Supernovae II
^^^^^^^^^^^^^

The supernovae rate is given by the number of stars massive enough that end their life at the required time:

.. math::
   \dot{N}_\textrm{SNII}(t) = \int_{M_l}^{M_u} \delta(t - \tau(m)) \frac{\phi(m)}{m} \mathrm{d}m \, ,

where :math:`M_l` and :math:`M_u` are the lower and upper mass limits for a star exploding in SNII, :math:`\delta` is the Dirac function and :math:`\phi` is the initial mass function (in mass).

The yields for SNII cannot be written in an analytical form; they depend on a few different tables that are based on the work of `Kobayashi et al. (2000) <https://ui.adsabs.harvard.edu/abs/2000ApJ...539...26K/abstract>`_ and `Tsujimoto et al. (1995) <https://ui.adsabs.harvard.edu/abs/1995MNRAS.277..945T/abstract>`_.

Supernovae Ia
^^^^^^^^^^^^^

The supernovae Ia are more complicated as they involve two different stars. The SNIa rate is:

.. math::
  \dot{N}_\textrm{SNIa}(t) = \left( \int_{M_{p,l}}^{M_{p,u}} \frac{\phi(m)}{m} \mathrm{d}m \right) \sum_i b_i \int_{M_{d,l,i}}^{M_{d,u,i}} \, ,
  \delta(t-\tau(m)) \frac{\phi_d(m)}{m}\mathrm{d}m \, ,

.. math::
   \phi_d(m) \propto m^{-0.35}

where :math:`M_{p,l}` and :math:`M_{p,u}` are the mass limits for a progenitor of a white dwarf, :math:`b_i` is the probability to have a companion and
:math:`M_{d,l,i}` and :math:`M_{d,u,i}` are the mass limits for each type of companion.
The first integral represents the number of white dwarfs progenitor, while the second part accounts the probability of forming a binary. The default parameters used for the SNIa rate are given in the following table. They can be set to different values via PyChem:

+------------------+--------------------+-------------------+------------------+
| Companion        |  :math:`M_{d,l,i}` | :math:`M_{d,u,i}` | :math:`b_i`      |
+==================+====================+===================+==================+
| Red giant        |   0.9              |    1.5            |    0.02          |
+------------------+--------------------+-------------------+------------------+
| Main sequence    |   1.8              |    2.5            |    0.05          |
+------------------+--------------------+-------------------+------------------+

The yields are based on the same papers as the SNII.

Stellar winds
^^^^^^^^^^^^^

The stellar wind model used in GEAR is based on the work of `Deng et al. (2024b) <https://arxiv.org/abs/2405.08869>`_.

The energy and mass ejected by a star (single type of stellar particle) in GEAR depends only on two parameters: its initial mass (:math:`M_{init}`) and its metallicity (:math:`Z`).
A set of power law allows to calculate two quantity, the mass ejected (:math:`\dot{M}`) and the terminal wind velocity (:math:`v_\infty`):

.. math::
  \log_{10}\mathcal{A} = \left\{ \begin{array}{rcl} b_0 + b_1 x + b_2 x^2 + b_3 x^3 + b_4 x^4\ , & x\leq x_0,\\
  \mathcal{A}_0 + \frac{d \log_{10}\mathcal{A}}{dx}\bigg|_{x=x_0}(x - x_0)\ , & x > x_0,
  \end{array} \right.

.. math::
  \mathcal{A}_0 = b_0 + b_1 x_0 + b_2 x^2_0 + b_3 x^3_0 + b_4 x^4_0,

.. math::
  \frac{d\log_{10}\mathcal{A}}{dx}\bigg|_{x=x_0} = (b_1 + 2b_2 x_0 + 3b_3 x^2_0 + 4b_4 x^3_0),

.. math::
  b_i = a_{i0} + a_{i1} y + a_{i2} y^2

with :math:`\mathcal{A}(M,Z)` being either :math:`\dot{M}` or/and :math:`v_\infty`, :math:`y\equiv \log_{10}(Z)`, :math:`x\equiv \log_{10}(M_{init})` where :math:`x_0` is a mass limit above which the behavior is better fitted by linear relation.
Then, the power ejected is simply:

.. math::
  \dot{E} = L = \frac{1}{2}\dot{M}v_\infty

The power and mass ejected by a stellar particle owning an entire or partial IMF is the sum of all the stars contained in the stellar particle, or more precisely, the integral over the part of the IMF that is not dead yet:

.. math::
  Q_{\star} = M_{\star}\int^{M_u}_{M_{min}}\frac{\phi(m)}{m}Q(m)dm

Where :math:`Q(m)` is the wanted quantity (either :math:`\dot{M}` or :math:`L`) for a single star of mass :math:`m`, :math:`Q_{\star}` the wanted quantity of the stellar particle, :math:`M_u` the upper mass limits for a star dying (or exploding), :math:`M_{min}` the minimal mass of the IMF, :math:`\delta` is the Dirac function and :math:`\phi` is the initial mass function (in mass).

Stellar evolution table
^^^^^^^^^^^^^^^^^^^^^^^

To generate the table, we use `pychem <https://www.astro.unige.ch/~revazy/PyChem/>`_ python module. Below, we provide an overview of the structure of the table. Please refer to the pychem documentation for more information.

.. graphviz:: feedback_table.dot

where the solid (dashed) squares represent a group (a dataset) with the name of the object underlined and the attributes written below. Everything is in solar mass or without units (e.g. mass fraction or unitless constant).

In ``Data``, the attribute ``elts`` is an array of string with the element names (the last should be ``Metals``, it corresponds to the sum of all the elements), ``MeanWDMass`` is the mass of the white dwarfs and ``SolarMassAbundances`` is an array of float containing the mass fraction of the different element in the sun.
simply
In ``IMF``, ``n + 1`` is the number of parts in the IMF, ``as`` are the exponent (``n+1`` elements), ``ms`` are the mass limits between each part (``n`` elements) and ``Mmin`` (``Mmax``) is the minimal (maximal) mass of a star.

In ``LifeTimes``, the coefficients are given in the form of a single table (``coeff_z`` with a 3x3 shape).

In ``SNIa``, ``a`` is the exponent of the distribution of binaries, ``bb1``  and ``bb2`` are the coefficients :math:`b_i` and the other attributes follow the same names as in the SNIa formulas.

The ``Metals`` group from the ``SNIa`` contains the name of each element (``elts``) and the metal mass fraction ejected by each supernova (``data``) in the same order. They must contain the same elements as in ``Data``.

Then, for the ``SNII``, the mass limits are given by ``Mmin`` and ``Mmax``. For the yields, the datasets required are ``Ej`` (mass fraction ejected [processed]), ``Ejnp`` (mass fraction ejected [non processed]) and one dataset for each element present in ``elts``. The datasets should all have the same size, be uniformly sampled in log and contain the attributes ``min`` (mass in log for the first element) and ``step`` (difference of mass in log between two elements).

Finally, in the ``SW`` group, under the subgroup ``MetallicityDependent``, there are 4 2D datasets: ``Energy`` and ``Mass`` for the power and mass ejected by the stellar winds for a single star, and ``Integrated_Energy`` and ``Integrated_Mass_Loss`` for the power and mass ejected by a continuous stellar particle. Each of these datasets contains 8 attributes: ``label`` (the label of the dataset), ``dims`` (the dimensions of the grid), ``m0`` and ``z0`` (the the mass and metallicity in log of the first element of the grid), ``dm`` and ``dz`` (the logarithmic spacing in mass and metallcity between two element of the grid), and ``nm`` and ``nz`` (the number of elements in mass and metallicity in the grid).

GEAR includes two types of tables, one for population II stars and one for population III. The tables are specified by ``GEARFeedback:yields_table`` and ``GEARFeedback:yields_table_first_stars``. The choice of the table depends on the metallicity [Fe/H] (``GEARFeedback:imf_transition_metallicity``). Below this metallicity, we use ``yields_table_first_stars`` ; above, we use ``yields_tables``. If we set ``imf_transition_metallicity`` to 0, we only use ``yields_tables``.

Star particles types
^^^^^^^^^^^^^^^^^^^^

In GEAR, we consider three types of star particles:

* **Single stellar population (SSP)** star particles: These particles represent a population of stars defined by the provided IMF. These particles are created by the star formation scheme.

* **Individual stars**: These particles represent individual stars. These particles are created by the sink particles scheme (see :ref:`sink_GEAR_model_summary`).

* **Continuous stars**: These particles represent the integrated low mass portion of the IMF (e.g. stars with :math:`M_\star < 8 \, M_\odot`), while the remaining mass is sampled as individual stars (e.g. stars with :math:`M_\star \geq 8 \, M_\odot`). These particles are created by the sink particles scheme (see :ref:`sink_GEAR_model_summary`). Note that by a careful choice of parameters, these stars may represent the full IMF and thus be treated like SSP stars, with the difference being that they are spawned by sink particles and not gas.

The stellar evolution is treated as follows for each star type:

- **Individual stars**: There is no IMF sampling or averaging needed. We know the star's mass and metallicity. Therefore, its stellar evolution properties are known, e.g. it will explode into exactly one SN II. There is no SNIa.

- **Single star population (SSP) and Continuous stars**: They are treated in the same manner, with the only difference being the IMF upper limit. For the SSP stars, this is ``Mmax`` defined in the stellar evolution tables; for the continuous stars, this is defined by ``GEARSink:minimal_discrete_mass_Msun`` for population II stars and ``GEARSink:minimal_discrete_mass_first_stars_Msun`` for population III stars. For these particles, we need to sample the IMF. This is explained in the next section.

IMF sampling
^^^^^^^^^^^^

To properly determine the number of supernovae and the mass of the ejected yields, we need to sample the IMF for SSP and continuous star particles. GEAR implements two methods `(Revaz et al. 2016) <https://ui.adsabs.harvard.edu/abs/2016A&A...588A..21R>`_:

- Continuous IMF sampling (CIMFS) (``GEARFeedback:discrete_yields: 0``): We integrate the quantities over the IMF and then explode a floating-point number of stars, which can be below 1 in some cases. This method works well for large stellar particle mass and supernovae rates, but not for low stellar particle mass or low supernovae rates. Note that SNIa often occur in the second regime, hence the method. The overall effect is similar to diluting the SN explosions over time. 


- Random discrete IMF sampling (RIMFS) (``GEARFeedback:discrete_yields: 1``): We avoid the issue of non-integer event numbers by taking the floor of the calculated SN count and stochastically adding an additional supernova based on the fractional part. We then compute the properties for a single stellar at a time.

Feedback energy, metals and momentum injection
----------------------------------------------

When a star goes into a supernova (type II and Ia), the stellar evolution determines how much energy, mass and metals are released during the explosion. For the SNII, the tables contain the value of the energy released for different stars' masses (usually :math:`10^{51}` erg, but can vary by model). For the SNIa, the energy is specified by the parameter ``GEARFeedback:supernovae_Ia_energy_erg``. In both cases, these are the energies *before* any supernova efficiency is applied.

The energy can be distributed as internal/thermal energy or as momentum. Thus, we need to distribute internal energy, momentum, mass and metals to the gas particles. We will group all these in the “fluxes” term. We have two models for the distribution of these fluxes and the subgrid modelling of the supernovae: GEAR model and GEAR mechanical model. We describe the two schemes in the  :ref:`gear_sn_feedback_models` page.


Model parameters
----------------

Here, we provide a summary of the feedback parameters contained in the YAML parameter file.

The first two parameters relate to the quantity of energy injected and are available for all feedback prescriptions.

* The energy released by a supernova: ``supernovae_Ia_energy_erg``. Note that the SNII energy is stored in the feedback tables ``yields_table_first_stars`` and ``yields_table``.

* The effective injected energy is multiplied by a factor ``supernovae_efficiency``. This factor is applied to SNII and SNIa.
  
* ``yields_table``: Table containing the stellar evolution parameters for population II stars, including the yields.
  
* ``yields_table_first_stars``: Similar to ``yields_table`` but for the first stars (population III).
  
* ``imf_transition_metallicity`` specifies which table to use based on [Fe/H]. If the gas metallicity is below ``imf_transition_metallicity``, we use the population III table; if it is above, we use the population table II.
  
* Two different models exist for the IMF sampling of supernovae (``GEARFeedback:discrete_yields``): the continuous mode (``GEARFeedback:discrete_yields: 1``) and the discrete mode (``GEARFeedback:discrete_yields: 0``).
  
* ``elements`` is the list of yields to read from the tables. The number of elements is specified at compile time  ``(--with-chemistry=GEAR_N)``

* ``discrete_star_minimal_gravity_mass_Msun``: Minimal gravity mass in solar masses after a discrete star completely explodes. Default: 0.1

* ``GEARSupernovaeII:interpolation_size`` is the number of elements to keep in the interpolation of the data.

* ``GEARStellar_wind:interpolation_size_mass``  Size of the mass array of the grid used in stellar winds yields.

* ``GEARStellar_wind:interpolation_size_metallicity`` Size of the metallicity array of the grid used in stellar winds yields.

Here is the whole feedback section:

.. code:: YAML

	  GEARFeedback:
	    supernovae_Ia_energy_erg: 1e51                           # Energy released by a single supernova.
	    supernovae_efficiency: 0.1                               # Supernovae energy efficiency, used for both SNIa and SNII. The energy released effectively is E_sn = supernovae_efficiency*E_sn
	    yields_table: chemistry-AGB+OMgSFeZnSrYBaEu-16072013.h5  # Table containing the yields.
	    yields_table_first_stars: chemistry-PopIII.hdf5          # Table containing the yields of the first stars (population III).
	    imf_transition_metallicity: -5                           # Maximal metallicity ([Fe/H]) for a first star (0 to deactivate).
	    discrete_yields: 0                                       # Should we use discrete yields or the IMF integrated one?
	    elements: [Fe, Mg, O, S, Zn, Sr, Y, Ba, Eu]              # Elements to read in the yields table. The number of elements should be one less than the number of elements (N) requested during the configuration (--with-chemistry=GEAR_N).
	    discrete_star_minimal_gravity_mass_Msun: 0.1             # Minimal gravity mass after a discrete star completely explodes. In M_sun. (Default: 0.1)

	  GEARSupernovaeII:
	    interpolation_size:  200                                 # Number of elements to keep in the interpolation of the data. (Default: 200)

    GEARStellar_wind:
      interpolation_size_mass:        200                      # Size of the mass array of the grid used in stellar winds yields.
      interpolation_size_metallicity: 110                      # Size of the metallicity array of the grid used in stellar winds yields.

References
----------

- `pychem <https://www.astro.unige.ch/~revazy/PyChem/>`_ : python module to generate GEAR tables

- `Tsujimoto et al. (1995) <https://ui.adsabs.harvard.edu/abs/1995MNRAS.277..945T/abstract>`_.

- `Kobayashi et al. (2000) <https://ui.adsabs.harvard.edu/abs/2000ApJ...539...26K/abstract>`_

- `Kroupa (2001) <https://ui.adsabs.harvard.edu/abs/2001MNRAS.322..231K/abstract>`_

- `Poirier (2004) <https://theses.fr/2004STR13003>`_

- `Revaz et al. (2016) <https://ui.adsabs.harvard.edu/abs/2016A&A...588A..21R>`_

- `Hausammann (2021) <https://infoscience.epfl.ch/entities/publication/3e6d2e54-a782-440a-86c3-05482e83794d>`_

- `Deng et al. (2024b) <https://arxiv.org/abs/2405.08869>`_
