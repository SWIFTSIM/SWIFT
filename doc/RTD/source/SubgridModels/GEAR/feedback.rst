.. GEAR sub-grid model feedback and stellar evolution
   Loic Hausammann, 17th April 2020
   Darwin Roduit, 20th January 2026

.. _gear_stellar_evolution_and_feedback:  

Stellar evolution and feedback
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The feedback prescription is composed of a few different models:
  - The initial mass function (IMF) defines the quantity of each type of star,
  - The lifetime of a star defines when a star will explode (or die),
  - The supernovae of type II (SNII) define the rates and yields,
  - The supernovae of type Ia (SNIa) define the rates and yields,
  - The energy injection that defines how to inject the energy / metals into the particles.

Most of the parameters are defined inside a table (``GEARFeedback:yields_table``). To generate the table, we use `pychem <https://www.astro.unige.ch/~revazy/PyChem/>`_ python module. You can get a table by clicking on `this link <https://virgodb.cosma.dur.ac.uk/swift-webstorage/FeedbackTables/POPIIsw.h5>`_. Some examples in ``swiftsim/examples/`` use this table, e.g. ``swiftsim/examples/GEAR``, ``swiftsim/examples/IsolatedGalaxy/IsolatedGalaxy_multi_component/GEAR/``


Stellar evolution
-----------------

Initial mass function
^^^^^^^^^^^^^^^^^^^^^

GEAR is using the IMF model from `Kroupa (2001) <https://ui.adsabs.harvard.edu/abs/2001MNRAS.322..231K/abstract>`_.
We also restrict the mass of the stars to be inside :math:`[\mathrm{Mmin}, \mathrm{Mmax}] M_\odot`, where ``Mmin`` and ``Mmax`` are read from the table. In the default model, the number of stars :math:`\xi` found in the mass range :math:`[m, m + \mathrm{d}m]` is:

.. math::
  \xi(m) \propto m^{-\alpha_i}\, \textrm{where}\,
  \begin{cases}
   \alpha_0 = 0.3,\, & 0.01 \leq m / M_\odot < 0.08, \\
   \alpha_1 = 1.3,\, & 0.08 \leq m / M_\odot < 0.50, \\
   \alpha_2 = 2.3,\, & 0.50 \leq m / M_\odot < 1.00, \\
   \alpha_3 = 2.3,\, & 1.00 \leq m / M_\odot,
  \end{cases}

Note that we have a difference of 1 in the exponent due to the use of IMF in mass rather than in number. The IMF is defined in mass through :math:`\phi \propto m^{1 - \alpha_i}`.

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
The first integral represents the number of white dwarfs progenitor, while the second part accounts the probability of forming a binary. The parameters used for the SNIa rate are given in the following table:

+------------------+--------------------+-------------------+------------------+
| Companion        |  :math:`M_{d,l,i}` | :math:`M_{d,u,i}` | :math:`b_i`      |
+==================+====================+===================+==================+
| Red giant        |   0.9              |    1.5            |    0.02          |
+------------------+--------------------+-------------------+------------------+
| Main sequence    |   1.8              |    2.5            |    0.05          |
+------------------+--------------------+-------------------+------------------+

The yields are based on the same papers as the SNII.

Stellar evolution table
^^^^^^^^^^^^^^^^^^^^^^^

To generate the table, we use `pychem <https://www.astro.unige.ch/~revazy/PyChem/>`_ python module. In summary, the feedback table is an HDF5 file with the following structure:

.. graphviz:: feedback_table.dot

where the solid (dashed) squares represent a group (a dataset) with the name of the object underlined and the attributes written below. Everything is in solar mass or without units (e.g. mass fraction or unitless constant).

In ``Data``, the attribute ``elts`` is an array of string with the element names (the last should be ``Metals``, it corresponds to the sum of all the elements), ``MeanWDMass`` is the mass of the white dwarfs and ``SolarMassAbundances`` is an array of float containing the mass fraction of the different element in the sun.
simply
In ``IMF``, ``n + 1`` is the number of parts in the IMF, ``as`` are the exponent (``n+1`` elements), ``ms`` are the mass limits between each part (``n`` elements) and ``Mmin`` (``Mmax``) is the minimal (maximal) mass of a star.

In ``LifeTimes``, the coefficients are given in the form of a single table (``coeff_z`` with a 3x3 shape).

In ``SNIa``, ``a`` is the exponent of the distribution of binaries, ``bb1``  and ``bb2`` are the coefficients :math:`b_i` and the other attributes follow the same names as in the SNIa formulas.

The ``Metals`` group from the ``SNIa`` contains the name of each element (``elts``) and the metal mass fraction ejected by each supernova (``data``) in the same order. They must contain the same elements as in ``Data``.

Finally, for the ``SNII``, the mass limits are given by ``Mmin`` and ``Mmax``. For the yields, the datasets required are ``Ej`` (mass fraction ejected [processed]), ``Ejnp`` (mass fraction ejected [non processed]) and one dataset for each element present in ``elts``. The datasets should all have the same size, be uniformly sampled in log and contain the attributes ``min`` (mass in log for the first element) and ``step`` (difference of mass in log between two elements).

GEAR includes two types of tables, one for population II stars and one for population III. The tables are specified by ``GEARFeedback:yields_table`` and ``GEARFeedback:yields_table_first_stars``. The choice of the table depends on the metallicity [Fe/H] (``GEARFeedback:imf_transition_metallicity``). Below this metallicity, we use ``yields_table_first_stars`` ; above, we use ``yields_tables``. If we set ``imf_transition_metallicity`` to 0, we only use ``yields_tables``.

Star particles types
^^^^^^^^^^^^^^^^^^^^

In GEAR, we consider three types of star particles:

* **Single star population (SSP)** star particles: These particles represent a population of stars defined by the provided IMF. These particles are created by the star formation scheme.

* **Individual stars**: These particles represent individual stars. These particles are created by the sink particles scheme (see :ref:`sink_GEAR_model_summary`).

* **Continuous stars**: These particles represent the integrated portion of the IMF (e.g. stars with :math:`M_\star < 8 \, M_\odot`), while the remaining mass is sampled as individual stars (e.g. stars with :math:`M_\star \geq 8 \, M_\odot`). These particles are created by the sink particles scheme (see :ref:`sink_GEAR_model_summary`). Note that by a careful choice of parameters, these stars may represent the full IMF and thus be treated like SSP stars, with the difference being that they are spawned by sink particles and not gas.

The stellar evolution is treated as follows for each star type:

- *Individual stars*: There is no IMF sampling or averaging needed. We know the star's mass and metallicity. Therefore, its stellar evolution properties are known, e.g. it will explode into exactly one SN II. There is no SNIa.

- *SSP and continuous stars*: They are treated in the same manner, with the only difference being the IMF upper limit. For the SSP stars, this is ``Mmax`` defined in the stellar evolution tables; for the continuous stars, this is defined by ``GEARSink:minimal_discrete_mass_Msun`` for population II stars and ``GEARSink:minimal_discrete_mass_first_stars_Msun`` for population III stars. For these particles, we need to sample the IMF. This is explained in the next section.

IMF sampling
^^^^^^^^^^^^

To properly determine the number of supernovae and the mass of the ejected yields, we need to sample the IMF for SSP and continuous star particles. GEAR implements two methods `(Revaz et al. 2016) <https://ui.adsabs.harvard.edu/abs/2016A&A...588A..21R>`_:

- Continuous IMF sampling (CIMFS) (``GEARFeedback:discrete_yields: 0``): We integrate the quantities over the IMF and then explode a floating-point number of stars, which can be below 1 in some cases. This method works well for large stellar particle mass and supernovae rates, but not for low stellar particle mass or low supernovae rates. Note that SNIa often occur in the second regime, hence the method. The overall effect is similar to diluting the SN explosions over time. 


- Random discrete IMF sampling (RIMFS) (``GEARFeedback:discrete_yields: 1``): We avoid the issue of non-integer event numbers by taking the floor of the calculated SN count and stochastically adding an additional supernova based on the fractional part. We then compute the properties for a single star at a time.

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
  
* ``yields_table_first_stars``: Similar to ``yields_table`` but for the first stars.
  
* ``imf_transition_metallicity`` specifies which table to use based on [Fe/H]. If the gas metallicity is below ``imf_transition_metallicity``, we use the population III table; if it is above, we use the population table II.
  
* Two different models exist for the IMF sampling of supernovae (``GEARFeedback:discrete_yields``): the continuous mode (``GEARFeedback:discrete_yields: 1``) and the discrete mode (``GEARFeedback:discrete_yields: 0``).
  
* ``elements`` is the list of yields to read from the tables. The number of elements is specified at compile time  ``(--with-chemistry=GEAR_N)``

* ``discrete_star_minimal_gravity_mass_Msun``: Minimal gravity mass in solar masses after a discrete star completely explodes. Default: 0.1

* ``GEARSupernovaeII:interpolation_size`` is the number of elements to keep in the interpolation of the data.

Here is the whole feedback section:

.. code:: YAML

	  GEARFeedback:
	    supernovae_Ia_energy_erg: 1e51                           # Energy released by a single supernova.
	    supernovae_efficiency: 0.1                               # Supernovae energy efficiency, used for both SNIa and SNII. The energy released effectively is E_sn = supernovae_efficiency*E_sn
	    yields_table: chemistry-AGB+OMgSFeZnSrYBaEu-16072013.h5  # Table containing the yields.
	    yields_table_first_stars: chemistry-PopIII.hdf5          # Table containing the yields of the first stars.
	    imf_transition_metallicity: -5                           # Maximal metallicity ([Fe/H]) for a first star (0 to deactivate).
	    discrete_yields: 0                                       # Should we use discrete yields or the IMF integrated one?
	    elements: [Fe, Mg, O, S, Zn, Sr, Y, Ba, Eu]              # Elements to read in the yields table. The number of elements should be one less than the number of elements (N) requested during the configuration (--with-chemistry=GEAR_N).
	    discrete_star_minimal_gravity_mass_Msun: 0.1             # Minimal gravity mass after a discrete star completely explodes. In M_sun. (Default: 0.1)

	  GEARSupernovaeII:
	    interpolation_size:  200                                 # Number of elements to keep in the interpolation of the data. (Default: 200)

References
----------

- `pychem <https://www.astro.unige.ch/~revazy/PyChem/>`_ : python module to generate GEAR tables

- `Tsujimoto et al. (1995) <https://ui.adsabs.harvard.edu/abs/1995MNRAS.277..945T/abstract>`_.

- `Kobayashi et al. (2000) <https://ui.adsabs.harvard.edu/abs/2000ApJ...539...26K/abstract>`_

- `Kroupa (2001) <https://ui.adsabs.harvard.edu/abs/2001MNRAS.322..231K/abstract>`_

- `Poirier (2004) <https://theses.fr/2004STR13003>`_

- `Revaz et al. (2016) <https://ui.adsabs.harvard.edu/abs/2016A&A...588A..21R>`_

- `Hausammann (2021) <https://infoscience.epfl.ch/entities/publication/3e6d2e54-a782-440a-86c3-05482e83794d>`_
