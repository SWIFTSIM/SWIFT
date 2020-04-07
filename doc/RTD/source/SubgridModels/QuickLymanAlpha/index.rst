.. Quick Lyman-alpha sub-grid model
   Matthieu Schaller, 7th March 2020


Quick Lyman-alpha (QLA) model
=============================

This section of the documentation gives a brief description of the different
components of the quick Lyman-alpha sub-grid model. We mostly focus on the
parameters and values output in the snapshots.

Given the nature of the model, no feedback or black holes are used. The star
formation model is minimalist and the chemistry/cooling models are limited to
primordial abundances.

.. _QLA_entropy_floors:

Gas entropy floor
~~~~~~~~~~~~~~~~~

The gas particles in the QLA model are prevented from cooling below a certain
temperature. The temperature limit depends on the density of the particles. The
floor is implemented as a polytropic "equation of state":math:`P = P_c
\left(\rho/\rho_c\right)^\gamma` (all done in physical coordinates), with the
constants derived from the user input given in terms of temperature and Hydrogen
number density. We use :math:`gamma=1` in this model. The code computing the
entropy floor is located in the directory ``src/entropy_floor/QLA/`` and the
floor is applied in the drift and kick operations of the hydro scheme.

An additional over-density criterion above the mean baryonic density is applied
to prevent gas not collapsed into structures from being affected. To be precise,
this criterion demands that the floor is applied only if :math:`\rho_{\rm com} >
\Delta_{\rm floor}\bar{\rho_b} = \Delta_{\rm floor} \Omega_b \rho_{\rm crit,0}`,
with :math:`\Delta_{\rm floor}` specified by the user, :math:`\rho_{\rm crit,0}
= 3H_0/8\pi G` the critical density at redshift zero [#f1]_, and
:math:`\rho_{\rm com}` the gas co-moving density. Typical values for
:math:`\Delta_{\rm floor}` are of order 10.

The model is governed by 3 parameters for each of the two limits. These are
given in the ``QLAEntropyFloor`` section of the YAML file. The parameters are
the Hydrogen number density (in :math:`cm^{-3}`) and temperature (in :math:`K`)
of the anchor point of the floor as well as the minimal over-density required to
apply the limit. To simplify things, all constants are converted to the internal
system of units upon reading the parameter file.

For a normal quick Lyman-alpha run, that section of the parameter file reads:

.. code:: YAML

  QLAEntropyFloor:
    density_threshold_H_p_cm3: 0.1       # Physical density above which the entropy floor kicks in expressed in Hydrogen atoms per cm^3.
    over_density_threshold:    10.       # Over-density above which the entropy floor can kick in.
    temperature_norm_K:        8000      # Temperature of the entropy floor at the density threshold expressed in Kelvin.


SWIFT will convert the temperature normalisations and Hydrogen number density
thresholds into internal energies and densities respectively assuming a neutral
gas with primordial abundance pattern. This implies that the floor may not be
exactly at the position given in the YAML file if the gas has different
properties. This is especially the case for the temperature limit which will
often be lower than the imposed floor by a factor :math:`\frac{\mu_{\rm
neutral}}{\mu_{ionised}} \approx \frac{1.22}{0.59} \approx 2` due to the
different ionisation states of the gas.

Recall that we additionally impose an absolute minimum temperature at all
densities with a value provided in the :ref:`Parameters_SPH` section of the
parameter file. This minimal temperature is typically set to 100 Kelvin.


.. _QLA_cooling:
     
Gas cooling: Wiersma+2009a with fixed primordial metallicity
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The gas cooling is based on the redshift-dependent tables of `Wiersma et
al. (2009)a <http://adsabs.harvard.edu/abs/2009MNRAS.393...99W>`_ that include
element-by-element cooling rates for the 11 elements (`H`, `He`, `C`, `N`, `O`,
`Ne`, `Mg`, `Si`, `S`, `Ca` and `Fe`) that dominate the total rates. The tables
assume that the gas is in ionization equilibrium with the cosmic microwave
background (CMB) as well as with the evolving X-ray and UV background from
galaxies and quasars described by the model of `Haardt & Madau (2001)
<http://adsabs.harvard.edu/abs/2001cghr.confE..64H>`_. Note that this model
ignores *local* sources of ionization, self-shielding and non-equilibrium
cooling/heating. The tables can be obtained from this `link
<http://virgodb.cosma.dur.ac.uk/swift-webstorage/CoolingTables/EAGLE/coolingtables.tar.gz>`_
which is a re-packaged version of the `original tables
<http://www.strw.leidenuniv.nl/WSS08/>`_. The code reading and interpolating the
table is located in the directory ``src/cooling/QLA/``.

The Wiersma tables containing the cooling rates as a function of redshift,
Hydrogen number density, Helium fraction (:math:`X_{He} / (X_{He} + X_{H})`) and
element abundance relative to the solar abundance pattern assumed by the tables
(see equation 4 in the original paper). Since the quick Lyman-alpha model is
only of interest for gas outside of haloes, we can make use of primordial gas
only. This means that the particles do not need to carry a metallicity array or
any individual element arrays. Another optimization is to ignore the cooling
rates of the metals in the tables.

Above the redshift of Hydrogen re-ionization we use the extra table containing
net cooling rates for gas exposed to the CMB and a UV + X-ray background at
redshift nine truncated above 1 Rydberg. At the redshift or re-ionization, we
additionally inject a fixed user-defined amount of energy per unit mass to all
the gas particles.

In addition to the tables we inject extra energy from Helium II re-ionization
using a Gaussian model with a user-defined redshift for the centre, width and
total amount of energy injected per unit mass. Additional energy is also
injected instantaneously for Hydrogen re-ionisation to all particles (active and
inactive) to make sure the whole Universe reaches the expected temperature
quickly (i.e not just via the interaction with the now much stronger UV
background).

The cooling itself is performed using an implicit scheme (see the theory
documents) which for small values of the cooling rates is solved explicitly. For
larger values we use a bisection scheme.  The cooling rate is added to the
calculated change in energy over time from the other dynamical equations. This
is different from other commonly used codes in the literature where the cooling
is done instantaneously.

We note that the QLA cooling model does not impose any restriction on the
particles' individual time-steps. The cooling takes place over the time span
given by the other conditions (e.g the Courant condition).

Finally, the cooling module also provides a function to compute the temperature
of a given gas particle based on its density, internal energy, abundances and
the current redshift. This temperature is the one used to compute the cooling
rate from the tables and similarly to the cooling rates, they assume that the
gas is in collisional equilibrium with the background radiation. The
temperatures are, in particular, computed every time a snapshot is written and
they are listed for every gas particle:

+---------------------+-------------------------------------+-----------+-------------------------------------+
| Name                | Description                         | Units     | Comments                            |
+=====================+=====================================+===========+=====================================+
| ``Temperatures``    | | Temperature of the gas as         | [U_T]     | | The calculation is performed      |
|                     | | computed from the tables.         |           | | using quantities at the last      |
|                     |                                     |           | | time-step the particle was active |
+---------------------+-------------------------------------+-----------+-------------------------------------+

Note that if one is running without cooling switched on at runtime, the
temperatures can be computed by passing the ``--temperature`` runtime flag (see
:ref:`cmdline-options`). Note that the tables then have to be available as in
the case with cooling switched on. 

The cooling model is driven by a small number of parameter files in the
`QLACooling` section of the YAML file. These are the re-ionization parameters
and the path to the tables. A valid section of the YAML file looks like:

.. code:: YAML

   QLACooling:
     dir_name:     /path/to/the/Wiersma/tables/directory # Absolute or relative path
     H_reion_z:            11.5      # Redshift of Hydrogen re-ionization
     H_reion_ev_p_H:        2.0      # Energy injected in eV per Hydrogen atom for Hydrogen re-ionization.
     He_reion_z_centre:     3.5      # Centre of the Gaussian used for Helium re-ionization
     He_reion_z_sigma:      0.5      # Width of the Gaussian used for Helium re-ionization
     He_reion_ev_p_H:       2.0      # Energy injected in eV per Hydrogen atom for Helium II re-ionization.

.. _QLA_star_formation:

Star formation
~~~~~~~~~~~~~~

The star formation in the Quick Lyman-alpha model is very simple. Any gas
particle with a density larger than a multiple of the critical density for
closure is directly turned into a star. The idea is to rapidly eliminate any gas
that is found within bound structures since we are only interested in what
happens in the inter-galactic medium. The over-density multiple is the only
parameter of this model.

The code applying this star formation law is located in the directory
``src/star_formation/QLA/``. 

For a normal Quick Lyman-alpha run, that section of the parameter file reads:

.. code:: YAML

   # Quick Lyman-alpha star formation parameters
   QLAStarFormation:
     over_density:              1000      # The over-density above which gas particles turn into stars.                                                                                          


.. [#f1] Recall that in a non-cosmological run the critical density is
	 set to 0, effectively removing the over-density
	 constraint of the floors.
