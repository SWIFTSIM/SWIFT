.. Basic sub-grid model
   Matthieu Schaller, 20th December 2018


Basic model (others)
====================

Sinks: Simple Bondi-Hoyle accretion
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

The ``Basic`` sink model provides a foundation on which new sink implementations could be built. It includes a prescription for Bondi-Hoyle gas accretion, and a method for sink-sink mergers that is a slightly simplified version of the implementation used in GEAR.

No other physics is implemented for this model. Sinks cannot form - to use this model, sink particles must already be present in the initial conditions. They also cannot spawn stars from the gas they accrete.

Bondi-Hoyle accretion can be done in one of two ways:

 * Gas particles within the sink's kernel are stochastically swallowed entirely, with a probability set by the Bondi-Hoyle rate. Specifically, the probability is set by the current difference between the sink's subgrid mass (determined by the accretion rate) and its dynamical mass (which tracks the number of particles/sinks actually swallowed). This mode is equivalent to the EAGLE black hole accretion model.
 * Gas particles within the sink's kernel are "nibbled" down to some minimal mass, which can be specified by the user. This method is equivalent to the black hole accretion model of Bahe et al. 2022.

This model has only two parameters that must be specified in your parameter ``yml`` file:

 * ``BasicSink:use_nibbling``: determines whether accretion is done by "nibbling" or by swallowing outright.
 * ``BasicSink:min_gas_mass_for_nibbling_Msun``: if using "nibbling", the minimum mass to which gas particles can be nibbled. A good default is half the original particle mass.

For an even more bare-bones starting point, the ``Default`` sink model contains no physics at all, and is a totally blank canvas on which to build your sink model.


Cooling: Analytic models
~~~~~~~~~~~~~~~~~~~~~~~~

Currently, we have 3 different simple cooling models (const-lambda, const-du
and Compton). These are all based on analytic formulas and can be used
to quickly understand how the cooling interacts with the rest of the
code before moving to more complex models.

Equations
---------

The first table compares the different analytical cooling while the next ones
are specific to a given cooling.  The quantities are the internal energy (\\( u
\\)), the density \\( rho \\), the element mass fraction (\\( X_i \\)), the
cooling function (\\(\\Lambda\\), the proton mass (\\( m_H \\)) and the time
step condition (\\( t\_\\text{step}\\)).  If not specified otherwise, all
cooling contains a temperature floor avoiding negative temperature.

.. csv-table:: Analytical Cooling
   :header: "Variable", "Const-Lambda", "Const-du"

   "\\( \\frac{ \\mathrm{d}u }{ \\mathrm{d}t } \\)", "\\( -\\Lambda \\frac{\\rho^2 X_H^2}{\\rho m_H^2} \\)", "const"
   "\\( \\Delta t\_\\text{max} \\)", "\\( t\_\\text{step} \\frac{u}{\\left|\\frac{ \\mathrm{d}u }{ \\mathrm{d}t }\\right|} \\)", "\\( t\_\\text{step} \\frac{u}{\\ \\left| \\frac{ \\mathrm{d}u }{ \\mathrm{d}t }\\right|} \\)"

TODO: Add description of the parameters and units.

TODO: Add Compton cooling model

Cooling: TREECOOL
~~~~~~~~~~~~~~~~~

The ``TREECOOL`` model implements the primordial cooling function of `Katz,
Weinberg & Hernquist (1996)
<https://ui.adsabs.harvard.edu/abs/1996ApJS..105...19K>`_ (hereafter KWH96).
The gas is assumed to be primordial (Hydrogen and Helium only, with the mass
fractions taken from the primordial Helium fraction in the physical constants)
and to be in ionization equilibrium with both the collisional processes and an
optically thin, spatially uniform UV background.

The abundances of HI, HII, HeI, HeII, HeIII and of the electrons are obtained
by solving eq. 33-38 of KWH96, and the following processes contribute to the
net cooling rate (KWH96, Table 1):

 * collisional excitation of HI and HeII,
 * collisional ionization of HI, HeI and HeII,
 * radiative recombination of HII, HeII and HeIII, plus the dielectronic
   recombination of HeII,
 * free-free emission (Bremsstrahlung),
 * inverse Compton cooling off the CMB,
 * photo-heating by the UV background.

The rate coefficients are the fits of KWH96, Table 2 (mostly from `Cen 1992
<https://ui.adsabs.harvard.edu/abs/1992ApJS...78..341C>`_). They are tabulated
on a regular grid in :math:`\log_{10}(T)` at start-up, so no external table of
rates is needed. Outside the tabulated range the gas is assumed to be entirely
neutral (below :math:`T_{\rm min}`) or entirely ionized (above
:math:`T_{\rm max}`, where only free-free and Compton cooling remain).

This is deliberately the basic model of KWH96 and of the ``COOLING`` option of
Gadget/Arepo. It does **not** include metal-line cooling, molecular cooling,
self-shielding of the UV background, the HeII reionization heating boost, or
any local radiation sources. The composition is fixed for the whole run: the
Hydrogen and Helium mass fractions are derived from the primordial Helium
fraction of the physical constants, and the element abundances tracked by the
chemistry module, if any, are ignored.

The UV background is read from a ``TREECOOL`` file, the plain-text format used
by most of the widely distributed UV background models. It has one row per
redshift and seven columns:

.. code-block:: none

   log10(1 + z)   Gamma_HI  Gamma_HeI  Gamma_HeII   eps_HI  eps_HeI  eps_HeII

where the :math:`\Gamma` are photo-ionization rates in
:math:`\rm s^{-1}` and the :math:`\epsilon` photo-heating rates in
:math:`\rm erg\,s^{-1}`, both per ion of the corresponding species. The
rates are interpolated logarithmically in :math:`\log_{10}(1+z)` once per
time-step. Above the highest redshift covered by the file the UV background is
switched off entirely, and the gas cools by collisional processes alone. An
example file, ``TREECOOL_UV_background.txt`` (Faucher-Giguère et al. 2009,
December 2011 update), is provided with the
``examples/Cooling/ConstantCosmoTempEvolution`` example, which runs a uniform
cosmological box with this model and compares the resulting thermal history of
the IGM to observations. Files in this format for the commonly used UV
background models (Faucher-Giguère et al., Haardt & Madau, Puchwein et al.)
are distributed with Gadget-4 and Arepo and by the authors of the models. The
rows must be in increasing order of redshift; many files are padded with rows
of zeros beyond the end of the model, and reading stops at the first such row.
At least two usable rows are required.

Note that the UV background is switched on instantaneously when the redshift
first drops below the end of the table (or below
``UV_background_start_redshift``), which produces a sudden jump in the
photo-heating rate at that time.

The energy is integrated implicitly: if the change in energy over the
time-step is small the explicit solution is used, and otherwise the equation
:math:`u_{\rm new} = u_{\rm old} + \Lambda(u_{\rm new})\,\mathrm{d}t` is
solved by bisection. The model therefore imposes no cooling time-step
criterion of its own.

Two limitations of the model are worth keeping in mind. First, the gas is
assumed to be entirely neutral below :math:`T_{\rm min}` while, in the
presence of the UV background, low-density gas just above :math:`T_{\rm min}`
is highly ionized. The mean molecular weight, and hence the relation between
internal energy and temperature, is therefore discontinuous at
:math:`T_{\rm min}`, and gas whose internal energy falls in the gap is
assigned :math:`T = T_{\rm min}`. To stay clear of this regime,
``log10_T_min`` should be set well below the lowest temperature the gas can
reach in the run (the reference implementation in Arepo uses one tenth of the
minimal gas temperature). Second, in non-cosmological runs the redshift is
zero, so the :math:`z = 0` entry of the ``TREECOOL`` table is used and the
Compton term cools the gas towards the present-day CMB temperature. The Arepo
implementation switches Compton cooling off in that case; set
``with_Compton_cooling: 0`` to reproduce this behaviour.

To use this model, configure with ``--with-cooling=TREECOOL`` and give the
following parameters:

.. code-block:: yaml

   TREECOOLCooling:
     TREECOOL_file:                ./TREECOOL_UV_background.txt
     rapid_cooling_threshold:      0.333333   # (Optional)
     UV_background_start_redshift: 1e30       # (Optional)
     with_Compton_cooling:         1          # (Optional)
     log10_T_min:                  1.0        # (Optional)
     log10_T_max:                  9.0        # (Optional)

Only ``TREECOOL_file`` is compulsory. ``rapid_cooling_threshold`` sets the
value of :math:`\mathrm{d}t / t_{\rm cool}` above which the new energy is
applied directly rather than through the time derivative of the internal
energy; a negative value always uses the latter.
``UV_background_start_redshift`` can be used to delay the onset of the UV
background beyond the range of the table itself, and ``with_Compton_cooling``
switches the inverse Compton cooling off the CMB on (default) or off. Finally,
``log10_T_min`` and ``log10_T_max`` set the range over which the rate
coefficients are tabulated.

In addition to the radiated energies, this model writes two fields per gas
particle to the snapshots. ``Temperatures`` are computed from the internal
energies using the mean molecular weight of the equilibrium ionization state,
so they differ from a conversion at a fixed mean molecular weight, in
particular for partially ionized gas. ``ElectronFractions`` are the electron
number densities in units of the Hydrogen number densities, again assuming
ionization equilibrium; they are set for all particles at the start of the run
and updated each time a particle is cooled. The parameters of the model are
recorded in the ``SubgridScheme`` group of the snapshots.

How to Implement a New Cooling
------------------------------

The developer should provide at least one function for:
 * writing the cooling name in HDF5
 * cooling a particle
 * the maximal time step possible
 * initializing a particle
 * computing the total energy radiated by a particle
 * initializing the cooling parameters
 * printing the cooling type

For implementation details, see ``src/cooling/none/cooling.h``

See :ref:`new_option` for the full list of changes required.
