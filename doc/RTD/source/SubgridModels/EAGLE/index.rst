.. EAGLE sub-grid model
   Matthieu Schaller, 20th December 2018


EAGLE model
===========

This section of the documentation gives a brief description of the
different components of the EAGLE sub-grid model. We mostly focus on
the parameters and values output in the snapshots.

Chemical tracers
~~~~~~~~~~~~~~~~

The gas particles in the EAGLE model carry metal abundance information
in the form of metal mass fractions. We follow the following 9
elements: `H`, `He`, `C`, `N`, `O`, `Ne`, `Mg`, `Si` and `Fe`. We
additionally follow the total metal mass fraction (i.e. absolute
metallicity) `Z`. This is typically larger than the sum of the 7
metals that are individually traced since this will also contain the
contribution of all the elements that are not individually followed.
We note that all of definitions are independent of any definition of
solar the solar metallicity :math:`Z_\odot` or of any solar abundance
pattern.

As part of the diagnostics, we additionally trace the elements coming
from the different stellar evolution channels. We store for each
particle the total mass coming from all the SNIa that enriched that
particle and the metal mass fraction from SNIa. This is the fraction
of the *total* gas mass that is in the form of metals originating from
SNIa stars. By construction this fraction will be smaller than the
total metal mass fraction. The same tracers exist for the SNII and AGB
channels. Finally, we also compute the iron gas fraction from
SNIa. This it the fraction of the *total* gas mass that is made of
iron originating from SNIa explosions. 

We finally also compute the smoothed versions of the individual
element mass fractions, of the total metal mass fractions, and of the
iron gas fraction from SNIa.

The chemistry module in ``src/chemistry/EAGLE`` includes all the arrays
that are added to the particles and the functions used to compute the
smoothed elements.

When a star is formed (see below), it inherits all the chemical
tracers of its parent gas particle.

In the snapshots, we output for each gas and star particle:

+----------------------------------+-------------------------------------+-----------+-----------------------------+
| Name                             | Description                         | Units     | Comments                    |
+==================================+=====================================+===========+=============================+
| ``ElementAbundance``             | | Fraction of the gas/star mass     | [-]       | | Array of length           |
|                                  | | in the different elements         |           | | 9 for each particle       |
+----------------------------------+-------------------------------------+-----------+-----------------------------+
| ``SmoothedElementAbundance``     | | Fraction of the gas/star mass     | [-]       | | Array of length           |
|                                  | | in the different elements         |           | | 9 for each particle       |
|                                  | | smoothed over SPH neighbours      |           |                             |
+----------------------------------+-------------------------------------+-----------+-----------------------------+
| ``Metallicity``                  | | Fraction of the gas/star mass     | [-]       |                             |
|                                  | | in *all* metals                   |           |                             |
+----------------------------------+-------------------------------------+-----------+-----------------------------+
| ``SmoothedMetallicity``          | | Fraction of the gas/star mass     | [-]       |                             |
|                                  | | in *all* metals                   |           |                             |
|                                  | | smoothed over SPH neighbours      |           |                             |
+----------------------------------+-------------------------------------+-----------+-----------------------------+
| ``TotalMassFromSNIa``            | | Total mass of the gas/star        | [U_M]     |                             |
|                                  | | that was produced by enrichment   |           |                             |
|                                  | | from SNIa stars                   |           |                             |
+----------------------------------+-------------------------------------+-----------+-----------------------------+
| ``MetalMassFracFromSNIa``        | | Fraction of the *total* gas/star  | [-]       |                             |
|                                  | | mass that is in metals produced   |           |                             |
|                                  | | by enrichment from SNIa stars     |           |                             |
+----------------------------------+-------------------------------------+-----------+-----------------------------+
| ``TotalMassFromAGB``             | | Total mass of the gas/star        | [U_M]     |                             |
|                                  | | that was produced by enrichment   |           |                             |
|                                  | | from AGB stars                    |           |                             |
+----------------------------------+-------------------------------------+-----------+-----------------------------+
| ``MetalMassFracFromAGB``         | | Fraction of the *total* gas/star  | [-]       |                             |
|                                  | | mass that is in metals produced   |           |                             |
|                                  | | by enrichment from AGB star       |           |                             |
+----------------------------------+-------------------------------------+-----------+-----------------------------+
| ``TotalMassFromSNII``            | | Total mass of the gas/star        | [U_M]     |                             |
|                                  | | that was produced by enrichment   |           |                             |
|                                  | | from SNII stars                   |           |                             |
+----------------------------------+-------------------------------------+-----------+-----------------------------+
| ``MetalMassFracFromSNII``        | | Fraction of the gas/star mass     | [-]       |                             |
|                                  | | that is in metals produced by     |           |                             |
|                                  | | enrichment from SNII stars        |           |                             |
+----------------------------------+-------------------------------------+-----------+-----------------------------+
| ``IronMassFracFromSNIa``         | | Fraction of the *total* gas/star  | [-]       |                             |
|                                  | | mass in *iron* produced produced  |           |                             |
|                                  | | by enrichment from SNIa stars     |           |                             |
+----------------------------------+-------------------------------------+-----------+-----------------------------+
| ``SmoothedIronMassFracFromSNIa`` | | Fraction of the *total* gas/star  | [-]       |                             |
|                                  | | mass in *iron* produced produced  |           |                             |
|                                  | | by enrichment from SNIa stars     |           |                             |
|                                  | | smoothed over SPH neighbours      |           |                             |
+----------------------------------+-------------------------------------+-----------+-----------------------------+

The stars will lose mass over their lifetime (up to ~45%). The
fractions will remain unchanged but if one is interested in computing
an absolute metal mass (say) for a star, the ``InitialMass`` (see
below) of the star must be used.

The chemistry model only requires a small number of parameters to be
specified in the YAML file. These are the initial values of the
metallicity and element mass fractions. These are then applied at the
start of a simulation to *all* the *gas* particles. All 9 elements have
to be specified An example section, for primordial abundances (typical
for a cosmological run), is:

.. code:: YAML

   EAGLEChemistry:
     Metallicity:                0.    # Mass fraction in all metals
     InitAbundance_Hydrogen:     0.755 # Mass fraction in Hydrogen
     InitAbundance_Helium:       0.245 # Mass fraction in Helium
     InitAbundance_Carbon:       0.    # Mass fraction in Carbon
     InitAbundance_Nitrogen:     0.    # Mass fraction in Nitrogen
     InitAbundance_Oxygen:       0.    # Mass fraction in Oxygen
     InitAbundance_Neon:         0.    # Mass fraction in Neon
     InitAbundance_Magnesium:    0.    # Mass fraction in Magnesium
     InitAbundance_Silicon:      0.    # Mass fraction in Silicon
     InitAbundance_Iron:         0.    # Mass fraction in Iron

Whilst one would use the following values for solar abundances
(typical for an idealised low-redshift run):

.. code:: YAML

   EAGLEChemistry:
     Metallicity:                0.014        # Mass fraction in all metals
     InitAbundance_Hydrogen:     0.70649785   # Mass fraction in Hydrogen
     InitAbundance_Helium:       0.28055534   # Mass fraction in Helium
     InitAbundance_Carbon:       2.0665436e-3 # Mass fraction in Carbon
     InitAbundance_Nitrogen:     8.3562563e-4 # Mass fraction in Nitrogen
     InitAbundance_Oxygen:       5.4926244e-3 # Mass fraction in Oxygen
     InitAbundance_Neon:         1.4144605e-3 # Mass fraction in Neon
     InitAbundance_Magnesium:    5.907064e-4  # Mass fraction in Magnesium
     InitAbundance_Silicon:      6.825874e-4  # Mass fraction in Silicon
     InitAbundance_Iron:         1.1032152e-3 # Mass fraction in Iron


     
Cooling: Wiersma+2008a
~~~~~~~~~~~~~~~~~~~~~~

Particle tracers
~~~~~~~~~~~~~~~~

Star formation: Schaye+2008
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Stellar enrichment: Wiersma+2008b
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Supernova feedback: Schaye+2012
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Black-hole creation
~~~~~~~~~~~~~~~~~~~

Black-hole accretion
~~~~~~~~~~~~~~~~~~~~

AGN feedback
~~~~~~~~~~~~
