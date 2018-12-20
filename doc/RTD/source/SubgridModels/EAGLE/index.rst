.. EAGLE sub-grid model
   Matthieu Schaller, 20th December 2018


EAGLE model
===========

Chemical tracers
~~~~~~~~~~~~~~~~

The gas particles in the EAGLE model carry metal abundance information
in the form of metal mass fractions. We follow the following 9
elements: `H`, `He`, `C`, `N`, `O`, `Ne`, `Mg`, `Si` and `Fe`. We
additionally follow the total metal mass fraction `Z`. This is
typically larger than the sum of the 7 metals that are individually
traced since this will also contain the contribution of all the
elements that are not individually followed.

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

In the snapshots, we output:

+------------------------------+--------------------------------+-----------+-----------------------------+
| Name                         | Description                    | Units     | Comments                    | 
+==============================+================================+===========+=============================+
| ``ElementAbundance``         | | Fraction of the gas mass     | [-]       | | Array of length           |
|                              | | in the different elements    |           | | 9 for each particle       |
+------------------------------+--------------------------------+-----------+-----------------------------+
| ``SmoothedElementAbundance`` | | Fraction of the gas mass     | [-]       | | Array of length           |
|                              | | in the different elements    |           | | 9 for each particle       |
|                              | | smoothed over SPH neighbours |           |                             |
+------------------------------+--------------------------------+-----------+-----------------------------+
| ``Metallicity``              | | Fraction of the gas mass     | [-]       |                             |
|                              | | in *all* metals              |           |                             |
+------------------------------+--------------------------------+-----------+-----------------------------+
| ``SmoothedMetallicity``      | | Fraction of the gas mass     | [-]       |                             |
|                              | | in *all* metals              |           |                             |
|                              | | smoothed over SPH neighbours |           |                             |
+------------------------------+--------------------------------+-----------+-----------------------------+

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
