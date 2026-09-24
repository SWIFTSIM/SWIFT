.. GEAR sub-grid model radiation pressure
   Darwin Roduit, 24th September 2026

.. _gear_radiation_pressure:

Radiation pressure
=====================

A star's bolometric light pushes on the dusty gas around it. GEAR models this with the LEBRON momentum-coupling scheme of Hopkins, Quataert and Murray (2012, MNRAS 421, 3488, Sec 2.1) and Hopkins et al. (2014, MNRAS 445, 581, Appendix A):

.. math::
   \dot{p} = (1 - e^{-\tau_\mathrm{NUV}}) (1 + \tau_\mathrm{IR}) \frac{L_\mathrm{bol}}{c}

:math:`L_\mathrm{bol}` is the star's bolometric luminosity, read from the radiation table and scaled by ``GEARFeedback:radiation_pressure_efficiency``. :math:`\tau_\mathrm{NUV}` and :math:`\tau_\mathrm{IR}` are the optical depths of the gas around the star to the non-ionizing continuum (912 Å to 3 :math:`\mu`\ m) and to the dust-reprocessed infrared, computed from a Sobolev approximation to the local gas column density (the star's own SPH density and density gradient, capped at the kernel support radius) and metallicity-scaled opacities (:math:`\kappa_\mathrm{NUV} = 1800\,\mathrm{cm}^2/\mathrm{g} \times Z/Z_\odot`, :math:`\kappa_\mathrm{IR} = 10\,\mathrm{cm}^2/\mathrm{g} \times Z/Z_\odot`). Both opacities are population-averaged (calibrated against STARBURST99 spectra), not single-star values.

The resulting momentum is deposited on the star's SPH gas neighbours, kernel-weighted and directed radially outward from the star.

Model parameters
------------------

* ``radiation_pressure_efficiency``: dimensionless factor multiplying the star's bolometric luminosity before it is used above. 0 (default) switches the channel off entirely.

.. code:: YAML

   GEARFeedback:
     radiation_pressure_efficiency: 0   # Boost factor on L_bol for the radiation-pressure momentum injection. 0 = off (Default: 0)

Snapshot outputs
------------------

These gas fields need ``--with-tracers=GEAR``:

.. list-table::
   :header-rows: 1
   :widths: 35 45 20

   * - Name
     - Description
     - Units
   * - ``CumulativeMomentumFromRadiationPressure``
     - Norm of the momentum received from radiation pressure, summed over events (scalar sum of the norms, so isotropic kicks do not cancel)
     - [U_M U_L U_T\ :sup:`-1`\ ]
   * - ``MaxKickVelocityFromRadiationPressure``
     - Largest single-event velocity kick this particle received from radiation pressure
     - [U_L U_T\ :sup:`-1`\ ]

References
-----------

- Hopkins, Quataert and Murray (2012), MNRAS 421, 3488, Section 2.1.
- Hopkins et al. (2014), MNRAS 445, 581, Appendix A.
