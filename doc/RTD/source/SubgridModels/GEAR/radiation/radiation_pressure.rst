.. GEAR sub-grid model radiation pressure
   Darwin Roduit, 24th September 2026

.. _gear_radiation_pressure:

Radiation pressure
=====================

A star's bolometric light pushes on the dusty gas around it. GEAR models this with LEBRON, a local momentum-coupling scheme (Hopkins, Quataert and Murray 2012, MNRAS 421, 3488, Sec 2.1; Hopkins et al. 2014, MNRAS 445, 581, Appendix A) that treats the star's light as locally absorbed and re-radiated by the surrounding dust rather than solving radiative transfer in full:

.. math::
   \dot{p} = (1 - e^{-\tau_\mathrm{NUV}}) (1 + \tau_\mathrm{IR}) \frac{L_\mathrm{bol}}{c}

:math:`L_\mathrm{bol}` is the star's bolometric luminosity, read from the radiation table and scaled by ``GEARFeedback:radiation_pressure_efficiency``. :math:`\tau_\mathrm{NUV}` and :math:`\tau_\mathrm{IR}` are the optical depths of the gas around the star to the non-ionizing continuum (912 Å to 3 :math:`\mu`\ m) and to the dust-reprocessed infrared, computed from a Sobolev approximation to the local gas column density (the star's own SPH density and density gradient, capped at the kernel support radius) and metallicity-scaled opacities (:math:`\kappa_\mathrm{NUV} = 1800\,\mathrm{cm}^2/\mathrm{g} \times Z/Z_\odot`, :math:`\kappa_\mathrm{IR} = 10\,\mathrm{cm}^2/\mathrm{g} \times Z/Z_\odot`). Both opacities are population-averaged (calibrated against STARBURST99 spectra), not single-star values.

The resulting momentum is deposited on the star's SPH gas neighbours, kernel-weighted and directed radially outward from the star.

Configuring and compiling
--------------------------

See :ref:`gear_radiation` for how to configure and build this channel, including the ``--with-tracers=GEAR`` requirement shared by every radiation channel.

Model parameters
------------------

* ``radiation_pressure_efficiency``: dimensionless factor the code multiplies the star's bolometric luminosity by before computing :math:`\dot{p}` above. 0 (default) switches the channel off entirely; 1 uses the star's bolometric luminosity as read from the table, unboosted. A value above 1 scales the effective luminosity used for the momentum injection beyond the star's tabulated output.

.. code:: YAML

   GEARFeedback:
     radiation_pressure_efficiency: 0   # Multiplies L_bol before the momentum injection. 0 = off, 1 = unboosted (Default: 0)

Snapshot outputs
------------------

The gas fields this channel writes are documented on the :ref:`gear_output_radiation_pressure` section of the :ref:`gear_output_fields` page.

References
-----------

- Hopkins, Quataert and Murray (2012), MNRAS 421, 3488, Section 2.1.
- Hopkins et al. (2014), MNRAS 445, 581, Appendix A.
