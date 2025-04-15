.. Supernova feedback in GEAR model
   Darwin Roduit, 30 March 2025

.. gear_sn_feedback_models:

.. _gear_sn_feedback_models:

GEAR supernova feedback
=======================

When a star goes into a supernova, we compute the amount of internal energy, mass and metals the explosion transfers to the star's neighbouring gas particles. We will group all these in the “fluxes” term.  
We have two models for the distribution of these fluxes and the subgrid modelling of the supernovae: GEAR model and GEAR mechanical model.

.. note::
   We may sometimes refer to GEAR feedback as GEAR thermal feedback to clearly distinguish it from GEAR mechanical feedback.


.. _gear_sn_feedback_gear_thermal:

GEAR model
----------

In the GEAR (thermal) model, the fluxes are distributed by weighing with the SPH kernel:

.. math::

   w_{{sj}} = W_i(\| \vec{{x}}_{{sj}} \|, \, h_s) \frac{{m_j}}{{\rho_s}}

for :math:`s` the star and :math:`j` the gas (`Revaz and Jablonka 2012 <https://ui.adsabs.harvard.edu/abs/2012A%26A...538A..82R/abstract>`_).

In the GEAR model, we do not inject momentum, only *internal energy*. Then, internal energy conversion to kinetic energy is left to the hydrodynamic solver, which will compute appropriately the gas density, temperature and velocity.  
However, if the cooling radius :math:`R_{\text{cool}}` of the explosion is unresolved, i.e. the cooling radius is smaller than our simulation resolution, the cooling radiates away the internal energy.

To understand why this happens, let us remind the main phases of an SN explosion in a homogeneous medium. We provide a simple picture that is more complicated than the one explained here. See `Haid et al. 2016 <https://ui.adsabs.harvard.edu/abs/2016MNRAS.460.2962H/abstract>`_ or `Thornton et al. 1998 <https://iopscience.iop.org/article/10.1086/305704>`_ for further details.

* The first stage of the SN explosion is the **free expansion**. In this momentum-conserving regime, the ejecta of the stars sweeps the ISM. At the end of this phase, 72% of the initial SN energy has been converted to thermal energy.
* Once the SN ejecta has swept an ISM mass of comparable mass, the blast wave enters the **energy-conserving Sedov-Taylor phase**. It continues with an adiabatic expansion, performing some :math:`P \, \mathrm{d}V` work on the gas. In this phase, the internal energy is converted into kinetic energy as the ejecta continues sweeping the ISM gas. This phase continues until radiative losses become significant after some radius :math:`R_{\text{cool}}`.
* At this point, the blast wave enters the **momentum-conserving snowplough phase** and forms a thin shell. In this regime, efficient cooling radiates away the internal energy, and thus, the blast wave slows down.

Now, we better understand why the internal energy is radiated away. It is a consequence of efficient cooling in the snowplough phase. When this happens, the feedback is unresolved and its energy does not affect the ISM, apart from the mass and metal injection. To circumvent this problem, GEAR thermal feedback implements a **fixed delayed cooling mechanism**. The cooling of the particles affected by feedback is deactivated during some mega year, usually 5 Myr in our simulations. The time is controlled by the ``GrackleCooling:thermal_time_myr`` parameter. This mechanism allows the internal energy to transform into kinetic energy without immediately being radiated away. However, such an approach poses the question of the time required to prevent gas from cooling in the simulations.
