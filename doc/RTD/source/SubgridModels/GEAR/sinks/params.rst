.. Sink particles in GEAR model
   Darwin Roduit, 15 March 2024

.. sink_GEAR_model:

.. _sink_GEAR_parameters:

Model parameters
----------------

The parameters of the GEAR sink model are grouped into the ``GEARSink`` section of the parameter file. 

The first two parameters are:

* The sink cut-off radius for gas and sink accretion: ``cut_off_radius``,
* The sink inner accretion radius fraction in terms of the cut-off radius: ``f_acc``.

The ``f_acc`` parameter is optional. Its default value is :math:`0.8`. Its value must respect :math:`0 \leq f_\text{acc} \leq 1` . It describes the inner radius :math:`f_{\text{acc}} \cdot r_{\text{cut-off}}` in which gas particles are swallowed without any further checks, as explained below. 

The next two mandatory parameters are:

* the gas maximal temperature to form a sink:``maximal_temperature``,
* the gas threshold density to form a sink:``density_threshold_g_per_cm3``.

These two parameters govern the first two criteria of the sink formation scheme. If these criteria are not passed, sink particles are not created. If they are passed, the code performs further checks to form sink particles. Some of those criteria checks can be disabled, as explained below.

The next set of parameters deals with the sampling of the IMF and the representation of star particles:

* minimal mass of stars represented by discrete particles: ``minimal_discrete_mass_Msun``,
* mass of the stellar particle representing the continuous part of the IMF: ``stellar_particle_mass_Msun``,
* minimal mass of the first stars represented by discrete particles: ``minimal_discrete_mass_first_stars_Msun``,
* mass of the stellar particle (first stars) representing the continuous part of the IMF:: ``stellar_particle_mass_first_stars_Msun``.

With sink particles, star particles can represent either a single star or a population of stars in the low mass part of the IMF (continuous IMF sampling). The stars in the continuous part of the IMF are put together in a particle of mass ``stellar_particle_mass_Msun`` or ``stellar_particle_mass_first_stars_Msun``, while individual stars in the discrete part have their mass sampled from the IMF. The limit between the continuous and discrete sampling of the IMF is controlled by  ``minimal_discrete_mass_Msun`` and ``minimal_discrete_mass_first_stars_Msun``.

The next set of parameters controls the sink formation scheme. More details are provided in the GEAR documentation. Here is a brief overview:

* whether or not the gas must be contracting: ``sink_formation_contracting_gas_criterion`` (default: 1), 
* whether or not the gas smoothing length must be small enough: ``sink_formation_smoothing_length_criterion`` (default: 1),
* whether or not the gas must be in a Jeans unstable state: ``sink_formation_jeans_instability_criterion`` (default: 1),
* whether or not the gas must be in a bound state: ``sink_formation_bound_state_criterion`` (default: 1),
* whether or not a new sink can be formed in a region where its ``cut_off_radius`` and the one of an existing sink overlap: ``sink_formation_overlapping_sink_criterion`` (default: 1).

Those criteria are checked if the density and temperature criteria are successfully passed. They control the behaviour of the sink formation scheme. By default, they are all activated and set to ``1``.

The last parameter is ``disable_sink_formation`` (default: 0). It controls whether sinks are formed or not in the simulation. The main purpose is when sinks are put in initial conditions and sinks are not wanted to be added during the run. This parameter is set to ``0`` by default, i.e. sink formation is *enabled*. 

The full section is:

.. code:: YAML
	  
   GEARSink:
     cut_off_radius:        1e-3                 # Cut off radius of the sink particles (in internal units).
     f_acc: 0.8                                  # (Optional) Fraction of the cut_off_radius that determines if a gas particle should be swallowed wihtout additional check. (Default: 0.8)
     maximal_temperature:        3e3             # Maximal gas temperature for forming a star (in K)
     density_threshold_g_per_cm3: 1.67e-21       # Minimal gas density for forming a star (in g/cm3 (1.67e-24 =1acc))
     stellar_particle_mass_Msun:      20              # Mass of the stellar particle representing the low mass stars (continuous IMF sampling) (in solar mass)
     minimal_discrete_mass_Msun:      8               # Minimal mass of stars represented by discrete particles (in solar mass)
     stellar_particle_mass_first_stars_Msun: 20       # Mass of the stellar particle representing the low mass stars (continuous IMF sampling) (in solar mass). First stars
     minimal_discrete_mass_first_stars_Msun: 8        # Minimal mass of stars represented by discrete particles (in solar mass). First stars
     star_spawning_sigma_factor: 0.2             # Factor to rescale the velocity dispersion of the stars when they are spawned. (Default: 0.2)
     sink_formation_contracting_gas_criterion: 1     # (Optional) Activate the contracting gas criterion for sink formation. (Default: 1)
     sink_formation_smoothing_length_criterion: 1    # (Optional) Activate the smoothing length criterion for sink formation. (Default: 1)
     sink_formation_jeans_instability_criterion: 1   # (Optional) Activate the two Jeans instability criteria for sink formation. (Default: 1)
     sink_formation_bound_state_criterion: 1         # (Optional) Activate the bound state criterion for sink formation. (Default: 1)
     sink_formation_overlapping_sink_criterion: 1    # (Optional) Activate the overlapping sink criterion for sink formation. (Default: 1)
     disable_sink_formation: 0                   # (Optional) Disable sink formation. (Default: 0)

.. warning::
   Some parameter choices can greatly impact the outcome of your simulations. Think twice when choosing them.

Sink accretion radius
~~~~~~~~~~~~~~~~~~~~~

The most critical parameter is ``cut_off_radius``. As explained in the theory, to form a sink, the gas smoothing kernel edge :math:`\gamma_k h` (:math:`\gamma_k` is a kernel dependent constant) must be smaller than ``cut_off_radius`` (if this criterion is enabled). Therefore, the cut-off radius strongly depends on the resolution of your simulations. Moreover, if you use a minimal gas smoothing length `h`, and plan to use sink particles, consider whether the cut-off radius will meet the smoothing length criterion. If `h` never meets the aforementioned criterion, you will never form sinks and thus never have stars.

On the contrary, if you set a too high cut-off radius, then sinks will accrete a lot of gas particles and spawn a lot of stars in the same cell, which the code might not like and crash with the error:

``runner_others.c:runner_do_star_formation_sink():274: Too many stars in the cell tree leaf! The sorting task will not be able to perform its duties. Possible solutions: (1) The code need to be run with different star formation parameters to reduce the number of star particles created. OR (2) The size of the sorting stack must be increased in runner_sort.c.``

This problem can be mitigated by choosing a higher value of ``stellar_particle_mass_Msun`` and ``stellar_particle_mass_first_stars_Msun``, or higher values of ``minimal_discrete_mass_Msun`` and ``minimal_discrete_mass_first_stars_Msun``. Of course, this comes at the price of having fewer individual stars. Finally, all parameters will depend on your needs.

*If you do not want to change your parameters*, you can increase the ``sort_stack_size`` variable at the beginning ``runner_sort.c``. The default value is 10 in powers of 2 (so the stack size is 1024 particles). Increase it to the desired value. Be careful to not overestimate this.

Guide to choose the the accretion radius or the density threshold
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We provide some advice to help you set up the sink accretion radius or the threshold density appropriately.

First, you must choose either the sink accretion radius or the threshold density. Choosing the density might be easier based on your previous work or if you have an expected star formation density. Once you fix the density or the accretion radius, you can use the following formula to estimate the remaining parameter. In the code, the gas smoothing length is determined with:

.. math::
   h = \eta \left( \frac{X_{\text{H}} m_B}{m_{\text{H}} n_{\text{H}}} \right)^{1/3} \, ,

where :math:`\eta` is a constant related to the number of neighbours in the kernel, :math:`X_{\text{H}}` is the hydrogen mass fraction, :math:`m_B` the gas particle's mass, :math:`m_{\text{H}}` the hydrogen particle mass and :math:`n_{\text{H}}` the hydrogen number density.

Let us provide an example. In GEAR, we do not model physical processes below the parsec scale. Hence, let us take :math:`h \sim 1` pc. In zoom-in simulations we have :math:`m_B \simeq 95 \; M_{\odot}`. The remaining parameters are :math:`\eta = 1.2348` and :math:`X_{\text{H}} = 0.76`. So, after inverting the formula, we find :math:`n_H \simeq 5500 \text{ hydrogen atoms/cm}^3`. In practice, we use :math:`n_H = 1000 \text{ hydrogen atoms/cm}^3`, close to the estimation, and an accretion radius :math:`r_{\text{acc}} = 10` pc. These values are slightly different for safety reasons, but they are consistent.

Remember that this was a way, among others, to determine good accretion radius and threshold density. It can help you with your first runs with sink particles.

Comment on star formation efficiency
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Notice that this model does not have parameters to control the star formation rate of the sink. The SFR is self-regulated by the gas/sink accretion and other feedback mechanisms. Supernovae tend to create bubbles of lower density at the site of star formation, removing the gas and preventing further gas accretion. However, the sink might run into this stack size problem by the time the first supernovae explode. Other pre-stellar feedback mechanisms could do the job earlier, though they are not implemented in GEAR.

.. note:: 
   We provide a piece of general advice: do some calibration on low-resolution simulations. This will help to see what works and what does not work. Keep in mind that you might want to put a higher ``stellar_particle_mass_X_Msun`` at the beginning to avoid spawning too many stars. For the high-resolution simulations, you then can lower the particle's mass.
