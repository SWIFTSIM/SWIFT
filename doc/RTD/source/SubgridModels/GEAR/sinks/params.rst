.. Sink particles in GEAR model
   Darwin Roduit, 15 March 2024

.. sink_GEAR_model:

.. warning::
  This page is under construction. It may lack information. 

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
* the gas threshold density to form a sink:``density_threshold``.

These two parameters govern the first two criteria of the sink formation scheme. If these criteria are not passed, sink particles are not created. If they are passed, the code performs further checks to form sink particles. Some of those criteria checks can be disabled, as explained below.

The next set of parameters deals with the sampling of the IMF and the representation of star particles:

* minimal mass of stars represented by discrete particles: ``minimal_discrete_mass``,
* mass of the stellar particle representing the continuous part of the IMF: ``stellar_particle_mass``,
* minimal mass of the first stars represented by discrete particles: ``minimal_discrete_mass_first_stars``,
* mass of the stellar particle (first stars) representing the continuous part of the IMF:: ``stellar_particle_mass_first_stars``.

With sink particles, star particles can represent either a single star or a population of stars in the low mass part of the IMF (continuous IMF sampling). The stars in the continuous part of the IMF are put together in a particle of mass ``stellar_particle_mass`` or ``stellar_particle_mass_first_stars``, while individual stars in the discrete part have their mass sampled from the IMF. The limit between the continuous and discrete sampling of the IMF is controlled by  ``minimal_discrete_mass`` and ``minimal_discrete_mass_first_stars``.

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
     density_threshold:          1.67e-21        # Minimal gas density for forming a star (in g/cm3 (1.67e-24 =1acc))
     stellar_particle_mass:      20              # Mass of the stellar particle representing the low mass stars (continuous IMF sampling) (in solar mass)
     minimal_discrete_mass:      8               # Minimal mass of stars represented by discrete particles (in solar mass)
     stellar_particle_mass_first_stars: 20       # Mass of the stellar particle representing the low mass stars (continuous IMF sampling) (in solar mass). First stars
     minimal_discrete_mass_first_stars: 8        # Minimal mass of stars represented by discrete particles (in solar mass). First stars
     sink_formation_contracting_gas_criterion: 1     # (Optional) Activate the contracting gas criterion for sink formation. (Default: 1)
     sink_formation_smoothing_length_criterion: 1    # (Optional) Activate the smoothing length criterion for sink formation. (Default: 1)
     sink_formation_jeans_instability_criterion: 1   # (Optional) Activate the two Jeans instability criteria for sink formation. (Default: 1)
     sink_formation_bound_state_criterion: 1         # (Optional) Activate the bound state criterion for sink formation. (Default: 1)
     sink_formation_overlapping_sink_criterion: 1    # (Optional) Activate the overlapping sink criterion for sink formation. (Default: 1)
     disable_sink_formation: 0                   # (Optional) Disable sink formation. (Default: 0)

.. warning::
   Some parameter choices can greatly impact the outcome of your simulations. Think twice when choosing them.

The most critical parameter is ``cut_off_radius``. As explained in the theory, to form a sink, the gas smoothing length `h` must be smaller than ``cut_off_radius / 2`` (if this criterion is enabled). Therefore, the cut-off radius strongly depends on the resolution of your simulations. Moreover, if you use a minimal gas smoothing length `h`, and plan to use sink particles, consider whether the cut-off radius will meet the smoothing length criterion. If `h` never meets the aforementioned criterion, you will never form sinks and thus never have stars.

On the contrary, if you set a too high cut-off radius, then sinks will accrete a lot of gas particles and spawn a lot of stars in the same cell, which the code might not like and crash with the error:

``runner_others.c:runner_do_star_formation_sink():274: Too many stars in the cell tree leaf! The sorting task will not be able to perform its duties. Possible solutions: (1) The code need to be run with different star formation parameters to reduce the number of star particles created. OR (2) The size of the sorting stack must be increased in runner_sort.c.``

This problem can be mitigated by choosing a higher value of ``stellar_particle_mass`` and ``stellar_particle_mass_first_stars``, or higher values of ``minimal_discrete_mass`` and ``minimal_discrete_mass_first_stars``. Of course, this comes at the price of having fewer individual stars. Finally, all parameters will depend on your needs.

*If you do not want to change your parameters*, you can increase the ``sort_stack_size`` variable at the beginning ``runner_sort.c``. The default value is 10 in powers of 2 (so the stack size is 1024 particles). Increase it to the desired value. Be careful to not overestimate this.

Notice that this model does not have parameters to control the star formation rate of the sink. The SFR is self-regulated by the gas/sink accretion and other feedback mechanisms. Supernovae tend to create bubbles of lower density at the site of star formation, removing the gas and preventing further gas accretion. However, the sink might run into this stack size problem by the time the first supernovae explode. Other pre-stellar feedback mechanisms could do the job earlier, though they are not implemented in GEAR.

.. note:: 
   We provide a piece of general advice: do some calibration on low-resolution simulations. This will help to see what works and what does not work. Keep in mind that you might want to put a higher ``stellar_particle_mass_X`` at the beginning to avoid spawning too many stars. For the high-resolution simulations, you then can lower the particle's mass.
