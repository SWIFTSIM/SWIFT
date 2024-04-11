.. Sink particles in GEAR model
   Darwin Roduit, 15 March 2024

.. sink_GEAR_model:

.. warning::
  This page is under construction. It may lack information. 

Model summary
-------------

Here, we provide a comprehensive summary of the model. Sink particles are an alternative to the current model of star formation that transforms gas particles into sink particles under some criteria explained below. Then, the sink can accrete gas and spawn stars. 

To spawn stars, an IMF is sampled. Details explanation of the IMF sampling are explained below. In short, the IMF is split into two parts In the lower part, star particles represent continuous stellar population in a similar way than what is currently implemented in common models. In the second upper part, star particles represent individual stars. Then, the feedback is improved to take into account both type of stars. Currently, only supernovae feedback is implemented. The sink particle method allows thus to track the effect of the supernovae of individual stars in the simulation.

The current model includes sink formation, gas accretion, sink merging, IMF sampling and star spawning and finally supernove feedback (type Ia and II).

Our main references are the following papers `Bate et al. <https://ui.adsabs.harvard.edu/abs/1995MNRAS.277..362B/abstract>`_, `Price et al. <https://ui.adsabs.harvard.edu/abs/2018PASA...35...31P/abstract>`_ and `Federrath et al. <https://ui.adsabs.harvard.edu/abs/2010ApJ...713..269F/abstract>`_

.. figure:: ../../../Task/sink.png
    :width: 400px
    :align: center
    :figclass: align-center
    :alt: Task dependencies for the sink scheme.

    This figure shows the task dependencies for the sink scheme.
    The first rectangle groups the tasks that determine if sink particles will swallow other
    sink particles or gas particles.
    In the second one, the gas particles tagged as "to be swallowed" are effectively swallowed.
    In the third one, the sink particles tagged as "to be swallowed" are effectively swallowed.
    This was done with SWIFT v0.9.0.


Conversion from comoving to physical space
------------------------------------------

In the following, we always refer to physical quantities. In non-cosmological quantities, there is no ambiguity since there is non expansion of the universe to account for and thus no scale factor :math:`a(t)`. However, in cosmological simulation, we need to take care to convert from comoving quantities to physical ones when needed, e.g. to compute energies. Here is a recap:

* :math:`\mathbf{x}_p = \mathbf{x}_ca`
* :math:`\mathbf{v}_p = \mathbf{v}_c/a + a H \mathbf{x}_c`
* :math:`\rho_p = \rho_c/a^3`
* :math:`\Phi_p = \Phi_c/a + c(a)`

Here, the subscript `p` stands for physical and `c` for comoving. 

Notice that the potential normalization constant has been chosen to be :math:`c(a) = 0`. 


Sink formation
--------------

At the core of the sink particle method is the sink formation algorithm. This is critical to form sink in regions adequate to star formation. Failling to can produce spurious sinks and stars, which is not desirable. However, there is not easy answer to the question. We made the choice to implement a simple and efficient algorithm.
The primary criteria required to transform a gas particle into a sink are:

1. the density of a given particle :math:`i` is bigger than a user-defined threshold density: :math:`\rho_i > \rho_{\text{threshold}}` ;
2. the temperature of a given particle is smaller than a user-defined threshold temperature: :math:`T_i < T_{\text{threshold}}`. 

The first critetion is common but not the second one. This is checked to ensure that sink particles, and thus stars, are not generated in hot regions. The parameters for those threshold quantities are respectively called ``density_threshold`` and ``maximal_temperature``.

Then, further criteria are checked. They are always checked for gas particles within the accretion radius :math:`r_{\text{acc}}` (also called cut-off radius) of a given gas particle :math:`i`. Such gas particles are called *neighbours*.

.. note::
   Notice that in the current implementation, the accretion radius is kept *fixed and the same* for all sink. However, for the sake of generality, the mathematical expressions are given as if the accretion radii could be different. 

So, the other criteria are the following:

3. The gas particle is at local potential minimum: :math:`\Phi_i = \min_j \Phi_j`.
4. Gas surrounding the particle is at rest or collapsing: :math:`\nabla \cdot \mathbf{v}_p \leq 0`.
5. The smoothing lenght of the particle is less than half the accretion radius: :math:`h_i < r_{\text{acc}} / 2`.
6. All neighbours are currently active.
7. The sum of thermal of the neighbours satisfies: :math:`E_{\text{therm}} < |E_{\text{pot}}|/2`.
8. The sum of thermal energy and rotational energy satisies: :math:`E_{\text{therm}} + E_{\text{rot}} < | E_{\text{pot}}|`.
9. The total energy of the neihbours si negative, i.e. the clump is bound to the sink: :math:`E_{\text{tot}} < 0`.
10. Forming a sink here will not overlap an existing sink :math:`s`: :math:`\left| \mathbf{x}_i - \mathbf{x}_s \right| > r_{\text{acc}, i} + r_{\text{acc}, s}`.


The different energies are computed as follow:

* :math:`E_{\text{therm}} = \displaystyle \sum_j m_j u_{j, p}`
* :math:`E_{\text{kin}} = \displaystyle \frac{1}{2} \sum_j m_j v_{j, p}^2`
* :math:`E_{\text{pot}} = \displaystyle \sum_j m_j * pi->sink_data.potential_p`
* :math:`E_{\text{rot}} = \displaystyle \sqrt{E_{\text{rot}, x}^2 + E_{\text{rot}, y}^2 + E_{\text{rot}, z}^2}`
* :math:`E_{\text{rot}, x} = \displaystyle \frac{1}{2} \sum_j xxx`
* :math:`E_{\text{rot}, y} = \displaystyle \frac{1}{2} \sum_j xxx`
* :math:`E_{\text{rot}, z} = \displaystyle \frac{1}{2} \sum_j xxx`
* ANGULAR MOMENTUM
* :math:`E_{\text{mag}} = \displaystyle E_{\text{mag}, j}`
* :math:`E_{\text{tot}} = E_{\text{kin}} + E_{\text{pot}} +  E_{\text{therm}} + E_{\text{mag}}`

.. note::
   Currently, magnetic energy is not included in the total energy, since the MHD scheme is in progress. However, the necessary modifications have already been taken care of.

   The :math:`p` subscript is to recall that we are using physical quantities to compute energies.


Some comments about the criteria:

The third criteria is mainly here to prevent two sink particles to form at a distance smaller than the sink accretion radius. Since we allow sinks merging, such situation raises the question of which sink should swallow the other one? This can depend on the order of the task, which is not desirable. As a result, this criterion is enforced.

The last criterion prevents the formation of spurious sinks. Experiences have shown that removing gas within the accretion radius biases the hydro density estimates: the gas feel a force toward the sink. At some point, there is an equilibrium and gas particles accumulate at the edge of the accretion radius, which can then spawn sink particles that do not fall onto the primary sink and thus never merges. This criterion can be disabled. 

.. note::
  Notice however than contrary to  `Bate et al. <https://ui.adsabs.harvard.edu/abs/1995MNRAS.277..362B/abstract>`_, no boundary conditions for sink particles are introduced in the hydrodynamics calculations.



Gas accretion
-------------




Sink merging
------------


IMF sampling
------------

Star spawning
-------------


Stellar feedback
----------------

Stellar feedback per se is not in the sink module, but in the feedback one. However, if one uses sink particles with individual stars, the feedback implementation must be adapted. Here is a recap of the GEAR feedback with sink particles. 
