.. Sink particles in GEAR model
   Darwin Roduit, 15 March 2024

.. sink_GEAR_model:

Introduction
------------

GEAR sink particles provide an alternative model to star formation. Instead of stochastically transforming gas particles into stars as is done in GEAR star formation scheme under some conditions, we transform a gas into a sink particle. The main property of the sink particle is its accretion radius. When gas particles within this accretion radius are eligible to be swallowed by the sink, we remove them and transfer their mass, momentum, angular momentum, chemistry properties, etc to the sink particle.

With the sink particles, the IMF splits into two parts: the continuous part and the discrete part. Those parts will correspond to two kinds of stars. Particles in the discrete part of the IMF represent individual stars. It means that discrete IMF-sampled stars have different masses. Particles in the continuous part represent a population of stars, all with the same mass.

The sink particle will randomly choose a target mass, accrete gas until it reaches this target mass and finally spawn a star. Then, the sink chooses a new target mass and repeats the same procedure. When stars are spawned, they are given a new position and velocity as well as chemical properties.

In ``theory.rst`` we outline all of the theory which is implemented as part of the model. This includes how sink particles are formed, how the gas is accreted, how sinks are merged and how stars are spawned. In ``params.rst`` we list and discuss all parameters used by the model. Below we outline how to configure and run the model.

Compiling and running the model
-------------------------------

You can configure the model with ``--with-sink=GEAR`` in combination with other configure options of the GEAR model. The model will then be used when the ``--sinks`` flag is among the runtime options.

Then, you do not need to do anything special. Sink particles will be created during your runs. If you want, you can have sink particles in your ICs. At the moment, sink particles do not have any special fields to be set.

A full list of all relevant parameters of the model is in :ref:`sink_GEAR_parameters`. We also briefly describe the most important parameters which need to be set to run the model, as well as how to run it in different configurations.

