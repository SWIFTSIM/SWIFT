.. Adding new schemes
   Darwin Roduit, 16 Ocotber 2024

.. _new_option_sink:

How to add your sink scheme
-------------------------------

Here, we provide comprehensive information to guide you in adding your sink scheme into Swift. To better understand how to add new schemes within Swift, read the general information provided on :ref:`new_option` page. 

The default sink scheme is empty and gives you an idea of the minimum required fields and functions for the code to compile. The GEAR sink module has the base functions plus some extra ones for its operations. It can provide you with a working example. However, it can only work with the GEAR feedback module since it relies on IMF properties that are only located there. 

As discussed in the GEAR sink :ref:`sink_GEAR_model_summary`, the physics relies on the following tasks: sink formation, gas and sink particle flagging, gas swallowing, sink swallowing and star formation. You do not need to care about the tasks, only the core functions within the sink module. However, you may need to locate where the code calls these functions. The file ``src/runner_others.c`` contains the ``runner_do_star_formation_sink()`` and ``runner_do_sink_formation()``. These functions are responsible for generating stars out of sinks and sinks from gas particles. The other general task-related functions are in ``src/runner_sinks.c``.

The following presents the most essential functions you need to implement. This will give you an idea of the workload. 


Sink formation
~~~~~~~~~~~~~~

Before forming a sink, the code loops over all gas and sink particles to gather data about its neighbourhood. This is performed in ``sink_prepare_part_sink_formation_gas_criteria()`` and ``sink_prepare_part_sink_formation_sink_criteria()`` functions. For instance, in GEAR, we compute the total potential energy, thermal energy, etc. 

Then, to decide if we can turn a gas particle into a sink particle, the function ``sink_is_forming()`` is called. Before forming a sink particle, there is a call to ``sink_should_convert_to_sink()``. This function determines whether the gas particle must transform into a sink. Both functions return either 0 or 1.

Gas and sink flagging: finding whom to eat
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Before accreting the gas/sink particles, the sink needs to look for eligible particles. The gas swallow interactions are performed within ``runner_iact_nonsym_sinks_gas_swallow()`` and the sink swallow in ``runner_iact_nonsym_sinks_sink_swallow()``.


Gas and sink swallowing: updating the sink properties
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When the sink swallows gas particles, it updates its internal properties based on the gas particles' properties. The ``sink_swallow_part()`` function takes care of this.

Similarly, when the sink swallows sink particles, it updates its properties from the to-be-swallowed sink particles. The ``sink_swallow_sink()`` performs the update.

There is no more than that. The code properly removes the swallowed particles. 

Star formation
~~~~~~~~~~~~~~

The most important function is ``sink_spawn_star()``. It controls whether the code should continue creating stars and must return 0 or 1.

The ``sink_copy_properties_to_star()`` does what it says. This function is also responsible for adequately initialising the stars' properties. In GEAR, we give new positions and velocities within this function. 

The following three functions allow you to update your sink particles (e.g. their masses) before, during and after the star formation loop: ``sink_update_sink_properties_before_star_formation()``, ``sink_update_sink_properties_during_star_formation()`` and ``sink_update_sink_properties_after_star_formation()``.

These functions are located in ``sink/Default/sink.h``
