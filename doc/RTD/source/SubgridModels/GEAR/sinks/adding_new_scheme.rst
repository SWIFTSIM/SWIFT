.. Sink particles in GEAR model
   Darwin Roduit, 14 July 2024

.. sink_GEAR_model:

How to add your own sink scheme
-------------------------------

Here, we provide comprehensive information to guide you in adding your sink scheme into Swift. To better understand how to add new schemes within Swift, read the general information provided on :ref:`new_option` page. 

The default sink scheme is empty and gives you an idea of the minimum required fields and functions for the code to compile. The GEAR sink module has the base functions plus some extra ones for its operations within the module. However, it can only work with the GEAR feedback module because it relies on IMF properties that are only located there. 

As discussed in the model summary, the physics relies on the following tasks: sink formation, gas and sink particle flagging, gas swallowing, sink swallowing and star formation. You do not need to care about the tasks, only the core functions within the sink module. However, you may need to locate where the code calls these functions. The file ``src/runner_others.c`` contains the ``runner_do_star_formation_sink()``. This function is responsible for generating stars out of sinks. All other functions are in ``src/runner_sinks.c``.
