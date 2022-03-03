.. Light Cones
   John Helly 29th April 2021

.. _lightcone_adding_outputs_label:

Adding New Types of Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~

New particle properties can be added to the particle light cones as follows:

* Add a field to the ``lightcone_<type>_data`` struct in ``lightcone_particle_io.h`` to store the new quantity
* Modify the ``lightcone_store_<type>`` function in ``lightcone_particle_io.c`` to set the new struct field from the particle data
* in ``lightcone_io_make_output_fields()``, add a call to ``lightcone_io_make_output_field()`` to define the new output

Here, <type> is the particle type: gas, dark_matter, stars, black_hole or neutrino.

To add a new type of HEALPIX map:

* Add a function to compute the quantity in ``lightcone_map_types.c``. See ``lightcone_map_total_mass()`` for an example.
* Add a new entry to the ``lightcone_map_types`` array in lightcone_map_types.h. This should specify the name of the new map type, a pointer to the function to compute the quantity, and the units of the quantity. The last entry in the array is not used and must have a NULL function pointer to act as an end marker.
