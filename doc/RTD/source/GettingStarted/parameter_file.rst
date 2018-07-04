.. Parameter File
   Loic Hausammann, 1 june 2018

.. _Parameter_File_label:

Parameter File
==============

To run SWIFT, you will need to provide a ``yaml`` parameter file.  An example is
given in ``examples/parameter_file.yml`` which should contain all possible
parameters.  Each section in this file corresponds to a different option in
SWIFT and are not always required depending on the configuration options and
the run time parameters.


Output Selection
~~~~~~~~~~~~~~~~

With SWIFT, you can select the particle fields to output in snapshot using the parameter file.
In section ``SelectOutput``, you can remove a field by adding a parameter formatted in the
following way ``field_parttype`` where ``field`` is the name of the field that you
want to remove (e.g. ``Masses``) and ``parttype`` is the type of particles that
contains this field (e.g. ``Gas``, ``DM`` or ``Star``).  For a parameter, the only
values accepted are 0 (skip this field when writing) or 1 (default, do not skip
this field when writing).

You can generate a ``yaml`` file containing all the possible fields with ``./swift -o output.yml``. By default, all the fields are written.
