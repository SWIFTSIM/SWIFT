.. Parameter File
   Loic Hausammann, 1 June 2018

.. _Output_list_label:

Output List
~~~~~~~~~~~

In the sections ``Snapshots`` and ``Statistics``, you can specify the
options ``output_list_on`` and ``output_list`` which receive an int
and a filename.  The ``output_list_on`` enable or not the output list
and ``output_list`` is the filename containing the output times.  With
the file header, you can choose between writing redshifts, scale
factors or times.

Example of file containing with times (in internal units)::

  # Time
  0.5
  1.5
  3.0
  12.5

Example of file with scale factors::

  # Scale Factor
  0.1
  0.2
  0.3

Example of file with redshift::

  # Redshift
  20
  15
  10
  5

If an output list is specified, the basic values for the first
snapshot (``time_first``, ``scale_factor_first``) and difference
(``delta_time``) are ignored.
  
.. _Output_selection_label:

Output Selection
~~~~~~~~~~~~~~~~

With SWIFT, you can select the particle fields to output in snapshot
using the parameter file.  In section ``SelectOutput``, you can remove
a field by adding a parameter formatted in the following way
``field_parttype`` where ``field`` is the name of the field that you
want to remove (e.g. ``Masses``) and ``parttype`` is the type of
particles that contains this field (e.g. ``Gas``, ``DM`` or ``Star``).
For a parameter, the only values accepted are 0 (skip this field when
writing) or 1 (default, do not skip this field when writing). By
default all fields are written.

This field is mostly used to remove unnecessary output by listing them
with 0's. A classic use-case for this feature is a DM-only simulation
(pure n-body) where all particles have the same mass. Outputting the
mass field in the snapshots results in extra i/o time and unnecessary
waste of disk space. The corresponding section of the ``yaml``
parameter file would look like this::

  SelectOutput:
    Masses_DM:   0

You can generate a ``yaml`` file containing all the possible fields
available for a given configuration of SWIFT by running ``./swift --output-params output.yml``.
