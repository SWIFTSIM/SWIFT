.. Light Cones
   John Helly 29th April 2021

.. _lightcone_running_label:

Running SWIFT with Light Cone Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To produce light cone particle output swift must be configured
with ``--enable-lightcone``. Additionally, making HEALPix maps
requires the HEALPix C library. If using MPI then parallel HDF5
is also required.

Light cone parameters are specified in the ``Lightcone`` section
of the SWIFT parameter file. See :ref:`Parameters_light_cone` for details.

SWIFT must be run with the ``--lightcone`` flag to activate light
cone outputs, otherwise the Lightcone section in the parameter file
is ignored.



