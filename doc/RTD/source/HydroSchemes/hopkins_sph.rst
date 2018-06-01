.. 'Hopkins'-SPH
   Josh Borrow 5th April 2018

Pressure-Entropy SPH
====================

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Contents:

A pressure-entropy SPH scheme is available in SWIFT, inspired by Hopkins 2013.
This includes a Monaghan AV scheme and a Balsara switch.


.. code-block:: bash
   
   ./configure --with-hydro="pressure-entropy"


Pressure-Energy SPH
===================

Pressure-energy SPH is now implemented in SWIFT, and like the pressure-entropy
scheme it includes a Monaghan AV scheme and a Balsara switch.


.. code-block:: bash
   
   ./configure --with-hydro="pressure-energy"
