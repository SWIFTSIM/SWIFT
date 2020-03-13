.. Minimal SPH
   Josh Borrow 4th April 2018

.. _minimal:

Minimal (Density-Energy) SPH
============================

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Contents:

This scheme is a textbook implementation of Density-Energy SPH, and can be used
as a pedagogical example. It also implements a Monaghan AV scheme with a
Balsara switch, like the GADGET-2 scheme. It uses very similar equations, but
differs in implementation details; namely it tracks the internal energy \(u\)
as the thermodynamic variable, rather than entropy \(A\). To use the minimal
scheme, use

.. code-block:: bash

    ./configure --with-hydro="minimal"

As it uses a very simple, fixed artificial viscosity, only the
``SPH:viscosity_alpha`` parameter has any effect for this scheme. This will
change the strength of the artificial viscosity throughout the simulation,
and has a default of 0.8.
