.. GIZMO (MFV)
   Josh Borrow, 5th April 2018

GIZMO-Like Scheme
=================

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Contents:


There is a meshless finite volume (MFV) GIZMO-like scheme implemented in SWIFT
(see Hopkins 2015 for more information). You will need a Riemann solver to run
this, and configure as follows:

.. code-block:: bash
   
   ./configure --with-hydro="gizmo" --with-riemann-solver="hllc" --disable-vec

