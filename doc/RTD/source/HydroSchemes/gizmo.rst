.. GIZMO (MFV)
   Josh Borrow, 5th April 2018

GIZMO-Like Scheme
=================

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Contents:


There is a mesh-less finite volume (MFV) GIZMO-like scheme implemented in SWIFT
(see Hopkins 2015 for more information). You will need a Riemann solver to run
this, and configure as follows:

.. code-block:: bash
   
   ./configure --with-hydro="gizmo-mfv" --with-riemann-solver="hllc"


We also have the mesh-less finite mass (MFM) GIZMO-like scheme. You can select
this at compile-time with the following configuration flags:

.. code-block:: bash
   
   ./configure --with-hydro="gizmo-mfm" --with-riemann-solver="hllc"
