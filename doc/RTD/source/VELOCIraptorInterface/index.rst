.. VELOCIraptor Interface
   Folkert Nobels, 8th October 2018

.. _VELOCIraptor_interface:
   
VELOCIraptor Interface
======================

This section includes information on the VELOCIraptor interface implemented in
SWIFT. There are mainly four subsections; the first section explains shortly 
how VELOCIraptor works, the second subsection explains how to configure SWIFT
with VELOCIraptor, the third subsection explains how to configure a standalone
version of VELOCIraptor and the fourth subsection explains how the output format
of VELOCIraptor works. There is an additional fifth section which describes the
optional output of halo most bound particles for use with semi-analytic galaxy
formation models.

.. warning::
   The VELOCIraptor interface is only tested for gravity-only
   simulations and the VELOCIraptor code itself is not as advanced as
   other packages for hydrodynamical simulations. The `HBT-HERONS
   <ui.adsabs.harvard.edu/abs/2025MNRAS.543.1339F>`_ package [#f1]_ is
   the recommended alternative developed by the SWIFT team and used
   for the FLAMINGO and COLIBRE simulation campaigns. Note that this finder
   directly uses the FOF outputs and does not require a code interface
   within SWIFT.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   whatis
   stfwithswift
   stfalone
   output
   orphanparticles


.. rubric:: Footnotes

.. [#f1] See https://hbt-herons.strw.leidenuniv.nl/
