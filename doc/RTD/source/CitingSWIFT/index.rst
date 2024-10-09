.. Citing SWIFT
   Matthieu Schaller, 5th April 2020

Disclaimer, Citing SWIFT & Giving Credit
========================================

First of all, thank you for using SWIFT for your projects!

Licensing
~~~~~~~~~

SWIFT is licensed under the LGPL version 3.0. A copy of the license is provided
in the file ``COPYING`` at the root of the repository. It can also be found
online `here <https://www.gnu.org/licenses/lgpl-3.0-standalone.html>`_ as an
extension to the `GPL <https://www.gnu.org/licenses/gpl-3.0-standalone.html>`_.

.. figure:: https://www.gnu.org/graphics/lgplv3-with-text-154x68.png
    :width: 100px
    :alt: LGPL v3.0 logo

Disclaimer
~~~~~~~~~~

We would like to emphasise that SWIFT comes without any warranty of accuracy,
correctness or efficiency. As mentioned in the license, the software comes
`as-is` and the onus is on the user to get meaningful results. Whilst the
authors will endeavour to answer questions related to using the code, we
recommend users build and maintain their own copies. This documentation contains
the most basic information to get started. Reading it and possibly also the
source code is the best way to start running simulations.

The users are responsible to understand what the code is doing and for the
results of their simulation runs.

Note also that the values of the parameters given in the examples are only
indicative. We recommend users experiment by themselves and a campaign of
experimentation with various values is highly encouraged. Each problem will
likely require different values and the sensitivity to the details of the
physical model is something left to the users to explore.

Acknowledgment & Citation
~~~~~~~~~~~~~~~~~~~~~~~~~

The SWIFT code was last described in this `publication
<https://ui.adsabs.harvard.edu/abs/2023arXiv230513380S/abstract>`_, which
introduced the core solver, the numerical methods used as well as many
extensions. We ask users running SWIFT for their research to please cite this
paper when they present their results. This paper is best referenced by the
following bibtex citation block:

.. code-block:: bibtex

  @ARTICLE{2023arXiv230513380S,
    author = {{Schaller}, Matthieu and others},
    title = "{SWIFT: A modern highly-parallel gravity and smoothed particle hydrodynamics solver for astrophysical and cosmological applications}",
    journal = {\mnras},
    keywords = {software: simulations, methods: numerical, software: public release, Astrophysics - Instrumentation and Methods for Astrophysics, Astrophysics - Cosmology and Nongalactic Astrophysics, Astrophysics - Earth and Planetary Astrophysics, Astrophysics - Astrophysics of Galaxies, Computer Science - Distributed, Parallel, and Cluster Computing},
    year = 2024,
    month = may,
    volume = {530},
    number = {2},
    pages = {2378-2419},
    doi = {10.1093/mnras/stae922},
    archivePrefix = {arXiv},
    eprint = {2305.13380},
    primaryClass = {astro-ph.IM},
    adsurl = {https://ui.adsabs.harvard.edu/abs/2024MNRAS.530.2378S},
    adsnote = {Provided by the SAO/NASA Astrophysics Data System}
  }


In order to keep track of usage and measure the impact of the software, we
kindly ask users publishing scientific results using SWIFT to add the following
sentence to the acknowledgment section of their papers:

.. code-block:: text
		
   The research in this paper made use of the SWIFT open-source simulation code
   (http://www.swiftsim.com, Schaller et al. 2018) version X.Y.Z.
   
with the version number set to the version used for the simulations and the
reference pointing to the `ASCL entry <https://ascl.net/1805.020>`_ of the
code. This corresponds to the following bibtex citation block:

.. code-block:: bibtex

   @MISC{2018ascl.soft05020S,
     author = {{Schaller}, Matthieu and others},
     title = "{SWIFT: SPH With Inter-dependent Fine-grained Tasking}",
     keywords = {Software },
     howpublished = {Astrophysics Source Code Library},
     year = 2018,
     month = may,
     eid = {ascl:1805.020},
     pages = {ascl:1805.020},
     archivePrefix = {ascl},
     eprint = {1805.020},
     adsurl = {https://ui.adsabs.harvard.edu/abs/2018ascl.soft05020S},
     adsnote = {Provided by the SAO/NASA Astrophysics Data System}
   }

When using models or parts of the code whose details were introduced in other
papers, we kindly ask that the relevant work is properly acknowledge and
cited. This includes the :ref:`subgrid`, the :ref:`planetary` extensions, the
:ref:`hydro` and :ref:`rt`, or the particle-based :ref:`neutrinos`.
