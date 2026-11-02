.. FLAMINGO sub-grid model
   Matthieu Schaller, 24th March 2026


FLAMINGO model
==============

This section of the documentation gives a brief description of the
different components of the FLAMINGO sub-grid model. The model
is described by `Schaye et al. (2023) <https://ui.adsabs.harvard.edu/abs/2023MNRAS.526.4978S/abstract>`_
and the procedure used to calibrate the free parameters was presented
by `Kugel et al. (2023) <https://ui.adsabs.harvard.edu/abs/2023MNRAS.526.6103K/abstract>`_.


At the level of the code, all the components of the FLAMINGO model are
directly taken from the :ref:`EAGLE` model as implemented in
SWIFT. The differences lie in the choice of runtime parameters
used. The parameter files resulting from the calibration effort are
given in the examples directory.
