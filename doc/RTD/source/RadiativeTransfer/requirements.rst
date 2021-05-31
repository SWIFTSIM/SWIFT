.. Radiative Transfer Scheme Requirements
    Mladen Ivkovic 05.2021

.. _rt_requirements:
   
Requirements
------------

.. warning::
    The radiative transfer schemes are still in development and are not useable
    at this moment. This page is currently a placeholder to document new
    features and requirements as the code grows.


To be able to run with radiative transfer, you'll need to run swift with the
following flags:

.. code::

    swift --radiation --hydro --feedback --stars --self-gravity [and/or] --external-gravity


Some notes:

- The radiation data is coupled to the gas data, so you can't run without 
  ``--hydro``. (Also the whole point of these schemes is the interaction between
  gas and radiation, so why would you want to?)

- Currently the only source of radiation are stars, so you need ``--stars``. 
  If you want to run without radiative sources, still run the code with
  ``--stars``, even if you don't have any stars in your initial conditions.

- Running with ``--stars`` requires some form of gravity, be it self-gravity or
  external gravity. Since we need stars, we inherit this requirement. If you want
  no gravity, run with ``--external-gravity`` and set the external potential to
  zero.

- We need ``--feedback`` in order to have meaningful smoothing lengths for
  stars. However, you don't need any specific feedback model; It'll work even if
  you configured ``--with-feedback=none``.


