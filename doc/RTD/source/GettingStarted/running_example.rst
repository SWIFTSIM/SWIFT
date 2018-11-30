.. Running an Example
   Josh Borrow, 5th April 2018

Running an Example
==================

Now that you have built the code, you will want to run an example! To do that,
you need to follow the following instructions (requires ``python2`` or
``python3`` with the ``h5py`` and other standard scientific packages, as well
as ``wget`` for grabbing the glass).

.. code-block:: bash
   
   cd examples/SodShock_3D
   ./getGlass.sh
   python makeIC.py
   ../swift --hydro --threads=4 sodShock.yml
   python plotSolution.py 1


This will run the 'SodShock' in 3D and produce a nice plot that shows you
how the density has varied. Try running with GIZMO-MFV (this will take
_significantly_ longer than with SPH) to see the difference. For that, you
will need to reconfigure with the following options:

.. code-block:: bash
   
   ./configure \
   --with-hydro=gizmo-mfv \
   --with-riemann-solver=hllc


To see the results that you should get, you should check out our developer
wiki at https://gitlab.cosma.dur.ac.uk/swift/swiftsim/wikis/Sod_3D.

If you don't get these results, please contact us on our GitHub page at
https://github.com/SWIFTSIM/swiftsim/issues.
