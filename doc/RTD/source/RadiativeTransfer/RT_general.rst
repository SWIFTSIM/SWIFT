.. Radiative Transfer Scheme Requirements
    Mladen Ivkovic 05.2021

.. _rt_general:


.. warning::
    The radiative transfer schemes are still in development and are not useable
    at this moment. This page is currently a placeholder to document new
    features and requirements as the code grows.




Radiative Transfer in SWIFT
---------------------------

We currently support two methods for radiative transfer in SWIFT:
:ref:`GEAR-RT <rt_GEAR>` and :ref:`SPHM1RT <rt_SPHM1>`. They are both
moment-based methods. For documentation on each of these methods, please
refer to their respective pages.

This page contains some general documentation on running SWIFT with RT, which
is valid irrespective of the exact method used.





Requirements
~~~~~~~~~~~~


To be able to run with radiative transfer, you'll need to run swift with the
following flags:

.. code::

    swift --radiation --hydro --feedback --stars --self-gravity [and/or] --external-gravity



Some notes on these runtime flags:


- The radiation data is coupled to the gas data, so you can't run without
  ``--hydro``. (Also, the whole point of these schemes is the interaction between
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





.. _rt_subcycling:

RT Sub-Cycling
~~~~~~~~~~~~~~

SWIFT allows to sub-cycle the solution of radiative transfer steps (both
photon propagation and thermochemistry) with respect to the hydrodynamics
time steps. Basically you can tell SWIFT to run up to X radiative transfer
steps during a single hydrodynamics step for all particles in the simulation.
The aim is to not waste time doing unnecessary hydrodynamics updates, which
typically allow for much higher time steps compared to radiation due to the
propagation speed of the respective advected quantity.

You will need to provide an upper limit on how many RT sub-cycles per hydro
step you want to allow. That is governed by the

.. code:: yaml

   TimeIntegration:
       max_nr_rt_subcycles: 128         # maximal number of RT sub-cycles per hydro step

parameter, which is mandatory for any RT runs. To turn off sub-cycling and
couple the radiative transfer and the hydrodynamics time steps one-to-one,
set this parameter to either 0 or 1.

Due to the discretization of individual particle time steps in time bins
with a factor of 2 difference in time step size from a lower to a higher
time bin, the ``max_nr_rt_subcycles`` parameter itself is required to be
a power of 2 as well.

Note that this parameter will set an upper limit to the number of sub-cycles
per hydro step. If the ratio of hydro-to-RT time step is greater than what
``max_nr_rt_subcycles`` allows for, then the hydro time step will be reduced
to fit the maximal threshold. If it is smaller, the particle will simply do
fewer sub-cycles.







Reading the output of time-step information
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


When running with :ref:`sub-cycling enabled <rt_subcycling>`, additional time
step information is written.

Firstly, whenever you run with a SWIFT compiled with any RT method, it will also
create a file called ``rtsubcycles.txt``, which  contains analogous data to the
``timesteps.txt`` file. However, if you run without sub-cycling (i.e. with
``TimeIntegration:max_nr_rt_subcycles`` set to ``0`` or ``1``), **then the file
will remain empty** (save for the header information). That's because its contents
would be redundant - all the data will be identical to what will be written in
``timesteps.txt``.



Secondly, when running with RT and sub-cycling enabled, the output to ``STDOUT``
contains additional information. More concretely, it changes from something like this:


.. code-block::

    [00003.0] engine_dump_snapshot: Dumping snapshot at t=0.000000e+00
    [00003.1] engine_print_stats: Saving statistics at t=0.000000e+00
    #   Step           Time Scale-factor     Redshift      Time-step Time-bins      Updates    g-Updates    s-Updates sink-Updates    b-Updates  Wall-clock time [ms]  Props    Dead time [ms]
           0   0.000000e+00    1.0000000    0.0000000   0.000000e+00    1   56        18000        31000        13000            0            0              2610.609    281           101.971
           1   1.220703e-05    1.0000000    0.0000000   1.220703e-05   43   43           11           11            0            0            0                61.686      1             1.324
           2   2.441406e-05    1.0000000    0.0000000   1.220703e-05   43   44        12685        12685            0            0            0              1043.433      0            35.461
           3   3.662109e-05    1.0000000    0.0000000   1.220703e-05   43   43           11           11            0            0            0                51.340      1             1.628
           4   4.882813e-05    1.0000000    0.0000000   1.220703e-05   43   45        18000        18000            0            0            0              1342.531      0            36.831
           5   6.103516e-05    1.0000000    0.0000000   1.220703e-05   43   43           11           11            0            0            0                48.412      1             1.325
           6   7.324219e-05    1.0000000    0.0000000   1.220703e-05   43   44        12685        12685            0            0            0              1037.307      0            34.718
           7   8.544922e-05    1.0000000    0.0000000   1.220703e-05   43   43           11           11            0            0            0                47.791      1             1.362
           8   9.765625e-05    1.0000000    0.0000000   1.220703e-05   43   46        18000        18004            4            0            0              1410.851      0            35.005
           9   1.098633e-04    1.0000000    0.0000000   1.220703e-05   43   43           11           11            0            0            0                48.322      1             1.327
          10   1.220703e-04    1.0000000    0.0000000   1.220703e-05   43   44        12685        12685            0            0            0              1109.944      0            33.691


To something like this:


.. code-block::

     [rt-sc] 0    0.000000e+00    1.000000    0.000000  1.220703e-05    1   56        18000            -            -            -            -                     -      -                 -
     [rt-sc] 1    1.220703e-05    1.000000    0.000000  1.220703e-05   43   43           11            -            -            -            -                 0.683      -             0.170
     [rt-sc] 2    2.441406e-05    1.000000    0.000000  1.220703e-05   43   44        12685            -            -            -            -               100.543      -             5.070
     [rt-sc] 3    3.662109e-05    1.000000    0.000000  1.220703e-05   43   43           11            -            -            -            -                 0.762      -             0.208
     [rt-sc] 4    4.882813e-05    1.000000    0.000000  1.220703e-05   43   45        18000            -            -            -            -               124.011      -             6.396
     [rt-sc] 5    6.103516e-05    1.000000    0.000000  1.220703e-05   43   43           11            -            -            -            -                 0.831      -             0.254
     [rt-sc] 6    7.324219e-05    1.000000    0.000000  1.220703e-05   43   44        12685            -            -            -            -               107.172      -             5.572
     [rt-sc] 7    8.544922e-05    1.000000    0.000000  1.220703e-05   43   43           11            -            -            -            -                 0.759      -             0.227
    [00003.0] engine_dump_snapshot: Dumping snapshot at t=0.000000e+00
    [00003.1] engine_print_stats: Saving statistics at t=0.000000e+00
    #   Step           Time Scale-factor     Redshift      Time-step Time-bins      Updates    g-Updates    s-Updates sink-Updates    b-Updates  Wall-clock time [ms]  Props    Dead time [ms]
           0   0.000000e+00    1.0000000    0.0000000   0.000000e+00    1   56        18000        31000        13000            0            0              2941.254    281           120.261
     [rt-sc] 0    9.765625e-05    1.000000    0.000000  1.220703e-05   43   46        18000            -            -            -            -                     -      -                 -
     [rt-sc] 1    1.098633e-04    1.000000    0.000000  1.220703e-05   43   43           11            -            -            -            -                 0.990      -             0.417
     [rt-sc] 2    1.220703e-04    1.000000    0.000000  1.220703e-05   43   44        12685            -            -            -            -               104.155      -             5.744
     [rt-sc] 3    1.342773e-04    1.000000    0.000000  1.220703e-05   43   43           11            -            -            -            -                 0.765      -             0.176
     [rt-sc] 4    1.464844e-04    1.000000    0.000000  1.220703e-05   43   45        18000            -            -            -            -               125.237      -             5.605
     [rt-sc] 5    1.586914e-04    1.000000    0.000000  1.220703e-05   43   43           11            -            -            -            -                 0.856      -             0.282
     [rt-sc] 6    1.708984e-04    1.000000    0.000000  1.220703e-05   43   44        12685            -            -            -            -               112.171      -             5.251
     [rt-sc] 7    1.831055e-04    1.000000    0.000000  1.220703e-05   43   43           11            -            -            -            -                 0.861      -             0.241
           1   9.765625e-05    1.0000000    0.0000000   9.765625e-05   46   46            4            8            4            0            0               546.225      1            24.648
     [rt-sc] 0    1.953125e-04    1.000000    0.000000  1.220703e-05   43   47        18000            -            -            -            -                     -      -                 -
     [rt-sc] 1    2.075195e-04    1.000000    0.000000  1.220703e-05   43   43           11            -            -            -            -                 0.842      -             0.212
     [rt-sc] 2    2.197266e-04    1.000000    0.000000  1.220703e-05   43   44        12685            -            -            -            -               126.674      -             6.295
     [rt-sc] 3    2.319336e-04    1.000000    0.000000  1.220703e-05   43   43           11            -            -            -            -                 0.797      -             0.289
     [rt-sc] 4    2.441406e-04    1.000000    0.000000  1.220703e-05   43   45        18000            -            -            -            -               142.086      -             5.511
     [rt-sc] 5    2.563477e-04    1.000000    0.000000  1.220703e-05   43   43           11            -            -            -            -                 0.919      -             0.196
     [rt-sc] 6    2.685547e-04    1.000000    0.000000  1.220703e-05   43   44        12685            -            -            -            -               131.550      -             5.896
     [rt-sc] 7    2.807617e-04    1.000000    0.000000  1.220703e-05   43   43           11            -            -            -            -                 0.809      -             0.186
           2   1.953125e-04    1.0000000    0.0000000   9.765625e-05   46   47           27           43           16            0            0               558.226      0            27.711
     [rt-sc] 0    2.929688e-04    1.000000    0.000000  1.220703e-05   43   46        18000            -            -            -            -                     -      -                 -
     [rt-sc] 1    3.051758e-04    1.000000    0.000000  1.220703e-05   43   43           11            -            -            -            -                 0.738      -             0.207
     [rt-sc] 2    3.173828e-04    1.000000    0.000000  1.220703e-05   43   44        12685            -            -            -            -               122.572      -             5.170
     [rt-sc] 3    3.295898e-04    1.000000    0.000000  1.220703e-05   43   43           11            -            -            -            -                 1.063      -             0.345
     [rt-sc] 4    3.417969e-04    1.000000    0.000000  1.220703e-05   43   45        18000            -            -            -            -               147.110      -             5.409
     [rt-sc] 5    3.540039e-04    1.000000    0.000000  1.220703e-05   43   43           11            -            -            -            -                 1.091      -             0.350
     [rt-sc] 6    3.662109e-04    1.000000    0.000000  1.220703e-05   43   44        12685            -            -            -            -               134.273      -             6.561
     [rt-sc] 7    3.784180e-04    1.000000    0.000000  1.220703e-05   43   43           11            -            -            -            -                 0.825      -             0.298
           3   2.929688e-04    1.0000000    0.0000000   9.765625e-05   46   46            4            8            4            0            0               557.164      0            24.760


Here's what's going on here:

- All lines beginning with the prefix ``[rt-sc]`` are time step data of RT
  sub-cycling steps, i.e. of sub-cycles.
- The sub-cycling index follows the prefix ``[rt-sc]``. For example, a line
  beginning with ``[rt-sc] 2`` is the sub-cycle with index ``2`` of some time step.
- The "sub-cycle" with index ``0`` is the one performed alongside all other tasks
  (e.g. hydro, gravity, stellar feedback, etc.) All other sub-cycle indices
  indicate actual sub-cycles, i.e. actual intermediate steps where only radiative
  transfer is being solved (or put differently: where only RT tasks are being launched).
- The sub-cycling lines are written to ``STDOUT`` *before* the line of the full
  time step data. More precisely, in the above example, time step ``1`` with all its
  RT sub-cycles is written to ``STDOUT`` as this block:

.. code-block::

    #   Step           Time Scale-factor     Redshift      Time-step Time-bins      Updates    g-Updates    s-Updates sink-Updates    b-Updates  Wall-clock time [ms]  Props    Dead time [ms]

    [ ... some lines omitted ... ]

     [rt-sc] 0    9.765625e-05    1.000000    0.000000  1.220703e-05   43   46        18000            -            -            -            -                     -      -                 -
     [rt-sc] 1    1.098633e-04    1.000000    0.000000  1.220703e-05   43   43           11            -            -            -            -                 0.990      -             0.417
     [rt-sc] 2    1.220703e-04    1.000000    0.000000  1.220703e-05   43   44        12685            -            -            -            -               104.155      -             5.744
     [rt-sc] 3    1.342773e-04    1.000000    0.000000  1.220703e-05   43   43           11            -            -            -            -                 0.765      -             0.176
     [rt-sc] 4    1.464844e-04    1.000000    0.000000  1.220703e-05   43   45        18000            -            -            -            -               125.237      -             5.605
     [rt-sc] 5    1.586914e-04    1.000000    0.000000  1.220703e-05   43   43           11            -            -            -            -                 0.856      -             0.282
     [rt-sc] 6    1.708984e-04    1.000000    0.000000  1.220703e-05   43   44        12685            -            -            -            -               112.171      -             5.251
     [rt-sc] 7    1.831055e-04    1.000000    0.000000  1.220703e-05   43   43           11            -            -            -            -                 0.861      -             0.241
           1   9.765625e-05    1.0000000    0.0000000   9.765625e-05   46   46            4            8            4            0            0               546.225      1            24.648

- Let's have a closer look at the written data:

  - In step ``1``, 4 hydro particles, 8 gravity particles, and 4 star particles
    were updated. (You can see that in the last line.)
  - The integration of the full step was performed over a time step with
    size ``9.765625e-05``. **That is valid for all physics except radiative
    transfer.** The RT was integrated 8 times with a time step size of ``1.220703e-05``.
    You can see this in the sub-cycling output lines.
  - Each RT sub-cycling line also tells you the minimal and maximal time bin
    size that was worked on, as well as how many hydro particles underwent RT
    updates.
  - RT sub-cycles only ever update radiative transfer. There will never be any
    gravity, star, sink, or black hole particle updates in it.
  - Since the sub-cycle with index ``0`` is performed alongside all other physics
    during the main step, the isolated wall-clock time and dead time fields are
    not available for it.
  - The wall-clock time and dead time fields in the full step line (the one starting
    *without* ``[rt-sc]``) include the data of the sub-cycles for this step as well.





