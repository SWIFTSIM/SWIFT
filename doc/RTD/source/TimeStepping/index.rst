.. Details of the time-stepping scheme
   Matthieu Schaller 9th November 2019

.. _time_stepping:

Time integration details
========================

This section of the documentation includes information on the time
stepping and time integration used in SWIFT. It starts with the integer
time-line and the time-bins, then describes how the time-step of a particle
is chosen and how the kick-drift-kick scheme runs as tasks. The time-step
limiter and the synchronisation are described next. The last page explains the
table of the steps that SWIFT prints, and what the time-step parameters cost.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   integer_time_line
   time_step_criteria
   kick_drift_kick
   timestep_limiter
   timestep_synchronization
   tuning

The figures
-----------

The figures of these pages are not stored in the repository. They are made by
the Python scripts in the directory ``TimeStepping/figures`` every time the
documentation is built. The build (``conf.py``) runs all the ``*.py`` files
under ``doc/RTD/source`` before Sphinx reads the pages, and it stops with an
error if one of the scripts fails. Every script writes its figure next to
itself and does nothing if the figure already exists, so a build in a
directory that already has the figures does not draw them again. To draw a
figure again, delete it and build the documentation, or run its script from
that directory, for instance ``python3 plot_time_bins.py``. The scripts need
``numpy`` and ``matplotlib``.

* ``timeline_helpers.py`` contains the functions of ``src/timeline.h`` and the
  rules of ``make_integer_timestep()`` written in Python, so that the figures
  use the definitions of the code. It draws nothing.
* ``plot_time_bins.py``, ``plot_active_bins.py``, ``plot_kdk_timeline.py`` and
  ``plot_limiter_interrupt.py`` draw the figures that explain the time-line, the
  scheme and the limiter. They need no input.
* ``plot_sedov_steps.py`` draws the two figures of the real run. It reads
  ``sedov_3d_timesteps.txt``, which is the ``timesteps.txt`` of the
  ``SedovBlast_3D`` example (:math:`64^3` glass, default parameter file,
  ``./swift --hydro --limiter --threads=4 sedov.yml``). This file is kept in the
  repository, so that the figures can be drawn without running SWIFT.
