Continuous Simulation Data Stream (CSDS)
========================================

   .. warning::

      CSDS is an emerging feature and is subject to ongoing development. We encourage users to report issues and ask questions on the SWIFT Slack for support."

The Continuous Simulation Data Stream (CSDS) is a new output strategy to replace or supplement the traditional snapshots.

Traditional snapshot work by writing particles, either all of them or a subset, at fixed user-specified time intervals. The drawback of this approach is the following. In regions where particles evolve a little, we redundantly sample these particles. In regions where particles evolve a lot due to the local dynamics (gravity, hydrodynamics, feedback, etc), we under-sample them and cannot accurately track their evolution. This phenomenon is particularly evident in cosmological simulations where particles evolve on a wide range of different timescales.

To circumvent this problem, the CSDS was developed. Instead of writing particles at fixed times, each particle is written according to its *own* timescale, e.g. every ten timesteps update. Crucially, particles are written independently in their own record, without requiring global synchronisation. Synchronisation is only required by the reader and is achieved through interpolation. Using interpolation, we can get positions, velocities, internal energies, ... at *any time* between the simulation start and end.

Implementation details can be found in the following paper, `Hausammann et al (2022) <https://ui.adsabs.harvard.edu/abs/2022A%26C....4100659H/abstract>`_.

Data writing
------------

The file where the CSDS writes is called the *logfile*. The logfile begins with a header that contains parameters for file handling and data identification. All available masks are listed by name and data size, followed by the specific masks designated for each particle type, ordered by their writing sequence.

Every time the CSDS logs data, it writes an individual *record*. Each record is composed of a small record header followed by the particle data. The record header stores two primary pieces of information:

* A *mask*: A bit mask that identifies which specific fields (coordinates, velocities, internal energies, etc.) are present in the record's data section. This allows a record to be identified as containing either a timestamp (see below) or specific particle fields.
* An *offset*: The distance in the file to the previous record belonging to the same particle. This offset is key to efficiently tracking a particle's history backwards through the file.

At the beginning of the simulation, the CSDS first records the time and initial conditions (positions, velocities, ...). During the simulation execution, the CSDS writes the current timestep at the start of each step and includes data for only the required active particles. Upon completion of the run, the CSDS records the particle data at the final time. The file is then concluded with a second write of the timestep, which acts as a sentinel to reliably indicate the end of the data file.

In the current implementation, the CSDS allows writing a few particle fields through masks. For instance, we write the positions, velocities and accelerations. Two particularly important masks are:

* ``Timestamp``: this masks logs the current timestep;
* ``SpecialFlags``: this masks logs when particles are created, deleted, changed their type (e.g. in star formation, a gas particle is converted into a star) and when particles enter or leave the MPI domain.

Data reading
------------

When retrieving particle data at a specific time, the CSDS-reader accesses the particles in the logfile and performs a time interpolation between two records.

The reader provides functionality to retrieve and track a defined set of particles based on their IDs across the simulation history.

To simplify data analysis, the CSDS reader provides provides a Python wrapper around its C++ code. You can enable the wrapper through the configuration option ``--with-python``. The CSDS Python wrapper requires the Boost library and this can be specified with the configuration option ``--with-boost``. An example of python file to read the data is provided in ``csds/examples/reader_example.py``.

.. note::
   The accuracy is governed by the frequency of recording for each particle. For hydrodynamic particles, it is advised to use a small frequency (< 10 steps) since these can undergo shocks that change their properties nonlinearly.

Running a simulation with CSDS
------------------------------

First, you need to download the reader within the Swift repos. This is done with the following command:

.. code::
   bash

   git submodule update --init --recursive

To run with the CSDS, you will need to use the configuration option ``--enable-csds`` and the run time argument ``--csds``.

Parameters
***********

The CSDS has the following parameters:

* The frequency to write a particle record: ``delta_step``. This is the same for all particles.
* The basename of the CSDS files: ``basename``.
* The initial size in gigabytes of the CSDS logfile: ``initial_buffer_size``.
* The factor used to automatically increase the logfile's buffer size when the current space is exhausted: ``buffer_scale``.

Currently, there is no way to predict the size of the CSDS output. Hence, you should carefully choose ``initial_buffer_size`` and ``buffer_scale``. A practical initial size estimation is often a factor (e.g. 10 times) the IC's size. This is to ensure the CSDS has enough space to write the particles for the initial timestep record.

.. code::
   yaml

   # Parameters governing the CSDS snapshot system
   CSDS:
     delta_step:           10     # Update the particle log every this many updates
     basename:             index  # Common part of the filenames
     initial_buffer_size:  1      # (Optional) Buffer size in GB
     buffer_scale:	   10     # (Optional) When buffer size is too small, update it with required memory times buffer_scale

Examples
--------

Currently, two examples in Swift use the CSDS.

We recommend starting with the ``examples/HydroTests/SedovBlast_3D`` with the CSDS. You only need to compile with CSDS and add ``--csds`` in the ``run.sh`` file. Once the simulation is completed, you can use the example ``csds/examples/reader_example.py``. This file is kept up to date with the most recent changes and includes a call to all the existing functions.
If you need some extra information, a doc string is provided for the class ``csds.Reader`` and all its methods.

The second example is a dark matter only cosmological simulation ``examples/SmallCosmoVolume/SmallCosmoVolume_DM``. The examples use Velociraptor to find halos and particle IDs and track the dark matter position evolution.

If you wish to obtain a snapshot from the CSDS, a script is available in ``csds/examples/create_snapshot.py``.


Current implementation state
----------------------------

The CSDS is currently (November 2025) only implemented for a few modules:

* hydro: ``SPHENIX`` and ``Gadget2``
* Stars: ``GEAR`` and ``Basic``
* Star formation: ``GEAR`` and ``none``
* Chemistry: ``GEAR``, ``AGORA`` and ``none``
* Gravity: ``MultiSoftening``

Extending the CSDS to other modules is easy. You only need to add the CSDS structure to the particles and implement the IO functions (see ``src/hydro/Gadget2/hydro_part.h``, ``src/hydro/Gadget2/hydro_csds.c`` and ``src/hydro/Gadget2/hydro_csds.h``).

Current limitations
-------------------

In the current implementation, the number of fields that can be written to the logfile is limited by the mask size (``CSDS_MASK_SIZE`` in ``csds/src/definition.hpp``). The default mask size is 2, hence we have a maximum of 16 fields. As we need a timestamp and the special flag masks, there are in total 14 different fields that can be written in the logfile. Fields can be grouped, e.g. in the ``SPHENIX`` model with the ``SPHENIXSecondaryFields``.

One may want to increase the default mask size to increase the number of fields. However, as this mask header is written for each chunk, this also increases the size of the logfile. Hence, we do not recommend increasing it. A better idea is to implement a particle-type mask instead of a field-based mask.

On the particle side, sink particles, black holes and neutrinos are not yet supported. On the subgrid side, cooling is not supported.

In the current implementation, we log all fields. However, we could imagine logging different fields with different frequencies.
