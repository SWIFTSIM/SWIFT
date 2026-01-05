.. GEAR sub-grid model
   Darwin Roduit, 5h January 2025

.. _gear_dev_notes:

GEAR developer notes
~~~~~~~~~~~~~~~~~~~~

This page groups some developer notes about the GEAR subgrid models.

Feedback order of operations
----------------------------

In GEAR, feedback operations are executed in a different conceptual order compared to EAGLE, although the task ordering itself is unchanged. The difference lies only in how feedback-related operations are staged conceptually.

The feedback sequence in GEAR is as follows:

#. **Drift**

#. **Hydrodynamics and gravity**

#. **Kick2**

#. **Star formation / star formation sink**:
   A star particle is formed.

#. **Stellar density and preparation tasks (``prep1`` to ``prep4``)**:
   Feedback-related quantities are computed.
   If the star particle has just been formed, there is nothing to do at this stage.

#. **Feedback application**:
   Feedback quantities are distributed to neighbouring gas particles, but the gas
   properties are *not updated yet*.
   Gas particles are flagged for synchronisation so that they become active in the next timestep.
   If the star particle has just been formed, there is nothing to do at this stage.

#. **Timestep update**:
   Stellar evolution quantities for the *next timestep* are computed (within
   ``feedback_will_do_feedback()``).
   If the star particle has just been formed, there is nothing to do at this stage.

The current timestep then ends, and a new one begins with:

#. **Drift** (``feedback_update_part()``):
   Gas particle properties are updated at this point, *before* the hydrodynamics loop.

Note that the *task order is not modified* by this scheme. For example, drift operations still always occur before hydrodynamics and gravity.

During the feedback loop, the quantities to be injected into gas particles, namely *internal energy, momentum, and metal mass* are computed but not applied immediately. Instead, they are injected during the drift step via ``feedback_update_part()``. Applying these updates directly within the feedback loop would make the result dependent on the order of particle interactions, since updates to the internal energy depend on the particle mass, which itself is modified by feedback.

Finally, stellar evolution is computed one timestep ahead of feedback. This allows star particles that will not produce any feedback to be skipped efficiently.
