.. Time-step synchronization
   Matthieu Schaller 9th November 2019

.. _time_step_sync:
   
Time-step synchronization
=========================

Enabled with command-line option ``--sync``.

We also implement a synchronization step to change the time-step 
of particles that have been directly affected by external source 
terms, typically feedback events. Durier & Dalla Vecchia (2012) 
showed that the Saitoh & Makino (2009) mechanism was not sufficient 
in scenarios where particles receive energy in the middle of their 
regular time-step. When particles are affected by feedback (see 
Sections 8.1, 8.2, and 8.3 of Schaller et al. (2024), MNRAS 530:2), 
we flag them for synchronization. A final pass over the particles, 
implemented as a task acting on any cell which was drifted to the 
current time, takes these flagged particles, interrupts their current 
step to terminate it at the current time and forces them back onto 
the timeline (Section 2.4, ibid) at the current step. They then recompute 
their time-step and get integrated forward in time as if they were 
on a short time-step all along. This guarantees a correct propagation 
of energy and hence an efficient implementation of feedback. The use 
of this mechanism is always recommended in simulations with external source terms.
