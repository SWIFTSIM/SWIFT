.. Time-step limiter
   Matthieu Schaller 9th November 2019

.. _time_step_limiter:
   
Time-step limiter
=================

Enabled with command-line option ``--limiter``.

The first limit we impose is to limit the time-step of active particles
(section 3.4, Schaller et al. (2024), MNRAS 530:2). 
When a particle computes the size of its next time-step, typically using the CFL condition, 
it also additionally considers the time-step size of all the particles 
it interacted within the loop computing accelerations. We then demand 
that the particle of interest’s time-step size is not larger than a 
factor :math:`\Delta` of the minimum of all the neighbours’ values. We typically 
use :math:`\Delta = 4` which fits naturally within the binary structure of the 
time-steps in the code. This first mechanism is always activated in 
Swift and does not require any additional loops or tasks; it is, however, 
not sufficient to ensure energy conservation in all cases.

The time-step limiter proposed by Saitoh & Makino (2009) is
also implemented in SWIFT and is a recommended option for all
simulations not using a fixed time-step size for all particles. This
extends the simple mechanism described above, 
by also considering inactive particles and waking them up 
if one of their active neighbours uses a much smaller time-step size. 

This is implemented by means
of an additional loop over the neighbours at the end of the regular
sequence. Once an active particle has computed its time-step length for the next step, 
we perform an additional loop over its
neighbours and activate any particles whose time-step length differs
by more than a factor :math:`\Delta` (usually also set to 4). 

As shown by Saitoh & Makino (2009), this is necessary to conserve energy and hence
yield the correct solution even in purely hydrodynamics problems
such as a Sedov–Taylor blast wave. The additional loop over the
neighbours is implemented by duplicating the already existing tasks
and changing the content of the particle interactions to activate the
requested neighbours.
