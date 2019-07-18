.. Random number generator
   Folkert Nobels, 11th of July 2019

Random number generator
=======================

Often subgrid models require random numbers, for this purpose 
SWIFT has a random number generator. The random number generator
of SWIFT is based on a combination of the standard `rand_r` and `erand48`
random number generators. Since for some applications in cosmology
we want to be able to sample random numbers with a probability lower than 
:math:`10^{-8}`, we could not simply use the 32-bit `rand_r` due to its cut-off
and spacing of around :math:`2^{-32} \approx 2 \times 10^{-10}`.
For the `erand48` algorithm with 48 bits the spacing and cutoff are 
significantly smaller and around :math:`2^{-48} \approx 3.5 \times 10^{-15}`,
so this is very suitable for our applications. 

Reproducible random numbers
~~~~~~~~~~~~~~~~~~~~~~~~~~~

In our simulations we want to be able to reproduce the exact same random 
numbers if we do exactly the same simulation twice. Because of this we 
want to seed the random number generator in a reproducible way. In order to do this
we seed with the particle ID of the current particle and the current 
integer simulation time. 

To produce several different types of random numbers we have an additional
argument called the type of random number which is basically the nth random
number for the specified seed, which is added to the particle ID, thus providing
a distinct state per random number type.

Implementation
~~~~~~~~~~~~~~

Our random number generator packs the particle ID (plus the random number type) and
the current simulation time as two 64-bit values, plus a constant 16-bit value,
into a 144-bit buffer. This buffer is treated as an array of 9 `uint16` values.

In a first pass we initialize the seed to 0 and run through the 9 `uint16` values,
xor-ing them with the seed and calling `rand_r` on the seed to advance it. Using
`rand_r` with the thus-generated seed, we generate a sequence of 9 16-bit values
and xor them with the original 144-bit buffer.

The 9 bytes of this modified buffer are then used for three passes of `erand48`,
xor-ing the state in the same way as before. `erand48` is then called one final
time with the last state, producing a random double-precision value with a
48-bit mantissa.

What to do if we break the random number generator?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The most likely case is that the RNG is not strong enough for our application,
in this case we could simply do multiple passes of both shuffling the state and
generating the final value from the state. This increases the computational cost but
will make the random number generator stronger. 

An other case is that we need probabilities that are lower than :math:`1 \times 10^{-17}`,
in this case we simply cannot use our random number generator and for example
need to generate two random numbers in order to probe these low probabilities. 

