.. Random number generator
   Folkert Nobels, 11th of July 2019

Random number generator
=======================

Often subgrid models require random numbers, for this purpose 
SWIFT has a random number generator. The random number generator
of SWIFT is based on Random123 from `Salmon et al. (2011) <https://dl.acm.org/citation.cfm?doid=2063405>`_ .
Random123 is a sophisticated random number generator that is a 
collection of counter-based random number generators. From Random123
SWIFT uses the Threefry algorithm to generate random numbers, Threefry
is a RNG that is based on the Treefish block cipher, in our case we 
use the 64 bit implementation because for some applications in cosmology
we want to be able to sample random numbers with a probability lower than 
:math:`10^{-8}`. If we want to do this properly and have enough statical
random numbers below this limit we are forced to use a 64 bit RNG, because
the 32 bit RNG has a cut-off and spacing of around 
:math:`2^{-32} \approx 2 \times 10^{-10}`.
For the Treefry algorithm with 64 bit the spacing and cut-off are 
significantly smaller and around :math:`2^{-64} \approx 5 \times 10^{-20}`,
so this is very suitable for our applications. 

Optimal parameters
~~~~~~~~~~~~~~~~~~

The Treefish algorithm has several parameters that can be tweaked. The main
two parameters that can be tweaked are the number of rounds that we want to 
use and if we want to use 2x64 or 4x64. In which 2x64 is in general faster 
for the same number of rounds but performs less for the same number of rounds.
After testing we found that for our practical purpose the fastest configuration
which still produces good enough random numbers is the 4x64 algorithm with 
12 cycles. Also according to `Salmon et al. (2011) <https://dl.acm.org/citation.cfm?doid=2063405>`_
this is still a number of cycles for which they are unaware of any statical 
flaws of Threefry. 

Reproducible random numbers
~~~~~~~~~~~~~~~~~~~~~~~~~~~

In our simulations we want to be able to reproduce the exact same random 
numbers if we do exactly the same simulation twice. Because of this we 
want to seed the Threefry algorithm in a reproducible way. In order to do this
we seed with the particle ID of the current particle and the current 
integer simulation time. Threefry is known to be able to generate random
numbers for very highly correlated seed values and this algorithm therefore
also passes all our tests in which we find no correlation between the 
random numbers. 

To produce several different types of random numbers we have an additional
argument called the type of random number which is basically the nth random
number for the specified seed, Threefry however does not require to calculate 
the random numbers before this and can immediately calculate the nth random
number without more calculation than the first.

What to do if we break the random number generator?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The most likely case is that the RNG is not strong enough for our application,
in this case we could simply increase the number of rounds to a higher value
for example 14 (maximum is 72). This will increase the computational cost but
will make the random number generator stronger. 

An other case is that we need probabilities that are lower than :math:`1 \times 10^{-17}`,
in this case we simply cannot use our random number generator and for example
need to generate two random numbers in order to probe these low probabilities. 

