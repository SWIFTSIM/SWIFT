Blob Test
=========

This test is very similar to the one presented in Agertz+ 2007
(https://ui.adsabs.harvard.edu/abs/2007MNRAS.380..963A/abstract) that tests
un-seeded fluid mixing.

You can use the makeIC.py script to generate the initial conditions, and
then the makeMovie.py script to look at the results. Note that this is a very
dynamic test, and similarly to the kelvin-helmholtz it does _not_ have an explicit
answer and is best viewed as a movie with more than a grain of salt.

The expected results are:

+ For schemes with surface tension terms, the blob stays together
+ For schemes without surface tension terms, the blob breaks up.

This is at much lower resolution than the original test in Agertz+ 2007.
To change the resolution, you will need to edit the `makeIC.py` file.
A similar resolution to the GAD\_10M presented in that paper can be gained
by using `num_on_side=128`.
