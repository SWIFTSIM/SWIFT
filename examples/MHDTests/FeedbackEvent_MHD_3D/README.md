Feedback Event
==============

This is a modified version of the SedovBlast_3D test, that is supposed to be used
with the sub-grid model, i.e. with

+ Hydrodynamics
+ Cooling
+ Star formation

This is here to ensure that we are handling the feedback events correctly in the code.

We should not form any stars here, but ensure that flag is switched on just in case.

This test emulates what the EAGLE model does to particles for feedback, i.e.

+ Heats a single particle to 10^7.5 K
+ Does _not_ switch off cooling
+ Runs to completion.


Running Multiple Tests
----------------------

If you would like to run a suite of tests, try the runs.sh script. You'll
need to set the directories in the parameter file to be one higher, i.e.
../coolingtables rather than ./coolingtables.
