.. main.c
   Mladen Ivkovic



Main Function
------------------

located at ``swiftsim/examples/main.c``

+ Initialize MPI
+ read command line arguments and options
+ write output parameter file if necessary (i.e. if it was given.)
+ checking, announcing, and setting up stuff : 
    + is parameter file given
    + check whether you chose appropriate flags
    + set up some other minor stuff
    + make some announcements for the log/screen
    + read and broadcast parameter file
+




Questions
------------------
+ What is output_parameter_file for? (main.c/359)
