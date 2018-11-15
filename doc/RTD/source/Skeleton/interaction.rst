.. How to print messages; error.h
   November 2018
   Mladen Ivkovic



Interacting with the User
----------------------------------


Some macros for the interaction with the user are available in ``swiftsim/src/error.h``

To get standardized output, you should use these defined macros:

+ ``error(msg)`` 
    print error message ``msg`` and abort.

+ ``message(msg)`` 
    print message ``msg`` to the log.
