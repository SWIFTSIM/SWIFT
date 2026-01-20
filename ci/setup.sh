#!/bin/bash

#  Common setup and functions for jobs.
#
#  Peter W. Draper 20-JAN-2026.

# When exiting in error report current configuration options.
function ONEXIT {
   if test "$?" != 0; then
      echo "Current configuration: $(grep "\./configure" config.log)"
   fi
}
trap ONEXIT EXIT

# Wrap calls to make so we can use -j <something> easily.
# Also suppress output from libtool. Those logs are too big for gitlab.
function do_make {
    make -j 4 V=0 $*
}

# Wrap ./configure calls so we can echo that line to the log.
function do_configure {
    echo "## CONFIGURE: $*"
    ./configure $*
}

# More chat from the scripts. Re-enable if desperate.
#set -e
#set -x
