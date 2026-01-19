#!/bin/bash -l

#  Tests to be ran when changes are committed to a merge request.

# When exiting in error report current configuration options.
function ONEXIT {
   if test "$?" != 0; then
      echo "Current configuration: $(grep "\./configure" config.log)"
   fi
}
trap ONEXIT EXIT

#  Lots of output.
set -e
set -x

#  Build toolchain.
source ci/intel-modules.sh

#  Clean sources.
git clean -fdx

#  And off we go.
./autogen.sh

echo
echo "non-MPI build"
echo "-------------"
./configure --disable-optimization
make
make clean

echo
echo "MPI build"
echo "---------"
./configure --with-parmetis --disable-optimization
make
make clean

#  Keep simple, may have a number of these happening.
exit
