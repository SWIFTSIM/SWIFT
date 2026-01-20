#!/bin/bash -l

#  SWIFT Intel toolchain complete set of builds and tests.
#
#  Peter W. Draper 20-JAN-2026.

source ci/intel-modules.sh
source ci/setup.sh

#  Start from clean sources.
git clean -fdx

echo "XXX check only XXX"

#  And off we go.
./autogen.sh

#  Checks without parallel HDF5.
source ci/swift-checks.sh

#  Checks with parallel HDF5.
source ci/intel-modules-parallel.sh
source ci/swift-checks-parallel.sh

#  Checks that distribution still works.
source ci/swift-checks-distribution.sh

exit
