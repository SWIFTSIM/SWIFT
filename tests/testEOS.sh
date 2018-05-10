#!/bin/bash

echo ""

rm -f testEOS_rho_u_P.txt

echo "Running  ./testEOS  21  1  0"

./testEOS  21  1  0

exit $?
