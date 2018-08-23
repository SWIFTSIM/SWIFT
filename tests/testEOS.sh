#!/bin/bash

echo ""

rm -f testEOS_rho_u_P_*.txt

echo "Running testEOS for each planetary material"

A1_mat_id=(
    100
    101
    102
    200
    201
    202
    300
    301
    302
    303
)

for mat_id in "${A1_mat_id[@]}"
do
    ./testEOS "$mat_id" 1
done

exit $?
