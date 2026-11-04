#!/bin/bash

snapshot="nifty_0303"
out_name="halo_0303"
vr_loc="./stf"
config_file="vrconfig_swift.cfg"

export OMP_NUM_THREADS=16

$vr_loc -i $snapshot -I 2 -o $out_name -C $config_file 
