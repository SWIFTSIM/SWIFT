#!/bin/bash

# Remove all output files from a previous run
rm -fv statistics.txt task_level*.txt *used_parameters.yml *.csv *.xmf
rm -fv timesteps.txt output*.txt
rm -fv *_0???.hdf5 snapshots/* restart/*
rm -fv *.png
