#!/bin/bash

# Run SWIFT
../../swift --hydro --cooling --limiter --threads=4 feedback.yml 2>&1 | tee output.log

# Plot the solution
python plotSolution.py 5
