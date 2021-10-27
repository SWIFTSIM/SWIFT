#!/bin/bash

if [ ! -f "eigenvals.txt" ]; then
    if [ ! -f "get_eigenvals.o" ]; then
        make
    fi

    ./get_eigenvalues.o
fi

python3 ./plot_eigenvalues.py
