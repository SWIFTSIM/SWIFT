#!/bin/bash
set -eo pipefail

# Run SWIFT
../../../swift --cosmology --hydro --self-gravity --threads=28 santa_barbara.yml

