#!/bin/bash

echo "Fetching initial conditions for the zoom in example..."
./getIC.sh

../../swift --cooling --feedback --cosmology  --limiter --sync --self-gravity --hydro --stars --star-formation --threads=8 zoom_in.yml 2>&1 | tee output.log
