#!/bin/bash

set -eo pipefail
# Convergence study: one run per resolution in run_N<N>/, then the L1 errors against N
# compared with the Arepo curves of Pakmor et al. 2016 (figure 4).
# Usage: [SWIFT_ARGS="-P SPH:resolution_eta:1.5 -P SPH:CFL_condition:0.2"] ./run_convergence.sh [N ...]
#        (default 32 64 128 256; 256 needs about ten minutes on 4 threads; SWIFT_ARGS
#        are extra swift arguments, e.g. MAGMA's recommended eta and CFL)
set -e
E=$(cd "$(dirname "$0")" && pwd)
[ -e "$E/glassPlane_128.hdf5" ] && [ -e "$E/glassPlane_32.hdf5" ] || ( cd "$E" && ./getGlass.sh )
for N in ${@:-32 64 128 256}; do
    D=$E/run_N$N
    mkdir -p "$D"
    cd "$D"
    [ -e yeeVortex.hdf5 ] || ( ln -sf "$E"/glassPlane_*.hdf5 . && python3 "$E/makeIC.py" $N )
    [ -e yeeVortex_0004.hdf5 ] || "$E/../../../swift" --hydro --threads=4 $SWIFT_ARGS "$E/yeeVortex.yml" 2>&1 | tee output.log
done
cd "$E"
python3 plotConvergence.py run_N*/yeeVortex_0004.hdf5
