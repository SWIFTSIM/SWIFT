#!/bin/bash

# Generate the initial conditions if they are not present.
if [ ! -e multiTypes.hdf5 ]
then
    echo "Generating initial conditions for the multitype box example..."
    python makeIC.py 9 13 7 1
fi

../swift --hydro --external-gravity --stars --threads=1 multiTypes.yml 2>&1 | tee output.log
