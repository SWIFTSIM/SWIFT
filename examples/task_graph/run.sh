#!/bin/bash

function test {
    $1
    if [[ $? -ne 0 ]]; then
	echo "Swift failed, please verify if you have compiled with debugging checks (--enable-debugging-checks)"
	exit 1
    fi
}

 # Generate the initial conditions if they are not present.
if [ ! -e EAGLE_ICs_6.hdf5 ]
then
    echo "Fetching initial conditions for the EAGLE 6Mpc example..."
    ./getIC.sh
fi

test "../swift -x -s task_graph.yml"

dot -Tpng dependency_graph.dot -o task_graph.png

echo "Task dependency graph written in task_graph.png"
