#!/bin/bash

# Check that the required file is present
if [ ! -e dependency_graph.dot ]
then
    echo "Missing task-graph output! Can generate figure."
else 
    dot -Tpng dependency_graph.dot -o task_graph.png
fi

