#!/bin/bash

#  Creates a graphic from the task graph file dependency_graph.dot.
#  Requires the graphviz command "dot".

if [ ! -e dependency_graph.dot ]; then
    echo "Missing task-graph output 'dependency_graph.dot'! Cannot generate figure."
else 
    dot -Tpng dependency_graph.dot -o task_graph.png
    echo "Output written to task_graph.png"
fi

exit
