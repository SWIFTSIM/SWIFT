#!/bin/bash

output=cell_hierarchy.html

# merge all mpi ranks together
csv_output=cell_hierarchy.csv
if [ -f $csv_output ]
then
   rm $csv_output
fi

for filename in ./cell_hierarchy_*.csv;
do
    cat $filename >> cell_hierarchy.csv
done

# copy HTML page to the repository
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
cp $DIR/data/cell_hierarchy.html $output

echo $output has been generated
