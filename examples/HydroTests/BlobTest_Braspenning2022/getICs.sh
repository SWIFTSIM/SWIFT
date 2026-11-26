#!/bin/bash
set -eo pipefail

# Only download the files that are not already present
for f in blob_16_100.hdf5 blob_16_10.hdf5 blob_32_100.hdf5 blob_32_10.hdf5 blob_64_100.hdf5 blob_64_10.hdf5 blob_128_100.hdf5 blob_128_10.hdf5
do
    if [ ! -e $f ]
    then
        wget https://virgodb.cosma.dur.ac.uk/swift-webstorage/ICs/BlobTest_Braspenning22/$f
    fi
done
