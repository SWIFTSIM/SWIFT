#!/bin/bash

# Make a simple animation of the snapshots
out="earth_impact.mp4"
ffmpeg -framerate 5 -i earth_impact_%?%?%?%?.png $out -y

echo Saved $out
