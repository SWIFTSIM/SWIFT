#!/bin/bash

d1=gizmo/png
d2=anarchy/png
d3=pu/png

for i in {0..1040}
do
    echo $i
    tmp=`printf "%04i" $i`
    convert $d1/rayleigh_taylor_$tmp.png $d2/rayleigh_taylor_$tmp.png $d3/rayleigh_taylor_$tmp.png +append movie/rayleigh_taylor_$tmp.png
done

ffmpeg -framerate 10 -i movie/rayleigh_taylor_%04d.png -c:v libx264 -r 30 -pix_fmt yuv420p rayleigh_taylor.mp4
