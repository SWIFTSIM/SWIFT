#!/bin/bash

ffmpeg -framerate 10 -i output_%04d.png -c:v libx264 -r 30 -pix_fmt yuv420p stromgren_sphere-2D.mp4
