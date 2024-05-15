#!/bin/bash

for snap in feedback_*.hdf5
do
        python3 talkPlot.py "$snap"
done

convert -delay 20 -loop 0 *.png Feedback.gif
