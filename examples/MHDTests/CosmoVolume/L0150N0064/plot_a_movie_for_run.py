#!/usr/bin/env python3

# Libraries
import subprocess
import numpy as np
import pandas as pd
from glob import glob
from PIL import Image
import os
import sys

path_to_snapshots = sys.argv[1]
path_to_images = sys.argv[2]
movie_name = sys.argv[3]
slice_at = sys.argv[4]
plane = sys.argv[5]


def make_mp4(frame_folder, name_mp4):
    this_directory = os.getcwd()
    go_to_frames_directory = " cd " + frame_folder
    command = "ffmpeg -framerate 1 -i %d.png " + name_mp4 + ".mp4"
    return_back = " cd " + this_directory
    command_sandwich = (
        " (" + go_to_frames_directory + " &&" + command + " &&" + return_back + " ) "
    )
    subprocess.call(command_sandwich, shell=True)


def plot_one_snapshot(snap_addr, img_addr):
    command = (
        "python "
        + "plot_snapshot_slice.py "
        + snap_addr
        + " "
        + img_addr
        + " "
        + slice_at
        + " "
        + plane
    )
    subprocess.call(command, shell=True)


def make_images_directory(path_to_run):
    this_directory = os.getcwd()
    go_to_run_directory = " cd " + path_to_run
    command = " mkdir frames "
    return_back = " cd " + this_directory
    command_sandwich = (
        " (" + go_to_run_directory + " &&" + command + " &&" + return_back + " ) "
    )

    print(command_sandwich)
    subprocess.call(command_sandwich, shell=True)


def make_movie_for_one_run():
    # make_images_directory(path_to_run)
    addr_book = sorted(glob(path_to_snapshots + "/snap_*.hdf5"))
    max_snapshot = int(len(addr_book))
    for i in range(min(len(addr_book), max_snapshot)):
        plot_one_snapshot(addr_book[i], path_to_images + "/" + str(i) + ".png")
    make_mp4(path_to_images, movie_name)


# make_movie_for_one_run('test_results/ec158')
# make_gif('test_results/ec157/frames', 'ec157',82)
# make_mp4('test_results/ec157/frames', 'ec157')
# data = data[176:182]
# print(data)
# make_movie_for_all_runs(data)
make_movie_for_one_run()
