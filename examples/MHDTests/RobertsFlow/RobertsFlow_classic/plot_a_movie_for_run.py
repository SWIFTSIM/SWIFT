#!/usr/bin/env python3

# Libraries
import subprocess
import numpy as np
import pandas as pd
from glob import glob
from PIL import Image
import os

#Upload test_run_parameters
data=pd.read_csv('test_run_parameters.csv', sep='\t')

#Mask complete runs
mask = data['Status']=='done'
data=data[mask]
test_results_folder = 'test_results/'

def sort_key(x):
    y = int(os.path.basename(x).split('.')[0])
    return y

def make_gif(frame_folder, name_gif, duration):
    addr_book = sorted(glob(frame_folder+'/*.png'),key=sort_key)
    frames = [Image.open(image) for image in addr_book]
    print(addr_book)
    frame_one = frames[0]
    frame_one.save(frame_folder+'/'+name_gif+'.gif', format="GIF", append_images=frames, save_all=True, loop=0)

def make_mp4(frame_folder, name_mp4):
    this_directory = os.getcwd()
    go_to_frames_directory = ' cd '+frame_folder
    command = 'ffmpeg -framerate 1 -i %d.png '+name_mp4+'.mp4'
    return_back = ' cd '+this_directory
    command_sandwich = ' ('+go_to_frames_directory+' &&'+command+' &&'+return_back+' ) '
    subprocess.call(command_sandwich, shell=True)

def plot_one_snapshot(snap_addr,img_addr):
    command = 'python '+'plot_snapshot_slice.py '+snap_addr+' '+img_addr
    subprocess.call(command, shell=True)

def make_images_directory(path_to_run):
    this_directory = os.getcwd()
    go_to_run_directory = ' cd '+path_to_run
    command = ' mkdir frames '
    return_back = ' cd '+this_directory
    command_sandwich = ' ('+go_to_run_directory+' &&'+command+' &&'+return_back+' ) '

    print(command_sandwich)
    subprocess.call(command_sandwich, shell=True)


def make_movie_for_one_run(path_to_run):
    make_images_directory(path_to_run)
    addr_book = sorted(glob(path_to_run+'/RobertsFlow_*.hdf5'))
    name_of_run = os.path.basename(path_to_run)
    path_to_frames = path_to_run + '/frames'
    max_snapshot=int(len(addr_book)/2)
    for i in range(min(len(addr_book),max_snapshot)):
        plot_one_snapshot(addr_book[i], path_to_frames+'/'+str(i)+'.png')
    make_mp4(path_to_frames, name_of_run)


def make_movie_for_all_runs(the_data):
    for i in range(len(the_data)):
        run_name = the_data['Run #'].iloc[i]
        path_to_snapshots = test_results_folder+run_name
        print(path_to_snapshots)
        make_movie_for_one_run(path_to_snapshots)


#make_movie_for_one_run('test_results/ec158')
#make_gif('test_results/ec157/frames', 'ec157',82)        
#make_mp4('test_results/ec157/frames', 'ec157')  
data = data[176:182]
print(data)
make_movie_for_all_runs(data)

