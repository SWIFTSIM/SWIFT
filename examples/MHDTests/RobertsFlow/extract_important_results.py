#!/usr/bin/env python3

# Libraries
import subprocess
import glob
import pandas as pd

test_results_folder = './test_results/'

export_storage = 'data_for_export'

def get_important_file_addresses():
    statistics_files_addr = glob.glob(test_results_folder+'**/statistics.txt')
    df = pd.DataFrame({'Addr':statistics_files_addr})
    df['folder_name'] = [str(addr).split('/')[2] for addr in statistics_files_addr]
    return df

def create_results_directory(res_dirname):
    # construct command to create results directory
    mkdir_command = "mkdir " + res_dirname
    # show command
    print(mkdir_command)
    # execute
    subprocess.call(mkdir_command, shell=True)

def copy_results(results_data, res_dirname):
    path_to_new_dir = res_dirname + "/" + str(results_data['folder_name'].values[0])
    mkdir_command = " mkdir " + path_to_new_dir

    file_addr = str(results_data['Addr'].values[0])

    # construct commands to move each type of files
    cp_all_statistics_files_command = " cp " + file_addr+ " " + path_to_new_dir

    # construct command to copy everything
    command_sandwich = (
        " ("
        + mkdir_command
        + " &&"
        + cp_all_statistics_files_command
        + " ) "
    )
    # show command
    print(command_sandwich)
    # execute
    subprocess.call(command_sandwich, shell=True)


create_results_directory(export_storage)

df = get_important_file_addresses()

# export statistics files
for i in range(len(df)):
    results_data=df.iloc[[i]]
    copy_results(results_data,export_storage)

# copy test run parameters file
copy_parameters_command = " cp " + "./test_run_parameters.csv " + "./"+ export_storage+"/"
subprocess.call(copy_parameters_command, shell=True)
