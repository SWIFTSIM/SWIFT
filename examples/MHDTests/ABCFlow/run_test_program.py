#!/usr/bin/env python3

# Libraries
import subprocess
import numpy as np
import pandas as pd
from glob import glob
import os

# Example set up run parameters here
# parameters = {'Run #':['test_run'],'v0':[5.0],'Vz_factor':[1.0],'eta':[0.05],'kv':[1],'kb':[1],'Lbox':[1],'Flow_kind':[0], 'Scheme':['FDI'], 'IAfile':['g16'],'monopole_subtraction':[None],'artificial_diffusion':[None],'hyperbolic_dedner':[1.0],'hyperbolic_dedner_divv':[0.5],'parabolic_dedner':[1.0]}
# parameter_data = pd.DataFrame(data = parameters)
# parameter_data.to_csv('test_run_parameters.csv', sep = '\t',index=False)

# some parameters and dictionaries
###################################################################################################

# read parameters table
all_parameter_data = pd.read_csv("test_run_parameters.csv", sep="\t")
all_parameter_data = all_parameter_data.replace({np.nan: None})

# mask all runs to do
mask = all_parameter_data["Status"] != "done"

# Dictionaries for shortcuts
# full names of schemes
schemes_dict = {
    "ODI": "direct-induction",
    "FDI": "direct-induction-fede",
    "VP": "vector-potential",
}

# full names of initial arrangement files (generate the files before running!)
IA_dict = {
    "g16": "glassCube_16.hdf5",
    "g32": "glassCube_32.hdf5",
    "g64": "glassCube_64.hdf5",
    "s16": "stacked_16.hdf5",
    "s24": "stacked_24.hdf5",
    "s32": "stacked_32.hdf5",
    "s40": "stacked_40.hdf5",
    "s48": "stacked_48.hdf5",
    "s56": "stacked_56.hdf5",
    "s64": "stacked_64.hdf5",
    "s72": "stacked_72.hdf5",
    "s80": "stacked_80.hdf5",
    "s88": "stacked_88.hdf5",
    "s96": "stacked_96.hdf5",
    "s104": "stacked_104.hdf5",
    "s112": "stacked_112.hdf5",
    "s120": "stacked_120.hdf5",
    "s128": "stacked_128.hdf5",
    "c32": "c32.hdf5",
    "c64": "c64.hdf5",
}

# full names of forcing types
# forcing_dict = {"v": "roberts-flow", "a": "roberts-flow-acceleration"}

# Number of threads
threads = 18  # 60

# Where to store results
results_directory_name = "test_results"

# funcions
###################################################################################################

# Function that reads parameter dataframe
def read_parameter_csv(addr_to_csv):
    # read parameters table
    all_parameter_data = pd.read_csv(addr_to_csv, sep="\t")
    all_parameter_data = all_parameter_data.replace({np.nan: None})

    # mask all runs to do
    mask = all_parameter_data["Status"] != "done"
    return all_parameter_data


def update_parameter_csv(the_parameters, i):
    # set status to done
    the_parameters.at[i, "Status"] = "done"
    # write to parameter table
    the_parameters.to_csv("test_run_parameters.csv", sep="\t", index=False)


# Function for configuring simulation
def configure_simulation(scheme, forcing, spline, eos, path_to_lib=False):

    path_to_configure = " ../../../"
    scheme_opt = " --with-spmhd=" + scheme
    kernel_opt = " --with-kernel=" + spline
    forcing_opt = " --with-forcing=" + forcing
    eos_opt = " --with-equation-of-state=" + eos
    hand_vec = " --disable-hand-vec"
    compiler_warnings = " --disable-compiler-warnings"
    doxygen = " --disable-doxygen-doc"
    hydro_opt = ""  #' --with-hydro=minimal'

    # if specific path to some libraries is needed, set this up here
    if path_to_lib == True:
        path_to_libraries = " --with-fftw=/opt/homebrew/Cellar/fftw/3.3.10_1/ --with-gsl=/opt/homebrew/Cellar/gsl/2.7.1/ --with-hdf5=/opt/homebrew/Cellar/hdf5/1.14.2/bin/h5cc"
    else:
        path_to_libraries = ""

    # get full path of the run directory
    this_directory = os.getcwd()
    # go to the main swift directory
    goto_swift_directory_command = " cd " + path_to_configure
    # configure simulation with the options (scheme, forcing, other things...)
    configure_command = (
        " ./configure "
        + kernel_opt
        + scheme_opt
        + forcing_opt
        + eos_opt
        + hand_vec
        + compiler_warnings
        + doxygen
        + path_to_libraries
        + hydro_opt
    )
    # compile code (maybe consider make clean before that)
    make_command = " make -j"
    # go back to the run directory
    goto_this_directory_command = " cd " + this_directory
    # constuct command sandwich
    command_sandwich = (
        "( "
        + goto_swift_directory_command
        + " &&"
        + configure_command
        + " &&"
        + make_command
        + " &&"
        + goto_this_directory_command
        + " )"
    )
    # show sandwich
    print("Script sends command: " + command_sandwich)
    # execute
    try:
        subprocess.run(command_sandwich, shell=True, check=True)
        print("Configuring swift successful")
    except Exception as e:
        print("Error: " + str(e))
        raise Exception("Execution was stopped due to an error")
    # maybe shoud addd some try except here


# Funciton for making ICfile
def make_IC(phys_parameters, IAfile):
    if len(phys_parameters) != 0:
        # get configuration parameters from phys_parameters row of parameters dataframe
        v0 = " --rms_velocity=" + str(phys_parameters["v0"].values[0])
        Vz_factor = " --Vz_factor=" + str(phys_parameters["Vz_factor"].values[0])
        kv = " --velocity_wavevector=" + str(phys_parameters["kv"].values[0])
        kb = " --magnetic_wavevector=" + str(phys_parameters["kb"].values[0])
        Lbox = " --boxsize=" + str(phys_parameters["Lbox"].values[0])
        a = " --A=" + str(phys_parameters["A"].values[0])
        b = " --B=" + str(phys_parameters["B"].values[0])
        c = " --C=" + str(phys_parameters["C"].values[0])
        path = " --IA_path=" + "./IAfiles/" + IAfile

        # Construct command to make ICs with selected parameters
        command = (
            " python3 "
            + "make_IC.py"
            + v0
            + path
            + kv
            + kb
            + Lbox
            + Vz_factor
            + a
            + b
            + c
        )
    else:
        command = " python3 " + "make_IC.py"
    # show command
    print("Script sends command: " + command)
    # execute
    try:
        subprocess.run(command, shell=True, check=True)
        print("created ICs")
    except Exception as e:
        print("Error: " + str(e))
        raise Exception("Execution was stopped due to an error")


# Function for running simulation
def run_simulation(phys_parameters, threads):
    # get configuration parameters from phys_parameters row of parameters dataframe
    v0 = phys_parameters["v0"].values[0]
    Vz_factor = phys_parameters["Vz_factor"].values[0]
    eta = phys_parameters["eta"].values[0]
    kv = phys_parameters["kv"].values[0]
    Lbox = phys_parameters["Lbox"].values[0]
    # viscosity_alpha = phys_parameters["viscosity_alpha"].values[0]
    a = phys_parameters["A"].values[0]
    b = phys_parameters["B"].values[0]
    c = phys_parameters["C"].values[0]

    # get parameters for the scheme
    monopole_subtraction = phys_parameters["monopole_subtraction"].values[0]
    artificial_diffusion = phys_parameters["artificial_diffusion"].values[0]
    hyperbolic_dedner = phys_parameters["hyperbolic_dedner"].values[0]
    hyperbolic_dedner_divv = phys_parameters["hyperbolic_dedner_divv"].values[0]
    parabolic_dedner = phys_parameters["parabolic_dedner"].values[0]

    # get timestep parameters
    tau_max = phys_parameters["tau_max"].values[0]
    min_steps = phys_parameters["min_steps"].values[0]

    with_hydro = " --hydro"

    threads_opt = " --threads=" + str(threads)

    # Construct string command to set up endtime and max timestep options
    timeI_pref = " -P TimeIntegration:"
    set_time_end = ""
    set_dt_max = ""
    if tau_max != None:
        t_c = Lbox / (2 * np.pi * kv * v0)
        time_end = tau_max * t_c
        set_time_end = timeI_pref + "time_end:" + str(time_end)
        if min_steps != None:
            dt_max = time_end / min_steps
            set_dt_max = timeI_pref + "dt_max:" + str(dt_max)

    set_timeI_par = set_time_end + set_dt_max

    # Construct command to set up SPH parameters
    # set_av = ""
    # set_av_max = ""
    # set_av_min = ""
    # sph_pref = " -P SPH:"
    # if Forcing_kind == "a":
    #    set_av = sph_pref + "viscosity_alpha:" + str(viscosity_alpha)
    #    set_av_max = sph_pref + "viscosity_alpha_max:" + str(viscosity_alpha)
    #    set_av_min = sph_pref + "viscosity_alpha_min:" + str(viscosity_alpha)
    # set_sph_par = set_av + set_av_max + set_av_min

    # Construct string command to set up forcing options

    f_pref = " -P ABC_Flow_Forcing:"

    set_u0 = ""
    if v0 != None:
        set_u0 = f_pref + "u0:" + str(v0)

    set_kv = ""
    if kv != None:
        set_kv = f_pref + "kv:" + str(kv)

    set_Vz_factor = ""
    if Vz_factor != None:
        set_Vz_factor = f_pref + "Vz_factor:" + str(int(Vz_factor))
    set_a = ""
    if a != None:
        set_a = f_pref + "A:" + str(a)
    set_b = ""
    if b != None:
        set_b = f_pref + "B:" + str(b)
    set_c = ""
    if c != None:
        set_c = f_pref + "C:" + str(c)

    set_forcing_par = set_u0 + set_kv + set_Vz_factor + set_a + set_b + set_c

    # Construct string command to set up MHD parameters
    MHD_pref = " -P MHD:"
    set_eta = ""
    if eta != None:
        set_eta = MHD_pref + "resistive_eta:" + str(eta)

    set_monopole_subtraction = ""
    if monopole_subtraction != None:
        set_monopole_subtraction = (
            MHD_pref + "monopole_subtraction:" + str(monopole_subtraction)
        )

    set_artificial_diffusion = ""
    if artificial_diffusion != None:
        set_artificial_diffusion = (
            MHD_pref + "artificial_diffusion:" + str(artificial_diffusion)
        )

    set_hyperbolic_dedner = ""
    if hyperbolic_dedner != None:
        set_hyperbolic_dedner = MHD_pref + "hyperbolic_dedner:" + str(hyperbolic_dedner)

    set_hyperbolic_dedner_divv = ""
    if hyperbolic_dedner_divv != None:
        set_hyperbolic_dedner_divv = (
            MHD_pref + "hyperbolic_dedner_divv:" + str(hyperbolic_dedner_divv)
        )

    set_parabolic_dedner = ""
    if parabolic_dedner != None:
        set_parabolic_dedner = MHD_pref + "parabolic_dedner:" + str(parabolic_dedner)

    set_MHD_par = (
        set_eta
        + set_monopole_subtraction
        + set_artificial_diffusion
        + set_hyperbolic_dedner
        + set_hyperbolic_dedner_divv
        + set_parabolic_dedner
    )

    # Add all options
    set_all_par = set_timeI_par + set_forcing_par + set_MHD_par  # + set_sph_par

    # path to parameter file
    parameter_file = " ./ABCFlow.yml"

    # where to write output
    write_output_file = " | tee output_log.txt"
    # construct command to start simulation
    command = (
        " ../../../swift"
        + with_hydro
        + threads_opt
        + set_all_par
        + parameter_file
        + write_output_file
    )
    # show command
    print("Script sends command: " + command)
    # execute
    try:
        subprocess.run(command, shell=True, check=True)
        print("run complete")
    except Exception as e:
        print("Error: " + str(e))
        raise Exception("Execution was stopped due to an error")


# a function that creates the directory for all run data storage
def create_results_directory(res_dirname):
    # construct command to create results directory
    mkdir_command = "mkdir " + res_dirname
    # show command
    print("Creating the results directory if it is not there: " + mkdir_command)
    # execute
    subprocess.call(mkdir_command, shell=True)


def move_results(phys_parameters, res_dirname):
    run_nr = phys_parameters["Run #"].values[0]
    path_to_new_dir = res_dirname + "/" + str(run_nr)
    mkdir_command = " mkdir " + path_to_new_dir

    from_folder = "."

    # construct commands to move each type of files
    mv_all_txt_files_command = " mv " + from_folder + "/*.txt " + path_to_new_dir
    mv_snapshots_command = " mv " + from_folder + "/*.hdf5 " + path_to_new_dir
    mv_all_ymls_command = (
        " mv " + from_folder + "/u*.yml " + path_to_new_dir
    )  # all .yml but not RobertsFlow.yml
    mv_all_xmfs_command = " mv " + from_folder + "/*.xmf " + path_to_new_dir
    mv_all_csvs_command = (
        " mv " + from_folder + "/d*.csv " + path_to_new_dir
    )  # all .csv but not test_run_parameters

    # construct command to move everything
    command_sandwich = (
        " ("
        + mkdir_command
        + " &&"
        + mv_all_txt_files_command
        + " &&"
        + mv_snapshots_command
        + " &&"
        + mv_all_ymls_command
        + " &&"
        + mv_all_xmfs_command
        + " &&"
        + mv_all_csvs_command
        + " ) "
    )
    # show command
    print("Script sends command: " + command_sandwich)
    # execute
    subprocess.run(command_sandwich, shell=True, check=True)
    print("Moved files to results folder")


# main program that takes the_parameters from table, creates ICs, configures code and executes each row and stores data
def run_all():
    the_parameters = read_parameter_csv("test_run_parameters.csv")
    create_results_directory(results_directory_name)

    # make_IC([],'')

    # prepare_glass()
    for i in range(len(the_parameters)):
        if mask[i]:
            parameters_for_the_run = the_parameters.iloc[[i]]
            print(parameters_for_the_run)
            scheme = schemes_dict[parameters_for_the_run["Scheme"].values[0]]
            forcing = "abc-flow"
            configure_simulation(scheme, forcing, "quintic-spline", "isothermal-gas")

            IAfile = IA_dict[parameters_for_the_run["IAfile"].values[0]]
            make_IC(parameters_for_the_run, IAfile)

            run_simulation(parameters_for_the_run, threads)

            move_results(parameters_for_the_run, results_directory_name)

            update_parameter_csv(the_parameters, i)


# code execution
###################################################################################################
try:
    run_all()
except Exception as e:
    print("Error: " + str(e))
