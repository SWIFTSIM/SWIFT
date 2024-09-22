import pandas as pd
import os
from glob import glob
import numpy as np
import matplotlib.pyplot as plt

results_directory_name = "test_results/"
tau_max = 70 #30 #70
take_last = int(40 / (5e-2)) #int(10 / (5e-2)) #int(40 / (5e-2))


def load_test_run_parameters():
    run_data = pd.read_csv("test_run_parameters.csv", sep="\t")
    run_data["Rm"] = round(
        run_data["v0"]
        * run_data["Lbox"]
        / (2 * np.pi * run_data["kv"] * run_data["eta"]),
        2,
    )
    mask = run_data["Status"] == "done"
    run_data = run_data[mask]
    return run_data


def select_runs_to_plot(run_data, the_key):
    if the_key["Scheme"] != "All":
        mask = run_data["Scheme"] == the_key["Scheme"]
        run_data = run_data[mask]
    if the_key["IAfile"] != "All":
        mask = run_data["IAfile"] == the_key["IAfile"]
        run_data = run_data[mask]
    if the_key["v0"] != "All":
        mask = run_data["v0"] == float(the_key["v0"])
        run_data = run_data[mask]
    if the_key["eta"] != "All":
        mask = run_data["eta"] == float(the_key["eta"])
        run_data = run_data[mask]
    if the_key["Rm"] != "All":
        mask = run_data["Rm"] == float(the_key["Rm"])
        run_data = run_data[mask]
    return run_data


def find_growth_rate(the_time, B_field, nlast=take_last):
    l = len(B_field)
    B_field_cut = B_field[l - 1 - nlast : -1]
    time_cut = the_time[l - 1 - nlast : -1]
    # print(time_cut,B_field_cut)
    res = np.polyfit(time_cut, B_field_cut, 1)[0]
    return np.round(res, 3)


def process_info(run_data, the_key, toplot=True, tosave=True):
    if toplot:
        fig, ax = plt.subplots(1, 2,  figsize=(10, 5),constrained_layout=True)
        name_of_the_plot = (
            "Scheme="
            + the_key["Scheme"]
            + "_"
            + "IAfile="
            + the_key["IAfile"]
            + "_"
            + "v0="
            + the_key["v0"]
            + "_"
            + "eta="
            + the_key["eta"]
            + "_"
            + "Rm="
            + the_key["Rm"]
        )
    growth_rates = []
    growth_rate_errors = []
    divergence_errors = []
    divergence_errors_err = []
    for i in range(len(run_data)):
        run_data_slice = run_data.iloc[[i]]
        # print(run_data_slice)
        Rm = run_data_slice["Rm"].values[0]
        # print(run_data_slice)
        # print(t_c)
        # print(i,run_data_slice['Run #'])
        the_addr = (
            results_directory_name
            + str(run_data_slice["Run #"].values[0])
            + "/statistics.txt"
        )
        the_statistics = np.transpose(np.loadtxt(the_addr))
        Time = np.array(the_statistics[1])
        B = np.array(the_statistics[38])
        B = B / B[0]
        divB = np.abs(np.array(the_statistics[35]))

        mask = Time <= tau_max
        Time = Time[mask]
        B = B[mask]
        divB = divB[mask]


        dBdt = np.diff(B)/np.diff(Time)
        cut_B = B[1:].copy()
        local_growth_rate = dBdt[-take_last:]/cut_B[-take_last:]
        cut_timestamps = Time[-take_last:]

        growth_rate = np.mean(local_growth_rate)
        growth_rate_err_2sigma = 2*np.std(local_growth_rate)
        growth_rate_str = str(np.round(growth_rate,5))+"$\pm$"+str(np.round(growth_rate_err_2sigma,5))
        
        the_name = (
            "#"
            + str(run_data_slice["Run #"].values[0])
            + "_"
            + str(run_data_slice["Scheme"].values[0])
            + "_"
            + str(run_data_slice["IAfile"].values[0])
            + "_$\eta$="
            + str(run_data_slice["eta"].values[0])
            + "_<s>="
            + growth_rate_str
        )
        print(the_name)
        if toplot:
            ax[0].plot(Time, B)
            #ax[1].plot(Time, divB, label=the_name)
            ax[1].plot(cut_timestamps, local_growth_rate,label=the_name)

        growth_rates.append(growth_rate)
        growth_rate_errors.append(growth_rate_err_2sigma)
        divergence_errors.append(np.mean(divB[-take_last:]))
        divergence_errors_err.append(2*np.std(divB[-take_last:]))
        
    if toplot:
        ax[0].set_xlabel("time", fontsize=8)
        #ax[1].set_xlabel("t$", fontsize=8)
        ax[1].set_xlabel("time", fontsize=8)
        ax[0].set_ylabel("$B_{rms}$/$B_{rms}$(0)", fontsize=8)
        #ax[1].set_ylabel("<divB*h/B>", fontsize=8)
        ax[1].set_ylabel("d(LnB)/dt", fontsize=8)
        ax[1].legend(loc="best", fontsize=8)
        ax[0].set_yscale("log")
        #ax[1].set_yscale("log")
        #ax[0].set_ylim(1e-2,1e3)
        #ax[1].set_ylim(-0.15,0.15)
        ax[0].grid()
        ax[1].grid()
        # ax[0].set_xlim(0,10)
        # ax[1].set_ylim(1e-4,1)
        #ax[0].set_aspect('equal')
        #ax[1].set_aspect('equal')
        plt.savefig(results_directory_name + name_of_the_plot + ".png", dpi=100)
    if tosave:
        run_data['growth_rate'] = growth_rates
        run_data['growth_rate_err'] = growth_rate_errors
        run_data['divB_err'] = divergence_errors
        run_data['divB_err_err'] = divergence_errors_err
        run_data.to_csv(the_key["Scheme"]+'_results.csv',sep="\t")

def sort_and_plot(run_data, the_key):
    selected_runs = select_runs_to_plot(run_data, the_key)
    print(selected_runs)
    process_info(selected_runs, the_key,toplot=True,tosave=True)


sort_key1 = {"Scheme": "All", "IAfile": "All", "v0": "All", "eta": "All", "Rm": "All"}
#sort_key2 = {"Scheme": "FDI", "IAfile": "All", "v0": "All", "eta": "All", "Rm": "All"}
#sort_key3 = {"Scheme": "ODI", "IAfile": "All", "v0": "All", "eta": "All", "Rm": "All"}
#sort_key4 = {"Scheme": "VP", "IAfile": "All", "v0": "All", "eta": "All", "Rm": "All"}

run_data = load_test_run_parameters()
# type selected runs to plot
#selected_runs = [
#'16r1',
#'16r2',
#'16r3',
#'16r4',
#'16r1t2',
#'16r2t2',
#'16r3t2',
#'16r4t2',
#]
#mask = run_data['Run #'].isin(selected_runs)
#run_data = run_data[mask]

# fill selected range of runs to plot
run_data = run_data[:]#36:126]
print(run_data)
sort_and_plot(run_data, sort_key1)
#sort_and_plot(run_data, sort_key2)
#sort_and_plot(run_data, sort_key3)
#sort_and_plot(run_data, sort_key4)
# sort_and_plot(run_data, sort_key5)
# sort_and_plot(run_data, sort_key6)
# sort_and_plot(run_data, sort_key7)
# sort_and_plot(run_data, sort_key8)
# sort_and_plot(run_data, sort_key9)
# sort_and_plot(run_data, sort_key10)
