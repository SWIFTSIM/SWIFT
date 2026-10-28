import pandas as pd
import os
from glob import glob
import numpy as np
import matplotlib.pyplot as plt

results_directory_name = "test_results/"
tau_max = 1000.0  # 30 #70
take_last = int(900 / (5e-2))  # int(10 / (5e-2)) #int(40 / (5e-2))


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


def process_info(toplot=True, tosave=True):
    if toplot:
        fig, ax = plt.subplots(1, 2, figsize=(10, 5), constrained_layout=True)
    growth_rates = []
    growth_rate_errors = []
    divergence_errors = []
    divergence_errors_err = []
    the_addr = ("./statistics.txt")
    #the_addr = ("./test_results/Planar_graddiv_Brookshaw_eta0p0_Ded_01/statistics.txt")
    the_name = "the run"
    the_statistics = np.transpose(np.loadtxt(the_addr))
    Time = np.array(the_statistics[1])
    B = np.array(the_statistics[39])
    B = B / B[0]
    divB = np.abs(np.array(the_statistics[35]))

    mask = Time <= tau_max
    Time = Time[mask]
    B = B[mask]
    divB = divB[mask]
                                     

    dBdt = np.diff(B) / np.diff(Time)
    cut_B = B[1:].copy()
    local_growth_rate = dBdt[-take_last:] / cut_B[-take_last:]
    cut_timestamps = Time[-take_last:]
    if len(local_growth_rate)>len(cut_timestamps):
        local_growth_rate = local_growth_rate[:-1]
    elif len(local_growth_rate)<len(cut_timestamps):
        cut_timestamps = cut_timestamps[:-1]

    print(len(cut_timestamps),len(local_growth_rate))

    growth_rate = np.mean(local_growth_rate)
    growth_rate_err_2sigma = 2 * np.std(local_growth_rate)
    growth_rate_str = (
        str(np.round(growth_rate, 5))
        + "$\pm$"
        + str(np.round(growth_rate_err_2sigma, 5))
    )

    x = np.linspace(10,100,10)
    Bquad = x**2
    B4 = x**4
    Blin = x

    if toplot:
        ax[0].plot(Time, B)
        #ax[0].plot(x, Bquad/10, color="black", linestyle='dashed')
        #ax[0].plot(x, B4/100, color="black", linestyle='-.')
        #ax[0].plot(x, Blin, color="black", linestyle='dotted')
        ax[1].plot(cut_timestamps, local_growth_rate, label=the_name)
        #ax[0].set_xscale('log')

    growth_rates.append(growth_rate)
    growth_rate_errors.append(growth_rate_err_2sigma)
    divergence_errors.append(np.mean(divB[-take_last:]))
    divergence_errors_err.append(2 * np.std(divB[-take_last:]))

    if toplot:
        ax[0].set_xlabel("time", fontsize=18)
        ax[1].set_xlabel("time", fontsize=18)
        ax[0].set_ylabel(r"$B_{rms}/B_{rms}(0)$", fontsize=18)
        ax[1].set_ylabel("d(LnB)/dt", fontsize=18)
        ax[1].legend(loc="best", fontsize=18)
        ax[0].set_yscale("log")
        ax[0].tick_params(labelsize=18)
        ax[1].tick_params(labelsize=18)
        ax[0].grid()
        ax[1].grid()
        plt.savefig("statistics.png", dpi=100)
    if tosave:
        run_data["growth_rate"] = growth_rates
        run_data["growth_rate_err"] = growth_rate_errors
        run_data["divB_err"] = divergence_errors
        run_data["divB_err_err"] = divergence_errors_err
        run_data.to_csv(the_key["Scheme"] + "_results.csv", sep="\t")


process_info(toplot=True, tosave=False)
                                                    
