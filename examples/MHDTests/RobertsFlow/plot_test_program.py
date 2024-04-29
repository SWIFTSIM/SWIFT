import pandas as pd
import os
from glob import glob
import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit

results_directory_name = "test_results/"
tau_max = 100
take_last = int(30/ (5e-2))


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

# Old function for calculating growth rates
def find_growth_rate(the_time, B_field, nlast=take_last):
    l = len(B_field)
    B_field_cut = B_field[l - 1 - nlast : -1]
    time_cut = the_time[l - 1 - nlast : -1]
    # print(time_cut,B_field_cut)
    res = np.polyfit(time_cut, B_field_cut, 1)[0]
    return np.round(res, 3)

def growth_rate_fit(time,growth_rate_1,growth_rate_2,amprel):
    return (growth_rate_1+amprel*growth_rate_2*np.exp((growth_rate_2-growth_rate_1)*time))/(1+amprel*np.exp((growth_rate_2-growth_rate_1)*time))

# Not used by default. The function calculates the growth rates for non-oscillating solutions and acconts for presence of subleading modes that can induce systematic uncertainty in growth determination
def find_growth_rate_with_subleading_mode(timestamps,B,nlast=take_last):
    dBdt = np.diff(B)/np.diff(timestamps)
    local_growth_rate = dBdt[-nlast:]/B[-nlast:]
    timestamps_cut = timestamps[-nlast:]

    mean_growth_rate = np.mean(local_growth_rate)
    std_growth_rate = np.std(local_growth_rate)
      
    fit_res = curve_fit(f=growth_rate_fit,xdata=timestamps_cut,ydata=local_growth_rate,p0=[mean_growth_rate,0.5*mean_growth_rate,std_growth_rate/mean_growth_rate])
    parameters = fit_res[0]
    covmatrix = fit_res[1]
    return mean_growth_rate, std_growth_rate, parameters, covmatrix

def plot_info(run_data, the_key):
    fig, ax = plt.subplots(1, 3, figsize=(15, 5))
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
    for i in range(len(run_data)):
        run_data_slice = run_data.iloc[[i]]
        # print(run_data_slice)
        v0 = run_data_slice["v0"].values[0]
        Lbox = run_data_slice["Lbox"].values[0]
        kv = run_data_slice["kv"].values[0]
        kv0 = 2 * np.pi * kv / Lbox
        Rm = run_data_slice["Rm"].values[0]
        t_c = 1 / (v0 * kv0)
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
        Time = Time / t_c
        B = np.array(the_statistics[38])
        B = B / B[0]
        Ekin = np.array(the_statistics[13])
        Mh = np.abs(np.array(the_statistics[37]))
        divB = np.abs(np.array(the_statistics[35]))

        mask = Time <= tau_max
        Time = Time[mask]
        B = B[mask]
        Mh = Mh[mask]
        divB = divB[mask]

        dBdt = np.diff(B)/np.diff(Time)
        local_growth_rate = dBdt[-take_last:]/B[-take_last:]
        cut_timestamps = Time[-take_last:]

#       growth_rate = find_growth_rate(Time, np.log(B))
#       growth_rate_err = 0

#       mean_growth_rate,mean_growth_rate_err,fit_parameters, fit_parameter_error = find_growth_rate_with_subleading_mode(Time,B)
#       growth_rate = round(fit_parameters[0],6)
#       growth_rate_err = round(2*np.sqrt(fit_parameter_error[0,0]),6)

        growth_rate = round(np.mean(local_growth_rate),6)
        growth_rate_err = round(2*np.std(local_growth_rate),6)

        the_name = (
            "#"
            + str(run_data_slice["Run #"].values[0])
            + "_"
            + str(run_data_slice["Scheme"].values[0])
            + "_"
            + str(run_data_slice["IAfile"].values[0])
#            + "_$v_0$="
#            + str(v0)
            + "_$\eta$="
            + str(run_data_slice["eta"].values[0])
            + "_<s>="
            + str(growth_rate)
            + "$\pm$"
            + str(growth_rate_err)
        )
        ax[0].plot(Time, B)
        ax[1].plot(Time, divB, label=the_name)
        ax[2].plot(cut_timestamps, local_growth_rate)
#        ax[2].plot(cut_timestamps, growth_rate_fit(cut_timestamps,fit_parameters[0],fit_parameters[1],fit_parameters[2]))
        # ax[2].plot(Time, Mh)
    ax[0].set_xlabel("t/$t_c$", fontsize=8)
    ax[1].set_xlabel("t/$t_c$", fontsize=8)
    ax[2].set_xlabel("t/$t_c$", fontsize=8)
    ax[0].set_ylabel("<$B_{rms}$>/<$B_{rms}$(0)>", fontsize=8)
    ax[1].set_ylabel("<divB*h/B>", fontsize=8)
    ax[2].set_ylabel("$d(Ln(B))/dt$", fontsize=8)
    ax[1].legend(loc="best", fontsize=8)
    ax[0].set_yscale("log")
    # ax[1].set_yscale("log")
    ax[1].set_yscale("log")
    # ax[0].set_ylim(1,1e4)
    # ax[0].set_xlim(0,10)
    # ax[1].set_ylim(1e-4,1)
    plt.savefig(results_directory_name + name_of_the_plot + ".png", dpi=100)


def sort_and_plot(run_data, the_key):
    selected_runs = select_runs_to_plot(run_data, the_key)
    print(selected_runs)
    plot_info(selected_runs, the_key)


sort_key1 = {"Scheme": "All", "IAfile": "All", "v0": "All", "eta": "All", "Rm": "All"}
sort_key2 = {"Scheme": "FDI", "IAfile": "All", "v0": "All", "eta": "All", "Rm": "All"}
sort_key3 = {"Scheme": "ODI", "IAfile": "All", "v0": "All", "eta": "All", "Rm": "All"}
sort_key4 = {"Scheme": "VP", "IAfile": "All", "v0": "All", "eta": "All", "Rm": "All"}

run_data = load_test_run_parameters()
# type selected runs to plot
# selected_runs = ['ec359',
#'ec384',
#'ec385',
# ]
# mask = run_data['Run #'].isin(selected_runs)
# run_data = run_data[mask]

# fill selected range of runs to plot
run_data = run_data[-1:]
print(run_data)
sort_and_plot(run_data, sort_key1)
sort_and_plot(run_data, sort_key2)
sort_and_plot(run_data, sort_key3)
sort_and_plot(run_data, sort_key4)
# sort_and_plot(run_data, sort_key5)
# sort_and_plot(run_data, sort_key6)
# sort_and_plot(run_data, sort_key7)
# sort_and_plot(run_data, sort_key8)
# sort_and_plot(run_data, sort_key9)
# sort_and_plot(run_data, sort_key10)
