import pandas as pd
import os
from glob import glob
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update(
    {
        "font.family": "serif",
        "mathtext.fontset": "cm",  # Computer Modern
    }
)

results_directory_name = "test_results/"
Time_max = 70  # maximal time
Time_growth = 40  # growth interval measurement


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


def find_growth_rate(time, B_field):
    res, cov = np.polyfit(time, np.log(B_field), 1, cov=True)
    std_err = np.sqrt(np.diag(cov))
    return res, cov, std_err


def process_info(run_data, the_key, toplot=True, tosave=True):
    if toplot:
        fig, ax = plt.subplots(1, 2, figsize=(10, 5), constrained_layout=True)
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
    legend_font_size = min([14, int(14 * (5 / len(run_data)))])
    for i in range(len(run_data)):
        run_data_slice = run_data.iloc[[i]]
        # print(run_data_slice)
        Rm = run_data_slice["Rm"].values[0]
        the_addr = (
            results_directory_name
            + str(run_data_slice["Run #"].values[0])
            + "/statistics.txt"
        )
        the_statistics = np.transpose(np.loadtxt(the_addr))
        Time = np.array(the_statistics[1])
        B = np.array(the_statistics[39])
        B /= B[0]
        divB = np.abs(np.array(the_statistics[35]))

        mask_B = Time <= Time_max
        mask_growth = (Time <= min([Time_max, Time[-1]])) & (
            Time >= min([Time_max, Time[-1]]) - Time_growth
        )

        Time_cut = Time[mask_growth]
        dBdt = np.diff(B[mask_growth]) / np.diff(Time_cut)
        cut_B = 0.5 * (B[mask_growth][1:] + B[mask_growth][:-1])

        local_growth_rate = dBdt / cut_B
        T_mean_cut = 0.5 * (Time[mask_growth][1:] + Time[mask_growth][:-1])
        divB = divB[mask_growth]

        res, cov, std_err = find_growth_rate(Time[mask_growth], B[mask_growth])

        # one can choose growth rates as mean of the local growth rate or as a linear fit in log space
        growth_rate = res[0]  # np.mean(local_growth_rate)
        growth_rate_err_2sigma = std_err[0]  # 2 * np.std(local_growth_rate)
        growth_rate_str = (
            str(np.round(growth_rate, 5))
            + "$\pm$"
            + str(np.round(growth_rate_err_2sigma, 5))
        )

        the_name = str(run_data_slice["Run #"].values[0])

        if toplot:
            ax[0].plot(Time[mask_B], B[mask_B])
            ax[1].plot(T_mean_cut, local_growth_rate, label=the_name)

        growth_rates.append(growth_rate)
        growth_rate_errors.append(growth_rate_err_2sigma)
        divergence_errors.append(np.mean(divB))
        divergence_errors_err.append(2 * np.std(divB))

    from matplotlib.ticker import AutoMinorLocator

    for a in ax.flat:
        # Ticks on all four sides
        a.tick_params(
            axis="both",
            which="major",
            direction="in",
            top=True,
            right=True,
            length=7,
            width=1.2,
            labelsize=14,
        )

        a.tick_params(
            axis="both",
            which="minor",
            direction="in",
            top=True,
            right=True,
            length=4,
            width=1.0,
        )

        # Add minor ticks
        a.xaxis.set_minor_locator(AutoMinorLocator())
        a.yaxis.set_minor_locator(AutoMinorLocator())

    if toplot:
        ax[0].set_xlabel("time", fontsize=16)
        ax[1].set_xlabel("time", fontsize=16)
        ax[0].set_ylabel(r"$ B_{\rm rms}^{\rm norm}$", fontsize=16)
        ax[1].set_ylabel(r"$d(\log B_{\rm rms}^{\rm norm})/dt$", fontsize=16)
        ax[1].legend(loc="best", fontsize=legend_font_size)
        ax[0].set_yscale("log")
        ax[0].grid()
        ax[1].grid()
        # ax[0].set_aspect('equal')
        # ax[1].set_aspect('equal')
        plt.savefig(results_directory_name + name_of_the_plot + ".png", dpi=100)
    if tosave:
        run_data["growth_rate"] = growth_rates
        run_data["growth_rate_err"] = growth_rate_errors
        run_data["divB_err"] = divergence_errors
        run_data["divB_err_err"] = divergence_errors_err
        run_data.to_csv(the_key["Scheme"] + "_results.csv", sep="\t")


def sort_and_plot(run_data, the_key):
    selected_runs = select_runs_to_plot(run_data, the_key)
    print(selected_runs)
    process_info(selected_runs, the_key, toplot=True, tosave=True)


sort_key1 = {"Scheme": "All", "IAfile": "All", "v0": "All", "eta": "All", "Rm": "All"}

run_data = load_test_run_parameters()
# type selected runs to plot
# selected_runs = [
#'test1',
#'test2',
# ]
# mask = run_data['Run #'].isin(selected_runs)
# run_data = run_data[mask]

# fill selected range of runs to plot
run_data = run_data[:]  # 36:126
print(run_data)
sort_and_plot(run_data, sort_key1)
