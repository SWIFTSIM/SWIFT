import os

if (
    os.path.exists("efficiencies.png")
    and os.path.exists("modes.png")
    and os.path.exists("spec_ang_mom.png")
    and os.path.exists("spinup.png")
):
    # no need to rerun script
    exit()

import numpy as np
from scipy.optimize import fsolve


def Z1(x):
    return 1 + (1 - x ** 2) ** 0.3333 * (
        (1 + np.absolute(x)) ** 0.3333 + (1 - np.absolute(x)) ** 0.3333
    )


def Z2(x):
    return np.sqrt(3 * x ** 2 + Z1(x) ** 2)


def r_hor(x):
    return 1 + np.sqrt(1 - x ** 2)


def r_isco(x):
    return 3 + Z2(x) - np.sign(x) * np.sqrt((3 - Z1(x)) * (3 + Z1(x) + 2 * Z2(x)))


def eps_NT(x):
    return 1 - np.sqrt(1 - 2 / 3 * 1 / r_isco(x))


def eff_sd(m, a, limit):
    return limit - (
        1
        / 10
        * 1
        / m
        * (
            0.985 / ((4.627 - 4.445 * a) ** -0.5524 + 1.6 * 1 / m)
            + 0.015 / ((827.3 - 718.1 * a) ** -0.706 + 1.6 * 1 / m)
        )
        * (0.9663 - 0.9292 * a) ** -0.5693
    ) / (1 - np.sqrt(1 - 2 / (3 * r_isco(a))))


def find_root(a, limit):
    return fsolve(eff_sd, x0=1, args=(a, limit))[0]


def A(x):
    return (0.9663 - 0.9292 * x) ** -0.5693


def B(x):
    return (4.627 - 4.445 * x) ** -0.5524


def C(x):
    return (827.3 - 718.1 * x) ** -0.706


def m_dot_crit1(x, limit):
    C1 = C(x) / B(x)
    eps_1 = 16 * limit / (0.985 * A(x)) * eps_NT(x)
    N1 = 0.015 / 0.985
    res1 = (
        1.6
        / B(x)
        * 1
        / (2 * C1 * eps_1)
        * (
            np.sqrt(
                (C1 * (1 - eps_1) + N1 - eps_1) ** 2 + 4 * eps_1 * C1 * (N1 - eps_1 + 1)
            )
            + C1 * (1 - eps_1)
            + N1
            - eps_1
        )
    )
    return res1


def m_dot_crit2(x, limit):
    C1 = C(x) / B(x)
    eps_1 = 16 * limit / (0.985 * A(x)) * eps_NT(x)
    N1 = 0.015 / 0.985
    res1 = (
        1.6
        / B(x)
        * 1
        / (2 * C1 * eps_1)
        * (
            -1
            * np.sqrt(
                (C1 * (1 - eps_1) + N1 - eps_1) ** 2 + 4 * eps_1 * C1 * (N1 - eps_1 + 1)
            )
            + C1 * (eps_1 - 1)
            + N1
            - eps_1
        )
    )
    return res1


def beta(alfa):
    return 1 / (1 + 2 * alfa)


def gamma(alpha):
    return (8 - 3 * beta(alpha)) / (6 - 3 * beta(alpha))


def t1(alpha):
    return -0.2703 * gamma(alpha) + 1.3603


def t2(alpha):
    return -0.94 + 4.475 * (gamma(alpha) - 1.444) - 5.1402 * (gamma(alpha) - 1.444) ** 2


def t3(alpha):
    return -0.84 * np.log10(alpha) - 0.919 + 0.643 * np.exp(-0.209 / alpha)


def t4(x):
    return (0.6365 * r_isco(x) - 0.4828) * (1 + 11.9 * np.exp(-0.838 * r_isco(x) ** 4))


def t5(x):
    return 1.444 * np.exp(-1.01 * r_isco(x) ** 0.86) + 0.1


def T(x, alpha):
    return (
        0.31
        * (1 + (t4(x) / r_hor(x)) ** 0.9) ** (t2(alpha) + t3(alpha))
        * 1
        / (r_hor(x) - t5(x)) ** (t1(alpha))
    )


def eta(x, alpha):
    return 1 + gamma(alpha) / (gamma(alpha) - 1) * T(x, alpha)


def f1(x):
    return 0.0871 * r_isco(x) - 0.1082


def f2(alpha):
    return 0.5 - 7.798 * (gamma(alpha) - 1.333) ** 1.26


def f3(x):
    return 0.153 * (r_isco(x) - 0.6) ** 0.3 + 0.105


def f4(x, alpha):
    return (
        f3(x)
        * (0.9 * gamma(alpha) - 0.2996)
        * (1.202 - 0.08 * (np.log10(alpha) + 2.5) ** 2.6)
    )


def f5(alpha):
    return -1.8 * gamma(alpha) + 4.299 - 0.018 + 0.018 * (np.log10(alpha) + 2) ** 3.571


def f6(x, alpha):
    return (
        f4(x, alpha)
        * (((0.14 * np.log10(r_hor(x) ** f5(alpha)) + 0.23) / f4(x, alpha)) ** 10 + 1)
        ** 0.1
    )


def L_adv(x, alpha):
    return (
        f2(alpha)
        + (f1(x) + 10 ** f6(x, alpha)) * (1.15 - 0.03 * (np.log10(alpha) + 3) ** 2.37)
    ) / eta(x, alpha)


a = np.arange(-1, 1, 0.0001)
mdotcrit1 = m_dot_crit1(a, 0.5)
mdotcrit2 = m_dot_crit2(a, 0.5)
m_a_90 = [find_root(x, 0.9) for x in a]
m_a_75 = [find_root(x, 0.75) for x in a]
m_a_50 = [find_root(x, 0.5) for x in a]

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

fig = plt.figure(figsize=(8, 6))

plt.style.use("classic")
plt.fill_between(a, [0.0001 for x in a], [0.028 for x in a], color="blue", alpha=0.2)
plt.fill_between(a, [0.028 for x in a], mdotcrit1, color="red", alpha=0.2)
plt.fill_between(a, mdotcrit1, [375 for x in a], color="orange", alpha=0.2)
plt.ylabel("$\dot{m}$", fontsize=24, usetex=True)
plt.xlabel("$a$", fontsize=24, usetex=True)
plt.tick_params(axis="y", right=True, direction="in")
plt.yscale("log")
plt.axis([-1, 1, 0.001, 100])
plt.text(-0.22, 0.004, "Thick disc", fontsize=20)
plt.text(-0.2, 0.33, "Thin disc", fontsize=20)
plt.text(-0.2, 18, "Slim disc", fontsize=20)
plt.minorticks_on()
plt.tick_params(
    axis="x",
    direction="in",
    bottom=True,
    top=True,
    length=8,
    width=1.2,
    which="major",
    labelsize=16,
)
plt.tick_params(
    axis="y",
    direction="in",
    left=True,
    right=True,
    length=8,
    width=1.2,
    which="major",
    labelsize=16,
)
plt.tick_params(
    axis="x",
    direction="in",
    bottom=True,
    top=True,
    length=4,
    width=0.9,
    which="minor",
    labelsize=16,
)
plt.tick_params(
    axis="y",
    direction="in",
    left=True,
    right=True,
    length=4,
    width=0.9,
    which="minor",
    labelsize=16,
)

plt.savefig("modes.png", bbox_inches="tight")
plt.close()

phi = -20.2 * a ** 3 - 14.9 * a ** 2 + 34.0 * a + 52.6
horizon_ang_vel = a / (2 * (1 + np.sqrt(1 - a ** 2)))
jet_factor = (
    0.05
    / (4.0 * np.pi)
    * 1
    / 0.3
    * phi ** 2
    * horizon_ang_vel ** 2
    * (1.0 + 1.38 * horizon_ang_vel ** 2 - 9.2 * horizon_ang_vel ** 4)
)
Z1_j = np.array(
    [
        1 + (1 - x ** 2) ** 0.333 * ((1 + abs(x)) ** 0.333 + (1 - abs(x)) ** 0.333)
        for x in a
    ]
)
Z2_j = np.array(np.sqrt(3 * a ** 2 + Z1_j ** 2))
r_iso = 3 + Z2_j - np.sign(np.array(a)) * np.sqrt((3 - Z1_j) * (3 + Z1_j + 2 * Z2_j))
eps_TD = 1 - np.sqrt(1 - 2 / (3 * r_iso))
eps_ADAF1 = 0.144 * (6 / r_iso) * eps_TD * min(1, 0.028 / 0.0044)
eps_ADAF2 = 0.144 * (6 / r_iso) * eps_TD * min(1, 0.001 / 0.0044)
Jet_ADAF = jet_factor * 0.3
Jet_SD = 0.22 * jet_factor
Jet_TD1 = 10 ** -3 * 0.1 ** (-0.1) * 100 ** 0.2 * 10 ** (2 * 0.1) * jet_factor
Jet_TD2 = 10 ** -3 * 0.1 ** (-0.1) * 10 ** (-1 * 0.1) * jet_factor
eps_SD1 = (
    1
    / 10
    * 1
    / 1
    * (
        0.985 / ((4.627 - 4.445 * a) ** -0.5524 + 1.6 * 1 / 1)
        + 0.015 / ((827.3 - 718.1 * a) ** -0.706 + 1.6 * 1 / 1)
    )
    * (0.9663 - 0.9292 * a) ** -0.5693
)
eps_SD2 = (
    1
    / 10
    * 1
    / 50
    * (
        0.985 / ((4.627 - 4.445 * a) ** -0.5524 + 1.6 * 1 / 50)
        + 0.015 / ((827.3 - 718.1 * a) ** -0.706 + 1.6 * 1 / 50)
    )
    * (0.9663 - 0.9292 * a) ** -0.5693
)

mdot_bh_ADAF1 = (1 - Jet_ADAF / 4.447) * (1 - eps_ADAF1 - Jet_ADAF)
mdot_bh_ADAF2 = (1 - Jet_ADAF / 4.447) * (1 - eps_ADAF2 - Jet_ADAF)
mdot_bh_ADAF3 = (1 - Jet_ADAF / 11.56) * (1 - eps_ADAF1 - Jet_ADAF)
mdot_bh_ADAF4 = (1 - Jet_ADAF / 11.56) * (1 - eps_ADAF2 - Jet_ADAF)
mdot_bh_TD1 = (1 - Jet_TD1 / 4.447) * (1 - eps_TD - Jet_TD1)
mdot_bh_TD2 = (1 - Jet_TD2 / 4.447) * (1 - eps_TD - Jet_TD2)


fig = plt.figure(figsize=(18, 4))
fig.subplots_adjust(wspace=0, hspace=0, top=1, bottom=0)
gs = gridspec.GridSpec(1, 3, width_ratios=[1, 1, 1])
plt.style.use("classic")

plt.subplot(gs[0])
plt.plot(
    a,
    eps_ADAF2,
    linewidth=2,
    label="$\epsilon_\mathrm{rad}(\dot{m}<0.0044)$",
    color="red",
)
plt.plot(
    a,
    eps_ADAF1,
    linewidth=2,
    label="$\epsilon_\mathrm{rad}(\dot{m}=0.028)$",
    color="red",
    linestyle="--",
)
plt.plot(a, 0.97 * Jet_ADAF, linewidth=2, label="$\epsilon_\mathrm{jet}$", color="blue")
plt.fill_between(a, eps_ADAF1, eps_ADAF2, color="red", alpha=0.2)
plt.ylabel("$\epsilon_\mathrm{feedback}$", fontsize=24, usetex=True)
plt.xlabel("$a$", fontsize=24, usetex=True)
plt.tick_params(axis="y", right=True, direction="in")
plt.legend(loc="upper left", prop={"size": 13})
plt.xticks([-1, -0.5, 0, 0.5, 1], [-1, -0.5, 0, 0.5, 1])
plt.yticks(
    [0.0001, 0.001, 0.01, 0.1, 1, 10, 100],
    ["", 10 ** (-3), 10 ** (-2), 10 ** (-1), 10 ** (-0), 10 ** (1), 10 ** (2)],
)
plt.minorticks_on()
plt.axis([-1, 1, 0.0001, 10])
plt.yscale("log")
plt.tick_params(
    axis="x",
    direction="in",
    bottom=True,
    top=True,
    length=8,
    width=1.2,
    which="major",
    labelsize=16,
)
plt.tick_params(
    axis="y",
    direction="in",
    left=True,
    right=True,
    length=8,
    width=1.2,
    which="major",
    labelsize=16,
)
plt.tick_params(
    axis="x",
    direction="in",
    bottom=True,
    top=True,
    length=4,
    width=0.9,
    which="minor",
    labelsize=16,
)
plt.tick_params(
    axis="y",
    direction="in",
    left=True,
    right=True,
    length=4,
    width=0.9,
    which="minor",
    labelsize=16,
)
plt.title("Thick disc", fontsize=16)

plt.subplot(gs[1])
plt.plot(a, 0.97 * eps_TD, linewidth=2, label="$\epsilon_\mathrm{rad}$", color="red")
plt.plot(
    a,
    Jet_TD2,
    linewidth=2,
    label="$\epsilon_\mathrm{jet}(\dot{m}=0.028,M_\mathrm{BH}=10^9 \mathrm{M}_\odot)$",
    color="blue",
)
plt.plot(
    a,
    Jet_TD1,
    linewidth=2,
    label="$\epsilon_\mathrm{jet}(\dot{m}=1,M_\mathrm{BH}=10^6 \mathrm{M}_\odot)$",
    color="blue",
    linestyle="--",
)
plt.fill_between(a, Jet_TD1, Jet_TD2, color="blue", alpha=0.2)
plt.xlabel("$a$", fontsize=24, usetex=True)
plt.tick_params(axis="y", right=True, direction="in")
plt.yscale("log")
plt.legend(loc="upper left", prop={"size": 13})
plt.xticks([-1, -0.5, 0, 0.5, 1], ["", -0.5, 0, 0.5, 1])
plt.axis([-1, 1, 0.0001, 10])
plt.yticks([0.001, 0.01, 0.1, 1, 10], ["", "", "", "", ""])
plt.minorticks_on()
plt.tick_params(
    axis="x",
    direction="in",
    bottom=True,
    top=True,
    length=8,
    width=1.2,
    which="major",
    labelsize=16,
)
plt.tick_params(
    axis="y",
    direction="in",
    left=True,
    right=True,
    length=8,
    width=1.2,
    which="major",
    labelsize=16,
)
plt.tick_params(
    axis="x",
    direction="in",
    bottom=True,
    top=True,
    length=4,
    width=0.9,
    which="minor",
    labelsize=16,
)
plt.tick_params(
    axis="y",
    direction="in",
    left=True,
    right=True,
    length=4,
    width=0.9,
    which="minor",
    labelsize=16,
)
plt.title("Thin disc", fontsize=16)

plt.subplot(gs[2])
plt.plot(
    a, eps_SD1, linewidth=2, label="$\epsilon_\mathrm{rad}(\dot{m}=1)$", color="red"
)
plt.plot(
    a,
    eps_SD2,
    linewidth=2,
    label="$\epsilon_\mathrm{rad}(\dot{m}=50)$",
    color="red",
    linestyle="--",
)
plt.plot(a, Jet_SD, linewidth=2, label="$\epsilon_\mathrm{jet}$", color="blue")
plt.fill_between(a, eps_SD1, eps_SD2, color="red", alpha=0.2)
plt.xlabel("$a$", fontsize=24, usetex=True)
plt.tick_params(axis="y", right=True, direction="in")
plt.yscale("log")
plt.legend(loc="upper left", prop={"size": 13})
plt.xticks([-1, -0.5, 0, 0.5, 1], ["", -0.5, 0, 0.5, 1])
plt.axis([-1, 1, 0.0001, 10])
plt.yticks([0.001, 0.01, 0.1, 1, 10], ["", "", "", "", ""])
plt.minorticks_on()
plt.tick_params(
    axis="x",
    direction="in",
    bottom=True,
    top=True,
    length=8,
    width=1.2,
    which="major",
    labelsize=16,
)
plt.tick_params(
    axis="y",
    direction="in",
    left=True,
    right=True,
    length=8,
    width=1.2,
    which="major",
    labelsize=16,
)
plt.tick_params(
    axis="x",
    direction="in",
    bottom=True,
    top=True,
    length=4,
    width=0.9,
    which="minor",
    labelsize=16,
)
plt.tick_params(
    axis="y",
    direction="in",
    left=True,
    right=True,
    length=4,
    width=0.9,
    which="minor",
    labelsize=16,
)
plt.title("Slim disc", fontsize=16)

plt.savefig("efficiencies.png", bbox_inches="tight")

L_isco1 = [2 / 3 * 1 / np.sqrt(3) * (1 + 2 * np.sqrt(3 * r_isco(x) - 2)) for x in a]

plt.style.use("classic")
fig = plt.figure(figsize=(8, 6), linewidth=4)
plt.plot(a, L_isco1, linewidth=2, label="$\ell_\mathrm{ISCO}$", color="red")
plt.plot(
    a,
    0.45 * np.array(L_isco1),
    linewidth=3,
    linestyle="--",
    label=r"$0.45\ell_\mathrm{ISCO}$",
    color="red",
)
plt.plot(
    a,
    [L_adv(x, 0.1) for x in a],
    linewidth=2,
    label=r"$\ell_\mathrm{adv},\alpha=0.1$",
    color="green",
)
plt.plot(
    a,
    [L_adv(x, 0.2) for x in a],
    linewidth=2,
    label=r"$\ell_\mathrm{adv},\alpha=0.2$",
    color="teal",
)
plt.plot(
    a,
    [L_adv(x, 0.3) for x in a],
    linewidth=2,
    label=r"$\ell_\mathrm{adv},\alpha=0.3$",
    color="purple",
)
plt.ylabel("$\ell_\mathrm{in}$", fontsize=24, usetex=True)
plt.xlabel("$a$", fontsize=24, usetex=True)
plt.tick_params(axis="y", right=True, direction="in")
plt.legend(loc="upper right", prop={"size": 14})
plt.minorticks_on()
plt.axis([-1, 1, 0, 5])
plt.tick_params(
    axis="x",
    direction="in",
    bottom=True,
    top=True,
    length=8,
    width=1.2,
    which="major",
    labelsize=16,
)
plt.tick_params(
    axis="y",
    direction="in",
    left=True,
    right=True,
    length=8,
    width=1.2,
    which="major",
    labelsize=16,
)
plt.tick_params(
    axis="x",
    direction="in",
    bottom=True,
    top=True,
    length=4,
    width=0.9,
    which="minor",
    labelsize=16,
)
plt.tick_params(
    axis="y",
    direction="in",
    left=True,
    right=True,
    length=4,
    width=0.9,
    which="minor",
    labelsize=16,
)

plt.savefig("spec_ang_mom.png", bbox_inches="tight")
plt.close()

z1 = np.array(
    [
        1 + (1 - x ** 2) ** 0.333 * ((1 + abs(x)) ** 0.333 + (1 - abs(x)) ** 0.333)
        for x in a
    ]
)
z2 = np.array(np.sqrt(3 * a ** 2 + z1 ** 2))
r_iso = 3 + z2 - np.sign(np.array(a)) * np.sqrt((3 - z1) * (3 + z1 + 2 * z2))

phi = -20.2 * a ** 3 - 14.9 * a ** 2 + 34.0 * a + 52.6
horizon_ang_vel = a / (2 * (1 + np.sqrt(1 - a ** 2)))
jet_factor = (
    0.05
    / (4.0 * np.pi)
    * 1
    / 0.3
    * phi ** 2
    * horizon_ang_vel ** 2
    * (1.0 + 1.38 * horizon_ang_vel ** 2 - 9.2 * horizon_ang_vel ** 4)
)

da_TD_acc_only = 2 / 3 * 1 / np.sqrt(3) * (
    1 + 2 * np.sqrt(3 * r_iso - 2)
) - 2 * a * np.sqrt(1 - 2 / (3 * r_iso))
da_TD_Benson = (
    2 / 3 * 1 / np.sqrt(3) * (1 + 2 * np.sqrt(3 * r_iso - 2))
    - 2 * a * np.sqrt(1 - 2 / (3 * r_iso))
    - (1.25 * 10 ** -3 * 0.1 ** (-0.1) * 100 ** 0.2 * 10 ** (2 * 0.1) * jet_factor)
    * 2
    / a
    * (np.sqrt(1 - a ** 2))
    * (1 + np.sqrt(1 - a ** 2))
)
da_ADAF_acc_only = L_adv(a, 0.1) - 2 * a
da_ADAF_Benson = (
    L_adv(a, 0.1)
    - 2 * a
    - (jet_factor * 0.3) * 2 / a * (np.sqrt(1 - a ** 2)) * (1 + np.sqrt(1 - a ** 2))
)
eps_SD = (
    1
    / 10
    * 1
    / 10
    * (
        0.985 / ((4.627 - 4.445 * a) ** -0.5524 + 1.6 * 1 / 10)
        + 0.015 / ((827.3 - 718.1 * a) ** -0.706 + 1.6 * 1 / 10)
    )
    * (0.9663 - 0.9292 * a) ** -0.5693
)
da_SD_acc_only = L_adv(a, 0.1) - 2 * a * (1 - eps_SD)
da_SD_Benson = (
    L_adv(a, 0.1)
    - 2 * a * (1 - eps_SD)
    - (jet_factor * 0.22) * 2 / a * (np.sqrt(1 - a ** 2)) * (1 + np.sqrt(1 - a ** 2))
)


fig = plt.figure(figsize=(18, 4))
fig.subplots_adjust(wspace=0, hspace=0, top=1, bottom=0)
gs = gridspec.GridSpec(1, 3, width_ratios=[1, 1, 1])
plt.style.use("classic")

plt.subplot(gs[0])
plt.plot(a, da_TD_acc_only, linewidth=2, label="Accretion only", color="black")
plt.plot(a, da_TD_Benson, linewidth=1.5, label="Jet spindown included", color="blue")
plt.plot(a, [0 for x in a], linewidth=1.5, color="black", linestyle="--")
plt.ylabel("$\mathrm{d}a/\mathrm{d}\ln M_\mathrm{BH,0}$", fontsize=24, usetex=True)
plt.xlabel("$a$", fontsize=24, usetex=True)
plt.tick_params(axis="y", right=True, direction="in")
plt.legend(loc="lower left", prop={"size": 15})
plt.minorticks_on()
plt.axis([-1, 1, -4, 7])
plt.tick_params(
    axis="x",
    direction="in",
    bottom=True,
    top=True,
    length=8,
    width=1.2,
    which="major",
    labelsize=16,
)
plt.tick_params(
    axis="y",
    direction="in",
    left=True,
    right=True,
    length=8,
    width=1.2,
    which="major",
    labelsize=16,
)
plt.tick_params(
    axis="x",
    direction="in",
    bottom=True,
    top=True,
    length=4,
    width=0.9,
    which="minor",
    labelsize=16,
)
plt.tick_params(
    axis="y",
    direction="in",
    left=True,
    right=True,
    length=4,
    width=0.9,
    which="minor",
    labelsize=16,
)
plt.title("Thin disc", fontsize=16)

plt.subplot(gs[1])
plt.plot(a, da_ADAF_acc_only, linewidth=2, label="Accretion only", color="black")
plt.plot(a, da_ADAF_Benson, linewidth=1.5, label="Jet spindown included", color="blue")
plt.plot(a, [0 for x in a], linewidth=1.5, color="black", linestyle="--")
plt.xlabel("$a$", fontsize=24, usetex=True)
plt.tick_params(axis="y", right=True, direction="in")
plt.xticks([-1.0, -0.5, 0.0, 0.5, 1.0], ["", -0.5, 0.0, 0.5, 1.0])
plt.yticks([-8, -6, -4, -2, 0, 2, 4, 6, 8], ["", "", "", "", "", "", "", "", ""])
plt.minorticks_on()
plt.axis([-1, 1, -4, 7])
plt.tick_params(
    axis="x",
    direction="in",
    bottom=True,
    top=True,
    length=8,
    width=1.2,
    which="major",
    labelsize=16,
)
plt.tick_params(
    axis="y",
    direction="in",
    left=True,
    right=True,
    length=8,
    width=1.2,
    which="major",
    labelsize=16,
)
plt.tick_params(
    axis="x",
    direction="in",
    bottom=True,
    top=True,
    length=4,
    width=0.9,
    which="minor",
    labelsize=16,
)
plt.tick_params(
    axis="y",
    direction="in",
    left=True,
    right=True,
    length=4,
    width=0.9,
    which="minor",
    labelsize=16,
)
plt.title("Thick disc", fontsize=16)

plt.subplot(gs[2])
plt.plot(a, da_SD_acc_only, linewidth=2, label="Accretion only", color="black")
plt.plot(a, da_SD_Benson, linewidth=1.5, label="Jet spindown included", color="blue")
plt.plot(a, [0 for x in a], linewidth=1.5, color="black", linestyle="--")
plt.xlabel("$a$", fontsize=24, usetex=True)
plt.tick_params(axis="y", right=True, direction="in")
plt.xticks([-1.0, -0.5, 0.0, 0.5, 1.0], ["", -0.5, 0.0, 0.5, 1.0])
plt.yticks([-8, -6, -4, -2, 0, 2, 4, 6, 8], ["", "", "", "", "", "", "", "", ""])
plt.minorticks_on()
plt.axis([-1, 1, -4, 7])
plt.tick_params(
    axis="x",
    direction="in",
    bottom=True,
    top=True,
    length=8,
    width=1.2,
    which="major",
    labelsize=16,
)
plt.tick_params(
    axis="y",
    direction="in",
    left=True,
    right=True,
    length=8,
    width=1.2,
    which="major",
    labelsize=16,
)
plt.tick_params(
    axis="x",
    direction="in",
    bottom=True,
    top=True,
    length=4,
    width=0.9,
    which="minor",
    labelsize=16,
)
plt.tick_params(
    axis="y",
    direction="in",
    left=True,
    right=True,
    length=4,
    width=0.9,
    which="minor",
    labelsize=16,
)
plt.title("Slim disc", fontsize=16)

plt.savefig("spinup.png", bbox_inches="tight")
