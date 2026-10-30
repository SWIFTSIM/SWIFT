###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2026 Matthieu Schaller (schaller@strw.leidenuniv.nl)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
##############################################################################
"""Mode amplitude, phase speed and damping of the sound wave from all snapshots."""
import glob
import sys

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ic = sys.argv[1] if len(sys.argv) > 1 else "soundWave.hdf5"
with h5py.File(ic, "r") as f:
    A = float(f["/Header"].attrs["SoundWaveAmplitude"])
    n = int(f["/Header"].attrs["SoundWaveNWave"])
    standing = bool(f["/Header"].attrs["SoundWaveStanding"])
    boost = float(f["/Header"].attrs["SoundWaveBoost"])
    cs = float(f["/Header"].attrs["SoundWaveSoundSpeed"])
    L = float(np.atleast_1d(f["/Header"].attrs["BoxSize"])[0])
k = 2.0 * np.pi * n / L
omega = cs * k

t, Zv, Zr = [], [], []
for fn in sorted(glob.glob("soundWave_[0-9][0-9][0-9][0-9].hdf5")):
    with h5py.File(fn, "r") as f:
        x = f["/PartType0/Coordinates"][:, 0] - boost * float(f["/Header"].attrs["Time"])
        vx = f["/PartType0/Velocities"][:, 0] - boost
        m = f["/PartType0/Masses"][:]
        rho = f["/PartType0/Densities"][:]
        t.append(float(f["/Header"].attrs["Time"]))
    w = m / m.sum()
    e = np.exp(-1j * k * x)
    Zv.append(2.0 * np.sum(w * vx * e))
    Zr.append(2.0 * np.sum(w * (rho / rho.mean() - 1.0) * e))
t, Zv, Zr = np.array(t), np.array(Zv), np.array(Zr)

if standing:
    E = np.abs(Zv) ** 2 + cs**2 * np.abs(Zr) ** 2
    damping = -0.5 * np.polyfit(t, np.log(E / (A * cs) ** 2), 1)[0]
    a = np.real(Zr * np.exp(1j * np.angle(Zr[0])) / np.abs(Zr[0]))  # ~ cos(omega t)
    crossings = np.where(np.diff(np.sign(a)) != 0)[0]
    omega_meas = np.pi / np.mean(np.diff(t[crossings])) if len(crossings) > 2 else float("nan")
    print(f"standing wave: omega / (cs k) = {omega_meas / omega:.4f}, "
          f"damping rate {damping:.3e} = {damping * 2 * np.pi / omega:.3e} per period")
    plt.semilogy(t, E / (A * cs) ** 2, ".-")
    plt.xlabel("t")
    plt.ylabel(r"$E_{\rm mode} / (A c_s)^2$")
else:
    phase = np.unwrap(np.angle(Zv))
    omega_meas = -np.polyfit(t, phase, 1)[0]
    damping = -np.polyfit(t, np.log(np.abs(Zv) / (A * cs)), 1)[0]
    print(f"travelling wave: phase speed / cs = {omega_meas / omega:.4f}, "
          f"damping rate {damping:.3e} = {damping * 2 * np.pi / omega:.3e} per period")
    fig, ax = plt.subplots(1, 2, figsize=(8, 3.2))
    ax[0].plot(t, np.abs(Zv) / (A * cs), ".-")
    ax[0].set_xlabel("t")
    ax[0].set_ylabel(r"$|Z_v| / (A c_s)$")
    ax[1].plot(t, phase + omega * t, ".-")
    ax[1].set_xlabel("t")
    ax[1].set_ylabel(r"phase $-$ exact [rad]")
plt.tight_layout()
plt.savefig("soundWave.png", dpi=200)
