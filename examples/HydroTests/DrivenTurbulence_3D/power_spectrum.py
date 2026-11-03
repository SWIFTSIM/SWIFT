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

"""Velocity power spectrum of a SWIFT snapshot the way Bauer & Springel (2012)
measure it: nearest-particle sampling of the velocity on a Cartesian grid
twice finer than the particle resolution, FFT, spherical shells in fine
logarithmic k bins merged to >= 12 modes, E(k) = 4 pi k^2 <P(k)> (L/2pi)^3 so that int E(k) dk = <|v|^2> over the grid.
Writes k, E(k), k E(k) and the rms Mach number.
Usage: python3 -I power_spectrum.py snapshot.hdf5 out.txt [grid] [cs]
"""

import sys, numpy as np, h5py
from scipy.spatial import cKDTree

fn, out = sys.argv[1], sys.argv[2]
with h5py.File(fn) as f:
    L = float(np.ravel(f["Header"].attrs["BoxSize"])[0])
    t = float(np.ravel(f["Header"].attrs["Time"])[0])
    x = f["PartType0/Coordinates"][:].astype(np.float64)
    v = f["PartType0/Velocities"][:].astype(np.float64)
    m = f["PartType0/Masses"][:].astype(np.float64)
n = len(x)
N = int(round(n ** (1 / 3)))
G = int(sys.argv[3]) if len(sys.argv) > 3 else 2 * N
cs = float(sys.argv[4]) if len(sys.argv) > 4 else 1.0
mach = np.sqrt((m * (v * v).sum(1)).sum() / m.sum()) / cs
# nearest-particle sampling on the grid (periodic)
tree = cKDTree(np.mod(x, L), boxsize=L)
g = (np.arange(G) + 0.5) * L / G
gx, gy, gz = np.meshgrid(g, g, g, indexing="ij")
pts = np.stack([gx.ravel(), gy.ravel(), gz.ravel()], axis=1)
_, idx = tree.query(pts, workers=4)
vg = v[idx].reshape(G, G, G, 3)
del pts, idx
sigma2 = (vg * vg).sum() / G**3  # <|v|^2> on the grid
# FFT; power per mode normalised so that sum over modes = sigma2
P = np.zeros((G, G, G // 2 + 1))
for c in range(3):
    vh = np.fft.rfftn(vg[..., c]) / G**3
    P += np.abs(vh) ** 2
# account for the missing half of the rfft (modes with kz > 0 appear once)
w = np.full(G // 2 + 1, 2.0)
w[0] = 1.0
if G % 2 == 0:
    w[-1] = 1.0
P *= w[None, None, :]
kf = np.fft.fftfreq(G, d=L / G) * 2 * np.pi
kz = np.fft.rfftfreq(G, d=L / G) * 2 * np.pi
kk = np.sqrt(kf[:, None, None] ** 2 + kf[None, :, None] ** 2 + kz[None, None, :] ** 2)
kmin = 2 * np.pi / L
kmax = kk.max()
edges = kmin * 2 ** (np.arange(0, np.log2(kmax / kmin) + 0.05, 0.05))
edges = np.concatenate([[0.5 * kmin], edges])
which = np.digitize(kk.ravel(), edges) - 1
Pw = P.ravel()
mask = (which >= 0) & (which < len(edges) - 1) & (kk.ravel() > 0)
# per fine bin: total power, number of modes of the full k-space and mean k. An rfft entry with kz > 0 stands for
# two modes of the full k-space, so the count must carry the same weights as the power; counting the rfft entries
# instead makes the mean power per mode (and E) too high by a factor ~2.
Wm = np.broadcast_to(w[None, None, :], P.shape).ravel()[mask]
power = np.bincount(which[mask], weights=Pw[mask], minlength=len(edges) - 1)
count = np.bincount(which[mask], weights=Wm, minlength=len(edges) - 1)
kc = np.bincount(
    which[mask], weights=kk.ravel()[mask] * Wm, minlength=len(edges) - 1
) / np.maximum(count, 1)
# merge bins until each has at least 12 modes (first bins: 18, 14, 24, 24, 12, 30, ... modes)
rows = []
acc_p = 0
acc_k = 0
acc_n = 0
lo = edges[0]
for i in range(len(power)):
    acc_p += power[i]
    acc_k += kc[i] * count[i]
    acc_n += count[i]
    if acc_n >= 12:
        # mean power per mode times the continuum shell volume 4 pi k^2 and the mode density (L/2pi)^3
        # (Bauer & Springel eq. 13); this removes the sawtooth that the discrete mode count per bin gives
        k = acc_k / acc_n
        E = 4 * np.pi * k**2 * (acc_p / acc_n) * (L / (2 * np.pi)) ** 3
        rows.append((k, E, k * E, acc_n))
        acc_p = acc_k = acc_n = 0
        lo = edges[i + 1]
rows = np.array(rows)
np.savetxt(
    out,
    rows,
    header=f"t={t:g} N={N} grid={G} sigma2={sigma2:.6e} mach_rms_massweighted={mach:.5f}\nk  E(k)  k*E(k)  n_modes",
)
print(
    f"{fn}: t={t:g} Mach={mach:.4f} sigma2={sigma2:.4e} int E dk={np.trapezoid(rows[:, 1], rows[:, 0]):.4e}"
)
