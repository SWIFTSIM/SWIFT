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
"""Kinetic-energy decay of the Taylor-Green vortex and the effective numerical
viscosity nu_eff defined by E_kin(t) = E_kin(0) exp(-16 pi^2 nu_eff t)."""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

data = np.loadtxt("statistics.txt", comments="#", ndmin=2)
t = data[:, 1]
E = data[:, 13]
ok = E > 0
slope = np.polyfit(t[ok], np.log(E[ok] / E[0]), 1)[0]
nu_eff = -slope / (16.0 * np.pi**2)
print(f"E_kin(t_end) / E_kin(0) = {E[-1] / E[0]:.4f}, nu_eff = {nu_eff:.3e}")

plt.figure(figsize=(5, 4))
plt.semilogy(t, E / E[0], "k-", label="SWIFT")
plt.semilogy(t, np.exp(-16.0 * np.pi**2 * nu_eff * t), "r--", label=rf"$\nu_{{\rm eff}}$ = {nu_eff:.2e}")
plt.xlabel("t")
plt.ylabel(r"$E_{\rm kin}(t) / E_{\rm kin}(0)$")
plt.legend()
plt.tight_layout()
plt.savefig("taylorGreen.png", dpi=150)
