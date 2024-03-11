###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2023 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
################################################################################

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import sys

plt.style.use("../../../tools/stylesheets/mnras.mplstyle")

data = np.loadtxt("statistics.txt")
t = data[:, 1]
E_kin = data[:, 13]
E_int = data[:, 14]
E_pot = data[:, 15]

E_tot = E_kin + E_int + E_pot

plt.figure(figsize=(5, 5))

plt.subplot(211)
plt.plot(t, E_tot, "k", label="total")
plt.plot(t, E_int, label="internal")
plt.plot(t, E_kin, label="kinetic")
plt.plot(t, E_pot, label="potential")

plt.legend(loc="upper left")

plt.ylabel("Energy")
plt.xlabel("Time")

plt.subplot(212)
plt.plot(t, E_tot / E_tot[0], "k")

plt.ylim(0.98, 1.02)
plt.ylabel("Total energy w.r.t. t=0")
plt.xlabel("Time")

plt.savefig("Energy.png")
