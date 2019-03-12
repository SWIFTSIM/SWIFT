################################################################################
# This file is part of SWIFT.
# Copyright (c) 2018 Ashley Kelly ()
#                    Folkert Nobels (nobels@strw.leidenuniv.nl)
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
from galpy.potential import NFWPotential
from galpy.orbit import Orbit
import numpy as np
import matplotlib.pyplot as plt
from astropy import units
import h5py as h5

C = 8.0
M_200 = 2.0


def read_data():
    R = np.array([])
    z = np.array([])
    for frame in range(0, 599, 1):
        try:
            sim = h5.File("output_%04d.hdf5" % frame, "r")
        except IOError:
            break

        boxSize = sim["/Header"].attrs["BoxSize"][0]
        pos = sim["/PartType1/Coordinates"][:, :] - boxSize / 2.0
        R = np.append(R, np.sqrt(pos[0, 0] ** 2 + pos[0, 1] ** 2))
        z = np.append(z, pos[0, 2])
    return (R, z)


def galpy_nfw_orbit():
    # Setting up the potential
    nfw = NFWPotential(conc=C, mvir=M_200, H=70.0, wrtcrit=True, overdens=200)
    nfw.turn_physical_on()
    vxvv = [
        8.0 * units.kpc,
        0.0 * units.km / units.s,
        240.0 * units.km / units.s,
        0.0 * units.pc,
        5.0 * units.km / units.s,
    ]

    # Calculating the orbit
    ts = np.linspace(0.0, 0.58, 1000) * units.Gyr
    o = Orbit(vxvv=vxvv)
    o.integrate(ts, nfw, method="odeint")

    return o


o = galpy_nfw_orbit()
(R, z) = read_data()

o.plot()
plt.scatter(R, z, s=1, color="black", marker="x")
plt.savefig("comparison.png")
plt.close()
