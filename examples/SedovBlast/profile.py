###############################################################################
 # This file is part of SWIFT.
 # Copyright (c) 2015 Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
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

import numpy as np
import h5py
import matplotlib
matplotlib.use('Agg')
import pylab as pl
import glob
import sedov

# plot the radial density profile for all snapshots in the current directory,
# as well as the analytical solution of Sedov

for f in sorted(glob.glob("output*.hdf5")):
    file = h5py.File(f, "r")
    t = file["/Header"].attrs["Time"]
    coords = np.array(file["/PartType0/Coordinates"])
    rho = np.array(file["/PartType0/Density"])

    radius = np.sqrt( (coords[:,0]-5.)**2 + (coords[:,1]-5.)**2 + \
                      (coords[:,2]-5.)**2 )
    
    r_theory, rho_theory = sedov.get_analytical_solution(100., 5./3., 3, t,
                                                         max(radius))
    
    pl.plot(r_theory, rho_theory, "r-")
    pl.plot(radius, rho, "k.")
    pl.savefig("{name}.png".format(name = f[:-5]))
    pl.close()
