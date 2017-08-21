###############################################################################
 # This file is part of SWIFT.
 # Copyright (c) 2016 Stefan Arridge (stefan.arridge@durham.ac.uk)
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
import h5py as h5
import matplotlib
matplotlib.use("Agg")
from pylab import *
import sys

def do_binning(x,y,x_bin_edges):

    #x and y are arrays, where y = f(x)
    #returns number of elements of x in each bin, and the total of the y elements corresponding to those x values

    n_bins = x_bin_edges.size - 1
    count = np.zeros(n_bins)
    y_totals = np.zeros(n_bins)
    
    for i in range(n_bins):
        ind = np.intersect1d(np.where(x > bin_edges[i])[0],np.where(x <= bin_edges[i+1])[0])
        count[i] = ind.size
        binned_y = y[ind]
        y_totals[i] = np.sum(binned_y)

    return(count,y_totals)


#for the plotting
max_r = float(sys.argv[1])
n_radial_bins = int(sys.argv[2])
n_snaps = int(sys.argv[3])

#some constants
OMEGA = 0.3 # Cosmological matter fraction at z = 0
PARSEC_IN_CGS = 3.0856776e18
KM_PER_SEC_IN_CGS = 1.0e5
CONST_G_CGS = 6.672e-8
CONST_m_H_CGS = 1.67e-24
h = 0.67777 # hubble parameter
gamma = 5./3.
eta = 1.2349
H_0_cgs = 100. * h * KM_PER_SEC_IN_CGS / (1.0e6 * PARSEC_IN_CGS)

#read some header/parameter information from the first snapshot

filename = "Hydrostatic_0000.hdf5"
f = h5.File(filename,'r')
params = f["Parameters"]
unit_mass_cgs = float(params.attrs["InternalUnitSystem:UnitMass_in_cgs"])
unit_length_cgs = float(params.attrs["InternalUnitSystem:UnitLength_in_cgs"])
unit_velocity_cgs = float(params.attrs["InternalUnitSystem:UnitVelocity_in_cgs"])
unit_time_cgs = unit_length_cgs / unit_velocity_cgs
v_c = float(params.attrs["IsothermalPotential:vrot"])
v_c_cgs = v_c * unit_velocity_cgs
header = f["Header"]
N = header.attrs["NumPart_Total"][0]
box_centre = np.array(header.attrs["BoxSize"])

#calculate r_vir and M_vir from v_c
r_vir_cgs = v_c_cgs / (10. * H_0_cgs * np.sqrt(OMEGA))
M_vir_cgs = r_vir_cgs * v_c_cgs**2 / CONST_G_CGS

for i in range(n_snaps):

    filename = "Hydrostatic_%04d.hdf5" %i
    f = h5.File(filename,'r')
    coords_dset = f["PartType0/Coordinates"]
    coords = np.array(coords_dset)

    #translate coords by centre of box
    header = f["Header"]
    snap_time = header.attrs["Time"]
    snap_time_cgs = snap_time * unit_time_cgs
    coords[:,0] -= box_centre[0]/2.
    coords[:,1] -= box_centre[1]/2.
    coords[:,2] -= box_centre[2]/2.
    radius = np.sqrt(coords[:,0]**2 + coords[:,1]**2 + coords[:,2]**2)
    radius_cgs = radius*unit_length_cgs
    radius_over_virial_radius = radius_cgs / r_vir_cgs

    #get the internal energies
    vel_dset = f["PartType0/Velocities"]
    vel = np.array(vel_dset)

    #make dimensionless
    vel /= v_c
    r = radius_over_virial_radius

    #find radial component of velocity
    v_r = np.zeros(r.size)
    for j in range(r.size):
        v_r[j] = -np.dot(coords[j,:],vel[j,:])/radius[j]

    bin_edges = np.linspace(0,max_r,n_radial_bins + 1)
    (hist,v_r_totals) = do_binning(r,v_r,bin_edges)
    
    bin_widths = bin_edges[1] - bin_edges[0]
    radial_bin_mids = np.linspace(bin_widths / 2. , max_r - bin_widths / 2. , n_radial_bins) 
    binned_v_r = v_r_totals / hist

    figure()
    plot(radial_bin_mids,binned_v_r,'ko',label = "Average radial velocity in shell")
    legend(loc = "upper right")
    xlabel(r"$r / r_{vir}$")
    ylabel(r"$v_r / v_c$")
    title(r"$\mathrm{Time}= %.3g \, s \, , \, %d \, \, \mathrm{particles} \,,\, v_c = %.1f \, \mathrm{km / s}$" %(snap_time_cgs,N,v_c))
    ylim((-1,1))
    plot_filename = "./plots/radial_velocity_profile/velocity_profile_%03d.png" %i
    savefig(plot_filename,format = "png")
    close()
