#!/usr/bin/env python3
################################################################################
# This file is part of SWIFT.
# Copyright (c) 2024 Yves Revaz (yves.revaz@epfl.ch)
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

from pNbody import *
from pNbody import ic
from astropy import units as u
from astropy import constants as c
from pNbody import thermodyn, ctes

#################################################################
# Units
#################################################################

UnitLength_in_cm = 3.085e21
UnitMass_in_g = 4.435693e44
UnitVelocity_in_cm_per_s = 97824708.2699

k = (
    c.k_B.to(u.g * u.cm ** 2 / u.s ** 2 / u.K)
    / UnitVelocity_in_cm_per_s ** 2
    / UnitMass_in_g
).value  # Boltzman constant   in code units
m_p = (c.m_p.to(u.g) / UnitMass_in_g).value  # proton mass         in code units

kmstoCodeUnits = (1 / (UnitVelocity_in_cm_per_s * u.cm / u.s).to(u.km / u.s)).value
MsoltoCodeUnits = 1 / (UnitMass_in_g * u.g).to(u.M_sun).value

#################################################################
# model parameters
#################################################################


m_ref = 200000.0  # mass of gas particles, in solar mass
nf = 1  # particle number multiplicative factor (used to reduce noise)
hydro = 1  # 1=gas is treated with SPH, 0=gas will be collision-less


# mass ratio between components
fm_gas = 1.0
fm_disk = 1.0
fm_bulge = 1.0
fm_halo = 5.0

# Stellar disk
if fm_disk == 0:
    M_disk = 0
else:
    M_disk = 2e10 * MsoltoCodeUnits  # Stellar disk total mass
Hr_disk = 2.0
Hz_disk = 0.3
fr_disk = 10.0
fz_disk = 10.0
toomre_disk = 2.0

# Dark halo
if fm_halo == 0:
    M_halo = 0
else:
    M_halo = 100e10 * MsoltoCodeUnits  # Dark halo total mass
Hr_halo = 50.0  # for plummer=30, for nfw, determined rmax
Rs_halo = 15.0  # for nfw
dR_halo = 0.25
fr_halo = 3.0

# Bulge
if fm_bulge == 0:
    M_bulge = 0
else:
    M_bulge = 0.4e10 * MsoltoCodeUnits  # Bulge total mass
Hr_bulge = 1.0
fr_bulge = 5.0

# Gas disk
if fm_gas == 0:
    M_gas = 0
else:
    M_gas = 0.5e10 * MsoltoCodeUnits  # Gas disk total mass
Hz_gas = 0.3
Hr_gas = 4.0
Rf = 40.0
Hr_gas = Hr_gas - Hz_gas
rmax_gas = 10 * Hr_gas
zmax_gas = 3 * Hz_gas
sigmavel_gas = 10 * kmstoCodeUnits  # km/s to code units
T_gas = 100  # K


# estimation of the gravitational softening length
eps = 0.05 / (m_ref / 200000) ** (1 / 3.0)

m_ref = m_ref * MsoltoCodeUnits  # to code units


#################################################################
# parameters for the velocities
#################################################################

ErrTolTheta = 0.5
AdaptativeSoftenning = False


###################################
# spherical components
###################################


# grid parameters halo
stats_name_halo = "stats_halo.dmp"
grmin_halo = 0  # grid minimal radius
grmax_halo = Hr_halo * fr_halo * 1.05  # grid maximal radius
nr_halo = 64  # number of radial bins
eps_halo = eps
# grid bins functions
rc_halo = Hr_halo / 5.0


def g_halo(r):
    return np.log(r / rc_halo + 1.0)


def gm_halo(r):
    return rc_halo * (np.exp(r) - 1.0)


# grid parameters bulge
stats_name_bulge = "stats_bulge.dmp"
grmin_bulge = 0  # grid minimal radius
grmax_bulge = Hr_bulge * fr_bulge * 1.05  # grid maximal radius
nr_bulge = 64  # number of radial bins
eps_bulge = eps
# grid bins functions
rc_bulge = Hr_halo


def g_bulge(r):
    return np.log(r / rc_bulge + 1.0)


def gm_bulge(r):
    return rc_bulge * (np.exp(r) - 1.0)


###################################
# cylindrical components
###################################

# grid parameters disk
stats_name_disk = "stats_disk.dmp"
grmin_disk = 0.0  # minimal grid radius
grmax_disk = Hr_disk * fr_disk  # maximal grid radius
gzmin_disk = -Hz_disk * fz_disk  # minimal grid z
gzmax_disk = Hz_disk * fz_disk  # maximal grid z
nr_disk = 32  # number of bins in r
nt_disk = 2  # number of bins in t
nz_disk = 64 + 1  # number of bins in z
# for an even value of nz, the potential is computed at z=0
# for an odd  value of nz, the density   is computed at z=0
eps_disk = eps
# grid bins functions
rc_disk = 3.0


def g_disk(r):
    return np.log(r / rc_disk + 1.0)


def gm_disk(r):
    return rc_disk * (np.exp(r) - 1.0)


mode_sigma_z = {"name": "jeans", "param": None}
mode_sigma_r = {"name": "isothropic", "param": 2}
mode_sigma_p = {"name": "epicyclic_approximation", "param": None}
params_disk = [mode_sigma_z, mode_sigma_r, mode_sigma_p]


# grid parameters gas
stats_name_gas = "stats_gas.dmp"
grmin_gas = 0.0  # minimal grid radius
grmax_gas = rmax_gas * 1.05  # maximal grid radius
gzmin_gas = -zmax_gas * 1.05  # minimal grid z
gzmax_gas = zmax_gas * 1.05  # maximal grid z
nr_gas = 32  # number of bins in r
nt_gas = 2  # number of bins in t
nz_gas = 64 + 1  # number of bins in z
# for an even value of nz, the potential is computed at z=0
# for an odd  value of nz, the density   is computed at z=0
eps_gas = eps
rc_gas = 3.0
# grid bins functions
def g_gas(r):
    return np.log(r / rc_gas + 1.0)


def gm_gas(r):
    return rc_gas * (np.exp(r) - 1.0)


mode_sigma_z = {"name": "jeans", "param": None}
mode_sigma_r = {"name": "constant", "param": sigmavel_gas}
mode_sigma_p = {"name": "epicyclic_approximation", "param": None}
params_gas = [mode_sigma_z, mode_sigma_r, mode_sigma_p]


#################################################################
# compute the number of particles for each component
#################################################################


# here we give explicitly the mass of the gas particles
m = m_ref

if fm_gas == 0:
    N_gas = 0
else:
    N_gas = int(M_gas / (m * fm_gas))

if fm_disk == 0:
    N_disk = 0
else:
    N_disk = int(M_disk / (m * fm_disk))

if fm_bulge == 0:
    N_bulge = 0
else:
    N_bulge = int(M_bulge / (m * fm_bulge))

if fm_halo == 0:
    N_halo = 0
else:
    N_halo = int(M_halo / (m * fm_halo))

print("N_gas   = %d" % N_gas)
print("N_disk  = %d" % N_disk)
print("N_bulge = %d" % N_bulge)
print("N_halo  = %d" % N_halo)
print("----------------------------")
print("N_tot   = %d" % (N_gas + N_disk + N_bulge + N_halo))
print("----------------------------")


print()
if N_gas > 0:
    print("m_gas   = %g Msol" % ((M_gas / N_gas) / MsoltoCodeUnits))

if N_disk > 0:
    print("m_disk  = %g Msol" % ((M_disk / N_disk) / MsoltoCodeUnits))

if N_bulge > 0:
    print("m_bulge = %g Msol" % ((M_bulge / N_bulge) / MsoltoCodeUnits))

if N_halo > 0:
    print("m_halo  = %g Msol" % ((M_halo / N_halo) / MsoltoCodeUnits))
print()


if nf > 1:
    N_gas = int(nf * N_gas)
    N_disk = int(nf * N_disk)
    N_bulge = int(nf * N_bulge)
    N_halo = int(nf * N_halo)


#################################################################
# generate models
#################################################################


#####################
# exponnential disk
#####################

nb_disk = None
if M_disk != 0.0:
    print("generating disk...")
    nb_disk = ic.expd(
        N_disk,
        Hr_disk,
        Hz_disk,
        fr_disk * Hr_disk,
        fz_disk * Hz_disk,
        irand=0,
        ftype="gh5",
    )
    nb_disk.set_tpe("disk")
    nb_disk.mass = (M_disk / N_disk) * np.ones(nb_disk.nbody).astype(np.float32)
    nb_disk.rename("disk.dat")
    nb_disk.write()


#####################
# halo
#####################


nb_halo = None
if M_halo != 0.0:
    print("generating halo...")
    nb_halo = ic.nfw(N_halo, Rs_halo, fr_halo * Hr_halo, dR_halo, ftype="gh5")
    nb_halo.set_tpe("halo")
    nb_halo.mass = (M_halo / N_halo) * np.ones(nb_halo.nbody).astype(np.float32)
    nb_halo.rename("halo.dat")
    nb_halo.write()


#####################
# bulge
#####################

nb_bulge = None
if M_bulge != 0.0:
    print("generating bulge...")
    nb_bulge = ic.plummer(
        N_bulge, 1, 1, 1, Hr_bulge, fr_bulge * Hr_bulge, vel="no", ftype="gh5"
    )
    nb_bulge.set_tpe("bulge")
    nb_bulge.mass = (M_bulge / N_bulge) * np.ones(nb_bulge.nbody).astype(np.float32)
    nb_bulge.rename("bulge.dat")
    nb_bulge.write()

#####################
# gas disk
#####################

nb_gas = None
if M_gas != 0.0:
    print("generating gas...")
    nb_gas = ic.miyamoto_nagai(
        N_gas, Hr_gas, Hz_gas, rmax_gas, zmax_gas, irand=-2, ftype="gh5"
    )
    nb_gas.set_tpe("gas")
    nb_gas.mass = (M_gas / N_gas) * np.ones(nb_gas.nbody).astype(np.float32)
    nb_gas.rename("gas.dat")
    nb_gas.write()


###############################################################
# merge all components
###############################################################

nb = None

if nb_disk is not None:
    if nb is None:
        nb = nb_disk
    else:
        nb = nb + nb_disk

if nb_halo is not None:
    if nb is None:
        nb = nb_halo
    else:
        nb = nb + nb_halo

if nb_bulge is not None:
    if nb is None:
        nb = nb_bulge
    else:
        nb = nb + nb_bulge

if nb_gas is not None:
    if nb is None:
        nb = nb_gas
    else:
        nb = nb + nb_gas

# save particles without velocities
nb.write("snapnf.hdf5")

###############################################################
# compute velocities
###############################################################


if nb_disk is not None:
    print("------------------------")
    print("disk velocities...")
    print("------------------------")

    nb_disk, phi, stats_disk = nb.Get_Velocities_From_Cylindrical_Grid(
        select="disk",
        disk=("disk", "gas"),
        eps=eps_disk,
        nR=nr_disk,
        nz=nz_disk,
        nt=nt_disk,
        Rmax=grmax_disk,
        zmin=gzmin_disk,
        zmax=gzmax_disk,
        params=params_disk,
        Phi=None,
        g=g_disk,
        gm=gm_disk,
        ErrTolTheta=ErrTolTheta,
        AdaptativeSoftenning=AdaptativeSoftenning,
    )
    iofunc.write_dmp(stats_name_disk, stats_disk)

    r = stats_disk["R"]
    z = stats_disk["z"]
    dr = r[1] - r[0]
    dz = z[nz_disk // 2 + 1] - z[nz_disk // 2]

    print("disk : Delta R :", dr, "=", dr // eps_disk, "eps")
    print("disk : Delta z :", dz, "=", dz // eps_disk, "eps")

    # reduce if needed
    if nf > 1:
        nb_disk = nb_disk.reduc(nf, mass=True)


if nb_gas is not None:
    print("------------------------")
    print("gas velocities...")
    print("------------------------")
    nb_gas, phi, stats_gas = nb.Get_Velocities_From_Cylindrical_Grid(
        select="gas",
        disk=("disk", "gas"),
        eps=eps_gas,
        nR=nr_gas,
        nz=nz_gas,
        nt=nt_gas,
        Rmax=grmax_gas,
        zmin=gzmin_gas,
        zmax=gzmax_gas,
        params=params_gas,
        Phi=None,
        g=g_gas,
        gm=gm_gas,
        ErrTolTheta=ErrTolTheta,
        AdaptativeSoftenning=AdaptativeSoftenning,
    )
    iofunc.write_dmp(stats_name_gas, stats_gas)

    r = stats_gas["R"]
    z = stats_gas["z"]
    dr = r[1] - r[0]
    dz = z[nz_gas // 2 + 1] - z[nz_gas // 2]
    print("gas   : Delta R :", dr, "=", dr / eps_gas, "eps")
    print("gas   : Delta z :", dz, "=", dz / eps_gas, "eps")

    # reduce if needed
    if nf > 1:
        nb_gas = nb_gas.reduc(nf, mass=True)


if nb_bulge is not None:
    print("------------------------")
    print("bulge velocities...")
    print("------------------------")
    nb_bulge, phi, stats_bulge = nb.Get_Velocities_From_Spherical_Grid(
        select="bulge",
        eps=eps_bulge,
        nr=nr_bulge,
        rmax=grmax_bulge,
        phi=None,
        g=g_bulge,
        gm=gm_bulge,
        UseTree=True,
        ErrTolTheta=ErrTolTheta,
    )
    iofunc.write_dmp(stats_name_bulge, stats_bulge)

    r = stats_bulge["r"]
    dr = r[1] - r[0]
    print("bulge : Delta r :", dr, "=", dr / eps_bulge, "eps")

    # reduce if needed
    if nf > 1:
        nb_bulge = nb_bulge.reduc(nf, mass=True)


if nb_halo is not None:
    print("------------------------")
    print("halo velocities...")
    print("------------------------")
    nb_halo, phi, stats_halo = nb.Get_Velocities_From_Spherical_Grid(
        select="halo",
        eps=eps_halo,
        nr=nr_halo,
        rmax=grmax_halo,
        phi=None,
        g=g_halo,
        gm=gm_halo,
        UseTree=True,
        ErrTolTheta=ErrTolTheta,
    )
    iofunc.write_dmp(stats_name_halo, stats_halo)

    r = stats_halo["r"]
    dr = r[1] - r[0]
    print("halo : Delta r :", dr, "=", dr / eps_halo, "eps")

    # reduce if needed
    if nf > 1:
        nb_halo = nb_halo.reduc(nf, mass=True)


###############################################################
# sum the different components and save the final model
###############################################################

nb = None

if nb_disk is not None:
    if nb is None:
        nb = nb_disk
    else:
        nb = nb + nb_disk

if nb_halo is not None:
    if nb is None:
        nb = nb_halo
    else:
        nb = nb + nb_halo

if nb_bulge is not None:
    if nb is None:
        nb = nb_bulge
    else:
        nb = nb + nb_bulge

if nb_gas is not None:
    if nb is None:
        nb = nb_gas
    else:
        nb = nb + nb_gas


# reorganize components
nb1 = nb.select("halo")
nb1.set_tpe("halo_1")  # position expected by swift

nb2 = nb.select("disk")
nb2.set_tpe("stars_1")  # position expected by swift

nb3 = nb.select("bulge")
nb3.set_tpe("stars_1")  # position expected by swift

nb4 = nb.select("gas")

if hydro == 0:
    nb4.set_tpe("bndry")

nb = nb1 + nb2 + nb3 + nb4

# convert to swift
nb = nb.set_ftype(ftype="swift")


# add units
np.UnitLength_in_cm = UnitLength_in_cm
nb.UnitMass_in_g = UnitMass_in_g
nb.UnitVelocity_in_cm_per_s = UnitVelocity_in_cm_per_s
nb.Unit_time_in_cgs = UnitLength_in_cm / UnitVelocity_in_cm_per_s

nb.boxsize = 150 * 10
nb.rsp_init = np.ones(nb.nbody) * eps
nb.birth_time_init = -1 * np.ones(nb.nbody)

if hydro:
    gamma = 5 / 3.0
    xi = 0.76
    ionisation = 0
    mu = thermodyn.MeanWeight(xi, ionisation)
    mumh = m_p * mu
    nb.u_init = T_gas / (gamma - 1.0) * k / mumh * np.ones(nb.nbody)

# save model
nb.rename("galaxy_multi_component.hdf5")
nb.write()

#%%
#Add the StellarParticleType attribute to the dataset
import h5py as h5
import numpy as np

nb_star = nb.select("stars")
N_star = np.sum(nb_star.npart)
star_tpe = 2 # Single population stars
star_type = np.ones(N_star)*star_tpe

with  h5.File("galaxy_multi_component.hdf5", "r+") as f:
    f["PartType4"].create_dataset("StellarParticleType", data=star_type)
