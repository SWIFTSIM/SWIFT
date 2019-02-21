#!/usr/bin/env python
###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2018 Folkert Nobels (nobels@strw.leidenuniv.nl)
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
import numpy as np
from hmf import MassFunction
import hmf
from astropy.cosmology import FlatLambdaCDM



def getHMFz(z, H0=70.3, Om0=0.276, Ob0=0.0455, Tcmb0=2.725, Mmin=1e10, Mmax=1e15):
    """ Fast function to call the HMF from hmf, this function only has 
        7 variables and will return the dn/d(log10 M) and M array.
        z: redshift
        H0: Hubble constant
        Om0: Matter density
        Ob0: Baryon density
        Tcmb0: CMB temperature at z=0
        Mmin: minimum mass (solar masses)
        Mmax: Maximum mass (solar masses) 
    """
    new_model = FlatLambdaCDM(H0=H0, Om0=Om0, Ob0=Ob0, Tcmb0=Tcmb0)
    hmff = MassFunction(
        cosmo_model=new_model,
        Mmax=np.log10(Mmax),
        Mmin=np.log10(Mmin),
        z=z,
        hmf_model="ST",
    )
    return hmff.m, hmff.dndlog10m


def getHMFztinker(z, H0=70.3, Om0=0.276, Ob0=0.0455, Tcmb0=2.725, Mmin=1e10, Mmax=1e15):
    """ Fast function to call the HMF from hmf, this function only has 
        6 variables and will return the dn/d(log10 M) and M array.
        H0: Hubble constant
        Om0: Matter density
        Ob0: Baryon density
        Tcmb0: CMB temperature at z=0
        Mmin: minimum mass (solar masses)
        Mmax: Maximum mass (solar masses) 
    """
    new_model = FlatLambdaCDM(H0=H0, Om0=Om0, Ob0=Ob0, Tcmb0=Tcmb0)
    hmff = MassFunction(
        cosmo_model=new_model, Mmax=np.log10(Mmax), Mmin=np.log10(Mmin), z=z
    )
    return hmff.m, hmff.dndlog10m


