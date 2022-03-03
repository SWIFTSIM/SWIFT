#!/bin/env python

import h5py
import numpy as np
import healpy as hp

import read_lightcone as rl

plt.subplot(2,2,1)
totalmass_map = rl.read_map("./lightcones/", "lightcone0", shell_nr=0, map_name="TotalMass")
hp.mollview(totalmass_map+1, norm="log", title="Projected mass", hold=True)

plt.subplot(2,2,2)
gasmass_map = rl.read_map("./lightcones/", "lightcone0", shell_nr=0, map_name="SmoothedGasMass")
hp.mollview(gasmass_map+1, norm="log", title="Gas mass", hold=True)

plt.subplot(2,2,3)
stellarmass_map = rl.read_map("./lightcones/", "lightcone0", shell_nr=0, map_name="StellarMass")
hp.mollview(stellarmass_map+1, norm="log", title="Stellar mass", hold=True)

plt.subplot(2,2,4)
xray_map = rl.read_map("./lightcones/", "lightcone0", shell_nr=0, map_name="XrayROSATIntrinsicPhotons")
hp.mollview(xray_map+1e50, norm="log", title="ROSAT photons", hold=True)

plt.suptitle("SmallCosmoVolume lightcone")
