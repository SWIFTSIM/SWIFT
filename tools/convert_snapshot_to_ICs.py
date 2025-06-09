#!/usr/bin/env python3

"""
Script to convert the NIFTY ICs to those that are compatible with SWIFT.
Note that this leaves h-factors as-is to be fixed in-place by SWIFT.

You will need:

    + swiftsimio
"""

from swiftsimio import Writer, load
import numpy as np
import unyt
import os

# Which file to read
#  filename = "./eagle_0000.hdf5"
filename = "output_0001.hdf5"
# What to call the output IC file
output_filename = "ICs.hdf5"


# -----------------------------------------------------------------------
# -----------------------------------------------------------------------


if not os.path.exists(filename):
    raise FileNotFoundError(filename)

warning = """
SWIFT snapshots are not written in a way that should be used as initial
conditions. By using this script, you're trying to do exactly that. Are
you sure you know what you're doing? [y/n]
"""

answer = input(warning)
if not answer.startswith("y") or answer.startswith("Y"):
    print("Quitting.")
    quit()

# Load data from snapshot.

snap = load(filename)

# Get units and create unyt.unitsystem.
length = snap.units.length
mass = snap.units.mass
time = snap.units.time
temperature = snap.units.temperature

my_units = unyt.UnitSystem("copiedFromSnapshot", length, mass, time, temperature)


# Store additional info about the snapshot in the ICs. Can't hurt to have the data lying around.
extra_header = {}

if snap.metadata.cosmology is not None:

    #  Accessing data through metadata.cosmology_raw works on older snapshots as
    #  well, so let's do it this way.

    extra_header[
        "snapshot cosmology Critical density [internal_units]"
    ] = snap.metadata.cosmology_raw["Critical density [internal units]"]
    extra_header["snapshot cosmology H [internal units]"] = snap.metadata.cosmology_raw[
        "H [internal units]"
    ]
    extra_header[
        "snapshot cosmology H0 [internal units]"
    ] = snap.metadata.cosmology_raw["H0 [internal units]"]
    extra_header[
        "snapshot cosmology Hubble time [internal units]"
    ] = snap.metadata.cosmology_raw["Hubble time [internal units]"]
    extra_header[
        "snapshot cosmology Lookback time [internal units]"
    ] = snap.metadata.cosmology_raw["Lookback time [internal units]"]
    extra_header["snapshot cosmology N_eff"] = snap.metadata.cosmology_raw["N_eff"]
    extra_header["snapshot cosmology N_nu"] = snap.metadata.cosmology_raw["N_nu"]
    extra_header["snapshot cosmology N_ur"] = snap.metadata.cosmology_raw["N_ur"]
    extra_header["snapshot cosmology Omega_b"] = snap.metadata.cosmology_raw["Omega_b"]
    extra_header["snapshot cosmology Omega_cdm"] = snap.metadata.cosmology_raw[
        "Omega_cdm"
    ]
    extra_header["snapshot cosmology Omega_g"] = snap.metadata.cosmology_raw["Omega_g"]
    extra_header["snapshot cosmology Omega_k"] = snap.metadata.cosmology_raw["Omega_k"]
    extra_header["snapshot cosmology Omega_lambda"] = snap.metadata.cosmology_raw[
        "Omega_lambda"
    ]
    extra_header["snapshot cosmology Omega_m"] = snap.metadata.cosmology_raw["Omega_m"]
    extra_header["snapshot cosmology Omega_nu"] = snap.metadata.cosmology_raw[
        "Omega_nu"
    ]
    extra_header["snapshot cosmology Omega_nu_0"] = snap.metadata.cosmology_raw[
        "Omega_nu_0"
    ]
    extra_header["snapshot cosmology Omega_r"] = snap.metadata.cosmology_raw["Omega_r"]
    extra_header["snapshot cosmology Omega_ur"] = snap.metadata.cosmology_raw[
        "Omega_ur"
    ]
    extra_header["snapshot cosmology T_CMB_0 [K]"] = snap.metadata.cosmology_raw[
        "T_CMB_0 [K]"
    ]
    extra_header[
        "snapshot cosmology T_CMB_0 [internal_units]"
    ] = snap.metadata.cosmology_raw["T_CMB_0 [internal units]"]
    extra_header["snapshot cosmology T_nu_0 [eV]"] = snap.metadata.cosmology_raw[
        "T_nu_0 [eV]"
    ]
    extra_header[
        "snapshot cosmology T_nu_0 [internal_units]"
    ] = snap.metadata.cosmology_raw["T_nu_0 [internal units]"]
    extra_header["snapshot cosmology h"] = snap.metadata.cosmology_raw["h"]
    extra_header["snapshot cosmology w"] = snap.metadata.cosmology_raw["w"]
    extra_header["snapshot cosmology w_0"] = snap.metadata.cosmology_raw["w_0"]
    extra_header["snapshot cosmology w_a"] = snap.metadata.cosmology_raw["w_a"]

else:
    print("Found no cosmology metadata in snapshot. Continuing without.")


extra_header["snapshot a"] = snap.metadata.a
extra_header["snapshot date"] = snap.metadata.snapshot_date
extra_header["snapshot dimension"] = snap.metadata.dimension
extra_header["snapshot filename"] = snap.metadata.filename
extra_header["snapshot gas gamma"] = snap.metadata.gas_gamma
extra_header["snapshot gravity scheme"] = snap.metadata.gravity_scheme
extra_header["snapshot header"] = snap.metadata.header
extra_header["snapshot hydro info"] = snap.metadata.hydro_info
extra_header["snapshot hydro scheme"] = snap.metadata.hydro_scheme
extra_header["snapshot initial mass table"] = snap.metadata.initial_mass_table
extra_header["snapshot library info"] = snap.metadata.library_info
extra_header["snapshot mass table"] = snap.metadata.mass_table
extra_header["snapshot redshift"] = snap.metadata.redshift
extra_header["snapshot run name"] = snap.metadata.run_name
extra_header["snapshot scale_factor"] = snap.metadata.scale_factor
extra_header["snapshot subgrid scheme"] = snap.metadata.subgrid_scheme
extra_header["snapshot stars properties"] = snap.metadata.stars_properties
extra_header["snapshot stars scheme"] = snap.metadata.stars_scheme
extra_header["snapshot system name"] = snap.metadata.system_name
extra_header["snapshot time"] = snap.metadata.time


# Initialize the Writer
writer = Writer(
    unit_system=my_units,
    box_size=snap.metadata.boxsize,
    dimension=snap.metadata.dimension,
    compress=True,
    extra_header=extra_header,
    scale_factor=snap.metadata.scale_factor,
)

#  /PartType0/ - Gas
#  /PartType1/ - Dark Matter
#  /PartType2/ - Background Dark Matter
#  /PartType3/ - Sinks
#  /PartType4/ - Stars
#  /PartType5/ - Black Holes
#  /PartType6/ - Neutrino Dark Matter

if snap.metadata.has_type[0]:
    # Get and write gas
    print("Adding gas data to ICs.")
    writer.gas.coordinates = snap.gas.coordinates
    writer.gas.velocities = snap.gas.velocities
    writer.gas.masses = snap.gas.masses
    writer.gas.internal_energy = snap.gas.internal_energies
    writer.gas.smoothing_length = snap.gas.smoothing_lengths
    writer.gas.particle_ids = snap.gas.particle_ids
else:
    print("Found no gas data in snapshot. Continuing without.")

if snap.metadata.has_type[1]:
    # Get and write dark matter
    print("Adding dark matter data to ICs.")
    writer.dark_matter.coordinates = snap.dark_matter.coordinates
    writer.dark_matter.velocities = snap.dark_matter.velocities
    writer.dark_matter.masses = snap.dark_matter.masses
    writer.dark_matter.particle_ids = snap.dark_matter.particle_ids
else:
    print("Found no dark matter data in snapshot. Continuing without.")

if snap.metadata.has_type[3]:
    # Get and write sinks
    print("Adding sinks data to ICs.")
    writer.sinks.coordinates = snap.sinks.coordinates
    writer.sinks.velocities = snap.sinks.velocities
    writer.sinks.masses = snap.sinks.masses
    writer.sinks.particle_ids = snap.sinks.particle_ids
else:
    print("Found no sinks data in snapshot. Continuing without.")

if snap.metadata.has_type[4]:
    # Get and write stars
    print("Adding stars data to ICs.")
    writer.stars.coordinates = snap.stars.coordinates
    writer.stars.velocities = snap.stars.velocities
    writer.stars.masses = snap.stars.masses
    writer.stars.smoothing_length = snap.stars.smoothing_lengths
    writer.stars.particle_ids = snap.stars.particle_ids
else:
    print("Found no star data in snapshot. Continuing without.")


if snap.metadata.has_type[5]:
    # Get and write black holes
    print("Adding black holes data to ICs.")
    writer.black_holes.coordinates = snap.black_holes.coordinates
    writer.black_holes.velocities = snap.black_holes.velocities
    writer.black_holes.masses = snap.black_holes.masses
    writer.black_holes.particle_ids = snap.black_holes.particle_ids
    writer.black_holes.smoothing_length = snap.black_holes.smoothing_lengths
else:
    print("Found no black hole data in snapshot. Continuing without.")


if snap.metadata.has_type[6]:
    # Get and write neutrinos
    # TODO
    raise NotImplementedError("Can't write neutrino ICs yet.")


# Finally, write everything.
writer.write(output_filename)
