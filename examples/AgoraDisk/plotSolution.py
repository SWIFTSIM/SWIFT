#!/usr/bin/env python2
#######################################################################
#
#  UNIFIED ANALYSIS SCRIPT FOR DISK SIMULATION FOR THE AGORA PROJECT
#
#  FOR SCRIPT HISTORY SEE VERSION CONTROL CHANGELOG
#
#  Note: This script requires yt-3.2 or yt-dev. Older versions may
#        yield incorrect results.
#
#######################################################################
# This script is a copy of the AGORA project (https://bitbucket.org/mornkr/agora-analysis-script/)
# modified in order to take into account SWIFT
import matplotlib
matplotlib.use('Agg')
import sys
import math
import copy
import csv
import numpy as np
import matplotlib.colorbar as cb
import matplotlib.lines as ln
import yt.utilities.physical_constants as constants
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from yt.mods import *
from yt.units.yt_array import YTQuantity
from mpl_toolkits.axes_grid1 import AxesGrid
from matplotlib.offsetbox import AnchoredText
from matplotlib.ticker import MaxNLocator
from yt.fields.particle_fields import add_volume_weighted_smoothed_field
from yt.fields.particle_fields import add_nearest_neighbor_field
from yt.analysis_modules.star_analysis.api import StarFormationRate
from yt.analysis_modules.halo_analysis.api import *
from yt.data_objects.particle_filters import add_particle_filter
from scipy.stats import kde
from subprocess import call
#mylog.setLevel(1)

import sys

if "-C" in sys.argv:
    with_cooling = 1
else:
    with_cooling = 0
    
draw_density_map                 = 1#1           # 0/1         = OFF/ON
draw_temperature_map             = 1#1           # 0/1         = OFF/ON
draw_cellsize_map                = 2#2           # 0/1/2/3     = OFF/ON/ON now with resolution map where particle_size is defined as [M/rho]^(1/3) for SPH/ON with both cell_size and resolution maps
draw_elevation_map               = 1#1           # 0/1         = OFF/ON
draw_metal_map                   = 0#0           # 0/1/2       = OFF/ON/ON with total metal mass in the simulation annotated (this will be inaccurate for SPH)
draw_zvel_map                    = 0#0           # 0/1         = OFF/ON
draw_star_map                    = 0#0           # 0/1         = OFF/ON
draw_star_clump_stats            = 0#0           # 0/1/2/3     = OFF/ON/ON with additional star_map with annotated clumps/ON with addtional star_map + extra clump mass function with GIZMO-ps2
draw_SFR_map                     = 0###0         # 0/1/2       = OFF/ON/ON with local K-S plot using patches
draw_PDF                         = 0###0         # 0/1/2/3     = OFF/ON/ON with constant pressure/entropy lines/ON with additional annotations such as 1D profile from a specific code (CHANGA)
draw_pos_vel_PDF                 = 0##0          # 0/1/2/3/4   = OFF/ON/ON with 1D profile/ON also with 1D dispersion profile/ON also with separate 1D vertical dispersion profiles
draw_star_pos_vel_PDF            = 0##0          # 0/1/2/3/4   = OFF/ON/ON with 1D profile/ON also with 1D dispersion profile/ON also with separate 1D vertical dispersion profiles
draw_rad_height_PDF              = 0##0          # 0/1/2/3     = OFF/ON/ON with 1D profile/ON with analytic ftn subtracted
draw_metal_PDF                   = 0##0          # 0/1         = OFF/ON
draw_density_DF                  = 0#0           # 0/1/2       = OFF/ON/ON with difference plot between 1st and 2nd datasets (when 2, dataset_num should be set to 2)
draw_radius_DF                   = 0#0           # 0/1         = OFF/ON
draw_star_radius_DF              = 0#0           # 0/1/2       = OFF/ON/ON with SFR profile and K-S plot (when 2, this automatically turns on draw_radius_DF)
draw_height_DF                   = 0#0           # 0/1         = OFF/ON
draw_SFR                         = 0#0           # 0/1/2       = OFF/ON/ON with extra line with GIZMO-ps2
draw_cut_through                 = 0#0           # 0/1         = OFF/ON
add_nametag                      = 1#1           # 0/1         = OFF/ON
add_mean_fractional_dispersion   = 0#0           # 0/1         = OFF/ON

dataset_num                      = 2             # 1/2         = 1st dataset(Grackle+noSF)/2nd dataset(Grackle+SF+ThermalFbck) for AGORA Paper 4
yt_version_pre_3_2_3             = 0             # 0/1         = NO/YES to "Is the yt version being used pre yt-dev-3.2.3?"
times                            = [0, 500]      # in Myr
figure_width                     = 30            # in kpc
n_ref                            = 64            # 8; for SPH codes
over_refine_factor               = 1             # 2; for SPH codes
aperture_size_SFR_map            = 750           # in pc       = Used if draw_SFR_map = 1, 750 matches Bigiel et al. (2008)
young_star_cutoff_SFR_map        = 20            # in Myr      = Used if draw_SFR_map = 1
young_star_cutoff_star_radius_DF = 20            # in Myr      = Used if draw_star_radius_DF = 1
mean_dispersion_radius_range     = [2, 10]       # in kpc      = Used if add_mean_fractional_dispersion = 1, for draw_pos_vel_PDF, draw_star_pos_vel_PDF, draw_rad_height_PDF, draw_radius_DF, draw_star_radius_DF
mean_dispersion_height_range     = [0, 0.6]      # in kpc      = Used if add_mean_fractional_dispersion = 1, for draw_height_DF
mean_dispersion_time_range       = [50, 500]     # in Myr      = Used if add_mean_fractional_dispersion = 1, for draw_SFR
mean_dispersion_density_range    = [1e-25, 1e-22]# in g/cm3    = Used if add_mean_fractional_dispersion = 1, for draw_density_DF
disk_normal_vector               = [0., 0., 1.]

gadget_default_unit_base = {'UnitLength_in_cm'         : 3.08567758e+21,
                            'UnitMass_in_g'            :    1.98848e+43,
                            'UnitVelocity_in_cm_per_s' :      100000}
color_names              = ['red', 'magenta', 'orange', 'gold', 'green', 'cyan', 'blue', 'blueviolet', 'black']
linestyle_names          = ['-']
marker_names             = ['s', 'o', 'p', 'v', '^', '<', '>', 'h', '*']

# file_location = ['/lustre/ki/pfs/mornkr/080112_CHaRGe/pfs-hyades/AGORA-DISK-repository-for-use/Grackle+noSF/',
#                '/lustre/ki/pfs/mornkr/080112_CHaRGe/pfs-hyades/AGORA-DISK-repository-for-use/Grackle+SF+ThermalFbck/']
# file_location = ['/global/project/projectdirs/agora/AGORA-DISK-repository-for-use/Grackle+noSF/',
#                '/global/project/projectdirs/agora/AGORA-DISK-repository-for-use/Grackle+SF+ThermalFbck/']
# codes = ['ART-I', 'ART-II', 'ENZO', 'RAMSES', 'CHANGA', 'GASOLINE', 'GADGET-3', 'GEAR', 'GIZMO']
# filenames = [[[file_location[0]+'ART-I/IC/AGORA_Galaxy_LOW.d', file_location[0]+'ART-I/t0.5Gyr/10MpcBox_csf512_02350.d'],
#             [file_location[0]+'ART-II/noSF_def_2p/OUT/AGORA_LOW_000000.art', file_location[0]+'ART-II/noSF_def_2p/OUT/AGORA_LOW_000098.art'],
#             [file_location[0]+'ENZO/DD0000/DD0000', file_location[0]+'ENZO/DD0100/DD0100'],
#             [file_location[0]+'RAMSES/output_00001/info_00001.txt', file_location[0]+'RAMSES/output_00068/info_00068.txt'],
#             [file_location[0]+'CHANGA/disklow/disklow.000000', file_location[0]+'CHANGA/disklow/disklow.000500'],
#             [file_location[0]+'GASOLINE/LOW_dataset1.00001', file_location[0]+'GASOLINE/LOW_dataset1.00335'],
#             [file_location[0]+'GADGET-3/AGORA_ISO_LOW_DRY/snap_iso_dry_000.hdf5', file_location[0]+'GADGET-3/AGORA_ISO_LOW_DRY/snap_iso_dry_010.hdf5'],
#             [file_location[0]+'GEAR/snapshot_0000', file_location[0]+'GEAR/snapshot_0500'],
#             [file_location[0]+'GIZMO/snapshot_temp_000', file_location[0]+'GIZMO/snapshot_temp_100']],
#            [[file_location[1]+'ART-I/IC/AGORA_Galaxy_LOW.d', file_location[1]+'ART-I/t0.5Gyr/10MpcBox_csf512_02350.d'],
#             [file_location[1]+'ART-II/SF_FBth_def_2p/OUT/AGORA_LOW_000000.art', file_location[1]+'ART-II/SF_FBth_def_2p/OUT/AGORA_LOW_000311.art'],
#             [file_location[1]+'ENZO/DD0000/DD0000', file_location[1]+'ENZO/DD0050/DD0050'],
#             [file_location[1]+'RAMSES/output_00001/info_00001.txt', file_location[1]+'RAMSES/output_00216/info_00216.txt'],
#             [file_location[1]+'CHANGA/disklow/disklow.000000', file_location[1]+'CHANGA/disklow/disklow.000500'],
#             [file_location[1]+'GASOLINE/LOW_dataset1.00001', file_location[1]+'GASOLINE/LOW_dataset2.00335'],
#             [file_location[1]+'GADGET-3/AGORA_ISO_LOW_SF_SNII_Thermal_Chevalier_SFT10/snap_iso_sf_000.hdf5', file_location[1]+'GADGET-3/AGORA_ISO_LOW_SF_SNII_Thermal_Chevalier_SFT10/snap_iso_sf_010.hdf5'],
#             [file_location[1]+'GEAR/snapshot_0000', file_location[1]+'GEAR/snapshot_0500'],
#             [file_location[1]+'GIZMO/snapshot_temp_000', file_location[1]+'GIZMO/snapshot_temp_100']]]
codes = ['SWIFT', 'GEAR']
filenames = [[["./agora_disk_IC.hdf5", "./agora_disk_500Myr.hdf5"],
              ["./snapshot_0000.hdf5", "./snapshot_0500.hdf5"]],
             [["./agora_disk_IC.hdf5", "./agora_disk_500Myr.hdf5"],
              ["./snapshot_0000.hdf5", "./snapshot_0500.hdf5"]]]

# codes = ["SWIFT"]
# filenames = [[["./agora_disk_0000.hdf5", "./agora_disk_0050.hdf5"]],
#             [["./agora_disk_0000.hdf5", "./agora_disk_0050.hdf5"]]]
# codes = ['ART-I']
# filenames = [[[file_location[0]+'ART-I/IC/AGORA_Galaxy_LOW.d', file_location[0]+'ART-I/t0.5Gyr/10MpcBox_csf512_02350.d']],
#            [[file_location[1]+'ART-I/IC/AGORA_Galaxy_LOW.d', file_location[1]+'ART-I/t0.5Gyr/10MpcBox_csf512_02350.d']]]
# codes = ['ART-II']
# filenames = [[[file_location[0]+'ART-II/noSF_def_2p/OUT/AGORA_LOW_000000.art', file_location[0]+'ART-II/noSF_def_2p/OUT/AGORA_LOW_000098.art']],
#            [[file_location[1]+'ART-II/SF_FBth_def_2p/OUT/AGORA_LOW_000000.art', file_location[1]+'ART-II/SF_FBth_def_2p/OUT/AGORA_LOW_000311.art']]]
# codes = ['ENZO']
# filenames = [[[file_location[0]+'ENZO/DD0000/DD0000', file_location[0]+'ENZO/DD0100/DD0100']],
#            [[file_location[1]+'ENZO/DD0000/DD0000', file_location[1]+'ENZO/DD0050/DD0050']]]
# codes = ['RAMSES']
# filenames = [[[file_location[0]+'RAMSES/output_00001/info_00001.txt', file_location[0]+'RAMSES/output_00068/info_00068.txt']],
#            [[file_location[1]+'RAMSES/output_00001/info_00001.txt', file_location[1]+'RAMSES/output_00216/info_00216.txt']]]
# codes = ['CHANGA']
# filenames = [[[file_location[0]+'CHANGA/disklow/disklow.000000', file_location[0]+'CHANGA/disklow/disklow.000500']],
#            [[file_location[1]+'CHANGA/disklow/disklow.000000', file_location[1]+'CHANGA/disklow/disklow.000500']]]
# codes = ['GASOLINE']
# filenames = [[[file_location[0]+'GASOLINE/LOW_dataset1.00001', file_location[0]+'GASOLINE/LOW_dataset1.00335']],
#            [[file_location[1]+'GASOLINE/LOW_dataset1.00001', file_location[1]+'GASOLINE/LOW_dataset2.00335']]]
# codes = ['GADGET-3']
# filenames = [[[file_location[0]+'GADGET-3/AGORA_ISO_LOW_DRY/snap_iso_dry_000.hdf5', file_location[0]+'GADGET-3/AGORA_ISO_LOW_DRY/snap_iso_dry_010.hdf5']],
#            [[file_location[1]+'GADGET-3/AGORA_ISO_LOW_SF_SNII_Thermal_Chevalier_SFT10/snap_iso_sf_000.hdf5', file_location[1]+'GADGET-3/AGORA_ISO_LOW_SF_SNII_Thermal_Chevalier_SFT10/snap_iso_sf_010.hdf5']]]
# codes = ['GEAR']
# filenames = [[['snapshot_0000', 'snapshot_0500']],
#            [['snapshot_0000', 'snapshot_0500']]]
# codes = ['GIZMO']
# filenames = [[[file_location[0]+'GIZMO/snapshot_temp_000', file_location[0]+'GIZMO/snapshot_temp_100']],
#            [[file_location[1]+'GIZMO/snapshot_temp_000', file_location[1]+'GIZMO/snapshot_temp_100']]]

def load_dataset(dataset_num, time, code, codes, filenames_entry):
        if codes[code] == 'ART-I': # ART frontend doesn't find the accompanying files, so we specify them; see http://yt-project.org/docs/dev/examining/loading_data.html#art-data
                dirnames = filenames_entry[code][time][:filenames_entry[code][time].rfind('/')+1]
                if time == 0:
                        timestamp = ''
                else:
                        timestamp = filenames_entry[code][time][filenames_entry[code][time].rfind('_'):filenames_entry[code][time].rfind('.')]
                pf = load(filenames_entry[code][time], file_particle_header=dirnames+'PMcrd'+timestamp+'.DAT', file_particle_data=dirnames+'PMcrs0'+timestamp+'.DAT', file_particle_stars=dirnames+'stars'+timestamp+'.dat')
        elif codes[code] == 'CHANGA' or codes[code] == 'GASOLINE': # For TIPSY frontend, always make sure to place your parameter file in the same directory as your datasets
                pf = load(filenames_entry[code][time], n_ref=n_ref, over_refine_factor=over_refine_factor)
        elif codes[code] == 'GADGET-3': # For GADGET-3 2nd dataset, there somehow exist particles very far from the center; so we use [-2000, 2000] for a bounding_box
                pf = load(filenames_entry[code][time], unit_base = gadget_default_unit_base, bounding_box=[[-2000.0, 2000.0], [-2000.0, 2000.0], [-2000.0, 2000.0]], n_ref=n_ref, over_refine_factor=over_refine_factor)
        elif codes[code] == "SWIFT":
                pf = load(filenames_entry[code][time], unit_base = gadget_default_unit_base, bounding_box=[[-2000.0, 2000.0], [-2000.0, 2000.0], [-2000.0, 2000.0]], n_ref=n_ref, over_refine_factor=over_refine_factor)
        elif codes[code] == 'GEAR':
                pf = load(filenames_entry[code][time], unit_base = gadget_default_unit_base, bounding_box=[[-2000.0, 2000.0], [-2000.0, 2000.0], [-2000.0, 2000.0]], n_ref=n_ref, over_refine_factor=over_refine_factor) 
               # from yt.frontends.gadget.definitions import gadget_header_specs
               #  from yt.frontends.gadget.definitions import gadget_ptype_specs
               #  from yt.frontends.gadget.definitions import gadget_field_specs
               #  if with_cooling:
               #          gadget_header_specs["chemistry"] = (('h1',  4, 'c'),('h2',  4, 'c'),('empty',  256, 'c'),)
               #          header_spec = "default+chemistry"
               #  else:
               #          header_spec = "default"
               #  gear_ptype_specs = ("Gas", "Stars", "Halo", "Disk", "Bulge", "Bndry")
               #  if dataset_num == 1:
               #          pf = GadgetDataset(filenames_entry[code][time], unit_base = gadget_default_unit_base, bounding_box=[[-1000.0, 1000.0], [-1000.0, 1000.0], [-1000.0, 1000.0]], header_spec=header_spec,
               #                             ptype_spec=gear_ptype_specs, n_ref=n_ref, over_refine_factor=over_refine_factor)
               #  elif dataset_num == 2:
               #          # For GEAR 2nd dataset (Grackle+SF+ThermalFbck), "Metals" is acutally 10-species field; check how metallicity field is added below.
               #          if with_cooling:
               #                  agora_gear  = ( "Coordinates",
               #                                  "Velocities",
               #                                  "ParticleIDs",
               #                                  "Mass",
               #                                  ("InternalEnergy", "Gas"),
               #                                  ("Density", "Gas"),
               #                                  ("SmoothingLength", "Gas"),
               #                                  ("Metals", "Gas"),
               #                                  ("StellarFormationTime", "Stars"),
               #                                  ("StellarInitMass", "Stars"),
               #                                  ("StellarIDs", "Stars"),
               #                                  ("StellarDensity", "Stars"),
               #                                  ("StellarSmoothingLength", "Stars"),
               #                                  ("StellarMetals", "Stars"),
               #                                  ("Opt1", "Stars"),
               #                                  ("Opt2", "Stars"),
               #                  )
               #          else:
               #                  agora_gear  = ( "Coordinates",
               #                                  "Velocities",
               #                                  "ParticleIDs",
               #                                  "Mass",
               #                                  ("InternalEnergy", "Gas"),
               #                                  ("Density", "Gas"),
               #                                  ("SmoothingLength", "Gas"),
               #                  )
               #          gadget_field_specs["agora_gear"] = agora_gear
               #          pf = GadgetDataset(filenames_entry[code][time], unit_base = gadget_default_unit_base, bounding_box=[[-1000.0, 1000.0], [-1000.0, 1000.0], [-1000.0, 1000.0]], header_spec=header_spec,
               #                             ptype_spec=gear_ptype_specs, field_spec="agora_gear", n_ref=n_ref, over_refine_factor=over_refine_factor)
        elif codes[code] == 'GIZMO':
                from yt.frontends.gadget.definitions import gadget_field_specs
                if dataset_num == 1:
                        agora_gizmo = ( "Coordinates",
                                        "Velocities",
                                        "ParticleIDs",
                                        "Mass",
                                        ("Temperature", "Gas"),
                                        ("Density", "Gas"),
                                        ("Electron_Number_Density", "Gas"),
                                        ("HI_NumberDensity", "Gas"),
                                        ("SmoothingLength", "Gas"),
                                        )
                elif dataset_num == 2:
                        agora_gizmo = ( "Coordinates",
                                        "Velocities",
                                        "ParticleIDs",
                                        "Mass",
                                        ("Temperature", "Gas"),
                                        ("Density", "Gas"),
                                        ("ElectronAbundance", "Gas"),
                                        ("NeutralHydrogenAbundance", "Gas"),
                                        ("SmoothingLength", "Gas"),
                                        ("StarFormationRate", "Gas"),
                                        ("StellarFormationTime", "Stars"),
                                        ("Metallicity", ("Gas", "Stars")),
                                        ("DelayTime", "Stars"),
                                        ("StellarInitMass", "Stars"),
                                        )
                gadget_field_specs["agora_gizmo"] = agora_gizmo
                pf = load(filenames_entry[code][time], unit_base = gadget_default_unit_base, bounding_box=[[-1000.0, 1000.0], [-1000.0, 1000.0], [-1000.0,1000.0]], field_spec="agora_gizmo",
                          n_ref=n_ref, over_refine_factor=over_refine_factor)
        else:
                pf = load(filenames_entry[code][time])
        return pf

fig_density_map             = []
fig_temperature_map         = []
fig_cellsize_map            = []
fig_cellsize_map_2          = []
fig_elevation_map           = []
fig_metal_map               = []
fig_zvel_map                = []
fig_star_map                = []
fig_star_map_2              = []
fig_star_map_3              = []
fig_degr_sfr_map            = []
fig_degr_density_map        = []
fig_PDF                     = []
fig_pos_vel_PDF             = []
fig_star_pos_vel_PDF        = []
fig_rad_height_PDF          = []
fig_metal_PDF               = []
grid_density_map            = []
grid_temperature_map        = []
grid_cellsize_map           = []
grid_cellsize_map_2         = []
grid_elevation_map          = []
grid_metal_map              = []
grid_zvel_map               = []
grid_star_map               = []
grid_star_map_2             = []
grid_star_map_3             = []
grid_degr_sfr_map           = []
grid_degr_density_map       = []
grid_PDF                    = []
grid_pos_vel_PDF            = []
grid_star_pos_vel_PDF       = []
grid_rad_height_PDF         = []
grid_metal_PDF              = []
star_clump_masses_hop       = []
star_clump_masses_fof       = []
star_clump_masses_hop_ref   = []
star_clump_masses_fof_ref   = []
pos_vel_xs                  = []
pos_vel_profiles            = []
pos_disp_xs                 = []
pos_disp_profiles           = []
pos_disp_vert_xs            = []
pos_disp_vert_profiles      = []
star_pos_vel_xs             = []
star_pos_vel_profiles       = []
star_pos_disp_xs            = []
star_pos_disp_profiles      = []
star_pos_disp_vert_xs       = []
star_pos_disp_vert_profiles = []
rad_height_xs               = []
rad_height_profiles         = []
density_DF_xs               = []
density_DF_profiles         = []
density_DF_1st_xs           = []
density_DF_1st_profiles     = []
radius_DF_xs                = []
radius_DF_profiles          = []
surface_density             = []
star_radius_DF_xs           = []
star_radius_DF_profiles     = []
star_surface_density        = []
sfr_radius_DF_xs            = []
sfr_radius_DF_profiles      = []
sfr_surface_density         = []
height_DF_xs                = []
height_DF_profiles          = []
height_surface_density      = []
sfr_ts                      = []
sfr_cum_masses              = []
sfr_sfrs                    = []
surf_dens_SFR_map           = []
sfr_surf_dens_SFR_map       = []
cut_through_zs              = []
cut_through_zvalues         = []
cut_through_xs              = []
cut_through_xvalues         = []

for time in range(len(times)):
        if draw_density_map == 1:
                fig_density_map      += [plt.figure(figsize=(100,20))]
                grid_density_map     += [AxesGrid(fig_density_map[time], (0.01,0.01,0.99,0.99), nrows_ncols = (2, len(codes)), axes_pad = 0.02, add_all = True, share_all = True,
                                                  label_mode = "1", cbar_mode = "single", cbar_location = "right", cbar_size = "2%", cbar_pad = 0.02)]
        if draw_temperature_map == 1:
                fig_temperature_map  += [plt.figure(figsize=(100,20))]
                grid_temperature_map += [AxesGrid(fig_temperature_map[time], (0.01,0.01,0.99,0.99), nrows_ncols = (2, len(codes)), axes_pad = 0.02, add_all = True, share_all = True,
                                                  label_mode = "1", cbar_mode = "single", cbar_location = "right", cbar_size = "2%", cbar_pad = 0.02)]
        if draw_cellsize_map >= 1:
                fig_cellsize_map     += [plt.figure(figsize=(100,20))]
                grid_cellsize_map    += [AxesGrid(fig_cellsize_map[time], (0.01,0.01,0.99,0.99), nrows_ncols = (2, len(codes)), axes_pad = 0.02, add_all = True, share_all = True,
                                                  label_mode = "1", cbar_mode = "single", cbar_location = "right", cbar_size = "2%", cbar_pad = 0.02)]
                fig_cellsize_map_2   += [plt.figure(figsize=(100,20))]
                grid_cellsize_map_2  += [AxesGrid(fig_cellsize_map_2[time], (0.01,0.01,0.99,0.99), nrows_ncols = (2, len(codes)), axes_pad = 0.02, add_all = True, share_all = True,
                                                  label_mode = "1", cbar_mode = "single", cbar_location = "right", cbar_size = "2%", cbar_pad = 0.02)]
        if draw_elevation_map == 1:
                fig_elevation_map    += [plt.figure(figsize=(100,20))]
                grid_elevation_map   += [AxesGrid(fig_elevation_map[time], (0.01,0.01,0.99,0.99), nrows_ncols = (1, len(codes)), axes_pad = 0.02, add_all = True, share_all = True,
                                                  label_mode = "1", cbar_mode = "single", cbar_location = "right", cbar_size = "6%", cbar_pad = 0.02)]
        if draw_metal_map >= 1:
                fig_metal_map        += [plt.figure(figsize=(100,20))]
                grid_metal_map       += [AxesGrid(fig_metal_map[time], (0.01,0.01,0.99,0.99), nrows_ncols = (2, len(codes)), axes_pad = 0.02, add_all = True, share_all = True,
                                                  label_mode = "1", cbar_mode = "single", cbar_location = "right", cbar_size = "2%", cbar_pad = 0.02)]
        if draw_zvel_map == 1:
                fig_zvel_map         += [plt.figure(figsize=(100,20))]
                grid_zvel_map        += [AxesGrid(fig_zvel_map[time], (0.01,0.01,0.99,0.99), nrows_ncols = (2, len(codes)), axes_pad = 0.02, add_all = True, share_all = True,
                                                  label_mode = "1", cbar_mode = "single", cbar_location = "right", cbar_size = "2%", cbar_pad = 0.02)]
        if draw_star_map == 1:
                fig_star_map         += [plt.figure(figsize=(100,20))]
                grid_star_map        += [AxesGrid(fig_star_map[time], (0.01,0.01,0.99,0.99), nrows_ncols = (2, len(codes)), axes_pad = 0.02, add_all = True, share_all = True,
                                                  label_mode = "1", cbar_mode = "single", cbar_location = "right", cbar_size = "2%", cbar_pad = 0.02)]
        if draw_star_clump_stats >= 1:
                star_clump_masses_hop.append([])
                star_clump_masses_fof.append([])
                fig_star_map_2       += [plt.figure(figsize=(100,20))]
                grid_star_map_2      += [AxesGrid(fig_star_map_2[time], (0.01,0.01,0.99,0.99), nrows_ncols = (2, len(codes)), axes_pad = 0.02, add_all = True, share_all = True,
                                                  label_mode = "1", cbar_mode = "single", cbar_location = "right", cbar_size = "2%", cbar_pad = 0.02)]
                fig_star_map_3       += [plt.figure(figsize=(100,20))]
                grid_star_map_3      += [AxesGrid(fig_star_map_3[time], (0.01,0.01,0.99,0.99), nrows_ncols = (2, len(codes)), axes_pad = 0.02, add_all = True, share_all = True,
                                                  label_mode = "1", cbar_mode = "single", cbar_location = "right", cbar_size = "2%", cbar_pad = 0.02)]
        if draw_SFR_map >= 1:
                fig_degr_density_map += [plt.figure(figsize=(100,20))]
                grid_degr_density_map+= [AxesGrid(fig_degr_density_map[time], (0.01,0.01,0.99,0.99), nrows_ncols = (1, len(codes)), axes_pad = 0.02, add_all = True, share_all = True,
                                                  label_mode = "1", cbar_mode = "single", cbar_location = "right", cbar_size = "6%", cbar_pad = 0.02)]
                fig_degr_sfr_map     += [plt.figure(figsize=(100,20))]
                grid_degr_sfr_map    += [AxesGrid(fig_degr_sfr_map[time], (0.01,0.01,0.99,0.99), nrows_ncols = (1, len(codes)), axes_pad = 0.02, add_all = True, share_all = True,
                                                  label_mode = "1", cbar_mode = "single", cbar_location = "right", cbar_size = "6%", cbar_pad = 0.02)]
                surf_dens_SFR_map.append([])
                sfr_surf_dens_SFR_map.append([])
        if draw_SFR_map >= 2 or draw_star_radius_DF >= 2:
                def import_text(filename, separator):
                        for line in csv.reader(open(filename), delimiter=separator, skipinitialspace=True):
                                if line:
                                        yield line
        if draw_PDF >= 1:
                fig_PDF              += [plt.figure(figsize=(50, 80))]
                grid_PDF             += [AxesGrid(fig_PDF[time], (0.01,0.01,0.99,0.99), nrows_ncols = (3, int(math.ceil(len(codes)/3.0))), axes_pad = 0.05, add_all = True, share_all = True,
                                                  label_mode = "1", cbar_mode = "single", cbar_location = "right", cbar_size = "2%", cbar_pad = 0.05, aspect = False)]
        if draw_pos_vel_PDF >= 1:
                fig_pos_vel_PDF      += [plt.figure(figsize=(50, 80))]
                grid_pos_vel_PDF     += [AxesGrid(fig_pos_vel_PDF[time], (0.01,0.01,0.99,0.99), nrows_ncols = (3, int(math.ceil(len(codes)/3.0))), axes_pad = 0.05, add_all = True, share_all = True,
                                                  label_mode = "1", cbar_mode = "single", cbar_location = "right", cbar_size = "2%", cbar_pad = 0.05, aspect = False)]
                pos_vel_xs.append([])
                pos_vel_profiles.append([])
                pos_disp_xs.append([])
                pos_disp_profiles.append([])
                pos_disp_vert_xs.append([])
                pos_disp_vert_profiles.append([])
        if draw_star_pos_vel_PDF >= 1:
                fig_star_pos_vel_PDF += [plt.figure(figsize=(50, 80))]
                grid_star_pos_vel_PDF+= [AxesGrid(fig_star_pos_vel_PDF[time], (0.01,0.01,0.99,0.99), nrows_ncols = (3, int(math.ceil(len(codes)/3.0))), axes_pad = 0.05, add_all = True, share_all = True,
                                                  label_mode = "1", cbar_mode = "single", cbar_location = "right", cbar_size = "2%", cbar_pad = 0.05, aspect = False)]
                star_pos_vel_xs.append([])
                star_pos_vel_profiles.append([])
                star_pos_disp_xs.append([])
                star_pos_disp_profiles.append([])
                star_pos_disp_vert_xs.append([])
                star_pos_disp_vert_profiles.append([])
        if draw_rad_height_PDF >= 1:
                fig_rad_height_PDF   += [plt.figure(figsize=(50, 80))]
                grid_rad_height_PDF  += [AxesGrid(fig_rad_height_PDF[time], (0.01,0.01,0.99,0.99), nrows_ncols = (3, int(math.ceil(len(codes)/3.0))), axes_pad = 0.05, add_all = True, share_all = True,
                                                  label_mode = "1", cbar_mode = "single", cbar_location = "right", cbar_size = "2%", cbar_pad = 0.05, aspect = False)]
                rad_height_xs.append([])
                rad_height_profiles.append([])
        if draw_metal_PDF == 1:
                fig_metal_PDF        += [plt.figure(figsize=(50, 80))]
                grid_metal_PDF       += [AxesGrid(fig_metal_PDF[time], (0.01,0.01,0.99,0.99), nrows_ncols = (3, int(math.ceil(len(codes)/3.0))), axes_pad = 0.05, add_all = True, share_all = True,
                                                  label_mode = "1", cbar_mode = "single", cbar_location = "right", cbar_size = "2%", cbar_pad = 0.05, aspect = False)]
        if draw_density_DF >= 1:
                density_DF_xs.append([])
                density_DF_profiles.append([])
                density_DF_1st_xs.append([])
                density_DF_1st_profiles.append([])
                if draw_density_DF == 2 and dataset_num == 1:
                        print("This won't work; resetting draw_density_DF to 1...")
                        draw_density_DF == 1
        if draw_star_radius_DF >= 1:
                star_radius_DF_xs.append([])
                star_radius_DF_profiles.append([])
                star_surface_density.append([])
                sfr_radius_DF_xs.append([])
                sfr_radius_DF_profiles.append([])
                sfr_surface_density.append([])
        if draw_star_radius_DF == 2:
                draw_radius_DF = 1
        if draw_radius_DF == 1:
                radius_DF_xs.append([])
                radius_DF_profiles.append([])
                surface_density.append([])
        if draw_height_DF == 1:
                height_DF_xs.append([])
                height_DF_profiles.append([])
                height_surface_density.append([])
        if draw_SFR >= 1:
                sfr_ts.append([])
                sfr_cum_masses.append([])
                sfr_sfrs.append([])
        if draw_cut_through == 1:
                cut_through_zs.append([])
                cut_through_zvalues.append([])
                cut_through_xs.append([])
                cut_through_xvalues.append([])

for time in range(len(times)):

        for code in range(len(codes)):

                if (dataset_num != 1 and dataset_num != 2) or (filenames[dataset_num-1][code] == []):
                        continue

                ####################################
                #        PRE-ANALYSIS STEPS        #
                ####################################

                # LOAD DATASETS
                pf = load_dataset(dataset_num, time, code, codes, filenames[dataset_num-1])

                # PARTICLE FILED NAMES FOR SPH CODES, AND STELLAR PARTICLE FILTERS
                PartType_Gas_to_use = "Gas"
                PartType_Star_to_use = "Stars"
                PartType_StarBeforeFiltered_to_use = "Stars"
                MassType_to_use = "Mass"
                MetallicityType_to_use = "Metallicity"
                FormationTimeType_to_use = "StellarFormationTime" # for GADGET/GEAR/GIZMO, this field has to be added in frontends/sph/fields.py, in which only "FormationTime" can be recognized
                SmoothingLengthType_to_use = "SmoothingLength"

                if codes[code] == 'CHANGA' or codes[code] == 'GASOLINE':
                        MetallicityType_to_use = "Metals"
                        PartType_Star_to_use = "NewStars"
                        PartType_StarBeforeFiltered_to_use = "Stars"
                        FormationTimeType_to_use = "FormationTime"
                        SmoothingLengthType_to_use = "smoothing_length"
                        def NewStars(pfilter, data): # see http://yt-project.org/docs/dev/analyzing/filtering.html#filtering-particle-fields
                                return (data[(pfilter.filtered_type, FormationTimeType_to_use)] > 0)
                        add_particle_filter(PartType_Star_to_use, function=NewStars, filtered_type=PartType_StarBeforeFiltered_to_use, requires=[FormationTimeType_to_use])
                        pf.add_particle_filter(PartType_Star_to_use)
                        pf.periodicity = (True, True, True) # this is needed especially when bPeriodic = 0 in GASOLINE, to avoid RuntimeError in geometry/selection_routines.pyx:855
                elif codes[code] == 'ART-I':
                        PartType_StarBeforeFiltered_to_use = "stars"
                        FormationTimeType_to_use = "particle_creation_time"
                        def Stars(pfilter, data):
                                return (data[(pfilter.filtered_type, FormationTimeType_to_use)] > 0)
#                               return ((data[(pfilter.filtered_type, FormationTimeType_to_use)] > 0) & (data[(pfilter.filtered_type, "particle_index")] >= 212500)) # without the fix metioned above, the above line won't work because all "stars"="specie1" have the same wrong particle_creation_time of 6.851 Gyr; so you will have to use this quick and dirty patch
                        add_particle_filter(PartType_Star_to_use, function=Stars, filtered_type=PartType_StarBeforeFiltered_to_use, requires=[FormationTimeType_to_use, "particle_index"])
                        pf.add_particle_filter(PartType_Star_to_use)
                elif codes[code] == 'ART-II':
                        PartType_StarBeforeFiltered_to_use = "STAR"
                        if time != 0: # BIRTH_TIME field exists in ART-II but as a dimensionless quantity for some reason in frontends/artio/fields.py; so we create StellarFormationTime field
                                def _FormationTime(field, data):
                                        return pf.arr(data["STAR", "BIRTH_TIME"].d, 'code_time')
                                pf.add_field(("STAR", FormationTimeType_to_use), function=_FormationTime, particle_type=True, take_log=False, units="code_time")
                        def Stars(pfilter, data):
                                return (data[(pfilter.filtered_type, FormationTimeType_to_use)] > 0)
                        add_particle_filter(PartType_Star_to_use, function=Stars, filtered_type=PartType_StarBeforeFiltered_to_use, requires=[FormationTimeType_to_use])
                        pf.add_particle_filter(PartType_Star_to_use)
                elif codes[code] == 'ENZO':
                        PartType_StarBeforeFiltered_to_use = "all"
                        FormationTimeType_to_use = "creation_time"
                        def Stars(pfilter, data):
                                return ((data[(pfilter.filtered_type, "particle_type")] == 2) & (data[(pfilter.filtered_type, FormationTimeType_to_use)] > 0))
                        add_particle_filter(PartType_Star_to_use, function=Stars, filtered_type=PartType_StarBeforeFiltered_to_use, requires=["particle_type", FormationTimeType_to_use])
                        pf.add_particle_filter(PartType_Star_to_use)
                elif codes[code] == "GADGET-3":
                        PartType_Gas_to_use = "PartType0"
                        PartType_Star_to_use = "PartType4"
                        PartType_StarBeforeFiltered_to_use = "PartType4"
                        MassType_to_use = "Masses"
                elif codes[code] == "SWIFT":
                        PartType_Gas_to_use = "PartType0"
                        PartType_Star_to_use = "PartType2"
                        PartType_StarBeforeFiltered_to_use = "PartType2"
                        MassType_to_use = "Masses"
                elif codes[code] == "GEAR":
                        PartType_Gas_to_use = "PartType0"
                        PartType_Star_to_use = "PartType1"
                        PartType_StarBeforeFiltered_to_use = "PartType1"
                        MassType_to_use = "Masses"
                elif codes[code] == 'RAMSES':
                        PartType_StarBeforeFiltered_to_use = "all"
                        pf.current_time = pf.arr(pf.parameters['time'], 'code_time') # reset pf.current_time because it is incorrectly set up in frontends/ramses/data_structure.py, and I don't wish to mess with units there
                        FormationTimeType_to_use = "particle_age" # particle_age field actually means particle creation time, at least for this particular dataset, so the new field below is not needed
                        # if time != 0: # Only particle_age field exists in RAMSES (only for new stars + IC stars), so we create StellarFormationTime field
                        #       def _FormationTime(field, data):
                        #               return pf.current_time - data["all", "particle_age"].in_units("s")
                        #       pf.add_field(("all", FormationTimeType_to_use), function=_FormationTime, particle_type=True, take_log=False, units="code_time")
                        def Stars(pfilter, data):
                                return (data[(pfilter.filtered_type, "particle_age")] > 0)
                        add_particle_filter(PartType_Star_to_use, function=Stars, filtered_type=PartType_StarBeforeFiltered_to_use, requires=["particle_age"])
                        pf.add_particle_filter(PartType_Star_to_use)

                # AXIS SWAP FOR PLOT COLLECTION
                pf.coordinates.x_axis[1] = 0
                pf.coordinates.y_axis[1] = 2
                pf.coordinates.x_axis['y'] = 0
                pf.coordinates.y_axis['y'] = 2

                # ADDITIONAL FIELDS I: GLOBALLY USED FIELDS
                def _density_squared(field, data):
                        return data[("gas", "density")]**2
                pf.add_field(("gas", "density_squared"), function=_density_squared, units="g**2/cm**6")

                # ADDITIONAL FIELDS II: FOR CELL SIZE AND RESOLUTION MAPS
                def _CellSizepc(field,data):
                        return (data[("index", "cell_volume")].in_units('pc**3'))**(1/3.)
                pf.add_field(("index", "cell_size"), function=_CellSizepc, units='pc', display_name="Resolution $\Delta$ x", take_log=True )
                def _Inv2CellVolumeCode(field,data):
                        return data[("index", "cell_volume")]**-2
                pf.add_field(("index", "cell_volume_inv2"), function=_Inv2CellVolumeCode, units='code_length**(-6)', display_name="Inv2CellVolumeCode", take_log=True)
                if draw_cellsize_map >= 2:
                        if codes[code] == 'CHANGA' or codes[code] == 'GEAR' or codes[code] == 'GADGET-3' or codes[code] == 'GASOLINE' or codes[code] == 'GIZMO' or codes[code] == "SWIFT":
                                def _ParticleSizepc(field, data):
                                        return (data[(PartType_Gas_to_use, MassType_to_use)]/data[(PartType_Gas_to_use, "Density")])**(1./3.)
                                pf.add_field((PartType_Gas_to_use, "particle_size"), function=_ParticleSizepc, units="pc", display_name="Resolution $\Delta$ x", particle_type=True, take_log=True)
                                def _Inv2ParticleVolumepc(field, data):
                                        return (data[(PartType_Gas_to_use, MassType_to_use)]/data[(PartType_Gas_to_use, "Density")])**(-2.)
                                pf.add_field((PartType_Gas_to_use, "particle_volume_inv2"), function=_Inv2ParticleVolumepc, units="pc**(-6)", display_name="Inv2ParticleVolumepc", particle_type=True, take_log=True)
                                # Also creating smoothed field following an example in yt-project.org/docs/dev/cookbook/calculating_information.html; use hardcoded num_neighbors as in frontends/gadget/fields.py
                                fn = add_volume_weighted_smoothed_field(PartType_Gas_to_use, "Coordinates", MassType_to_use, SmoothingLengthType_to_use, "Density", "particle_size", pf.field_info, nneighbors=64)
                                fn = add_volume_weighted_smoothed_field(PartType_Gas_to_use, "Coordinates", MassType_to_use, SmoothingLengthType_to_use, "Density", "particle_volume_inv2", pf.field_info, nneighbors=64)

                                # Alias doesn't work -- e.g. pf.field_info.alias(("gas", "particle_size"), fn[0]) -- check alias below; so I simply add ("gas", "particle_size")
                                def _ParticleSizepc_2(field, data):
                                        return data["deposit", PartType_Gas_to_use+"_smoothed_"+"particle_size"]
                                pf.add_field(("gas", "particle_size"), function=_ParticleSizepc_2, units="pc", force_override=True, display_name="Resolution $\Delta$ x", particle_type=False, take_log=True)
                                def _Inv2ParticleVolumepc_2(field, data):
                                        return data["deposit", PartType_Gas_to_use+"_smoothed_"+"particle_volume_inv2"]
                                pf.add_field(("gas", "particle_volume_inv2"), function=_Inv2ParticleVolumepc_2, units="pc**(-6)", force_override=True, display_name="Inv2ParticleVolumepc", particle_type=False, take_log=True)

                # ADDITIONAL FIELDS III: TEMPERATURE
                if codes[code] == 'GEAR' or codes[code] == 'GADGET-3' or codes[code] == 'RAMSES' or codes[code] == "SWIFT":
                        # From grackle/src/python/utilities/convenience.py: Calculate a tabulated approximation to mean molecular weight (valid for data that used Grackle 2.0 or below)
                        def calc_mu_table_local(temperature):
                                tt = np.array([1.0e+01, 1.0e+02, 1.0e+03, 1.0e+04, 1.3e+04, 2.1e+04, 3.4e+04, 6.3e+04, 1.0e+05, 1.0e+09])
                                mt = np.array([1.18701555, 1.15484424, 1.09603514, 0.9981496, 0.96346395, 0.65175895, 0.6142901, 0.6056833, 0.5897776, 0.58822635])
                                logttt= np.log(temperature)
                                logmu = np.interp(logttt,np.log(tt),np.log(mt)) # linear interpolation in log-log space
                                return np.exp(logmu)
                        temperature_values = []
                        mu_values = []
                        T_over_mu_values = []
                        current_temperature = 1e1
                        final_temperature = 1e7
                        dlogT = 0.1
                        while current_temperature < final_temperature:
                                temperature_values.append(current_temperature)
                                current_mu = calc_mu_table_local(current_temperature)
                                mu_values.append(current_mu)
                                T_over_mu_values.append(current_temperature/current_mu)
                                current_temperature = np.exp(np.log(current_temperature)+dlogT)
                        def convert_T_over_mu_to_T(T_over_mu):
                                logT_over_mu = np.log(T_over_mu)
                                logT = np.interp(logT_over_mu, np.log(T_over_mu_values), np.log(temperature_values)) # linear interpolation in log-log space
                                return np.exp(logT)
                        if codes[code] == 'GEAR' or codes[code] == 'GADGET-3' or codes[code] == "SWIFT":
                                def _Temperature_3(field, data):
                                        gamma = 5.0/3.0
                                        T_over_mu = (data[PartType_Gas_to_use, "InternalEnergy"] * (gamma-1) * constants.mass_hydrogen_cgs / constants.boltzmann_constant_cgs).in_units('K').d # T/mu
                                        return YTArray(convert_T_over_mu_to_T(T_over_mu), 'K') # now T
                                pf.add_field((PartType_Gas_to_use, "Temperature"), function=_Temperature_3, particle_type=True, force_override=True, units="K")
                        elif codes[code] == 'RAMSES':
                                # The pressure field includes the artificial pressure support term, so one needs to be careful (compare with the exsiting frontends/ramses/fields.py)
                                def _temperature_3(field, data):
                                        T_J = 1800.0  # in K
                                        n_J = 8.0     # in H/cc
                                        gamma_0 = 2.0
                                        x_H = 0.76
                                        mH = 1.66e-24      # from pymses/utils/constants/__init__.py  (vs. in yt, mass_hydrogen_cgs = 1.007947*amu_cgs = 1.007947*1.660538921e-24 = 1.6737352e-24)
                                        kB = 1.3806504e-16 # from pymses/utils/constants/__init__.py  (vs. in yt, boltzmann_constant_cgs = 1.3806488e-16)
                                        if time != 0:
                                                T_over_mu = data["gas", "pressure"].d/data["gas", "density"].d * mH / kB - T_J * (data["gas", "density"].d * x_H / mH / n_J)**(gamma_0 - 1.0) # T/mu = T2 in Ramses
                                        else:
                                                T_over_mu = data["gas", "pressure"].d/data["gas", "density"].d * mH / kB # IC: no pressure support
                                        return YTArray(convert_T_over_mu_to_T(T_over_mu), 'K') # now T
                                pf.add_field(("gas", "temperature"), function=_temperature_3, force_override=True, units="K")

                # ADDITIONAL FIELDS IV: METALLICITY (IN MASS FRACTION, NOT IN ZSUN)
                if draw_metal_map >= 1 or draw_metal_PDF == 1:
                        if codes[code] == 'ART-I': # metallicity field in ART-I has a different meaning (see frontends/art/fields.py), and metallicity field in ART-II is missing
                                def _metallicity_2(field, data):
                                        return (data["gas", "metal_ii_density"] + data["gas", "metal_ia_density"]) / data["gas", "density"] # metal_ia_density needs to be added to account for initial 0.5 Zsun, even though we don't have SNe Ia
                                pf.add_field(("gas", "metallicity"), function=_metallicity_2, force_override=True, display_name="Metallicity", take_log=True, units="")
                        elif codes[code] == 'ART-II':
                                def _metallicity_2(field, data):
                                        return data["gas", "metal_ii_density"] / data["gas", "density"]
                                pf.add_field(("gas", "metallicity"), function=_metallicity_2, force_override=True, display_name="Metallicity", take_log=True, units="")
                        elif codes[code] == 'ENZO': # metallicity field in ENZO is in Zsun, so we create a new field
                                def _metallicity_2(field, data):
                                        return data["gas", "metal_density"] / data["gas", "density"]
                                pf.add_field(("gas", "metallicity"), function=_metallicity_2, force_override=True, display_name="Metallicity", take_log=True, units="")
                        elif codes[code] == 'GEAR': # "Metals" in GEAR is 10-species field ([:,9] is the total metal fraction), so requires a change in _vector_fields in frontends/gadget/io.py: added ("Metals", 10)
                                def _metallicity_2(field, data):
                                        if len(data[PartType_Gas_to_use, "Metals"].shape) == 1:
                                                return data[PartType_Gas_to_use, "Metals"]
                                        else:
                                                return data[PartType_Gas_to_use, "Metals"][:,9].in_units("") # in_units("") turned out to be crucial!; otherwise code_metallicity will be used and it will mess things up
                                # We are creating ("Gas", "Metallicity") here, different from ("Gas", "metallicity") which is auto-generated by yt but doesn't work properly
                                pf.add_field((PartType_Gas_to_use, MetallicityType_to_use), function=_metallicity_2, display_name="Metallicity", particle_type=True, take_log=True, units="")
                                # Also creating smoothed field following an example in yt-project.org/docs/dev/cookbook/calculating_information.html; use hardcoded num_neighbors as in frontends/gadget/fields.py
                                fn = add_volume_weighted_smoothed_field(PartType_Gas_to_use, "Coordinates", MassType_to_use, "SmoothingLength", "Density", MetallicityType_to_use, pf.field_info, nneighbors=64)
                                # Alias doesn't work -- e.g. pf.field_info.alias(("gas", "metallicity"), fn[0]) -- probably because pf=GadgetDataset()?, not load()?; so I add and replace existing ("gas", "metallicity")
                                def _metallicity_3(field, data):
                                        return data["deposit", PartType_Gas_to_use+"_smoothed_"+MetallicityType_to_use]
                                pf.add_field(("gas", "metallicity"), function=_metallicity_3, force_override=True, display_name="Metallicity", particle_type=False, take_log=True, units="")
                        if draw_metal_map >= 2:
                                def _metal_mass(field, data):
                                        return data["gas", "cell_mass"] * data["gas", "metallicity"]
                                pf.add_field(("gas", "metal_mass"), function=_metal_mass, force_override=True, display_name="Metal Mass", take_log=True, units="Msun")

                # ADDITIONAL FIELDS V: FAKE PARTICLE FIELDS, CYLINDRICAL COORDINATES, etc.
                def rho_agora_disk(r, z):
                        r_d = YTArray(3.432, 'kpc')
                        z_d = 0.1*r_d
                        M_d = YTArray(4.297e10, 'Msun')
                        f_gas = 0.2
                        r_0 = (f_gas*M_d / (4.0*np.pi*r_d**2*z_d)).in_units('g/cm**3')
                        return r_0*numpy.exp(-r/r_d)*numpy.exp(-z/z_d)

                if codes[code] == "ART-I" or codes[code] == "ART-II" or codes[code] == "ENZO"  or codes[code] == "RAMSES":
                        def _cylindrical_z_abs(field, data):
                                return numpy.abs(data[("index", "cylindrical_z")])
                        pf.add_field(("index", "cylindrical_z_abs"), function=_cylindrical_z_abs, take_log=False, particle_type=False, units="cm")
                        def _density_minus_analytic(field, data):
                                return data[("gas", "density")] - rho_agora_disk(data[("index", "cylindrical_r")], data[("index", "cylindrical_z_abs")])
                        pf.add_field(("gas", "density_minus_analytic"), function=_density_minus_analytic, take_log=False, particle_type=False, display_name="density residual abs", units="g/cm**3")
                else:
                        def _particle_position_cylindrical_z_abs(field, data):
                                return numpy.abs(data[(PartType_Gas_to_use, "particle_position_cylindrical_z")])
                        pf.add_field((PartType_Gas_to_use, "particle_position_cylindrical_z_abs"), function=_particle_position_cylindrical_z_abs, take_log=False, particle_type=True, units="cm")
                        # particle_type=False doesn't make sense, but is critical for PhasePlot/ProfilePlot to work for particle fields
                        # requires a change in data_objects/data_container.py: remove raise YTFieldTypeNotFound(ftype)
                        def _Density_2(field, data):
                                return data[(PartType_Gas_to_use, "Density")].in_units('g/cm**3')
                        pf.add_field((PartType_Gas_to_use, "Density_2"), function=_Density_2, take_log=True, particle_type=False, display_name="Density", units="g/cm**3")
                        def _Temperature_2(field, data):
                                return data[(PartType_Gas_to_use, "Temperature")].in_units('K')
                        pf.add_field((PartType_Gas_to_use, "Temperature_2"), function=_Temperature_2, take_log=True, particle_type=False, display_name="Temperature", units="K")
                        def _Mass_2(field, data):
                                return data[(PartType_Gas_to_use, MassType_to_use)].in_units('Msun')
                        pf.add_field((PartType_Gas_to_use, "Mass_2"), function=_Mass_2, take_log=True, particle_type=False, display_name="Mass", units="Msun")
                        if draw_metal_map >= 1 or draw_metal_PDF == 1:
                                def _Metallicity_2(field, data):
                                        return data[(PartType_Gas_to_use, MetallicityType_to_use)]
                                pf.add_field((PartType_Gas_to_use, "Metallicity_2"), function=_Metallicity_2, take_log=True, particle_type=False, display_name="Metallicity", units="")
                        def _Density_2_minus_analytic(field, data):
                                return data[(PartType_Gas_to_use, "density")] - rho_agora_disk(data[(PartType_Gas_to_use, "particle_position_cylindrical_radius")], data[(PartType_Gas_to_use, "particle_position_cylindrical_z_abs")])
                        pf.add_field((PartType_Gas_to_use, "Density_2_minus_analytic"), function=_Density_2_minus_analytic, take_log=False, particle_type=False, display_name="Density Residual", units="g/cm**3")

                ####################################
                #      MAIN ANALYSIS ROUTINES      #
                ####################################

                # FIND CENTER AND PROJ_REGION
                v, cen = pf.h.find_max(("gas", "density")) # find the center to keep the galaxy at the center of all the images; here we assume that the gas disk is no bigger than 30 kpc in radius
                sp = pf.sphere(cen, (30.0, "kpc"))
                cen2 = sp.quantities.center_of_mass(use_gas=True, use_particles=False).in_units("kpc")
                sp2 = pf.sphere(cen2, (1.0, "kpc"))
                cen3 = sp2.quantities.max_location(("gas", "density"))
                center = pf.arr([cen3[1].d, cen3[2].d, cen3[3].d], 'code_length') # naive usage such as YTArray([cen3[1], cen3[2], cen3[3]]) doesn't work somehow for ART-II data
                if yt_version_pre_3_2_3 == 1:
                        center = pf.arr([cen3[2].d, cen3[3].d, cen3[4].d], 'code_length') # for yt-3.2.3 or before
                if codes[code] == "GASOLINE" and dataset_num == 2 and time == 1:
                        #center = pf.arr([2.7912903206399826e+20, 1.5205303149849894e+21, 1.5398968883245956e+21], 'cm') # a temporary hack for GASOLINE (center of the most massive stellar clump)
                        center = pf.arr([ 0.09045956,  0.49277032,  0.49904659], 'code_length')
                # if codes[code] == "GADGET-3" and dataset_num == 2 and time == 1:
                #       #center = pf.arr([6.184520935812296e+21, 4.972678132728082e+21, 6.559067311284074e+21], 'cm') # a temporary hack for GADGET-3 (center of the most massive stellar clump)
                #       center = pf.arr([ 2.00426673674,  1.61153523084,  2.12564895041], 'code_length')
                if codes[code] == "ART-I" or codes[code] == "ART-II" or codes[code] == "ENZO"  or codes[code] == "RAMSES":
                        proj_region = pf.box(center - YTArray([figure_width, figure_width, figure_width], 'kpc'),
                                             center + YTArray([figure_width, figure_width, figure_width], 'kpc')) # projected images made using a (2*figure_width)^3 box for AMR codes
                else:
                        proj_region = pf.all_data()

                # DENSITY MAPS
                if draw_density_map == 1:
                        my_cmap = copy.copy(matplotlib.cm.get_cmap('algae'))
                        my_cmap.set_bad(my_cmap(0)) # cmap range [0, 1)
                        my_cmap.set_under(my_cmap(0))
                        for ax in range(1, 3):
                                p = ProjectionPlot(pf, ax, ("gas", "density"), center = center, data_source=proj_region, width = (figure_width, 'kpc'), weight_field = None, fontsize=9)
                                p.set_zlim(("gas", "density"), 1e-5, 1e-1)
                                p.set_cmap(("gas", "density"), my_cmap)
                                plot = p.plots[("gas", "density")]

                                plot.figure = fig_density_map[time]
                                plot.axes = grid_density_map[time][(ax-1)*len(codes)+code].axes
                                if code == 0: plot.cax = grid_density_map[time].cbar_axes[0]
                                p._setup_plots()

                        if add_nametag == 1:
                                at = AnchoredText("%s" % codes[code], loc=2, prop=dict(size=6), frameon=True)
                                grid_density_map[time][code].axes.add_artist(at)

                # TEMPERATURE MAPS
                if draw_temperature_map == 1:
                        my_cmap = copy.copy(matplotlib.cm.get_cmap('algae'))
                        my_cmap.set_bad(my_cmap(0.6)) # (log(1e4) - log(1e1)) / (log(1e6) - log(1e1)) = 0.6
                        my_cmap.set_under(my_cmap(0))
                        for ax in range(1, 3):
                                p2 = ProjectionPlot(pf, ax, ("gas", "temperature"), center = center, data_source=proj_region, width = (figure_width, 'kpc'), weight_field = ("gas", "density_squared"), fontsize=9)
                                p2.set_zlim(("gas", "temperature"), 1e1, 1e6)
                                p2.set_cmap(("gas", "temperature"), my_cmap)
                                plot2 = p2.plots[("gas", "temperature")]

                                plot2.figure = fig_temperature_map[time]
                                plot2.axes = grid_temperature_map[time][(ax-1)*len(codes)+code].axes
                                if code == 0: plot2.cax = grid_temperature_map[time].cbar_axes[0]
                                p2._setup_plots()

                        if add_nametag == 1:
                                at = AnchoredText("%s" % codes[code], loc=2, prop=dict(size=6), frameon=True)
                                grid_temperature_map[time][code].axes.add_artist(at)

                # CELL-SIZE MAPS
                if draw_cellsize_map >= 1:
                        my_cmap = copy.copy(matplotlib.cm.get_cmap('algae'))
                        my_cmap.set_bad(my_cmap(0.99999)) # 1 doesn't work!
                        my_cmap.set_under(my_cmap(0))
                        if draw_cellsize_map == 1 or draw_cellsize_map == 3:
                                for ax in range(1, 3):
                                        p25 = ProjectionPlot(pf, ax, ("index", "cell_size"), center = center, data_source=proj_region, width = (figure_width, 'kpc'), weight_field = ("index", "cell_volume_inv2"), fontsize=9)
                                        p25.set_zlim(("index", "cell_size"), 10, 1e3)
                                        p25.set_cmap(("index", "cell_size"), my_cmap)
                                        plot25 = p25.plots[("index", "cell_size")]

                                        plot25.figure = fig_cellsize_map[time]
                                        plot25.axes = grid_cellsize_map[time][(ax-1)*len(codes)+code].axes
                                        if code == 0: plot25.cax = grid_cellsize_map[time].cbar_axes[0]
                                        p25._setup_plots()

                        # Create another map with a different resolution definition if requested
                        if draw_cellsize_map == 2 or draw_cellsize_map == 3:
                                for ax in range(1, 3):
                                        if codes[code] == "ART-I" or codes[code] == "ART-II" or codes[code] == "ENZO"  or codes[code] == "RAMSES":
                                                p251 = ProjectionPlot(pf, ax, ("index", "cell_size"), center = center, data_source=proj_region, width = (figure_width, 'kpc'), \
                                                                             weight_field = ("index", "cell_volume_inv2"), fontsize=9)
                                                p251.set_zlim(("index", "cell_size"), 10, 1e3)
                                                p251.set_cmap(("index", "cell_size"), my_cmap)
                                                plot251 = p251.plots[("index", "cell_size")]
                                        else:
                                                p251 = ProjectionPlot(pf, ax, ("gas", "particle_size"), center = center, data_source=proj_region, width = (figure_width, 'kpc'), \
                                                                             weight_field = ("gas", "particle_volume_inv2"), fontsize=9)
                                                p251.set_zlim(("gas", "particle_size"), 10, 1e3)
                                                p251.set_cmap(("gas", "particle_size"), my_cmap)
                                                plot251 = p251.plots[("gas", "particle_size")]

                                        plot251.figure = fig_cellsize_map_2[time]
                                        plot251.axes = grid_cellsize_map_2[time][(ax-1)*len(codes)+code].axes
                                        if code == 0: plot251.cax = grid_cellsize_map_2[time].cbar_axes[0]
                                        p251._setup_plots()

                        if add_nametag == 1:
                                if draw_cellsize_map == 1 or draw_cellsize_map == 3:
                                        at = AnchoredText("%s" % codes[code], loc=2, prop=dict(size=6), frameon=True)
                                        grid_cellsize_map[time][code].axes.add_artist(at)
                                if draw_cellsize_map == 2 or draw_cellsize_map == 3:
                                        at = AnchoredText("%s" % codes[code], loc=2, prop=dict(size=6), frameon=True)
                                        grid_cellsize_map_2[time][code].axes.add_artist(at)

                # ELEVATION MAPS
                if draw_elevation_map == 1:
                        def _CellzElevationpc(field,data):
                                return ((data[("index", "z")] - center[2]).in_units('pc'))
                        pf.add_field(("index", "z_elevation"), function=_CellzElevationpc, units='kpc', display_name="$z$ Elevation", take_log=False)
                        my_cmap = copy.copy(matplotlib.cm.get_cmap('algae'))
                        my_cmap.set_bad(my_cmap(0.5))
                        my_cmap.set_under(my_cmap(0))
                        for ax in range(2, 3):
                                p255 = ProjectionPlot(pf, ax, ("index", "z_elevation"), center = center, data_source=proj_region, width = (figure_width, 'kpc'), weight_field = ("gas", "density"), fontsize=9)
                                p255.set_zlim(("index", "z_elevation"), -1, 1)
                                p255.set_cmap(("index", "z_elevation"), my_cmap)
                                plot255 = p255.plots[("index", "z_elevation")]

                                plot255.figure = fig_elevation_map[time]
                                plot255.axes = grid_elevation_map[time][(ax-2)*len(codes)+code].axes
                                if code == 0: plot255.cax = grid_elevation_map[time].cbar_axes[0]
                                p255._setup_plots()

                        if add_nametag == 1:
                                at = AnchoredText("%s" % codes[code], loc=2, prop=dict(size=6), frameon=True)
                                grid_elevation_map[time][code].axes.add_artist(at)

                # Z-VELOCITY MAPS
                if draw_zvel_map == 1:
                        my_cmap = copy.copy(matplotlib.cm.get_cmap('algae'))
                        my_cmap.set_bad(my_cmap(0.5))
                        my_cmap.set_under(my_cmap(0))
                        for ax in range(1, 3):
                                p26 = ProjectionPlot(pf, ax, ("gas", "velocity_z"), center = center, data_source=proj_region, width = (figure_width, 'kpc'), weight_field = ("gas", "density_squared"), fontsize=9)
                                p26.set_zlim(("gas", "velocity_z"), -1e7, 1e7)
                                p26.set_cmap(("gas", "velocity_z"), my_cmap)
                                plot26 = p26.plots[("gas", "velocity_z")]

                                plot26.figure = fig_zvel_map[time]
                                plot26.axes = grid_zvel_map[time][(ax-1)*len(codes)+code].axes
                                if code == 0: plot26.cax = grid_zvel_map[time].cbar_axes[0]
                                p26._setup_plots()

                        if add_nametag == 1:
                                at = AnchoredText("%s" % codes[code], loc=2, prop=dict(size=6), frameon=True)
                                grid_zvel_map[time][code].axes.add_artist(at)

                # METAL MAPS
                if draw_metal_map >= 1:
                        my_cmap = copy.copy(matplotlib.cm.get_cmap('algae'))
                        my_cmap.set_bad(my_cmap(0.347)) # (0.02041 - 0.01) / (0.04 - 0.01) = 0.347
                        my_cmap.set_under(my_cmap(0))
                        for ax in range(1, 3):
                                pf.field_info[("gas", "metallicity")].take_log = False
                                p26 = ProjectionPlot(pf, ax, ("gas", "metallicity"), center = center, data_source=proj_region, width = (figure_width, 'kpc'), weight_field = ("gas", "density_squared"), fontsize=9)
                                p26.set_zlim(("gas", "metallicity"), 0.01, 0.04)
                                p26.set_cmap(("gas", "metallicity"), my_cmap)
                                plot26 = p26.plots[("gas", "metallicity")]

                                plot26.figure = fig_metal_map[time]
                                plot26.axes = grid_metal_map[time][(ax-1)*len(codes)+code].axes
                                if code == 0: plot26.cax = grid_metal_map[time].cbar_axes[0]
                                p26._setup_plots()

                        if add_nametag == 1:
                                if draw_metal_map == 2:
                                        sp = pf.sphere(center, (1, "Mpc"))
                                        total_metal_mass = sp.quantities.total_quantity([("gas", "metal_mass")])
                                        at = AnchoredText("%s - %e" % (codes[code], total_metal_mass), loc=2, prop=dict(size=6), frameon=True)
                                else:
                                        at = AnchoredText("%s" % codes[code], loc=2, prop=dict(size=6), frameon=True)
                                grid_metal_map[time][code].axes.add_artist(at)

                # STELLAR MAPS
                if draw_star_map == 1 and time != 0:
                        for ax in range(1, 3):
                                p27 = ParticleProjectionPlot(pf, ax, (PartType_Star_to_use, "particle_mass"), center = center, data_source=proj_region, width = (figure_width, 'kpc'), weight_field = None, fontsize=9)
                                p27.set_unit((PartType_Star_to_use, "particle_mass"), 'Msun')
                                p27.set_zlim((PartType_Star_to_use, "particle_mass"), 1e4, 1e7)
                                p27.set_cmap((PartType_Star_to_use, "particle_mass"), 'algae')
                                p27.set_buff_size(400) # default is 800
                                p27.set_colorbar_label((PartType_Star_to_use, "particle_mass"), "Stellar Mass Per Pixel ($\mathrm{M}_{\odot}$)")
                                plot27 = p27.plots[(PartType_Star_to_use, "particle_mass")]
                                # p27 = ParticleProjectionPlot(pf, ax, (PartType_Star_to_use, "particle_velocity_cylindrical_theta"), center = center, data_source=proj_region, width = (figure_width, 'kpc'), weight_field = (PartType_Star_to_use, "particle_mass"), fontsize=9)
                                # p27.set_unit((PartType_Star_to_use, "particle_velocity_cylindrical_theta"), 'km/s')
                                # p27.set_buff_size(400)
                                # p27.set_colorbar_label((PartType_Star_to_use, "particle_velocity_cylindrical_theta"), "Rotational Velocity (km/s)")
                                # plot27 = p27.plots[(PartType_Star_to_use, "particle_velocity_cylindrical_theta")]

                                plot27.figure = fig_star_map[time]
                                plot27.axes = grid_star_map[time][(ax-1)*len(codes)+code].axes
                                if code == 0: plot27.cax = grid_star_map[time].cbar_axes[0]
                                p27._setup_plots()

                        if add_nametag == 1:
                                at = AnchoredText("%s" % codes[code], loc=2, prop=dict(size=6), frameon=True)
                                grid_star_map[time][code].axes.add_artist(at)

                # STELLAR CLUMP STATISTICS AND ANNOTATED STELLAR MAPS
                if draw_star_clump_stats >= 1 and time != 0:
                        # For make yt's HaloCatalog to work with non-cosmological dataset, a fix needed to be applied to analysis_modules/halo_finding/halo_objects.py: self.period = ds.arr([3.254, 3.254, 3.254], 'Mpc')
                        pf.hubble_constant = 0.71; pf.omega_lambda = 0.73; pf.omega_matter = 0.27; pf.omega_curvature = 0.0; pf.current_redshift = 0.0 # another trick to make HaloCatalog work especially for ART-I dataset
                        if os.path.exists("./halo_catalogs/hop_%s_%05d/hop_%s_%05d.0.h5" % (codes[code], times[time], codes[code], times[time])) == False:
                                hc = HaloCatalog(data_ds=pf, finder_method='hop', output_dir="./halo_catalogs/hop_%s_%05d" % (codes[code], times[time]), \
                                                         finder_kwargs={'threshold': 2e8, 'dm_only': False, 'ptype': PartType_Star_to_use})
#                                                        finder_kwargs={'threshold': 2e5, 'dm_only': False, 'ptype': "all"})
                                hc.add_filter('quantity_value', 'particle_mass', '>', 2.6e6, 'Msun') # exclude halos with less than 30 particles
                                hc.add_filter('quantity_value', 'particle_mass', '<', 2.6e8, 'Msun') # exclude the most massive halo (threshold 1e8.4 is hand-picked, so one needs to be careful!)
                                hc.create()
                        halo_ds = load("./halo_catalogs/hop_%s_%05d/hop_%s_%05d.0.h5" % (codes[code], times[time], codes[code], times[time]))
                        hc = HaloCatalog(halos_ds=halo_ds, output_dir="./halo_catalogs/hop_%s_%05d" % (codes[code], times[time]))
                        hc.load()

                        halo_ad = hc.halos_ds.all_data()
                        star_clump_masses_hop[time].append(np.log10(halo_ad['particle_mass'][:].in_units("Msun")))

                        if os.path.exists("./halo_catalogs/fof_%s_%05d/fof_%s_%05d.0.h5" % (codes[code], times[time], codes[code], times[time])) == False:
                                hc2 = HaloCatalog(data_ds=pf, finder_method='fof', output_dir="./halo_catalogs/fof_%s_%05d" % (codes[code], times[time]), \
                                                          finder_kwargs={'link': 0.0025, 'dm_only': False, 'ptype': PartType_Star_to_use})
#                                                         finder_kwargs={'link': 0.02, 'dm_only': False, 'ptype': "all"})
                                hc2.add_filter('quantity_value', 'particle_mass', '>', 2.6e6, 'Msun') # exclude halos with less than 30 particles
                                hc2.add_filter('quantity_value', 'particle_mass', '<', 2.6e8, 'Msun') # exclude the most massive halo (threshold 1e8.4 is hand-picked, so one needs to be careful!)
                                hc2.create()
                        halo_ds2 = load("./halo_catalogs/fof_%s_%05d/fof_%s_%05d.0.h5" % (codes[code], times[time], codes[code], times[time]))
                        hc2 = HaloCatalog(halos_ds=halo_ds2, output_dir="./halo_catalogs/fof_%s_%05d" % (codes[code], times[time]))
                        hc2.load()

                        halo_ad2 = hc2.halos_ds.all_data()
                        star_clump_masses_fof[time].append(np.log10(halo_ad2['particle_mass'][:].in_units("Msun")))
                        # most_massive = halo_ad['particle_mass'].max().in_units('Msun')
                        # cen4 = [halo_ad['particle_position_x'][halo_ad['particle_mass'] == most_massive][0].in_units('code_length').d,
                        #       halo_ad['particle_position_y'][halo_ad['particle_mass'] == most_massive][0].in_units('code_length').d,
                        #       halo_ad['particle_position_z'][halo_ad['particle_mass'] == most_massive][0].in_units('code_length').d]


                        # Add additional star_map with annotated clumps if requested
                        if draw_star_clump_stats >= 2:
                                for ax in range(1, 3):
                                        p271 = ParticleProjectionPlot(pf, ax, (PartType_Star_to_use, "particle_mass"), center = center, data_source=proj_region, width = (figure_width, 'kpc'), weight_field = None, fontsize=9)
                                        p271.set_unit((PartType_Star_to_use, "particle_mass"), 'Msun')
                                        p271.set_zlim((PartType_Star_to_use, "particle_mass"), 1e4, 1e7)
                                        p271.set_cmap((PartType_Star_to_use, "particle_mass"), 'algae')
                                        p271.set_buff_size(400) # default is 800
                                        p271.set_colorbar_label((PartType_Star_to_use, "particle_mass"), "Stellar Mass Per Pixel ($\mathrm{M}_{\odot}$)")
                                        plot271 = p271.plots[(PartType_Star_to_use, "particle_mass")]
                                        if ax == 2:
                                                p271.annotate_halos(hc, factor=1.0, circle_args={'linewidth':0.7, 'alpha':0.8, 'facecolor':'none', 'edgecolor':'k'})

                                        plot271.figure = fig_star_map_2[time]
                                        plot271.axes = grid_star_map_2[time][(ax-1)*len(codes)+code].axes
                                        if code == 0: plot271.cax = grid_star_map_2[time].cbar_axes[0]
                                        p271._setup_plots()

                                if add_nametag == 1:
                                        at = AnchoredText("%s" % codes[code], loc=2, prop=dict(size=6), frameon=True)
                                        grid_star_map_2[time][code].axes.add_artist(at)

                                for ax in range(1, 3):
                                        p272 = ParticleProjectionPlot(pf, ax, (PartType_Star_to_use, "particle_mass"), center = center, data_source=proj_region, width = (figure_width, 'kpc'), weight_field = None, fontsize=9)
                                        p272.set_unit((PartType_Star_to_use, "particle_mass"), 'Msun')
                                        p272.set_zlim((PartType_Star_to_use, "particle_mass"), 1e4, 1e7)
                                        p272.set_cmap((PartType_Star_to_use, "particle_mass"), 'algae')
                                        p272.set_buff_size(400) # default is 800
                                        p272.set_colorbar_label((PartType_Star_to_use, "particle_mass"), "Stellar Mass Per Pixel ($\mathrm{M}_{\odot}$)")
                                        plot272 = p272.plots[(PartType_Star_to_use, "particle_mass")]
                                        # p272 = ParticleProjectionPlot(pf, ax, ("all", "particle_mass"), center = center, data_source=proj_region, width = (figure_width, 'kpc'), weight_field = None, fontsize=9)
                                        # p272.set_unit(("all", "particle_mass"), 'Msun')
                                        # p272.set_zlim(("all", "particle_mass"), 1e4, 1e9)
                                        # p272.set_cmap(("all", "particle_mass"), 'algae')
                                        # p272.set_buff_size(400) # default is 800
                                        # p272.set_colorbar_label(("all", "particle_mass"), "Mass Per Pixel ($\mathrm{M}_{\odot}$)")
                                        # plot272 = p272.plots[("all", "particle_mass")]
                                        if ax == 2:
                                                p272.annotate_halos(hc2, factor=1.0, circle_args={'linewidth':0.7, 'alpha':0.8, 'facecolor':'none', 'edgecolor':'k'})#, annotate_field='particle_mass')

                                        plot272.figure = fig_star_map_3[time]
                                        plot272.axes = grid_star_map_3[time][(ax-1)*len(codes)+code].axes
                                        if code == 0: plot272.cax = grid_star_map_3[time].cbar_axes[0]
                                        p272._setup_plots()

                                if add_nametag == 1:
                                        at = AnchoredText("%s" % codes[code], loc=2, prop=dict(size=6), frameon=True)
                                        grid_star_map_3[time][code].axes.add_artist(at)

                # SFR MAPS
                if draw_SFR_map >= 1 and time != 0:
                        sp = pf.sphere(center, (1.0*figure_width, "kpc"))
                        xy_bins = int(float(figure_width) / (float(aperture_size_SFR_map/1000.)))
                        my_cmap = copy.copy(matplotlib.cm.get_cmap('algae'))
                        my_cmap.set_bad(my_cmap(0))
                        my_cmap.set_under(my_cmap(0))
                        if codes[code] == "ART-I" or codes[code] == "ART-II" or codes[code] == "ENZO"  or codes[code] == "RAMSES":
                                def _position_x_from_center(field, data):
                                        return (data[("index", "x")] - center[0]).in_units('kpc')
                                pf.add_field(("index", "position_x_from_center"), function=_position_x_from_center, take_log=False, particle_type=False, units="kpc")
                                def _position_y_from_center(field, data):
                                        return (data[("index", "y")] - center[1]).in_units('kpc')
                                pf.add_field(("index", "position_y_from_center"), function=_position_y_from_center, take_log=False, particle_type=False, units="kpc")
                                def _cell_mass_surface_density_SFR_map(field, data):
                                        trans = np.zeros(data[("gas", "cell_mass")].shape)
                                        ind = np.where(data[("gas", "cell_mass")] > 0)
                                        trans[ind] = data[("gas", "cell_mass")][ind].in_units('Msun').d/(float(aperture_size_SFR_map))**2
                                        return data.ds.arr(trans, "Msun/pc**2").in_base(data.ds.unit_system.name)
                                pf.add_field(("gas", "cell_mass_surface_density_SFR_map"), function=_cell_mass_surface_density_SFR_map, \
                                                     display_name="$\mathrm{\Sigma}_\mathrm{gas}$", units='Msun/pc**2', particle_type=False, take_log=True)

                                p28 = PhasePlot(sp, ("index", "position_x_from_center"), ("index", "position_y_from_center"), \
                                                        ("gas", "cell_mass_surface_density_SFR_map"), weight_field=None, fontsize=9, x_bins=xy_bins, y_bins=xy_bins)
                                p28.set_zlim(("gas", "cell_mass_surface_density_SFR_map"), 1e0, 1e3)
                                p28.set_cmap(("gas", "cell_mass_surface_density_SFR_map"), my_cmap)
                                plot28 = p28.plots[("gas", "cell_mass_surface_density_SFR_map")]
                        else:
                                def _particle_position_x_from_center(field, data):
                                        return (data[(PartType_Gas_to_use, "particle_position_x")] - center[0]).in_units('kpc')
                                pf.add_field((PartType_Gas_to_use, "particle_position_x_from_center"), function=_particle_position_x_from_center, take_log=False, particle_type=True, units="kpc")
                                def _particle_position_y_from_center(field, data):
                                        return (data[(PartType_Gas_to_use, "particle_position_y")] - center[0]).in_units('kpc')
                                pf.add_field((PartType_Gas_to_use, "particle_position_y_from_center"), function=_particle_position_y_from_center, take_log=False, particle_type=True, units="kpc")
                                def _particle_mass_surface_density_SFR_map(field, data):
                                        trans = np.zeros(data[(PartType_Gas_to_use, "particle_mass")].shape)
                                        ind = np.where(data[(PartType_Gas_to_use, "particle_mass")] > 0)
                                        trans[ind] = data[(PartType_Gas_to_use, "particle_mass")][ind].in_units('Msun').d/(float(aperture_size_SFR_map))**2
                                        return data.ds.arr(trans, "Msun/pc**2").in_base(data.ds.unit_system.name)
                                pf.add_field((PartType_Gas_to_use, "particle_mass_surface_density_SFR_map"), function=_particle_mass_surface_density_SFR_map, \
                                                     display_name="$\mathrm{\Sigma}_\mathrm{gas}$", units='Msun/pc**2', particle_type=True, take_log=True)

                                p28 = ParticlePhasePlot(sp, (PartType_Gas_to_use, "particle_position_x_from_center"), (PartType_Gas_to_use, "particle_position_y_from_center"), \
                                                                (PartType_Gas_to_use, "particle_mass_surface_density_SFR_map"), weight_field=None, deposition='cic', fontsize=9, x_bins=xy_bins, y_bins=xy_bins)
                                p28.set_zlim((PartType_Gas_to_use, "particle_mass_surface_density_SFR_map"), 1e0, 1e3)
                                p28.set_cmap((PartType_Gas_to_use, "particle_mass_surface_density_SFR_map"), my_cmap)
                                plot28 = p28.plots[(PartType_Gas_to_use, "particle_mass_surface_density_SFR_map")]

                        p28.set_xlabel("x (kpc)")
                        p28.set_ylabel("y (kpc)")
                        p28.set_xlim(-15, 15)
                        p28.set_ylim(-15, 15)
                        plot28.figure = fig_degr_density_map[time]
                        plot28.axes = grid_degr_density_map[time][code].axes
                        if code == 0: plot28.cax = grid_degr_density_map[time].cbar_axes[0]
                        p28._setup_plots()

                        def _particle_position_x_from_center_2(field, data):
                                return (data[(PartType_StarBeforeFiltered_to_use, "particle_position_x")] - center[0]).in_units('kpc')
                        pf.add_field((PartType_StarBeforeFiltered_to_use, "particle_position_x_from_center"), function=_particle_position_x_from_center_2, take_log=False, particle_type=True, units="kpc")
                        pf.add_particle_filter(PartType_Star_to_use) # This is needed for a filtered particle type PartType_Star_to_use to work, because we have just created new particle fields.
                        def _particle_position_y_from_center_2(field, data):
                                return (data[(PartType_StarBeforeFiltered_to_use, "particle_position_y")] - center[0]).in_units('kpc')
                        pf.add_field((PartType_StarBeforeFiltered_to_use, "particle_position_y_from_center"), function=_particle_position_y_from_center_2, take_log=False, particle_type=True, units="kpc")
                        pf.add_particle_filter(PartType_Star_to_use)
                        def _particle_mass_young_stars_SFR_surface_density_SFR_map(field, data):
                                trans = np.zeros(data[(PartType_StarBeforeFiltered_to_use, "particle_mass")].shape)
                                ind = np.where(data[(PartType_StarBeforeFiltered_to_use, FormationTimeType_to_use)].in_units('Myr') > (pf.current_time.in_units('Myr').d - young_star_cutoff_SFR_map)) # mass for young stars only
                                trans[ind] = data[(PartType_StarBeforeFiltered_to_use, "particle_mass")][ind].in_units('Msun').d/(young_star_cutoff_SFR_map*1e6)/(float(aperture_size_SFR_map)/1000.)**2
                                return data.ds.arr(trans, "Msun/yr/kpc**2").in_base(data.ds.unit_system.name)
                        pf.add_field((PartType_StarBeforeFiltered_to_use, "particle_mass_young_stars_SFR_surface_density_SFR_map"), function=_particle_mass_young_stars_SFR_surface_density_SFR_map, \
                                             display_name="$\mathrm{\Sigma}_\mathrm{SFR}$", units='Msun/yr/kpc**2', particle_type=True, take_log=True)
                        pf.add_particle_filter(PartType_Star_to_use)

                        p281 = ParticlePhasePlot(sp, (PartType_Star_to_use, "particle_position_x_from_center"), (PartType_Star_to_use, "particle_position_y_from_center"), \
                                                         (PartType_Star_to_use, "particle_mass_young_stars_SFR_surface_density_SFR_map"), deposition='cic', weight_field=None, fontsize=9, x_bins=xy_bins, y_bins=xy_bins)
                        p281.set_zlim((PartType_Star_to_use, "particle_mass_young_stars_SFR_surface_density_SFR_map"), 3e-4, 3e-1)
                        p281.set_cmap((PartType_Star_to_use, "particle_mass_young_stars_SFR_surface_density_SFR_map"), my_cmap)
                        plot281 = p281.plots[(PartType_Star_to_use, "particle_mass_young_stars_SFR_surface_density_SFR_map")]

                        p281.set_xlabel("x (kpc)")
                        p281.set_ylabel("y (kpc)")
                        p281.set_xlim(-15, 15)
                        p281.set_ylim(-15, 15)
                        plot281.figure = fig_degr_sfr_map[time]
                        plot281.axes = grid_degr_sfr_map[time][code].axes
                        if code == 0: plot281.cax = grid_degr_sfr_map[time].cbar_axes[0]
                        p281._setup_plots()

                        if add_nametag == 1:
                                at = AnchoredText("%s" % codes[code], loc=2, prop=dict(size=6), frameon=True)
                                grid_degr_density_map[time][code].axes.add_artist(at)
                                at2 = AnchoredText("%s" % codes[code], loc=2, prop=dict(size=6), frameon=True)
                                grid_degr_sfr_map[time][code].axes.add_artist(at2)

                        # Record degraded maps for the later use in a local K-S plot
                        if draw_SFR_map == 2 and time != 0:
                                fn = p28.profile.save_as_dataset() # see http://yt-project.org/docs/dev/reference/api/generated/yt.data_objects.profiles.ProfileND.save_as_dataset.html
                                phaseplot_ds = load(fn)
                                if codes[code] == "ART-I" or codes[code] == "ART-II" or codes[code] == "ENZO"  or codes[code] == "RAMSES":
                                        phaseplot_surf_dens = phaseplot_ds.data[('data', "cell_mass_surface_density_SFR_map")]
                                else:
                                        phaseplot_surf_dens = phaseplot_ds.data[('data', "particle_mass_surface_density_SFR_map")]
                                fn_1 = p281.profile.save_as_dataset()
                                phaseplot_ds_1 = load(fn_1)
                                phaseplot_sfr_surf_dens = phaseplot_ds_1.data[('data', "particle_mass_young_stars_SFR_surface_density_SFR_map")]
                                xy_shape = phaseplot_surf_dens.shape[0]
                                surf_dens_SFR_map[time].append(phaseplot_surf_dens.reshape(1, xy_shape**2)[0])
                                sfr_surf_dens_SFR_map[time].append(phaseplot_sfr_surf_dens.reshape(1, xy_shape**2)[0])
                                call(["rm", "%s_%s.h5" % (str(pf), 'ParticleProfile')]) # Remove unwanted .h5 files saved
                                if codes[code] == "ART-I" or codes[code] == "ART-II" or codes[code] == "ENZO"  or codes[code] == "RAMSES":
                                        call(["rm", "%s_%s.h5" % (str(pf), 'Profile2D')])

                # DENSITY-TEMPERATURE PDF
                if draw_PDF >= 1:
                        sp = pf.sphere(center, (0.5*figure_width, "kpc"))
                        my_cmap2 = copy.copy(matplotlib.cm.get_cmap('algae'))
                        my_cmap2.set_under(my_cmap2(0))
                        if codes[code] == "ART-I" or codes[code] == "ART-II" or codes[code] == "ENZO"  or codes[code] == "RAMSES":
                                p3 = PhasePlot(sp, ("gas", "density"), ("gas", "temperature"), ("gas", "cell_mass"), weight_field=None, fontsize=12, x_bins=300, y_bins=300)
                                p3.set_unit("cell_mass", 'Msun')
                                p3.set_zlim(("gas", "cell_mass"), 1e3, 1e8)
                                p3.set_cmap(("gas", "cell_mass"), my_cmap2)
                                p3.set_colorbar_label(("gas", "cell_mass"), "Mass ($\mathrm{M}_{\odot}$)")
                                plot3 = p3.plots[("gas", "cell_mass")]
                        else:
                                # Because ParticlePhasePlot doesn't yet work for a log-log PDF for some reason, I will do the following trick.
                                p3 = PhasePlot(sp, (PartType_Gas_to_use, "Density_2"), (PartType_Gas_to_use, "Temperature_2"), (PartType_Gas_to_use, "Mass_2"), weight_field=None, fontsize=12, x_bins=300, y_bins=300)
                                p3.set_zlim((PartType_Gas_to_use, "Mass_2"), 1e3, 1e8)
                                p3.set_cmap((PartType_Gas_to_use, "Mass_2"), my_cmap2)
                                plot3 = p3.plots[(PartType_Gas_to_use, "Mass_2")]

                        p3.set_xlim(1e-29, 1e-21)
                        p3.set_ylim(10, 1e7)

                        plot3.figure = fig_PDF[time]
                        plot3.axes = grid_PDF[time][code].axes
                        if code == 0: plot3.cax = grid_PDF[time].cbar_axes[0]
                        p3._setup_plots()

                        # Add constant pressure and constant entropy lines to guide the eyes
                        if draw_PDF >= 2:
                                dummy_density = np.logspace(-29, -21, num=4)
                                for constant_i in range(1, 100):
                                        constant = 1e-100*100**constant_i
                                        line = ln.Line2D(dummy_density, constant/dummy_density, linestyle=":", linewidth=0.6, color='k', alpha=0.7) # constant pressure P ~ n*T
                                        grid_PDF[time][code].axes.add_line(line)
                                        line2 = ln.Line2D(dummy_density, constant*dummy_density**(2./3.), linestyle="-.", linewidth=0.6, color='k', alpha=0.7) # constant entropy S ~ n**(-2/3)*T
                                        grid_PDF[time][code].axes.add_line(line2)

                        if add_nametag == 1:
                                at = AnchoredText("%s" % codes[code], loc=3, prop=dict(size=10), frameon=True)
                                grid_PDF[time][code].axes.add_artist(at)

                # POSITION-VELOCITY PDF FOR GAS
                if draw_pos_vel_PDF >= 1:
                        sp = pf.sphere(center, (0.5*figure_width, "kpc"))
                        sp.set_field_parameter("normal", disk_normal_vector)
                        if codes[code] == "ART-I" or codes[code] == "ART-II" or codes[code] == "ENZO"  or codes[code] == "RAMSES":
                                sp_dense = sp.cut_region(["obj['gas', 'density'].in_units('g/cm**3') > 1.e-25"]) # For AMR codes, consider only the cells that are dense enough
                                pf.field_info[("index", "cylindrical_r")].take_log = False
                                pf.field_info[("gas", "cylindrical_tangential_velocity")].take_log = False
                                p4 = PhasePlot(sp_dense, ("index", "cylindrical_r"), ("gas", "cylindrical_tangential_velocity"), ("gas", "cell_mass"), weight_field=None, fontsize=12, x_bins=300, y_bins=300)
                                p4.set_unit("cylindrical_r", 'kpc')
                                p4.set_unit("cylindrical_tangential_velocity", 'km/s')
                                p4.set_unit("cell_mass", 'Msun')
                                p4.set_zlim(("gas", "cell_mass"), 1e3, 1e8)
                                p4.set_cmap(("gas", "cell_mass"), 'algae')
                                p4.set_colorbar_label(("gas", "cell_mass"), "Mass ($\mathrm{M}_{\odot}$)")
                                plot4 = p4.plots[("gas", "cell_mass")]
                        else:
                                pf.field_info[(PartType_Gas_to_use, "particle_position_cylindrical_radius")].take_log = False
                                pf.field_info[(PartType_Gas_to_use, "particle_velocity_cylindrical_theta")].take_log = False
                                p4 = ParticlePhasePlot(sp, (PartType_Gas_to_use, "particle_position_cylindrical_radius"), (PartType_Gas_to_use, "particle_velocity_cylindrical_theta"), \
                                                               (PartType_Gas_to_use, MassType_to_use), deposition='ngp', weight_field=None, fontsize=12, x_bins=300, y_bins=300)
                                p4.set_unit("particle_position_cylindrical_radius", 'kpc')
                                p4.set_unit("particle_velocity_cylindrical_theta", 'km/s')
                                p4.set_unit(MassType_to_use, 'Msun')
                                p4.set_zlim((PartType_Gas_to_use, MassType_to_use), 1e3, 1e8)
                                p4.set_cmap((PartType_Gas_to_use, MassType_to_use), 'algae')
                                p4.set_colorbar_label((PartType_Gas_to_use, MassType_to_use), "Mass ($\mathrm{M}_{\odot}$)")
                                plot4 = p4.plots[(PartType_Gas_to_use, MassType_to_use)]

                        p4.set_xlabel("Cylindrical Radius (kpc)")
                        p4.set_ylabel("Rotational Velocity (km/s)")
                        p4.set_xlim(0, 14)
                        p4.set_ylim(-100, 350)

                        plot4.figure = fig_pos_vel_PDF[time]
                        plot4.axes = grid_pos_vel_PDF[time][code].axes
                        if code == 0: plot4.cax = grid_pos_vel_PDF[time].cbar_axes[0]
                        p4._setup_plots()

                        # Add 1D profile line if requested
                        if draw_pos_vel_PDF >= 2 and time != 0:
                                if codes[code] == "ART-I" or codes[code] == "ART-II" or codes[code] == "ENZO"  or codes[code] == "RAMSES":
                                        p5 = ProfilePlot(sp_dense, ("index", "cylindrical_r"),  ("gas", "cylindrical_tangential_velocity"), \
                                                                 weight_field=("gas", "cell_mass"), n_bins=50, x_log=False)
                                        p5.set_log("cylindrical_tangential_velocity", False)
                                        p5.set_log("cylindrical_r", False)
                                        p5.set_unit("cylindrical_r", 'kpc')
                                        p5.set_xlim(1e-3, 14)
                                        p5.set_ylim("cylindrical_tangential_velocity", -100, 350)
                                        line = ln.Line2D(p5.profiles[0].x.in_units('kpc'), p5.profiles[0]["cylindrical_tangential_velocity"].in_units('km/s'), linestyle="-", linewidth=2, color='k', alpha=0.7)
                                        pos_vel_xs[time].append(p5.profiles[0].x.in_units('kpc').d)
                                        pos_vel_profiles[time].append(p5.profiles[0]["cylindrical_tangential_velocity"].in_units('km/s').d)
                                else:
                                        p5 = ProfilePlot(sp, (PartType_Gas_to_use, "particle_position_cylindrical_radius"), (PartType_Gas_to_use, "particle_velocity_cylindrical_theta"), \
                                                                 weight_field=(PartType_Gas_to_use, MassType_to_use), n_bins=50, x_log=False)
                                        p5.set_log("particle_velocity_cylindrical_theta", False)
                                        p5.set_log("particle_position_cylindrical_radius", False)
                                        p5.set_unit("particle_position_cylindrical_radius", 'kpc')
                                        p5.set_xlim(1e-3, 14)
                                        p5.set_ylim("particle_velocity_cylindrical_theta", -100, 350)
                                        line = ln.Line2D(p5.profiles[0].x.in_units('kpc'), p5.profiles[0]["particle_velocity_cylindrical_theta"].in_units('km/s'), linestyle="-", linewidth=2, color='k', alpha=0.7)
                                        pos_vel_xs[time].append(p5.profiles[0].x.in_units('kpc').d)
                                        pos_vel_profiles[time].append(p5.profiles[0]["particle_velocity_cylindrical_theta"].in_units('km/s').d)
                                grid_pos_vel_PDF[time][code].axes.add_line(line)

                        # Add velocity dispersion profile if requested
                        if draw_pos_vel_PDF >= 3 and time != 0:
                                if codes[code] == "ART-I" or codes[code] == "ART-II" or codes[code] == "ENZO"  or codes[code] == "RAMSES":
                                        # To calculate velocity dispersion -- residual velocity components other than the rotational velocity found above -- we will need a mass-weighted average of [v_i - v_rot_local]**2
                                        # Here we assumed that the disk is on x-y plane; i.e. the variable disk_normal_vector above needs to be [0.0, 0.0, 1.0]
                                        def _local_rotational_velocity_x(field, data):
                                                trans = np.zeros(data[("gas", "velocity_x")].shape)
                                                dr = 0.5*(pos_vel_xs[time][code][1] - pos_vel_xs[time][code][0])
                                                for radius, v_rot_local in zip(pos_vel_xs[time][code], pos_vel_profiles[time][code]):
                                                        ind = np.where((data[("index", "cylindrical_r")].in_units("kpc") >= (radius - dr)) & (data[("index", "cylindrical_r")].in_units("kpc") < (radius + dr)))
                                                        trans[ind] = -np.sin(data["index", 'cylindrical_theta'][ind]) * v_rot_local * 1e5 # in cm/s
                                                return data.ds.arr(trans, "cm/s").in_base(data.ds.unit_system.name)
                                        pf.add_field(("gas", "local_rotational_velocity_x"), function=_local_rotational_velocity_x, take_log=False, particle_type=False, units="cm/s")
                                        def _local_rotational_velocity_y(field, data):
                                                trans = np.zeros(data[("gas", "velocity_y")].shape)
                                                dr = 0.5*(pos_vel_xs[time][code][1] - pos_vel_xs[time][code][0])
                                                for radius, v_rot_local in zip(pos_vel_xs[time][code], pos_vel_profiles[time][code]):
                                                        ind = np.where((data[("index", "cylindrical_r")].in_units("kpc") >= (radius - dr)) & (data[("index", "cylindrical_r")].in_units("kpc") < (radius + dr)))
                                                        trans[ind] =  np.cos(data["index", 'cylindrical_theta'][ind]) * v_rot_local * 1e5
                                                return data.ds.arr(trans, "cm/s").in_base(data.ds.unit_system.name)
                                        pf.add_field(("gas", "local_rotational_velocity_y"), function=_local_rotational_velocity_y, take_log=False, particle_type=False, units="cm/s")
                                        def _velocity_minus_local_rotational_velocity_squared(field, data):
                                                return (data[("gas", "velocity_x")] - data[("gas", "local_rotational_velocity_x")])**2 + \
                                                    (data[("gas", "velocity_y")] - data[("gas", "local_rotational_velocity_y")])**2 + \
                                                    (data[("gas", "velocity_z")])**2
                                        pf.add_field(("gas", "velocity_minus_local_rotational_velocity_squared"), function=_velocity_minus_local_rotational_velocity_squared,
                                                     take_log=False, particle_type=False, units="cm**2/s**2")
#					slc = SlicePlot(pf, 'z', ("gas", "local_rotational_velocity_y"), center = center, width = (figure_width, 'kpc'))
#					slc = SlicePlot(pf, 'z', ("gas", "velocity_minus_local_rotational_velocity_squared"), center = center, width = (figure_width, 'kpc')) # this should give zeros everywhere for IC at t=0
#					slc.save()

                                        p55 = ProfilePlot(sp_dense, ("index", "cylindrical_r"),  ("gas", "velocity_minus_local_rotational_velocity_squared"), \
                                                                  weight_field=("gas", "cell_mass"), n_bins=50, x_log=False)
                                        p55.set_log("cylindrical_r", False)
                                        p55.set_unit("cylindrical_r", 'kpc')
                                        p55.set_xlim(1e-3, 14)
                                        pos_disp_xs[time].append(p55.profiles[0].x.in_units('kpc').d)
                                        pos_disp_profiles[time].append(np.sqrt(p55.profiles[0]["velocity_minus_local_rotational_velocity_squared"]).in_units('km/s').d)
                                else:
                                        def _particle_local_rotational_velocity_x(field, data):
                                                trans = np.zeros(data[(PartType_Gas_to_use, "particle_velocity_x")].shape)
                                                dr = 0.5*(pos_vel_xs[time][code][1] - pos_vel_xs[time][code][0])
                                                for radius, v_rot_local in zip(pos_vel_xs[time][code], pos_vel_profiles[time][code]):
                                                        ind = np.where((data[(PartType_Gas_to_use, "particle_position_cylindrical_radius")].in_units("kpc") >= (radius - dr)) & \
                                                                               (data[(PartType_Gas_to_use, "particle_position_cylindrical_radius")].in_units("kpc") < (radius + dr)))
                                                        trans[ind] = -np.sin(data[(PartType_Gas_to_use, "particle_position_cylindrical_theta")][ind]) * v_rot_local * 1e5
                                                return data.ds.arr(trans, "cm/s").in_base(data.ds.unit_system.name)
                                        pf.add_field((PartType_Gas_to_use, "particle_local_rotational_velocity_x"), function=_particle_local_rotational_velocity_x, take_log=False, particle_type=True, units="cm/s")
                                        def _particle_local_rotational_velocity_y(field, data):
                                                trans = np.zeros(data[(PartType_Gas_to_use, "particle_velocity_y")].shape)
                                                dr = 0.5*(pos_vel_xs[time][code][1] - pos_vel_xs[time][code][0])
                                                for radius, v_rot_local in zip(pos_vel_xs[time][code], pos_vel_profiles[time][code]):
                                                        ind = np.where((data[(PartType_Gas_to_use, "particle_position_cylindrical_radius")].in_units("kpc") >= (radius - dr)) & \
                                                                               (data[(PartType_Gas_to_use, "particle_position_cylindrical_radius")].in_units("kpc") < (radius + dr)))
                                                        trans[ind] = np.cos(data[(PartType_Gas_to_use, "particle_position_cylindrical_theta")][ind]) * v_rot_local * 1e5
                                                return data.ds.arr(trans, "cm/s").in_base(data.ds.unit_system.name)
                                        pf.add_field((PartType_Gas_to_use, "particle_local_rotational_velocity_y"), function=_particle_local_rotational_velocity_y, take_log=False, particle_type=True, units="cm/s")
                                        def _particle_velocity_minus_local_rotational_velocity_squared(field, data):
                                                return (data[(PartType_Gas_to_use, "particle_velocity_x")] - data[(PartType_Gas_to_use, "particle_local_rotational_velocity_x")])**2 + \
                                                    (data[(PartType_Gas_to_use, "particle_velocity_y")] - data[(PartType_Gas_to_use, "particle_local_rotational_velocity_y")])**2 + \
                                                    (data[(PartType_Gas_to_use, "particle_velocity_z")])**2
                                        pf.add_field((PartType_Gas_to_use, "particle_velocity_minus_local_rotational_velocity_squared"), function=_particle_velocity_minus_local_rotational_velocity_squared,
                                                     take_log=False, particle_type=True, units="cm**2/s**2")

                                        p55 = ProfilePlot(sp, (PartType_Gas_to_use, "particle_position_cylindrical_radius"), (PartType_Gas_to_use, "particle_velocity_minus_local_rotational_velocity_squared"), \
                                                                 weight_field=(PartType_Gas_to_use, MassType_to_use), n_bins=50, x_log=False)
                                        p55.set_log("particle_position_cylindrical_radius", False)
                                        p55.set_unit("particle_position_cylindrical_radius", 'kpc')
                                        p55.set_xlim(1e-3, 14)
                                        pos_disp_xs[time].append(p55.profiles[0].x.in_units('kpc').d)
                                        pos_disp_profiles[time].append(np.sqrt(p55.profiles[0]["particle_velocity_minus_local_rotational_velocity_squared"]).in_units('km/s').d)

                                # Add vertical velocity dispersion profile if requested
                                if codes[code] == "ART-I" or codes[code] == "ART-II" or codes[code] == "ENZO"  or codes[code] == "RAMSES":
                                        def _velocity_z_squared(field, data):
                                                return (data[("gas", "velocity_z")])**2
                                        pf.add_field(("gas", "velocity_z_squared"), function=_velocity_z_squared, take_log=False, particle_type=False, units="cm**2/s**2")

                                        p551 = ProfilePlot(sp_dense, ("index", "cylindrical_r"),  ("gas", "velocity_z_squared"), weight_field=("gas", "cell_mass"), n_bins=50, x_log=False)
                                        p551.set_log("cylindrical_r", False)
                                        p551.set_unit("cylindrical_r", 'kpc')
                                        p551.set_xlim(1e-3, 14)
                                        pos_disp_vert_xs[time].append(p551.profiles[0].x.in_units('kpc').d)
                                        pos_disp_vert_profiles[time].append(np.sqrt(p551.profiles[0]["velocity_z_squared"]).in_units('km/s').d)
                                else:
                                        def _particle_velocity_z_squared(field, data):
                                                return (data[(PartType_Gas_to_use, "particle_velocity_z")])**2
                                        pf.add_field((PartType_Gas_to_use, "particle_velocity_z_squared"), function=_particle_velocity_z_squared, take_log=False, particle_type=True, units="cm**2/s**2")

                                        p551 = ProfilePlot(sp, (PartType_Gas_to_use, "particle_position_cylindrical_radius"), (PartType_Gas_to_use, "particle_velocity_z_squared"), \
                                                                 weight_field=(PartType_Gas_to_use, MassType_to_use), n_bins=50, x_log=False)
                                        p551.set_log("particle_position_cylindrical_radius", False)
                                        p551.set_unit("particle_position_cylindrical_radius", 'kpc')
                                        p551.set_xlim(1e-3, 14)
                                        pos_disp_vert_xs[time].append(p551.profiles[0].x.in_units('kpc').d)
                                        pos_disp_vert_profiles[time].append(np.sqrt(p551.profiles[0]["particle_velocity_z_squared"]).in_units('km/s').d)

                        if add_nametag == 1:
                                at = AnchoredText("%s" % codes[code], loc=4, prop=dict(size=10), frameon=True)
                                grid_pos_vel_PDF[time][code].axes.add_artist(at)

                # POSITION-VELOCITY PDF FOR NEW STARS
                if draw_star_pos_vel_PDF >= 1 and time != 0:
                        sp = pf.sphere(center, (0.5*figure_width, "kpc"))
                        sp.set_field_parameter("normal", disk_normal_vector)
                        pf.field_info[(PartType_Star_to_use, "particle_position_cylindrical_radius")].take_log = False
                        pf.field_info[(PartType_Star_to_use, "particle_velocity_cylindrical_theta")].take_log = False
                        pf.field_info[(PartType_Star_to_use, "particle_mass")].output_units = 'code_mass' # this turned out to be crucial!; otherwise wrong output_unit 'g' is assumed in ParticlePhasePlot->create_profile in visualizaiton/particle_plots.py for ART-I/ENZO/RAMSES
                        p41 = ParticlePhasePlot(sp, (PartType_Star_to_use, "particle_position_cylindrical_radius"), (PartType_Star_to_use, "particle_velocity_cylindrical_theta"), \
                                                       (PartType_Star_to_use, "particle_mass"), deposition='ngp', weight_field=None, fontsize=12, x_bins=300, y_bins=300)
                        p41.set_unit("particle_position_cylindrical_radius", 'kpc')
                        p41.set_unit("particle_velocity_cylindrical_theta", 'km/s')
                        p41.set_unit("particle_mass", 'Msun') # requires a change in set_unit in visualization/profile_plotter.py: remove self.plots[field].zmin, self.plots[field].zmax = (None, None)
#			p41.set_unit((PartType_Star_to_use, "particle_mass"), 'Msun') # Neither this nor above works without such change
                        p41.set_zlim((PartType_Star_to_use, "particle_mass"), 1e3, 1e7)
                        p41.set_cmap((PartType_Star_to_use, "particle_mass"), 'algae')

                        p41.set_colorbar_label((PartType_Star_to_use, "particle_mass"), "Newly Formed Stellar Mass ($\mathrm{M}_{\odot}$)")
                        plot41 = p41.plots[(PartType_Star_to_use, "particle_mass")]

                        p41.set_xlabel("Cylindrical Radius (kpc)")
                        p41.set_ylabel("Rotational Velocity (km/s)")
                        p41.set_xlim(0, 14)
                        p41.set_ylim(-100, 350)

                        plot41.figure = fig_star_pos_vel_PDF[time]
                        plot41.axes = grid_star_pos_vel_PDF[time][code].axes
                        if code == 0: plot41.cax = grid_star_pos_vel_PDF[time].cbar_axes[0]
                        p41._setup_plots()

                        # Add 1D profile line if requested
                        if draw_star_pos_vel_PDF >= 2 and time != 0:
                                p51 = ProfilePlot(sp, (PartType_Star_to_use, "particle_position_cylindrical_radius"), (PartType_Star_to_use, "particle_velocity_cylindrical_theta"), \
                                                         weight_field=(PartType_Star_to_use, "particle_mass"), n_bins=50, x_log=False)
                                p51.set_log((PartType_Star_to_use, "particle_velocity_cylindrical_theta"), False)
                                p51.set_log((PartType_Star_to_use, "particle_position_cylindrical_radius"), False)
                                p51.set_unit("particle_position_cylindrical_radius", 'kpc')
                                p51.set_xlim(1e-3, 14)
                                p51.set_ylim("particle_velocity_cylindrical_theta", -100, 350)
                                line = ln.Line2D(p51.profiles[0].x.in_units('kpc'), p51.profiles[0]["particle_velocity_cylindrical_theta"].in_units('km/s'), linestyle="-", linewidth=2, color='k', alpha=0.7)
                                star_pos_vel_xs[time].append(p51.profiles[0].x.in_units('kpc').d)
                                star_pos_vel_profiles[time].append(p51.profiles[0]["particle_velocity_cylindrical_theta"].in_units('km/s').d)
                                grid_star_pos_vel_PDF[time][code].axes.add_line(line)

                        # Add velocity dispersion profile if requested
                        if draw_star_pos_vel_PDF >= 3 and time != 0:
                                def _particle_local_rotational_velocity_x(field, data):
                                        trans = np.zeros(data[(PartType_StarBeforeFiltered_to_use, "particle_velocity_x")].shape)
                                        dr = 0.5*(star_pos_vel_xs[time][code][1] - star_pos_vel_xs[time][code][0])
                                        for radius, v_rot_local in zip(star_pos_vel_xs[time][code], star_pos_vel_profiles[time][code]):
                                                ind = np.where((data[(PartType_StarBeforeFiltered_to_use, "particle_position_cylindrical_radius")].in_units("kpc") >= (radius - dr)) & \
                                                                       (data[(PartType_StarBeforeFiltered_to_use, "particle_position_cylindrical_radius")].in_units("kpc") < (radius + dr)))
                                                trans[ind] = -np.sin(data[(PartType_StarBeforeFiltered_to_use, "particle_position_cylindrical_theta")][ind]) * v_rot_local * 1e5
                                        return data.ds.arr(trans, "cm/s").in_base(data.ds.unit_system.name)
                                pf.add_field((PartType_StarBeforeFiltered_to_use, "particle_local_rotational_velocity_x"), function=_particle_local_rotational_velocity_x, take_log=False, particle_type=True, units="cm/s")
                                def _particle_local_rotational_velocity_y(field, data):
                                        trans = np.zeros(data[(PartType_StarBeforeFiltered_to_use, "particle_velocity_y")].shape)
                                        dr = 0.5*(star_pos_vel_xs[time][code][1] - star_pos_vel_xs[time][code][0])
                                        for radius, v_rot_local in zip(star_pos_vel_xs[time][code], star_pos_vel_profiles[time][code]):
                                                ind = np.where((data[(PartType_StarBeforeFiltered_to_use, "particle_position_cylindrical_radius")].in_units("kpc") >= (radius - dr)) & \
                                                                       (data[(PartType_StarBeforeFiltered_to_use, "particle_position_cylindrical_radius")].in_units("kpc") < (radius + dr)))
                                                trans[ind] = np.cos(data[(PartType_StarBeforeFiltered_to_use, "particle_position_cylindrical_theta")][ind]) * v_rot_local * 1e5
                                        return data.ds.arr(trans, "cm/s").in_base(data.ds.unit_system.name)
                                pf.add_field((PartType_StarBeforeFiltered_to_use, "particle_local_rotational_velocity_y"), function=_particle_local_rotational_velocity_y, take_log=False, particle_type=True, units="cm/s")
                                def _particle_velocity_minus_local_rotational_velocity_squared(field, data):
                                        return (data[(PartType_StarBeforeFiltered_to_use, "particle_velocity_x")] - data[(PartType_StarBeforeFiltered_to_use, "particle_local_rotational_velocity_x")])**2 + \
                                            (data[(PartType_StarBeforeFiltered_to_use, "particle_velocity_y")] - data[(PartType_StarBeforeFiltered_to_use, "particle_local_rotational_velocity_y")])**2 + \
                                            (data[(PartType_StarBeforeFiltered_to_use, "particle_velocity_z")])**2
                                pf.add_field((PartType_StarBeforeFiltered_to_use, "particle_velocity_minus_local_rotational_velocity_squared"), function=_particle_velocity_minus_local_rotational_velocity_squared,
                                             take_log=False, particle_type=True, units="cm**2/s**2")
                                pf.add_particle_filter(PartType_Star_to_use) # This is needed for a filtered particle type PartType_Star_to_use to work, because we have just created new particle fields.

                                p515 = ProfilePlot(sp, (PartType_Star_to_use, "particle_position_cylindrical_radius"), (PartType_Star_to_use, "particle_velocity_minus_local_rotational_velocity_squared"), \
                                                          weight_field=(PartType_Star_to_use, "particle_mass"), n_bins=50, x_log=False)
                                p515.set_log((PartType_Star_to_use, "particle_position_cylindrical_radius"), False)
                                p515.set_unit("particle_position_cylindrical_radius", 'kpc')
                                p515.set_xlim(1e-3, 14)
                                star_pos_disp_xs[time].append(p515.profiles[0].x.in_units('kpc').d)
                                star_pos_disp_profiles[time].append(np.sqrt(p515.profiles[0]["particle_velocity_minus_local_rotational_velocity_squared"]).in_units('km/s').d)

                        # Add vertical velocity dispersion profile if requested
                        if draw_star_pos_vel_PDF >= 4 and time != 0:
                                def _particle_velocity_z_squared(field, data):
                                        return (data[(PartType_StarBeforeFiltered_to_use, "particle_velocity_z")])**2
                                pf.add_field((PartType_StarBeforeFiltered_to_use, "particle_velocity_z_squared"), function=_particle_velocity_z_squared,
                                             take_log=False, particle_type=True, units="cm**2/s**2")
                                pf.add_particle_filter(PartType_Star_to_use) # This is needed for a filtered particle type PartType_Star_to_use to work, because we have just created new particle fields.

                                p5151 = ProfilePlot(sp, (PartType_Star_to_use, "particle_position_cylindrical_radius"), (PartType_Star_to_use, "particle_velocity_z_squared"), \
                                                          weight_field=(PartType_Star_to_use, "particle_mass"), n_bins=50, x_log=False)
                                p5151.set_log((PartType_Star_to_use, "particle_position_cylindrical_radius"), False)
                                p5151.set_unit("particle_position_cylindrical_radius", 'kpc')
                                p5151.set_xlim(1e-3, 14)
                                star_pos_disp_vert_xs[time].append(p5151.profiles[0].x.in_units('kpc').d)
                                star_pos_disp_vert_profiles[time].append(np.sqrt(p5151.profiles[0]["particle_velocity_z_squared"]).in_units('km/s').d)

                        if add_nametag == 1:
                                at = AnchoredText("%s" % codes[code], loc=4, prop=dict(size=10), frameon=True)
                                grid_star_pos_vel_PDF[time][code].axes.add_artist(at)

                # RADIUS-HEIGHT PDF
                if draw_rad_height_PDF >= 1:
                        sp = pf.sphere(center, (0.5*figure_width, "kpc"))
                        if codes[code] == "ART-I" or codes[code] == "ART-II" or codes[code] == "ENZO"  or codes[code] == "RAMSES":
                                if draw_rad_height_PDF == 1 or draw_rad_height_PDF == 2:
                                        p55 = PhasePlot(sp, ("index", "cylindrical_r"), ("index", "cylindrical_z_abs"), ("gas", "density"), weight_field=("gas", "cell_mass"), fontsize=12, x_bins=200, y_bins=200)
                                        p55.set_zlim(("gas", "density"), 1e-26, 1e-21)
                                        p55.set_cmap(("gas", "density"), 'algae')
                                        p55.set_log("cylindrical_r", False)
                                        p55.set_log("cylindrical_z_abs", False)
                                        p55.set_unit("cylindrical_r", 'kpc')
                                        p55.set_unit("cylindrical_z_abs", 'kpc')
                                        plot55 = p55.plots[("gas", "density")]
                                elif draw_rad_height_PDF == 3:
                                        p55 = PhasePlot(sp, ("index", "cylindrical_r"), ("index", "cylindrical_z_abs"), ("gas", "density_minus_analytic"), weight_field=("gas", "cell_mass"), fontsize=12, x_bins=200, y_bins=200)
                                        p55.set_zlim(("gas", "density_minus_analytic"), -1e-24, 1e-24)
                                        p55.set_cmap(("gas", "density_minus_analytic"), 'algae')
                                        p55.set_log("cylindrical_r", False)
                                        p55.set_log("cylindrical_z_abs", False)
                                        p55.set_unit("cylindrical_r", 'kpc')
                                        p55.set_unit("cylindrical_z_abs", 'kpc')
                                        plot55 = p55.plots[("gas", "density_minus_analytic")]
                        else:
                                # Because ParticlePhasePlot doesn't work for these newly created fields for some reason, I will do the following trick.
                                if draw_rad_height_PDF == 1 or draw_rad_height_PDF == 2:
                                        p55 = PhasePlot(sp, (PartType_Gas_to_use, "particle_position_cylindrical_radius"), (PartType_Gas_to_use, "particle_position_cylindrical_z_abs"), (PartType_Gas_to_use, "Density_2"), weight_field=(PartType_Gas_to_use, "Mass_2"), fontsize=12, x_bins=200, y_bins=200)
                                        p55.set_zlim((PartType_Gas_to_use, "Density_2"), 1e-26, 1e-21)
                                        p55.set_cmap((PartType_Gas_to_use, "Density_2"), 'algae')
                                        p55.set_log("particle_position_cylindrical_radius", False)
                                        p55.set_log("particle_position_cylindrical_z_abs", False)
                                        p55.set_unit("particle_position_cylindrical_radius", 'kpc')
                                        p55.set_unit("particle_position_cylindrical_z_abs", 'kpc')
                                        plot55 = p55.plots[(PartType_Gas_to_use, "Density_2")]
                                elif draw_rad_height_PDF == 3:
                                        p55 = PhasePlot(sp, (PartType_Gas_to_use, "particle_position_cylindrical_radius"), (PartType_Gas_to_use, "particle_position_cylindrical_z_abs"), (PartType_Gas_to_use, "Density_2_minus_analytic"), weight_field=(PartType_Gas_to_use, "Mass_2"), fontsize=12, x_bins=200, y_bins=200)
                                        p55.set_zlim((PartType_Gas_to_use, "Density_2_minus_analytic"), -1e-24, 1e-24)
                                        p55.set_cmap((PartType_Gas_to_use, "Density_2_minus_analytic"), 'algae')
                                        p55.set_log("particle_position_cylindrical_radius", False)
                                        p55.set_log("particle_position_cylindrical_z_abs", False)
                                        p55.set_unit("particle_position_cylindrical_radius", 'kpc')
                                        p55.set_unit("particle_position_cylindrical_z_abs", 'kpc')
                                        plot55 = p55.plots[(PartType_Gas_to_use, "Density_2_minus_analytic")]

                        p55.set_xlabel("Cylindrical Radius (kpc)")
                        p55.set_ylabel("Vertical Height (kpc)")
                        p55.set_xlim(0, 14)
                        p55.set_ylim(0, 1.4)

                        plot55.figure = fig_rad_height_PDF[time]
                        plot55.axes = grid_rad_height_PDF[time][code].axes
                        if code == 0: plot55.cax = grid_rad_height_PDF[time].cbar_axes[0]
                        p55._setup_plots()

                        # Add 1D profile line if requested
                        if draw_rad_height_PDF == 2:
                                if codes[code] == "ART-I" or codes[code] == "ART-II" or codes[code] == "ENZO"  or codes[code] == "RAMSES":
                                        p56 = ProfilePlot(sp, ("index", "cylindrical_r"),  ("index", "cylindrical_z_abs"), \
                                                                  weight_field=("gas", "cell_mass"), n_bins=50, x_log=False)
                                        p56.set_log("cylindrical_r", False)
                                        p56.set_log("cylindrical_z_abs", False)
                                        p56.set_unit("cylindrical_r", 'kpc')
                                        p56.set_unit("cylindrical_z_abs", 'kpc')
                                        p56.set_xlim(0, 14)
                                        p56.set_ylim("cylindrical_z_abs", 0, 1.4)
                                        line = ln.Line2D(p56.profiles[0].x.in_units('kpc'), p56.profiles[0]["cylindrical_z_abs"].in_units('kpc'), linestyle="-", linewidth=2, color='k', alpha=0.7)
                                        rad_height_xs[time].append(p56.profiles[0].x.in_units('kpc').d)
                                        rad_height_profiles[time].append(p56.profiles[0]["cylindrical_z_abs"].in_units('kpc').d)
                                else:
                                        p56 = ProfilePlot(sp, (PartType_Gas_to_use, "particle_position_cylindrical_radius"), (PartType_Gas_to_use, "particle_position_cylindrical_z_abs"), \
                                                                  weight_field=(PartType_Gas_to_use, MassType_to_use), n_bins=50, x_log=False)
                                        p56.set_log("particle_position_cylindrical_radius", False)
                                        p56.set_log("particle_position_cylindrical_z_abs", False)
                                        p56.set_unit("particle_position_cylindrical_radius", 'kpc')
                                        p56.set_unit("particle_position_cylindrical_z_abs", 'kpc')
                                        p56.set_xlim(0, 14)
                                        p56.set_ylim("particle_position_cylindrical_z_abs", 0, 1.4)
                                        line = ln.Line2D(p56.profiles[0].x.in_units('kpc'), p56.profiles[0]["particle_position_cylindrical_z_abs"].in_units('kpc'), linestyle="-", linewidth=2, color='k', alpha=0.7)
                                        rad_height_xs[time].append(p56.profiles[0].x.in_units('kpc').d)
                                        rad_height_profiles[time].append(p56.profiles[0]["particle_position_cylindrical_z_abs"].in_units('kpc').d)
                                grid_rad_height_PDF[time][code].axes.add_line(line)

                        if add_nametag == 1:
                                at = AnchoredText("%s" % codes[code], loc=1, prop=dict(size=10), frameon=True)
                                grid_rad_height_PDF[time][code].axes.add_artist(at)

                # DENSITY-TEMPERATURE-METALLICITY PDF
                if draw_metal_PDF == 1:
                        sp = pf.sphere(center, (0.5*figure_width, "kpc"))
                        my_cmap2 = copy.copy(matplotlib.cm.get_cmap('algae'))
                        my_cmap2.set_under('w')
                        if codes[code] == "ART-I" or codes[code] == "ART-II" or codes[code] == "ENZO"  or codes[code] == "RAMSES":
                                pf.field_info[("gas", "metallicity")].take_log = False
                                p3 = PhasePlot(sp, ("gas", "density"), ("gas", "temperature"), ("gas", "metallicity"), weight_field=("gas", "cell_mass"), fontsize=12, x_bins=300, y_bins=300)
                                p3.set_zlim(("gas", "metallicity"), 0.01, 0.04)
                                p3.set_cmap(("gas", "metallicity"), my_cmap2)
                                p3.set_colorbar_label(("gas", "metallicity"), "Metallicity (Mass-weighted average of mass fraction)")
                                plot3 = p3.plots[("gas", "metallicity")]
                        else:
                                # Because ParticlePhasePlot doesn't yet work for a log-log PDF for some reason, I will do the following trick.
                                pf.field_info[(PartType_Gas_to_use, "Metallicity_2")].take_log = False
                                p3 = PhasePlot(sp, (PartType_Gas_to_use, "Density_2"), (PartType_Gas_to_use, "Temperature_2"), (PartType_Gas_to_use, "Metallicity_2"), weight_field=(PartType_Gas_to_use, "Mass_2"), fontsize=12, x_bins=300, y_bins=300)
                                p3.set_zlim((PartType_Gas_to_use, "Metallicity_2"), 0.01, 0.04)
                                p3.set_cmap((PartType_Gas_to_use, "Metallicity_2"), my_cmap2)
                                plot3 = p3.plots[(PartType_Gas_to_use, "Metallicity_2")]

                        p3.set_xlim(1e-29, 1e-21)
                        p3.set_ylim(10, 1e7)

                        plot3.figure = fig_metal_PDF[time]
                        plot3.axes = grid_metal_PDF[time][code].axes
                        if code == 0: plot3.cax = grid_metal_PDF[time].cbar_axes[0]
                        p3._setup_plots()

                        if add_nametag == 1:
                                at = AnchoredText("%s" % codes[code], loc=3, prop=dict(size=10), frameon=True)
                                grid_metal_PDF[time][code].axes.add_artist(at)

                # DENSITY DF (DISTRIBUTION FUNCTION)
                if draw_density_DF >= 1:
                        sp = pf.sphere(center, (0.5*figure_width, "kpc"))
                        if codes[code] == "ART-I" or codes[code] == "ART-II" or codes[code] == "ENZO"  or codes[code] == "RAMSES":
                                p6 = ProfilePlot(sp, ("gas", "density"),  ("gas", "cell_mass"), weight_field=None, n_bins=50, x_log=True, accumulation=False)
#				p6 = ProfilePlot(sp, ("gas", "density"),  ("gas", "cell_mass"), weight_field=None, n_bins=50, x_log=True, accumulation=True)
                                p6.set_log("cell_mass", True)
                                p6.set_xlim(1e-29, 1e-21)
                                density_DF_xs[time].append(p6.profiles[0].x.in_units('g/cm**3').d)
                                density_DF_profiles[time].append(p6.profiles[0]["cell_mass"].in_units('Msun').d)
                        else:
                                # Because ParticleProfilePlot doesn't exist, I will do the following trick.
                                p6 = ProfilePlot(sp, (PartType_Gas_to_use, "Density_2"),  (PartType_Gas_to_use, "Mass_2"), weight_field=None, n_bins=50, x_log=True, accumulation=False)
#				p6 = ProfilePlot(sp, (PartType_Gas_to_use, "Density_2"),  (PartType_Gas_to_use, "Mass_2"), weight_field=None, n_bins=50, x_log=True, accumulation=True)
                                p6.set_log("Mass_2", True)
                                p6.set_xlim(1e-29, 1e-21)
                                density_DF_xs[time].append(p6.profiles[0].x.in_units('g/cm**3').d)
                                density_DF_profiles[time].append(p6.profiles[0]["Mass_2"].in_units('Msun').d)

                        # Add difference plot between 1st and 2nd datasets, if requested
                        if draw_density_DF == 2 and time != 0:
                                if dataset_num == 2:
                                        pf_1st = load_dataset(1, time, code, codes, filenames[0]) # load 1st datasets
                                        v, cen = pf_1st.h.find_max(("gas", "density"))
                                        sp = pf_1st.sphere(cen, (30.0, "kpc"))
                                        cen2 = sp.quantities.center_of_mass(use_gas=True, use_particles=False).in_units("kpc")
                                        sp2 = pf_1st.sphere(cen2, (1.0, "kpc"))
                                        cen3 = sp2.quantities.max_location(("gas", "density"))
                                        center_1st = pf_1st.arr([cen3[1].d, cen3[2].d, cen3[3].d], 'code_length')
                                        if yt_version_pre_3_2_3 == 1:
                                                center_1st = pf_1st.arr([cen3[2].d, cen3[3].d, cen3[4].d], 'code_length') # for yt-3.2.3 or before
                                        sp_1st = pf_1st.sphere(center_1st, (0.5*figure_width, "kpc"))
                                        if codes[code] == "ART-I" or codes[code] == "ART-II" or codes[code] == "ENZO"  or codes[code] == "RAMSES":
                                                p61 = ProfilePlot(sp_1st, ("gas", "density"),  ("gas", "cell_mass"), weight_field=None, n_bins=50, x_log=True, accumulation=False)
                                                p61.set_log("cell_mass", True)
                                                p61.set_xlim(1e-29, 1e-21)
                                                density_DF_1st_xs[time].append(p61.profiles[0].x.in_units('g/cm**3').d)
                                                density_DF_1st_profiles[time].append(p61.profiles[0]["cell_mass"].in_units('Msun').d)
                                        else:
                                                def _Density_2(field, data):
                                                        return data[(PartType_Gas_to_use, "Density")].in_units('g/cm**3')
                                                pf_1st.add_field((PartType_Gas_to_use, "Density_2"), function=_Density_2, take_log=True, particle_type=False, display_name="Density", units="g/cm**3")
                                                def _Mass_2(field, data):
                                                        return data[(PartType_Gas_to_use, MassType_to_use)].in_units('Msun')
                                                pf_1st.add_field((PartType_Gas_to_use, "Mass_2"), function=_Mass_2, take_log=True, particle_type=False, display_name="Mass", units="Msun")
                                                p61 = ProfilePlot(sp_1st, (PartType_Gas_to_use, "Density_2"),  (PartType_Gas_to_use, "Mass_2"), weight_field=None, n_bins=50, x_log=True, accumulation=False)
                                                p61.set_log("Mass_2", True)
                                                p61.set_xlim(1e-29, 1e-21)
                                                density_DF_1st_xs[time].append(p61.profiles[0].x.in_units('g/cm**3').d)
                                                density_DF_1st_profiles[time].append(p61.profiles[0]["Mass_2"].in_units('Msun').d)
                                else:
                                        print("This won't work; consider setting dataset_num to 2...")
                                        continue

                # CYLINDRICAL RADIUS DF + RADIALLY-BINNED GAS SURFACE DENSITY
                if draw_radius_DF == 1:
                        sp = pf.sphere(center, (0.5*figure_width, "kpc"))
                        sp.set_field_parameter("normal", disk_normal_vector)
                        if codes[code] == "ART-I" or codes[code] == "ART-II" or codes[code] == "ENZO"  or codes[code] == "RAMSES":
                                p7 = ProfilePlot(sp, ("index", "cylindrical_r"),  ("gas", "cell_mass"), weight_field=None, n_bins=50, x_log=False, accumulation=False)
                                p7.set_log("cell_mass", True)
                                p7.set_log("cylindrical_r", False)
                                p7.set_unit("cylindrical_r", 'kpc')
                                p7.set_xlim(1e-3, 15)
                                radius_DF_xs[time].append(p7.profiles[0].x.in_units('kpc').d)
                                radius_DF_profiles[time].append(p7.profiles[0]["cell_mass"].in_units('Msun').d)
                        else:
                                p7 = ProfilePlot(sp, (PartType_Gas_to_use, "particle_position_cylindrical_radius"),  (PartType_Gas_to_use, "Mass_2"), weight_field=None, n_bins=50, x_log=False, accumulation=False)
                                p7.set_log("Mass_2", True)
                                p7.set_log("particle_position_cylindrical_radius", False)
                                p7.set_unit("particle_position_cylindrical_radius", 'kpc')
                                p7.set_xlim(1e-3, 15)
                                radius_DF_xs[time].append(p7.profiles[0].x.in_units('kpc').d)
                                radius_DF_profiles[time].append(p7.profiles[0]["Mass_2"].in_units('Msun').d)

                # CYLINDRICAL RADIUS DF + RADIALLY-BINNED SURFACE DENSITY FOR NEW STARS
                if draw_star_radius_DF >= 1 and time != 0:
                        sp = pf.sphere(center, (0.5*figure_width, "kpc"))
                        sp.set_field_parameter("normal", disk_normal_vector)
                        pf.field_info[(PartType_Star_to_use, "particle_mass")].take_log = True
                        pf.field_info[(PartType_Star_to_use, "particle_mass")].output_units = 'code_mass' # this turned out to be crucial!; check output_units above
                        pf.field_info[(PartType_Star_to_use, "particle_position_cylindrical_radius")].take_log = False
                        p71 = ProfilePlot(sp, (PartType_Star_to_use, "particle_position_cylindrical_radius"),  (PartType_Star_to_use, "particle_mass"), weight_field=None, n_bins=50, x_log=False, accumulation=False)
                        p71.set_unit("particle_position_cylindrical_radius", 'kpc')
                        p71.set_xlim(1e-3, 15)
                        star_radius_DF_xs[time].append(p71.profiles[0].x.in_units('kpc').d)
                        star_radius_DF_profiles[time].append(p71.profiles[0]["particle_mass"].in_units('Msun').d)

                        # Add RADIALLY-BINNED SFR SURFACE DENSITY PROFILE if requested (SFR estimated using stars younger than 20 Myrs old)
                        if draw_star_radius_DF == 2 and time != 0:
                                def _particle_mass_young_stars_star_radius_DF(field, data):
                                        trans = np.zeros(data[(PartType_StarBeforeFiltered_to_use, "particle_mass")].shape)
                                        ind = np.where(data[(PartType_StarBeforeFiltered_to_use, FormationTimeType_to_use)].in_units('Myr') > (pf.current_time.in_units('Myr').d - young_star_cutoff_star_radius_DF))
                                        trans[ind] = data[(PartType_StarBeforeFiltered_to_use, "particle_mass")][ind].in_units('code_mass')
                                        return data.ds.arr(trans, "code_mass").in_base(data.ds.unit_system.name)
                                pf.add_field((PartType_StarBeforeFiltered_to_use, "particle_mass_young_stars_star_radius_DF"), function=_particle_mass_young_stars_star_radius_DF, \
                                                     display_name="Young Stellar Mass", units='code_mass', particle_type=True, take_log=True)
                                pf.add_particle_filter(PartType_Star_to_use)

                                pf.field_info[(PartType_Star_to_use, "particle_mass_young_stars_star_radius_DF")].take_log = True
                                pf.field_info[(PartType_Star_to_use, "particle_mass_young_stars_star_radius_DF")].output_units = 'code_mass' # this turned out to be crucial!; check output_units above
                                pf.field_info[(PartType_Star_to_use, "particle_position_cylindrical_radius")].take_log = False
                                p72 = ProfilePlot(sp, (PartType_Star_to_use, "particle_position_cylindrical_radius"), (PartType_Star_to_use, "particle_mass_young_stars_star_radius_DF"), \
                                                          weight_field=None, n_bins=50, x_log=False, accumulation=False)
                                p72.set_unit("particle_position_cylindrical_radius", 'kpc')
                                p72.set_xlim(1e-3, 15)
                                sfr_radius_DF_xs[time].append(p72.profiles[0].x.in_units('kpc').d)
                                sfr_radius_DF_profiles[time].append(p72.profiles[0]["particle_mass_young_stars_star_radius_DF"].in_units('Msun').d/young_star_cutoff_star_radius_DF/1e6) # in Msun/yr

                # VERTICAL HEIGHT DF + VERTICALLY-BINNED GAS SURFACE DENSITY
                if draw_height_DF == 1:
                        sp = pf.sphere(center, (0.5*figure_width, "kpc"))
                        sp.set_field_parameter("normal", disk_normal_vector)
                        if codes[code] == "ART-I" or codes[code] == "ART-II" or codes[code] == "ENZO"  or codes[code] == "RAMSES":
                                p8 = ProfilePlot(sp, ("index", "cylindrical_z_abs"),  ("gas", "cell_mass"), weight_field=None, n_bins=10, x_log=False, accumulation=False)
                                p8.set_log("cell_mass", True)
                                p8.set_log("cylindrical_z_abs", False)
                                p8.set_unit("cylindrical_z_abs", 'kpc')
                                p8.set_xlim(1e-3, 1.4)
                                height_DF_xs[time].append(p8.profiles[0].x.in_units('kpc').d)
                                height_DF_profiles[time].append(p8.profiles[0]["cell_mass"].in_units('Msun').d)
                        else:
                                p8 = ProfilePlot(sp, (PartType_Gas_to_use, "particle_position_cylindrical_z_abs"),  (PartType_Gas_to_use, "Mass_2"), weight_field=None, n_bins=10, x_log=False, accumulation=False)
                                p8.set_log("Mass_2", True)
                                p8.set_log("particle_position_cylindrical_z_abs", False)
                                p8.set_unit("particle_position_cylindrical_z_abs", 'kpc')
                                p8.set_xlim(1e-3, 1.4)
                                height_DF_xs[time].append(p8.profiles[0].x.in_units('kpc').d)
                                height_DF_profiles[time].append(p8.profiles[0]["Mass_2"].in_units('Msun').d)

                # STAR FORMATION RATE + CUMULATIVE STELLAR MASS GROWTH IN TIME
                if draw_SFR >= 1 and time != 0:
                        from yt.units.dimensions import length # Below are tricks to make StarFormationRate() work, particularly with "volume" argument, as it currently works only with comoving datasets
                        pf.unit_registry.add('pccm', pf.unit_registry.lut['pc'][0], length, "\\rm{pc}/(1+z)")
                        pf.hubble_constant = 0.71; pf.omega_lambda = 0.73; pf.omega_matter = 0.27; pf.omega_curvature = 0.0

                        sp = pf.sphere(center, (0.5*figure_width, "kpc"))
                        draw_SFR_mass = sp[(PartType_Star_to_use, "particle_mass")].in_units('Msun')
                        draw_SFR_ct   = sp[(PartType_Star_to_use, FormationTimeType_to_use)].in_units('Myr')
                        sfr = StarFormationRate(pf, star_mass = draw_SFR_mass, star_creation_time = draw_SFR_ct,
                                                volume = sp.volume(), bins = 25) # 25 bins hardcoded; see: http://yt-project.org/docs/dev/analyzing/analysis_modules/star_analysis.html

                        sfr_ts[time].append(sfr.time.in_units('Myr')) # in Myr
                        sfr_cum_masses[time].append(sfr.Msol_cumulative) # in Msun
                        sfr_sfrs[time].append(sfr.Msol_yr) # in Msun/yr

                # DENSITY ALONG THE ORTHO-RAY OBJECT CUTTING THROUGH THE CENTER
                if draw_cut_through == 1:
                        ray = pf.ortho_ray(2, (center[0].in_units('code_length'), center[1].in_units('code_length'))) # see: http://yt-project.org/doc/visualizing/manual_plotting.html#line-plots
                        srt = np.argsort(ray['z'])
                        cut_through_zs[time].append(np.array(ray['z'][srt].in_units('kpc').d - center[2].in_units('kpc').d))
                        cut_through_zvalues[time].append(np.array(ray[("gas", "density")][srt].in_units('g/cm**3').d))

                        ray = pf.ortho_ray(0, (center[1].in_units('code_length'), center[2].in_units('code_length')))
                        srt = np.argsort(ray['x'])
                        cut_through_xs[time].append(np.array(ray['x'][srt].in_units('kpc').d - center[2].in_units('kpc').d))
                        cut_through_xvalues[time].append(np.array(ray[("gas", "density")][srt].in_units('g/cm**3').d))

        ####################################
        #        POST-ANALYSIS STEPS       #
        ####################################

        # SAVE FIGURES
        if draw_density_map == 1:
                fig_density_map[time].savefig("Sigma_%dMyr" % times[time], bbox_inches='tight', pad_inches=0.03, dpi=300)
        if draw_temperature_map == 1:
                fig_temperature_map[time].savefig("Temp_%dMyr" % times[time], bbox_inches='tight', pad_inches=0.03, dpi=300)
        if draw_cellsize_map == 1 or draw_cellsize_map == 3:
                fig_cellsize_map[time].savefig("Cell_%dMyr" % times[time], bbox_inches='tight', pad_inches=0.03, dpi=300)
        if draw_cellsize_map == 2 or draw_cellsize_map == 3:
                fig_cellsize_map_2[time].savefig("Resolution_%dMyr" % times[time], bbox_inches='tight', pad_inches=0.03, dpi=300)
        if draw_elevation_map == 1:
                fig_elevation_map[time].savefig("Elevation_%dMyr" % times[time], bbox_inches='tight', pad_inches=0.03, dpi=300)
        if draw_metal_map >= 1:
                fig_metal_map[time].savefig("Metal_%dMyr" % times[time], bbox_inches='tight', pad_inches=0.03, dpi=300)
        if draw_zvel_map == 1:
                fig_zvel_map[time].savefig("z-vel_%dMyr" % times[time], bbox_inches='tight', pad_inches=0.03, dpi=300)
        if draw_star_map == 1 and time != 0:
                fig_star_map[time].savefig("Star_%dMyr" % times[time], bbox_inches='tight', pad_inches=0.03, dpi=300)
        if draw_star_clump_stats >= 1 and time != 0:
                if draw_star_clump_stats >= 2:
                        fig_star_map_2[time].savefig("Star_with_clumps_HOP_%dMyr" % times[time], bbox_inches='tight', pad_inches=0.03, dpi=300)
                        fig_star_map_3[time].savefig("Star_with_clumps_FOF_%dMyr" % times[time], bbox_inches='tight', pad_inches=0.03, dpi=300)
                if draw_star_clump_stats == 3:
                        pf_clump_ref = load_dataset(2, 1, 0, ['GIZMO'],
                                                    [['dummy', '/lustre/ki/pfs/mornkr/080112_CHaRGe/pfs-hyades/AGORA-DISK-repository/Grackle+SF+ThermalFbck/GIZMO/AGORA_disk_second_ps1.5/snapshot_temp_100_last']])
                        PartType_Star_to_use = "Stars"
                        if os.path.exists("./halo_catalogs/hop_%s_%05d_high_floor/hop_%s_%05d_high_floor.0.h5" % ('GIZMO', 500, 'GIZMO', 500)) == False:
                                hc_clump_ref = HaloCatalog(data_ds=pf_clump_ref, finder_method='hop', output_dir="./halo_catalogs/hop_%s_%05d_high_floor" % ('GIZMO', 500), \
                                                         finder_kwargs={'threshold': 2e8, 'dm_only': False, 'ptype': PartType_Star_to_use})
                                hc_clump_ref.add_filter('quantity_value', 'particle_mass', '>', 2.6e6, 'Msun')
                                hc_clump_ref.add_filter('quantity_value', 'particle_mass', '<', 2.6e8, 'Msun')
                                hc_clump_ref.create()
                        halo_ds_clump_ref = load("./halo_catalogs/hop_%s_%05d_high_floor/hop_%s_%05d_high_floor.0.h5" % ('GIZMO', 500, 'GIZMO', 500))
                        hc_clump_ref = HaloCatalog(halos_ds=halo_ds_clump_ref, output_dir="./halo_catalogs/hop_%s_%05d_high_floor" % ('GIZMO', 500))
                        hc_clump_ref.load()

                        halo_ad_clump_ref = hc_clump_ref.halos_ds.all_data()
                        star_clump_masses_hop_ref.append(np.log10(halo_ad_clump_ref['particle_mass'][:].in_units("Msun")))

                        if os.path.exists("./halo_catalogs/fof_%s_%05d_high_floor/fof_%s_%05d_high_floor.0.h5" % ('GIZMO', 500, 'GIZMO', 500)) == False:
                                hc2_clump_ref = HaloCatalog(data_ds=pf_clump_ref, finder_method='fof', output_dir="./halo_catalogs/fof_%s_%05d_high_floor" %  ('GIZMO', 500), \
                                                          finder_kwargs={'link': 0.0025, 'dm_only': False, 'ptype': PartType_Star_to_use})
                                hc2_clump_ref.add_filter('quantity_value', 'particle_mass', '>', 2.6e6, 'Msun')
                                hc2_clump_ref.add_filter('quantity_value', 'particle_mass', '<', 2.6e8, 'Msun')
                                hc2_clump_ref.create()
                        halo_ds2_clump_ref = load("./halo_catalogs/fof_%s_%05d_high_floor/fof_%s_%05d_high_floor.0.h5" % ('GIZMO', 500, 'GIZMO', 500))
                        hc2_clump_ref = HaloCatalog(halos_ds=halo_ds2_clump_ref, output_dir="./halo_catalogs/fof_%s_%05d_high_floor" % ('GIZMO', 500))
                        hc2_clump_ref.load()

                        halo_ad2_clump_ref = hc2_clump_ref.halos_ds.all_data()
                        star_clump_masses_fof_ref.append(np.log10(halo_ad2_clump_ref['particle_mass'][:].in_units("Msun")))
                plt.clf()
                fig = plt.figure(figsize=(8, 4))
                gridspec.GridSpec(1, 2)
                for star_clump_stats_i in range(1,3,1):
                        codes_plotted = []
#			plt.subplot(1,2,star_clump_stats_i, aspect=0.25)
                        plt.subplot2grid((1, 2), (0, star_clump_stats_i-1))
                        plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=0.5)
                        for code in range(len(codes)):
                                if (star_clump_stats_i == 1 and (codes[code] == "ART-I" or codes[code] == "ART-II"  or codes[code] == "ENZO" or codes[code] == "RAMSES")) or \
                                   (star_clump_stats_i == 2 and (codes[code] == "CHANGA" or codes[code] == "GASOLINE" or codes[code] == "GADGET-3" or codes[code] == "GEAR" or codes[code] == "GIZMO" or codes[code] == "SWIFT")):
                                        hist = np.histogram(star_clump_masses_hop[time][code], bins=10, range=(6., 8.5))
                                        dbin = 0.5*(hist[1][1] - hist[1][0])
                                        lines = plt.plot(hist[1][:-1]+dbin, hist[0], color=color_names[code], linestyle=linestyle_names[np.mod(code, len(linestyle_names))], marker=marker_names[code], markeredgecolor='none', linewidth=1.2, alpha=0.8)
                                        codes_plotted.append(codes[code])
                        plt.xlim(6, 8.5)
                        plt.ylim(-0.1, 15)
                        plt.grid(True)
                        plt.xlabel("$\mathrm{log[Newly\ Formed\ Stellar\ Clump\ Mass\ (M_{\odot})]}$")
                        if star_clump_stats_i == 1:
                                plt.ylabel("$\mathrm{Stellar\ Clump\ Counts, \ \ N_{clump}(M)}$")
                        plt.legend(codes_plotted, loc=1, frameon=True, ncol=1, fancybox=True)
                        leg = plt.gca().get_legend()
                        ltext = leg.get_texts()
                        plt.setp(ltext, fontsize='xx-small')
                plt.savefig("star_clump_stats_HOP_%dMyr" % times[time], bbox_inches='tight', pad_inches=0.03, dpi=300)
                plt.clf()
                # Reiterate for cumulative plots
                fig = plt.figure(figsize=(8, 4))
                gridspec.GridSpec(1, 2)
                for star_clump_stats_i in range(1,3,1):
                        codes_plotted = []
#			plt.subplot(1,2,star_clump_stats_i, aspect=0.15)
                        plt.subplot2grid((1, 2), (0, star_clump_stats_i-1))
                        plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=0.5)
                        for code in range(len(codes)):
                                if (star_clump_stats_i == 1 and (codes[code] == "ART-I" or codes[code] == "ART-II"  or codes[code] == "ENZO" or codes[code] == "RAMSES")) or \
                                   (star_clump_stats_i == 2 and (codes[code] == "CHANGA" or codes[code] == "GASOLINE" or codes[code] == "GADGET-3" or codes[code] == "GEAR" or codes[code] == "GIZMO" or codes[code] == "SWIFT")):
                                        hist = np.histogram(star_clump_masses_hop[time][code], bins=10, range=(6., 8.5))
                                        dbin = 0.5*(hist[1][1] - hist[1][0])
                                        lines = plt.plot(hist[1][:-1]+dbin, np.cumsum(hist[0][::-1])[::-1], color=color_names[code], linestyle=linestyle_names[np.mod(code, len(linestyle_names))], \
                                                                 marker=marker_names[code], markeredgecolor='none', linewidth=1.2, alpha=0.8)
                                        codes_plotted.append(codes[code])
                        if draw_star_clump_stats == 3 and star_clump_stats_i == 2:
                                hist_clump_ref = np.histogram(star_clump_masses_hop_ref, bins=10, range=(6., 8.5))
                                dbin = 0.5*(hist[1][1] - hist[1][0])
                                lines = plt.plot(hist_clump_ref[1][:-1]+dbin, np.cumsum(hist_clump_ref[0][::-1])[::-1], color='k', linestyle='--', marker='*', markeredgecolor='none', linewidth=1.2, alpha=0.8)
                                codes_plotted.append('GIZMO-ps2')
                        plt.xlim(6, 8.5)
#			plt.ylim(-0.1, 30)
                        plt.ylim(0.9, 50)
                        plt.semilogy()
                        plt.grid(True)
                        plt.xlabel("$\mathrm{log[Newly\ Formed\ Stellar\ Clump\ Mass\ (M_{\odot})]}$")
                        if star_clump_stats_i == 1:
                                plt.ylabel("$\mathrm{Cumulative\ Clump\ Counts, \ \ N_{clump}(>M)}}$")
                        plt.legend(codes_plotted, loc=1, frameon=True, ncol=1, fancybox=True, labelspacing=0.15, borderpad=0.3)
                        leg = plt.gca().get_legend()
                        ltext = leg.get_texts()
                        plt.setp(ltext, fontsize='xx-small')
                plt.savefig("star_clump_stats_HOP_cumulative_%dMyr" % times[time], bbox_inches='tight', pad_inches=0.03, dpi=300)
                plt.clf()

                fig = plt.figure(figsize=(8, 4))
                gridspec.GridSpec(1, 2)
                for star_clump_stats_i in range(1,3,1):
                        codes_plotted = []
#			plt.subplot(1,2,star_clump_stats_i, aspect=0.25)
                        plt.subplot2grid((1, 2), (0, star_clump_stats_i-1))
                        plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=0.5)
                        for code in range(len(codes)):
                                if (star_clump_stats_i == 1 and (codes[code] == "ART-I" or codes[code] == "ART-II"  or codes[code] == "ENZO" or codes[code] == "RAMSES")) or \
                                   (star_clump_stats_i == 2 and (codes[code] == "CHANGA" or codes[code] == "GASOLINE" or codes[code] == "GADGET-3" or codes[code] == "GEAR" or codes[code] == "GIZMO" or codes[code] == "SWIFT")):
                                        hist = np.histogram(star_clump_masses_fof[time][code], bins=10, range=(6., 8.5))
                                        dbin = 0.5*(hist[1][1] - hist[1][0])
                                        lines = plt.plot(hist[1][:-1]+dbin, hist[0], color=color_names[code], linestyle=linestyle_names[np.mod(code, len(linestyle_names))], marker=marker_names[code], markeredgecolor='none', linewidth=1.2, alpha=0.8)
                                        codes_plotted.append(codes[code])
                        plt.xlim(6, 8.5)
                        plt.ylim(-0.1, 15)
                        plt.grid(True)
                        plt.xlabel("$\mathrm{log[Newly\ Formed\ Stellar\ Clump\ Mass\ (M_{\odot})]}$")
                        if star_clump_stats_i == 1:
                                plt.ylabel("$\mathrm{Stellar\ Clump\ Counts, \ \ N_{clump}(M)}$")
                        plt.legend(codes_plotted, loc=1, frameon=True, ncol=1, fancybox=True)
                        leg = plt.gca().get_legend()
                        ltext = leg.get_texts()
                        plt.setp(ltext, fontsize='xx-small')
                plt.savefig("star_clump_stats_FOF_%dMyr" % times[time], bbox_inches='tight', pad_inches=0.03, dpi=300)
                plt.clf()
                # Reiterate for cumulative plots
                fig = plt.figure(figsize=(8, 4))
                gridspec.GridSpec(1, 2)
                for star_clump_stats_i in range(1,3,1):
                        codes_plotted = []
#			plt.subplot(1,2,star_clump_stats_i, aspect=0.15)
                        plt.subplot2grid((1, 2), (0, star_clump_stats_i-1))
                        plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=0.5)
                        for code in range(len(codes)):
                                if (star_clump_stats_i == 1 and (codes[code] == "ART-I" or codes[code] == "ART-II"  or codes[code] == "ENZO" or codes[code] == "RAMSES")) or \
                                   (star_clump_stats_i == 2 and (codes[code] == "CHANGA" or codes[code] == "GASOLINE" or codes[code] == "GADGET-3" or codes[code] == "GEAR" or codes[code] == "GIZMO" or codes[code] == "SWIFT")):
                                        hist = np.histogram(star_clump_masses_fof[time][code], bins=10, range=(6., 8.5))
                                        dbin = 0.5*(hist[1][1] - hist[1][0])
                                        lines = plt.plot(hist[1][:-1]+dbin, np.cumsum(hist[0][::-1])[::-1], color=color_names[code], linestyle=linestyle_names[np.mod(code, len(linestyle_names))], \
                                                                 marker=marker_names[code], markeredgecolor='none', linewidth=1.2, alpha=0.8)
                                        codes_plotted.append(codes[code])
                        if draw_star_clump_stats == 3 and star_clump_stats_i == 2:
                                hist_clump_ref = np.histogram(star_clump_masses_fof_ref, bins=10, range=(6., 8.5))
                                dbin = 0.5*(hist[1][1] - hist[1][0])
                                lines = plt.plot(hist_clump_ref[1][:-1]+dbin, np.cumsum(hist_clump_ref[0][::-1])[::-1], color='k', linestyle='--', marker='*', markeredgecolor='none', linewidth=1.2, alpha=0.8)
                                codes_plotted.append('GIZMO-ps2')
                        plt.xlim(6, 8.5)
#			plt.ylim(-0.1, 30)
                        plt.ylim(0.9, 50)
                        plt.semilogy()
                        plt.grid(True)
                        plt.xlabel("$\mathrm{log[Newly\ Formed\ Stellar\ Clump\ Mass\ (M_{\odot})]}$")
                        if star_clump_stats_i == 1:
                                plt.ylabel("$\mathrm{Cumulative\ Clump\ Counts, \ \ N_{clump}(>M)}}$")
                        plt.legend(codes_plotted, loc=1, frameon=True, ncol=1, fancybox=True, labelspacing=0.15, borderpad=0.3)
                        leg = plt.gca().get_legend()
                        ltext = leg.get_texts()
                        plt.setp(ltext, fontsize='xx-small')
                plt.savefig("star_clump_stats_FOF_cumulative_%dMyr" % times[time], bbox_inches='tight', pad_inches=0.03, dpi=300)
                plt.clf()
        if draw_SFR_map >= 1 and time != 0:
                fig_degr_density_map[time].savefig("degraded_Sigma_%dMyr" % times[time], bbox_inches='tight', pad_inches=0.03, dpi=300)
                fig_degr_sfr_map[time].savefig("degraded_SFR_map_%dMyr" % times[time], bbox_inches='tight', pad_inches=0.03, dpi=300)
                if draw_SFR_map == 2 and time != 0:
                        # If requested, draw the local K-S plot using the degraded maps created above
                        plt.clf()
#			plt.subplot(111)
                        fig = plt.figure(figsize=(8, 7))
                        Bigiel_surface_density = []
                        Bigiel_sfr_surface_density = []
                        for code in range(len(codes)):
                                # Remove bins where SFR surface density is zero
                                KS_x = np.array(surf_dens_SFR_map[time][code])
                                KS_x[np.where(np.array(sfr_surf_dens_SFR_map[time][code]) < 1e-10)] = 0
                                KS_x = np.log10(KS_x)
                                KS_y = np.log10(np.array(sfr_surf_dens_SFR_map[time][code]))
                                KS_x = KS_x[~np.isinf(KS_x)]
                                KS_y = KS_y[~np.isinf(KS_y)]
                                if codes[code] == "CHANGA":
                                        plt.scatter(KS_x, KS_y, color='k', edgecolor='k', s=20, linewidth=0.7, marker=marker_names[code], alpha=0.1)
                                # Draw a 80% contour rather than scattering all the datapoints; see http://stackoverflow.com/questions/19390320/scatterplot-contours-in-matplotlib
                                Gaussian_density_estimation_nbins = 20
                                kernel = kde.gaussian_kde(np.vstack([KS_x, KS_y]))
                                xi, yi = np.mgrid[KS_x.min():KS_x.max():Gaussian_density_estimation_nbins*1j, KS_y.min():KS_y.max():Gaussian_density_estimation_nbins*1j]
                                zi = np.reshape(kernel(np.vstack([xi.flatten(), yi.flatten()])), xi.shape)
                                contours = plt.contour(xi, yi, zi, np.array([0.2]), colors=color_names[code], linewidths=1.2, alpha=1.0)    # 80%    percentile contour
#				contours = plt.contour(xi, yi, zi, np.array([0.3173]), colors=color_names[code], linewidths=1.2, alpha=1.0) # 68.27% percentile contour (1-sigma)
                                contours.collections[0].set_label(codes[code]) # setting names for legend
                                plt.clabel(contours, fmt=codes[code], inline=True, fontsize=7) # setting labels
                        plt.xlim(0, 3)
                        plt.ylim(-4, 1)
                        plt.grid(True)
                        plt.xlabel("$\mathrm{log[Gas\ Surface\ Density\ (M_{\odot}/pc^2)]}$")
                        plt.ylabel("$\mathrm{log[Star\ Formation\ Rate\ Surface\ Density\ (M_{\odot}/yr/kpc^2)]}$")
                        plt.legend(loc=2, frameon=True, ncol=2, fancybox=True) # note difference from others since we want the legends to list contours, not scattered points
                        leg = plt.gca().get_legend()
                        ltext = leg.get_texts()
                        plt.setp(ltext, fontsize='small')
                        t = np.arange(-2, 5, 0.01)
                        for line in import_text("./Bigieletal2008_Fig8_Contour.txt", " "):
                                Bigiel_surface_density.append(float(line[0]) + np.log10(1.36)) # for the factor 1.36, see Section 2.3.1 of Bigiel et al. 2008
                                Bigiel_sfr_surface_density.append(float(line[1]))
                        plt.fill(Bigiel_surface_density, Bigiel_sfr_surface_density, fill=True, color='b', alpha = 0.1, hatch='\\') # contour by Bigiel et al. 2008 (Fig 8), or Feldmann et al. 2012 (Fig 1)
                        plt.plot(t, 1.37*t - 3.78, 'k--', linewidth = 2, alpha = 0.7) # observational fit by Kennicutt et al. 2007
#			plt.axhline(y = np.log10(8.593e4/young_star_cutoff_SFR_map/1e6/(float(aperture_size_SFR_map)/1000.)**2),
#                                   color='k', linestyle ='-.', linewidth=1, alpha=0.7) # SFR surface density cutoff due to limited star particle mass resolution
                        plt.savefig("K-S_local_%dMyr" % times[time], bbox_inches='tight', pad_inches=0.03, dpi=300)
                        plt.clf()
        if draw_PDF >= 1:
                if draw_PDF == 3 and time != 0:
                        # If requested, find a 1D profile line from a specific snapshot as a reference (in this case ENZO or CHANGA noSF run at 500 Myr), then we add it to all panels
                        # pf_profile_ref = load_dataset(1, 1, 0, ['ENZO'], [['dummy', file_location[0]+'ENZO/DD0100/DD0100']])
                        pf_profile_ref = load_dataset(1, 1, 0, ['CHANGA'], [['dummy', file_location[0]+'CHANGA/disklow/disklow.000500']])
                        def _Density_2(field, data):
                                return data[(PartType_Gas_to_use, "Density")].in_units('g/cm**3')
                        pf_profile_ref.add_field((PartType_Gas_to_use, "Density_2"), function=_Density_2, take_log=True, particle_type=False, display_name="Density", units="g/cm**3")
                        def _Temperature_2(field, data):
                                return data[(PartType_Gas_to_use, "Temperature")].in_units('K')
                        pf_profile_ref.add_field((PartType_Gas_to_use, "Temperature_2"), function=_Temperature_2, take_log=True, particle_type=False, display_name="Temperature", units="K")
                        def _Mass_2(field, data):
                                return data[(PartType_Gas_to_use, MassType_to_use)].in_units('Msun')
                        pf_profile_ref.add_field((PartType_Gas_to_use, "Mass_2"), function=_Mass_2, take_log=True, particle_type=False, display_name="Mass", units="Msun")

                        v, cen = pf_profile_ref.h.find_max(("gas", "density"))
                        sp = pf_profile_ref.sphere(cen, (30.0, "kpc"))
                        cen2 = sp.quantities.center_of_mass(use_gas=True, use_particles=False).in_units("kpc")
                        sp2 = pf_profile_ref.sphere(cen2, (1.0, "kpc"))
                        cen3 = sp2.quantities.max_location(("gas", "density"))
                        center_profile_ref = pf_profile_ref.arr([cen3[1].d, cen3[2].d, cen3[3].d], 'code_length')
                        if yt_version_pre_3_2_3 == 1:
                                center_profile_ref = pf_profile_ref.arr([cen3[2].d, cen3[3].d, cen3[4].d], 'code_length')
                        sp_profile_ref = pf_profile_ref.sphere(center_profile_ref, (0.5*figure_width, "kpc"))

                        # p35 = ProfilePlot(sp_profile_ref, ("gas", "density"),  ("gas", "temperature"), weight_field=("gas", "cell_mass"), n_bins=30)
                        # p35.set_xlim(1e-29, 1e-21)
                        # p35.set_ylim("temperature", 10, 1e7)
                        p35 = ProfilePlot(sp_profile_ref, (PartType_Gas_to_use, "Density_2"), (PartType_Gas_to_use, "Temperature_2"), weight_field=(PartType_Gas_to_use, "Mass_2"), n_bins=30)
                        p35.set_xlim(1e-26, 1e-21)
                        p35.set_ylim("Temperature_2", 10, 1e7)
                        for code in range(len(codes)):
#				line_profile_ref = ln.Line2D(p35.profiles[0].x.in_units('g/cm**3'), p35.profiles[0]["temperature"].in_units('K'), linestyle="--", linewidth=2, color='k', alpha=0.7)
                                line_profile_ref = ln.Line2D(p35.profiles[0].x.in_units('g/cm**3'), p35.profiles[0]["Temperature_2"].in_units('K'), linestyle="--", linewidth=2, color='k', alpha=0.7)
                                grid_PDF[time][code].axes.add_line(line_profile_ref)
                fig_PDF[time].savefig("PDF_%dMyr" % times[time], bbox_inches='tight', pad_inches=0.03, dpi=300)
        if draw_pos_vel_PDF >= 1:
                fig_pos_vel_PDF[time].savefig("pos_vel_PDF_%dMyr" % times[time], bbox_inches='tight', pad_inches=0.03, dpi=300)
                if draw_pos_vel_PDF >= 2 and time != 0:
                        plt.clf()
#			plt.subplot(111, aspect=0.04)
                        fig = plt.figure(figsize=(8, 8))
                        gridspec.GridSpec(4, 1)
                        plt.subplot2grid((4, 1), (0, 0), rowspan=3)
                        plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=0.5)
                        for code in range(len(codes)):
                                lines = plt.plot(pos_vel_xs[time][code], pos_vel_profiles[time][code], color=color_names[code], linestyle=linestyle_names[np.mod(code, len(linestyle_names))], marker=marker_names[code], markeredgecolor='none', linewidth=1.2, alpha=0.8)
                        plt.xlim(0, 14)
                        plt.ylim(0, 270)
                        plt.grid(True)
                        plt.xlabel("$\mathrm{Cylindrical\ Radius\ (kpc)}$", fontsize='large')
                        plt.ylabel("$\mathrm{Rotational\ Velocity\ (km/s)}$", fontsize='large')
                        plt.legend(codes, loc=4, frameon=True, ncol=2, fancybox=True)
                        leg = plt.gca().get_legend()
                        ltext = leg.get_texts()
                        plt.setp(ltext, fontsize='small')

                        plt.subplot2grid((4, 1), (3, 0), rowspan=1)
                        ave_profiles = np.mean(np.array(pos_vel_profiles[time]), axis=0)
                        for code in range(len(codes)):
                                lines = plt.plot(pos_vel_xs[time][code], (pos_vel_profiles[time][code] - ave_profiles)/ave_profiles, color=color_names[code], linestyle='none', marker=marker_names[code], markeredgecolor='none', linewidth=1.2, alpha=0.8)
                        plt.xlim(0, 14)
                        plt.ylim(-0.5, 0.5)
                        plt.grid(True)
                        plt.xlabel("$\mathrm{Cylindrical\ Radius\ (kpc)}$", fontsize='large')
                        plt.ylabel(r"$\mathrm{Residual,}\/(\mathrm{v} - \mathrm{\bar v})/\mathrm{\bar v}$", fontsize='large')
                        if add_mean_fractional_dispersion == 1:
                                mean_fractional_dispersion = (np.std(np.array(pos_vel_profiles[time]), axis=0)/ave_profiles)[(mean_dispersion_radius_range[0] < pos_vel_xs[time][0]) & (pos_vel_xs[time][0] < mean_dispersion_radius_range[1])].mean()
                                mean_fractional_dispersion_in_dex = np.log10(1.+mean_fractional_dispersion)
                                plt.annotate("Mean fractional dispersion for %d < r < %d kpc = %.3f (%.3f dex) " % (mean_dispersion_radius_range[0], mean_dispersion_radius_range[1], mean_fractional_dispersion, mean_fractional_dispersion_in_dex), xy=(0.35, 0.08), size=10, xycoords='axes fraction')
#				plt.text(12, 0.3, "mean dispersion = %.3f or %.3f dex " % (mean_fractional_dispersion, mean_fractional_dispersion_in_dex), ha="center", va="center", size=8, bbox=dict(boxstyle="round", fc="w", ec="0.5", alpha=0.9))
                        plt.savefig("pos_vel_%dMyr" % times[time], bbox_inches='tight', pad_inches=0.03, dpi=300)
                        plt.clf()
                if draw_pos_vel_PDF >= 3 and time != 0:
                        plt.clf()
#			plt.subplot(111, aspect=0.064)
                        fig = plt.figure(figsize=(8, 10))
                        gridspec.GridSpec(5, 1)
                        plt.subplot2grid((5, 1), (0, 0), rowspan=3)
                        plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=0.5)
                        for code in range(len(codes)):
                                lines = plt.plot(pos_disp_xs[time][code], pos_disp_profiles[time][code], color=color_names[code], linestyle=linestyle_names[np.mod(code, len(linestyle_names))], marker=marker_names[code], markeredgecolor='none', linewidth=1.2, alpha=0.8)
                        plt.xlim(0, 14)
                        plt.ylim(0, 170)
                        plt.grid(True)
                        plt.xlabel("$\mathrm{Cylindrical\ Radius\ (kpc)}$", fontsize='large')
                        plt.ylabel("$\mathrm{Velocity\ Dispersion\ (km/s)}$", fontsize='large')
                        plt.legend(codes, loc=1, frameon=True, ncol=2, fancybox=True)
                        leg = plt.gca().get_legend()
                        ltext = leg.get_texts()
                        plt.setp(ltext, fontsize='small')

                        plt.subplot2grid((5, 1), (3, 0), rowspan=1)
                        ave_profiles = np.mean(np.array(pos_disp_profiles[time]), axis=0)
                        for code in range(len(codes)):
                                lines = plt.plot(pos_disp_xs[time][code], (pos_disp_profiles[time][code] - ave_profiles)/ave_profiles, color=color_names[code], linestyle='none', marker=marker_names[code], markeredgecolor='none', linewidth=1.2, alpha=0.8)
                        plt.xlim(0, 14)
                        plt.ylim(-1.0, 1.0)
                        plt.grid(True)
                        plt.xlabel("$\mathrm{Cylindrical\ Radius\ (kpc)}$", fontsize='large')
                        plt.ylabel(r"$\mathrm{Residual,}\/(\mathrm{\sigma} - \mathrm{\bar \sigma})/\mathrm{\bar \sigma}$", fontsize='large')
                        if add_mean_fractional_dispersion == 1:
                                mean_fractional_dispersion = (np.std(np.array(pos_disp_profiles[time]), axis=0)/ave_profiles)[(mean_dispersion_radius_range[0] < pos_disp_xs[time][0]) & (pos_disp_xs[time][0] < mean_dispersion_radius_range[1])].mean()
                                mean_fractional_dispersion_in_dex = np.log10(1.+mean_fractional_dispersion)
                                plt.annotate("Mean fractional dispersion for %d < r < %d kpc = %.3f (%.3f dex) " % (mean_dispersion_radius_range[0], mean_dispersion_radius_range[1], mean_fractional_dispersion, mean_fractional_dispersion_in_dex), xy=(0.35, 0.08), size=10, xycoords='axes fraction')
                        plt.savefig("pos_disp_%dMyr" % times[time], bbox_inches='tight', pad_inches=0.03, dpi=300)
                        plt.subplot2grid((5, 1), (4, 0), rowspan=1)
                        for code in range(len(codes)):
                                lines = plt.plot(pos_disp_xs[time][code], pos_disp_vert_profiles[time][code]/pos_disp_profiles[time][code], color=color_names[code], linestyle='none', marker=marker_names[code], markeredgecolor='none', linewidth=1.2, alpha=0.8)
                        plt.xlim(0, 14)
                        plt.ylim(0.0, 1.0)
                        plt.grid(True)
                        plt.xlabel("$\mathrm{Cylindrical\ Radius\ (kpc)}$", fontsize='large')
                        plt.ylabel(r"$\mathrm{Ratio,}\/\mathrm{\sigma}_{\perp}/\mathrm{\sigma}$", fontsize='large')
                        plt.savefig("pos_disp_%dMyr" % times[time], bbox_inches='tight', pad_inches=0.03, dpi=300)
                        plt.clf()
                if draw_pos_vel_PDF == 4 and time != 0:
                        plt.clf()
#			plt.subplot(111, aspect=0.064)
                        fig = plt.figure(figsize=(8, 8))
                        gridspec.GridSpec(4, 1)
                        plt.subplot2grid((4, 1), (0, 0), rowspan=3)
                        plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=0.5)
                        for code in range(len(codes)):
                                lines = plt.plot(pos_disp_vert_xs[time][code], pos_disp_vert_profiles[time][code], color=color_names[code], linestyle=linestyle_names[np.mod(code, len(linestyle_names))], marker=marker_names[code], markeredgecolor='none', linewidth=1.2, alpha=0.8)
                        plt.xlim(0, 14)
                        plt.ylim(0, 170)
                        plt.grid(True)
                        plt.xlabel("$\mathrm{Cylindrical\ Radius\ (kpc)}$", fontsize='large')
                        plt.ylabel("$\mathrm{Vertical\ Velocity\ Dispersion\ (km/s)}$", fontsize='large')
                        plt.legend(codes, loc=1, frameon=True, ncol=2, fancybox=True)
                        leg = plt.gca().get_legend()
                        ltext = leg.get_texts()
                        plt.setp(ltext, fontsize='small')

                        plt.subplot2grid((4, 1), (3, 0), rowspan=1)
                        ave_profiles = np.mean(np.array(pos_disp_vert_profiles[time]), axis=0)
                        for code in range(len(codes)):
                                lines = plt.plot(pos_disp_vert_xs[time][code], (pos_disp_vert_profiles[time][code] - ave_profiles)/ave_profiles, color=color_names[code], linestyle='none', marker=marker_names[code], markeredgecolor='none', linewidth=1.2, alpha=0.8)
                        plt.xlim(0, 14)
                        plt.ylim(-1.0, 1.0)
                        plt.grid(True)
                        plt.xlabel("$\mathrm{Cylindrical\ Radius\ (kpc)}$", fontsize='large')
                        plt.ylabel(r"$\mathrm{Residual,}\/(\mathrm{\sigma} - \mathrm{\bar \sigma})/\mathrm{\bar \sigma}$", fontsize='large')
                        if add_mean_fractional_dispersion == 1:
                                mean_fractional_dispersion = (np.std(np.array(pos_disp_vert_profiles[time]), axis=0)/ave_profiles)[(mean_dispersion_radius_range[0] < pos_disp_vert_xs[time][0]) & (pos_disp_vert_xs[time][0] < mean_dispersion_radius_range[1])].mean()
                                mean_fractional_dispersion_in_dex = np.log10(1.+mean_fractional_dispersion)
                                plt.annotate("Mean fractional dispersion for %d < r < %d kpc = %.3f (%.3f dex) " % (mean_dispersion_radius_range[0], mean_dispersion_radius_range[1], mean_fractional_dispersion, mean_fractional_dispersion_in_dex), xy=(0.35, 0.08), size=10, xycoords='axes fraction')
                        plt.savefig("pos_disp_vert_%dMyr" % times[time], bbox_inches='tight', pad_inches=0.03, dpi=300)
                        plt.clf()
        if draw_star_pos_vel_PDF >= 1 and time != 0:
                fig_star_pos_vel_PDF[time].savefig("star_pos_vel_PDF_%dMyr" % times[time], bbox_inches='tight', pad_inches=0.03, dpi=300)
                if draw_star_pos_vel_PDF >= 2 and time != 0:
                        plt.clf()
#			plt.subplot(111, aspect=0.04)
                        fig = plt.figure(figsize=(8, 8))
                        gridspec.GridSpec(4, 1)
                        plt.subplot2grid((4, 1), (0, 0), rowspan=3)
                        plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=0.5)
                        for code in range(len(codes)):
                                lines = plt.plot(star_pos_vel_xs[time][code], star_pos_vel_profiles[time][code], color=color_names[code], linestyle=linestyle_names[np.mod(code, len(linestyle_names))], marker=marker_names[code], markeredgecolor='none', linewidth=1.2, alpha=0.8)
                        plt.xlim(0, 14)
                        plt.ylim(0, 270)
                        plt.grid(True)
                        plt.xlabel("$\mathrm{Cylindrical\ Radius\ (kpc)}$", fontsize='large')
                        plt.ylabel("$\mathrm{Rotational\ Velocity\ (km/s)}$", fontsize='large')
                        plt.legend(codes, loc=4, frameon=True, ncol=2, fancybox=True)
                        leg = plt.gca().get_legend()
                        ltext = leg.get_texts()
                        plt.setp(ltext, fontsize='small')

                        plt.subplot2grid((4, 1), (3, 0), rowspan=1)
                        ave_profiles = np.mean(np.array(star_pos_vel_profiles[time]), axis=0)
                        for code in range(len(codes)):
                                lines = plt.plot(star_pos_vel_xs[time][code], (star_pos_vel_profiles[time][code] - ave_profiles)/ave_profiles, color=color_names[code], linestyle='none', marker=marker_names[code], markeredgecolor='none', linewidth=1.2, alpha=0.8)
                        plt.xlim(0, 14)
                        plt.ylim(-0.5, 0.5)
                        plt.grid(True)
                        plt.xlabel("$\mathrm{Cylindrical\ Radius\ (kpc)}$", fontsize='large')
                        plt.ylabel(r"$\mathrm{Residual,}\/(\mathrm{v} - \mathrm{\bar v})/\mathrm{\bar v}$", fontsize='large')
                        if add_mean_fractional_dispersion == 1:
                                mean_fractional_dispersion = (np.std(np.array(star_pos_vel_profiles[time]), axis=0)/ave_profiles)[(mean_dispersion_radius_range[0] < star_pos_vel_xs[time][0]) & (star_pos_vel_xs[time][0] < mean_dispersion_radius_range[1])].mean()
                                mean_fractional_dispersion_in_dex = np.log10(1.+mean_fractional_dispersion)
                                plt.annotate("Mean fractional dispersion for %d < r < %d kpc = %.3f (%.3f dex) " % (mean_dispersion_radius_range[0], mean_dispersion_radius_range[1], mean_fractional_dispersion, mean_fractional_dispersion_in_dex), xy=(0.35, 0.08), size=10, xycoords='axes fraction')
                        plt.savefig("star_pos_vel_%dMyr" % times[time], bbox_inches='tight', pad_inches=0.03, dpi=300)
                        plt.clf()
                if draw_star_pos_vel_PDF >= 3 and time != 0:
                        plt.clf()
#			plt.subplot(111, aspect=0.064)
                        fig = plt.figure(figsize=(8, 10))
                        gridspec.GridSpec(5, 1)
                        plt.subplot2grid((5, 1), (0, 0), rowspan=3)
                        plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=0.5)
                        for code in range(len(codes)):
                                lines = plt.plot(star_pos_disp_xs[time][code], star_pos_disp_profiles[time][code], color=color_names[code], linestyle=linestyle_names[np.mod(code, len(linestyle_names))], marker=marker_names[code], markeredgecolor='none', linewidth=1.2, alpha=0.8)
                        plt.xlim(0, 14)
                        plt.ylim(0, 170)
                        plt.grid(True)
                        plt.xlabel("$\mathrm{Cylindrical\ Radius\ (kpc)}$", fontsize='large')
                        plt.ylabel("$\mathrm{Velocity\ Dispersion\ (km/s)}$", fontsize='large')
                        plt.legend(codes, loc=1, frameon=True, ncol=2, fancybox=True)
                        leg = plt.gca().get_legend()
                        ltext = leg.get_texts()
                        plt.setp(ltext, fontsize='small')

                        plt.subplot2grid((5, 1), (3, 0), rowspan=1)
                        ave_profiles = np.mean(np.array(star_pos_disp_profiles[time]), axis=0)
                        for code in range(len(codes)):
                                lines = plt.plot(star_pos_disp_xs[time][code], (star_pos_disp_profiles[time][code] - ave_profiles)/ave_profiles, color=color_names[code], linestyle='none', marker=marker_names[code], markeredgecolor='none', linewidth=1.2, alpha=0.8)
                        plt.xlim(0, 14)
                        plt.ylim(-1.0, 1.0)
                        plt.grid(True)
                        plt.xlabel("$\mathrm{Cylindrical\ Radius\ (kpc)}$", fontsize='large')
                        plt.ylabel(r"$\mathrm{Residual,}\/(\mathrm{\sigma} - \mathrm{\bar \sigma})/\mathrm{\bar \sigma}$", fontsize='large')
                        if add_mean_fractional_dispersion == 1:
                                mean_fractional_dispersion = (np.std(np.array(star_pos_disp_profiles[time]), axis=0)/ave_profiles)[(mean_dispersion_radius_range[0] < star_pos_disp_xs[time][0]) & (star_pos_disp_xs[time][0] < mean_dispersion_radius_range[1])].mean()
                                mean_fractional_dispersion_in_dex = np.log10(1.+mean_fractional_dispersion)
                                plt.annotate("Mean fractional dispersion for %d < r < %d kpc = %.3f (%.3f dex) " % (mean_dispersion_radius_range[0], mean_dispersion_radius_range[1], mean_fractional_dispersion, mean_fractional_dispersion_in_dex), xy=(0.35, 0.08), size=10, xycoords='axes fraction')
                        plt.savefig("star_pos_disp_%dMyr" % times[time], bbox_inches='tight', pad_inches=0.03, dpi=300)
                        plt.subplot2grid((5, 1), (4, 0), rowspan=1)
                        for code in range(len(codes)):
                                lines = plt.plot(star_pos_disp_xs[time][code], star_pos_disp_vert_profiles[time][code]/star_pos_disp_profiles[time][code], color=color_names[code], linestyle='none', marker=marker_names[code], markeredgecolor='none', linewidth=1.2, alpha=0.8)
                        plt.xlim(0, 14)
                        plt.ylim(0.0, 1.0)
                        plt.grid(True)
                        plt.xlabel("$\mathrm{Cylindrical\ Radius\ (kpc)}$", fontsize='large')
                        plt.ylabel(r"$\mathrm{Ratio,}\/\mathrm{\sigma}_{\perp}/\mathrm{\sigma}$", fontsize='large')
                        plt.savefig("star_pos_disp_%dMyr" % times[time], bbox_inches='tight', pad_inches=0.03, dpi=300)
                        plt.clf()
                if draw_star_pos_vel_PDF == 4 and time != 0:
                        plt.clf()
#			plt.subplot(111, aspect=0.064)
                        fig = plt.figure(figsize=(8, 8))
                        gridspec.GridSpec(4, 1)
                        plt.subplot2grid((4, 1), (0, 0), rowspan=3)
                        plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=0.5)
                        for code in range(len(codes)):
                                lines = plt.plot(star_pos_disp_vert_xs[time][code], star_pos_disp_vert_profiles[time][code], color=color_names[code], linestyle=linestyle_names[np.mod(code, len(linestyle_names))], marker=marker_names[code], markeredgecolor='none', linewidth=1.2, alpha=0.8)
                        plt.xlim(0, 14)
                        plt.ylim(0, 170)
                        plt.grid(True)
                        plt.xlabel("$\mathrm{Cylindrical\ Radius\ (kpc)}$", fontsize='large')
                        plt.ylabel("$\mathrm{Vertical\ Velocity\ Dispersion\ (km/s)}$", fontsize='large')
                        plt.legend(codes, loc=1, frameon=True, ncol=2, fancybox=True)
                        leg = plt.gca().get_legend()
                        ltext = leg.get_texts()
                        plt.setp(ltext, fontsize='small')

                        plt.subplot2grid((4, 1), (3, 0), rowspan=1)
                        ave_profiles = np.mean(np.array(star_pos_disp_vert_profiles[time]), axis=0)
                        for code in range(len(codes)):
                                lines = plt.plot(star_pos_disp_vert_xs[time][code], (star_pos_disp_vert_profiles[time][code] - ave_profiles)/ave_profiles, color=color_names[code], linestyle='none', marker=marker_names[code], markeredgecolor='none', linewidth=1.2, alpha=0.8)
                        plt.xlim(0, 14)
                        plt.ylim(-1.0, 1.0)
                        plt.grid(True)
                        plt.xlabel("$\mathrm{Cylindrical\ Radius\ (kpc)}$", fontsize='large')
                        plt.ylabel(r"$\mathrm{Residual,}\/(\mathrm{\sigma} - \mathrm{\bar \sigma})/\mathrm{\bar \sigma}$", fontsize='large')
                        if add_mean_fractional_dispersion == 1:
                                mean_fractional_dispersion = (np.std(np.array(star_pos_disp_vert_profiles[time]), axis=0)/ave_profiles)[(mean_dispersion_radius_range[0] < star_pos_disp_vert_xs[time][0]) & (star_pos_disp_vert_xs[time][0] < mean_dispersion_radius_range[1])].mean()
                                mean_fractional_dispersion_in_dex = np.log10(1.+mean_fractional_dispersion)
                                plt.annotate("Mean fractional dispersion for %d < r < %d kpc = %.3f (%.3f dex) " % (mean_dispersion_radius_range[0], mean_dispersion_radius_range[1], mean_fractional_dispersion, mean_fractional_dispersion_in_dex), xy=(0.35, 0.08), size=10, xycoords='axes fraction')
                        plt.savefig("star_pos_disp_vert_%dMyr" % times[time], bbox_inches='tight', pad_inches=0.03, dpi=300)
                        plt.clf()
        if draw_rad_height_PDF >= 1:
                fig_rad_height_PDF[time].savefig("rad_height_PDF_%dMyr" % times[time], bbox_inches='tight', pad_inches=0.03, dpi=300)
                if draw_rad_height_PDF == 2 and time != 0:
                        plt.clf()
#			plt.subplot(111, aspect=18)
                        fig = plt.figure(figsize=(8, 8))
                        gridspec.GridSpec(4, 1)
                        plt.subplot2grid((4, 1), (0, 0), rowspan=3)
                        plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=0.5)
                        for code in range(len(codes)):
                                lines = plt.plot(rad_height_xs[time][code], rad_height_profiles[time][code], color=color_names[code], linestyle=linestyle_names[np.mod(code, len(linestyle_names))], marker=marker_names[code], markeredgecolor='none', linewidth=1.2, alpha=0.8)
                        plt.xlim(0, 14)
                        plt.ylim(0, 0.45)
                        plt.grid(True)
                        plt.xlabel("$\mathrm{Cylindrical\ Radius\ (kpc)}$", fontsize='large')
                        plt.ylabel("$\mathrm{Average\ Vertical\ Height\ (kpc)}$", fontsize='large')
                        plt.legend(codes, loc=2, frameon=True, ncol=2, fancybox=True)
                        leg = plt.gca().get_legend()
                        ltext = leg.get_texts()
                        plt.setp(ltext, fontsize='small')

                        plt.subplot2grid((4, 1), (3, 0), rowspan=1)
                        ave_profiles = np.mean(np.array(rad_height_profiles[time]), axis=0)
                        for code in range(len(codes)):
                                lines = plt.plot(rad_height_xs[time][code], (rad_height_profiles[time][code] - ave_profiles)/ave_profiles, color=color_names[code], linestyle='none', marker=marker_names[code], markeredgecolor='none', linewidth=1.2, alpha=0.8)
                        plt.xlim(0, 14)
                        plt.ylim(-1.0, 1.0)
                        plt.grid(True)
                        plt.xlabel("$\mathrm{Cylindrical\ Radius\ (kpc)}$", fontsize='large')
                        plt.ylabel(r"$\mathrm{Residual,}\/(\mathrm{h} - \mathrm{\bar h})/\mathrm{\bar h}$", fontsize='large')
                        if add_mean_fractional_dispersion == 1:
                                mean_fractional_dispersion = (np.std(np.array(rad_height_profiles[time]), axis=0)/ave_profiles)[(mean_dispersion_radius_range[0] < rad_height_xs[time][0]) & (rad_height_xs[time][0] < mean_dispersion_radius_range[1])].mean()
                                mean_fractional_dispersion_in_dex = np.log10(1.+mean_fractional_dispersion)
                                plt.annotate("Mean fractional dispersion for %d < r < %d kpc = %.3f (%.3f dex) " % (mean_dispersion_radius_range[0], mean_dispersion_radius_range[1], mean_fractional_dispersion, mean_fractional_dispersion_in_dex), xy=(0.35, 0.08), size=10, xycoords='axes fraction')
                        plt.savefig("rad_height_%dMyr" % times[time], bbox_inches='tight', pad_inches=0.03, dpi=300)
                        plt.clf()
        if draw_metal_PDF == 1:
                fig_metal_PDF[time].savefig("metal_PDF_%dMyr" % times[time], bbox_inches='tight', pad_inches=0.03, dpi=300)
        if draw_density_DF >= 1:
                plt.clf()
#		plt.subplot(111, aspect=1)
                fig = plt.figure(figsize=(8, 8))
                gridspec.GridSpec(4, 1)
                plt.subplot2grid((4, 1), (0, 0), rowspan=3)
                plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=0.5)
                for code in range(len(codes)):
                        lines = plt.plot(density_DF_xs[time][code], density_DF_profiles[time][code], color=color_names[code], linestyle=linestyle_names[np.mod(code, len(linestyle_names))], marker=marker_names[code], markeredgecolor='none', linewidth=1.2, alpha=0.8)
                plt.semilogx()
                plt.semilogy()
                plt.xlim(1e-28, 1e-21)
                plt.ylim(1e4, 1e9) # accumulation=False
                #plt.ylim(1e5, 2e10) # accumulation=True
                if time != 0 and dataset_num == 2:
                        plt.axvline(x = 1.67e-24*10, color='k', linestyle ='--', linewidth=2, alpha=0.7) # SF threshold density
                plt.grid(True)
                plt.xlabel("$\mathrm{Density\ (g/cm^3)}$", fontsize='large')
                plt.ylabel(r"$\mathrm{Mass,}\/\mathrm{d}M\mathrm{/dlog}\/\mathrm{\rho}\/\mathrm{(M_{\odot})}$", fontsize='large')
                plt.legend(codes, loc=4, frameon=True, ncol=2, fancybox=True)
                leg = plt.gca().get_legend()
                ltext = leg.get_texts()
                plt.setp(ltext, fontsize='small')

                plt.subplot2grid((4, 1), (3, 0), rowspan=1)
                ave_profiles = np.mean(np.array(density_DF_profiles[time]), axis=0)
                for code in range(len(codes)):
                        lines = plt.plot(density_DF_xs[time][code], (density_DF_profiles[time][code] - ave_profiles)/ave_profiles, color=color_names[code], linestyle='none', marker=marker_names[code], markeredgecolor='none', linewidth=1.2, alpha=0.8)
                plt.semilogx()
                plt.xlim(1e-28, 1e-21)
                plt.ylim(-1.0, 3.0)
                plt.grid(True)
                plt.locator_params(axis='y',nbins=4)
                plt.xlabel("$\mathrm{Density\ (g/cm^3)}$", fontsize='large')
                plt.ylabel(r"$\mathrm{Residual,}\/(\mathrm{m} - \mathrm{\bar m})/\mathrm{\bar m}$", fontsize='large')
                if add_mean_fractional_dispersion == 1:
                        mean_fractional_dispersion = (np.std(np.array(density_DF_profiles[time]), axis=0)/ave_profiles)[(mean_dispersion_density_range[0] < density_DF_xs[time][0]) & (density_DF_xs[time][0] < mean_dispersion_density_range[1])].mean()
                        mean_fractional_dispersion_in_dex = np.log10(1.+mean_fractional_dispersion)
                        plt.annotate("Mean fractional dispersion for %.1e < rho < %.1e = %.3f (%.3f dex) " % (mean_dispersion_density_range[0], mean_dispersion_density_range[1], mean_fractional_dispersion, mean_fractional_dispersion_in_dex), xy=(0.27, 0.84), size=10, xycoords='axes fraction')
                plt.savefig("density_DF_%dMyr" % times[time], bbox_inches='tight', pad_inches=0.03, dpi=300)
                plt.clf()
                if draw_density_DF == 2 and time != 0 and dataset_num == 2:
                        fig = plt.figure(figsize=(8, 8))
                        gridspec.GridSpec(16, 1)
                        plt.subplot2grid((16, 1), (0, 0), rowspan=9)
                        plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=0.5)
                        for code in range(len(codes)):
                                lines = plt.plot(density_DF_xs[time][code], (density_DF_profiles[time][code] - density_DF_1st_profiles[time][code])/1e8, color=color_names[code], linestyle=linestyle_names[np.mod(code, len(linestyle_names))], marker=marker_names[code], markeredgecolor='none', linewidth=1.2, alpha=0.8)
                        plt.xscale('log')
                        plt.xlim(1e-28, 1e-21)
                        plt.ylim(-7, 3)
                        plt.axvline(x = 1.67e-24*10, color='k', linestyle ='--', linewidth=2, alpha=0.7) # SF threshold density
                        plt.grid(True)
                        plt.xlabel("$\mathrm{Density\ (g/cm^3)}$", fontsize='large')
                        plt.ylabel(r"$\mathrm{Mass\ Change,}\/\mathrm{\Delta}(\mathrm{d}M\mathrm{/dlog}\/\mathrm{\rho})\/\mathrm{(10^8 M_{\odot})}$", fontsize='large')
                        plt.legend(codes, loc=3, frameon=True, ncol=2, fancybox=True)
                        leg = plt.gca().get_legend()
                        ltext = leg.get_texts()
                        plt.setp(ltext, fontsize='small')

                        plt.subplot2grid((16, 1), (9, 0), rowspan=7)
                        for code in range(len(codes)):
                                lines = plt.plot(density_DF_xs[time][code], (density_DF_profiles[time][code] - density_DF_1st_profiles[time][code])/1e8, color=color_names[code], linestyle='none', marker=marker_names[code], markeredgecolor='none', linewidth=1.2, alpha=0.8)
                        plt.xlabel("$\mathrm{Density\ (g/cm^3)}$", fontsize='large')
                        plt.ylabel(r"$\mathtt{symlog}[\/\mathrm{\Delta}(\mathrm{d}M\mathrm{/dlog}\/\mathrm{\rho})\/\mathrm{(10^8 M_{\odot})}\/]$", fontsize='large')
                        plt.xscale('log')
                        plt.yscale('symlog', linthreshy=0.01)
                        plt.xlim(1e-28, 1e-21)
                        plt.ylim(-10, 10)
                        plt.axvline(x = 1.67e-24*10, color='k', linestyle ='--', linewidth=2, alpha=0.7) # SF threshold density
                        plt.grid(True)
                        plt.savefig("density_DF_change_%dMyr" % times[time], bbox_inches='tight', pad_inches=0.03, dpi=300)
                        plt.clf()
        if draw_radius_DF == 1:
                plt.clf()
#		plt.subplot(111, aspect=1)
                fig = plt.figure(figsize=(8, 6))
                for code in range(len(codes)):
                        lines = plt.plot(radius_DF_xs[time][code], np.add.accumulate(radius_DF_profiles[time][code]), color=color_names[code], linestyle=linestyle_names[np.mod(code, len(linestyle_names))], marker=marker_names[code], markeredgecolor='none', linewidth=1.2, alpha=0.8)
                plt.semilogy()
                plt.xlim(0, 14)
                plt.ylim(1e7, 2e10)
                plt.grid(True)
                plt.xlabel("$\mathrm{Cylindrical\ Radius\ (kpc)}$", fontsize='large')
                plt.ylabel("$\mathrm{Mass\ (M_{\odot})}$", fontsize='large')
                plt.legend(codes, loc=4, frameon=True, ncol=2, fancybox=True)
                leg = plt.gca().get_legend()
                ltext = leg.get_texts()
                plt.setp(ltext, fontsize='small')
                plt.savefig("radius_DF_%dMyr" % times[time], bbox_inches='tight', pad_inches=0.03, dpi=300)

                plt.clf()
#		plt.subplot(111, aspect=1)
                fig = plt.figure(figsize=(8, 8))
                gridspec.GridSpec(4, 1)
                plt.subplot2grid((4, 1), (0, 0), rowspan=3)
                plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=0.5)
                for code in range(len(codes)):
                        temp = []
                        dr = 0.5*(radius_DF_xs[time][code][1] - radius_DF_xs[time][code][0]) # Here we assume that ProfilePlot was made with linearly binned radius_DF_xs
                        for radius in range(len(radius_DF_profiles[time][code])):
                                surface_area = np.pi*(((radius_DF_xs[time][code][radius]+dr)*1e3)**2 - ((radius_DF_xs[time][code][radius]-dr)*1e3)**2)
                                temp.append(radius_DF_profiles[time][code][radius] / surface_area)
                        surface_density[time].append(temp)
                        lines = plt.plot(radius_DF_xs[time][code], surface_density[time][code], color=color_names[code], linestyle=linestyle_names[np.mod(code, len(linestyle_names))], marker=marker_names[code], markeredgecolor='none', linewidth=1.2, alpha=0.8)
                plt.semilogy()
                plt.xlim(0, 14)
                plt.ylim(1e-1, 2e3)
                plt.grid(True)
                plt.xlabel("$\mathrm{Cylindrical\ Radius\ (kpc)}$", fontsize='large')
                plt.ylabel("$\mathrm{Surface\ Density\ (M_{\odot}/pc^2)}$", fontsize='large')
                plt.legend(codes, loc=1, frameon=True, ncol=2, fancybox=True)
                leg = plt.gca().get_legend()
                ltext = leg.get_texts()
                plt.setp(ltext, fontsize='small')

                plt.subplot2grid((4, 1), (3, 0), rowspan=1)
                ave_profiles = np.mean(np.array(surface_density[time]), axis=0)
                for code in range(len(codes)):
                        lines = plt.plot(radius_DF_xs[time][code], (surface_density[time][code] - ave_profiles)/ave_profiles, color=color_names[code], linestyle='none', marker=marker_names[code], markeredgecolor='none', linewidth=1.2, alpha=0.8)
                plt.xlim(0, 14)
                plt.ylim(-1.0, 1.0)
                plt.grid(True)
                plt.xlabel("$\mathrm{Cylindrical\ Radius\ (kpc)}$", fontsize='large')
                plt.ylabel(r"$\mathrm{Residual,}\/(\mathrm{\Sigma} - \mathrm{\bar \Sigma})/\mathrm{\bar \Sigma}$", fontsize='large')
                if add_mean_fractional_dispersion == 1:
                        mean_fractional_dispersion = (np.std(np.array(surface_density[time]), axis=0)/ave_profiles)[(mean_dispersion_radius_range[0] < radius_DF_xs[time][0]) & (radius_DF_xs[time][0] < mean_dispersion_radius_range[1])].mean()
                        mean_fractional_dispersion_in_dex = np.log10(1.+mean_fractional_dispersion)
                        plt.annotate("Mean fractional dispersion for %d < r < %d kpc = %.3f (%.3f dex) " % (mean_dispersion_radius_range[0], mean_dispersion_radius_range[1], mean_fractional_dispersion, mean_fractional_dispersion_in_dex), xy=(0.35, 0.08), size=10, xycoords='axes fraction')
                plt.savefig("gas_surface_density_radial_%dMyr" % times[time], bbox_inches='tight', pad_inches=0.03, dpi=300)
                plt.clf()
        if draw_star_radius_DF >= 1 and time != 0:
                plt.clf()
#		plt.subplot(111, aspect=1)
                fig = plt.figure(figsize=(8, 6))
                for code in range(len(codes)):
                        lines = plt.plot(star_radius_DF_xs[time][code], np.add.accumulate(star_radius_DF_profiles[time][code]), color=color_names[code], linestyle=linestyle_names[np.mod(code, len(linestyle_names))], marker=marker_names[code], markeredgecolor='none', linewidth=1.2, alpha=0.8)
                plt.semilogy()
                plt.xlim(0, 14)
                plt.ylim(1e7, 2e10)
                plt.grid(True)
                plt.xlabel("$\mathrm{Cylindrical\ Radius\ (kpc)}$", fontsize='large')
                plt.ylabel("$\mathrm{Newly\ Formed\ Stellar\ Mass\ (M_{\odot})}$", fontsize='large')
                plt.legend(codes, loc=4, frameon=True, ncol=2, fancybox=True)
                leg = plt.gca().get_legend()
                ltext = leg.get_texts()
                plt.setp(ltext, fontsize='small')
                plt.savefig("star_radius_DF_%dMyr" % times[time], bbox_inches='tight', pad_inches=0.03, dpi=300)

                plt.clf()
#		plt.subplot(111, aspect=1)
                fig = plt.figure(figsize=(8, 8))
                gridspec.GridSpec(4, 1)
                plt.subplot2grid((4, 1), (0, 0), rowspan=3)
                plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=0.5)
                for code in range(len(codes)):
                        temp = []
                        dr = 0.5*(star_radius_DF_xs[time][code][1] - star_radius_DF_xs[time][code][0])
                        for radius in range(len(star_radius_DF_profiles[time][code])):
                                surface_area = np.pi*(((star_radius_DF_xs[time][code][radius]+dr)*1e3)**2 - ((star_radius_DF_xs[time][code][radius]-dr)*1e3)**2)
                                temp.append(star_radius_DF_profiles[time][code][radius] / surface_area)
                        star_surface_density[time].append(temp)
                        lines = plt.plot(star_radius_DF_xs[time][code], star_surface_density[time][code], color=color_names[code], linestyle=linestyle_names[np.mod(code, len(linestyle_names))], marker=marker_names[code], markeredgecolor='none', linewidth=1.2, alpha=0.8)
                plt.semilogy()
                plt.xlim(0, 14)
                plt.ylim(1e-1, 2e3)
                plt.grid(True)
                plt.xlabel("$\mathrm{Cylindrical\ Radius\ (kpc)}$", fontsize='large')
                plt.ylabel("$\mathrm{Newly\ Formed\ Stellar\ Surface\ Density\ (M_{\odot}/pc^2)}$", fontsize='large')
                plt.legend(codes, loc=1, frameon=True, ncol=2, fancybox=True)
                leg = plt.gca().get_legend()
                ltext = leg.get_texts()
                plt.setp(ltext, fontsize='small')

                plt.subplot2grid((4, 1), (3, 0), rowspan=1)
                ave_profiles = np.mean(np.array(star_surface_density[time]), axis=0)
                for code in range(len(codes)):
                        lines = plt.plot(star_radius_DF_xs[time][code], (star_surface_density[time][code] - ave_profiles)/ave_profiles, color=color_names[code], linestyle='none', marker=marker_names[code], markeredgecolor='none', linewidth=1.2, alpha=0.8)
                plt.xlim(0, 14)
                plt.ylim(-1.0, 3.0)
                plt.grid(True)
                plt.locator_params(axis='y',nbins=4)
                plt.xlabel("$\mathrm{Cylindrical\ Radius\ (kpc)}$", fontsize='large')
                plt.ylabel(r"$\mathrm{Residual,}\/(\mathrm{\Sigma} - \mathrm{\bar \Sigma})/\mathrm{\bar \Sigma}$", fontsize='large')
                if add_mean_fractional_dispersion == 1:
                        mean_fractional_dispersion = (np.std(np.array(star_surface_density[time]), axis=0)/ave_profiles)[(mean_dispersion_radius_range[0] < star_radius_DF_xs[time][0]) & (star_radius_DF_xs[time][0] < mean_dispersion_radius_range[1])].mean()
                        mean_fractional_dispersion_in_dex = np.log10(1.+mean_fractional_dispersion)
                        plt.annotate("Mean fractional dispersion for %d < r < %d kpc = %.3f (%.3f dex) " % (mean_dispersion_radius_range[0], mean_dispersion_radius_range[1], mean_fractional_dispersion, mean_fractional_dispersion_in_dex), xy=(0.02, 0.84), size=10, xycoords='axes fraction')
                plt.savefig("star_surface_density_radial_%dMyr" % times[time], bbox_inches='tight', pad_inches=0.03, dpi=300)
                plt.clf()

                if draw_star_radius_DF == 2 and time != 0:
#			plt.subplot(111, aspect=1)
                        fig = plt.figure(figsize=(8, 8))
                        gridspec.GridSpec(4, 1)
                        plt.subplot2grid((4, 1), (0, 0), rowspan=3)
                        plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=0.5)
                        for code in range(len(codes)):
                                temp = []
                                dr = 0.5*(sfr_radius_DF_xs[time][code][1] - sfr_radius_DF_xs[time][code][0])
                                for radius in range(len(sfr_radius_DF_profiles[time][code])):
                                        surface_area = np.pi*((sfr_radius_DF_xs[time][code][radius]+dr)**2 - (sfr_radius_DF_xs[time][code][radius]-dr)**2)
                                        temp.append(sfr_radius_DF_profiles[time][code][radius] / surface_area)
                                sfr_surface_density[time].append(temp)
                                lines = plt.plot(sfr_radius_DF_xs[time][code], sfr_surface_density[time][code], color=color_names[code], linestyle=linestyle_names[np.mod(code, len(linestyle_names))], marker=marker_names[code], markeredgecolor='none', linewidth=1.2, alpha=0.8)
                        plt.semilogy()
                        plt.xlim(0, 14)
                        plt.ylim(1e-4, 1e1)
                        plt.grid(True)
                        plt.xlabel("$\mathrm{Cylindrical\ Radius\ (kpc)}$", fontsize='large')
                        plt.ylabel("$\mathrm{Star\ Formation\ Rate\ Surface\ Density\ (M_{\odot}/yr/kpc^2)}$", fontsize='large')
                        plt.legend(codes, loc=1, frameon=True, ncol=2, fancybox=True)
                        leg = plt.gca().get_legend()
                        ltext = leg.get_texts()
                        plt.setp(ltext, fontsize='small')

                        plt.subplot2grid((4, 1), (3, 0), rowspan=1)
                        ave_profiles = np.mean(np.array(sfr_surface_density[time]), axis=0)
                        for code in range(len(codes)):
                                lines = plt.plot(sfr_radius_DF_xs[time][code], (sfr_surface_density[time][code] - ave_profiles)/ave_profiles, color=color_names[code], linestyle='none', marker=marker_names[code], markeredgecolor='none', linewidth=1.2, alpha=0.8)
                        plt.xlim(0, 14)
                        plt.ylim(-1.0, 3.0)
                        plt.grid(True)
                        plt.locator_params(axis='y',nbins=4)
                        plt.xlabel("$\mathrm{Cylindrical\ Radius\ (kpc)}$", fontsize='large')
                        plt.ylabel(r"$\mathrm{Residual,}\/(\mathrm{\Sigma} - \mathrm{\bar \Sigma})/\mathrm{\bar \Sigma}$", fontsize='large')
                        # if add_mean_fractional_dispersion == 1:
                        #       mean_fractional_dispersion = (np.std(np.array(sfr_surface_density[time]), axis=0)/ave_profiles)[(mean_dispersion_radius_range[0] < sfr_radius_DF_xs[time][0]) & (sfr_radius_DF_xs[time][0] < mean_dispersion_radius_range[1])].mean()
                        #       mean_fractional_dispersion_in_dex = np.log10(1.+mean_fractional_dispersion)
                        #       plt.annotate("Mean fractional dispersion for %d < r < %d kpc = %.3f (%.3f dex) " % (mean_dispersion_radius_range[0], mean_dispersion_radius_range[1], mean_fractional_dispersion, mean_fractional_dispersion_in_dex), xy=(0.02, 0.84), size=10, xycoords='axes fraction')
                        plt.savefig("sfr_surface_density_radial_%dMyr" % times[time], bbox_inches='tight', pad_inches=0.03, dpi=300)
                        plt.clf()

                        # Draw K-S plot; below assumes that surface_density (or radius_DF_profiles) and sfr_surface_density (or sfr_radius_DF_profiles) have the identical size (n_bins in ProfilePlot)
#			plt.subplot(111)
                        fig = plt.figure(figsize=(8, 7))
                        KS_fit_t1 = []
                        KS_fit_t2 = []
                        Bigiel_surface_density = []
                        Bigiel_sfr_surface_density = []
                        for code in range(len(codes)):
                                # Remove bins where SFR surface density is zero
                                KS_x = np.array(surface_density[time][code])
                                KS_x[np.where(np.array(sfr_surface_density[time][code]) < 1e-10)] = 0
                                KS_x = np.log10(KS_x)
                                KS_y = np.log10(np.array(sfr_surface_density[time][code]))
                                KS_x = KS_x[~np.isinf(KS_x)]
                                KS_y = KS_y[~np.isinf(KS_y)]
                                t1, t2 = np.polyfit(KS_x, KS_y, 1)
                                KS_fit_t1.append(t1)
                                KS_fit_t2.append(t2)
                                plt.scatter(KS_x, KS_y, color=color_names[code], edgecolor=color_names[code], s=30, linewidth=0.7, marker=marker_names[code], alpha=0.8)
                        plt.xlim(0, 3)
                        plt.ylim(-4, 1)
                        plt.grid(True)
                        plt.xlabel("$\mathrm{log[Gas\ Surface\ Density\ (M_{\odot}/pc^2)]}$")
                        plt.ylabel("$\mathrm{log[Star\ Formation\ Rate\ Surface\ Density\ (M_{\odot}/yr/kpc^2)]}$")
                        plt.legend(codes, loc=2, frameon=True, ncol=2, fancybox=True)
                        leg = plt.gca().get_legend()
                        ltext = leg.get_texts()
                        plt.setp(ltext, fontsize='small')
                        t = np.arange(-2, 5, 0.01)
                        for line in import_text("./Bigieletal2008_Fig8_Contour.txt", " "):
                                Bigiel_surface_density.append(float(line[0]) + np.log10(1.36)) # for the factor 1.36, see Section 2.3.1 of Bigiel et al. 2008
                                Bigiel_sfr_surface_density.append(float(line[1]))
                        plt.fill(Bigiel_surface_density, Bigiel_sfr_surface_density, fill=True, color='b', alpha = 0.1, hatch='\\') # contour by Bigiel et al. 2008 (Fig 8), or Feldmann et al. 2012 (Fig 1)
#			plt.plot(Bigiel_surface_density, Bigiel_sfr_surface_density, 'b-', linewidth = 0.7, alpha = 0.4)
                        plt.plot(t, 1.37*t - 3.78, 'k--', linewidth = 2, alpha = 0.7) # observational fit by Kennicutt et al. 2007
                        plt.savefig("K-S_%dMyr" % times[time], bbox_inches='tight', pad_inches=0.03, dpi=300)
                        for code in range(len(codes)):
                                plt.plot(t, np.polyval([KS_fit_t1[code], KS_fit_t2[code]], t), color=color_names[code], linestyle=linestyle_names[np.mod(code, len(linestyle_names))], alpha = 0.7) # linear fits
                        plt.savefig("K-S_with_fits_%dMyr" % times[time], bbox_inches='tight', pad_inches=0.03, dpi=300)
                        plt.clf()
        if draw_height_DF == 1:
                plt.clf()
#		plt.subplot(111, aspect=0.5)
                fig = plt.figure(figsize=(8, 6))
                for code in range(len(codes)):
                        lines = plt.plot(height_DF_xs[time][code], np.add.accumulate(height_DF_profiles[time][code]), color=color_names[code], linestyle=linestyle_names[np.mod(code, len(linestyle_names))], marker=marker_names[code], markeredgecolor='none', linewidth=1.2, alpha=0.8)
                plt.semilogy()
                plt.xlim(0, 1.4)
                plt.ylim(1e9, 2e10)
                plt.grid(True)
                plt.xlabel("$\mathrm{Vertical\ Height\ (kpc)}$", fontsize='large')
                plt.ylabel("$\mathrm{Mass\ (M_{\odot})}$", fontsize='large')
                plt.legend(codes, loc=4, frameon=True, ncol=2, fancybox=True)
                leg = plt.gca().get_legend()
                ltext = leg.get_texts()
                plt.setp(ltext, fontsize='small')
                plt.savefig("height_DF_%dMyr" % times[time], bbox_inches='tight', pad_inches=0.03, dpi=300)

                plt.clf()
#		plt.subplot(111, aspect=0.8)
                fig = plt.figure(figsize=(8, 8))
                gridspec.GridSpec(4, 1)
                plt.subplot2grid((4, 1), (0, 0), rowspan=3)
                plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=0.5)
                for code in range(len(codes)):
                        temp = []
                        dh = height_DF_xs[time][code][1] - height_DF_xs[time][code][0]
                        for height in range(len(height_DF_profiles[time][code])):
                                surface_area = 2 * dh*1e3 * figure_width*1e3 # surface_area = 2*d(height)*figure_width in pc^2
                                temp.append(height_DF_profiles[time][code][height] / surface_area)
                        height_surface_density[time].append(temp)
                        lines = plt.plot(height_DF_xs[time][code], height_surface_density[time][code], color=color_names[code], linestyle=linestyle_names[np.mod(code, len(linestyle_names))], marker=marker_names[code], markeredgecolor='none', linewidth=1.2, alpha=0.8)
                plt.semilogy()
                plt.xlim(0, 1.4)
                plt.ylim(1e-1, 2e3)
                plt.grid(True)
                plt.xlabel("$\mathrm{Vertical\ Height\ (kpc)}$", fontsize='large')
                plt.ylabel("$\mathrm{Surface\ Density\ (M_{\odot}/pc^2)}$", fontsize='large')
                plt.legend(codes, loc=1, frameon=True, ncol=2, fancybox=True)
                leg = plt.gca().get_legend()
                ltext = leg.get_texts()
                plt.setp(ltext, fontsize='small')

                plt.subplot2grid((4, 1), (3, 0), rowspan=1)
                ave_profiles = np.mean(np.array(height_surface_density[time]), axis=0)
                for code in range(len(codes)):
                        lines = plt.plot(height_DF_xs[time][code], (height_surface_density[time][code] - ave_profiles)/ave_profiles, color=color_names[code], linestyle='none', marker=marker_names[code], markeredgecolor='none', linewidth=1.2, alpha=0.8)
                plt.xlim(0, 1.4)
                plt.ylim(-1.0, 3.0)
                plt.grid(True)
                plt.locator_params(axis='y',nbins=4)
                plt.xlabel("$\mathrm{Vertical\ Height\ (kpc)}$", fontsize='large')
                plt.ylabel(r"$\mathrm{Residual,}\/(\mathrm{\Sigma} - \mathrm{\bar \Sigma})/\mathrm{\bar \Sigma}$", fontsize='large')
                if add_mean_fractional_dispersion == 1:
                        mean_fractional_dispersion = (np.std(np.array(height_surface_density[time]), axis=0)/ave_profiles)[(mean_dispersion_height_range[0] < height_DF_xs[time][0]) & (height_DF_xs[time][0] < mean_dispersion_height_range[1])].mean()
                        mean_fractional_dispersion_in_dex = np.log10(1.+mean_fractional_dispersion)
                        plt.annotate("Mean fractional dispersion for %.1f < z < %.1f kpc = %.3f (%.3f dex) " % (mean_dispersion_height_range[0], mean_dispersion_height_range[1], mean_fractional_dispersion, mean_fractional_dispersion_in_dex), xy=(0.02, 0.84), size=10, xycoords='axes fraction')
                plt.savefig("gas_surface_density_vertical_%dMyr" % times[time], bbox_inches='tight', pad_inches=0.03, dpi=300)
                plt.clf()
        if draw_SFR >= 1 and time != 0:
                plt.clf()
#		plt.subplot(111)#, aspect=1e-7)
                fig = plt.figure(figsize=(8, 6))
                for code in range(len(codes)):
                        lines = plt.plot(sfr_ts[time][code], sfr_cum_masses[time][code], color=color_names[code], linestyle=linestyle_names[np.mod(code, len(linestyle_names))], marker=marker_names[code], markeredgecolor='none', linewidth=1.2, alpha=0.8)
                plt.xlim(0, times[time])
                plt.ylim(0, 2.5e9)
                plt.grid(True)
                plt.xlabel("$\mathrm{Time\ (Myr)}$", fontsize='large')
                plt.ylabel("$\mathrm{Stellar\ Mass\ (M_{\odot})}$", fontsize='large')
                plt.legend(codes, loc=2, frameon=True, ncol=2, fancybox=True)
                leg = plt.gca().get_legend()
                ltext = leg.get_texts()
                plt.setp(ltext, fontsize='small')
                plt.savefig("Stellar_mass_evolution_%dMyr" % times[time], bbox_inches='tight', pad_inches=0.03, dpi=300)
                plt.clf()
#		plt.subplot(111, aspect=55)
                fig = plt.figure(figsize=(8, 8))
                gridspec.GridSpec(4, 1)
                plt.subplot2grid((4, 1), (0, 0), rowspan=3)
                plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=0.5)
                for code in range(len(codes)):
                        lines = plt.plot(sfr_ts[time][code], sfr_sfrs[time][code], color=color_names[code], linestyle=linestyle_names[np.mod(code, len(linestyle_names))], marker=marker_names[code], markeredgecolor='none', linewidth=1.2, alpha=0.8)
                plt.xlim(0, times[time])
                plt.ylim(0, 8)
                plt.grid(True)
                plt.xlabel("$\mathrm{Time\ (Myr)}$", fontsize='large')
                plt.ylabel("$\mathrm{Star\ Formation\ Rate\ (M_{\odot}/yr)}$", fontsize='large')
                if draw_SFR == 2:
                        pf_sfr_ref = load_dataset(2, 1, 0, ['GIZMO'],
                                                  [['dummy', '/lustre/ki/pfs/mornkr/080112_CHaRGe/pfs-hyades/AGORA-DISK-repository/Grackle+SF+ThermalFbck/GIZMO/AGORA_disk_second_ps1.5/snapshot_temp_100_last']])
                        pf_sfr_ref.unit_registry.add('pccm', pf_sfr_ref.unit_registry.lut['pc'][0], length, "\\rm{pc}/(1+z)")
                        pf_sfr_ref.hubble_constant = 0.71; pf_sfr_ref.omega_lambda = 0.73; pf_sfr_ref.omega_matter = 0.27; pf_sfr_ref.omega_curvature = 0.0
                        sp_sfr_ref = pf_sfr_ref.sphere(center, (0.5*figure_width, "kpc"))
                        draw_SFR_mass_ref = sp_sfr_ref[(PartType_Star_to_use, "particle_mass")].in_units('Msun')
                        draw_SFR_ct_ref   = sp_sfr_ref[(PartType_Star_to_use, FormationTimeType_to_use)].in_units('Myr')
                        sfr_ref = StarFormationRate(pf_sfr_ref, star_mass = draw_SFR_mass_ref, star_creation_time = draw_SFR_ct_ref,
                                                    volume = sp_sfr_ref.volume(), bins = 25)
                        lines = plt.plot(sfr_ref.time.in_units('Myr'), sfr_ref.Msol_yr, color='k', linestyle='--', marker='*', markeredgecolor='none', linewidth=1.2, alpha=0.8)
                plt.legend(codes+['GIZMO-ps2'], loc=2, frameon=True, ncol=2, fancybox=True)
                leg = plt.gca().get_legend()
                ltext = leg.get_texts()
                plt.setp(ltext, fontsize='small')
                plt.subplot2grid((4, 1), (3, 0), rowspan=1)
                ave_profiles = np.mean(np.array(sfr_sfrs[time]), axis=0)
                for code in range(len(codes)):
                        lines = plt.plot(sfr_ts[time][code], (sfr_sfrs[time][code].d - ave_profiles)/ave_profiles, color=color_names[code], linestyle='none', marker=marker_names[code], markeredgecolor='none', linewidth=1.2, alpha=0.8)
                plt.xlim(0, times[time])
                plt.ylim(-1.0, 1.0)
                plt.grid(True)
                plt.xlabel("$\mathrm{Time\ (Myr)}$", fontsize='large')
                plt.ylabel(r"$\mathrm{Residual,}\/(\mathrm{\dot{\rho}_\star} - \mathrm{\bar{\dot{\rho}_\star}})/\mathrm{\bar{\dot{\rho}_\star}}$", fontsize='large')
                if add_mean_fractional_dispersion == 1:
                        mean_fractional_dispersion = (np.std(np.array(sfr_sfrs[time]), axis=0)/ave_profiles)[(mean_dispersion_time_range[0] < sfr_ts[time][0]) & (sfr_ts[time][0] < mean_dispersion_time_range[1])].mean()
                        mean_fractional_dispersion_in_dex = np.log10(1.+mean_fractional_dispersion)
                        plt.annotate("Mean fractional dispersion for %d < t < %d Myr = %.3f (%.3f dex) " % (mean_dispersion_time_range[0], mean_dispersion_time_range[1], mean_fractional_dispersion, mean_fractional_dispersion_in_dex), xy=(0.35, 0.08), size=10, xycoords='axes fraction')
                plt.savefig("SFR_%dMyr" % times[time], bbox_inches='tight', pad_inches=0.03, dpi=300)
                plt.clf()
        if draw_cut_through == 1:
                plt.clf()
#		plt.subplot(111, aspect=0.5)
                fig = plt.figure(figsize=(8, 6))
                for code in range(len(codes)):
                        lines = plt.plot(cut_through_zs[time][code], cut_through_zvalues[time][code], color=color_names[code], linestyle=linestyle_names[np.mod(code, len(linestyle_names))], marker=marker_names[code], markeredgecolor='none', linewidth=1.2, alpha=0.8)
                plt.semilogy()
                plt.xlim(-1.4, 1.4)
                plt.ylim(1e-26, 1e-22)
                plt.grid(True)
                plt.xlabel("$\mathrm{Vertical\ Height\ (kpc)}$", fontsize='large')
                plt.ylabel("$\mathrm{Density\ (g/cm^3)}$", fontsize='large')
                plt.legend(codes, loc=1, frameon=True, ncol=2, fancybox=True)
                leg = plt.gca().get_legend()
                ltext = leg.get_texts()
                plt.setp(ltext, fontsize='small')
                z = np.arange(-1.4, 1.4, 0.05)
                plt.plot(z, rho_agora_disk(0, np.abs(z)), linestyle="--", linewidth=2, color='k', alpha=0.7)
                plt.savefig("cut_through_z_%dMyr" % times[time], bbox_inches='tight', pad_inches=0.03, dpi=300)
                plt.clf()
#		plt.subplot(111, aspect=0.5)
                fig = plt.figure(figsize=(8, 6))
                for code in range(len(codes)):
                        lines = plt.plot(cut_through_xs[time][code], cut_through_xvalues[time][code], color=color_names[code], linestyle=linestyle_names[np.mod(code, len(linestyle_names))], marker=marker_names[code], markeredgecolor='none', linewidth=1.2, alpha=0.8)
                plt.semilogy()
                plt.xlim(-14, 14)
                plt.ylim(1e-26, 1e-22)
                plt.grid(True)
                plt.xlabel("$\mathrm{Cylindrical\ Radius\ (kpc)}$", fontsize='large')
                plt.ylabel("$\mathrm{Density\ (g/cm^3)}$", fontsize='large')
                plt.legend(codes, loc=1, frameon=True, ncol=2, fancybox=True)
                leg = plt.gca().get_legend()
                ltext = leg.get_texts()
                plt.setp(ltext, fontsize='small')
                x = np.arange(-14, 14, 0.05)
                plt.plot(x, rho_agora_disk(np.abs(x), 0), linestyle="--", linewidth=2, color='k', alpha=0.7)
                plt.savefig("cut_through_x_%dMyr" % times[time], bbox_inches='tight', pad_inches=0.03, dpi=300)
                plt.clf()
