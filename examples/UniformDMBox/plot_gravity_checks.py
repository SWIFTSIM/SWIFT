#!/usr/bin/env python

import sys
import glob
import re
import numpy as np
import matplotlib.pyplot as plt

params = {'axes.labelsize': 14,
'axes.titlesize': 18,
'font.size': 12,
'legend.fontsize': 12,
'xtick.labelsize': 14,
'ytick.labelsize': 14,
'text.usetex': True,
'figure.figsize': (10, 10),
'figure.subplot.left'    : 0.06,
'figure.subplot.right'   : 0.99  ,
'figure.subplot.bottom'  : 0.06  ,
'figure.subplot.top'     : 0.985  ,
'figure.subplot.wspace'  : 0.14  ,
'figure.subplot.hspace'  : 0.14  ,
'lines.markersize' : 6,
'lines.linewidth' : 3.,
'text.latex.unicode': True
}
plt.rcParams.update(params)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Times']})

min_error = 1e-6
max_error = 1e-1
num_bins = 51
bin_edges = np.linspace(np.log10(min_error), np.log10(max_error), num_bins + 1)
bin_size = (np.log10(max_error) - np.log10(min_error)) / num_bins
bins = 0.5*(bin_edges[1:] + bin_edges[:-1])
bin_edges = 10**bin_edges
bins = 10**bins

# Time-step to plot
step = int(sys.argv[1])

# Find the files for the different expansion orders
order_list = glob.glob("gravity_checks_step%d_order*.dat"%step)
num_order = len(order_list)

# Get the multipole orders
order = np.zeros(num_order)
for i in range(num_order):
    order[i] = int(order_list[i][26])
    
# Start the plot
plt.figure()

# Get the Gadget-2 data if existing
gadget2_file_list = glob.glob("forcetest_gadget2.txt")
if len(gadget2_file_list) != 0:

    data = np.loadtxt(gadget2_file_list[0])
    pos = data[:,1:4]
    a_exact = data[:,4:7]
    a_grav = data[:, 7:10]

    # Compute the error norm
    diff = a_exact - a_grav

    norm_diff = np.sqrt(diff[:,0]**2 + diff[:,1]**2 + diff[:,2]**2)
    norm_a = np.sqrt(a_exact[:,0]**2 + a_exact[:,1]**2 + a_exact[:,2]**2)

    norm_error = norm_diff / norm_a
    error_x = diff[:,0] / norm_a
    error_y = diff[:,1] / norm_a
    error_z = diff[:,2] / norm_a
    
    # Bin the error
    norm_error_hist,_ = np.histogram(norm_error, bins=bin_edges, density=False) / (np.size(norm_error) * bin_size)
    error_x_hist,_ = np.histogram(error_x, bins=bin_edges, density=False) / (np.size(norm_error) * bin_size)
    error_y_hist,_ = np.histogram(error_y, bins=bin_edges, density=False) / (np.size(norm_error) * bin_size)
    error_z_hist,_ = np.histogram(error_z, bins=bin_edges, density=False) / (np.size(norm_error) * bin_size)

    plt.subplot(221)    
    plt.semilogx(bins, norm_error_hist, 'k--', label="Gadget-2")
    plt.subplot(222)    
    plt.semilogx(bins, error_x_hist, 'k--', label="Gadget-2")
    plt.subplot(223)    
    plt.semilogx(bins, error_y_hist, 'k--', label="Gadget-2")
    plt.subplot(224)    
    plt.semilogx(bins, error_z_hist, 'k--', label="Gadget-2")


# Plot the different histograms
for i in range(num_order):
    data = np.loadtxt(order_list[i])
    pos = data[:,1:4]
    a_exact = data[:,4:7]
    a_grav = data[:, 7:10]

    # Compute the error norm
    diff = a_exact - a_grav

    norm_diff = np.sqrt(diff[:,0]**2 + diff[:,1]**2 + diff[:,2]**2)
    norm_a = np.sqrt(a_exact[:,0]**2 + a_exact[:,1]**2 + a_exact[:,2]**2)

    norm_error = norm_diff / norm_a
    error_x = diff[:,0] / norm_a
    error_y = diff[:,1] / norm_a
    error_z = diff[:,2] / norm_a
    
    # Bin the error
    norm_error_hist,_ = np.histogram(norm_error, bins=bin_edges, density=False) / (np.size(norm_error) * bin_size)
    error_x_hist,_ = np.histogram(error_x, bins=bin_edges, density=False) / (np.size(norm_error) * bin_size)
    error_y_hist,_ = np.histogram(error_y, bins=bin_edges, density=False) / (np.size(norm_error) * bin_size)
    error_z_hist,_ = np.histogram(error_z, bins=bin_edges, density=False) / (np.size(norm_error) * bin_size)

    plt.subplot(221)
    plt.semilogx(bins, norm_error_hist, label="SWIFT Multipoles order %d"%order[i])
    plt.subplot(222)    
    plt.semilogx(bins, error_x_hist, label="SWIFT Multipoles order %d"%order[i])
    plt.subplot(223)    
    plt.semilogx(bins, error_y_hist, label="SWIFT Multipoles order %d"%order[i])
    plt.subplot(224)    
    plt.semilogx(bins, error_z_hist, label="SWIFT Multipoles order %d"%order[i])

    
plt.subplot(221)
plt.xlabel("$|\delta \overrightarrow{a}|/|\overrightarrow{a}_{exact}|$")
plt.ylabel("Density")
plt.xlim(min_error, 2*max_error)
plt.ylim(0,3)
plt.legend(loc="upper left")
plt.subplot(222)    
plt.xlabel("$\delta a_x/|\overrightarrow{a}_{exact}|$")
plt.ylabel("Density")
plt.xlim(min_error, 2*max_error)
plt.ylim(0,1)
plt.legend(loc="upper left")
plt.subplot(223)    
plt.xlabel("$\delta a_y/|\overrightarrow{a}_{exact}|$")
plt.ylabel("Density")
plt.xlim(min_error, 2*max_error)
plt.ylim(0,1)
plt.legend(loc="upper left")
plt.subplot(224)    
plt.xlabel("$\delta a_z/|\overrightarrow{a}_{exact}|$")
plt.ylabel("Density")
plt.xlim(min_error, 2*max_error)
plt.ylim(0,1)
plt.legend(loc="upper left")



plt.savefig("gravity_checks_step%d.png"%step)
