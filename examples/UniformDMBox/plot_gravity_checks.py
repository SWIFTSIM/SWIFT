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

# Construct the bins
bin_edges = np.linspace(np.log10(min_error), np.log10(max_error), num_bins + 1)
bin_size = (np.log10(max_error) - np.log10(min_error)) / num_bins
bins = 0.5*(bin_edges[1:] + bin_edges[:-1])
bin_edges = 10**bin_edges
bins = 10**bins

# Colours
cols = ['b', 'g', 'r', 'm']

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

    gadget2_data = np.loadtxt(gadget2_file_list[0])
    gadget2_ids = gadget2_data[:,0]
    gadget2_pos = gadget2_data[:,1:4]
    gadget2_a_exact = gadget2_data[:,4:7]
    gadget2_a_grav = gadget2_data[:, 7:10]

    # Sort stuff
    sort_index = np.argsort(gadget2_ids)
    gadget2_ids = gadget2_ids[sort_index]
    gadget2_pos = gadget2_pos[sort_index, :]
    gadget2_a_exact = gadget2_a_exact[sort_index, :]
    gadget2_a_grav = gadget2_a_grav[sort_index, :]
    
    # Compute the error norm
    diff = gadget2_a_exact - gadget2_a_grav

    norm_diff = np.sqrt(diff[:,0]**2 + diff[:,1]**2 + diff[:,2]**2)
    norm_a = np.sqrt(gadget2_a_exact[:,0]**2 + gadget2_a_exact[:,1]**2 + gadget2_a_exact[:,2]**2)

    norm_error = norm_diff / norm_a
    error_x = abs(diff[:,0]) / norm_a
    error_y = abs(diff[:,1]) / norm_a
    error_z = abs(diff[:,2]) / norm_a
    
    # Bin the error
    norm_error_hist,_ = np.histogram(norm_error, bins=bin_edges, density=False) / (np.size(norm_error) * bin_size)
    error_x_hist,_ = np.histogram(error_x, bins=bin_edges, density=False) / (np.size(norm_error) * bin_size)
    error_y_hist,_ = np.histogram(error_y, bins=bin_edges, density=False) / (np.size(norm_error) * bin_size)
    error_z_hist,_ = np.histogram(error_z, bins=bin_edges, density=False) / (np.size(norm_error) * bin_size)
    
    norm_median = np.median(norm_error)
    median_x = np.median(error_x)
    median_y = np.median(error_y)
    median_z = np.median(error_z)

    norm_per95 = np.percentile(norm_error,95)
    per95_x = np.percentile(error_x,95)
    per95_y = np.percentile(error_y,95)
    per95_z = np.percentile(error_z,95)

    plt.subplot(221)    
    plt.semilogx(bins, norm_error_hist, 'k--', label="Gadget-2")
    plt.plot([norm_median, norm_median], [2.7, 3], 'k-', lw=1)
    plt.plot([norm_per95, norm_per95], [2.7, 3], 'k:', lw=1)
    plt.subplot(222)
    plt.semilogx(bins, error_x_hist, 'k--', label="Gadget-2")
    plt.plot([median_x, median_x], [1.8, 2], 'k-', lw=1)
    plt.plot([per95_x, per95_x], [1.8, 2], 'k:', lw=1)
    plt.subplot(223)    
    plt.semilogx(bins, error_y_hist, 'k--', label="Gadget-2")
    plt.plot([median_y, median_y], [1.8, 2], 'k-', lw=1)
    plt.plot([per95_y, per95_y], [1.8, 2], 'k:', lw=1)
    plt.subplot(224)    
    plt.semilogx(bins, error_z_hist, 'k--', label="Gadget-2")
    plt.plot([median_z, median_z], [1.8, 2], 'k-', lw=1)
    plt.plot([per95_z, per95_z], [1.8, 2], 'k:', lw=1)


# Plot the different histograms
for i in range(num_order-1, -1, -1):
    data = np.loadtxt(order_list[i])
    ids = data[:,0]
    pos = data[:,1:4]
    a_exact = data[:,4:7]
    a_grav = data[:, 7:10]

    # Sort stuff
    sort_index = np.argsort(ids)
    ids = ids[sort_index]
    pos = pos[sort_index, :]
    a_exact = a_exact[sort_index, :]
    a_grav = a_grav[sort_index, :]

    # Cross-checks
    if not np.array_equal(ids, gadget2_ids):
        print "Comparing different IDs !"

    if not np.array_equal(pos, gadget2_pos):
        print "Comparing different positions ! max difference:", np.max(pos - gadget2_pos)

    if not np.array_equal(a_exact, gadget2_a_exact):
        print "Comparing different exact accelerations ! max difference:", np.max(a_exact - gadget2_a_exact)
        
        
    # Compute the error norm
    diff = a_exact - a_grav

    norm_diff = np.sqrt(diff[:,0]**2 + diff[:,1]**2 + diff[:,2]**2)
    norm_a = np.sqrt(a_exact[:,0]**2 + a_exact[:,1]**2 + a_exact[:,2]**2)

    norm_error = norm_diff / norm_a
    error_x = abs(diff[:,0]) / norm_a
    error_y = abs(diff[:,1]) / norm_a
    error_z = abs(diff[:,2]) / norm_a
    
    # Bin the error
    norm_error_hist,_ = np.histogram(norm_error, bins=bin_edges, density=False) / (np.size(norm_error) * bin_size)
    error_x_hist,_ = np.histogram(error_x, bins=bin_edges, density=False) / (np.size(norm_error) * bin_size)
    error_y_hist,_ = np.histogram(error_y, bins=bin_edges, density=False) / (np.size(norm_error) * bin_size)
    error_z_hist,_ = np.histogram(error_z, bins=bin_edges, density=False) / (np.size(norm_error) * bin_size)

    norm_median = np.median(norm_error)
    median_x = np.median(error_x)
    median_y = np.median(error_y)
    median_z = np.median(error_z)

    norm_per95 = np.percentile(norm_error,95)
    per95_x = np.percentile(error_x,95)
    per95_y = np.percentile(error_y,95)
    per95_z = np.percentile(error_z,95)
    
    plt.subplot(221)
    plt.semilogx(bins, norm_error_hist, color=cols[i],label="SWIFT m-poles order %d"%order[i])
    plt.plot([norm_median, norm_median], [2.7, 3],'-', color=cols[i], lw=1)
    plt.plot([norm_per95, norm_per95], [2.7, 3],':', color=cols[i], lw=1)
    plt.subplot(222)    
    plt.semilogx(bins, error_x_hist, color=cols[i],label="SWIFT m-poles order %d"%order[i])
    plt.plot([median_x, median_x], [1.8, 2],'-', color=cols[i], lw=1)
    plt.plot([per95_x, per95_x], [1.8, 2],':', color=cols[i], lw=1)
    plt.subplot(223)    
    plt.semilogx(bins, error_y_hist, color=cols[i],label="SWIFT m-poles order %d"%order[i])
    plt.plot([median_y, median_y], [1.8, 2],'-', color=cols[i], lw=1)
    plt.plot([per95_y, per95_y], [1.8, 2],':', color=cols[i], lw=1)
    plt.subplot(224)    
    plt.semilogx(bins, error_z_hist, color=cols[i],label="SWIFT m-poles order %d"%order[i])
    plt.plot([median_z, median_z], [1.8, 2],'-', color=cols[i], lw=1)
    plt.plot([per95_z, per95_z], [1.8, 2],':', color=cols[i], lw=1)

    
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
plt.ylim(0,2)
#plt.legend(loc="center left")
plt.subplot(223)    
plt.xlabel("$\delta a_y/|\overrightarrow{a}_{exact}|$")
plt.ylabel("Density")
plt.xlim(min_error, 2*max_error)
plt.ylim(0,2)
#plt.legend(loc="center left")
plt.subplot(224)    
plt.xlabel("$\delta a_z/|\overrightarrow{a}_{exact}|$")
plt.ylabel("Density")
plt.xlim(min_error, 2*max_error)
plt.ylim(0,2)
#plt.legend(loc="center left")



plt.savefig("gravity_checks_step%d.png"%step)
