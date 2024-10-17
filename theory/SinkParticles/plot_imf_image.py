#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  12 13:28 2024

This script creates an image of some IMF that we split into the continuous part and discrete part.
This image illustrates the star spawning algorithm with GEAR sink particles

@author: Darwin Roduit
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib

# For swift doc, use the following and not the styleheet
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Computer Modern Roman']

mmin = 0.5
mmax = 300

matplotlib.rcParams.update({'font.size': 16})

figsize = (6.4, 4.8)
fig, ax = plt.subplots(num=1, ncols=1, nrows=1,
                       figsize=figsize, layout="tight")
ax.set_xlim([0.3, 400])
ax.set_ylim([1, 5e4])

ax.set_xticks([1, 2, 4, 8, 20, 50, 100, 300])
ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

ax.set_xlabel("$M_{\star}$ $[M_\odot]$")
ax.set_ylabel("d$N/$d$M$ [arbitrary units]")
ax.set_xscale('log')
ax.set_yscale('log')

# theoretical imf
s = -1.3
bins = 10**np.linspace(np.log10(mmin), np.log10(mmax), 100)
n = 0.9*10000*bins**s
ax.plot(bins, n, "k--")

bins = 10**np.linspace(np.log10(mmin), np.log10(8), 100)
n = 0.9*10000*bins**s
ax.fill_between(bins, 0.1, n, color="red", alpha=0.1)

bins = 10**np.linspace(np.log10(8), np.log10(mmax), 100)
n = 0.9*10000*bins**s
ax.fill_between(bins, 0.1, n, color="blue", alpha=0.1)

ax.text(2, 1e2, r"$M_{\rm c}$", horizontalalignment='center')
ax.text(50, 2, r"$M_{\rm d}$", horizontalalignment='center')

# Add limit
ax.vlines(x=8, ymin=0, ymax=600, color='k', linestyle='-')

# Add text to the vertical line
ax.text(8, 800, r"$m_{t}$", horizontalalignment='center')

# fig.patch.set_facecolor('none')  # Remove figure background
# ax.set_facecolor('none')         # Remove axes background

plt.savefig('sink_imf.png', dpi=300, bbox_inches='tight')
plt.close()
