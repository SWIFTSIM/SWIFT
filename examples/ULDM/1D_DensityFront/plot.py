#!/usr/bin/python3
###########################################################################################
#  package:   pNbody
#  file:      plotSphericalProfile
#  brief:     
#  copyright: GPLv3
#             Copyright (C) 2019 EPFL (Ecole Polytechnique Federale de Lausanne)
#             LASTRO - Laboratory of Astrophysics of EPFL
#  author:    Yves Revaz <yves.revaz@epfl.ch>
#
# This file is part of pNbody.
###########################################################################################

import os
import numpy as np
import argparse

import matplotlib.pyplot as plt
from pNbody import plot
from pNbody import *
from pNbody import cosmo
from pNbody import libutil

from astropy import units as u
from astropy import constants as c

import copy

####################################################################
# option parser
####################################################################

description="""plot the laplacian density for each particle as a function of the radius."""
epilog     ="""
Examples:
--------


"""

parser = argparse.ArgumentParser(description=description,epilog=epilog,formatter_class=argparse.RawDescriptionHelpFormatter)

plot.add_files_options(parser)
plot.add_arguments_units(parser)
plot.add_arguments_reduc(parser)
plot.add_arguments_center(parser)
plot.add_arguments_select(parser)
plot.add_arguments_info(parser)
plot.add_arguments_legend(parser)
plot.add_arguments_icshift(parser)
plot.add_arguments_cmd(parser)

parser.add_argument(action="store", 
                    dest="files", 
                    metavar='FILE', 
                    type=str,
                    default=None,
                    nargs='*',
                    help='a list of files')                     


parser.add_argument("-o",
                    action="store",
                    type=str,
                    dest="outputfilename",
                    default=None,
                    help="Name of the output file")  



parser.add_argument('--xmin',
                    action="store", 
                    dest="xmin", 
                    metavar='FLOAT', 
                    type=float,
                    default=None,
                    help='x min')

parser.add_argument('--xmax',
                    action="store", 
                    dest="xmax", 
                    metavar='FLOAT', 
                    type=float,
                    default=None,
                    help='x max')
                    
parser.add_argument('--ymin',
                    action="store", 
                    dest="ymin", 
                    metavar='FLOAT', 
                    type=float,
                    default=None,
                    help='y min')

parser.add_argument('--ymax',
                    action="store", 
                    dest="ymax", 
                    metavar='FLOAT', 
                    type=float,
                    default=None,
                    help='y max')

                                        
parser.add_argument('--log',
                    action="store", 
                    dest="log", 
                    metavar='STR', 
                    type=str,
                    default=None,
                    help='log scale (None,x,y,xy)')


parser.add_argument('-y','--y',
                    action="store", 
                    dest="y", 
                    metavar='STR', 
                    type=str,
                    default='density',
                    help='quantity to plot in the y axis')

parser.add_argument("--rmax",
                  action="store",
                  dest="rmax",
                  type=float,
                  default=0.5,
                  help="max radius of bins",
                  metavar=" FLOAT")


parser.add_argument("--colormap",
                    action="store", 
                    default='jet',
                    help='matplotlib colormap name (e.g. mycmap, tab20c, Greys, jet, binary)') 
                    
                    



#######################################
# MakePlot
#######################################


def MakePlot(opt):
  
  params = {
    "axes.labelsize": 22,
    "axes.titlesize": 22,
    "font.size": 22,
    "legend.fontsize": 22,
    "xtick.labelsize": 22,
    "ytick.labelsize": 22,
    "text.usetex": True,
    "figure.subplot.left": 0.15,
    "figure.subplot.right": 0.95,
    "figure.subplot.bottom": 0.15,
    "figure.subplot.top": 0.95,
    "figure.subplot.wspace": 0.2,
    "figure.subplot.hspace": 0.2,
    "figure.figsize" : (8, 6),
    "lines.markersize": 6,
    "lines.linewidth": 2.0,
  }
  plt.rcParams.update(params)
  


  # get a list of color
  colors = plot.ColorList(n=len(opt.files),colormap=opt.colormap)
  
  # list of points
  datas = []


  fig, axs = plt.subplots(3, 1, figsize=(12, 16))
  
  
  #################################
  # loop over files
  
  for i,filename in enumerate(opt.files):
    
    nb = Nbody(filename, ftype=opt.ftype)
    
    ################
    # units
    ################
    
    # define local units
    unit_params = plot.apply_arguments_units(opt)
    nb.set_local_system_of_units(params=unit_params)
    
    ################
    # apply options
    ################
    nb = plot.apply_arguments_icshift(nb, opt)
    nb = plot.apply_arguments_reduc(nb, opt)
    nb = plot.apply_arguments_select(nb, opt)
    nb = plot.apply_arguments_center(nb, opt)
    nb = plot.apply_arguments_cmd(nb, opt)
    nb = plot.apply_arguments_info(nb, opt)
    nb = plot.apply_arguments_display(nb, opt)
    nb = plot.apply_arguments_verbose(nb, opt)


    ################
    # plot the density
    ################    

    x = nb.x()
    y = nb.rho

    axs[0].plot(x,y,label=filename,c='r',ls=':')


    ################
    # plot the gradient
    ################    

    x = nb.x()
    y = nb.grad_rho[:,0]/nb.rho

    axs[1].plot(x,y,label=filename,c='r',ls=':')
    
    
  # ground truth
    
  xmax  = 10000
  x     = np.linspace(0,xmax,1000) # kpc
  x0    = 5000                     # kpc
  rho0  = 1e-8/3                   # 1e10Msol/kpc^3
  c     = 1.0
  sigma = 500                      # kpc
   
  def fct_Rho(x):
    return rho0 * (c+1-np.tanh((x-x0)/sigma))
   
  def fct_GradRho(x):
    return -rho0/sigma * 1/np.cosh((x-x0)/sigma)**2


  
  # density
  rho  = fct_Rho(x)
  axs[0].plot(x,rho,c='k')
  axs[0].set_xlabel(r'$r\,\rm{[kpc]}$')
  axs[0].set_ylabel(r'$\rho [10^{10}\,\rm{M_{\odot}\cdot kpc^{-3}}$')
  axs[0].set_xlim(1000,9000)
  axs[0].set_ylim(0.3e-8,1.1e-8)

  # density gradient
  rho  = fct_GradRho(x)/fct_Rho(x)
  axs[1].plot(x,rho,c='k')
  axs[1].set_xlabel(r'$r\,\rm{[kpc]}$')
  axs[1].set_ylabel(r'$\nabla\rho /\rho$')
  axs[1].set_xlim(1000,9000)
  axs[1].set_ylim(-1.1e-3,0.1e-3)
  

  # save or display
  if opt.outputfilename:
    plt.savefig(opt.outputfilename)
  else:
    plt.show()    
      


#################################
# main
#################################

if __name__ == '__main__':  
  
  opt = parser.parse_args()
  MakePlot(opt)
  



