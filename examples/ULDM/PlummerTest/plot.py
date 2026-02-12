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

parser.add_argument("--nr",
                  action="store",
                  dest="nr",
                  type=int,
                  default=32,
                  help="number of bins in r",
                  metavar=" INT")   

parser.add_argument("--nlos",
                  action="store",
                  dest="nlos",
                  type=int,
                  default=1,
                  help="number of ligne of signt",
                  metavar=" INT")   

parser.add_argument("--eps",
                     action="store",
                     dest="eps",
                     type=float,
                     default=0.1,
                     help="smoothing length",
                     metavar=" FLOAT")


parser.add_argument("--colormap",
                    action="store", 
                    default='jet',
                    help='matplotlib colormap name (e.g. mycmap, tab20c, Greys, jet, binary)') 
                    
                    

def PlummerDensity(r,a,M):
  """
  return Plummer radial density
  """
  cte = 3*M/(4*np.pi*a**3)
  return cte * (r**2/a**2  + 1)**(-5/2)

def PlummerGradientDensity(r,a,M):
  """
  return Plummer norm of the density gradient
  """
  cte = 3*M/(4*np.pi*a**3)
  return cte*(5*r/a**2)*(r**2/a**2  + 1)**(-7/2)

def PlummerLapacianDensity(r,a,M):
  """
  return Plummer radial Laplacian density
  """
  cte = 3*M/(4*np.pi*a**3)
  return cte*(5/a**2)*(4*r**2/a**2 - 3)*(r**2/a**2  + 1)**(-9/2)


  

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
  
  for filename in opt.files:
    
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
    # some info
    ################
    print("---------------------------------------------------------")
    nb.localsystem_of_units.info()
    nb.HubbleFactorCorrectionInfo()
    nb.ComovingToProperConversionInfo()
    print("---------------------------------------------------------")


    ################
    # get data
    ################    

    x = nb.rxyz()

    ###############
    # density
    ###############
    y = nb.rho
    # grid division
    rc = 1
    def f(r): return np.log(r / rc + 1.)
    def fm(r): return rc * (np.exp(r) - 1.)    
    G = libgrid.Spherical_1d_Grid(rmin=0, rmax=opt.rmax, nr=opt.nr, g=f, gm=fm)
    x = G.get_r()      
    y = G.get_MeanValMap(nb, y)
    axs[0].plot(x,y)

    ###############
    # gradient
    ###############
    y = nb.norm_grad_rho
    # grid division
    rc = 1
    def f(r): return np.log(r / rc + 1.)
    def fm(r): return rc * (np.exp(r) - 1.)    
    G = libgrid.Spherical_1d_Grid(rmin=0, rmax=opt.rmax, nr=opt.nr, g=f, gm=fm)
    x = G.get_r()      
    y = G.get_MeanValMap(nb, y)
    axs[1].plot(x,y)
    
    ###############
    # density laplacian
    ###############
    y = nb.laplacian_rho
    # grid division
    rc = 1
    def f(r): return np.log(r / rc + 1.)
    def fm(r): return rc * (np.exp(r) - 1.)    
    G = libgrid.Spherical_1d_Grid(rmin=0, rmax=opt.rmax, nr=opt.nr, g=f, gm=fm)
    x = G.get_r()      
    y = G.get_MeanValMap(nb, y)
    axs[2].plot(x,y)
    


  # Plummer ground truth
  a = 0.1
  M = 1
  
  r = G.get_r()

  #################
  # density
  #################
  rho = PlummerDensity(r,a,M)
  axs[0].plot(r,rho)
  axs[0].set_xlabel('Radius')
  axs[0].set_ylabel(r'$\rho$')
  axs[0].set_xlim(0,opt.rmax)  

  #################
  # gradient
  #################
  norm_grad_rho = PlummerGradientDensity(r,a,M)
  axs[1].plot(r,norm_grad_rho)
  axs[1].set_xlabel('Radius')
  axs[1].set_ylabel(r'$|\vec{\nabla} \rho|$')
  axs[1].set_xlim(0,opt.rmax)
  
  #################
  # laplacian
  #################
  laplacian_rho = PlummerLapacianDensity(r,a,M)
  axs[2].plot(r,laplacian_rho)
  axs[2].set_xlabel('Radius')
  axs[2].set_ylabel(r'$\nabla^2 \rho$')
  axs[2].set_xlim(0,opt.rmax)

  

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
  



