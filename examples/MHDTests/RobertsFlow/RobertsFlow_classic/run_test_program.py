#!/usr/bin/env python3

# Libraries
import subprocess
import numpy as np
import pandas as pd
from glob import glob
import os

# Set up run parameters here
parameters = {'v0':[5.0,5.0,5.0,5.0],'eta':[0.2,0.5,1.0,1.5],'kv':[1,2,3,4],'kb':[1,2,3,4],'Lbox':[1,2,3,4],'Flow_kind':[0,1,2,3], 'Scheme':['ODI','FDI','VP','VP'], 'IAfile':['g64','g64','g64','g64'],'artificial_diffusion':[None,None,None,None],'hyperbolic_dedner':[None,None,None,None],'hyperbolic_dedner_divv':[None,None,None,None],'parabolic_dedner':[None,None,None,None]}
parameter_data = pd.DataFrame(data = parameters)
parameter_data.to_csv('test_run_parameters.csv', sep = ';',index=False)

parameter_data = pd.read_csv('test_run_parameters.csv', sep = ';')
parameter_data = parameter_data.replace({np.nan: None})

# Dictionaries for shortcuts
schemes_dict = {'ODI':'direct-induction','FDI':'direct-induction-fede','VP':'vector-potential'}
IA_dict = {'g16':'glassCube_16.hdf5','g32':'glassCube_32.hdf5','g64':'glassCube_64.hdf5'}

# Number of threads
threads = 18

# Where to store results
results_directory_name = 'test_results/'

# Function for configuring simulation
def configure_simulation(scheme, forcing, spline, eos, path_to_lib = True):

 path_to_configure = ' ../../../../'
 scheme_opt = ' --with-spmhd='+scheme
 kernel_opt = ' --with-kernel='+spline
 forcing_opt = ' --with-forcing='+forcing
 eos_opt = ' --with-equation-of-state='+eos 
 hand_vec = ' --disable-hand-vec'
 compiler_warnings = ' --disable-compiler-warnings'
 doxygen = ' --disable-doxygen-doc'
	
 if path_to_lib==True: 
  path_to_libraries = ' --with-fftw=/opt/homebrew/Cellar/fftw/3.3.10_1/ --with-gsl=/opt/homebrew/Cellar/gsl/2.7.1/ --with-hdf5=/opt/homebrew/Cellar/hdf5/1.14.2/bin/h5cc'
 else:
  path_to_libraries = ''

 this_directory = os.getcwd()

 goto_swift_directory_command = ' cd '+path_to_configure
 configure_command = ' ./configure'+kernel_opt+scheme_opt+forcing_opt+eos_opt+hand_vec+compiler_warnings+doxygen+path_to_libraries
 make_command = ' make -j'
 goto_this_directory_command = ' cd '+this_directory

 command_sandwich = '( '+goto_swift_directory_command + ' &&' + configure_command + ' &&' + make_command + ' &&' + goto_this_directory_command + ' )'  

 #subprocess.call(command_sandwich, shell=True)
 print(command_sandwich) 

 print('configuring swift complete')

# Funciton for making ICfile
def make_IC(phys_parameters, IAfile):

 v0 = ' '+str(phys_parameters['v0'].values[0])
 eta = ' '+str(phys_parameters['eta'].values[0])
 kv = ' '+str(phys_parameters['kv'].values[0])
 kb = ' '+str(phys_parameters['kb'].values[0])
 Lbox = ' '+str(phys_parameters['Lbox'].values[0])

 Flow_kind = str(int(phys_parameters['Flow_kind'].values[0]))

 IC_folder = './flow_'+Flow_kind
 
 this_directory = os.getcwd()
 go_to_flow_directory = ' cd '+IC_folder
 command = ' python '+'make_IC.py' + v0 + eta + ' ./../IAfiles/'+IAfile + kv + kb + Lbox
 return_back = ' cd '+this_directory
 command_sandwich = ' ('+go_to_flow_directory+' &&'+command+' &&'+return_back+' ) '

 #subprocess.call(command_sandwich, shell=True)
 print(command_sandwich) 

 print('Created ICs')

# Function for running simulation
def run_simulation(phys_parameters, threads):
 v0 = phys_parameters['v0'].values[0]
 eta = phys_parameters['eta'].values[0]
 kv = phys_parameters['kv'].values[0]
 Flow_kind = phys_parameters['Flow_kind'].values[0]

 artificial_diffusion = phys_parameters['artificial_diffusion'].values[0]
 hyperbolic_dedner = phys_parameters['hyperbolic_dedner'].values[0]
 hyperbolic_dedner_divv = phys_parameters['hyperbolic_dedner_divv'].values[0]
 parabolic_dedner = phys_parameters['parabolic_dedner'].values[0]

 with_hydro = ' --hydro'
 threads_opt = ' --threads='+str(threads)
 
 f_pref = ' -P RobertsFlowForcing:'
 set_u0 = ''
 if v0!=None:
  set_u0 = f_pref+'u0:'+str(v0)

 set_kv = ''
 if kv!=None:
  set_kv = f_pref+'kv:'+str(kv)

 set_Flow_kind = ''
 if Flow_kind!=None:
  set_Flow_kind = f_pref+'Flow_kind:'+str(int(Flow_kind))
 
 set_forcing_par = set_u0+set_kv+set_Flow_kind

 MHD_pref = ' -P MHD:'
 set_eta = ''
 if eta!=None:
  set_eta = MHD_pref+'resistive_eta:'+str(eta)

 set_artificial_diffusion = ''
 if artificial_diffusion!=None:
  set_artificial_diffusion = MHD_pref+'artificial_diffusion:'+str(artificial_diffusion)  
 
 set_hyperbolic_dedner = ''
 if hyperbolic_dedner!=None:
  set_hyperbolic_dedner = MHD_pref+'hyperbolic_dedner:'+str(hyperbolic_dedner)

 set_hyperbolic_dedner_divv = ''
 if hyperbolic_dedner_divv!=None:
  set_hyperbolic_dedner = MHD_pref+'hyperbolic_dedner_divv:'+str(hyperbolic_dedner_divv)

 set_parabolic_dedner = ''
 if parabolic_dedner!=None:
  set_parabolic_dedner = MHD_pref+'parabolic_dedner:'+str(parabolic_dedner)

 set_MHD_par = set_eta+set_artificial_diffusion+set_hyperbolic_dedner+set_hyperbolic_dedner_divv+set_parabolic_dedner

 set_all_par = set_forcing_par + set_MHD_par

 parameter_file = ' ./../RobertsFlow.yml' 

 write_output_file = ' | tee output_log.txt'

 IC_folder = './flow_'+str(int(Flow_kind))

 this_directory = os.getcwd()
 go_to_flow_directory = ' cd '+IC_folder
 command = ' ../../../../../swift'+with_hydro+threads_opt+set_all_par+parameter_file+write_output_file
 return_back = ' cd '+this_directory
 command_sandwich = ' ('+go_to_flow_directory+' &&'+command+' &&'+return_back+' ) '
 
 #subprocess.call(command_sandwich, shell=True)
 print(command_sandwich)

#def move_run_results_to_directory()
 
def run_all(the_parameters):
 for i in range(len(the_parameters)):
  parameters_for_the_run = the_parameters.iloc[[i]]
  print(parameters_for_the_run)

  scheme = schemes_dict[parameters_for_the_run['Scheme'].values[0]]
  configure_simulation(scheme,'roberts-flow','quintic-spline','isothermal-gas')

  IAfile = IA_dict[parameters_for_the_run['IAfile'].values[0]]
  make_IC(parameters_for_the_run,IAfile)

  run_simulation(parameters_for_the_run, threads=18)

  
  


run_all(parameter_data)

#configure_simulation('direct-induction', 'roberts-flow', 'quintic-spline', 'isothermal-gas')

#make_IC(parameter_data.head(1),'glassCube_16.hdf5')

#run_simulation(parameter_data.head(1),18)
