import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt


keys=[
"Step",
"Time",
"Scale-factor",
"Redshift",
"Total_mass",
"Gas_mass",
"DM_mass",
"Sink_mass",
"Star_mass",
"BH_mass",
"Gas_Z_mass",
"Star_Z_mass",
"BH_Z_mass",
"Kin_Energy",
"Int_Energy",
"Pot_energy",
"Rad_energy",
"Gas_Entropy",
"CoM_x",
"CoM_y",
"CoM_z",
"Mom_x",
"Mom_y",
"Mom_z",
"Ang_mom_x",
"Ang_mom_y",
"Ang_mom_z",
"BH_acc_rate",
"BH_acc_mass",
"BH_sub_mass",
"Gas_H_mass",
"Gas_H2_mass",
"Gas_HI_mass",
"Gas_He_mass",
"Mag_Energy",
"DivB_err",
"Cr_Helicity",
"Mag_Helicity",
"BH_bol_lum",
"BH_jet_power"] 



description="Plot the content of a statistics.txt Swift output file."
epilog     ="""
Examples:
--------
sw_plotStatistics statistics.txt
sw_plotStatistics.py statistics.txt --x Time --y Rel_Tot_Energy

"""

parser = argparse.ArgumentParser(description=description,epilog=epilog,formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument(action="store", 
                    dest="files", 
                    metavar='FILE', 
                    type=str,
                    nargs="*",
                    default=None,
                    help='lis of files') 

parser.add_argument("--mode",
                    action="store", 
                    dest="mode", 
                    metavar='STRING', 
                    required=False,
                    type=str,
                    default=None,
                    help='plotting mode') 

parser.add_argument("--x",
                    action="store", 
                    dest="x", 
                    metavar='STRING', 
                    required=False,
                    type=str,
                    default="Step",
                    help='x value %s'%keys)

parser.add_argument("--y",
                    action="store", 
                    dest="y", 
                    metavar='STRING', 
                    required=False,
                    type=str,
                    default="Rel_Tot_Energy",
                    help='y value %s'%keys)



parser.add_argument("-o",
                    action="store", 
                    dest="outputfile", 
                    metavar='STRING', 
                    required=False,
                    type=str,
                    default=None,
                    help='save the output to the given filename')




params = {
    "axes.labelsize": 14,
    "axes.titlesize": 18,
    "font.size": 12,
    "legend.fontsize": 12,
    "xtick.labelsize": 14,
    "ytick.labelsize": 14,
    "text.usetex": True,
    "figure.subplot.left": 0.19,
    "figure.subplot.right": 0.98,
    "figure.subplot.bottom": 0.12,
    "figure.subplot.top": 0.95,
    "figure.subplot.wspace": 0.14,
    "figure.subplot.hspace": 0.12,
    "figure.figsize" : (8,6),
    "lines.markersize": 6,
    "lines.linewidth": 3.0,
}
plt.rcParams.update(params)




def getValueAndLabel(value):
  
  if value == "Wall-clock_time":
    label = "Wall clock time [hour]"
    value = data[value] / 1000. / 3600.
    value = np.add.accumulate(value)
    return value,label  
  
  elif value == "Step":
    label = "Step"
    value = data[value]
    return value,label  

  elif value == "Rel_Tot_Energy":
    label = "Relative Total Energy [\%]"
    value = data[value]
    return value,label  
  
  else:
    label = value
    value = data[value]
    return value,label      
    

  




##########################################
#
#  M A I N
#
##########################################

opt = parser.parse_args()



############################
# load the data

datas = []

for fle in opt.files:

  # do not use loadtxt as it crashes when a line is cutted
  #data_tmp = np.loadtxt(f,comments="#")

  fle = open(fle,'r')
  lines = fle.readlines()
  fle.close()
  lines = lines[:-1]
  boo = list(map(lambda x:x[0]!="#", lines))
  # get the last commented lines to read the labels
  #tmp,count=np.unique(boo,return_counts=True)
  #idx = count[0]-1
  #keys = str.split(lines[idx-1])[1:]


  # keep non commented lines
  lines = np.compress(boo,lines)
  lines = list(map(str.split, lines))
  
  lines = np.array(list(map(np.array, lines)))
  data_tmp = lines.astype(np.double)
  

  # create a dictionary
  data = {}
  for i,key in enumerate(keys):
    data[key] = data_tmp[:,i]
  
  del data_tmp

  datas.append(data)


# add additional quantities
datas[0]['Tot_Energy']     = datas[0]["Kin_Energy"] + datas[0]["Int_Energy"] +  datas[0]["Pot_energy"] + datas[0]["Rad_energy"]
datas[0]['Rel_Tot_Energy'] = 100*(datas[0]['Tot_Energy'] - datas[0]['Tot_Energy'][0])/datas[0]['Tot_Energy']


  
############################
# do the plot

for data in datas:
  
  if opt.mode is None:
  
    x,xlabel = getValueAndLabel(opt.x)
    y,ylabel = getValueAndLabel(opt.y)
    
    ax = plt.gca()
    ax.plot(x,y)

    


ax = plt.gca()
ax.set_xlabel(xlabel)
ax.set_ylabel(ylabel)


if opt.outputfile is not None:
  plt.savefig(opt.outputfile)
else:
  plt.show()







