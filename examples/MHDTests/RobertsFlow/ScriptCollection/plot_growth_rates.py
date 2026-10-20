import pandas as pd
import os
from glob import glob
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({
    "font.family": "serif",
    "mathtext.fontset": "cm",    # Computer Modern
})

results_SWIFT = "SWIFT_results.csv"
results_Pencil = "Pencil_results.csv"

def load_results_csv(result_file = results_SWIFT):
    run_data = pd.read_csv(result_file, sep="\t")
    try:
        mask = run_data["Status"] == "done"
    except:
        mask = run_data["status"] == "done"
    run_data = run_data[mask]
    return run_data

data_SWIFT = load_results_csv(results_SWIFT)
data_Pencil = load_results_csv(results_Pencil)

mask_064_128 = np.array([data_Pencil['run_name'].iloc[i].split('_')[2]=='064' for i in range(len(data_Pencil))])
mask_128_256 = np.array([data_Pencil['run_name'].iloc[i].split('_')[2]=='128' for i in range(len(data_Pencil))])

#print(data_Pencil)

print(np.round(0.15 * 10**np.linspace(-1.5,0,15) ,5))


flow_kinds = np.unique([data_SWIFT['Flow_kind'].values])

symbols = {1:r'$\rm \,\,\,\,I$',2:r'$\rm \,\,\,II$',3:r'$\rm III$',4:r'$\rm IV$' }

Pencil_keys = {1:'Roberts-I',2:'Roberts-II',3:'Roberts-III',4:'Roberts-IVc' }

# Generate figure
nx = 2
ny = 2
fig, ax = plt.subplots(ny, nx, sharey=True,sharex=True, figsize=(5 * nx, 5 * ny))

fig.subplots_adjust(
    wspace=0,  # horizontal space
    hspace=0   # vertical space
)

fig.subplots_adjust(
    left=0,
    right=1,
    bottom=0,
    top=1,
    wspace=0,
    hspace=0
)

from matplotlib.ticker import AutoMinorLocator

for a in ax.flat:
    # Ticks on all four sides
    a.tick_params(axis='both',
                  which='major',
                  direction='in',
                  top=True,
                  right=True,
                  length=7,
                  width=1.2,
                  labelsize=16)

    a.tick_params(axis='both',
                  which='minor',
                  direction='in',
                  top=True,
                  right=True,
                  length=4,
                  width=1.0)

    # Add minor ticks
    a.xaxis.set_minor_locator(AutoMinorLocator())
    a.yaxis.set_minor_locator(AutoMinorLocator())

for i in range(nx):
    for j in range(ny):
        idx = i*ny + j
        if idx<len(flow_kinds):
            # read SWIFT
            flow_mask = data_SWIFT['Flow_kind']==flow_kinds[idx]
            data_flow_SWIFT= data_SWIFT[flow_mask].sort_values(by='eta')
            growth_rate = data_flow_SWIFT['growth_rate']
            growth_rate_err = data_flow_SWIFT['growth_rate_err']
            resistivity = data_flow_SWIFT['eta']
            #ax[i][j].errorbar(x=resistivity, y=growth_rate, yerr=growth_rate_err, label='SWIFT',capsize=2.5,marker='o', markersize=2.5,color='red', linestyle='')
            ax[i][j].plot(resistivity, growth_rate, marker='x', markersize=3.6,markeredgewidth=1,color='red', linestyle='-',label=r'SWIFT  $64^2\times128$')
            #ax[i][j].set_title('Flow '+symbols[flow_kinds[idx]],fontsize=18)

            # read Pencil
            file_name = 'Pencil_RF'+str(flow_kinds[idx])+'_o1.csv'
            #if ( (flow_kinds[idx]==4) | (flow_kinds[idx]==1) | (flow_kinds[idx]==2)):
            flow_mask = data_Pencil['run.in:kinematic_flow']==Pencil_keys[flow_kinds[idx]]
            print(Pencil_keys[flow_kinds[idx]])
   
            flow_mask_064_128 =( mask_064_128 & flow_mask )
            data_flow_Pencil = data_Pencil[flow_mask_064_128].sort_values(by='run.in:eta')
            resistivity = data_flow_Pencil['run.in:eta'].to_numpy()
            growth_rate = data_flow_Pencil['growth_rate'].to_numpy()
            ax[i][j].plot(resistivity,growth_rate, label=r'Pencil   $64^2\times128$',marker='+', markersize=5,markeredgewidth=1,color='black',linestyle='-',linewidth=1)
             #else:
            #    data_flow_Pencil = pd.read_csv(file_name,sep=',',header=None)
            #    data_flow_Pencil = data_flow_Pencil.sort_values(by=0).to_numpy()
            #    resistivity = data_flow_Pencil[:,0]
            #    growth_rate = data_flow_Pencil[:,1]

            flow_mask_128_256 = ( mask_128_256 & flow_mask )
            data_flow_Pencil = data_Pencil[flow_mask_128_256].sort_values(by='run.in:eta')
            resistivity = data_flow_Pencil['run.in:eta'].to_numpy()
            growth_rate = data_flow_Pencil['growth_rate'].to_numpy()
            ax[i][j].plot(resistivity, growth_rate, label=r'Pencil $128^2\times256$',marker='o', markersize=5,color='green',linestyle='-',fillstyle='none',markeredgewidth=1,linewidth=1)
 
 
            #ax[i][j].set_title('Flow '+symbols[flow_kinds[idx]],fontsize=20)
            ax[i][j].tick_params(labelsize=16)
            #ax[i][j].text(x=0.0015,y=0.13, s='Flow '+symbols[flow_kinds[idx]],fontsize=20)
            #ax[i][j].text(x=0.12,y=0.135, s='Flow '+symbols[flow_kinds[idx]],fontsize=16)
            ax[i][j].text(x=0.12,y=0.13, s='RF '+symbols[flow_kinds[idx]],fontsize=22)

ax[1][0].legend(fontsize=16, loc='lower left')
#ax[0][0].set_ylim(-0.01,0.16)
#ax[0][0].set_xlim(0.0,0.6)
ax[1][0].set_xlabel('$\eta\, k_0/ v_0 $', fontsize=20)
ax[1][1].set_xlabel('$\eta\, k_0/ v_0 $', fontsize=20)
ax[1][0].set_ylabel('$\lambda / \, (k_0 \,v_0)$', fontsize=20)
ax[0][0].set_ylabel('$\lambda / \, (k_0 \,v_0)$', fontsize=20)

ax[0][0].set_xscale('log')
#ax[0][0].set_yscale('log')
#ax[0][0].set_ylim(0.002,0.2)
#ax[0][0].set_ylim(5e-5,1e1)
ax[0][0].set_ylim(-0.05,0.15)
#ax[0][0].set_xlim(0.003,0.7)
ax[0][0].set_xlim(7e-4,0.7)
#ax[0][0].set_xlim(0.0001,1.0)
ax[0][0].axhline(y=0.0, linestyle='dotted', color='black')
ax[0][1].axhline(y=0.0, linestyle='dotted', color='black')
ax[1][0].axhline(y=0.0, linestyle='dotted', color='black')
ax[1][1].axhline(y=0.0, linestyle='dotted', color='black')

# add OW limit lines
Rm_max_pencil = 326
Rm_max_swift = 36.2
#ax[0][0].axvline(x=0.866/Rm_max_pencil, linestyle=(0, (1, 5)), color='black', ymax=0.85)
#ax[0][0].axvline(x=0.866/Rm_max_swift, linestyle=(0, (1, 10)), color='red')
#ax[0][0].text(x=0.866/Rm_max_pencil*1.2,y=-0.03, s=r'$R_M^{\rm max}(\lambda)$',fontsize=20,color='black')
#ax[0][0].text(x=0.866/Rm_max_swift*1.2,y=-0.03, s=r'$R_M^{\rm max}(3\lambda)$',fontsize=20,color='red')
#ax[0][1].axvline(x=0.866/Rm_max_pencil, linestyle=(0, (1, 5)), color='black', ymax=0.85)
#ax[0][1].axvline(x=0.866/Rm_max_swift, linestyle=(0, (1, 10)), color='red')
#ax[1][0].axvline(x=0.866/Rm_max_pencil, linestyle=(0, (1, 5)), color='black', ymax=0.85)
#ax[1][0].axvline(x=0.866/Rm_max_swift, linestyle=(0, (1, 10)), color='red')
#ax[1][1].axvline(x=1.0/Rm_max_pencil, linestyle=(0, (1, 5)), color='black', ymax=0.85)
#ax[1][1].axvline(x=1.0/Rm_max_swift, linestyle=(0, (1, 10)), color='red')

#ax[1][1].axvspan(2.8e-3,1.9e-2,color='gray', alpha=0.3,ymax = 0.85)

# gray overwinding band Pencil
ax[0][0].axvspan(7e-4,0.866/Rm_max_pencil,color='gray', alpha=0.25,ymax = 1.0)
ax[1][0].axvspan(7e-4,0.866/Rm_max_pencil,color='gray', alpha=0.25,ymax = 1.0)
ax[0][1].axvspan(7e-4,0.866/Rm_max_pencil,color='gray', alpha=0.25,ymax = 1.0)
ax[1][1].axvspan(7e-4,1.000/Rm_max_pencil,color='gray', alpha=0.25,ymax = 1.0)

# red overwidning band SWIFT
ax[0][0].axvspan(7e-4,0.866/Rm_max_swift,color='red', alpha=0.05,ymax = 1.0)
ax[1][0].axvspan(7e-4,0.866/Rm_max_swift,color='red', alpha=0.05,ymax = 1.0)
ax[0][1].axvspan(7e-4,0.866/Rm_max_swift,color='red', alpha=0.05,ymax = 1.0)
ax[1][1].axvspan(7e-4,1.000/Rm_max_swift,color='red', alpha=0.05,ymax = 1.0)

# test showing oscillatory regime
ax[1][1].text(x=6e-3,y=2e-2, s=r'osc.',fontsize=16,color='black')


# the limit for 128^2x256 Pencil run
#ax[0][1].axvline(x=0.866/Rm_max_pencil*1.5/4, linestyle=(0, (1, 5)), color='green', ymax=0.85)


# eta^0.06 line
xpts_00 = np.logspace(np.log10(1.1e-3),np.log10(1e-2),10)
ax[0][0].plot(xpts_00,0.95e-1*(xpts_00/1.1e-3)**(0.07), color='blue', linestyle='dotted')
ax[0][0].text(x=4e-3,y=9e-2, s='$\eta^{0.07}$',fontsize=20,color='blue')


fig.tight_layout()
plt.savefig('roberts_flows_growth_rates.png', dpi=220)

# save only flow II
nx = 1
ny = 1
fig, ax = plt.subplots(ny, nx, sharey=True,sharex=True, figsize=(5.5 * nx, 5 * ny))

from matplotlib.ticker import AutoMinorLocator

ax.tick_params(axis='both',
              which='major',
              direction='in',
              top=True,
              right=True,
              length=7,
              width=1.2,
              labelsize=16)

ax.tick_params(axis='both',
              which='minor',
              direction='in',
              top=True,
              right=True,
              length=4,
              width=1.0)

# Add minor ticks
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())

# read SWIFT
idx=1
flow_mask = data_SWIFT['Flow_kind']==flow_kinds[idx]
data_flow_SWIFT= data_SWIFT[flow_mask].sort_values(by='eta')
growth_rate = data_flow_SWIFT['growth_rate']
growth_rate_err = data_flow_SWIFT['growth_rate_err']
resistivity = data_flow_SWIFT['eta']
#ax[i][j].errorbar(x=resistivity, y=growth_rate, yerr=growth_rate_err, label='SWIFT',capsize=2.5,marker='o', markersize=2.5,color='red', linestyle='')
ax.plot(resistivity, growth_rate, marker='x', markersize=3.6,markeredgewidth=1,color='red', linestyle='-',label=r'SWIFT  $64^2\times128$')
#ax[i][j].set_title('Flow '+symbols[flow_kinds[idx]],fontsize=18)

# read Pencil
file_name = 'Pencil_RF'+str(flow_kinds[idx])+'_o1.csv'
#if ( (flow_kinds[idx]==4) | (flow_kinds[idx]==1) | (flow_kinds[idx]==2)):
flow_mask = data_Pencil['run.in:kinematic_flow']==Pencil_keys[flow_kinds[idx]]

flow_mask_064_128 =( mask_064_128 & flow_mask )
data_flow_Pencil = data_Pencil[flow_mask_064_128].sort_values(by='run.in:eta')
resistivity = data_flow_Pencil['run.in:eta'].to_numpy()
growth_rate = data_flow_Pencil['growth_rate'].to_numpy()
ax.plot(resistivity,growth_rate, label=r'Pencil   $64^2\times128$',marker='+', markersize=5,markeredgewidth=1,color='black',linestyle='-',linewidth=1)
 #else:
#    data_flow_Pencil = pd.read_csv(file_name,sep=',',header=None)
#    data_flow_Pencil = data_flow_Pencil.sort_values(by=0).to_numpy()
#    resistivity = data_flow_Pencil[:,0]
#    growth_rate = data_flow_Pencil[:,1]

flow_mask_128_256 = ( mask_128_256 & flow_mask )
data_flow_Pencil = data_Pencil[flow_mask_128_256].sort_values(by='run.in:eta')
resistivity = data_flow_Pencil['run.in:eta'].to_numpy()
growth_rate = data_flow_Pencil['growth_rate'].to_numpy()
ax.plot(resistivity, growth_rate, label=r'Pencil $128^2\times256$',marker='o', markersize=5,color='green',linestyle='-',fillstyle='none',markeredgewidth=1,linewidth=1)


ax.tick_params(labelsize=16)
ax.text(x=0.08,y=0.055, s='RF '+symbols[flow_kinds[idx]],fontsize=22)

ax.legend(fontsize=16, loc='lower left')
ax.set_xlabel('$\eta\, k_0/ v_0 $', fontsize=20)
ax.set_ylabel('$\lambda / \, (k_0 \,v_0) $', fontsize=20)

ax.set_xscale('log')
ax.set_ylim(-0.05,0.07)
ax.set_xlim(7e-4,0.7)
ax.axhline(y=0.0, linestyle='dotted', color='black')

# add OW limit lines
Rm_max_pencil = 326
Rm_max_swift = 36.2

# gray overwinding band Pencil
ax.axvspan(7e-4,0.866/Rm_max_pencil,color='gray', alpha=0.25,ymax = 1.0)
# red overwidning band SWIFT
ax.axvspan(7e-4,0.866/Rm_max_swift,color='red', alpha=0.05,ymax = 1.0)

# eta^0.06 line
#xpts_00 = np.logspace(np.log10(1.1e-3),np.log10(1e-2),10)
#ax.plot(xpts_00,0.95e-1*(xpts_00/1.1e-3)**(0.07), color='blue', linestyle='dotted')
#ax.text(x=4e-3,y=9e-2, s='$\eta^{0.07}$',fontsize=20,color='blue')


fig.tight_layout()
plt.savefig('roberts_flow_II_growth_rates.png', dpi=220)


