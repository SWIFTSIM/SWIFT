import pandas as pd
import os
from glob import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

plt.rcParams.update({
    "font.family": "serif",
    "mathtext.fontset": "cm",    # Computer Modern
})

def load_results_csv(result_file):
    run_data = pd.read_csv(result_file, sep="\t")
    try:
        mask = run_data["Status"] == "done"
    except:
        mask = run_data["status"] == "done"
    run_data = run_data[mask]
    return run_data

def load_results_txt(result_file):
    run_data = pd.read_csv(result_file, sep=r'\s+')
    return run_data

results_Pencil_2nd = "./data_for_export_Pencil_one_mode/5em2_2nd_rerun.txt" #"5em2_2nd.txt" #"./data_for_export_Pencil/5em2_2nd_rerun.txt" #"5em2_2nd.txt"
results_Pencil_6th = "./data_for_export_Pencil_one_mode/5em2_6th_rerun.txt" #"5em2_6th.txt" #"./data_for_export_Pencil/5em2_6th_rerun.txt"  #"5em2_6th.txt"
results_Pencil_10th = "./data_for_export_Pencil_one_mode/5em2_10th_rerun.txt" #"5em2_10th.txt" #"./data_for_export_Pencil/5em2_10th_rerun.txt" #"5em2_10th.txt"
results_SWIFT = "5em2_SWIFT.csv"

#data_SWIFT = load_results_csv(results_SWIFT)
data_Pencil_2nd = load_results_txt(results_Pencil_2nd)
data_Pencil_6th = load_results_txt(results_Pencil_6th)
data_Pencil_10th = load_results_txt(results_Pencil_10th)
print(data_Pencil_10th)
data_SWIFT = load_results_csv(results_SWIFT)

keys_Pencil = ['lam','Bm/Brms','-WL/B2', 'Jrms/Brms', 'JB/B2', 'JBrms/B2'] 
keys_SWIFT = ['growth_rate_Brms_stat','B_bar_rms','WL','Jrms','JdotBrms','JdotB_mean']

# rename these SWIFT columns
rename_dict = dict(zip(keys_SWIFT, keys_Pencil))
data_SWIFT = data_SWIFT.rename(columns=rename_dict)

exact_Pencil_10th = data_Pencil_10th.iloc[0]
data_Pencil_10th = data_Pencil_10th.iloc[1:]
#exact_Pencil_10th = data_Pencil_10th.iloc[-1]
#data_Pencil_10th = data_Pencil_10th.iloc[:-1]
exact_SWIFT = data_SWIFT.iloc[-1]
data_SWIFT = data_SWIFT.iloc[:-1]
# subtract exact result
for key in keys_Pencil:
    data_Pencil_10th[key]=np.abs(data_Pencil_10th[key] - exact_Pencil_10th[key])
    data_Pencil_2nd[key]=np.abs(data_Pencil_2nd[key] - exact_Pencil_10th[key])
    data_Pencil_6th[key]=np.abs(data_Pencil_6th[key] - exact_Pencil_10th[key])

for key in keys_Pencil:
    data_SWIFT[key]=np.abs(data_SWIFT[key] - exact_SWIFT[key])

#
#for i in range(nx):
#    for j in range(ny):
#        idx = i*ny + j
#        if idx<len(flow_kinds):
#            # read SWIFT
#            flow_mask = data_SWIFT['Flow_kind']==flow_kinds[idx]
#            data_flow_SWIFT= data_SWIFT[flow_mask].sort_values(by='eta')
#            growth_rate = data_flow_SWIFT['growth_rate']
#            growth_rate_err = data_flow_SWIFT['growth_rate_err']
#            resistivity = data_flow_SWIFT['eta']
#            #ax[i][j].errorbar(x=resistivity, y=growth_rate, yerr=growth_rate_err, label='SWIFT',capsize=2.5,marker='o', markersize=2.5,color='red', linestyle='')
#            ax[i][j].plot(resistivity, growth_rate, marker='x', markersize=7.1,markeredgewidth=2,color='red', linestyle='',label='SWIFT')
#            #ax[i][j].set_title('Flow '+symbols[flow_kinds[idx]],fontsize=18)
#
#            # read Pencil
#            file_name = 'Pencil_RF'+str(flow_kinds[idx])+'_o1.csv'
#            #if ( (flow_kinds[idx]==4) | (flow_kinds[idx]==1) | (flow_kinds[idx]==2)):
#            flow_mask = data_Pencil['run.in:kinematic_flow']==Pencil_keys[flow_kinds[idx]]
#            print(Pencil_keys[flow_kinds[idx]])
#            data_flow_Pencil = data_Pencil[flow_mask].sort_values(by='run.in:eta')
#            resistivity = data_flow_Pencil['run.in:eta'].to_numpy()
#            growth_rate = data_flow_Pencil['growth_rate'].to_numpy()
#            #else:
#            #    data_flow_Pencil = pd.read_csv(file_name,sep=',',header=None)
#            #    data_flow_Pencil = data_flow_Pencil.sort_values(by=0).to_numpy()
#            #    resistivity = data_flow_Pencil[:,0]
#            #    growth_rate = data_flow_Pencil[:,1]
#            ax[i][j].plot(resistivity,growth_rate, label='Pencil',marker='+', markersize=10,markeredgewidth=2,color='black',linestyle='')
#         
#            #ax[i][j].set_title('Flow '+symbols[flow_kinds[idx]],fontsize=20)
#            ax[i][j].tick_params(labelsize=16)
#            ax[i][j].text(x=0.14,y=0.13, s='Flow '+symbols[flow_kinds[idx]],fontsize=20)
#
#ax[1][0].legend(fontsize=16, loc='lower left')
##ax[0][0].set_ylim(-0.01,0.16)
##ax[0][0].set_xlim(0.0,0.6)
#ax[1][0].set_xlabel('$\eta\, k_0/ u_0 $', fontsize=20)
#ax[1][1].set_xlabel('$\eta\, k_0/ u_0 $', fontsize=20)
#ax[1][0].set_ylabel('$\lambda$', fontsize=20)
#ax[0][0].set_ylabel('$\lambda$', fontsize=20)
#
#ax[0][0].set_xscale('log')
#ax[0][0].set_yscale('log')
#ax[0][0].set_ylim(0.002,0.2)
#ax[0][0].set_xlim(0.003,0.7)
#
#
#fig.tight_layout()
#plt.savefig('roberts_flows_growth_rates.png', dpi=220)

def plot_results():

    #keys_to_plot = ['lam','Bm/Brms','-WL/B2', 'Jrms/Brms', 'JB/B2', 'JBrms/B2']
    #key_names = [r'$ \rm err \left( \lambda \right) $',r'$\rm err \left( \overline{B}_{rms}/B_{rms} \right)$',r'$\rm err \left( W_L/v_{rms}B_{rms}^2 \right)$',r'$ \rm err \left(  J_{rms}/B_{rms} \right)$',r'$\rm err \left( (\vec J \cdot \vec B)_{rms}/B_{rms}^2 \right)$',r'$\rm err \left( \langle\vec J \cdot \vec B \rangle \right)/B_{rms}^2$'] 
 
    keys_to_plot = ['lam','Bm/Brms','-WL/B2', 'Jrms/Brms', 'JBrms/B2', 'JB/B2']
    key_names = [r'$ \rm err \left( \lambda \right) $',r'$\rm err \left( \overline{B}_{rms}/B_{rms} \right)$',r'$\rm err \left( W_L/v_{rms}B_{rms}^2 \right)$',r'$ \rm err \left(  J_{rms}/B_{rms} \right)$',r'$\rm err \left( \langle\vec J \cdot \vec B \rangle \right)/B_{rms}^2$',r'$\rm err \left( (\vec J \cdot \vec B)_{rms}/B_{rms}^2 \right)$'] 
 

    Nkeys = len(keys_to_plot)
    ncol = 2
    nrow = int(Nkeys/ncol)
    fig, ax = plt.subplots(nrow, ncol, sharey=True,sharex=True, figsize=(5 * ncol, 4.5 * nrow))


    fig.subplots_adjust(
        wspace=0,  # horizontal space
        hspace=0   # vertical space
    )
    
    
    #delta = 0.1
    fig.subplots_adjust(
        left=0.1,
        right=0.95,
        bottom=0.05,
        top=1.0,
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
    
    for i in range(nrow):
        for j in range(ncol):
             index = i*ncol+j
             if index<Nkeys:
                 column_name = keys_to_plot[index]
                 col_Pencil_2nd = data_Pencil_2nd[column_name] 
                 Nx_2nd = data_Pencil_2nd['Nx'] 
                 col_Pencil_6th = data_Pencil_6th[column_name]
                 Nx_6th = data_Pencil_6th['Nx'] 
                 col_Pencil_10th = data_Pencil_10th[column_name]
                 Nx_10th = data_Pencil_10th['Nx'] 
                 col_SWIFT = data_SWIFT[column_name]

                 N = np.array([8,16,32,64,128]) #data_Pencil_10th['Nx']

                 ax[i][j].plot(N[1:],col_SWIFT,marker='x',markersize=7,linestyle='dashed',color='red', label='SWIFT',fillstyle='none',markeredgewidth=2)
                 ax[i][j].plot(Nx_2nd,col_Pencil_2nd,marker='s',markersize=6,linestyle='dashed',color='black', label='Pencil 2nd',fillstyle='none',markeredgewidth=2)
                 ax[i][j].plot(Nx_6th,col_Pencil_6th,marker='o',markersize=6,linestyle='dashdot',color='black', label='Pencil 6th',fillstyle='none',markeredgewidth=2)
                 ax[i][j].plot(Nx_10th,col_Pencil_10th,marker='+',markersize=10,linestyle='solid',color='black', label='Pencil 10th',markeredgewidth=2)

                 ax[i][j].set_yscale('log')
                 ax[i][j].set_xscale('log')
                 #ax[i, j].set_aspect("equal", adjustable="box")
                 ax[i][j].set_ylabel(key_names[index],fontsize=20)
                 if j==1:
                     ax[i][j].yaxis.set_label_position('right')
                 #ax[i][j].set_ylim([7e-9,2e0])
                 ax[i][j].set_ylim([7e-12,2e0])
                 ax[i][j].set_xticks(N,N,fontsize=16)
                 #ax[i][j].tick_params(labelsize=16)
                 ax[i][j].set_yticks([1e0,1e-1,1e-2,1e-3,1e-4,1e-5,1e-6,1e-7,1e-8,1e-9,1e-10,1e-11],['$10^0$','$10^{-1}$','$10^{-2}$','$10^{-3}$','$10^{-4}$','$10^{-5}$','$10^{-6}$','$10^{-7}$','$10^{-8}$','$10^{-9}$','$10^{-10}$','$10^{-11}$'],fontsize=20)
                 ax[i][j].xaxis.set_minor_locator(mticker.NullLocator())
                 ax[i][j].tick_params(labelsize=16)
        ax[-1][0].legend(fontsize=14, loc='lower left')
        ax[-1][0].set_xlabel('$ N_\perp $', fontsize=20)
        ax[-1][1].set_xlabel('$ N_\perp $', fontsize=20)
  
        # N^-2 lines
        xpts_00 = np.logspace(np.log10(64),np.log10(128),10)
        ax[0][0].plot(xpts_00,7e-4*(xpts_00/64)**(-2), color='blue', linestyle='dotted')

        xpts_01 = np.logspace(np.log10(64),np.log10(128),10)
        ax[0][1].plot(xpts_01,5e-4*(xpts_01/64)**(-2), color='blue', linestyle='dotted')

        xpts_10 = np.logspace(np.log10(64),np.log10(128),10)
        ax[1][0].plot(xpts_01,2e-2*(xpts_01/64)**(-2), color='blue', linestyle='dotted')

        xpts_11 = np.logspace(np.log10(64),np.log10(128),10)
        ax[1][1].plot(xpts_01,2e-2*(xpts_01/64)**(-2), color='blue', linestyle='dotted')

        xpts_20 = np.logspace(np.log10(64),np.log10(128),10)
        ax[2][0].plot(xpts_01,4e-3*(xpts_01/64)**(-2), color='blue', linestyle='dotted')

        xpts_21 = np.logspace(np.log10(64),np.log10(128),10)
        ax[2][1].plot(xpts_01,3e-3*(xpts_01/64)**(-2), color='blue', linestyle='dotted')

        # N^-6 lines
        xpts_00 = np.logspace(np.log10(64),np.log10(128),10)
        ax[0][0].plot(xpts_00,3e-6*(xpts_00/64)**(-6), color='blue', linestyle=(0, (1, 5)))

        xpts_01 = np.logspace(np.log10(64),np.log10(128),10)
        ax[0][1].plot(xpts_01,1e-5*(xpts_01/64)**(-6), color='blue', linestyle=(0, (1, 5)))

        xpts_10 = np.logspace(np.log10(64),np.log10(128),10)
        ax[1][0].plot(xpts_10,1e-4*(xpts_10/64)**(-6), color='blue', linestyle=(0, (1, 5)))

        xpts_11 = np.logspace(np.log10(64),np.log10(128),10)
        ax[1][1].plot(xpts_11,4e-4*(xpts_11/64)**(-6), color='blue', linestyle=(0, (1, 5)))

        xpts_20 = np.logspace(np.log10(64),np.log10(128),10)
        ax[2][0].plot(xpts_20,1e-4*(xpts_20/64)**(-6), color='blue', linestyle=(0, (1, 5)))

        xpts_21 = np.logspace(np.log10(64),np.log10(128),10)
        ax[2][1].plot(xpts_21,8e-6*(xpts_21/64)**(-6), color='blue', linestyle=(0, (1, 5)))

        # N^-10 lines
        xpts_00 = np.logspace(np.log10(64),np.log10(128),10)
        ax[0][0].plot(xpts_00,3e-8*(xpts_00/64)**(-10), color='blue', linestyle=(0, (1, 10)))

        xpts_01 = np.logspace(np.log10(64),np.log10(128),10)
        ax[0][1].plot(xpts_01,2e-7*(xpts_01/64)**(-10), color='blue', linestyle=(0, (1, 10)))

        xpts_10 = np.logspace(np.log10(64),np.log10(128),10)
        ax[1][0].plot(xpts_10,4e-6*(xpts_10/64)**(-10), color='blue', linestyle=(0, (1, 10)))

        xpts_11 = np.logspace(np.log10(64),np.log10(128),10)
        ax[1][1].plot(xpts_11,2e-5*(xpts_11/64)**(-10), color='blue', linestyle=(0, (1, 10)))

        xpts_20 = np.logspace(np.log10(64),np.log10(128),10)
        ax[2][0].plot(xpts_20,2e-6*(xpts_20/64)**(-10), color='blue', linestyle=(0, (1, 10)))

        xpts_21 = np.logspace(np.log10(64),np.log10(128),10)
        ax[2][1].plot(xpts_21,1e-7*(xpts_21/64)**(-10), color='blue', linestyle=(0, (1, 10)))

        # line names
        ax[0][0].text(x=100,y=1.5e-5, s='$N_\perp^{-2}$',fontsize=14,color='blue')
        ax[0][0].text(x=100,y=5e-7, s='$N_\perp^{-6}$',fontsize=14,color='blue')
        ax[0][0].text(x=100,y=1e-9, s='$N_\perp^{-10}$',fontsize=14,color='blue')

    plt.savefig('5em2.png',dpi=220)

def plot_results_flow_I():

    #keys_to_plot = ['lam','Bm/Brms','-WL/B2', 'Jrms/Brms', 'JB/B2', 'JBrms/B2']
    #key_names = [r'$ \rm err \left( \lambda \right) $',r'$\rm err \left( \overline{B}_{rms}/B_{rms} \right)$',r'$\rm err \left( W_L/v_{rms}B_{rms}^2 \right)$',r'$ \rm err \left(  J_{rms}/B_{rms} \right)$',r'$\rm err \left( (\vec J \cdot \vec B)_{rms}/B_{rms}^2 \right)$',r'$\rm err \left( \langle\vec J \cdot \vec B \rangle \right)/B_{rms}^2$'] 
 
    key_to_plot = 'lam'
    key_name = r'$ \rm err \left( \lambda \right) $'

    ncol = 1
    nrow = 1#int(Nkeys/ncol)
    fig, ax = plt.subplots(nrow, ncol, sharey=True,sharex=True, figsize=(5 * ncol, 4.5 * nrow))

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
    
    column_name = key_to_plot
    col_Pencil_2nd = data_Pencil_2nd[column_name] 
    Nx_2nd = data_Pencil_2nd['Nx'] 
    col_Pencil_6th = data_Pencil_6th[column_name]
    Nx_6th = data_Pencil_6th['Nx'] 
    col_Pencil_10th = data_Pencil_10th[column_name]
    Nx_10th = data_Pencil_10th['Nx'] 
    col_SWIFT = data_SWIFT[column_name]

    N = np.array([8,16,32,64,128]) #data_Pencil_10th['Nx']

    ax.plot(N[1:],col_SWIFT,marker='x',markersize=7,linestyle='dashed',color='red', label='SWIFT',fillstyle='none',markeredgewidth=2)
    ax.plot(Nx_2nd,col_Pencil_2nd,marker='s',markersize=6,linestyle='dashed',color='black', label='Pencil 2nd',fillstyle='none',markeredgewidth=2)
    ax.plot(Nx_6th,col_Pencil_6th,marker='o',markersize=6,linestyle='dashdot',color='black', label='Pencil 6th',fillstyle='none',markeredgewidth=2)
    ax.plot(Nx_10th,col_Pencil_10th,marker='+',markersize=10,linestyle='solid',color='black', label='Pencil 10th',markeredgewidth=2)

    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_ylabel(key_name,fontsize=20)
    ax.set_ylim([7e-12,2e0])
    ax.set_xticks(N,N,fontsize=16)
    ax.set_yticks([1e0,1e-1,1e-2,1e-3,1e-4,1e-5,1e-6,1e-7,1e-8,1e-9,1e-10,1e-11],['$10^0$','$10^{-1}$','$10^{-2}$','$10^{-3}$','$10^{-4}$','$10^{-5}$','$10^{-6}$','$10^{-7}$','$10^{-8}$','$10^{-9}$','$10^{-10}$','$10^{-11}$'],fontsize=20)
    ax.xaxis.set_minor_locator(mticker.NullLocator())
    ax.tick_params(labelsize=16)
    ax.legend(fontsize=14, loc='lower left')
    ax.set_xlabel('$ N_\perp $', fontsize=20)
  
    # N^-2 lines
    xpts_00 = np.logspace(np.log10(64),np.log10(128),10)
    ax.plot(xpts_00,7e-4*(xpts_00/64)**(-2), color='blue', linestyle='dotted')

    # N^-6 lines
    xpts_00 = np.logspace(np.log10(64),np.log10(128),10)
    ax.plot(xpts_00,3e-6*(xpts_00/64)**(-6), color='blue', linestyle=(0, (1, 5)))

    # N^-10 lines
    xpts_00 = np.logspace(np.log10(64),np.log10(128),10)
    ax.plot(xpts_00,3e-8*(xpts_00/64)**(-10), color='blue', linestyle=(0, (1, 10)))

    # line names
    ax.text(x=100,y=1.5e-5, s='$N_\perp^{-2}$',fontsize=14,color='blue')
    ax.text(x=100,y=5e-7, s='$N_\perp^{-6}$',fontsize=14,color='blue')
    ax.text(x=100,y=1e-9, s='$N_\perp^{-10}$',fontsize=14,color='blue')

    fig.tight_layout()

    plt.savefig('5em2_lam.png',dpi=220)


def find_growth_rate(N, error):
    a,b = np.polyfit(np.log(N), np.log(error), 1, cov=False)
    n = - round(a,1)
    N0 = round(np.exp(-b/a),1)
    return [n,N0]


def generate_fit_table():

 
    keys_to_plot = ['lam','Bm/Brms','-WL/B2', 'Jrms/Brms', 'JBrms/B2', 'JB/B2']
    key_names = [r'$ \rm err \left( \lambda \right) $',r'$\rm err \left( \overline{B}_{rms}/B_{rms} \right)$',r'$\rm err \left( W_L/v_{rms}B_{rms}^2 \right)$',r'$ \rm err \left(  J_{rms}/B_{rms} \right)$',r'$\rm err \left( \langle\vec J \cdot \vec B \rangle \right)/B_{rms}^2$',r'$\rm err \left( (\vec J \cdot \vec B)_{rms}/B_{rms}^2 \right)$'] 
 

    Nkeys = len(keys_to_plot)
    ncol = 2
    nrow = int(Nkeys/ncol)

    tab_res=[]
   
    for i in range(nrow):
        for j in range(ncol):
             index = i*ncol+j
             if index<Nkeys:
                 column_name = keys_to_plot[index]
                 col_Pencil_2nd = data_Pencil_2nd[column_name].to_numpy()
                 Nx_2nd = data_Pencil_2nd['Nx'].to_numpy()
                 idx_2nd = np.argsort(Nx_2nd)
                 col_Pencil_6th = data_Pencil_6th[column_name].to_numpy()
                 Nx_6th = data_Pencil_6th['Nx'].to_numpy()
                 idx_6th = np.argsort(Nx_6th)
                 col_Pencil_10th = data_Pencil_10th[column_name].to_numpy()
                 Nx_10th = data_Pencil_10th['Nx'].to_numpy() 
                 idx_10th = np.argsort(Nx_10th)
                 col_SWIFT = data_SWIFT[column_name].to_numpy()

                 N = np.array([8,16,32,64,128]) #data_Pencil_10th['Nx']

                 res_2nd = find_growth_rate(Nx_2nd[idx_2nd][-3:], col_Pencil_2nd[idx_2nd][-3:])
                 res_6th = find_growth_rate(Nx_6th[idx_6th][-2:], col_Pencil_6th[idx_6th][-2:])
                 res_10th = find_growth_rate(Nx_10th[idx_10th][-2:], col_Pencil_10th[idx_10th][-2:])
                 res_SWIFT = find_growth_rate(N[-3:], col_SWIFT[-3:])
 

                 # required precsision at Nperp = 256 10th order
                 err_128_10th = (128/res_10th[-1])**(-res_10th[0])
                 err_256_10th = (256/res_10th[-1])**(-res_10th[0])
                 print("reference run error order, Nperp = 128", err_128_10th)
                 print("reference run error order, Nperp = 256", err_256_10th)
 
                 tab_res.append({'scheme':'Pencil_2nd','quantity':column_name,'n_tilt':res_2nd[0], 'N0_shift':res_2nd[-1]})
                 tab_res.append({'scheme':'Pencil_6th','quantity':column_name,'n_tilt':res_6th[0], 'N0_shift':res_6th[-1]})
                 tab_res.append({'scheme':'Pencil_10th','quantity':column_name,'n_tilt':res_10th[0], 'N0_shift':res_10th[-1]})
                 tab_res.append({'scheme':'SWIFT','quantity':column_name,'n_tilt':res_SWIFT[0], 'N0_shift':res_SWIFT[-1]})

    tab_res = pd.DataFrame(tab_res)
    tab_res = tab_res.sort_values(by='scheme')
    return tab_res 

plot_results()
plot_results_flow_I()


# create fitting table
tab = generate_fit_table()
print(tab)
