
#from swiftsimio import load_statistics
from unyt import g, s, statA
import numpy as np
import matplotlib.pyplot as plt
import sys

# load reference
with_reference = False

folder='reference'
if with_reference:
    import pandas as pd
    Berr_vs_t_data = pd.read_csv('./'+folder+'/Berr_vs_t_Tricco_Price_2012.csv',header=None)

fig, ax = plt.subplots(2, 1, figsize=((2*4, 4*4)))

file = sys.argv[1]
output = sys.argv[2]

data = np.loadtxt(file)

t = data[:,1]
   
divB_err = data[:,35] 
#Brms = data[:,38]

#ln1 = ax1.semilogy(t, Brms/Brms[0],'k-', label=r'$B_{rms}$')
#ln2 = ax2.semilogy(t, divB_err, 'r:', label=r'Error')
ax[0].plot(t,divB_err, label = 'SWIFT', color='black')

if with_reference:
    ax[0].plot(Berr_vs_t_data[0],Berr_vs_t_data[1],label='Terrence S.Tricco, Daniel J. Price. 2012',color='red')


ax[0].set_yscale('log')
ax[0].set_ylim([1e-2,1e-0])
ax[0].set_xlim([0.0,0.2])

ax[0].set_xlabel(r'$t$')
ax[0].set_ylabel(r'$ h |\nabla \cdot B|/ B $')

ax[0].legend()


from swiftsimio import load
from swiftsimio.visualisation.slice import slice_gas
# Load snapshot
filename = './CurlB_0000.hdf5'
data = load(filename)
# Retrieve some information about the simulation run
artDiffusion = data.metadata.hydro_scheme["Artificial Diffusion Constant"]
dedHyp = data.metadata.hydro_scheme["Dedner Hyperbolic Constant"]
dedHypDivv = data.metadata.hydro_scheme["Dedner Hyperbolic div(v) Constant"]
dedPar = data.metadata.hydro_scheme["Dedner Parabolic Constant"]
eta = data.metadata.hydro_scheme["Resistive Eta"]
git = data.metadata.code["Git Revision"]
gitBranch = data.metadata.code["Git Branch"]
hydroScheme = data.metadata.hydro_scheme["Scheme"]
kernel = data.metadata.hydro_scheme["Kernel function"]
neighbours = data.metadata.hydro_scheme["Kernel target N_ngb"]
n_gas = data.metadata.n_gas

# Display mean divB and h
divB = data.gas.magnetic_divergences
print(np.max(np.abs(divB)),np.min(np.abs(divB)))
print(np.abs(divB))

# Add panel with infromation about the run
text_common_args = dict(
  fontsize=10, ha="center", va="center", transform=ax[1].transAxes
)
ax[1].text(
  0.5,
  0.8,
  "Monopole test 3D",
  **text_common_args,
)
ax[1].text(0.5, 0.7, "SWIFT %s" % git.decode("utf-8"), **text_common_args)
ax[1].text(0.5, 0.6, "Branch %s" % gitBranch.decode("utf-8"), **text_common_args)
ax[1].text(0.5, 0.5, hydroScheme.decode("utf-8"), **text_common_args)
ax[1].text(
  0.5,
  0.4,
  kernel.decode("utf-8") + " with $%.2f$ neighbours" % (neighbours),
  **text_common_args,
)
ax[1].text(
  0.5, 0.3, "Artificial diffusion: $%.2f$ " % (artDiffusion), **text_common_args
)
ax[1].text(
  0.5,
  0.2,
  "Dedner Hyp, Hyp_div(v), Par: $%.2f,%.2f,%.2f$ " % (dedHyp, dedHypDivv, dedPar),
  **text_common_args,
)
ax[1].text(
  0.5, 0.1, r"Physical resistivity: $\eta$: $%.2f$ " % (eta), **text_common_args
)
ax[1].text(
  0.5, 0.0, r"Number of particles: $N_p$: $%.0f$ " % (n_gas), **text_common_args
)
ax[1].axis("off")


plt.tight_layout()

#plt.savefig('DivergencePlot.png', dpi=100)
plt.savefig(output, dpi=100)

to_plot="hist"
########################################## Make histogram
if to_plot == "hist":
    quantity = (
        np.abs(divB).to_physical()
    ).value  # (rho.to_physical()/rho_mean).value#(abs_vec(v).to_physical()/vrms).value
    qname = "$divB$"  #'rho/<rho>' #'$|v|/v_{rms}$'
    Nbins = 100
    # bins = np.linspace(np.min(quantity),np.max(quantity),Nbins)
    bins = np.logspace(
        -1, 1, Nbins
    )
    fig, ax = plt.subplots()
    plt.hist(quantity, bins=bins, color="blue", label=qname, alpha=0.5, density=False)
    #plt.hist(
    #    quantity[mask],
    #    bins=bins,
    #    color="red",
    #    label=qname + "($R_0$>" + str(min_metric) + ")",
    #    alpha=0.5,
    #    density=False,
    #)
    plt.xscale("log")
    plt.yscale("log")
    plt.legend()
    ax.set_ylabel("$N_{particles}$/bin")
    ax.set_xlabel(qname)
    
    plt.savefig('hist', dpi=100)
