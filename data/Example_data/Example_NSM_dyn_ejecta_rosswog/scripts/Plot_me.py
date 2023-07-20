# Author: M. Reichert
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
# Import the WinNet python plot-routine
# Note that you can also just append the path to your .bashrc
import sys
sys.path.append('../../bin')
from class_files.winnet_class import winnet

# Create a first Figure.
fig,ax = plt.subplots(1,2,figsize=(10,4),sharex=True)
plt.subplots_adjust(wspace=0)

# Load data from the mainout
nsm_example = winnet('.')
nsm_example.read_mainout()
time   = nsm_example.get_mainout_time()
yheavy = nsm_example.get_mainout_yheavy()
yn     = nsm_example.get_mainout_yn()

# Plot the neutron to seed ratio
ax[0].plot(time,yn/yheavy,lw=2)
# Calculate freezeout point
index_freezout = np.argmin(abs(yn/yheavy-1))
time_freezout  = time[index_freezout]
# ... and plot it
ax[0].plot([time_freezout,time_freezout],[0,1],ls='--',color='grey',lw=2)
ax[0].text(time_freezout,1e-6,"t = "+str(round(time_freezout,3))+' s',ha='right',va='center',rotation=90,color='k',fontsize=9)
ax[0].plot([0,time_freezout],[1,1],ls='--',color='grey',lw=2)
ax[0].text(5e-2,1,r"Y$_n$/Y$_{\mathrm{Heavy}}$ = 1",ha='center',va='bottom',color='k',fontsize=9)
ax[0].scatter(time_freezout,1,marker='x',color='r',s=100)
ax[0].text(time_freezout+1,1.1,"Neutron freezout",ha='left',va='center',color='r',fontsize=9)

# Read the energy generation
energy_file = "generated_energy.dat"
df = pd.read_csv(energy_file,skiprows=2,header=None,sep='\s+')
with open(energy_file,"r") as f:
    lines = f.readlines()
    header = lines[1].split()[1:]
header = list(map(lambda x: x.strip(),header))
df.columns = header

time  = df["time[s]"].values
engen = df["Engen(Total)"].values
# ... and plot it
ax[1].plot(time,engen,label="Total",lw=2.5)
ax[1].plot(time,df["Engen(fiss)"].values,color="tab:orange" ,label="Fission",lw=1)
ax[1].plot(time,df["Engen(weak)"].values,color="tab:cyan"   ,label="Weak",lw=1)
ax[1].plot(time,df["Engen(n,g)"].values,color="tab:green"   ,label=r"($n,\gamma$)",lw=1)
ax[1].plot(time,df["Engen(a,g)"].values,color="saddlebrown" ,label=r"($\alpha,\gamma$)",lw=1)
index_freezout = np.argmin(np.abs(time-time_freezout))
time_freezout  = time[index_freezout]
ax[1].plot([time_freezout,time_freezout],[0,engen[index_freezout]],ls='--',color='grey',lw=2)
ax[1].scatter(time_freezout,engen[index_freezout],marker='x',color='r',s=100,zorder=100)
ax[1].text(time_freezout+1,1.1*engen[index_freezout],"Neutron freezout",ha='left',va='center',color='r',fontsize=9)
ax[1].text(5e-2,2.3*engen[index_freezout],r"$\propto$ constant",ha='center',va='center',color='k',fontsize=10)
ax[1].text(5e2,1e14,r"$\propto$ t$^{- \alpha}$",ha='center',va='center',rotation=-42,color='k',fontsize=10)

ax[1].legend()

# Set the range and labels of the plots
ax[0].set_xlim(3e-3,5e5)
ax[0].set_ylim(1e-9,1e5)
ax[1].set_ylim(1e9,1e20)
ax[0].loglog()
ax[1].loglog()
ax[0].set_ylabel(r"Y$_n$/Y$_{\mathrm{Heavy}}$")
ax[1].set_ylabel(r"Energy generation [$\mathrm{erg} \, \mathrm{g}^{-1} \, \mathrm{s}^{-1}$]")
fig.text(0.5,0.008,r"Time [s]")

# Set the ticks of the right plot to the right
ax[1].yaxis.set_label_position("right")
ax[1].yaxis.tick_right()



# Also plot the final abundances
fig2 = plt.figure()
ax   = fig2.gca()
# Read the finab
nsm_example.read_finab()
A,X = nsm_example.get_final_A_X()
# Plot the massfraction, summed over equal mass numbers
ax.plot(A,X)
ax.set_ylabel("Mass fraction")
ax.set_xlabel("Mass number")
ax.set_yscale("log")
ax.set_xlim(0,249)
ax.set_ylim(1e-8,1)

fig.savefig("Neutron_to_seed_and_heating_power_nsm_dyn_ejecta_rosswog.pdf",bbox_inches="tight")
fig2.savefig("Final_mass_fractions_nsm_dyn_ejecta_rosswog.pdf",bbox_inches="tight")

plt.show()
