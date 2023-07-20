# Author: M. Reichert
import numpy as np
import matplotlib.pyplot as plt
from matplotlib          import cm
# Import the WinNet python plot-routine
# Note that you can also just append the path to your .bashrc
import sys
sys.path.append('../../bin')
from class_files.winnet_class import winnet


# First, plot the final abundances
# Create figures to plot at
fig_finab = plt.figure(figsize=(5,3))
ax_finab  = fig_finab.gca()

# Plot the result with default reaclib decays
# The table has been created out of the reaclib data
path = "beta_decay_reaclib.dat/finabsum.dat"
A,Y,X = np.loadtxt(path,unpack=True)
ax_finab.plot(A,X,label='REACLIB')

# Plot the Moeller beta decay data
# The decay data was accessed at
# https://www.matthewmumpower.com/publications/paper/2019/moller/nuclear-properties-for-astrophysical-and-radioactive-ion-beam-applications-ii
# and converted to the WinNet format
# The publication of Moeller et al. 2019 to these rates can be found here:
# https://www.sciencedirect.com/science/article/pii/S0092640X18300263?via%3Dihub
path = "beta_decay_moeller.dat/finabsum.dat"
A,Y,X = np.loadtxt(path,unpack=True)
ax_finab.plot(A,X,label='Möller et al. 2019')

# Plot the result with Marketin beta decay data
# The data and paper of Marketin et al. 2016 can be
# accessed at:
# https://journals.aps.org/prc/abstract/10.1103/PhysRevC.93.025805
path = "beta_decay_marketin.dat/finabsum.dat"
A,Y,X = np.loadtxt(path,unpack=True)
ax_finab.plot(A,X,label='Marketin et al. 2016')

# Make the plot look nice
ax_finab.legend(fontsize=7)
ax_finab.set_yscale('log')
ax_finab.set_ylim(1e-5,1)
ax_finab.set_xlim(55,240)
ax_finab.set_ylabel("Mass fraction")
ax_finab.set_xlabel("Mass number")

# Also plot the mass fractions in the nuclear chart at freezeout
fig_chart,ax_chart = plt.subplots(1,3,figsize=(15,4),sharex=True,sharey=True)
plt.subplots_adjust(hspace=0,wspace=0)


# Setup for nuclear chart plot
lw    = 0.1        # use a small linewidth for the nuclei-rectangle
min_X = 1e-12      # Min of mass fraction
max_X = 1e-4       # Max of mass fraction
cmap  = cm.magma_r # Use the jet colormap

# Reaclib
w = winnet('./beta_decay_reaclib.dat/')
w.read_mainout()
t,yn,yh = w.get_mainout_time(), w.get_mainout_yn(), w.get_mainout_yheavy()
# Find the neutron freezeout time (yn/yh = 1)
t_freezeout = t[np.argmin(np.abs(yn/yh-1))]

# Plot the path in the nuclear chart
w.plot_nuclear_chart_at(t_freezeout,figure=ax_chart[0],axes_label=False,element_labels=False,fig_is_ax=True,
                        colorbar=False,min_X=min_X,max_X=max_X,cmap= cmap,nuclei_linewidths=lw)
# Label the plot
ax_chart[0].text(0.01,0.98,"Reaclib",transform=ax_chart[0].transAxes,va='top',ha='left',fontsize=13,color="grey")


# Moeller
w = winnet('./beta_decay_moeller.dat/')
w.read_mainout()
t,yn,yh = w.get_mainout_time(), w.get_mainout_yn(), w.get_mainout_yheavy()
# Find the neutron freezeout time (yn/yh = 1)
t_freezeout = t[np.argmin(np.abs(yn/yh-1))]
# Plot the path in the nuclear chart
w.plot_nuclear_chart_at(t_freezeout,figure=ax_chart[1],axes_label=False,element_labels=False,fig_is_ax=True,
                        colorbar=False,min_X=min_X,max_X=max_X,cmap= cmap,nuclei_linewidths=lw)
# Label the plot
ax_chart[1].text(0.01,0.98,"Möller et al. 2019",transform=ax_chart[1].transAxes,va='top',ha='left',fontsize=13,color="grey")


# Marketin
w = winnet('./beta_decay_marketin.dat/')
w.read_mainout()
t,yn,yh = w.get_mainout_time(), w.get_mainout_yn(), w.get_mainout_yheavy()
# Find the neutron freezeout time (yn/yh = 1)
t_freezeout = t[np.argmin(np.abs(yn/yh-1))]
# Plot the path in the nuclear chart
w.plot_nuclear_chart_at(t_freezeout,figure=ax_chart[2],axes_label=False,element_labels=False,fig_is_ax=True,colorbar=True,
                        colorbar_position=[0.27, 0.85, 0.5, 0.025],colorbar_inset=True,min_X=min_X,max_X=max_X
                        ,cmap= cmap,nuclei_linewidths=lw)
# Label the plot
ax_chart[2].text(0.01,0.98,"Marketin et al. 2016",transform=ax_chart[2].transAxes,va='top',ha='left',fontsize=13,color="grey")


# Set the labels
ax_chart[0].set_ylabel("Proton number")
ax_chart[1].set_xlabel("Neutron number")
# Set the correct limits
[a.set_xlim(50,150) for a in ax_chart]
[a.set_ylim(30,85) for a in ax_chart]

# Save and show the figures
fig_finab.savefig("final_massfractions.pdf",bbox_inches='tight')
fig_chart.savefig("nuclear_chart_beta_decay.pdf",bbox_inches='tight')
plt.show()