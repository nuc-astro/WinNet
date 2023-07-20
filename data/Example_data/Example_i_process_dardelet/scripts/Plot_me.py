# Author: M. Reichert
import numpy as np
import matplotlib.pyplot as plt
from matplotlib          import cm
# Import the WinNet python plot-routine
# Note that you can also just append the path to your .bashrc
import sys
sys.path.append('../../bin')
from class_files.winnet_class import winnet

# Create figures to plot at
fig_neutron_dens = plt.figure()
ax_neutron_dens  = fig_neutron_dens.gca()

fig_chart,ax_chart = plt.subplots(1,3,figsize=(15,4),sharex=True,sharey=True)
plt.subplots_adjust(hspace=0,wspace=0)


# Setup for nuclear chart plot
lw    = 0.1    # use a small linewidth for the nuclei-rectangle
min_X = 1e-10  # Min of mass fraction
max_X = 1e-4   # Max of mass fraction
cmap  = cm.jet # Use the jet colormap

# X(H) = 20%
w = winnet('./i_h20.dat/')
w.read_mainout()
t,n_dens = w.get_mainout_time(), w.get_mainout_yn()
n_dens = n_dens * 1.e4 * 6.022140857e23
ax_neutron_dens.plot(t,n_dens,label='X(H)=20%',c="b",ls='-.')
# Plot the path in the nuclear chart
w.plot_nuclear_chart_at(1e4,figure=ax_chart[0],axes_label=False,element_labels=False,fig_is_ax=True,
                        colorbar=False,min_X=min_X,max_X=max_X,cmap= cmap,nuclei_linewidths=lw)

# X(H) = 10%
w = winnet('./i_h10.dat/')
w.read_mainout()
t,n_dens = w.get_mainout_time(), w.get_mainout_yn()
n_dens = n_dens * 1.e4 * 6.022140857e23
ax_neutron_dens.plot(t,n_dens,label='X(H)=10%',c="g",ls='--')
# Plot the path in the nuclear chart
w.plot_nuclear_chart_at(1e4,figure=ax_chart[1],axes_label=False,element_labels=False,fig_is_ax=True,
                        colorbar=False,min_X=min_X,max_X=max_X,cmap= cmap,nuclei_linewidths=lw)



# X(H) = 5%
w = winnet('./i_h05.dat/')
w.read_mainout()
t,n_dens = w.get_mainout_time(), w.get_mainout_yn()
n_dens = n_dens * 1.e4 * 6.022140857e23
ax_neutron_dens.plot(t,n_dens,label='X(H)=5%',c="r",ls='-')
# Plot the path in the nuclear chart
w.plot_nuclear_chart_at(1e4,figure=ax_chart[2],axes_label=False,element_labels=False,fig_is_ax=True,colorbar=True,
                        colorbar_position=[0.27, 0.85, 0.5, 0.025],colorbar_inset=True,min_X=min_X,max_X=max_X
                        ,cmap= cmap,nuclei_linewidths=lw)


# General setup for the neutron density plot
ax_neutron_dens.set_title(r"Without $^{13}$N($p, \gamma $)$^{14}$O rate!")
ax_neutron_dens.set_yscale('log')
ax_neutron_dens.set_xscale('log')
ax_neutron_dens.set_xlim(1.e3,1.e6)
ax_neutron_dens.set_ylim(1.e11,1.e16)
ax_neutron_dens.set_xlabel('Time [s]')
ax_neutron_dens.set_ylabel(r'Neutron density [cm$^{-3}$]')
ax_neutron_dens.legend()
ax_neutron_dens.grid()

# General setup for the nuclear chart plot
ax_chart[0].set_xlim(16,129)
ax_chart[0].set_ylim(12,75)
ax_chart[0].text(0.22,0.92,"X(H)=20%",transform=ax_chart[0].transAxes,ha='right')
ax_chart[1].text(0.22,0.92,"X(H)=10%",transform=ax_chart[1].transAxes,ha='right')
ax_chart[2].text(0.22,0.92,"X(H)=5%",transform=ax_chart[2].transAxes ,ha='right')
fig_chart.text(0.5,0.1,"Neutron number",ha='center')
ax_chart[0].set_ylabel("Proton number")

# Save the neutron density plot
# Compare this plot with the one from Dardelet et al 2015, Fig. 1
fig_neutron_dens.savefig('Neutron_density_iprocess_dardelet.pdf',bbox_inches="tight")

# Save the path in the nuclear chart
fig_chart.savefig('Path_iprocess_dardelet.pdf',bbox_inches="tight")

plt.show()
