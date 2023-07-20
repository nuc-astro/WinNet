# Author: M. Reichert
# Date: 07.06.22
import numpy               as np
import matplotlib.pyplot   as plt
import matplotlib.gridspec as gridspec
from matplotlib            import cm
# Import the WinNet python plot-routine
# Note that you can also just append the path to your .bashrc
import sys
sys.path.append('../../bin')
from class_files.winnet_class import winnet



# Create Figure object
fig = plt.figure(figsize=(5,6))

# Make a grid for the plot
gs = fig.add_gridspec(nrows=3, ncols=6,
                      hspace=0, wspace=0.0)

# The mass fractions should have a large panel
ax0 = fig.add_subplot(gs[0, :])
# Panels splitted in the middle for the rest of the plots
ax1 = fig.add_subplot(gs[1, 0:3])
ax2 = fig.add_subplot(gs[1, 3:])
ax3 = fig.add_subplot(gs[2, 0:3])
ax4 = fig.add_subplot(gs[2, 3:])

# Reshape axes into different shape
ax = np.array([[ax1,ax2],[ax3,ax4]])


def plot_tr(path,color,ls2='--',ls3='-',label=""):
    """
      Function to plot mass fractions, temperature, density, radius, and electron
      fraction of one WinNet run.
    """
    # Get information of the run
    w = winnet(path)
    # Read the template and get the path to the trajectory
    w.read_template()
    path_trajectory = w.get_trajectory_path()
    # Read the mainout and the finab of the run
    w.read_mainout()
    w.read_finab()
    # Get mass numbers and mass fractions
    A,X=w.get_final_A_X()
    # Be sure that they are arrays
    A = np.array(A)
    X = np.array(X)

    # Plot the final mass fractions
    ax0.plot(A[A<220],X[A<220],color=color,ls=ls3,lw=2,label=label)
    # Make diamonds for the actinides
    ax0.scatter(A[A>220],X[A>220],color=color,marker="d",s=10)
    # Configure the axis
    ax0.set_yscale("log")
    ax0.set_ylim(2e-6,2)
    ax0.set_xlim(0,240)
    ax0.set_xlabel("Mass number")
    # Change the style of the plot, move the ticks to the top, etc.
    ax0.xaxis.tick_top()
    ax0.xaxis.set_label_position("top")
    ax0.spines["bottom"].set_linewidth(2.5)

    # Choose a linewidth of 2
    lw = 2
    # Get mainout quantities
    mainout_iteration, mainout_time, mainout_temp, mainout_dens, mainout_ye, mainout_rad, \
    mainout_yn, mainout_yp, mainout_ya, mainout_ylight, mainout_yheavy, mainout_zbar, \
    mainout_abar, mainout_entropy = w.get_mainout()
    # Read the trajectory file
    time,temp,dens,ye = np.loadtxt(path_trajectory,unpack=True,usecols=[0,1,2,4])
    # Get final and initial time
    tf = time[-1]
    t0 = mainout_time[0]
    # Create masks in order to be able to plot dashed lines for the expansion
    mask = mainout_time<tf
    mask2 = ~mask

    # Plot the temperature
    ax[0,0].plot(mainout_time[mask]-t0,mainout_temp[mask],lw=lw,color=color)
    ax[0,0].plot(mainout_time[mask2]-t0,mainout_temp[mask2],lw=lw,color=color,ls=ls2)
    ax[0,0].set_xlim(-0.1,2)

    # Plot the density
    ax[0,1].plot(mainout_time[mask]-t0,mainout_dens[mask],lw=lw,color=color)
    ax[0,1].plot(mainout_time[mask2]-t0,mainout_dens[mask2],lw=lw,color=color,ls=ls2)
    ax[0,1].set_yscale('log')
    ax[0,1].set_xlim(-0.1,2)
    ax[0,1].set_ylim(2e1,5e9)
    ax[0,1].yaxis.tick_right()
    ax[0,1].yaxis.set_label_position("right")
    ax[0,1].set_ylabel(r"Density [g cm$^{-3}$]")

    # Plot the electron fraction
    ax[1,0].plot(mainout_time[mask]-t0,mainout_ye[mask],lw=lw,color=color)
    ax[1,0].plot(mainout_time[mask2]-t0,mainout_ye[mask2],lw=lw,color=color,ls=ls2)
    ax[0,0].set_ylim(0.1,9)
    ax[1,0].set_xlim(-0.1,2)

    # Plot the radius
    ax[1,1].plot(mainout_time[mask]-t0,mainout_rad[mask],lw=lw,color=color)
    ax[1,1].plot(mainout_time[mask2]-t0,mainout_rad[mask2],lw=lw,color=color,ls=ls2)
    ax[1,1].set_yscale('log')
    ax[1,1].set_xlim(-0.1,2)
    ax[1,1].set_ylim(1e2,4e6)
    ax[1,1].yaxis.tick_right()
    ax[1,1].yaxis.set_label_position("right")
    ax[1,1].set_ylabel(r"Radius [km]")



# First plot the trajectory that get ejected promptly
path = "tr_Prompt"
plot_tr(path,"#005AA9",label="Prompt")
# Second, plot the trajectory that get ejected via a PNS shape change
path = "tr_PNS_shape"
plot_tr(path,"#009D81",label="PNS-shape")
# Third, plot the trajectory that has a high entropy
path = "tr_High_entropy"
plot_tr(path,"#EC6500",label="High entropy")

# Make the legend for the uppermost axis
ax0.legend()

# Write time on the lower end of the figure
fig.text(0.5,0.05,"Time [s]",ha="center")

# Make the xscale logarithmic and set proper time limits
for a in ax.flatten():
    a.set_xscale("log")
    a.set_xlim(5e-4,2e1)

# Write Mass fraction on the top axis
# Alternatively, one can just type ax0.set_ylabel("Mass fraction"),
# however, this will not put the y-labels on the same height
xpos = 0.02
pos = ax0.get_position()
totpos = (pos.y1+pos.y0)/2.
fig.text(xpos,totpos,"Mass fraction",rotation=90,ha="center",va="center")

# Rest of the y-labels
l = ["Temperature [GK]",r"Electron fraction"]
for i,a in enumerate(ax[:,0]):
    pos = a.get_position()
    totpos = (pos.y1+pos.y0)/2.
    fig.text(xpos,totpos,l[i],rotation=90,ha="center",va="center")

fig.savefig("comparison_rprocess.pdf",bbox_inches="tight")





def plot_nuclear_chart_at_freezeout(path,ax,label):
    """
      Plot the nuclear chart at neutron freezeout, i.e., Yn/Yh = 1.
    """
    w = winnet(path)
    w.read_mainout()
    time   = w.get_mainout_time()
    yheavy = w.get_mainout_yheavy()
    yn     = w.get_mainout_yn()
    # get index and time for neutron freezout
    index = np.argmin(abs(yn/yheavy-1))
    t_freezout = time[index]

    lw    = 0.1    # Linewidth of nuclei
    min_X = 1e-10  # Min of mass fraction
    max_X = 1e-1   # Max of mass fraction
    cmap  = cm.jet # Use the jet colormap
    # Plot the path in the nuclear chart
    w.plot_nuclear_chart_at(t_freezout,figure=ax,axes_label=False,element_labels=False,fig_is_ax=True,colorbar=True,
                            colorbar_position=[0.27, 0.85, 0.5, 0.025],colorbar_inset=True,min_X=min_X,max_X=max_X,
                            cmap= cmap,nuclei_linewidths=lw)

    ax.text(0.01,0.98,label,transform=ax.transAxes,va="top",ha="left")
    ax.set_xlim(60,170)
    ax.set_ylim(40,95)

fig_chart,ax_chart = plt.subplots(1,3,figsize=(15,4),sharex=True,sharey=True)

plt.subplots_adjust(wspace=0)

# First plot the trajectory that get ejected promptly
path = "tr_Prompt"
plot_nuclear_chart_at_freezeout(path,ax_chart[0],"Prompt")
# Second, plot the trajectory that get ejected via a PNS shape change
path = "tr_PNS_shape"
plot_nuclear_chart_at_freezeout(path,ax_chart[1],"PNS shape")
# Third, plot the trajectory that has a high entropy
path = "tr_High_entropy"
plot_nuclear_chart_at_freezeout(path,ax_chart[2],"High entropy")

fig_chart.savefig("nuclear_chart_at_freezeout.pdf",bbox_inches="tight")
plt.show()
