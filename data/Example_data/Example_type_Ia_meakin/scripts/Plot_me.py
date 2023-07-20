# Author: M. Reichert
import numpy as np
import matplotlib.pyplot as plt

# Within this script, we are recreating Fig. 13 and Fig. 14
# of Meakin et al 2009
# ( https://iopscience.iop.org/article/10.1088/0004-637X/693/2/1188/pdf )


# First lets plot the thermodynamic conditions. The temperature should be a
# power law, the electron-fraction starts at 0.5, and the density is
# created with the Timmes EOS according to adiabatic conditions.
# Initially, we set up a figure with two y-axis:
fig = plt.figure(figsize=(5,4))
ax  = fig.gca()
ax2 = ax.twinx()
# Now we read the thermodynamic conditions
time, temp, rho, ye = np.loadtxt("mainout.dat",unpack=True,usecols=[1,2,3,4])

# Plot temperature, density, electron-fraction
ax.plot(time,temp   ,color="k"      ,label="Temperature")
ax.plot(time,rho/1e8,color="tab:red",label="Density")
ax2.plot(time,ye                    ,label=r"Y$_e$")
# Make the plot more beautiful, put captions and so on...
ax.set_yscale("log")
ax2.set_ylim(0.4984,0.5006)
ax2.set_ylabel("Electron fraction Y$_e$")
ax.set_ylim(7e-3,11)
ax.set_xlim(0,0.7)
ax.legend()
ax2.legend(loc="center right")
ax.set_ylabel(r"Temperature [GK]     Density [$10^8$ g cm$^{-3}$]")
ax.set_xlabel("Time [s]")
# Save the plot, It is very similar to Fig. 13 in Meakin et al. 2009.
fig.savefig("trajectory.pdf",bbox_inches="tight")


# Now, lets recreate Fig. 14. For this, we want to plot
# H, He4, Fe52, Fe53, Fe54, Co55, Ni56, Ni57, Ni58 versus temperature
fig2 = plt.figure(figsize=(5,4))
ax   = fig2.gca()

# Give a legend NSE/Network
l1, = ax.plot(np.nan,np.nan,color="k")
l2, = ax.plot(np.nan,np.nan,color="k",ls="--")
legend1 = plt.legend([l1,l2], ["Network", "NSE"], loc="upper left",fontsize=8)
ax.add_artist(legend1)

time,Yh,Yhe4,Yfe52,Yfe53,Yfe54,Yco55,Yni56,Yni57,Yni58 = np.loadtxt("tracked_nuclei.dat",unpack=True)
ax.plot(temp,Yh      ,color="tab:red"    ,label=r"$^1$H"    )
ax.plot(temp,Yhe4*4  ,color="tab:orange" ,label=r"$^4$He"   )
ax.plot(temp,Yfe52*52,color="tab:cyan"   ,label=r"$^{52}$Fe")
ax.plot(temp,Yfe53*53,color="tab:green"  ,label=r"$^{53}$Fe")
ax.plot(temp,Yfe54*54,color="lightgreen" ,label=r"$^{54}$Fe")
ax.plot(temp,Yco55*55,color="navy"       ,label=r"$^{55}$Co")
ax.plot(temp,Yni56*56,color="k"          ,label=r"$^{56}$Ni")
ax.plot(temp,Yni57*57,color="magenta"    ,label=r"$^{57}$Ni")
ax.plot(temp,Yni58*58,color="saddlebrown",label=r"$^{58}$Ni")


# We calculated the NSE abundances beforehand. For this,
# we just set the NSE temperature to a very low temperature.
# Now we only have to read and plot it.
time,temp,Yh,Yhe4,Yfe52,Yfe53,Yfe54,Yco55,Yni56,Yni57,Yni58 = np.load("nse_nuclei.npy")
ax.plot(temp,Yh      ,color="tab:red"    ,ls="--")
ax.plot(temp,Yhe4*4  ,color="tab:orange" ,ls="--")
ax.plot(temp,Yfe52*52,color="tab:cyan"   ,ls="--")
ax.plot(temp,Yfe53*53,color="tab:green"  ,ls="--")
ax.plot(temp,Yfe54*54,color="lightgreen" ,ls="--")
ax.plot(temp,Yco55*55,color="navy"       ,ls="--")
ax.plot(temp,Yni56*56,color="k"          ,ls="--")
ax.plot(temp,Yni57*57,color="magenta"    ,ls="--")
ax.plot(temp,Yni58*58,color="saddlebrown",ls="--")

ax.set_ylabel("Mass fraction")
ax.set_xlabel("Temperature [GK]")

ax.legend(fontsize=8)
ax.set_ylim(1e-4,1)
ax.set_xlim(7,2.25)
ax.set_yscale("log")

# Save the figure. We note that the evolution looks different compared to
# the original publication. In our case, mass fractions of the iron group
# freeze out much earlier and alphas reach much lower mass fractions.
# The NSE mass fractions are very similar.
fig2.savefig("mass_fractions.pdf",bbox_inches="tight")


# Ultimately, lets also plot the final abundances
fig3 = plt.figure(figsize=(5,4))
ax   = fig3.gca()
# Read the mass fractions summed over equal A
A, Y, X = np.loadtxt("finabsum.dat",unpack=True)

ax.plot(A,X)
ax.set_yscale("log")
ax.set_xlabel("Mass number A")
ax.set_ylabel("Mass fraction X")
fig3.savefig("final_mafras.pdf",bbox_inches="tight")
plt.show()
