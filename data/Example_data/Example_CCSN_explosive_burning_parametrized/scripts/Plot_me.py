# Author: M. Reichert
import numpy as np
import matplotlib.pyplot as plt

fig, ax = plt.subplots(2,1,figsize=(5,5),sharex=True)
plt.subplots_adjust(hspace=0)
time,temp,rho = np.loadtxt("mainout.dat",unpack=True,usecols=[1,2,3])

# WinNet density and temperature
ax[0].plot(time-time[0],rho/1e6,color="tab:blue",lw=4,alpha=0.7,label="WinNet")
ax[1].plot(time-time[0],temp,color="tab:red",lw=4,alpha=0.7)


# Calculate the thermodynamic quantities and check if they are okay:
R0    = 0.2
rho_0 = 1e6
T9_analytic  = lambda x: 2.4*(R0)**(-3./4.)*np.exp(-x/ (3*(446/np.sqrt(7*rho_0))))
rho_analytic = lambda x: 7*rho_0 *np.exp(-x / (446/np.sqrt(7*rho_0)))
T9_gridpoint  = T9_analytic(time)
rho_gridpoint = rho_analytic(time)
ax[0].plot(time-time[0],rho_gridpoint/1e6,ls="--",color="k",label="Analytic")
ax[1].plot(time-time[0],T9_gridpoint,ls="--",color="k")


fig2 = plt.figure(figsize=(5,3))
ax2  = fig2.gca()
A,X  = np.loadtxt("finabsum.dat",unpack=True,usecols=[0,2])
ax2.plot(A,X)
ax2.set_xlim(0,80)
ax2.set_ylim(1e-8,1)
ax2.set_yscale("log")
ax2.set_title("Final mass fractions")
ax2.set_ylabel("Mass fraction X")
ax2.set_xlabel("Mass number A")


ax[0].set_ylabel(r"$\rho$ [10$^6$ g cm$^{-3}$]")
ax[1].set_ylabel("T[GK]")
ax[1].set_xlabel("Time [s]")
ax[0].set_xlim(0,3)
ax[0].legend()
fig2.savefig("final_massfractions.pdf",bbox_inches="tight")
plt.show()
