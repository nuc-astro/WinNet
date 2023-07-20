# Author: M. Reichert
import numpy as np
import matplotlib.pyplot as plt
# Import the WinNet python plot-routine
# Note that you can also just append the path to your .bashrc
import sys
sys.path.append('../../bin')
from class_files.winnet_class import winnet


fig  = plt.figure()
ax   = fig.gca()

# Load data from the mainout
mrsn_example = winnet('.')

# Read and plot the tracked nuclei
mrsn_example.read_tracked_nuclei()
time = mrsn_example.get_tracked_time()
nucs = ["eu151","eu153","th232","u236"]
label = [r"$^{151}$Eu","$^{153}$Eu","$^{232}$Th","$^{236}$U"]
for ind,n in enumerate(nucs):
    X = mrsn_example.get_tracked_nuclei(n,massfraction=True)
    ax.plot(time,X,label=label[ind],lw=2)
ax.legend(loc="upper left")
ax.set_yscale("log")
ax.set_xscale("log")
ax.set_ylim(1e-7,1e-2)
ax.set_xlim(1e-1,1e17)
ax.set_ylabel("Mass fraction")
ax.set_xlabel("Time [s]")


# Also plot the final abundances
fig2 = plt.figure()
ax2   = fig2.gca()
# Read the finab
mrsn_example.read_finab()
A,X = mrsn_example.get_final_A_X()
# Plot the massfraction, summed over equal mass numbers
ax2.plot(A,X)
ax2.set_ylabel("Mass fraction")
ax2.set_xlabel("Mass number")
ax2.set_yscale("log")
ax2.set_xlim(0,249)
ax2.set_ylim(1e-8,1)

# fig2.savefig("Final_mass_fractions.pdf",bbox_inches="tight")

plt.show()
