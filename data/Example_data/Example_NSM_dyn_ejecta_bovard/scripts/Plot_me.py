# Author: M. Reichert
import numpy             as np
import matplotlib.pyplot as plt
from matplotlib          import cm
# Import the WinNet python plot-routine
# Note that you can also just append the path to your .bashrc
import sys
sys.path.append('../../bin')
from class_files.winnet_class import winnet

# Create a first Figure.
fig = plt.figure()
ax  = fig.gca()
# Read the finab of a normal representative case
nsm_example = winnet('tracer_492.dat')
nsm_example.read_finab()
A,X = nsm_example.get_final_A_X()
# Plot the massfraction, summed over equal mass numbers
ax.plot(A,X,lw=2,color="tab:blue",label="Representative with high mass")

# Read the finab of a high entropy case (see Bovard et al. 2017)
nsm_example_2 = winnet('tracer_1131.dat')
nsm_example_2.read_finab()
A,X = nsm_example_2.get_final_A_X()
# Plot the massfraction, summed over equal mass numbers
ax.plot(A,X,lw=2,color="tab:orange",label="High entropy with low mass")




ax.set_ylabel("Mass fraction")
ax.set_xlabel("Mass number")
ax.set_yscale("log")
ax.set_xlim(100,219)
ax.set_ylim(1e-6,1)
ax.legend(loc="upper left")
plt.savefig("massfractions_nsm_dyn_ejecta_bovard.pdf",bbox_inches="tight")
plt.show()
