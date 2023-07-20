# Author: M. Reichert
import numpy as np
import matplotlib.pyplot as plt
# Import the WinNet python plot-routine
# Note that you can also just append the path to your .bashrc
import sys
sys.path.append('../../bin')
from class_files.winnet_class import winnet


# Read the mainout
nsm_example = winnet(".")
nsm_example.read_mainout()
# Plot the mainout
fig = nsm_example.plot_mainout()

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
ax.set_xlim(0,140)
ax.set_ylim(1e-8,1)

# Read the template file
nsm_example.read_template()
# To get an overview we can plot the considered nuclei in the calculation
fig3 = nsm_example.plot_sunet()


# Save the figures
fig.savefig("mainout_nsm_wind_martin.pdf",bbox_inches="tight")
fig2.savefig("Final_mass_fractions_nsm_wind_martin.pdf",bbox_inches="tight")
fig3.savefig("Included_nuclei_nsm_wind_martin.pdf",bbox_inches="tight")

plt.show()
