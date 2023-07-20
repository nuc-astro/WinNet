# Author: M. Reichert
import numpy as np
import matplotlib.pyplot as plt
# Import the WinNet python plot-routine
# Note that you can also just append the path to your .bashrc
import sys
sys.path.append('../../bin')
from class_files.winnet_class import winnet


fig  = plt.figure(figsize=(8,4))
ax   = fig.gca()

mrsn_example = winnet('.')

# Load data from the finab
mrsn_example.read_finab()
# Plot the isotopes from the finab as a chain versus mass number
mrsn_example.plot_final_isotopes(figure=fig,lower_limit=1e-5)
ax.set_ylim(1e-5,1)


# Setup a second plot
fig2  = plt.figure(figsize=(8,4))
ax2   = fig2.gca()
# Load data from the timescales
mrsn_example.read_timescales()
mrsn_example.plot_timescales(figure=fig2)

ax2.legend()

ax2.set_xlim(3e-2,1e2)
ax2.set_ylim(1e-11,1e13)

fig.savefig("final_mass_fractions_isotopes_weakr_obergaulinger.pdf",bbox_inches="tight")
fig2.savefig("timescales_weakr_obergaulinger.pdf",bbox_inches="tight")

plt.show()
