# Author: M. Reichert
import matplotlib.pyplot as plt
import numpy as np

# The here produced plot can be compared with the
# yields of the CPR2 group in
# Bliss et al. 2018 (https://ui.adsabs.harvard.edu/abs/2018ApJ...855..135B/abstract)

# Read the summed (over Z) abundances
Z,Y = np.loadtxt("finabelem.dat",unpack=True)

# Plot the final abundances versus atomic number
plt.plot(Z,Y)

# Set the plot range correctly
plt.ylim(1e-8,1)
plt.xlim(1,58)

# Other plot options
plt.yscale("log")
plt.ylabel("Abundance Y")
plt.xlabel("Atomic number Z")

# Finally save the plot
plt.savefig("CCSN_wind_bliss.pdf",bbox_inches="tight")
plt.show()
