# Author: M. Reichert
import h5py
import sys
import numpy as np
import matplotlib.pyplot as plt

# This script will only work in case the h5fc compiler was enabled in the makefile.
# Otherwise, no hdf5 will be generated.

# Read the hdf5 file
try:
    data = h5py.File("WinNet_data.h5","r")
except:
    print("Enable h5fc compiler in the Makefile! No hdf5 file was found!")
    sys.exit()

# Read the time
time = data["snapshots/time"][:]

# Read mass and proton number
A    = data["snapshots/A"][:]
Z    = data["snapshots/Z"][:]

# Read abundance
Y = data["snapshots/Y"][:,:]

# masks for some isotopes
ni56  = (A==56)  & (Z==28)
th232 = (A==232) & (Z==90)
u236  = (A==236) & (Z==92)
eu151 = (A==151) & (Z==63)
eu153 = (A==153) & (Z==63)

# Axis scaling
plt.yscale("log")
plt.xscale("log")

# Plotting
plt.plot(time,151.0*Y[:,eu151],label="$^{151}$Eu")
plt.plot(time,153.0*Y[:,eu153],label="$^{153}$Eu")
plt.plot(time,232.0*Y[:,th232],label="$^{232}$Th")
plt.plot(time,236.0*Y[:,u236] ,label="$^{236}$U")

# Make a legend
plt.legend()

# Axis labels
plt.ylabel("Mass fraction")
plt.xlabel("Time [s]")

# Axis limits
plt.ylim(1e-7,1e-2)
plt.xlim(1e-1,1e17)

plt.show()
