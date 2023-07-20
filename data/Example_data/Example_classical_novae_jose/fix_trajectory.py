# Author: M. Reichert
# Date  : 21.03.2023
import numpy as np


# The trajectory downloaded from https://zenodo.org/record/6474694
# DOI: 10.5281/zenodo.6474694 , Version 1.1.1
# includes the same time with multiple thermodynamic entries.
# This non-monotonic time is not supported by the code and gets
# fixed with the following script:
path = "trajectory_ONe1p25nova"
time,temp,dens = np.loadtxt(path, skiprows=7,unpack=True)
mask = time[1:]-time[:-1] != 0
mask = np.append(mask, True)
out = np.array([time[mask],temp[mask],dens[mask]]).T
np.savetxt('trajectory_ONe1p25nova_mon_time', out, header="time [yrs] temperature [GK] density [g/cm^3]", fmt="%.16e %.6e %.6e")