# Author: M. Reichert
import numpy as np
import matplotlib.pyplot as plt
# Import the WinNet python plot-routine
# Note that you can also just append the path to your .bashrc
import sys
import os
sys.path.append('../../bin')
from class_files.winnet_class import winnet

fig = plt.figure(figsize=(7,5))
ax  = fig.gca()


# Find all run folders
folders = np.array(os.listdir('.'))
folder_mask = list(map(lambda x: 'trajectory_' in x, folders))
folders = np.array(folders)[folder_mask]
single_run = winnet(folders[0])

# A large array that covers all possible nuclei
A_all = np.arange(260)
X_all = np.zeros(260)
# Loop over the folders
for f in folders:
    # Plot every single run and in the end an integrated abundance
    single_run = winnet(f)
    single_run.read_finab()
    A,X = single_run.get_final_A_X()
    X_all[A] += X
    ax.plot(A,X,color='lightgrey',zorder=1)

# Plot the integrated mass fractions
ax.plot(A_all,X_all/len(folders),color="tab:orange",zorder=3,label="Integrated mass fractions")
# Make the legend and a "fake" plot to get the label
ax.plot(np.nan,np.nan,color='lightgrey',zorder=1,label="Individual mass fractions")
ax.legend(loc="upper right")

# Set up the plot
ax.set_xlabel("Mass number")
ax.set_ylabel("Mass fraction")
ax.set_yscale("log")
ax.set_ylim(1e-6,1)
ax.set_xlim(49,245)
plt.savefig("massfractions_nsm_disc_wu.pdf",bbox_inches="tight")
plt.show()
