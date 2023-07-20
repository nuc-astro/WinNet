# Author: M. Reichert
import numpy as np
import matplotlib.pyplot as plt
# Import the WinNet python plot-routine
# Note that you can also just append the path to your .bashrc
import sys
sys.path.append('../../bin')
from class_files.winnet_class import winnet

# Create the figure objects
fig  = plt.figure(figsize=(8,4))
ax   = fig.gca()

path = "bin_"

# Plot all bins
for i in range(1,5):
    # Create the proper path
    p_tmp = path+str(i)+".dat"
    # Read the final abundances
    run   = winnet(p_tmp)
    run.read_finab()
    A,X   = run.get_final_A_X()

    plt.plot(A,X,label="Bin "+str(i))


plt.ylim(1e-8,1)
plt.legend()
plt.xlabel("Mass number")
plt.ylabel("Mass fraction")
plt.yscale("log")
plt.savefig("final_massfractions.pdf",bbox_inches="tight")
plt.show()
