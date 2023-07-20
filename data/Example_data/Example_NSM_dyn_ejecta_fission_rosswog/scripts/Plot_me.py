# Author: M. Reichert
import numpy as np
import matplotlib.pyplot as plt

# Set up the figure with a reasonable size
fig = plt.figure(figsize=(5,3))

# The runs, the numbers correspond to different fissflags
# Fissflag 1: Panov et al. 2001
# Fissflag 2: Kodama & Takahashi 1975
# Fissflag 3: Mumpower et al. 2020 for neutron-induced and beta-delayed fission,
#             Kodama & Takahashi 1975 for spontaneous fission
runs=["1","2","3"]
labels=["Panov et al. 2001","Kodama & Takahashi 1975", "Mumpower et al. 2020"]

# Loop over the runs and plot the data
for ind,r in enumerate(runs):
    path = r+"/finabsum.dat"
    A,Y,X = np.loadtxt(path,unpack=True)
    plt.plot(A,X,label=labels[ind])

# Make the plot look nice
plt.yscale("log")
plt.ylim(1e-8,1)
plt.xlim(40,220)
plt.legend(loc="lower right",framealpha=1)
plt.xlabel("Mass number")
plt.ylabel("Mass fraction")
# Save the plot
plt.savefig("different_fission_fragments.pdf",bbox_inches="tight")
# Show the plot
plt.show()