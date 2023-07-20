# Author: M. Reichert
# Date: 09.07.22
import numpy as np
import matplotlib.pyplot as plt
# Import the WinNet python plot-routine
# Note that you can also just append the path to your .bashrc
import sys
sys.path.append('../../bin')
from class_files.winnet_class import winnet

# This script will plot the mass fractions and energy generation
# of CO burning. The plot can be compared to the one on cococubed:
# https://cococubed.com/code_pages/burn_helium.shtml

# Set up the figure
fig = plt.figure()
ax  = fig.gca()
ax.set_title("Hydrostatic CO burning")


def create_nucleus_name(s):
    """
      Function to convert the element string that is given by WinNet
      to a latex style string.
      Example:
        create_nucleus_name("ne20")
      will output: r"$^{20}$Ne"
    """
    digits =""
    chars  =""
    for d in s:
        if d.isdigit():
            digits+=d
        else:
            chars+=d
    nuc_name = r"$^{"+digits+"}$"+chars[0].upper()+chars[1:]
    return nuc_name

# Read the mass fractions of the elements
w = winnet(".")
w.read_tracked_nuclei()
names = w.get_tracked_nuclei_names()
time  = w.get_tracked_time()

# Define the x-position of each element label
xpos_dic = {"he4":4e-7,"c12":3.8e-5,"o16":1e-7,"ne20":5e-3,"mg24":4e2,"si28":4e2,\
            "s32":4e2,"ar36":4e2,"ca40":5e9,"ti44":4e2,"cr48":5e6,"fe52":1e10,"ni56":1e10}

for n in names:
    # Plot the mass fractions of each element
    X = w.get_tracked_nuclei(n)
    line, = ax.plot(time,X)
    # Plot the element labels
    if n in xpos_dic:
        out_name = create_nucleus_name(n)
        idx = np.argmin(abs(time-xpos_dic[n]))
        ax.text(xpos_dic[n],X[idx],out_name,ha="left",va="bottom",color=line.get_color())

# Plot also the generated energy on a second y-axis
ax_energy = ax.twinx()
time_energy,energy = np.loadtxt("generated_energy.dat",unpack=True,usecols=[0,1])
ax_energy.plot(time_energy,energy,ls="--",lw=1,color="saddlebrown",zorder=-30,alpha=0.8)
ax_energy.text(5e-9,2e25,r"$\epsilon$",color="k")
ax_energy.set_yscale("log")
ax_energy.set_ylabel(r"Energy [erg g$^{-1}$ s$^{-1}$]")
ax_energy.set_ylim(1e11,1e27)

# Set the labels, limits, and scales
ax.set_xlabel("Time [s]")
ax.set_ylabel("Mass fraction")
ax.loglog()
ax.set_xlim(1e-11,1e12)
ax.set_ylim(1e-8,1e1)
plt.savefig("CO_burning.pdf",bbox_inches="tight")
plt.show()
