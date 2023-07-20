# Author: M. Reichert
import numpy as np
import matplotlib.pyplot as plt
# Import the WinNet python plot-routine
# Note that you can also just append the path to your .bashrc
import sys
import os
sys.path.append('../../bin')
from class_files.nucleus_multiple_class import nucleus_multiple

def sum_over(A,Y):
    """
      Function to sum up mass fractions or abundances over equal A.
    """
    A = A.astype(int)
    max_A = max(A)
    # second_dimension = len(Y[0,:])
    Y_new = np.zeros([max_A+1,])
    for i in range(len(A)):
        Y_new[int(A[i])] += Y[i]
    return np.array(range(max_A+1)),Y_new

# Set up two figures
fig = plt.figure(figsize=(5,2.5))
ax  = fig.gca()

fig2 = plt.figure(figsize=(5,2.5))
ax2  = fig2.gca()

# Read and plot initial mass fractions
A,X = np.loadtxt("../../data/Example_data/Example_AGB_nishimura/iniab1.4E-02As09.ppn",
                 unpack=True,usecols=[2,3])
Asum,Xsum = sum_over(A,X)
ax2.plot(Asum,Xsum,label="Initial")

# Read the result
A,Z,N,Y,X = np.loadtxt("finab.dat",unpack=True,usecols=[0,1,2,3,4])
# Convert to integers
A = A.astype(int)
Z = Z.astype(int)
# Create a nucleus multiple instance. This is a class to
# deal with abundances.
nm_winnet = nucleus_multiple(A=A,Z=Z,Y=Y)

Asum,Xsum = nm_winnet.A_X
ax2.plot(Asum,Xsum,label="Final")

# Read the solar abundances
Z,A,Y=np.loadtxt("solar_abu_lodders.txt",unpack=True)
Z = Z.astype(int)
A = A.astype(int)
# Also put them in a nucleus multiple class
nm_solar = nucleus_multiple(A=A,Z=Z,Y=Y)

Asum,Xsum = nm_solar.A_X
ax2.plot(Asum,Xsum,label="Solar",lw=0.5)

# Calculate overproduction
nm_overprod= nm_winnet/nm_solar

# Get unique amount of protons to loop through
Z = np.unique(nm_overprod.Z)

# List with the elementnames, at place Z is the name of the element, e.g., elementnames[1]="h"
elementnames = ('neutron','h','he','li','be','b','c','n','o','f','ne','na','mg','al','si','p','s','cl','ar','k','ca','sc','ti','v','cr','mn','fe',
  'co','ni','cu','zn','ga','ge','as','se','br','kr','rb','sr','y','zr','nb','mo','tc','ru','rh','pd','ag','cd','in','sn','sb',
  'te', 'i','xe','cs','ba','la','ce','pr','nd','pm','sm','eu','gd','tb','dy','ho','er','tm','yb','lu','hf','ta','w','re','os',
  'ir','pt','au','hg','tl','pb','bi','po','at','rn','fr','ra','ac','th','pa','u','np','pu','am','cm','bk','cf','es','fm','md',
  'no','lr','rf','db','sg','bh','hs','mt','ds','rg','ub','ut','uq','up','uh','us','uo')

# Loop through it, plot the overproduction chains
# and write the correct name at the correct place
for Ztmp in Z:
    Y = nm_overprod.Y[(nm_overprod.Z==Ztmp) & (nm_overprod.Y>=1e-2)]
    A = nm_overprod.A[(nm_overprod.Z==Ztmp) & (nm_overprod.Y>=1e-2)]
    # The "try" is necessary for empty sequences
    try:
        maxY_arg = np.argmax(Y)
    except:
        continue
    ion_name = elementnames[Ztmp][0].upper()+elementnames[Ztmp][1:]
    ax.scatter(A,Y,s=10)
    p = ax.plot(A,Y)
    ax.text(A[maxY_arg],Y[maxY_arg]*1.1,ion_name,clip_on=True,ha="center",color=p[0].get_color())

# Make lines as in the Nishimura paper
ax.axhline(1e-1,color="lightgrey",zorder=-1,lw=0.8,ls="--")
ax.axhline(1e0 ,color="lightgrey",zorder=-1,lw=0.8,ls="dotted")
ax.axhline(1e1 ,color="lightgrey",zorder=-1,lw=0.8,ls="--")

# Set the limits and scales
ax.set_ylim(1e-2,1e3)
ax.set_xlim(40,160)
ax.set_yscale("log")
# Write the labels
ax.set_ylabel(r"Overproduction X/X$_\odot$")
ax.set_xlabel("Mass number")

ax2.legend()
ax2.set_yscale("log")
ax2.set_ylim(1e-11,1)
ax2.set_xlim(0,220)
ax2.set_xlabel("Mass number")
ax2.set_ylabel("Mass fraction")

# Save and show the plots
fig.savefig("overproduction.pdf",bbox_inches="tight")
fig2.savefig("mass_fractions.pdf",bbox_inches="tight")
plt.show()
