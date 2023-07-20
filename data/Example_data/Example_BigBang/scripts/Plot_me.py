# Author: M. Reichert
import matplotlib.pyplot    as plt
import numpy                as np

# Plot the nuclei:    n    p    d    t  he3  he4  li6  li7  be7
# The mass numbers are:
A     = np.array([1,1,2,3,3,4,6,7,7])
names = ['n','p','d','t',r'$^3$He',r'$^4$He',r'$^6$Li',r'$^7$Li',r'$^7$Be']

# Load the data
data = np.loadtxt("tracked_nuclei.dat",unpack=True)
# Time [s]
time = data[0]
# Abundance
Y    = data[1:].T
# Mass fraction
X    = Y * A

# Set up the plot
fig = plt.figure()
ax  = fig.gca()
colors = ['g','b','c','k','y','m','r','navy','teal']

# Plot the mass fraction of every nucleus
for i in range(len(X.T)):
    ax.plot(time,X.T[i],label=names[i],color=colors[i])

# Set the style of the plot
ax.set_xlabel('Time [s]')
ax.set_ylabel('Mass fraction')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(time[0],time[-1])
ax.set_ylim(1e-14,1.1)
plt.legend(loc="upper center",bbox_to_anchor=(0.5, 1.17),ncol = 5,fancybox=True)

# Save the plot
fig.savefig("bbn.pdf",bbox_inches='tight')

# Show the plot
plt.show()
