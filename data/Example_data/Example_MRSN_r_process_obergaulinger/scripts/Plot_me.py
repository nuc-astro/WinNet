# Author: M. Reichert
import numpy as np
import matplotlib.pyplot as plt
# Import the WinNet python plot-routine
# Note that you can also just append the path to your .bashrc
import sys
sys.path.append('../../bin')
from class_files.winnet_class import winnet

fig = plt.figure(figsize=(10,5))
ax  = fig.gca()

# Load data from the mainout
mrsn_example = winnet('.')

mrsn_example.read_mainout()
time   = mrsn_example.get_mainout_time()
yheavy = mrsn_example.get_mainout_yheavy()
yn     = mrsn_example.get_mainout_yn()

# Calculate freezeout point#
# Assume that the output is frequent enough
index_freezout = np.argmin(abs(yn/yheavy-1))
time_freezout  = time[index_freezout]

# Plot the path of abundances in the nuclear chart at freezeout
# This only works if snapshots were enabled!!
anim = mrsn_example.animate_nuclear_chart(figure=fig,plot_magic=True,time_title=True,min_X=1e-8,element_labels=False)
# To make only one snapshot at freezeout you can plot it with:
# mrsn_example.plot_nuclear_chart_at(time_freezout,figure=fig,plot_magic=True,time_title=True,min_X=1e-8,element_labels=False)


ax.set_xlim(40,130)
ax.set_ylim(25,75)
plt.show()
