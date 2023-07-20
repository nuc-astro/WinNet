# Author: M. Reichert
import numpy as np
import matplotlib.pyplot as plt
# Import the WinNet python plot-routine
# Note that you can also just append the path to your .bashrc
import sys
sys.path.append('../../bin')
from class_files.winnet_class import winnet

fig,ax = plt.subplots(1,3,figsize=(10,4),sharex=True,sharey=True)

plt.subplots_adjust(hspace=0,wspace=0)



# Load data from the mainout
mrsn_example = winnet('.')

# We have to know the time of the snapshots
snapshot_time = [10,1.8E+03,8.64E+04] # In seconds. 10s, 30min, 1day

# Plot everything in the nuclear chart. Only the last ax object has the colorbar
mrsn_example.plot_nuclear_chart_at(snapshot_time[0],figure=ax[0],axes_label=False,element_labels=False,fig_is_ax=True,colorbar=False)
mrsn_example.plot_nuclear_chart_at(snapshot_time[1],figure=ax[1],axes_label=False,element_labels=False,fig_is_ax=True,colorbar=False)
mrsn_example.plot_nuclear_chart_at(snapshot_time[2],figure=ax[2],axes_label=False,element_labels=False,fig_is_ax=True,colorbar=True,
                                   colorbar_position=[0.27, 0.85, 0.5, 0.025],colorbar_inset=True)

ax[0].set_ylabel("Proton number")
fig.text(0.5,0.15,"Neutron number",ha="center")

ax[0].set_xlim(9,59)
ax[0].set_ylim(15,49)
plt.savefig("massfractions_chart_nup_obergaulinger.pdf",bbox_inches="tight")
plt.show()
