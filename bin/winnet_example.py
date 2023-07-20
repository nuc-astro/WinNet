# Example file for winnet class
# Author: Moritz Reichert
# Date  : 08.05.17

# Import winnet class
from class_files.winnet_class   import winnet
import matplotlib.pyplot        as plt
import pandas                   as pd

# The class is able to plot several things and return the data.

################
# INITIALIZING #
################

# Create an instance of the reaclib class. If you want the class to be quiet, add quiet=True.
w_testclass = winnet('testfiles/winnet_testfiles/example_run')
# There is a possibility to read a "First started" guide:
print(w_testclass.get_man())
# As told, you first have to read data. For simplicity our test folder only contains the finabs, mainout and the param.out, but in principle the class is able to read nearly all data.
# possible read values are given by:
# w_testclass.read_finab()
# w_testclass.read_mainout()
# w_testclass.read_snapshots()
# w_testclass.read_timescales()
# w_testclass.read_template()
# Or if you want to read the whole run w_testclass.read_run(). This might complain that some things are not there
w_testclass.read_finab()
w_testclass.read_template()
w_testclass.read_mainout()




#############################
# PLOTTING AND GETTING DATA #
#############################

# The path for the sunet is extracted from param.out
# We can plot the contained nuclei in the run
w_testclass.plot_sunet()
# Show the content
plt.show()

# or plot more normal things like mass number over mass fraction
w_testclass.plot_final_A_X()
# Show the content
plt.show()

# or plot the different isotopes
# This function has several options:
# lower_limit             : lowest massfraction that will be plotted (default 10^-10)
# isotopes                : isotopes to be plotted. If none (default), all isotopes of finab will get plotted. Type is list that contains the Z number
# isotope_label           : Plot the labels of the isotopes?
# ignore_isotope_names    : which name label should be ignored? (list contains Z number)
# text_clipping           : Text clipping yes or no? ( Should the text shown outside the plotting area? )
# ytext_pos_factor        : offset factor for the y position of the labels. Only has an effect when isotope_label is true
# xtext_pos_offset        : offset for the x position of the labels. Only has an effect when isotope_label is true
w_testclass.plot_final_isotopes()
plt.show()

# or just show different quantities of the mainout and look for crazy things
w_testclass.plot_mainout()
plt.show()

# Useful can be the time evolution of the nuclei from snapshots. Because there are no snapshots in the example folder, we can not show this here.
# in principle you have to read the snapshots first
# w_testclass.read_snapshots()
# and you will have access to the data afterwards through (e.g. for he4)
# w_testclass.get_time_evolution('he4')


# There is much more included. A list with all possible routines can be accessed by
print(w_testclass.get_methods())
