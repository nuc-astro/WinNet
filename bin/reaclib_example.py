# Example file for reaclib class
# Author: Moritz Reichert
# Date  : 08.05.17

# Import the reaclib class
from class_files.reaclib_class  import reaclib
import matplotlib.pyplot        as plt
import pandas                   as pd

# The class is able to make several tests with a reaclib file. It is also able to plot rates of a reaclib.
# Furthermore, you can update the rates with other reaclibs, specifying what you want to update.
# This should be done carefully and always checked afterwards!!

################
# INITIALIZING #
################

# Create an instance of the reaclib class. If you want the class to be quiet, add quiet=True.
reaclib_file_path = " INSERT REACLIB FILE PATH "
reaclib_file_update_path = " INSERT MORE RECENT REACLIB FILE PATH "

r_testclass = reaclib(reaclib_file_path)
# Read the file
r_testclass.read_reaclib()


#########
# TESTS #
#########

# Make several tests with the file
# Tests are:
# 1. Test for rate being in the correct chapter
# 2. Test for consistency of weak rates
# 3. Test for mass consistency
# 4. Test for overflows of the rate
r_testclass.test_reaclib()
# Save the result to "reaclib_errors.html"
r_testclass.get_rate_error_html('reaclib_errors.html')
# Remove all reactions that contain errors
r_testclass.drop_errors()


############
# PLOTTING #
############

reactants = ['he4','c12']
products  = ['o16']
# Also kwargs for plotting are possible
# Plot it (optional values for this function are figure=None (give a figure object to plot to the same axes ),axlabel=True (label the x and y axis))
figure = r_testclass.plot_rate(reactants,products)
# Get the axes object
# you could change the range by e.g. ax.set_ylim(1e-50,1e10)
ax     = figure.gca()
# Show it
plt.show()


#######################
# UPDATE AND ANALYSIS #
#######################

# Get a brief overview of the reactions. How much n-gamma, gamma-n, alpha-n,... reactions are included. Return value is a pandas series
print('Overview of contained reactions:')
print((r_testclass.get_statistics()))

# Update reaclib with other more recent one.
# Possible values are:
# reaclib_path  - path of the reaclib that contains new rates
# dataframe     - pandas dataframe that contains reactions. If this is given, you don't have to give a reaclib_path.
# reaction_type - Reaction type that will be updated (string or list of strings), for possible values see ".get_statistics()"
# possible reaction types are: 'n-gamma', 'gamma-n', 'p-gamma', 'gamma-p', 'a-gamma', 'gamma-a', 'n-p', 'p-n', 'n-a', 'a-n', 'p-a', 'a-p', 'weak', 'fission', 'other'
# chapter       - Chapter that will be updated (list of integer or integer)
# ignore_label  - ignore the label for updating and only look at the reactions itself
# ignore_reverse- ignore if the reaction is reverse
# ignore_type   - ignore the type
# r_testclass.update('testfiles/20161205ReaclibV2.0no910',reaction_type=['n-gamma','gamma-n'])
r_testclass.update(reaclib_file_update_path,ignore_label=False)
print((r_testclass.get_statistics()))
# Save the reaclib again
r_testclass.save_reaclib('updated_reaclib.dat')

# It's also possible to get all 'n-gamma' rates and play around with them to merge them in afterwards again
n_gamma = r_testclass.get_dataframe(reaction_type='n-gamma')
# You can also print them: (comment it in again)
# print(n_gamma)
# Or multiply factor 10 to all a0's in all n-gamma reactions
n_gamma['a0'] = n_gamma['a0'].apply(lambda x: x*10.)
# And merge it in again.
r_testclass.update(dataframe=n_gamma)

# Or get critical low temperature rates:
# with min_temperature range of temperature
# amount_points - points for grid to look if the rate diverges
# max_rate - at which value is the rate critical?
crit_low_temp_dataframe = r_testclass.get_critical_low_temperature_rates(min_temperature=1e-3,amount_points=20,max_rate=1.e100)
# And do whatever you want with it e.g. save it or print the amount of rates
print(('Amount critical low temperature rates : ' + str(crit_low_temp_dataframe.count()[0])))
r_testclass.save_reaclib('crit_lowtemp.dat',dataframe = crit_low_temp_dataframe)
