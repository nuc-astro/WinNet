# Example file for nucleus class
# Author: Moritz Reichert
# Date  : 04.05.17

# Import the nucleus class (Attention, the import has to be different for python 3)
# This in only working in python 2!
from class_files.nucleus_class       import nucleus

# Simple class to work with nuclei. The class is able to analyze a nucleus from a
# name of a nucleus and extract several properties from this name.



################
# INITIALIZING #
################

# Create instances of the nucleus class
o16 = nucleus('o16')
ne15= nucleus('ne15')
he4 = nucleus('he4')


##############
# PROPERTIES #
##############


# Get the properties of a nucleus
mass            = o16.get_A()
atomic_number   = o16.get_Z()
neutron_number  = o16.get_N()
name            = o16.get_name()
# Print it
print(name+' has a mass number of '+str(mass)+', with Z = '+str(atomic_number)+' and N = '+str(neutron_number)+'.')

# The nucleus object can also store an abundance
o16.set_Y(1e-10)
# or alternatively set the mass fraction with o16.set_X(1e-10)
# You can access it with o16.get_Y()



###########
# SORTING #
###########


# Put all nuclei in a list
nuclei_list = [ne15,he4,o16]

# Sort it corresponding to proton number
sort_Z = sorted(nuclei_list)
sort_Z_names = map(lambda x: x.get_name(),sort_Z)
print('Sorted for Z: ',sort_Z_names)

# Sort it for A
# First change the sort criteria for every nucleus in list
for n in nuclei_list:
    n.set_sortcriteria('A')

sort_A = sorted(nuclei_list)
sort_A_names = map(lambda x: x.get_name(),sort_A)
print('Sorted for A: ',sort_A_names)
