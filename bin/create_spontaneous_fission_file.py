# Author: Moritz Reichert
# Date: 15.03.2023
#
# This file will create a file that contains entries for
# spontaneous fission reactions. These entries have to be copied
# into the fission rates file (e.g., fissionrates_frdm).
# To create the file here we use the fission barriers from
# Moeller et al. 2015 https://journals.aps.org/prc/supplemental/10.1103/PhysRevC.91.024310 and
# the FRDM fit of Khuyagbaatar 2020 (https://ui.adsabs.harvard.edu/abs/2020NuPhA100221958K/abstract)
# In contrast to Khuyagbaatar 2020, we also include odd nuclei, not only even ones.
# Additionally we add experimentally determined rates that can be accessed from the ENDFS database
# nds.iaea.org/relnsd/v1/data?/fields=ground_states&nuclides=all.
# It is in principle possible to also create spontaneous fission rates for the ETFSI model, but the
# here used fit has to be changed and the correct fission barriers have to be used.
import numpy       as np
import pandas      as pd
from nucleus_class import nucleus


# Path to the fission barriers.
path_barriers      = 'fission_barriers_frdm.dat'

# Path to experimental spontaneous fission half-lives
path_exp_halflifes = 'sf_halflife.csv'

# Read the fission barriers of frdm. They can be accessed at
# https://journals.aps.org/prc/supplemental/10.1103/PhysRevC.91.024310
Z,N,A,B=np.loadtxt(path_barriers,unpack=True)


# Read the experimental spontaneous fission half-lives.
df = pd.read_csv(path_exp_halflifes)
df['sf_halflife'] = pd.to_numeric(df['sf_halflife'],errors='coerce')
df = df[np.isfinite(df['sf_halflife'])]
hl_exp = df['sf_halflife'].values


def T05(Z,N,B):
    """
      Calculate the half-life as in Khuyagbaatar 2020.
      The link to the paper is the following:
      https://ui.adsabs.harvard.edu/abs/2020NuPhA100221958K/abstract
      Here, we use the fit of the FRDM data (Eq. 3) in combination with
      Eq. 1 and 2.
    """
    A = Z+N
    hqueromega = 0.14025*Z**2/A-4.6335
    T = (1+np.exp(np.pi*2*B/hqueromega))**(-1)
    hquer = 6.582119514e-16
    n = hqueromega/hquer/(2*np.pi)
    # n=1e14
    Tsf = np.log(2)/(n*T)
    # Negative values can happen due to hqueromega being negative.
    # These cases have an infinite half-life.
    if Tsf<0:
        Tsf = np.inf
    return Tsf



# Create lists for the half-lifes and nuclei names
t05 = []
nuc = []

# Proton number range for spontaneous
# fission half-lives.
zmin = 90
zmax = 112

# Loop through the half-lives and create a new array
# containing them. Put experimental ones whenever possible.
for ind,ztmp in enumerate(Z):
    # Skip irrelevant nuclei
    if (ztmp>zmax) or (ztmp <zmin):
        continue
    # Get the name of the nucleus
    nuc.append(nucleus(Z=int(Z[ind]),N=int(N[ind])).get_name())
    # Check if experimental half-life is available
    if np.any((df["z"]==ztmp) & (df["n"]==N[ind])):
        ht = hl_exp[(df["z"]==ztmp) & (df["n"]==N[ind])][0]
        # Check if the value is not nan
        if not np.isnan(ht):
            t05.append(float(ht))
        else:
            t05.append(T05(ztmp,N[ind],B[ind]))
    else:
        t05.append(T05(ztmp,N[ind],B[ind]))


# Create the reaclib file
out = ''
# Maximum half-life to be considered
# 1e25s = 3.17e17 years
thalf_max = 1e25
# Minimum cut-off half-life for too short half-lives
# Note that too fast rates will mess up the network equations
# and lead to a crash
thalf_min = 1e-15

# Loop through the half-lives and create the reaclib file
for ind,n in enumerate(nuc):
    if t05[ind]>thalf_max:
        continue
    if t05[ind]<thalf_min:
        t05[ind]=thalf_min

    # Calculate reaclib parameter
    a0 = np.log(np.log(2)/t05[ind])
    # Create the entry in correct format
    out += '     '+n.rjust(5)+33*' '+'sfis'+4*' '+'{:13.6e}'.format(0)+'\n'
    out += '{:13.6e}'.format(a0)+3*'{:13.6e}'.format(0)+'\n'
    out += 3*'{:13.6e}'.format(0)+'\n'


# Write the file
with open("frdm_sfis.dat","w") as f:
    f.write(out)