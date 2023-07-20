# Author: M. Reichert
# Import packages
import numpy             as np
import matplotlib.pyplot as plt
# Determining the trajectory of standart Big Bang.
# We need to know:
# - Temperature evolution
# - Density evolution
# - Electron fraction at weak freezeout

# The trajectory is derived for a flat isotropic and homogeneous universe, with adiabatic expandion, following (Winteler 2013).
# The Friedmann equations for this are given by:
#
# (1):   d^2R/dt^2 = -4pi//(3c^2) * (rho_eps + 3P) GR + 1/(3c^2) Lambda R
# (2):   (dR/dt)^2 = H(t)^2 = 8piG/(3c^2) rho_eps - kc^2/(R^2) + 1/(3c^2) Lambda
# (3):   d(rho_eps R^3)/dt + P dR^3/dt = 0
# (for a flat universe k=0 and Lambda = 0)


# At aproximately kb T_weak = 0.8 MeV [Winteler 2013, Fields 2006] the neutron to proton ratio is frozen out. From Maxwell Boltzmann distribution and the
# equilibrium of chemical potentials one can derive n/p = 1/6 (with n+p = 1 => n=1/7 and p=1/6)

# Constants
mn          =  939.5654133              #  mass neutron    [ MeV / c**2 ]
mp          =  938.2720813              #  mass proton     [ MeV / c**2 ]
Delta       = mn - mp                   #  mass difference [ MeV / c**2 ]

# Assume weak freezout at T=0.8 MeV
kB_Tweak = 0.8                          # [MeV]
n_p_ratio = np.exp(-Delta/kB_Tweak)     # Followed from Maxwell Boltzmann expression for the chemical potentials

# From n+p = 1:
neutron_freezout = n_p_ratio/(1. + n_p_ratio)
proton_freezout  = 1. - neutron_freezout
# The initial electron fraction is therefore given as:
ye_freezout = proton_freezout/(neutron_freezout+proton_freezout) # Same as neutron abundance, because there are only nucleons

# Calculate the first time after freezeout.
# MeV -> GK: 11.604519
first_time = (13.336 * 1./(kB_Tweak * 11.604519))**2.
# Get the time
time = np.logspace(np.log10(first_time),5,num=200) # [s]

# From adiabatic expansion and from relativistic degrees of freedom g = 3.3626:
temperature = 13.336 * 1./(time**0.5) # [GK]

# The density can be calculated in terms of the photon to baryon ratio eta:
# (Change eta for creating the plot shown in Schramm, Turner 1998)
eta = 6.11e-10                  # Measurement of WMAP
# Define some other constants:
m_u   = 1.660540e-27 * 1e3      # Weight of a nucleon
pi    = 3.14159265359
kB    = 8.6173303e-2            # [MeV / GK]
hbarc = 197.3269718 * 1.e-13    # [MeV * cm]
# Calculate the density in [g/ccm]. (From adiabatic expansion and assuming ultra relativistic particle [see Winteler 2012])
dens  = 2.404 * eta/(pi**2.) * ((kB*temperature)/(hbarc))**3. * m_u # [g/ccm]

# Create ye column. Note that the network only uses the first value, so we can keep it constant
ye           = [ye_freezout for i in range(len(time))]

# Save the trajectory
out = np.array([time,temperature,dens,ye]).T
np.savetxt('bbn.dat',out,header='time[s], T9[GK], density[g/cm^3], Ye \n')
