# Author: M. Reichert
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# Script to create tabulated rates for WinNet

# The default talys temperature grid as used in WinNet
talys_t = np.array([1.0e-4,5.0e-4,1.0e-3,5.0e-3,1.0e-2,5.0e-2,1.0e-1,1.5e-1,2.0e-1,2.5e-1,
                    3.0e-1,4.0e-1,5.0e-1,6.0e-1,7.0e-1,8.0e-1,9.0e-1,1.0e+0,1.5e+0,2.0e+0,
                    2.5e+0,3.0e+0,3.5e+0,4.0e+0,5.0e+0,6.0e+0,7.0e+0,8.0e+0,9.0e+0,1.0e+1 ])

# Define the string that will contain the tabulated rates in correct format
# as defined by WinNet
rate_file = ""



# Interpolate the tabulation from the paper on the talys grid
# Here its Ne22(a,g)Mg26
T9,adopt,exp = np.loadtxt("ne22amg26_angulo1999",unpack=True,skiprows=1,)
rate = adopt*10**exp
f = interp1d(np.log10(T9),np.log10(rate),fill_value=-99,bounds_error=False,kind="cubic")
rates = 10**f(np.log10(talys_t))
# If you want you can plot:
# plt.plot(talys_t,rates)

# Write it into the string (first reaclib chapter, then rate)
rate_file += "4\n\n"
rate_file += "       he4 ne22 mg26                       an99      1.06148e+01          \n"
stringlist = ["{:.3e}".format(ytmp) for ytmp in rates ]
rate_file += " "+" ".join(stringlist)+"\n"

# The next rate, O17(a,g)Ne21
# This rate is already in reaclib! However, they multiply it by 0.1 in Nishimura et al. 2017
T9,adopt,exp = np.loadtxt("o17agne21_cf88",unpack=True,skiprows=1,)
rate = adopt*10**exp
rate[rate==0]=1e-98
f = interp1d(np.log10(T9),np.log10(rate),fill_value=-98,bounds_error=False,kind="cubic")
rates = 10**f(np.log10(talys_t))
# If you want you can plot:
# plt.plot(talys_t,rates)

# Write the rate into the string again
rate_file += "       he4  o17 ne21                       cf88n     7.35100e+00          \n"
stringlist = ["{:.3e}".format(ytmp*0.1) for ytmp in rates ]
rate_file += " "+" ".join(stringlist)+"\n"


# The Ne22(a,n)Mg25 is given in a slightly different format in Jaeger et al. 2001
a = np.array([4.04,2.302e-4,6900,1.881e7])
b = np.array([0,-0.6,3.19,0.358])
c = np.array([7.74,6.14,11.3,26.7])
func = lambda T: np.sum(a*T**b*np.exp(-c/T))
y = np.array([func(ttmp) for ttmp in talys_t])
y[y<1e-99]=1e-99
# If you want you can plot:
# plt.plot(talys_t,y)

# Write the rate into the string again
rate_file += "5\n\n"
rate_file += "       he4 ne22    n mg25                  ja01     -4.78290e-01          \n"
stringlist = ["{:.3e}".format(ytmp) for ytmp in y ]
rate_file += " "+" ".join(stringlist)+"\n"


# The last rate: O17(a,n)Ne20
T9,adopt,exp = np.loadtxt("017anne20_angulo1999",unpack=True,skiprows=1,)
rate = adopt*10**exp
f = interp1d(np.log10(T9),np.log10(rate),fill_value=-99,bounds_error=False,kind="cubic")
rates = 10**f(np.log10(talys_t))
rates[rates<1e-99]=1e-99
# If you want you can plot:
# plt.plot(talys_t,rates)

# Write again to the string
rate_file += "       he4  o17    n ne20                  an99      5.86000e-01          \n"
stringlist = ["{:.3e}".format(ytmp) for ytmp in rates ]
rate_file += " "+" ".join(stringlist)+"\n"

# Write the string to the file
with open("tab_rates.dat","w") as f:
    f.write(rate_file)


# In case you plotted one of the rates:
# plt.yscale("log")
# plt.xscale("log")
# plt.xlim(1e-2,1e1)
# plt.ylim(1e-50,1e11)
# plt.show()
