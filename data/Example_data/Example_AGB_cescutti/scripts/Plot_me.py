# Author: M. Reichert
import numpy as np
import matplotlib.pyplot as plt


# We note that small differences between the result here and
# in Cescutti et al. 2018 can originate from different reaction
# rate input.


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


# Get the result of WinNet and plot it
A,Y,X = np.loadtxt('trajectory_C13pocket/finabsum.dat',unpack=True)
plt.plot(A,X,label="WinNet")

# Read and plot the result from Cescutti et al. 2018
A,Yfinal,Yseed = np.loadtxt('final_abundances_C13pocket.out',unpack=True,usecols=[3,4,6])
Xfinal         = Yfinal*A
Asum,Xsum      = sum_over(A,Xfinal)
plt.plot(Asum,Xsum,label="Cescutti et al. 2018")

# Plot the initial mass fractions
Xseed          = Yseed*A
Asum,Xsum      = sum_over(A,Xseed)
plt.plot(Asum,Xsum,label="Initial mass fractions")

# Also plot the result of the thermal pulse
A,Y,X = np.loadtxt('trajectory_TP/finabsum.dat',unpack=True)
plt.plot(A,X,label="WinNet (after thermal pulse)")

# Give labels to the plot
plt.ylabel("Mass fractions")
plt.xlabel("Mass number")

# Make it more beautiful
plt.yscale("log")
plt.ylim(1e-11,1)
plt.xlim(0,215)
plt.legend()

# Save and show the plot
plt.savefig("main_s_process.pdf",bbox_inches="tight")
plt.show()
