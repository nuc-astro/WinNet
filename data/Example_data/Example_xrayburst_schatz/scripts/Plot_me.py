# Author: M. Reichert
import numpy as np
import matplotlib.pyplot as plt


def sum_over(A,X):
    """
     Sum mass fractions over equal mass numbers
    """
    max_A = max(A)
    X_new = np.zeros([max_A+1,])
    for i in range(len(A)):
        X_new[int(A[i])] += X[i]
    return np.array(range(max_A+1)),X_new

# Load the final abundances of SkyNet (without calculating detailed balance there)
A,Z,N,Y,X = np.loadtxt("finab_SkyNet.dat",unpack=True)

# Sum over equal A's
A_summed, X_summed = sum_over(A.astype(int),X)

# Plot SkyNet
plt.plot(A_summed,X_summed,label="SkyNet")

# Load WinNet
A,Z,N,Y,X = np.loadtxt("finab.dat",unpack=True)

# Sum over equal A's
A_summed, X_summed = sum_over(A.astype(int),X)

# Plot WinNet
plt.plot(A_summed,X_summed,label="WinNet")

# Show the legend
plt.legend()
plt.xlabel("Mass number A")
plt.ylabel("Mass fraction X")
plt.ylim(1e-10,1)
plt.yscale("log")
plt.savefig("comparison_SkyNet_WinNet.pdf",bbox_inches="tight")
plt.show()
