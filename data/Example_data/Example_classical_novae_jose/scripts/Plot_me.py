# Default skeleton for python files
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('../../bin')
from class_files.nucleus_multiple_class import nucleus_multiple

# Create a figure to determine the figure size
fig = plt.figure(figsize=(5,3))


# Read and plot reference final abundances
# A bit lengthy due to the file format
path = "final_abundances.txt"

abundances = []
nuclei     = []
with open(path,"r") as f:
    lines = f.readlines()

    mode = 'nuc'
    for l in lines[2:]:
        nuc = ""
        abu = ""
        for char in l:
            if mode == "nuc":
                if char == "=" and nuc.strip()!="":
                    mode = "abund"
                    nuc  = nuc.strip().replace(" ","")
                    # Take care of protons, deuterons, and tritium
                    if nuc.lower()=="prot":
                        nuc = "p"
                    elif nuc.lower()=="deut":
                        nuc = "d"
                    elif nuc.lower()=="trit":
                        nuc = "t"
                    elif (nuc.lower()[:4]=="al26"):
                        nuc = "al26"

                    # Take care of excited states of al26
                    if (nuc!="al26") or (not ("al26" in nuclei)):
                        nuclei.append(nuc)
                    else:
                        continue

                else:
                    nuc = nuc+char
            if mode== "abund":
                if ((char == " ") or (char == "\n")) and abu.strip()!="":
                    # Take care of excited states of al26
                    if nuc!="al26":
                        abundances.append(float(abu.replace("D","e")))
                    else:
                        idx = nuclei.index("al26")
                        if idx == len(abundances):
                            abundances.append(float(abu.replace("D","e")))
                        else:
                            abundances[idx] = abundances[idx] + float(abu.replace("D","e"))
                    abu = ""
                    nuc = ""
                    mode = "nuc"
                else:
                    if char!="=":
                        abu = abu+char

# Convert to arrays
abundances = np.array(abundances)
nuclei     = np.array(nuclei)

# Note that the abundances are mass fractions here
nm = nucleus_multiple(nuclei,X=abundances)
# Get the mass fractions summed over mass number
A,X = nm.A_X
plt.plot(A,X,label="Mean composition\n(Jose & Hernanz 1998, model ONe5)")

# Plot initial composition, again the file contains mass fractions
Z,A,X = np.loadtxt("iniab.txt",unpack=True,usecols=[0,2,3])
Z = Z.astype(int)
A = A.astype(int)
# Sum over equal A with the nucleus multiple class
nm = nucleus_multiple(Z=Z,A=A,X=X)
A,X = nm.A_X
plt.plot(A,X,label="Initial composition",ls="--")



A,Y,X = np.loadtxt("finabsum.dat",unpack=True)

plt.plot(A,X,label="WinNet, single trajectory")



plt.yscale("log")
plt.ylim(1e-8,1)
plt.xlim(0,55)

plt.legend(loc="upper right",bbox_to_anchor=(1.0,1.4))
plt.ylabel("Mass fraction")
plt.xlabel("Mass number")
plt.savefig('massfractions.pdf',bbox_inches='tight')
plt.show()