# Author: M. Reichert
import matplotlib.pyplot as plt
import numpy as np
import os

# Plot idea from https://wmap.gsfc.nasa.gov/universe/bb_tests_ele.html

# Create Figure
fig = plt.figure(figsize=(5,5))

# Find all run folders
folders = np.array(os.listdir('.'))
folder_mask = list(map(lambda x: 'bbn_' in x, folders))
folders = np.array(folders)[folder_mask]

# Create and sort different baryon to photon ratios
etas = np.array(list(map(lambda x: float(x.replace('bbn_','')), folders)))
key = np.argsort(etas)

folders = folders[key]
etas    = etas[key]

# Prepare lists to store the final abundances for each run
he_list = []
hydrogren_list = []
d_list  = []
he3_list  = []
li7_list  = []
for f in folders:
    path = f
    # 1:A 2:Z 3:N 4:Yi 5:Xi
    A,Z,N,Y,X = np.loadtxt(path+'/finab.dat',unpack=True)

    # Store the abundances
    hydrogen = Y[(A==1) & (Z==1)][0]
    he       = Y[(A==4) & (Z==2)][0]
    try:
        d    = Y[(A==2) & (Z==1)][0]
    except:
        d    = 0

    he3      = Y[(A==3) & (Z==2)][0]
    li7      = Y[(A==7) & (Z==3)][0]
    try:
        be7  = Y[(A==7) & (Z==4)][0]
    except:
        be7  = 0
    hydrogren_list.append(hydrogen)
    he_list.append(he/hydrogen)
    d_list.append(d/hydrogen)
    he3_list.append(he3/hydrogen)
    li7_list.append((li7+be7)/hydrogen)


# Plot the Planck Satellite
plt.axvspan(5.96e-10,6.22e-10,alpha=0.5,color='k')
plt.text(5.5e-10,1e-6,'Planck Satellite',color='k',ha='center',va='top',rotation=90)

# Plotting properties
alpha=0.3

# Deuterium
# Cooke, R.~J., Pettini, M., \& Steidel, C.~C.\ 2018, \apj, 855, 102. doi:10.3847/1538-4357/aaab53
plt.axhspan( (2.527- 0.03) * 10**(-5), (2.527+ 0.03) * 10**(-5),color = 'tab:orange',alpha=alpha)
plt.text(1.58e-10,1e-3,'D',color='tab:orange',ha='center')

# Helium 3
# Bania, T.~M., Rood, R.~T., \& Balser, D.~S.\ 2002, \nat, 415, 54. doi:10.1038/415054a
plt.axhspan((1.1-0.2)*10**(-5),(1.1+0.2)*10**(-5),color = 'tab:green',alpha=alpha)
plt.text(6e-12,3e-4,'He 3',color='tab:green',ha='center')


# Helium 4
# Aver, E., Olive, K. A., and Skillman, E. D. A new approach to systematic uncertainties and self-consistency in helium abundance determinations. JCAP 2010, 05 (2010), 003.
hydrogren_list = np.array(hydrogren_list)
ind = np.argmin(abs(etas - 6.275e-11))
h = hydrogren_list[ind]
plt.axhspan((0.2561-0.0108)/4./h,(0.2561+0.0108)/4./h,color = 'tab:blue',alpha=alpha)
plt.text(1.6e-8,0.24,'He 4',color='tab:blue',ha='center')

# Lithium + Berillium
plt.axhspan((1.23-0.32)*1e-10,(1.23+0.68)*1e-10,color = 'tab:red',alpha=alpha)
plt.text(2.8e-11,5e-9,'Li 7 + Be 7',color='r',ha='center')

# Plot the calculated final abundances
plt.plot(etas,he_list,lw=2)
plt.plot(etas,d_list,lw=2)
plt.plot(etas,he3_list,lw=2)
plt.plot(etas,li7_list,lw=2)

# Set the scale and limits of the axes
plt.xscale('log')
plt.yscale('log')
plt.ylim(1e-11,1)
plt.xlim(1e-12,1e-7)

# Set the labels
plt.xlabel('$\eta$')
plt.ylabel('Abundance relative to hydrogen')

# Save and show the figure
plt.savefig('different_etas.pdf',bbox_inches='tight')
plt.show()
