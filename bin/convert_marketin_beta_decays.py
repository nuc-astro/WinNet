# Author: M. Reichert
# Date  : 30.03.2023
# Convert beta decay file of Marketin into a WinNet readable format.
# The Marketing file can be accessed at :
# https://journals.aps.org/prc/supplemental/10.1103/PhysRevC.93.025805
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Names of the elements
elname = ['neutron','h','he','li','be','b','c','n','o','f','ne','na','mg','al','si','p','s','cl','ar','k','ca','sc','ti','v','cr','mn','fe',
          'co','ni','cu','zn','ga','ge','as','se','br','kr','rb','sr','y','zr','nb','mo','tc','ru','rh','pd','ag','cd','in','sn','sb',
          'te', 'i','xe','cs','ba','la','ce','pr','nd','pm','sm','eu','gd','tb','dy','ho','er','tm','yb','lu','hf','ta','w','re','os',
          'ir','pt','au','hg','tl','pb','bi','po','at','rn','fr','ra','ac','th','pa','u','np','pu','am','cm','bk','cf','es','fm','md',
          'no','lr','rf','db','sg','bh','hs','mt','ds','rg','ub','ut','uq','up','uh','us','uo']



# The file with the marketin rates
path = "betadecay_supplemental_material.dat"
# Get the column names
with open(path,"r") as f:
    lines = f.readlines()

    for line in lines:
        if len(line.split("  "))>2:
            linesplit = line.split("  ")[1:]
            linesplit = list(map(lambda x: x.strip(),linesplit))

            if "Z" in linesplit:
                linesplit = np.array(["name"]+linesplit)
                header = linesplit[linesplit!=""]
                break

df = pd.read_csv(path,names=header, comment='#',delim_whitespace=True)

out = ""
for ind, row in df.iterrows():
    ztmp = int(row["Z"])
    atmp = int(row["A"])
    name = elname[ztmp]+str(atmp)

    halflive = row["T_1/2 [s]"]
    Etot     = row[["<E_e> [MeV]","<E_nu> [MeV]", "<E_gamma> [MeV]"]].sum()
    Enu      = row["<E_nu> [MeV]"]
    if halflive == 1e100:
        print("Skipped "+str(name)+" due to large halflife.")
        continue
    print(Etot,Enu,Enu/Etot)
    out+=name.rjust(5)+" "*4+"{:12.6e}".format(halflive)+"  "+"{:10.6e}".format(Etot)+"  "+"{:10.6e}".format(Enu)+"\n"
    pns = np.zeros(11)
    for i in range(6):
        pns[i] = row["P"+str(i)+"n"]
    for p in pns:
        out+="{:1.4f}".format(p)+"   "
    out+="\n"


with open("marketin.dat","w") as f:
    f.write(out)


