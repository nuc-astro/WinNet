# Author: M. Reichert
# Date  : 03.06.2024
# Convert beta delayed fission file of Mumpower et al. 2022 into a WinNet readable format.
# The Mumpower et al. 2022 file can be accessed at :
# https://journals.aps.org/prc/supplemental/10.1103/PhysRevC.106.065805
import numpy as np
import pandas as pd
from class_files.nucleus_class import nucleus


path = "mumpower_beoh350_bdf_frdm2012_bstr12-gtff_frldm2015.dat"


# Columns of the file
columns = ["name", "Z", "N",  "A", "P0", "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9", "P10",
                   "P0f", "P1f", "P2f", "P3f", "P4f", "P5f", "P6f", "P7f", "P8f", "P9f", "P10f"]

total_data =[]
rate_counter = 0

# Open and parse through the file
with open(path,"r") as f:
    lines = f.readlines()
    for l in lines:
        if (l[0] == "#") or (l.strip()==""):
            continue
        else:
            data = np.array(l.split()).astype(float)
            if len(data) == 16:
                # Get the nucleus name
                nucleus_name = nucleus(Z=int(data[0]), N=int(data[1])).get_name()
                # Create the first part of the entry
                entry = [nucleus_name]+list(data[0:14])
            else:
                # Create the second part of the entry
                entry = entry+list(data)
                # Add it to the total data
                total_data.append(entry)
                rate_counter += 1

# Create a pandas dataframe
df = pd.DataFrame(total_data, columns=columns)

# Calculate the total fission probability
df["Total_pf"] = df["P0f"]+df["P1f"]+df["P2f"]+df["P3f"]+df["P4f"]+df["P5f"]+df["P6f"]+df["P7f"]+df["P8f"]+df["P9f"]+df["P10f"]
# Drop entries that have no fission probability
df = df[df["Total_pf"]>0]


# Create the file
fcont = ""

# Go through the dataframe
for it, row in df.iterrows():
    fcont = fcont + row["name"].rjust(5) + 4*" "+4*" "+"mp22\n"
    for i in range(10):
        fcont = fcont + f"{row[f'P{i}f']:.3f}"+4*" "
    fcont = fcont + "\n"

# Write to file
with open("fission_rates_bdf_mp22.dat","w") as f:
    f.write(fcont)

# Say something
print(f"File fission_rates_bdf_mp22.dat created with {rate_counter} entries.")