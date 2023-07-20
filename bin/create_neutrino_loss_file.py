# Author: M. Reichert
# Date  : 31.03.2023
# Script to download the neutrino loss data from the ENDSF database
# and to put it into a WinNet readable file format.
# The file can be used to give more precise information on how
# much energy is radiated away by neutrinos in case that
# nuclear heating is enabled.
# A description about the data API can be found here:
# https://www-nds.iaea.org/relnsd/vcharthtml/api_v0_guide.html
import numpy as np
from tqdm import tqdm
import pandas as pd
import os
import h5py

# Name of the output filename
output_file = "nu_loss_data.dat"
# If the neutrino loss is zero, do not include it in the file
ignore_zeros = True
# The nuclei names are read from the sunet_complete file
sunet_path ="sunet_complete"
nuclei_names = np.loadtxt(sunet_path,usecols=[0],unpack=True,dtype=str)


# Create folders for the data
if not os.path.exists('bp'):
   os.makedirs('bp')
if not os.path.exists('bm'):
   os.makedirs('bm')

# Loop through the data and download it
for ind,n in enumerate(tqdm(nuclei_names)):
    baselink = '"https://www-nds.iaea.org/relnsd/v0/data?fields=decay_rads&nuclides='
    nuc = n
    user_agent = ' -U "Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:77.0) Gecko/20100101 Firefox/77.0" '
    bp = baselink+n+'&rad_types=bp"'
    bm = baselink+n+'&rad_types=bm"'
    # In principle, one can also download the other
    # decay information
    # g  = baselink+n+'&rad_types=g"'
    # auger= baselink+n+'&rad_types=e"'
    # alpha= baselink+n+'&rad_types=a"'
    # xrays= baselink+n+'&rad_types=x"'

    if not os.path.exists("bp/"+nuc+".csv"):
        cmd = "wget -O bp/"+nuc+".csv "+user_agent+bp
        os.system(cmd)

    if not os.path.exists("bm/"+nuc+".csv"):
        cmd = "wget -O bm/"+nuc+".csv "+user_agent+bm
        os.system(cmd)

# Get the neutrino loss
nu_loss = {}

# After download, loop over them and extract the data
for ind,n in enumerate(tqdm(nuclei_names)):
    try:
        df_bm = pd.read_csv("bm/"+n+".csv",sep=',',na_values=" ")
    except:
        df_bm = pd.DataFrame()
    try:
        df_bp = pd.read_csv("bp/"+n+".csv",sep=',',na_values=" ")
    except:
        df_bp = pd.DataFrame()


    if not df_bm.empty:
        df_bm['intensity_beta']      = pd.to_numeric(df_bm['intensity_beta'], errors='coerce')
        df_bm['anti_nu_mean_energy'] = pd.to_numeric(df_bm['anti_nu_mean_energy'], errors='coerce')
        av_anu = np.nansum(df_bm["anti_nu_mean_energy"].values*df_bm["intensity_beta"].values/100.)
        nu_loss[n] = av_anu/1000
    if not df_bp.empty:
        df_bp['intensity_beta'] = pd.to_numeric(df_bp['intensity_beta'], errors='coerce')
        df_bp['nu_mean_energy'] = pd.to_numeric(df_bp['nu_mean_energy'], errors='coerce')
        df_bp['energy_ec']      = pd.to_numeric(df_bp['energy_ec'], errors='coerce')
        df_bp['intensity_ec']   = pd.to_numeric(df_bp['intensity_ec'], errors='coerce')
        av_ec = np.nansum(df_bp["energy_ec"].values*df_bp["intensity_ec"].values/100.)

        av_nu = np.nansum(df_bp["nu_mean_energy"].values*df_bp["intensity_beta"].values/100.)

        if n in nu_loss:
            # The energy of the electron capture is added to the neutrino loss
            nu_loss[n] += (av_nu+av_ec)/1000
        else:
            # The energy of the electron capture is added to the neutrino loss
            nu_loss[n] = (av_nu+av_ec)/1000


# Now create the file
out = ""
for n in tqdm(nuclei_names):
    if n in nu_loss:
        if ignore_zeros and nu_loss[n] == 0:
            continue
        out += n.rjust(5)+" "+"{:16.6E}".format(nu_loss[n])+"\n"

# Save the file
with open(output_file,"w") as f:
    f.write(out)