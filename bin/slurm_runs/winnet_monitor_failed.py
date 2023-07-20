#!/usr/bin/env python
import numpy             as np
import matplotlib.pyplot as plt
import pandas            as pd
from tqdm                import tqdm
from datetime            import datetime,time,timedelta
from winnet_class        import winnet
import sys
import subprocess
import os
import h5py

# Define user
user = "mreicher"

use_slurm_times = False


if len(sys.argv) < 2:
    path ="."
else:
    path = sys.argv[1]

# Get number of slurm processes
cmd        = "squeue | grep "+user+" | grep R | wc -l"
x          = subprocess.check_output(cmd, shell=True)
nr_process = int(x.strip())

allTotal = 0
allRem   = 0
allFin   = 0
allFail  = 0

eltime = []
all_folders = os.listdir(path)
faillist  = []
ecodelist = []
for p in all_folders:
    tot_path = os.path.join(path,p)
    if not os.path.isdir(tot_path):
        continue

    w_path = os.path.join(tot_path,"winnet")
    b_path = os.path.join(tot_path,"blocked")
    f_path = os.path.join(tot_path,"finab.dat")
    h_path = os.path.join(tot_path,"WinNet_data.h5")
    o_path = os.path.join(tot_path,"OUT")
    e_path = os.path.join(tot_path,"ERR")
    # Some other folder
    if not os.path.isfile(w_path) and not os.path.isfile(b_path):
        continue
    elif not os.path.isfile(w_path) and os.path.isfile(b_path):
        # Either Fail or finished
        if os.path.isfile(f_path):
            # Read elapsed time of finshed runs
            if not use_slurm_times:
                w = winnet(tot_path)
                w.read_OUT()
                eltime.append(w.elapsed_time)
            allFin += 1
        elif os.path.isfile(h_path):
            ftmp = h5py.File(h_path,"r")
            if "finab/" in ftmp:
                allFin += 1
            else:
                allFail += 1
                faillist.append(p)
                # Try to read the error code
                try:
                    with open(e_path,"r") as f:
                        lines = f.readlines()
                        ecode = lines[1].split()[1].strip()
                except:
                    ecode="-"
                ecodelist.append(ecode)

            ftmp.close()
        else:
            allFail += 1
            faillist.append(p)
            # Try to read the error code
            try:
                with open(e_path,"r") as f:
                    lines = f.readlines()
                    ecode = lines[1].split()[1].strip()
            except:
                ecode="-"
            ecodelist.append(ecode)
    elif os.path.isfile(w_path) and os.path.isfile(o_path):
        allRem += 1

    allTotal +=1



outstr = ""
outstr += "          WinNet monitoring     "+"\n"
outstr += "===================================="+"\n"
# outstr += "\n"
outstr += "| List of failed runs:             | \n"
outstr += "|----------------------------------| \n"
outstr += "| Name                | Error code | \n"
lll = ""
for ind,p in enumerate(faillist):
    lll += "| "+p[:19].ljust(19)+" | "+ecodelist[ind].ljust(10)+" | \n"

outstr += lll
outstr += "|"+"_"*34 +"| \n"

outstr += 'Different error codes: \n'
outstr += ' '.join(list(set(ecodelist)))


# Output
print(outstr)
