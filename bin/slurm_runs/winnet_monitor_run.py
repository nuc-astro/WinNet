#!/usr/bin/env python
# Author: M. Reichert
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

# Define slurm user here (not really necessary, but it may improve the estimate)
user = "User"

# Use slurm to estimate times?
# This will speed up the monitoring
use_slurm_times = False

# Get the input
if len(sys.argv) < 2:
    path ="."
else:
    path = sys.argv[1]

# Get number of slurm processes
cmd        = "squeue | grep "+user+" | grep R | wc -l"
x          = subprocess.check_output(cmd, shell=True)
nr_process = int(x.strip())

# Set variables initially to zero
allTotal = 0
allRem   = 0
allFin   = 0
allFail  = 0

eltime = []
all_folders = os.listdir(path)
for p in all_folders:
    tot_path = os.path.join(path,p)
    if not os.path.isdir(tot_path):
        continue

    w_path = os.path.join(tot_path,"winnet")
    b_path = os.path.join(tot_path,"blocked")
    f_path = os.path.join(tot_path,"finab.dat")
    h_path = os.path.join(tot_path,"WinNet_data.h5")
    o_path = os.path.join(tot_path,"OUT")
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
                if not use_slurm_times:
                    w = winnet(tot_path)
                    w.read_OUT()
                    eltime.append(w.elapsed_time)
                allFin += 1
            else:
                allFail += 1
        else:
            allFail += 1
    elif os.path.isfile(w_path) and os.path.isfile(b_path):
        allRem += 1

    allTotal +=1

if use_slurm_times:
    # Get the time running already
    cmd        = "squeue | grep "+user+" | grep R | awk '{print $6}' "
    x          = subprocess.check_output(cmd, shell=True)

    times = x.split("\n")
    ttt   = [] # time in hours
    for t in times:
        if t.strip() == "":
            continue
        days = t.split("-")
        if len(days)==2:
            tottime = float(days[0])*24
            days    = days[1]
        else:
            tottime = 0
            days    = days[0]
        days = days.split(":")

        if len(days)== 2:
            tottime += float(days[0])/60.+float(days[1])/60./60.
        elif len(days)== 3:
            tottime += float(days[0])+float(days[1])/60.+float(days[2])/60./60.

        ttt.append(tottime)
    # Average time [h] the processes were running
    av_time = np.average(ttt)
else:
    ttt = np.array(eltime)/60./60.

    if allRem!=0:
        av_per_h = 1./(np.nanmean(eltime)/60./60./allRem)
        av_time  = allFin/av_per_h
    elif nr_process!=0:
        av_per_h = 1./allFin/(np.nanmean(eltime)/60./60./nr_process)
        av_time  = allFin/av_per_h
    else:
        av_per_h = 1./allFin/(np.nanmean(eltime)/60./60./100.)
        av_time  = allFin/av_per_h

# Run per hour
if allFin != 0 and len(ttt)>0:
    finished = False
elif allFin+allFail == allTotal:
    finished =True
else:
    if allRem != 0:
        # Make rough estimate, 15min per run
        av_per_h = 1./((1.5/6.)/allRem)
    elif nr_process!=0:
        av_per_h = 1./((1.5/6.)/nr_process)
    else:
        av_per_h = 1./((1.5/6.)/100)
    finished = False


# Estimated rest time
rest_runs = allTotal-allFin-allFail
if not finished:
    est_time = rest_runs/av_per_h
    # Make it well readable
    est_time = np.round(est_time*60*60)
    secs = int(est_time % 60)
    est_time = (est_time-secs)/60.
    mins = int(est_time % 60)
    est_time = int((est_time-mins)/60.)
    hours = est_time
else:
    secs  = 0
    mins  = 0
    hours = 0

timestr = str(hours).zfill(2)+":"+str(mins).zfill(2)+":"+str(secs).zfill(2)

# Duration
date_object = timedelta(hours=hours,minutes=mins,seconds=secs)
# Current date
now = datetime.now()

# Estimated time of finish
finishtime   = now+date_object

# Average time per tracer
est_time = np.round(av_time*60*60)
secs2 = int(est_time % 60)
est_time = (est_time-secs)/60.
mins2 = int(est_time % 60)
est_time = int((est_time-mins)/60.)
hours2 = est_time
timestr2 = str(hours2).zfill(2)+":"+str(mins2).zfill(2)+":"+str(secs2).zfill(2)


# Create the output

outstr = ""
outstr += "          WinNet monitoring     "+"\n"
outstr += "===================================="+"\n"
# outstr += "\n"
outstr += "| Number of runs      | "+str(allTotal).rjust(10)+" | \n"
outstr += "| Running             | "+str(allRem).rjust(10)  +" | \n"
outstr += "| Finished            | "+str(allFin).rjust(10)  +" | \n"
outstr += "| Failed              | "+str(allFail).rjust(10) +" | \n"
outstr += "| --------------------------------"+" | \n"
# outstr += "\n"
outstr += "| Av. time / tracer   | "+timestr2.rjust(10) +" | \n"
outstr += "| Estimated duration  | "+timestr.rjust(10)  +" | \n"
outstr += "| Estimated finish    | "+finishtime.strftime("%H:%M:%S").rjust(10) +" | \n"
outstr += "|"+"_"*34 +"| \n"


# Output
print(outstr)
