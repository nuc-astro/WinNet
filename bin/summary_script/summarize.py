# Author: M. Reichert
# Date  : 25.09.2024
import numpy as np
import h5py
import sys
import os
import optparse
import logging
import re
from tqdm                             import tqdm
from src_files.wreader                import wreader
from src_files.template_class         import template
from src_files.nucleus_multiple_class import nucleus_multiple



#--- define options ----------------------------------------------------------
p = optparse.OptionParser()
p.add_option("-i","--input" , action="store", dest="rundir",  default='.',  \
  help="Simulation directory to summarize (default: current directory)")
p.add_option("-o","--output", action="store", dest="outdir",  default='./summary.hdf5',  \
  help="Output path (default: summary.hdf5)")
p.add_option("-b", "--buf"  , action="store", dest="buffersize",  default='500',  \
  help="Buffer size before writing to the file (default: 500)")
p.add_option("-f", "--force", action="store_true", dest="force",  default=False,  \
    help="Force overwrite of the output file (default: False)")
p.add_option("-v", "--verbose", action="store_true", dest="verbose",  default=False,  \
    help="Enable verbose output. If enabled, output is written to 'debug.log' (default: False)")
p.add_option("--time_file", action="store", dest="time_file",  default=None,  \
    help="File containing a time grid in seconds to map the data to (default: None)")
p.add_option("--time_final", action="store", dest="time_final",  default=None,  \
    help="Final time for the time grid. Only used if '--time_file' "+\
          "is not given (default: Read from template file)")
p.add_option("--time_initial", action="store", dest="time_initial",  default=None,  \
    help="Initial time for the time grid. Only used if '--time_file' "+\
          "is not given (default: 1e-5)")
p.add_option("--time_number", action="store", dest="time_number",  default=None,  \
    help="Number of time steps for the time grid. Only used if '--time_file' "+\
            "is not given (default: 200)")
p.add_option("--sunet_path", action="store", dest="sunet_path",  default=None,  \
    help="Path to the sunet file (default: Read from template)")
p.add_option("--disable_mainout", action="store_true", dest="disable_mainout",  default=False,  \
    help="Disable the summary of the mainout output (default: False)")
p.add_option("--disable_energy", action="store_true", dest="disable_energy",  default=False,  \
    help="Disable the summary of the energy output (default: False)")
p.add_option("--disable_timescales", action="store_true", dest="disable_timescales",  default=False,  \
    help="Disable the summary of the timescales output (default: False)")
p.add_option("--disable_tracked_nuclei", action="store_true", dest="disable_tracked_nuclei",  default=False,  \
    help="Disable the summary of the tracked_nuclei output (default: False)")
p.add_option("--disable_nuloss", action="store_true", dest="disable_nuloss",  default=False,  \
    help="Disable the summary of the nuloss output (default: False)")
p.add_option("--disable_snapshots", action="store_true", dest="disable_snapshots",  default=False,  \
    help="Disable the summary of the snapshots output (default: False)")
p.set_usage("""

  Usage:   ./summarize.py -i <rundir>
  Example: ./summarize.py -i runs/test""")


#--- parse options -----------------------------------------------------------
(options,args) = p.parse_args()
run_path = options.rundir


# Remove later and make it an parameter
buffsize = int(options.buffersize)


# Verbose mode or not?
if options.verbose:
    # Set up a logger to trace the progress
    logging.basicConfig(
        format="%(asctime)s - %(levelname)-10s - %(message)s",
        style="%",
        datefmt="%Y-%m-%d %H:%M",
        level=logging.DEBUG,
    )
    file_handler = logging.FileHandler("debug.log", mode="w", encoding="utf-8")
    logger = logging.getLogger(__name__)
    # Set the style also for the file handler
    file_handler.setFormatter(logging.Formatter(
        fmt="%(asctime)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M",
        style="%",
    ))
    logger.addHandler(file_handler)
else:
    # Set up a logger to trace the progress
    logging.basicConfig(level=logging.ERROR)
    logger = logging.getLogger(__name__)

# Disable the logging in the terminal
logger.propagate = False
# Say something
logger.info(f"Started summarizing run at {run_path}.")

# Get a list of all directories in the run_path. Ignore "network_data" directory
dirs = [d for d in os.listdir(run_path) if os.path.isdir(os.path.join(run_path, d)) and d != "network_data"]

# Say something
logger.info(f"Found {len(dirs)} directories in {run_path}.")


# Create output hdf5 file
output_file = options.outdir
# Check if the file already exists
if os.path.exists(output_file) and not options.force:
    logger.error(f"Output file {output_file} already exists. Exiting.")
    raise ValueError(f"Output file {output_file} already exists. Either delete or use -f option to overwrite.")
elif os.path.exists(output_file) and options.force:
    logger.warning(f"Output file {output_file} already exists. Overwriting it.")
f_hdf = h5py.File(output_file, 'w')


# Check already one run to see what has been outputted
# Find a run that didnt crash
found = False
for d in dirs:
    data = wreader(os.path.join(run_path, d))
    if not data.is_crashed:
        logger.info(f"Found {d} to look up the output.")
        found = True
        break
if not found:
    # Raise an error if all runs are crashed
    logger.error("All runs are crashed!")
    raise ValueError("All runs are crashed!")

# Get the template name (ends with .par)
template_name = [f for f in os.listdir(os.path.join(run_path, d)) if f.endswith('.par')][0]
t = template(os.path.join(run_path, d, template_name))

# Check if the sunet is given in the template
if options.sunet_path is not None:
    net_source = options.sunet_path
else:
    if (not "net_source" in t.entries):
        # Raise an error if the net_source is not given
        raise ValueError("net_source not given in the template file.")
    else:
        # Get the net_source
        net_source = t["net_source"]

# Read the net_source file
logger.info(f"Using sunet file from {net_source}.")
nuclei = np.loadtxt(net_source,dtype=str)
nuclei_data = nucleus_multiple(nuclei)


# Create the time grid
if options.time_file is not None:
    logger.info(f"Using time grid from {options.time_file}.")
    mainout_time = np.loadtxt(options.time_file, dtype=float, unpack=True)
else:
    if options.time_final is not None:
        final_time = float(options.time_final)
    else:
        final_time = 3.15e16
        # Make some reasonable time grid
        # Check if termination criterion is given
        if "termination_criterion" in t.entries:
            if t["termination_criterion"] == "1":
                if "final_time" in t.entries:
                    final_time = float(t["final_time"])

    if options.time_initial is not None:
        initial_time = float(options.time_initial)
    else:
        initial_time = 1e-5

    if options.time_number is not None:
        time_number = int(options.time_number)
    else:
        time_number = 200

    mainout_time = np.logspace(np.log10(initial_time), np.log10(final_time), time_number)


# Possible entries in the data
possible_entries = []
if not options.disable_mainout:
    possible_entries.append("mainout")
else:
    logger.info("Ignoring mainout output.")
if not options.disable_energy:
    possible_entries.append("energy")
else:
    logger.info("Ignoring energy output.")
if not options.disable_timescales:
    possible_entries.append("timescales")
else:
    logger.info("Ignoring timescales output.")
if not options.disable_tracked_nuclei:
    possible_entries.append("tracked_nuclei")
else:
    logger.info("Ignoring tracked_nuclei output.")
if not options.disable_nuloss:
    possible_entries.append("nuloss")
else:
    logger.info("Ignoring nuloss output.")


entry_dict = {}
for entry in possible_entries:
    # Check existence
    if data.check_existence(entry) != 0:
        entry_dict[entry] = {}
        for key in data[entry].keys():
            # Ignore iteration and time key
            if ((key == "iteration") or (key == "time") or (key == "A") or (key == "Z") or (key == "N")
               or (key == "names") or (key == "latex_names")):
                continue
            # Ignore temperature, density, and radius for nuloss
            if (entry == "nuloss") and ((key == "temp") or (key == "dens") or (key == "rad")):
                continue
            entry_dict[entry][key] = np.zeros((len(mainout_time),buffsize))

        # Write the time already
        f_hdf[entry+"/time"] = mainout_time

        # Say something
        logger.info(f"Found {entry} in the data.")

        # Check if it is the tracked nuclei and put A, Z, and N
        # if entry == "tracked_nuclei":
        #     f_hdf[entry+"/A"] = data[entry]["A"]
        #     f_hdf[entry+"/Z"] = data[entry]["Z"]
        #     f_hdf[entry+"/N"] = data[entry]["N"]
        #     f_hdf[entry+"/names"] = data[entry]["names"]


# Take care of snapshots, check if they are custom or not
if (data.check_existence("snapshot") != 0) and (not options.disable_snapshots):
    if ("custom_snapshots" in t.entries) or ("h_custom_snapshots" in t.entries):
        # Check if either ascii or hdf5 custom snapshots are given
        summarize_snapshots = False
        if ("custom_snapshots" in t.entries):
            summarize_snapshots = (t["custom_snapshots"].strip().lower() == "yes")
            # Debug statement
            if (t["custom_snapshots"].strip().lower() == "yes"):
                logger.info("Found custom snapshots in ascii format.")
        if ("h_custom_snapshots" in t.entries):
            summarize_snapshots = (summarize_snapshots or (t["h_custom_snapshots"].strip().lower() == "yes"))
            # Debug statement
            if (t["h_custom_snapshots"].strip().lower() == "yes"):
                logger.info("Found custom snapshots in hdf5 format.")

        # Read the time for the custom snapshots
        if summarize_snapshots:
            if "snapshot_file" not in t.entries:
                raise ValueError("Invalid template file. snapshot_file not given in the template file.")
            snapshot_time = np.loadtxt(os.path.join(t["snapshot_file"]),dtype=float)
            # Convert from days to seconds
            snapshot_time *= 24*3600
            # Write the time already
            f_hdf["snapshots/time"] = snapshot_time
            # Write the A and Z data to the hdf5 file
            f_hdf["snapshots/A"] = nuclei_data.A
            f_hdf["snapshots/Z"] = nuclei_data.Z
            f_hdf["snapshots/N"] = nuclei_data.N
            # Create an array to buffer the data
            snapshot_data = np.zeros((len(nuclei),len(snapshot_time),buffsize))
            logger.info(f"Summarize custom snapshots as well.")
    else:
        summarize_snapshots = False
else:
    summarize_snapshots = False
    if not options.disable_snapshots:
        logger.info("Ignoring snapshots output.")


# Finab should always be in
finab_data_Y = np.zeros((len(nuclei),buffsize))
finab_data_X = np.zeros((len(nuclei),buffsize))
# Write already the A and Z data to the hdf5 file, this is the same for
# all runs
f_hdf["finab/A"] = nuclei_data.A
f_hdf["finab/Z"] = nuclei_data.Z
f_hdf["finab/N"] = nuclei_data.N

# Create array that will contain the names of the runs
# Array of strings:
run_names = np.zeros(buffsize,dtype="S100")
run_ids   = np.zeros(buffsize,dtype=int)

# Loop over all directories
ind = -1
for counter, d in enumerate(tqdm(dirs)):
    # Load data
    data = wreader(os.path.join(run_path, d),silent=True)

    if data.is_crashed:
        logger.warning(f"Run {d} is crashed. Skipping it.")
        continue

    # Increase the index
    ind += 1

    #### Finab ####
    ###############
    # Put the data in the finab_data, Notice that it should be at the right A and Z position
    # A is contained in data.finab["A"] and Z in data.finab["Z"]. It should fit to the nuclei_data.A and nuclei_data.Z
    # All of them are 1D arrays
    # Getting indices where match occurs
    indices = [(np.where((nuclei_data.A.astype(int) == A) & (nuclei_data.Z.astype(int) == Z))[0][0])
               for A, Z in zip(data.finab["A"].astype(int), data.finab["Z"].astype(int))]

    finab_data_Y[indices,ind % buffsize] = data.finab["Y"][:]
    finab_data_X[indices,ind % buffsize] = data.finab["X"][:]

    # Check if the finab_data is full and write it to the hdf5 file
    if ind % buffsize == 0 and ind != 0:
        # Check if the dataset is already created and if not create it
        if "finab/Y" not in f_hdf:
            f_hdf.create_dataset("finab/Y", (len(nuclei),ind+1), maxshape=(len(nuclei),None))
            f_hdf.create_dataset("finab/X", (len(nuclei),ind+1), maxshape=(len(nuclei),None))
        # If necessary extend the dataset
        if ind > buffsize:
            f_hdf["finab/Y"].resize((len(nuclei),ind+1))
            f_hdf["finab/X"].resize((len(nuclei),ind+1))
        # Write the data to the hdf5 file
        f_hdf["finab/Y"][:,ind-buffsize:ind] = finab_data_Y
        f_hdf["finab/X"][:,ind-buffsize:ind] = finab_data_X


    #### Run name ####
    ##################
    # Save the run name
    run_names[ind % buffsize] = d

    # Check if the run_names is full and write it to the hdf5 file
    if ind % buffsize == 0 and ind != 0:
        # Check if the dataset is already created and if not create it
        if "run_names" not in f_hdf:
            f_hdf.create_dataset("run_names", (ind+1,), maxshape=(None,),dtype="S100")
        # If necessary extend the dataset
        if ind > buffsize:
            f_hdf["run_names"].resize((ind+1,))
        # Write the data to the hdf5 file
        f_hdf["run_names"][ind-buffsize:ind] = run_names


    # Get numbers out from the run name
    try:
        run_ids[ind % buffsize] = int(re.findall(r'\d+', d)[0])
    except:
        run_ids[ind % buffsize] = -1

    # Check if the run_ids is full and write it to the hdf5 file
    if ind % buffsize == 0 and ind != 0:
        # Check if the dataset is already created and if not create it
        if "run_ids" not in f_hdf:
            f_hdf.create_dataset("run_ids", (ind+1,), maxshape=(None,),dtype=int)
        # If necessary extend the dataset
        if ind > buffsize:
            f_hdf["run_ids"].resize((ind+1,))
        # Write the data to the hdf5 file
        f_hdf["run_ids"][ind-buffsize:ind] = run_ids


    #### Custom snapshots ####
    ##########################

    if summarize_snapshots:
        # Get the time of the snapshots
        snapstime = data.snapshot_time
        # Now get the indexes of the entries that agree with snapshot_time
        indexes = np.searchsorted(snapstime, snapshot_time)
        # Put the data in the snapshot_data
        indices_nuclei = [(np.where((nuclei_data.A.astype(int) == A) & (nuclei_data.Z.astype(int) == Z))[0][0])
                           for A, Z in zip(data.A.astype(int), data.Z.astype(int))]
        # Store it
        snapshot_data[indices_nuclei,:,ind % buffsize] = data.Y[indexes][:, :].T

        # Check if the snapshot_data is full and write it to the hdf5 file
        if ind % buffsize == 0 and ind != 0:
            # Check if the dataset is already created and if not create it
            if "snapshots/Y" not in f_hdf:
                f_hdf.create_dataset("snapshots/Y", (len(nuclei),len(snapshot_time),ind+1), maxshape=(len(nuclei),len(snapshot_time),None))
                f_hdf.create_dataset("snapshots/X", (len(nuclei),len(snapshot_time),ind+1), maxshape=(len(nuclei),len(snapshot_time),None))
            # If necessary extend the dataset
            if ind > buffsize:
                f_hdf["snapshots/Y"].resize((len(nuclei),len(snapshot_time),ind+1))
                f_hdf["snapshots/X"].resize((len(nuclei),len(snapshot_time),ind+1))
            # Write the data to the hdf5 file
            f_hdf["snapshots/Y"][:,:,ind-buffsize:ind] = snapshot_data
            f_hdf["snapshots/X"][:,:,ind-buffsize:ind] = snapshot_data*nuclei_data.A[:,np.newaxis,np.newaxis]



    #### Other entries ####
    #######################

    for entry in entry_dict.keys():
        # Put the data in the data_dict
        for key in entry_dict[entry].keys():
            entry_dict[entry][key][:,ind % buffsize] = np.interp(mainout_time,data[entry]["time"],data[entry][key],left=np.nan,right=np.nan)

        # Check if the data_dict is full and write it to the hdf5 file
        if ind % buffsize == 0 and ind != 0:
            for key in entry_dict[entry].keys():
                # Check if the dataset is already created and if not create it
                if entry+"/"+key not in f_hdf:
                    f_hdf.create_dataset(entry+"/"+key, (len(mainout_time),ind+1), maxshape=(len(mainout_time),None))
                # If necessary extend the dataset
                if ind > buffsize:
                    f_hdf[entry+"/"+key].resize((len(mainout_time),ind+1))
                # Write the data to the hdf5 file
                f_hdf[entry+"/"+key][:,ind-buffsize:ind] = entry_dict[entry][key]



# Say something
logger.info(f"Finished looping over all directories, writing the last data to the hdf5 file.")


# Write the last data to the hdf5 file

#### Finab ####
###############

if "finab/Y" not in f_hdf:
    f_hdf.create_dataset("finab/Y", (len(nuclei),ind+1), maxshape=(len(nuclei),None))
    f_hdf.create_dataset("finab/X", (len(nuclei),ind+1), maxshape=(len(nuclei),None))
else:
    f_hdf["finab/Y"].resize((len(nuclei),ind+1))
    f_hdf["finab/X"].resize((len(nuclei),ind+1))
# Write the missing entries
if ind>buffsize:
    f_hdf["finab/Y"][:,ind-buffsize+1:ind+1] = finab_data_Y[:,:buffsize]
    f_hdf["finab/X"][:,ind-buffsize+1:ind+1] = finab_data_X[:,:buffsize]
else:
    f_hdf["finab/Y"][:,:ind+1] = finab_data_Y[:,:ind+1]
    f_hdf["finab/X"][:,:ind+1] = finab_data_X[:,:ind+1]


#### Run name ####
##################

if "run_names" not in f_hdf:
    f_hdf.create_dataset("run_names", (ind+1,), maxshape=(None,),dtype="S100")
else:
    f_hdf["run_names"].resize((ind+1,))
# Write the missing entries
if ind>buffsize:
    f_hdf["run_names"][ind-buffsize+1:ind+1] = run_names[:buffsize]
else:
    f_hdf["run_names"][:ind+1] = run_names[:ind+1]


if "run_ids" not in f_hdf:
    f_hdf.create_dataset("run_ids", (ind+1,), maxshape=(None,),dtype=int)
else:
    f_hdf["run_ids"].resize((ind+1,))
# Write the missing entries
if ind>buffsize:
    f_hdf["run_ids"][ind-buffsize+1:ind+1] = run_ids[:buffsize]
else:
    f_hdf["run_ids"][:ind+1] = run_ids[:ind+1]


#### Custom snapshots ####
##########################

if summarize_snapshots:
    if "snapshots/Y" not in f_hdf:
        f_hdf.create_dataset("snapshots/Y", (len(nuclei),len(snapshot_time),ind+1), maxshape=(len(nuclei),len(snapshot_time),None))
        f_hdf.create_dataset("snapshots/X", (len(nuclei),len(snapshot_time),ind+1), maxshape=(len(nuclei),len(snapshot_time),None))
    else:
        f_hdf["snapshots/Y"].resize((len(nuclei),len(snapshot_time),ind+1))
        f_hdf["snapshots/X"].resize((len(nuclei),len(snapshot_time),ind+1))
    # Write the missing entries
    if ind>buffsize:
        f_hdf["snapshots/Y"][:,:,ind-buffsize+1:ind+1] = snapshot_data[:,:,:buffsize]
        f_hdf["snapshots/X"][:,:,ind-buffsize+1:ind+1] = snapshot_data[:,:,:buffsize]*nuclei_data.A[:,np.newaxis,np.newaxis]
    else:
        f_hdf["snapshots/Y"][:,:,:ind+1] = snapshot_data[:,:,:ind+1]
        f_hdf["snapshots/X"][:,:,:ind+1] = snapshot_data[:,:,:ind+1]*nuclei_data.A[:,np.newaxis,np.newaxis]


#### Other entries ####
#######################

for entry in entry_dict.keys():
    for key in entry_dict[entry].keys():
        if entry+"/"+key not in f_hdf:
            f_hdf.create_dataset(entry+"/"+key, (len(mainout_time),ind+1), maxshape=(len(mainout_time),None))
        else:
            f_hdf[entry+"/"+key].resize((len(mainout_time),ind+1))
        # Write the missing entries
        if ind>buffsize:
            f_hdf[entry+"/"+key][:,ind-buffsize+1:ind+1] = entry_dict[entry][key][:,:buffsize]
        else:
            f_hdf[entry+"/"+key][:,:ind+1] = entry_dict[entry][key][:,:ind+1]



# Say something
logger.info(f"Finished summarizing run at {run_path}.")

# Close the hdf5 file
f_hdf.close()

