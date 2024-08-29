import numpy             as np
import matplotlib.pyplot as plt
from   FlowAnimation     import FlowAnimation
from   mpi4py            import MPI
from   tqdm              import tqdm
import os
import sys
import pickle


# Load the dictionary from the pkl file
# Get folder location of this script
script_path = os.path.dirname(os.path.realpath(__file__))
option_dict_path = os.path.join(script_path, 'data/options.pkl')
with open(option_dict_path, 'rb') as f:
    kwargs = pickle.load(f)

# Also run ffmpeg per default unless told otherwise
create_movie = True

# Get the arguments from the command line
# Run path, Frame_min, Frame_max, interval
run_path  = sys.argv[1]
frame_min = int(sys.argv[2])
frame_max = int(sys.argv[3])
interval  = int(sys.argv[4])
if len(sys.argv) > 5:
    create_movie = False

# Set up the MPI environment
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# Set up the figure
fig = plt.figure(figsize=(15, 8))

# Set up the frame directory
frame_dir=f'{run_path}/frames'
if not os.path.exists(frame_dir):
    os.mkdir(frame_dir)

# Make barrier here
comm.Barrier()

# Set up the animation
anim = FlowAnimation(run_path, fig, **kwargs)

frames = np.arange(frame_min,frame_max)[rank::size]
if rank == 0:
    for ii in tqdm(frames):
        anim.save_frame(ii)
else:
    for ii in frames:
        anim.save_frame(ii)

# Wait again
comm.Barrier()

# Combine the frames
if ((rank == 0) and create_movie):
    os.system(f'ffmpeg -y -r {interval} -start_number {frame_min} -i {frame_dir}/frame_%d.png  {run_path}/flow.mp4')