#!/usr/bin/env python
# Authors: M. Jacobi, J. Kuske, M. Reichert
# Movie script inspired by Skynet (J. Lippuner)
import sys
import os
# Get the path of the script and add it to the path (necessary for the imports)
script_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(script_path,'src_files/'))
from src_files.FlowAnimation import FlowAnimation
from src_files.wreader       import wreader
import matplotlib.pyplot     as plt
import optparse



#--- define options ----------------------------------------------------------
p = optparse.OptionParser()
p.add_option("-i","--input", action="store", dest="rundir",  default='.',  \
  help="Simulation directory to visualize (default: current directory)")
p.add_option("--disable_flow", action="store_true", dest="plot_flow", default=False, \
  help="Whether or not to plot the flow arrows.")
p.add_option("--flow_min", action="store", dest="flow_min", default="", \
  help="Lower limit of the flow.")
p.add_option("--flow_max", action="store", dest="flow_max", default="", \
  help="Upper limit of the flow.")
p.add_option("--fix_flows", action="store_true", dest="fix_flows", default=False, \
  help="Whether or not the flows are adapted to the data or lie between flow_min and flow_max.")
p.add_option("--flow_range", action="store", dest="flow_range", default="", \
  help="Log range of the flows in case that they are not fixed.")
p.add_option("--fix_flow_arrow_width", action="store_true", dest="fix_flow_arrow_width", default=False, \
  help="Fix the width of the flow arrows to constant width.")
p.add_option("--flow_cmap", action="store", dest="flow_cmap", default="", \
  help="Colormap of the flows.")
p.add_option("--separate_fission", action="store_true", dest="separate_fission", default=False, \
  help="Whether or not to show arrows also for fission. If not present, hatched areas will be plotted.")
p.add_option("--fission_minflow", action="store", dest="fission_minflow", default="", \
  help="Minimum flow to get indicated as fission region in case the separate fission flag is not given.")
p.add_option("--x_min", action="store", dest="x_min", default="", \
  help="Lower limit of the mass fraction.")
p.add_option("--x_max", action="store", dest="x_max", default="", \
  help="Upper limit of the mass fraction.")
p.add_option("--x_cmap", action="store", dest="x_cmap", default="", \
  help="Colormap of the mass fractions.")
p.add_option("--disable_abar", action="store_true", dest="disable_abar", default=False, \
  help="Whether or not disabling the indication of Abar.")
p.add_option("--mass_bins_cmap", action="store", dest="mass_bins_cmap", default="", \
  help="Colormap of the background colors.")
p.add_option("--disable_magic", action="store_true", dest="disable_magic", default=False, \
  help="Whether or not disabling the indication for the magic number.")
p.add_option("--additional_plot", action="store", dest="additional_plot", default="", \
  help="Whether or not plotting average timescales.")
p.add_option("--tau_min", action="store", dest="tau_min", default="", \
  help="Lower limit of the average timescales.")
p.add_option("--tau_max", action="store", dest="tau_max", default="", \
  help="Upper limit of the average timescales.")
p.add_option("--engen_min", action="store", dest="engen_min", default="", \
  help="Lower limit of the Energy.")
p.add_option("--engen_max", action="store", dest="engen_max", default="", \
  help="Upper limit of the Energy.")
p.add_option("--tracked_min", action="store", dest="tracked_min", default="", \
  help="Lower limit of the tracked nuclei mass fractions.")
p.add_option("--tracked_max", action="store", dest="tracked_max", default="", \
    help="Upper limit of the tracked nuclei mass fractions.")
p.add_option("--time_min", action="store", dest="t_min", default="", \
  help="Lower limit of the time.")
p.add_option("--time_max", action="store", dest="t_max", default="", \
  help="Upper limit of the time.")
p.add_option("--disable_mainout", action="store_true", dest="disable_mainout", default=False, \
  help="Whether or not disabling the mainout plot.")
p.add_option("--density_min", action="store", dest="density_min", default="", \
  help="Lower limit of the density.")
p.add_option("--density_max", action="store", dest="density_max", default="", \
  help="Upper limit of the density.")
p.add_option("--temperature_min", action="store", dest="temperature_min", default="", \
  help="Lower limit of the temperatures.")
p.add_option("--temperature_max", action="store", dest="temperature_max", default="", \
  help="Upper limit of the temperature.")
p.add_option("--ye_min", action="store", dest="ye_min", default="", \
  help="Lower limit of the electron fraction.")
p.add_option("--ye_max", action="store", dest="ye_max", default="", \
  help="Upper limit of the electron fraction.")
p.add_option("--frame_min", action="store", dest="frame_min", default="", \
  help="Value of the first frame (default = 1).")
p.add_option("--frame_max", action="store", dest="frame_max", default="", \
  help="Value of the last frame (default = end of the simulation).")
p.add_option("--save", action="store_true", dest="save", default=False, \
  help="Whether or not saving the movie.")
p.add_option("--output", action="store", dest="output_name", default='flow_movie.mp4', \
  help="Output name of the movie.")
p.add_option('--parallel_save', action='store_true', dest='parallel_save', default=False, \
             help='Whether or not to save the movie in parallel.')
p.add_option('--parallel_cpus', action='store', dest='parallel_cpus', default='5', \
             help='Number of CPUs to use for parallel saving.')
p.add_option("--interval", action="store", dest="interval", default='10', \
  help="Interval of the movie (larger value equals slower).")
p.add_option("--mpirun_path", action="store", dest="mpirun_path", default='', \
  help="Path of the mpirun command to use for parallel saving.")
p.set_usage("""
  Visualize a WinNet simulation. Ensure that at least
  snapshot_every or h_snapshot_every parameter was enabled in the
  parameter file. To plot timescales, energy, tracked nuclei, mainout,
  or reaction flows, the necessary parameters have to be enabled in the
  parameter file.

  Usage:   ./winnet_movie.py -i <rundir>
  Example: ./winnet_movie.py -i runs/test""")


#--- parse options -----------------------------------------------------------
(options,args) = p.parse_args()
run_path = options.rundir

kwargs = {}
kwargs['timescalerange']   = (1e-12, 1e10)
kwargs['trackedrange']     = (1e-8, 1e0)
kwargs['energyrange']      = (1e10, 1e20)
kwargs['timerange']        = (1e-5 , 1e5)
kwargs['densityrange']     = (1e-5, 1e12)
kwargs['temperaturerange'] = (0, 10)
kwargs['yerange']          = (0.0, 0.55)

if options.flow_min:  kwargs['flow_min'] = float(options.flow_min)
if options.flow_max:  kwargs['flow_max'] = float(options.flow_max)
if options.plot_flow: kwargs['plot_flow'] = False if options.plot_flow else True
if options.separate_fission: kwargs['separate_fission'] = True
if options.fission_minflow: kwargs['fission_minflow'] = float(options.fission_minflow)
if options.fix_flows: kwargs['flow_adapt_prange'] = False
if options.flow_range: kwargs['flow_prange'] = float(options.flow_range)
if options.fix_flow_arrow_width: kwargs['flow_adapt_width'] = True
if options.flow_cmap: kwargs['cmapNameFlow'] = options.flow_cmap
if options.x_min: kwargs['X_min'] = float(options.x_min)
if options.x_max: kwargs['X_max'] = float(options.x_max)
if options.x_cmap: kwargs['cmapNameX'] = options.x_cmap
if options.disable_abar: kwargs['plot_abar'] = (not options.disable_abar)
if options.mass_bins_cmap: kwargs['cmapNameMassBins'] = options.mass_bins_cmap
if options.disable_magic: kwargs['plot_magic'] = (not options.disable_magic)
if options.additional_plot: kwargs['additional_plot'] = options.additional_plot.lower().strip()
if options.tau_min: kwargs['timescalerange'] = (float(options.tau_min),kwargs['timescalerange'][1])
if options.tau_max: kwargs['timescalerange'] = (kwargs['timescalerange'][0], float(options.tau_max))
if options.engen_min: kwargs['energyrange'] = (float(options.engen_min), kwargs['energyrange'][1])
if options.engen_max: kwargs['energyrange'] = (kwargs['energyrange'][0], float(options.engen_max))
if options.tracked_min: kwargs['trackedrange'] = (float(options.tracked_min), kwargs['trackedrange'][1])
if options.tracked_max: kwargs['trackedrange'] = (kwargs['trackedrange'][0], float(options.tracked_max))
if options.t_min: kwargs['timerange'] = (float(options.t_min), kwargs['timerange'][1])
if options.t_max: kwargs['timerange'] = (kwargs['timerange'][0], float(options.t_max))
if options.disable_mainout: kwargs['plot_mainout'] = (not options.disable_mainout)
if options.density_min: kwargs['densityrange'] = (float(options.density_min), kwargs['densityrange'][1])
if options.density_max: kwargs['densityrange'] = (kwargs['densityrange'][0], float(options.density_max))
if options.temperature_min: kwargs['temperaturerange'] = (float(options.temperature_min), kwargs['temperaturerange'][1])
if options.temperature_max: kwargs['temperaturerange'] = (kwargs['temperaturerange'][0], float(options.temperature_max))
if options.ye_min: kwargs['yerange'] = (float(options.ye_min), kwargs['yerange'][1])
if options.ye_max: kwargs['yerange'] = (kwargs['yerange'][0], float(options.ye_max))



# Sanity checks
w = wreader(run_path)

# Check if the run has snapshots
value = w.check_existence('snapshot')
if value == 0:
    raise ValueError('No snapshots found. Please enable snapshots in the parameter file.')

# Sanity for timescales and so on, disable if not found
if options.additional_plot:
    if options.additional_plot == 'timescales':
        value = w.check_existence('timescales')
        if value == 0:
            print('No timescales found. Disabling timescales. Remove --additional_plot to disable this message.')
            kwargs['additional_plot'] = 'none'
    elif options.additional_plot == 'energy':
        value = w.check_existence('energy')
        if value == 0:
            print('No energy found. Disabling energy. Remove --additional_plot to disable this message.')
            kwargs['additional_plot'] = 'none'
    elif options.additional_plot == 'tracked':
        value = w.check_existence('tracked_nuclei')
        if value == 0:
            print('No tracked nuclei found. Disabling tracked nuclei. Remove --additional_plot to disable this message.')
            kwargs['additional_plot'] = 'none'
if not options.disable_mainout:
    value = w.check_existence('mainout')
    if value == 0:
        print('No mainout found. Disabling mainout. Set --disable_mainout to disable this message.')
        kwargs['plot_mainout'] = False
if not options.plot_flow:
    value = w.check_existence('flows')
    if value == 0:
        print('No flow found. Disabling flow. Set --disable_flow to disable this message.')
        kwargs['plot_flow'] = False


if options.frame_min: frame_min = int(options.frame_min)
else: frame_min = 1
if options.frame_max: frame_max = int(options.frame_max)
else: frame_max = w.nr_of_snaps

# Check if things should be saved or shown
if options.save:
    if not options.parallel_save:
        # Funcanimation
        fig = plt.figure(figsize=(15, 8))

        # Animate the flows
        anim = FlowAnimation(run_path, fig, **kwargs)

        ani = anim.get_funcanimation(frames=range(frame_min, frame_max))
        ani.save(options.output_name, fps=int(options.interval))
    else: # Parallel saving
        # Sanity check, does mpi4py exist?
        try:
            from mpi4py import MPI
        except ImportError:
            raise ImportError('mpi4py not found. Please install it to use parallel saving.')

        # Next check, is ffmpeg installed?
        if os.system('ffmpeg -version > /dev/null') != 0:
            raise ImportError('ffmpeg not found. Please install it to use parallel saving.')

        # Get folder location of this script
        script_path = os.path.dirname(os.path.realpath(__file__))
        # The options have to be passed to the parallel_save.py script
        # Therefore save them
        # try to import pickle
        try:
            import pickle
        except ImportError:
            raise ImportError('pickle not found. Please install it to use parallel saving.')

        option_dict_path = os.path.join(script_path, 'src_files/data/options.pkl')
        with open(option_dict_path, 'wb') as f:
            pickle.dump(kwargs, f)

        # Get the path to the parallel_save.py script
        parallel_save_path = os.path.join(script_path, 'src_files', 'parallel_save.py')

        # Check if the mpirun path is given
        if not options.mpirun_path:
            # Get the correct mpirun
            mpirun = os.popen('whereis mpirun').read().strip().split()
            # Get the mpirun that has oneapi in the path
            mpirun = [m for m in mpirun if 'oneapi' in m]
            if not mpirun:
                raise ImportError('No mpirun with oneapi found. Please install it to use parallel saving.')
            mpirun = mpirun[0]
        else:
            # Take the given mpirun path
            mpirun = options.mpirun_path

        # Test if the mpirun path is correct
        if os.system(f'{mpirun} -version > /dev/null') != 0:
            raise ImportError('mpirun not found or wrong path. Please install it to use parallel saving.')

        # Run the parallel saving
        os.system(f'{mpirun} -n {options.parallel_cpus} python {parallel_save_path} {run_path} {frame_min} {frame_max} {options.interval}')

        # Remove the options file
        os.remove(option_dict_path)

    print('Finished saving movie!')
else:
    # Funcanimation
    fig = plt.figure(figsize=(15, 8))

    # Animate the flows
    anim = FlowAnimation(run_path, fig, **kwargs)
    ani = anim.get_funcanimation(interval=int(options.interval), frames=range(frame_min, frame_max))
    plt.show()
