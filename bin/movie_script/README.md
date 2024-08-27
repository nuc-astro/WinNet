## Folder with movie script

This folder contains a script to analyze a WinNet run. It will show or save a movie of the mass fractions over time in the nuclear chart, inspired by the movie of SkyNet (done by J. Lippuner).

To create a video, the run must have snapshot output enabled (parameters snapshots_every or h_snapshots_every). Also timescales, mainout
properties, or abundance flows can be plotted if the frequency of the output is **set to the same value as the one of the snapshots**. 

An example command to generate the video is:

```bash
    python winnet_movie.py -i ../../runs/Example_NSM_dyn_ejecta_rosswog
```    

This command will display the movie in a separate window. Note that this process may be slow as the movie is generated in real-time. For a smoother experience, you can save the video to a file using the --save option.

#### Options

- `-h`, `--help`  
  Show this help message and exit.

- `-i RUNDIR`, `--input=RUNDIR`  
  Simulation directory to visualize (default: current directory).

- `--disable_flow`  
  Whether or not to plot the flow arrows.

- `--flow_min=FLOW_MIN`  
  Lower limit of the flow.

- `--flow_max=FLOW_MAX`  
  Upper limit of the flow.

- `--fix_flows`  
  Whether or not the flows are adapted to the data or lie between `flow_min` and `flow_max`.

- `--flow_range=FLOW_RANGE`  
  Log range of the flows in case they are not fixed.

- `--fix_flow_arrow_width`  
  Fix the width of the flow arrows to a constant width.

- `--flow_cmap=FLOW_CMAP`  
  Colormap of the flows.

- `--separate_fission`  
  Whether or not to show arrows also for fission. If not present, hatched areas will be plotted.

- `--fission_minflow=FISSION_MINFLOW`  
  Minimum flow to get indicated as a fission region in case the separate fission flag is not given.

- `--x_min=X_MIN`  
  Lower limit of the mass fraction.

- `--x_max=X_MAX`  
  Upper limit of the mass fraction.

- `--x_cmap=X_CMAP`  
  Colormap of the mass fractions.

- `--disable_abar`  
  Whether or not disabling the indication of Abar.

- `--mass_bins_cmap=MASS_BINS_CMAP`  
  Colormap of the background colors.

- `--disable_magic`  
  Whether or not disabling the indication for the magic number.

- `--additional_plot=ADDITIONAL_PLOT`  
  Whether or not plotting average timescales.

- `--tau_min=TAU_MIN`  
  Lower limit of the average timescales.

- `--tau_max=TAU_MAX`  
  Upper limit of the average timescales.

- `--engen_min=ENGEN_MIN`  
  Lower limit of the Energy.

- `--engen_max=ENGEN_MAX`  
  Upper limit of the Energy.

- `--tracked_min=TRACKED_MIN`  
  Lower limit of the tracked nuclei mass fractions.

- `--tracked_max=TRACKED_MAX`  
  Upper limit of the tracked nuclei mass fractions.

- `--time_min=T_MIN`  
  Lower limit of the time.

- `--time_max=T_MAX`  
  Upper limit of the time.

- `--disable_mainout`  
  Whether or not disabling the mainout plot.

- `--density_min=DENSITY_MIN`  
  Lower limit of the density.

- `--density_max=DENSITY_MAX`  
  Upper limit of the density.

- `--temperature_min=TEMPERATURE_MIN`  
  Lower limit of the temperatures.

- `--temperature_max=TEMPERATURE_MAX`  
  Upper limit of the temperature.

- `--ye_min=YE_MIN`  
  Lower limit of the electron fraction.

- `--ye_max=YE_MAX`  
  Upper limit of the electron fraction.

- `--frame_min=FRAME_MIN`  
  Value of the first frame (default: 1).

- `--frame_max=FRAME_MAX`  
  Value of the last frame (default: end of the simulation).

- `--save`  
  Whether or not saving the movie.

- `--output=OUTPUT_NAME`  
  Output name of the movie.

- `--parallel_save`  
  Whether or not to save the movie in parallel.

- `--parallel_cpus=PARALLEL_CPUS`  
  Number of CPUs to use for parallel saving.

- `--interval=INTERVAL`  
  Interval of the movie (larger value equals slower).

- `--mpirun_path=MPIRUN_PATH`  
  Path of the `mpirun` command to use for parallel saving.
  
  
#### Example

An example output could look like the following:

![Simulation visualization](../../doc/doxygen/figures/winteler_mhd.gif)
