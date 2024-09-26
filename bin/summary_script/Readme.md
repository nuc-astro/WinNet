# Simulation Data Summary Tool

This Python script summarizes the data from simulation runs that have been run with the --many option. The result is saved in an HDF5 file. It offers flexible options for customization and supports a variety of output formats.
The output of the individual runs is (linearly) interpolated to one common time grid. If the trajectory starts later or ends earlier than this time grid, the data in the summary is filled with NaNs. 

## Usage

To run the script, use:

```bash
python summarize.py -i <rundir> [options]
```

## Options

- **`-i`, `--input`**  
  Specifies the simulation directory to summarize.  
  **Default**: `.` (current directory)

- **`-o`, `--output`**  
  Specifies the output path for the summary file.  
  **Default**: `./summary.hdf5`

- **`-b`, `--buf`**  
  Buffer size before writing data to the output file.  
  **Default**: `500`

- **`-f`, `--force`**  
  Forces overwriting the output file if it already exists.  
  **Default**: `False`

- **`-v`, `--verbose`**  
  Enables verbose output, which is logged to `debug.log`.  
  **Default**: `False`

- **`--time_file`**  
  Specifies the path to a file containing the time grid in seconds.  
  **Default**: `None`

- **`--time_final`**  
  Specifies the final time for the time grid. This is only used if `--time_file` is not given.  
  **Default**: Read from the template file.

- **`--time_initial`**  
  Specifies the initial time for the time grid. This is only used if `--time_file` is not given.  
  **Default**: `1e-5`

- **`--time_number`**  
  Specifies the number of time steps for the time grid. This is only used if `--time_file` is not given.  
  **Default**: `200`

- **`--sunet_path`**  
  Specifies the path to the `sunet` file.  
  **Default**: Read from the template.

### Output Disabling Options

You can choose to disable certain parts of the summary process using the following options:

- **`--disable_mainout`**  
  Disables summarizing the `mainout` output.  
  **Default**: `False`

- **`--disable_energy`**  
  Disables summarizing the `energy` output.  
  **Default**: `False`

- **`--disable_timescales`**  
  Disables summarizing the `timescales` output.  
  **Default**: `False`

- **`--disable_tracked_nuclei`**  
  Disables summarizing the `tracked_nuclei` output.  
  **Default**: `False`

- **`--disable_snapshots`**  
  Disables summarizing the `snapshots` output.  
  **Default**: `False`
  
## Hdf5 output

The HDF5 file contains several key datasets based on the simulation results. 
Note that there are tools like HDF compass to visualize the content of the Hdf5 file. 
Below is a summary of the main entries in the HDF5 file:

1. **`finab/`**:
   - **`A`**: Array of atomic mass numbers (A) for the nuclei.
   - **`Z`**: Array of proton numbers (Z) for the nuclei.
   - **`N`**: Array of neutron numbers (N) for the nuclei.
   - **`Y`**: Abundance of the nuclei at the end of the simulation.
   - **`X`**: Mass fraction of the nuclei at the end of the simulation.

2. **`run_names`**: Array of run names (strings) that were summarized during the process.

3. **`run_ids`**: Array of run IDs, extracted from run directory names, which typically contain numeric identifiers.

4. **`mainout/`** (if not disabled):
   - **`time`**: Time grid for the "mainout" data.
   - **Other Keys**: Various interpolated quantities from the "mainout" data, excluding iteration, time, and nucleosynthesis keys (like A, Z, N).

5. **`energy/`** (if not disabled):
   - **`time`**: Time grid for the "energy" data.
   - **Other Keys**: Energy-related quantities, interpolated over the time grid.

6. **`timescales/`** (if not disabled):
   - **`time`**: Time grid for the "timescales" data.
   - **Other Keys**: Timescale-related quantities, interpolated over the time grid.

7. **`tracked_nuclei/`** (if not disabled):
   - **`time`**: Time grid for the "tracked_nuclei" data.
   - **Other Keys**: Data related to the abundance and properties of tracked nuclei, interpolated over the time grid.

8. **`snapshot/`** (if applicable and not disabled):
   - **`time`**: Custom snapshot time grid (if using custom snapshots).
   - **`Y`**: Data corresponding to the abundances at snapshot times.

  

## Examples

### Example 1: Summarizing a simulation

```bash
python summarize.py -i runs/test
```

Summarizes the simulation in the directory `runs/test` and outputs the summary to the default path `./summary.hdf5`.

### Example 2: Force overwrite and set buffer size

```bash
python summarize.py -i runs/test -o ./output/summary.hdf5 -f -b 1000
```

Summarizes the simulation in `runs/test`, forces the overwrite of `./output/summary.hdf5`, and sets the buffer size to `1000`.

### Example 3: Verbose output and custom time file

```bash
python summarize.py -i runs/test -v --time_file custom_time.txt
```

Summarizes the simulation in `runs/test` with verbose logging enabled (output saved to `debug.log`) and uses `custom_time.txt` for the time grid.

## Dependencies

Make sure you have the following libraries installed:

- `numpy`
- `h5py`
- `tqdm`
