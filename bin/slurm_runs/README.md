## Helper scripts to run WinNet on a cluster

The folder contains scripts to run WinNet on a cluster. In the following we give a brief overview of the files that are contained in this folder.\
Use the scripts on your own risk!

-----

#### [submit_script.sh](submit_script.sh)

**Description**\
When you already have prepared a folder that contains one folder for every trajectory including the WinNet executable and a parameter file (for example with the 
[makerun.py](../../makerun.py.example)), the [submit_script.sh](submit_script.sh) can be used to start several slurm jobs to run WinNet. 
For this you have to modify the header of the script to, e.g., point to your directory in the line\
``` #SBATCH -D /path/to/run ```\
If all necessary modifications have been done, you can start your job by, e.g., running\
``` sbatch --array 1-100 ```\
which will in this case run WinNet with 100 cores. 

-----

#### [winnet_monitor_failed.py](winnet_monitor_failed.py)

**Description**\
In case you have started a run with the [submit_script.sh](submit_script.sh) script, you can monitor the failed runs. When running this script with\
``` python winnet_monitor_failed.py path_to_run```\
it will output a list with all failed runs together with the error code that is described [here](https://nuc-astro.github.io/WinNet/error_codes.html).
An example output could be:

```
          WinNet monitoring     
====================================
| List of failed runs:             | 
|----------------------------------| 
| Name                | Error code | 
| tracer_576724.dat   | W170006    | 
|__________________________________| 
Different error codes: 
W170006
```

-----

#### [winnet_monitor_run.py](winnet_monitor_run.py)

**Description**\
In case you have started a run with the [submit_script.sh](submit_script.sh) script, you can monitor the runs. When running this script with\
``` python winnet_monitor_run.py path_to_run```\
it will output how many runs are still running, how many failed, and an approximation how long it will still take. An example output could be:

```
          WinNet monitoring     
====================================
| Number of runs      |       2024 |
| Running             |          5 | 
| Finished            |       1982 |
| Failed              |          1 | 
| -------------------------------- | 
| Av. time / tracer   |   00:14:07 | 
| Estimated duration  |   01:20:00 | 
| Estimated finish    |   11:15:25 | 
|__________________________________|
```
