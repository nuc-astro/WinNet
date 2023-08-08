## Folder with useful scripts and binaries

This folder contains useful scripts to change input reaction rates and for running automatic test cases. Also the WinNet binary will be created here. 
In the following we briefly list the content of the folder together with a small description.\
Use all scripts on your own risk!

------

#### [convert_marketin_beta_decays.py](bin/convert_marketin_beta_decays.py)

**Description**\
Convert beta decay file of Marketin into a WinNet readable format.
The [Marketin et al. 2016](https://ui.adsabs.harvard.edu/abs/2016PhRvC..93b5805M/abstract) file 
can be accessed [here](https://journals.aps.org/prc/supplemental/10.1103/PhysRevC.93.025805).

------

#### [create_alpha_decay_file.py](bin/create_alpha_decay_file.py)

**Description**\
Creates a file with alpha decay half-lifes. The half-lifes
are calculated by the Viola-Seaborg formula. The fit parameters
are partly taken from [Dong & Ren 2005](https://ui.adsabs.harvard.edu/abs/2005EPJA...26...69D/abstract) 
and have been partly fitted to experimental available data of the reaclib.
To create this file, a winvn file is needed.

------

#### [create_neutrino_loss_file.py](bin/create_neutrino_loss_file.py)

**Description**\
Script to download the neutrino loss data from the ENDSF database and to put it into a WinNet readable file format.
The file can be used to give more precise information on how much energy is radiated away by neutrinos in case that
nuclear heating is enabled. A description about the data API can be found [here](https://www-nds.iaea.org/relnsd/vcharthtml/api_v0_guide.html).

------

#### [create_neutrino_nucleon_file.py](bin/create_neutrino_nucleon_file.py)

**Description**\
Script to calculate the electron neutrino and electron antineutrino cross sections on nucleons according to 
[Burrows et al. 2006](https://ui.adsabs.harvard.edu/abs/2006NuPhA.777..356B/abstract) with weak magnetism and recoil corrections as in 
[Horowitz et al. 2002](https://ui.adsabs.harvard.edu/abs/2002PhRvD..65d3001H/abstract). Furthermore, it calculates the average energy of
the absorbed neutrino to consider it in the case that nuclear heating is enabled.

------

#### [create_spontaneous_fission_file.py](bin/create_spontaneous_fission_file.py)

**Description**\
This file will create a file that contains entries for spontaneous fission reactions. These entries have to be copied
into the fission rates file (e.g., fissionrates_frdm). To create the file here we use the fission barriers from
[Moeller et al. 2015](https://ui.adsabs.harvard.edu/abs/2015PhRvC..91b4310M/abstract) that can be downloaded 
[here](https://journals.aps.org/prc/supplemental/10.1103/PhysRevC.91.024310) and
the FRDM fit of [Khuyagbaatar 2020](https://ui.adsabs.harvard.edu/abs/2020NuPhA100221958K/abstract).
In contrast to [Khuyagbaatar 2020](https://ui.adsabs.harvard.edu/abs/2020NuPhA100221958K/abstract), we also include odd nuclei, 
not only even ones. Additionally we add experimentally determined rates that can be accessed from the 
[ENDFS database](nds.iaea.org/relnsd/v1/data?/fields=ground_states&nuclides=all).
It is in principle possible to also create spontaneous fission rates for the ETFSI model, but the
here used fit has to be changed and the correct fission barriers have to be used.

------

#### [examplecase_class.py](bin/examplecase_class.py)

**Description**\
Script to deal with the example cases of WinNet. New example cases have to be added here. The script is used within the [makerun.py](makerun.py) 
when executing\
``` python makerun.py --example ```

------

#### [nucleus_example.py](bin/nucleus_example.py)

**Description**\
Example on how to use the helper script to deal with nuclei in python that is contained in [bin/class_files/nucleus_class.py](bin/class_files/nucleus_class.py). 
With this script it is possible to create a class that contains basic information of a nucleus such as neutron number, proton number, and mass number.

------

#### [reaclib_example.py](bin/reaclib_example.py)

**Description**\
Example on how to use the helper script to deal with reaclib reaction rates in python that is contained in [bin/class_files/reaclib_class.py](bin/class_files/reaclib_class.py). 
With this script it is possible to create a class that contains basic information of a reaclib file. Also analyzing the reaction rates in terms of bad rate fits or erroneous 
entries is possible. Additionally one can merge two different reaclib files.

------

#### [run_tests.py](bin/run_tests.py)

**Description**\
Script that is called when running the automatic test cases with\
``` make tests ```\
It contains the basic comparison operations that are applied for the tests that are contained in the [test](test) folder.

------

#### [testcase_class.py](bin/testcase_class.py)

**Description**\
Helper script for [run_tests.py](bin/run_tests.py) script.

------

#### [winnet_example.py](winnet_example.py)

**Description**\
Example on how to use the helper script to deal with a WinNet run in python that is contained in [bin/class_files/winnet_class.py](bin/class_files/winnet_class.py). 
With this script it is possible to create a class that contains basic information of a WinNet run.
