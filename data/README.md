## Data folder
This folder contains the data for reaction rates, considered nuclei, nuclear masses, fission fragments, ....

We list the file contents including the corresponding references below:

### alpha_decays.dat
#### Contains:

Name of the alpha-decaying nucleus as well as alpha-decay half-life. The values have been calculated according to:

$$ \log_{10}T_\alpha = (aZ + b)Q_\alpha^{-0.5} + (cZ + d)+h_{log}  $$

with certain fit parameters that are given in the according literature.

#### Relevant parameters
use_alpha_decay_file , alpha_decay_file
#### Literature
[Reichert et al. 2023](https://ui.adsabs.harvard.edu/abs/2023arXiv230507048R/abstract), [Dong & Ren 2005](https://ui.adsabs.harvard.edu/abs/2005EPJA...26...69D/abstract)
#### See also
[Viola & Seaborg 1966](https://www.sciencedirect.com/science/article/abs/pii/0022190266804128)
#### Useful scripts
[create_alpha_decay_file.py](../bin/create_alpha_decay_file.py)

------

### beta_decay_marketin.dat
#### Contains:

File with half lifes and beta-delayed neutron emission probabilities (P0,...P10n), average Q-value, and average energy of
released neutrinos of the decay.

#### Relevant parameters
use_beta_decay_file , beta_decay_file , heating_mode
#### Literature
[Marketin et al. 2016](https://ui.adsabs.harvard.edu/abs/2016PhRvC..93b5805M/abstract)
#### Useful scripts
[convert_marketin_beta_decays.py](../bin/convert_marketin_beta_decays.py)

------

### beta_decay_moeller.dat
#### Contains:

File with half lifes and beta-delayed neutron emission probabilities (P0,...P10n).

#### Relevant parameters
use_beta_decay_file , beta_decay_file
#### Literature
[Moeller et al. 2019](https://ui.adsabs.harvard.edu/abs/2019ADNDT.125....1M/abstract)

------

### beta_decay_reaclib.dat
#### Contains:

File with half lifes and beta-delayed neutron emission probabilities (P0,...P10n).

#### Relevant parameters
use_beta_decay_file , beta_decay_file
#### Literature
[Cyburt et al. 2010](https://ui.adsabs.harvard.edu/abs/2010ApJS..189..240C)

------

### chem_table.dat
#### Contains:

Tabulated chemical potential of electron-positron gas taken from the Helmholtz equation of state.

#### Relevant parameters
chem_pot_file , use_timmes_mue
#### Literature
[Timmes & Arnett 1999](https://ui.adsabs.harvard.edu/abs/1999ApJS..125..277T/abstract),
[Cococubed](https://cococubed.com/code_pages/chemical_potential.shtml)

------

### datafile2.txt
#### Contains:

Tabulated partition functions exceeding 10 GK for the FRDM mass model.

#### Relevant parameters
htpf_file , use_htpf
#### Literature
[Rauscher 2003](https://ui.adsabs.harvard.edu/abs/2003ApJS..147..403R/abstract)

------

### FISS_Mumpower
#### Contains:

Fission fragment distribution for neutron-induced and beta delayed fission.

#### Relevant parameters
fissflag , nfission_file
#### Literature
[Mumpower et al. 2020](https://ui.adsabs.harvard.edu/abs/2020PhRvC.101e4607M/abstract)

------

### fissionrates_frdm
#### Contains:

Fission rates for the FRDM mass model in Reaclib file format. Neutron-induced fission
is given in [Panov et al. 2010](https://ui.adsabs.harvard.edu/abs/2010A%26A...513A..61P/abstract),
beta-delayed fission in [Panov et al. 2005](https://ui.adsabs.harvard.edu/abs/2005NuPhA.747..633P/abstract),
and spontaneous fission was calculated using the semi-empirical formula given in
[Khuyagbaatar 2020](https://ui.adsabs.harvard.edu/abs/2020NuPhA100221958K/abstract) with the
fission barriers of [Moeller et al. 2015](https://journals.aps.org/prc/supplemental/10.1103/PhysRevC.91.024310).

#### Relevant parameters
fissflag , fission_rates
#### Literature
[Panov et al. 2005](https://ui.adsabs.harvard.edu/abs/2005NuPhA.747..633P/abstract),
[Panov et al. 2010](https://ui.adsabs.harvard.edu/abs/2010A%26A...513A..61P/abstract)
[Moeller et al. 2015](https://journals.aps.org/prc/supplemental/10.1103/PhysRevC.91.024310)
[Khuyagbaatar 2020](https://ui.adsabs.harvard.edu/abs/2020NuPhA100221958K/abstract)
#### Useful scripts
[create_spontaneous_fission_file.py](../bin/create_spontaneous_fission_file.py)

------

### frdm_sn.dat
#### Contains:

Table with neutron separation energies with the FRDM mass model.

#### Relevant parameters
calc_nsep_energy , nsep_energies_file

------

### stable_isotopes.txt
#### Contains:
List that contains mass number, atomic number, and neutron number of stable isotopes.

------

### neunucleons.dat
#### Contains:
Tabulated neutrino reactions as well as average energies of the absorped
neutrinos on neutrons and protons.

#### Relevant parameters
nuflag , nunucleo_rates_file
#### Literature
[Burrows et al 2006](https://ui.adsabs.harvard.edu/abs/2006NuPhA.777..356B/abstract)
[Horowitz 2002](https://ui.adsabs.harvard.edu/abs/2002PhRvD..65d3001H/abstract)
#### See also
[Burrows & Thompson 2002](https://ui.adsabs.harvard.edu/abs/2002astro.ph.11404B/abstract)
#### Useful scripts
[create_neutrino_nucleon_file.py](../bin/create_spontaneous_fission_file.py)

------

### nu_channels
#### Contains:
Channels of neutrino reactions. Both, charged current and neutral current channels are included.

#### Relevant parameters
nuflag , nuchannel_file
#### Literature
[Sieverding et al. 2018](https://ui.adsabs.harvard.edu/abs/2018ApJ...865..143S/abstract)

------

### nucross.dat
#### Contains:
Neutrino reactions on heavy nuclei. Both, charged current and neutral current channels are included.

#### Relevant parameters
nuflag , nurates_file
#### Literature
[Sieverding et al. 2018](https://ui.adsabs.harvard.edu/abs/2018ApJ...865..143S/abstract)

------

### nu_loss_data.dat
#### Contains:
Average energy of neutrinos in MeV that are produced in the beta-decay.

#### Relevant parameters
use_neutrino_loss_file , neutrino_loss_file , heating_mode
#### Literature
[Brown et al. 2018](https://ui.adsabs.harvard.edu/abs/2018NDS...148....1B/abstract)
#### Useful scripts
[create_neutrino_loss_file.py](../bin/create_neutrino_loss_file.py)

------

### rateseff.out
#### Contains:
Tabulated weak reactions.

#### Relevant parameters
iwformat , weak_rates_file, temp_reload_exp_weak_rates
#### Literature
[Langanke & Martinez-Pinedo 2001](https://ui.adsabs.harvard.edu/abs/2001ADNDT..79....1L/abstract)
#### See also
[Fuller et al. 1985](https://ui.adsabs.harvard.edu/abs/1985ApJ...293....1F/abstract),
[Oda et al. 1994](https://ui.adsabs.harvard.edu/abs/1994ADNDT..56..231O/abstract)

------

### Reaclib_18_9_20
#### Contains:
Nuclear reaction rates in reaclib file format.

#### Relevant parameters
reaclib_file
#### Literature
[Cyburt et al. 2010](https://ui.adsabs.harvard.edu/abs/2010ApJS..189..240C)
#### See also
https://reaclib.jinaweb.org/

------

### theoretical_weak_rates.dat
#### Contains:
Tabulated weak reactions (electron-/positron- captures and $\beta$-decays) from different literatures
The reaction rates have been accessed via the [Weak rate library](https://groups.nscl.msu.edu/charge_exchange/weakrates.html) and have been converted into
log \<ft\> format.

#### Relevant parameters
iwformat , weak_rates_file, temp_reload_exp_weak_rates
#### Literature
[Fuller et al. 1985](https://ui.adsabs.harvard.edu/abs/1985ApJ...293....1F/abstract),
[Oda et al. 1994](https://ui.adsabs.harvard.edu/abs/1994ADNDT..56..231O/abstract),
[Langanke & Martinez-Pinedo 2001](https://ui.adsabs.harvard.edu/abs/2001ADNDT..79....1L/abstract),
[Pruet & Fuller 2003](https://ui.adsabs.harvard.edu/abs/2003ApJS..149..189P/abstract),
[Suzuki, Toki & Nomoto 2016](https://ui.adsabs.harvard.edu/abs/2016ApJ...817..163S/abstract)
#### See also
[Weak rate library](https://groups.nscl.msu.edu/charge_exchange/weakrates.html),
[tw_rate_module.f90](../src/tw_rate_module.f90)

------

### winvne_v2.0.dat
#### Contains:
Nuclear properties, such as the spin of the ground state, the mass excess, neutron and proton numbers and tabulated partition functions.

#### Relevant parameters
isotopes_file
#### Literature
[Cyburt et al.2010](https://ui.adsabs.harvard.edu/abs/2010ApJS..189..240C)
#### See also
https://reaclib.jinaweb.org/