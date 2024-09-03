## Example data
This folder contains files to run the examples of WinNet.

A list of all examples can be obtained by running the following terminal command:

```python makerun.py --example```

This will output a list of all possible example cases in the terminal.
Below is the full list of example cases.

<br><br>
### Big Bang:

---

**Description:**

Big Bang nucleosynthesis. The corresponding trajectory is calculated as described in Winteler 2014. The example calculates the element production with a photon to baryon ratio of 6.11e-10.

**Run it with:**
```bash
python makerun.py -p Example_BigBang.par -r Example_BigBang
```

To calculate many runs with different photon to baryon ratios you can run:
```bash
python makerun.py --many --prepare -p Example_BigBang_many.par -r Example_BigBang_many
```

**Literature:**

`-`

---

<br><br>
### Neutron star merger:

---

**Description:**

Run a model with 30 trajectories. The simulation was performed using two neutron stars with 1 solar mass each. The trajectories only included the dynamic ejecta. When running the example, the trajectories are downloaded from [here](https://compact-merger.astro.su.se/downloads_fluid_trajectories.html).  
If you want to run only one trajectory, you can run the following example case:

**Run it with:**
```bash
python makerun.py -p Example_NSM_dyn_ejecta_rosswog.par -r Example_NSM_dyn_ejecta_rosswog
```

Or in case you want to run this trajectory with three different fission fragment distributions:
```bash
python makerun.py -p Example_NSM_dyn_ejecta_fission_rosswog.par -r Example_NSM_dyn_ejecta_fission_rosswog --many --val 1,2,3
```

In case you want to run all available trajectories (30), run the following command. Be aware that this will most likely take more than one hour:
```bash
python makerun.py --many --prepare -p Example_NSM_dyn_ejecta_ns10ns10_rosswog.par -r Example_NSM_dyn_ejecta_ns10ns10_rosswog
```

**Literature:**
- [Korobkin et al. 2012](https://ui.adsabs.harvard.edu/abs/2012MNRAS.426.1940K/abstract)
- [Rosswog et al. 2013](https://ui.adsabs.harvard.edu/abs/2013MNRAS.430.2585R/abstract)
- [Piran et al. 2013](https://ui.adsabs.harvard.edu/abs/2013MNRAS.430.2121P/abstract)

---

**Description:**

Run two trajectories of model LS220-M1.35. One trajectory has a low mass and one is representative for the whole model. The trajectories calculate the dynamic ejecta of a neutron star merger.

**Run it with:**
```bash
python makerun.py --many -p Example_NSM_dyn_ejecta_bovard.par -r Example_NSM_dyn_ejecta_bovard
```

**Literature:**
- [Bovard et al. 2017](https://ui.adsabs.harvard.edu/abs/2017PhRvD..96l4005B/abstract)

---

**Description:**

Representative trajectory (bin 4, 90ms) of the neutrino driven wind of a neutron star merger. 

**Run it with:**
```bash
python makerun.py -p Example_NSM_wind_martin.par -r Example_NSM_wind_martin
```

Or alternatively run all different directional bins (4 tracers).

**Run it with:**
```bash
python makerun.py --many --prepare -p Example_NSM_wind_bins_martin.par -r Example_NSM_wind_bins_martin
```

**Literature:**
- [Perego et al. 2014](https://ui.adsabs.harvard.edu/abs/2014MNRAS.443.3134P/abstract)
- [Martin et al. 2015](https://ui.adsabs.harvard.edu/abs/2015ApJ...813....2M/abstract)

---

**Description:**

Representative trajectories (10) of the viscous/disk ejecta (model S-def) of a neutron star merger. The model has also been investigated in Lippuner et al. 2017.

**Run it with:**
```bash
python makerun.py --many --prepare -p Example_NSM_disk_wu.par -r Example_NSM_disk_wu
```

**Literature:**
- [Wu et al. 2016](https://ui.adsabs.harvard.edu/abs/2016MNRAS.463.2323W/abstract)
- [Lippuner et al. 2017](https://ui.adsabs.harvard.edu/abs/2017MNRAS.472..904L/abstract)

---

<br><br>
### Neutron star black hole merger:

---

**Description:**

Run 20 trajectories of a neutron star black hole merger. The simulation was performed using one neutron star with 1.4 solar masses and a black-hole with 10 solar masses. The trajectories only included the dynamic ejecta. When running the example, the trajectories are downloaded from [here](https://compact-merger.astro.su.se/downloads_fluid_trajectories.html). Due to the amount of calculations, this example case may take a while.

**Run it with:**
```bash
python makerun.py --many --prepare -p Example_NSBH_rosswog.par -r Example_NSBH_rosswog
```

**Literature:**
- [Korobkin et al. 2012](https://ui.adsabs.harvard.edu/abs/2012MNRAS.426.1940K/abstract)
- [Rosswog et al. 2013](https://ui.adsabs.harvard.edu/abs/2013MNRAS.430.2585R/abstract)
- [Piran et al. 2013](https://ui.adsabs.harvard.edu/abs/2013MNRAS.430.2121P/abstract)

---

<br><br>
### Magneto-rotational driven supernova:

---

**Description:**

Neutron-rich trajectory of a magneto-rotational driven supernova. Within this tracer particle, heavy elements get synthesized via the r-process.

**Run it with:**
```bash
python makerun.py -p Example_MRSN_r_process_winteler.par -r Example_MRSN_r_process_winteler
```

Or in case you want to run this trajectory with three different beta-decay rates:
```bash
python makerun.py -p Example_MRSN_r_process_beta_winteler.par -r Example_MRSN_r_process_beta_winteler --many --val beta_decay_marketin.dat,beta_decay_moeller.dat,beta_decay_reaclib.dat
```

**Literature:**
- [Winteler et al. 2012](https://ui.adsabs.harvard.edu/abs/2012ApJ...750L..22W/abstract)

---

**Description:**

Neutron-rich trajectory of a magneto-rotational driven supernova (35OC-Rs). Within this tracer particle, heavy elements get synthesized via the r-process.

**Run it with:**
```bash
python makerun.py -p Example_MRSN_r_process_obergaulinger.par -r Example_MRSN_r_process_obergaulinger
```

To run a less neutron rich trajectory that synthesizes elements by the weak-r process from the same model run:
```bash
python makerun.py -p Example_MRSN_weakr_obergaulinger.par -r Example_MRSN_weakr_obergaulinger
```

A proton-rich trajectory that synthesizes elements via the nu-p process from model 35OC-RO can be run by:
```bash
python makerun.py -p Example_MRSN_nup_process_obergaulinger.par -r Example_MRSN_nup_process_obergaulinger
```

**Literature:**
- [Obergaulinger & Aloy 2017](https://ui.adsabs.harvard.edu/abs/2017MNRAS.469L..43O/abstract)
- [Reichert et al. 2021](https://ui.adsabs.harvard.edu/abs/2021MNRAS.501.5733R/abstract)

---

**Description:**

High entropy trajectory as well as two neutron-rich tracers of a magneto-rotational driven supernova. The high entropy tracer originates in model P, the neutron-rich tracers from model 35OC-Rs_N. In all cases heavy elements get synthesized via the r-process.

**Run it with:**
```bash
python makerun.py --many -p Example_MRSN_r_process_reichert.par -r Example_MRSN_r_process_reichert
```

**Literature:**
- [Aloy & Obergaulinger 2021](https://ui.adsabs.harvard.edu/abs/2021MNRAS.500.4365A/abstract)
- [Obergaulinger & Aloy 2021](https://ui.adsabs.harvard.edu/abs/2021MNRAS.503.4942O/abstract)
- [Reichert et al. 2023](https://ui.adsabs.harvard.edu/abs/2023MNRAS.518.1557R/abstract)

---

<br><br>
### Classical Novae:

---

**Description:**

Trajectory of a classical novae (model ONe5). The trajectory can be downloaded at [Zenodo](https://zenodo.org/record/6474694) (v 1.1.1). Also, the following description is given there: The model follows 1 outburst from the initiation of the accretion stage and through the explosion, expansion, and ejection. It relies on a 1.25 Msun, ONe WD hosting the explosion. For simplicity, it is assumed pre-enriched, accreted material with a composition corresponding to 50% solar (from Lodders 2009) and 50% outer layers of the WD substrate (from Ritossa et al. 1996). The initial luminosity of the WD is assumed to be 10^-2 Lsun, and the mass-accretion rate is 2x10^-10 Msun yr^-1.

**Run it with:**
```bash
python makerun.py -p Example_classical_novae_jose.par -r Example_classical_novae_jose
```

**Literature:**
- [Jose & Hernanz 1998](https://ui.adsabs.harvard.edu/abs/1998ApJ...494..680J/abstract)
- [Jose 2022](https://doi.org/10.5281/zenodo.6474694)

---

<br><br>
### Accreting neutron stars:

---

**Description:**

Trajectory of an X-ray burst. The trajectory originates from a simulation of an accreting neutron star. The nucleosynthesis is dominated by the *rp-process*.

**Run it with:**
```bash
python makerun.py -p Example_xrayburst_schatz.par -r Example_xrayburst_schatz
```

**Literature:**
- [Schatz et al. 2002](https://ui.adsabs.harvard.edu/abs/2001NuPhA.688..150S/abstract)

---

<br><br>
### Regular Core-collapse supernovae:

---

**Description:**

Simple parametric model to calculate *complete Si burning*. The trajectory starts out of NSE and no initial composition is required.

**Run it with:**
```bash
python makerun.py -p Example_CCSN_explosive_burning_parametrized.par -r Example_CCSN_explosive_burning_parametrized
```

**Literature:**

`-`

---

**Description:**

Trajectory assuming a steady state model of a *neutrino driven wind*. The trajectory (CPR2) is publicly available at [here](https://theorie.ikp.physik.tu-darmstadt.de/astro/resources.php).

**Run it with:**
```bash
python makerun.py -p Example_CCSN_wind_bliss.par -r Example_CCSN_wind_bliss
```

**Literature:**
- [Bliss et al. 2018](https://ui.adsabs.harvard.edu/abs/2018ApJ...855..135B/abstract)

---

<br><br>
### Type Ia Supernova:

---

**Description:**

Parametrization of the detonation phase of a Type Ia Supernova. The trajectory undergoes *explosive burning*.

**Run it with:**
```bash
python makerun.py -p Example_type_Ia_meakin.par -r Example_type_Ia_meakin
```

**Literature:**
- [Meakin et al. 2009](https://ui.adsabs.harvard.edu/abs/2009ApJ...693.1188M/abstract)

---

<br><br>
### AGB stars:

---

**Description:**

Carbon 13 and thermal pulse trajectories of a 3 Msun, Z = 0.014 (solar metallicity) star. Elements get synthesized within a *strong s-process*. The two trajectories were calculated with the 1D stellar evolution code MESA. They were accessed at [Zenodo](https://zenodo.org/record/6474686) (v.1.2.1).

**Run it with:**
```bash
python makerun.py --many -p Example_AGB_cescutti.par -r Example_AGB_cescutti
```

**Literature:**
- [Cescutti et al. 2018](https://ui.adsabs.harvard.edu/abs/2018MNRAS.478.4101C/abstract)
- [Cescutti 2022](https://doi.org/10.5281/zenodo.6474686)

---

**Description:**

Representative trajectory of a *weak s-process*. The trajectory was extracted from a 25 Msun, Z = 0.014 (solar metallicity) stellar evolution model. It was chosen because it roughly corresponds to the average weak s-process production in massive stars weighted over the initial mass function. This information as well as the trajectory itself were accessed at [Zenodo](https://zenodo.org/record/6474728) (v 1.1.1).

**Run it with:**
```bash
python makerun.py -p Example_AGB_nishimura.par -r Example_AGB_nishimura
```

**Literature:**
- [Hirschi et al. 2004](https://ui.adsabs.harvard.edu/abs/2004A%26A...425..649H/abstract)
- [Nishimura et al. 2017](https://ui.adsabs.harvard.edu/abs/2017MNRAS.469.1752N/abstract)
- [Pignatari & Hirschi 2022](https://zenodo.org/record/6474728)

---

<br><br>
### Hydrostatic burning:

**Description:**

Hydrostatic *hydrogen burning* for a constant density of 100 g/ccm and different temperatures (between 1e7 K and 4e7 K). To be able to plot the output, the hdf5 output has to be enabled and configured in the Makefile.

**Run it with:**
```bash
python makerun.py -p Example_hydrostatic_hydrogen_burning.par -r Example_hydrostatic_hydrogen_burning --many --prepare --val_min=1 --val_max=4 --val_it=0.1
```

**Literature:**


`-`

---

**Description:**

Hydrostatic *carbon-oxygen burning* for a constant density of 1e9 g/ccm and a temperature of 3 GK.

**Run it with:**
```bash
python makerun.py -p Example_CO_burning.par -r Example_CO_burning
```

**Literature:**

`-`

---

**Description:**

A simple model simulating an *i-process* within a one-zone nuclear reaction network.

**Run it with:**
```bash
python makerun.py --many -p Example_i_process_dardelet.par -r Example_i_process_dardelet
```

**Literature:**
- [Dardelet et al. 2015](https://ui.adsabs.harvard.edu/abs/2015arXiv150505500D/abstract)

---
