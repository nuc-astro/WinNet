### Main s-process example
When using the data in this folder please read and cite the following literature.
#### Literature
- [Cescutti G., Hirschi R., Nishimura N., Hartogh J. W. den., Rauscher T., Murphy A. S. J., Cristallo S., 2018, MNRAS, 478, 4101. doi:10.1093/mnras/sty1185](https://ui.adsabs.harvard.edu/abs/2018MNRAS.478.4101C/abstract)
- [Cescutti G. (2022). Main s-process (1.2.1) [Data set]. Zenodo. https://doi.org/10.5281/zenodo.6474686](https://zenodo.org/record/6474686)

#### Access link
https://zenodo.org/record/6474686

#### Additional information from the providers of the trajectory:
The trajectory was used in Cescutti et al. (2018), and was extracted from a 3 Msun, Z = 0.014 (solar metallicity) stellar evolution model, calculate with 1D stellar evolution code MESA, revision number 6208 (Paxton et al. 2011).
The trajectory was taken from the 13C-pocket following the 6th Thermal-Pulse (TP)

The temperature and density profiles of the trajectory are in the corresponding file (trajectory_C13pocket: columns  time [s]   temperature [T9]  density [g/cm^3]). The abundances are in the attached file “iniab_C13pocket.

In the paper we run other 2 models in which we simply multiplied (c13 seed x 2) and  (c13 seed x 0.5).

The trajectory for TP conditions is also provided. The trajectory lasts for 1 yr, with a constant temperature of 0.245 GK and a constant density of 5 10^3 g cm−3. The initial abundances are summarized in Table 2 for light elements; for the other elements, the final abundances of
the standard 13C-pocket were diluted by a factor of 20 to take into account the diluting effect of the PDCZ.
The temperature and density profiles of the trajectory are in the corresponding file (trajectory_TP: columns  time [s]   temperature [T9]  density [g/cm^3]) The abundances are in iniab_TP).

These trajectories does not follow exactly what happens in real stars, but provide  the conditions that lead to an s-process production similar to that predicted  using full stellar models. Indeed it covers a single 13C-pocket and thermal-pulse contrary to the up to thousand computed in full stellar models. Therefore, these trajectories can be used for nuclear sensitivity studies, but they cannot be used to compute stellar yields.

Concerning the final abundances, please note that the output is the number density of nuclei whereas the initial composition is in mass fraction.
