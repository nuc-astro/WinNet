<img src="/doc/doxygen/doxygen_logo/winnet_logo_light_theme.png#gh-light-mode-only" width="40%">
<img src="/doc/doxygen/doxygen_logo/winnet_logo_dark_theme.png#gh-dark-mode-only"  width="40%">

A flexible, multi-purpose, single-zone nuclear reaction network.


## Documentation
See [WinNet-documentation](https://nuc-astro.github.io/WinNet/) for documentation and further information.


## Literature
[![arXiv](https://img.shields.io/badge/arxiv-2305.07048-b31b1b)](https://arxiv.org/abs/2305.07048) 
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8220549.svg)](https://doi.org/10.5281/zenodo.8220549)



Please read and cite [Reichert et al. 2023](https://ui.adsabs.harvard.edu/abs/2023arXiv230507048R/abstract) when you use this reaction network in a publication. 
Furthermore give credit to the data that you are using (reaction rates, thermodynamic conditions, equation of state, etc.). 
For this you can be guided by the example cases that contain the relevant literature in the header of their parameter files in the [par/](par/) folder. 
Additionally, the [Readme](data/README.md) in the [data/](data/) folder contains relevant information about the references of input files. The origin of
the thermodynamic trajectories is given in the list of examples when running
```python makerun.py --example```
or in the Readme of the data folders of the examples (e.g., [Readme](data/Example_data/Example_MRSN_r_process_winteler/README.md)).


## License
WinNet is available as open source under the terms of the revised BSD 3-Clause License. See the LICENSE file for more details.


## Acknowledgments
This work makes use of the following packages:
- The quadpack integration modules, accessed at https://netlib.org/quadpack/index.html
- The MINPACK-I Fortran 90 package accessed at https://people.sc.fsu.edu/~jburkardt/f_src/fsolve/fsolve.html. The code is distributed under the GNU LGPL license. The license file can be accessed [here](src/external_tools/LICENSE).
- The Timmes Equation of State accessed at https://cococubed.com/code_pages/eos.shtml that was presented in [Timmes & Swesty 2000](https://ui.adsabs.harvard.edu/abs/2000ApJS..126..501T/abstract).
- An interpolation routine for electron chemical potentials of F. Timmes accessed at https://cococubed.com/code_pages/chemical_potential.shtml.
- Fortran routines from https://cococubed.com/code_pages/nuloss.shtml that incorporate the analytic expressions of the thermal neutrino energies of [Itoh et al 1996](https://ui.adsabs.harvard.edu/abs/1996ApJS..102..411I/abstract).
- Other used data such as reaction rates and thermodynamic conditions for the example files are stated in the parameter files in the [par/](par/) folder.
