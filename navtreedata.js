/*
 @licstart  The following is the entire license notice for the JavaScript code in this file.

 The MIT License (MIT)

 Copyright (C) 1997-2020 by Dimitri van Heesch

 Permission is hereby granted, free of charge, to any person obtaining a copy of this software
 and associated documentation files (the "Software"), to deal in the Software without restriction,
 including without limitation the rights to use, copy, modify, merge, publish, distribute,
 sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in all copies or
 substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING
 BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
 DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

 @licend  The above is the entire license notice for the JavaScript code in this file
*/
var NAVTREE =
[
  [ "WinNet", "index.html", [
    [ "First steps", "index.html", [
      [ "Getting started", "index.html#getting_started", [
        [ "Prerequisits", "index.html#prerequisits", null ]
      ] ],
      [ "Starting to run the code", "index.html#start_running", [
        [ "The parameter file", "index.html#parameter", null ],
        [ "Example cases", "index.html#example_cases", null ],
        [ "Monitoring simulations", "index.html#monitoring_simulations", null ],
        [ "Compiling a local version of the documentation", "index.html#compile_docs", null ],
        [ "List of publications", "index.html#publications", null ]
      ] ]
    ] ],
    [ "Coding guidelines", "coding_conventions.html", null ],
    [ "GitHub Packages", "DockerUsage.html", [
      [ "Prerequisites", "DockerUsage.html#prerequisites", null ],
      [ "Steps to Use the Dockerfile", "DockerUsage.html#steps", [
        [ "Coding conventions", "coding_conventions.html#autotoc_md0", null ],
        [ "Usage of git-hooks", "coding_conventions.html#autotoc_md1", null ],
        [ "Doxygen styling conventions", "coding_conventions.html#autotoc_md2", null ],
        [ "Pushing changes to master branch", "coding_conventions.html#autotoc_md3", null ],
        [ "Authenticate Docker to GitHub Packages", "DockerUsage.html#authenticate", null ],
        [ "Pull the Docker Image", "DockerUsage.html#pull", null ],
        [ "Run the Docker Container", "DockerUsage.html#run", null ]
      ] ],
      [ "Example", "DockerUsage.html#example", null ],
      [ "Additional Notes", "DockerUsage.html#additional_notes", null ]
    ] ],
    [ "Error codes", "error_codes.html", null ],
    [ "How it works", "how_it_works.html", null ],
    [ "Input files", "input_files.html", null ],
    [ "Linux installation", "linux_install.html", null ],
    [ "Mac installation", "mac_install.html", null ],
    [ "Windows installation", "windows_install.html", null ],
    [ "Known issues", "known_issues.html", null ],
    [ "Parameters", "parameters.html", [
      [ "adapt_stepsize_maxcount", "parameters.html#adapt_stepsize_maxcount", null ],
      [ "alpha_decay_file", "parameters.html#alpha_decay_file", null ],
      [ "alpha_decay_ignore_all", "parameters.html#alpha_decay_ignore_all", null ],
      [ "alpha_decay_src_ignore", "parameters.html#alpha_decay_src_ignore", null ],
      [ "alpha_decay_zmax", "parameters.html#alpha_decay_zmax", null ],
      [ "alpha_decay_zmin", "parameters.html#alpha_decay_zmin", null ],
      [ "beta_decay_file", "parameters.html#beta_decay_file", null ],
      [ "beta_decay_src_ignore", "parameters.html#beta_decay_src_ignore", null ],
      [ "bfission_file", "parameters.html#bfission_file", null ],
      [ "calc_nsep_energy", "parameters.html#calc_nsep_energy", null ],
      [ "chem_pot_file", "parameters.html#chem_pot_file", null ],
      [ "custom_snapshots", "parameters.html#custom_snapshots", null ],
      [ "detailed_balance_src_ignore", "parameters.html#detailed_balance_src_ignore", null ],
      [ "detailed_balance_src_q_reac", "parameters.html#detailed_balance_src_q_reac", null ],
      [ "detailed_balance_src_q_winvn", "parameters.html#detailed_balance_src_q_winvn", null ],
      [ "Enue", "parameters.html#Enue", null ],
      [ "Enuebar", "parameters.html#Enuebar", null ],
      [ "Enux", "parameters.html#Enux", null ],
      [ "Enuxbar", "parameters.html#Enuxbar", null ],
      [ "expansiontype", "parameters.html#expansiontype", null ],
      [ "extrapolation_width", "parameters.html#extrapolation_width", null ],
      [ "final_dens", "parameters.html#final_dens", null ],
      [ "final_temp", "parameters.html#final_temp", null ],
      [ "final_time", "parameters.html#final_time", null ],
      [ "fissflag", "parameters.html#fissflag", null ],
      [ "fission_format_beta_delayed", "parameters.html#fission_format_beta_delayed", null ],
      [ "fission_format_n_induced", "parameters.html#fission_format_n_induced", null ],
      [ "fission_format_spontaneous", "parameters.html#fission_format_spontaneous", null ],
      [ "fission_frag_beta_delayed", "parameters.html#fission_frag_beta_delayed", null ],
      [ "fission_frag_missing", "parameters.html#fission_frag_missing", null ],
      [ "fission_frag_n_induced", "parameters.html#fission_frag_n_induced", null ],
      [ "fission_frag_spontaneous", "parameters.html#fission_frag_spontaneous", null ],
      [ "fission_rates_beta_delayed", "parameters.html#fission_rates_beta_delayed", null ],
      [ "fission_rates_n_induced", "parameters.html#fission_rates_n_induced", null ],
      [ "fission_rates_spontaneous", "parameters.html#fission_rates_spontaneous", null ],
      [ "flow_every", "parameters.html#flow_every", null ],
      [ "freeze_rate_temp", "parameters.html#freeze_rate_temp", null ],
      [ "gear_cFactor", "parameters.html#gear_cFactor", null ],
      [ "gear_eps", "parameters.html#gear_eps", null ],
      [ "gear_escale", "parameters.html#gear_escale", null ],
      [ "gear_ignore_adapt_stepsize", "parameters.html#gear_ignore_adapt_stepsize", null ],
      [ "gear_nr_eps", "parameters.html#gear_nr_eps", null ],
      [ "gear_nr_maxcount", "parameters.html#gear_nr_maxcount", null ],
      [ "gear_nr_mincount", "parameters.html#gear_nr_mincount", null ],
      [ "gear_timestep_max", "parameters.html#gear_timestep_max", null ],
      [ "h_custom_snapshots", "parameters.html#h_custom_snapshots", null ],
      [ "h_engen_every", "parameters.html#h_engen_every", null ],
      [ "h_engen_detailed", "parameters.html#h_engen_detailed", null ],
      [ "h_finab", "parameters.html#h_finab", null ],
      [ "h_flow_every", "parameters.html#h_flow_every", null ],
      [ "h_mainout_every", "parameters.html#h_mainout_every", null ],
      [ "h_nu_loss_every", "parameters.html#h_nu_loss_every", null ],
      [ "h_snapshot_every", "parameters.html#h_snapshot_every", null ],
      [ "h_timescales_every", "parameters.html#h_timescales_every", null ],
      [ "h_track_nuclei_every", "parameters.html#h_track_nuclei_every", null ],
      [ "heating_density", "parameters.html#heating_density", null ],
      [ "heating_mode", "parameters.html#heating_mode", null ],
      [ "heating_frac", "parameters.html#heating_frac", null ],
      [ "htpf_file", "parameters.html#htpf_file", null ],
      [ "initemp_cold", "parameters.html#initemp_cold", null ],
      [ "initemp_hot", "parameters.html#initemp_hot", null ],
      [ "initial_stepsize", "parameters.html#initial_stepsize", null ],
      [ "interp_mode", "parameters.html#interp_mode", null ],
      [ "isotopes_file", "parameters.html#isotopes_file", null ],
      [ "iwformat", "parameters.html#iwformat", null ],
      [ "iwinterp", "parameters.html#iwinterp", null ],
      [ "Le", "parameters.html#Le", null ],
      [ "Lebar", "parameters.html#Lebar", null ],
      [ "Lx", "parameters.html#Lx", null ],
      [ "Lxbar", "parameters.html#Lxbar", null ],
      [ "mainout_every", "parameters.html#mainout_every", null ],
      [ "net_source", "parameters.html#net_source", null ],
      [ "neutrino_loss_file", "parameters.html#neutrino_loss_file", null ],
      [ "neutrino_mode", "parameters.html#neutrino_mode", null ],
      [ "nfission_file", "parameters.html#nfission_file", null ],
      [ "nr_maxcount", "parameters.html#nr_maxcount", null ],
      [ "nr_mincount", "parameters.html#nr_mincount", null ],
      [ "nr_tol", "parameters.html#nr_tol", null ],
      [ "nrdiag_every", "parameters.html#nrdiag_every", null ],
      [ "nse_calc_every", "parameters.html#nse_calc_every", null ],
      [ "nse_delt_t9min", "parameters.html#nse_delt_t9min", null ],
      [ "nse_descend_t9start", "parameters.html#nse_descend_t9start", null ],
      [ "nse_max_it", "parameters.html#nse_max_it", null ],
      [ "nse_nr_tol", "parameters.html#nse_nr_tol", null ],
      [ "nse_solver", "parameters.html#nse_solver", null ],
      [ "nsep_energies_file", "parameters.html#nsep_energies_file", null ],
      [ "nsetemp_cold", "parameters.html#nsetemp_cold", null ],
      [ "nsetemp_hot", "parameters.html#nsetemp_hot", null ],
      [ "nuchannel_file", "parameters.html#nuchannel_file", null ],
      [ "nuflag", "parameters.html#nuflag", null ],
      [ "nu_loss_every", "parameters.html#nu_loss_every", null ],
      [ "nunucleo_rates_file", "parameters.html#nunucleo_rates_file", null ],
      [ "nurates_file", "parameters.html#nurates_file", null ],
      [ "out_every", "parameters.html#out_every", null ],
      [ "prepared_network_path", "parameters.html#prepared_network_path", null ],
      [ "reaclib_file", "parameters.html#reaclib_file", null ],
      [ "read_initial_composition", "parameters.html#read_initial_composition", null ],
      [ "rho_analytic", "parameters.html#rho_analytic", null ],
      [ "Rkm_analytic", "parameters.html#Rkm_analytic", null ],
      [ "screening_mode", "parameters.html#screening_mode", null ],
      [ "seed_file", "parameters.html#seed_file", null ],
      [ "seed_format", "parameters.html#seed_format", null ],
      [ "sfission_file", "parameters.html#sfission_file", null ],
      [ "snapshot_every", "parameters.html#snapshot_every", null ],
      [ "snapshot_file", "parameters.html#snapshot_file", null ],
      [ "solver", "parameters.html#solver", null ],
      [ "t_analytic", "parameters.html#t_analytic", null ],
      [ "T9_analytic", "parameters.html#T9_analytic", null ],
      [ "tabulated_rates_file", "parameters.html#tabulated_rates_file", null ],
      [ "tabulated_temperature_file", "parameters.html#tabulated_temperature_file", null ],
      [ "temp_reload_exp_weak_rates", "parameters.html#temp_reload_exp_weak_rates", null ],
      [ "termination_criterion", "parameters.html#termination_criterion", null ],
      [ "timescales_every", "parameters.html#timescales_every", null ],
      [ "timestep_hydro_factor", "parameters.html#timestep_hydro_factor", null ],
      [ "timestep_max", "parameters.html#timestep_max", null ],
      [ "timestep_factor", "parameters.html#timestep_factor", null ],
      [ "timestep_Ymin", "parameters.html#timestep_Ymin", null ],
      [ "timestep_traj_limit", "parameters.html#timestep_traj_limit", null ],
      [ "top_engen_every", "parameters.html#top_engen_every", null ],
      [ "track_nuclei_every", "parameters.html#track_nuclei_every", null ],
      [ "track_nuclei_file", "parameters.html#track_nuclei_file", null ],
      [ "trajectory_file", "parameters.html#trajectory_file", null ],
      [ "trajectory_format", "parameters.html#trajectory_format", null ],
      [ "trajectory_mode", "parameters.html#trajectory_mode", null ],
      [ "use_alpha_decay_file", "parameters.html#use_alpha_decay_file", null ],
      [ "use_beta_decay_file", "parameters.html#use_beta_decay_file", null ],
      [ "use_detailed_balance", "parameters.html#use_detailed_balance", null ],
      [ "use_detailed_balance_q_reac", "parameters.html#use_detailed_balance_q_reac", null ],
      [ "use_htpf", "parameters.html#use_htpf", null ],
      [ "use_neutrino_loss_file", "parameters.html#use_neutrino_loss_file", null ],
      [ "use_prepared_network", "parameters.html#use_prepared_network", null ],
      [ "use_tabulated_rates", "parameters.html#use_tabulated_rates", null ],
      [ "use_thermal_nu_loss", "parameters.html#use_thermal_nu_loss", null ],
      [ "use_timmes_mue", "parameters.html#use_timmes_mue", null ],
      [ "weak_rates_file", "parameters.html#weak_rates_file", null ],
      [ "Ye_analytic", "parameters.html#Ye_analytic", null ]
    ] ],
    [ "Automatic testcases", "atests.html", [
      [ "Neutron decay", "atests.html#neutron_decay", null ],
      [ "The decay chain of Ni56", "atests.html#ni_decay", null ],
      [ "(n,gamma)(gamma,n) equilibrium", "atests.html#n_gamma", null ],
      [ "Alpha network", "atests.html#alphanetwork", null ],
      [ "Big Bang nucleosynthesis", "atests.html#bigbang", null ],
      [ "Nuclear statistical equilibrium", "atests.html#nse", null ],
      [ "Neutrino reactions", "atests.html#neutrino_reactions", null ],
      [ "Beta-delayed fission", "atests.html#bdelayed_fission", null ],
      [ "Beta-decay format", "atests.html#bdecay_format", null ],
      [ "Analytic trajectory mode", "atests.html#parametric", null ],
      [ "Switch evolution mode", "atests.html#switch_evolution", null ],
      [ "Tabulated rates", "atests.html#tabulated_rates", null ],
      [ "Theoretical weak rates", "atests.html#twr", null ],
      [ "Screening corrections", "atests.html#screening", null ],
      [ "Adiabatic expansion", "atests.html#expansion", null ]
    ] ],
    [ "Todo List", "todo.html", null ],
    [ "Modules", "namespaces.html", [
      [ "Modules List", "namespaces.html", "namespaces_dup" ],
      [ "Module Members", "namespacemembers.html", [
        [ "All", "namespacemembers.html", "namespacemembers_dup" ],
        [ "Functions/Subroutines", "namespacemembers_func.html", "namespacemembers_func" ],
        [ "Variables", "namespacemembers_vars.html", "namespacemembers_vars" ]
      ] ]
    ] ],
    [ "Data Types List", "annotated.html", [
      [ "Data Types List", "annotated.html", "annotated_dup" ],
      [ "Class Hierarchy", "hierarchy.html", "hierarchy" ],
      [ "Data Fields", "functions.html", [
        [ "All", "functions.html", "functions_dup" ],
        [ "Functions/Subroutines", "functions_func.html", null ],
        [ "Variables", "functions_vars.html", "functions_vars" ]
      ] ]
    ] ],
    [ "Files", "files.html", [
      [ "File List", "files.html", "files_dup" ],
      [ "File Members", "globals.html", [
        [ "All", "globals.html", null ],
        [ "Functions/Subroutines", "globals_func.html", null ],
        [ "Macros", "globals_defs.html", null ]
      ] ]
    ] ]
  ] ]
];

var NAVTREEINDEX =
[
"DockerUsage.html",
"Example__hydrostatic__hydrogen__burning_2scripts_2Plot__me_8py.html#aae6215e6246e894c74949ae707c2e7b3",
"class__files_2____init_____8py.html",
"classbin_1_1class__files_1_1winnet__class_1_1winnet.html#afb5370fc07b3c71215e7c088b5619e69",
"detailed__balance_8f90.html#a86bdc4751e8e1e5dfe0be6da31f94424",
"fission__rate__module_8f90.html#ac474826c9bf98253a575b32fef774bce",
"hdf5__module_8f90.html#a243c791a91c11655e7897b99d17c3e90",
"macros_8h.html#a1c1bb39c350a467f8e28fe128e8b40b7",
"namespacetabulated__rates.html",
"parameter__class_8f90.html#a80ce8d36d8505de32c4ee5baaee7cb19",
"parser__module_8f90.html#af58518a18fa64d4dcf9263533ec457ca",
"structfission__rate__module_1_1fissionrate__type.html#a9f3d7260208e8901cd9b9d6553bfc698",
"tw__rate__module_8f90.html#a73952e6c67e0868cf5fe2446878eb4ee"
];

var SYNCONMSG = 'click to disable panel synchronisation';
var SYNCOFFMSG = 'click to enable panel synchronisation';