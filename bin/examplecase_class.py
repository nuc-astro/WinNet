import os
import shutil
import sys
import numpy as np

class bcolors:
    """
      Class to deal with colors in the terminal
    """
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'



class examplecase(object):
    def __init__(self):
        """
           Initialize the examplecase class. This class prepares example runs,
           e.g., creates trajectories from other scripts and copies plotting routines
           to the example run folder.
        """
        # Define the data path of the example cases
        self.__ex_script_dir = os.path.dirname(os.path.abspath(__file__))
        self.__path_examples = os.path.join(self.__ex_script_dir,"../data/Example_data")
        self.__script_folder = "scripts/"

        self.__example_pars = ['Example_BigBang.par',\
                               'Example_BigBang_many.par',\
                               'Example_NSM_dyn_ejecta_rosswog.par',\
                               'Example_NSM_dyn_ejecta_fission_rosswog.par',\
                               'Example_NSM_wind_martin.par',\
                               'Example_NSM_wind_bins_martin.par',\
                               'Example_NSM_dyn_ejecta_ns10ns10_rosswog.par',\
                               "Example_NSM_dyn_ejecta_bovard.par",\
                               "Example_NSM_disk_wu.par",\
                               "Example_NSBH_rosswog.par",\
                               'Example_MRSN_r_process_beta_winteler.par',\
                               'Example_MRSN_r_process_winteler.par',\
                               'Example_MRSN_r_process_obergaulinger.par',\
                               'Example_MRSN_weakr_obergaulinger.par',\
                               'Example_MRSN_nup_process_obergaulinger.par',\
                               "Example_MRSN_r_process_reichert.par",\
                               "Example_CCSN_wind_bliss.par",\
                               "Example_CCSN_explosive_burning_parametrized.par",\
                               "Example_i_process_dardelet.par",\
                               "Example_hydrostatic_hydrogen_burning.par",\
                               "Example_xrayburst_schatz.par",\
                               "Example_CO_burning.par",\
                               "Example_type_Ia_meakin.par",\
                               "Example_classical_novae_jose.par",\
                               "Example_AGB_cescutti.par",\
                               "Example_AGB_nishimura.par"
                               ]


    def prepare_examples(self,parname):
        """
          Prepare the simulation before any runfolder exists.
        """
        # Check if the run is an example case
        if (not (parname in self.__example_pars)):
            return
        # Get the current working directory
        current_dir = os.getcwd()
        # Define the data path of the example case
        example_dir = os.path.join(self.__path_examples,parname.replace('.par',''))
        # Special operations for the individual examples
        # Create the trajectory by executing a script
        if (parname == 'Example_BigBang.par'):
            os.chdir(example_dir)
            os.system('python create_trajectory.py')
        # Create the trajectories by executing a script
        if (parname == 'Example_BigBang_many.par'):
            os.chdir(example_dir)
            os.system('python create_trajectories.py')
        # Neutron star merger example, download from webpage
        if (parname == 'Example_NSM_dyn_ejecta_ns10ns10_rosswog.par'):
            os.chdir(example_dir)
            if not os.path.isfile("trajectory_00001.dat"):
                try:
                    # Download the model
                    url = 'https://compact-merger.astro.su.se/data/trajectories_chk30_ns10_ns10.zip'
                    os.system("wget -q "+url)
                    zipname = os.path.basename(url)
                except:
                    print("Failed to download the example from "+url+" . Exiting.")
                    sys.exit(1)
                try:
                    # Unzip the model
                    os.system("unzip -qq "+zipname)
                    os.system("mv "+zipname.replace(".zip","")+"/* .")
                    os.system("rm -r "+zipname.replace(".zip",""))
                    os.system("rm "+zipname)
                    if not os.path.isfile("trajectory_00001.dat"):
                        raise Exception
                except:
                    print("Failed to unzip the trajectories. Exiting.")
                    sys.exit(1)
                # Ensure that they have a "#" infront for reading
                all_trajs = os.listdir(example_dir)
                for tr in all_trajs:
                    if not "trajectory_" in tr:
                        continue
                    command ="sed -i 's%t%# t%g' "+str(tr)
                    os.system(command)
        # Neutron star blackhole merger example, download from webpage
        if (parname == 'Example_NSBH_rosswog.par'):
            os.chdir(example_dir)
            if not os.path.isfile("trajectory_ns14_BH10_irrot_1_checked.dat"):
                try:
                    # Download the model
                    url = 'https://compact-merger.astro.su.se/data/trajecs_BH10_checked.zip'
                    os.system("wget -q "+url)
                    zipname = os.path.basename(url)
                except:
                    print("Failed to download the example from "+url+" . Exiting.")
                    sys.exit(1)
                try:
                    # Unzip the model
                    os.system("unzip -qq "+zipname)
                    os.system("mv "+zipname.replace(".zip","")+"/* .")
                    os.system("rm -r "+zipname.replace(".zip",""))
                    os.system("rm "+zipname)
                    if not os.path.isfile("trajectory_ns14_BH10_irrot_1_checked.dat"):
                        raise Exception
                except:
                    print("Failed to unzip the trajectories. Exiting.")
                    sys.exit(1)
                # Ensure that they have a "#" infront for reading
                all_trajs = os.listdir(example_dir)
                for tr in all_trajs:
                    if not "trajectory_" in tr:
                        continue
                    command ="sed -i 's%t%# t%g' "+str(tr)
                    os.system(command)

        # Modify reaclib
        if (parname == 'Example_i_process_dardelet.par'):
            os.chdir(example_dir)

            if not os.path.isfile("modified_reaclib"):
                reaclib_path = os.path.join(example_dir,"../../Reaclib_18_9_20")

                with open(reaclib_path,"r") as f:
                    lines = f.readlines()
                    mask  = np.ones(len(lines),dtype=bool)
                    for ind,line in enumerate(lines):
                        # Remove all occurences n13(p,gamma)o14
                        if "p  n13  o14" in line:
                            mask[ind:ind+3] = False
                        # Skip behind chapter 5
                        if "5"+50*" " in line:
                            break
                    lines_without_rate = np.array(lines)[mask]

                with open("modified_reaclib","w") as f:
                    f.write(''.join(lines_without_rate))

        if (parname == 'Example_xrayburst_schatz.par'):
            os.chdir(example_dir)

            if not os.path.isfile("sunet"):
                try:
                    # Unzip the model
                    os.system("unzip -qq example_data.zip")
                    # os.system("mv "+zipname.replace(".zip","")+"/* .")
                    # os.system("rm -r "+zipname.replace(".zip",""))
                    # os.system("rm "+zipname)
                    if not os.path.isfile("sunet"):
                        raise Exception
                except:
                    print("Failed to unzip the example data. Exiting.")
                    sys.exit(1)

        os.chdir(current_dir)


    def copy_scripts(self,parname,rundir):
        """
          Copy plotting routines to the run folders.
        """
        # Check if the run is an example case
        if (not (parname in self.__example_pars)):
            return

        # Get the current working directory
        current_dir = os.getcwd()
        # Define the data path of the example case
        example_dir = os.path.join(self.__path_examples,parname.replace('.par',''))

        # Path that contains all scripts
        s_dir       = os.path.join(example_dir,self.__script_folder)
        all_scripts = os.listdir(s_dir)

        try:
            # Copy all scripts
            for f in all_scripts:
                path = os.path.join(s_dir,f)
                shutil.copy(path, rundir)
        except:
            pass

    def __link_lit(self,string,link):
        """
          Create a colored hyperlink in the terminal
        """
        out = "\e]8;;"+link+"\a"+bcolors.OKGREEN+string+bcolors.ENDC+"\e]8;;\a"
        return out

    def __repr__(self):
        """
          Representative text for the examples
        """
        # Header
        s = "\n"
        s+= "WinNet contains example cases from hydrodynamical simulations in different astrophysical environments. "\
            "The parameter files that act as an interface between the user and the code for these cases are contained in the "\
            +bcolors.WARNING+"par/"+bcolors.ENDC+" directory. The trajectories (containing temperature, density, ...) are contained "\
            "in "+bcolors.WARNING+"data/Example_data/"+bcolors.ENDC+". Each example will contain a file called "\
            +bcolors.WARNING+"Plot_me.py "+bcolors.ENDC+"in the "+bcolors.WARNING+"runs/"+bcolors.ENDC+" folder, "\
            "which plots different outputs of the nucleosynthesis calculation. "\
            "In case you want to use trajectories for a publication, "\
            "please cite the corresponding papers that are marked in green below. In contrast to the original works, "\
            "the example cases use different reaction rates from a recent compilation of JINA (Cyburt et al. 2010) and possibly different input parameters. "\
            "Furthermore, all example cases are only representative trajectories out of a much larger set used in the publications. "\
            "The results can therefore differ quantitatively with respect to the original publication.\n\n"
        # s = bcolors.HEADER+"Example cases:\n"+bcolors.ENDC
        # s+= bcolors.HEADER+"--------------\n"+bcolors.ENDC
        s+= "\n"

        #BigBang
        s+= bcolors.UNDERLINE+"Big Bang:\n"+bcolors.ENDC
        s+= '\n'

        s+= "Big Bang nucleosynthesis. The corresponding trajectory is calculated as described in "+"Winteler 2014. "+\
            "The example calculates the element production "\
            "with a"+bcolors.OKCYAN+" photon to baryon ratio of 6.11e-10.\n"+bcolors.ENDC
        s+= bcolors.BOLD  +"Run it with:\n"+bcolors.ENDC
        s+= bcolors.OKBLUE+"python makerun.py -p Example_BigBang.par -r Example_BigBang\n"+bcolors.ENDC
        s+= "To calculate many runs with"+bcolors.OKCYAN+" different photon to baryon ratios"+bcolors.ENDC+" you can run:\n"
        s+= bcolors.OKBLUE+"python makerun.py --many -p Example_BigBang_many.par -r Example_BigBang_many\n"+bcolors.ENDC
        s+= bcolors.BOLD  +"Literature:\n"+bcolors.ENDC
        s+= "-\n"
        s+= '\n\n\n'

        # NSM
        s+= bcolors.UNDERLINE+"Neutron star merger:\n"+bcolors.ENDC
        s+= '\n'


        # Rosswog/Korobkin/Piran et al.
        s+= "Run a model with 30 trajectories. The simulation was performed using two neutron stars with 1 solar mass each. "\
            "The trajectories only included the "+bcolors.OKCYAN+"dynamic ejecta"+bcolors.ENDC+". "\
            "When running the example, the trajectories are downloaded from "\
            "https://compact-merger.astro.su.se/downloads_fluid_trajectories.html.\n"
        s+= "If you want to run only one trajectory, you can run the following example case:\n"
        s+= bcolors.BOLD  +"Run it with:\n"+bcolors.ENDC
        s+= bcolors.OKBLUE+"python makerun.py -p Example_NSM_dyn_ejecta_rosswog.par -r Example_NSM_dyn_ejecta_rosswog\n"+bcolors.ENDC
        s+= "Or in case you want to run this trajectory with three different fission fragment distributions:\n"
        s+= bcolors.OKBLUE+"python makerun.py -p Example_NSM_dyn_ejecta_fission_rosswog.par -r Example_NSM_dyn_ejecta_fission_rosswog --many --val 1,2,3\n"+bcolors.ENDC
        s+= "In case you want to run all available trajectories (30), run the following command. Be aware that this will most likely take more than one hour:\n"
        s+= bcolors.OKBLUE+"python makerun.py --many -p Example_NSM_dyn_ejecta_ns10ns10_rosswog.par -r Example_NSM_dyn_ejecta_ns10ns10_rosswog"+bcolors.ENDC+"\n"
        s+= bcolors.BOLD  +"Literature:\n"+bcolors.ENDC
        s+= self.__link_lit("- Korobkin et al. 2012","https://ui.adsabs.harvard.edu/abs/2012MNRAS.426.1940K/abstract")+"\n"
        s+= self.__link_lit("- Rosswog et al. 2013","https://ui.adsabs.harvard.edu/abs/2013MNRAS.430.2585R/abstract")+"\n"
        s+= self.__link_lit("- Piran et al. 2013","https://ui.adsabs.harvard.edu/abs/2013MNRAS.430.2121P/abstract")+"\n"
        s+= '\n\n'

        # Bovard et al.
        s+= "Run two trajectories of model LS220-M1.35. One trajectory has a low mass and one is representative for the whole model. "\
            "The trajectories calculate the "+bcolors.OKCYAN+"dynamic ejecta"+bcolors.ENDC+" of a neutron star merger.\n"
        s+= bcolors.BOLD  +"Run it with:\n"+bcolors.ENDC
        s+= bcolors.OKBLUE+"python makerun.py --many -p Example_NSM_dyn_ejecta_bovard.par -r Example_NSM_dyn_ejecta_bovard"+bcolors.ENDC
        s+= '\n'
        s+= bcolors.BOLD  +"Literature:\n"+bcolors.ENDC
        s+= self.__link_lit("- Bovard et al. 2017","https://ui.adsabs.harvard.edu/abs/2017PhRvD..96l4005B/abstract")+"\n"
        s+= '\n\n'

        # Martin et al.
        s+= "Representative trajectory (bin 4, 90ms) of the "+bcolors.OKCYAN+"neutrino driven wind"+bcolors.ENDC+" of a neutron star merger. \n"
        s+= bcolors.BOLD  +"Run it with:\n"+bcolors.ENDC
        s+= bcolors.OKBLUE+"python makerun.py -p Example_NSM_wind_martin.par -r Example_NSM_wind_martin\n"+bcolors.ENDC
        s+= "Or alternatively run all different directional bins (4 tracers). \n"
        s+= bcolors.BOLD  +"Run it with:\n"+bcolors.ENDC
        s+= bcolors.OKBLUE+"python makerun.py --many -p Example_NSM_wind_bins_martin.par -r Example_NSM_wind_bins_martin\n"+bcolors.ENDC
        s+= bcolors.BOLD  +"Literature:\n"+bcolors.ENDC
        s+= self.__link_lit("- Perego et al. 2014","https://ui.adsabs.harvard.edu/abs/2014MNRAS.443.3134P/abstract")+"\n"
        s+= self.__link_lit("- Martin et al. 2015","https://ui.adsabs.harvard.edu/abs/2015ApJ...813....2M/abstract")+"\n"
        s+= '\n\n'

        # Wu et al.
        s+= "Representative trajectories (10) of the "+bcolors.OKCYAN+"viscous/disk ejecta"+bcolors.ENDC+" (model S-def) of a neutron star merger. "\
            "The model has also been investigated in Lippuner et al. 2017.\n"
        s+= bcolors.BOLD  +"Run it with:\n"+bcolors.ENDC
        s+= bcolors.OKBLUE+"python makerun.py --many -p Example_NSM_disk_wu.par -r Example_NSM_disk_wu\n"+bcolors.ENDC
        s+= bcolors.BOLD  +"Literature:\n"+bcolors.ENDC
        s+= self.__link_lit("- Wu et al. 2016","https://ui.adsabs.harvard.edu/abs/2016MNRAS.463.2323W/abstract")+"\n"
        s+= self.__link_lit("- Lippuner et al. 2017","https://ui.adsabs.harvard.edu/abs/2017MNRAS.472..904L/abstract")+"\n"
        s+= '\n\n\n'


        # NSBH
        s+= bcolors.UNDERLINE+"Neutron star black hole merger:\n"+bcolors.ENDC
        s+= '\n'
        # Rosswog/Korobkin/Piran et al.
        s+= "Run 20 trajectories of a neutron star black hole merger. "\
            "The simulation was performed using one neutron star with 1.4 solar masses and a black-hole with 10 solar masses. "\
            "The trajectories only included the "+bcolors.OKCYAN+"dynamic ejecta"+bcolors.ENDC+". "\
            "When running the example, the trajectories are downloaded from "\
            "https://compact-merger.astro.su.se/downloads_fluid_trajectories.html. "\
            "Due to the amount of calculations, this example case may take a while.\n"
        s+= bcolors.BOLD  +"Run it with:\n"+bcolors.ENDC
        s+= bcolors.OKBLUE+"python makerun.py --many -p Example_NSBH_rosswog.par -r Example_NSBH_rosswog"+bcolors.ENDC
        s+= '\n'
        s+= bcolors.BOLD  +"Literature:\n"+bcolors.ENDC
        s+= self.__link_lit("- Korobkin et al. 2012","https://ui.adsabs.harvard.edu/abs/2012MNRAS.426.1940K/abstract")+"\n"
        s+= self.__link_lit("- Rosswog et al. 2013","https://ui.adsabs.harvard.edu/abs/2013MNRAS.430.2585R/abstract")+"\n"
        s+= self.__link_lit("- Piran et al. 2013","https://ui.adsabs.harvard.edu/abs/2013MNRAS.430.2121P/abstract")+"\n"
        s+= '\n\n\n'


        # MR-supernova
        s+= bcolors.UNDERLINE+"Magneto-rotational driven supernova:\n"+bcolors.ENDC
        s+= '\n'


        # Winteler et al.
        s+="Neutron-rich trajectory of a magneto-rotational driven supernova. Within this tracer particle, "\
           "heavy elements get synthesized via the "+bcolors.OKCYAN+"r-process"+bcolors.ENDC+".\n"
        s+= bcolors.BOLD  +"Run it with:\n"+bcolors.ENDC
        s+= bcolors.OKBLUE+"python makerun.py -p Example_MRSN_r_process_winteler.par -r Example_MRSN_r_process_winteler\n"+bcolors.ENDC
        s+= "Or in case you want to run this trajectory with three different beta-decay rates:\n"
        s+= bcolors.OKBLUE+"python makerun.py -p Example_MRSN_r_process_beta_winteler.par -r Example_MRSN_r_process_beta_winteler --many --val beta_decay_marketin.dat,beta_decay_moeller.dat,beta_decay_reaclib.dat\n"+bcolors.ENDC
        s+= bcolors.BOLD  +"Literature:\n"+bcolors.ENDC
        s+= self.__link_lit("- Winteler et al. 2012","https://ui.adsabs.harvard.edu/abs/2012ApJ...750L..22W/abstract")+"\n"
        s+= '\n\n'

        # Obergaulinger/Reichert et al.
        s+="Neutron-rich trajectory of a magneto-rotational driven supernova (35OC-Rs). Within this tracer particle, "\
           "heavy elements get synthesized via the "+bcolors.OKCYAN+"r-process"+bcolors.ENDC+".\n"
        s+= bcolors.BOLD  +"Run it with:\n"+bcolors.ENDC
        s+= bcolors.OKBLUE+"python makerun.py -p Example_MRSN_r_process_obergaulinger.par -r Example_MRSN_r_process_obergaulinger\n"+bcolors.ENDC
        s+="To run a less neutron rich trajectory that synthesizes elements by the"+bcolors.OKCYAN+" weak-r process "+bcolors.ENDC+\
            "from the same model run:\n"
        s+= bcolors.OKBLUE+"python makerun.py -p Example_MRSN_weakr_obergaulinger.par -r Example_MRSN_weakr_obergaulinger\n"+bcolors.ENDC
        s+="A proton-rich trajectory that synthesizes elements via the"+bcolors.OKCYAN+" nu-p process "+bcolors.ENDC+\
            "from model 35OC-RO can be run by:\n"
        s+= bcolors.OKBLUE+"python makerun.py -p Example_MRSN_nup_process_obergaulinger.par -r Example_MRSN_nup_process_obergaulinger\n"+bcolors.ENDC
        s+= bcolors.BOLD  +"Literature:\n"+bcolors.ENDC
        s+= self.__link_lit("- Obergaulinger & Aloy 2017","https://ui.adsabs.harvard.edu/abs/2017MNRAS.469L..43O/abstract")+"\n"
        s+= self.__link_lit("- Reichert et al. 2021","https://ui.adsabs.harvard.edu/abs/2021MNRAS.501.5733R/abstract")+"\n"
        s+= '\n\n'

        # Aloy, Obergaulinger, Reichert
        s+="High entropy trajectory as well as two neutron-rich tracers of a magneto-rotational driven supernova. "\
           "The high entropy tracer originates in model P, the neutron-rich tracers from model 35OC-Rs_N. "\
           "In all cases heavy elements get synthesized via the "+bcolors.OKCYAN+"r-process"+bcolors.ENDC+".\n"
        s+= bcolors.BOLD  +"Run it with:\n"+bcolors.ENDC
        s+= bcolors.OKBLUE+"python makerun.py --many -p Example_MRSN_r_process_reichert.par -r Example_MRSN_r_process_reichert\n"+bcolors.ENDC
        s+= bcolors.BOLD  +"Literature:\n"+bcolors.ENDC
        s+= self.__link_lit("- Aloy & Obergaulinger 2021","https://ui.adsabs.harvard.edu/abs/2021MNRAS.500.4365A/abstract")+"\n"
        s+= self.__link_lit("- Obergaulinger & Aloy 2021","https://ui.adsabs.harvard.edu/abs/2021MNRAS.503.4942O/abstract")+"\n"
        s+= self.__link_lit("- Reichert et al. 2023","https://ui.adsabs.harvard.edu/abs/2023MNRAS.518.1557R/abstract")+"\n"
        s+= '\n\n\n'

        s+= bcolors.UNDERLINE+"Classical Novae:\n"+bcolors.ENDC
        s+= '\n'

        # Jose & Hernandez
        s+="Trajectory of a classical novae (model ONe5)."\
           "The trajectory can be downloaded at https://zenodo.org/record/6474694 (v 1.1.1). Also the following description is given there: "\
           "The model follows 1 outburst from the initiation of the accretion stage and through the explosion, expansion and ejection. "\
           "It relies on a 1.25 Msun, ONe WD hosting the explosion. For simplicity, it is assumed pre-enriched, accreted material with "\
           "a composition corresponding to 50% solar (from Lodders 2009) and 50% outer layers of the WD substrate (from Ritossa et al. 1996). "\
           "The initial luminosity of the WD is assumed to be 10^-2 Lsun, and the mass-accretion rate is 2x10^-10 Msun yr^-1.\n"
        s+= bcolors.BOLD  +"Run it with:\n"+bcolors.ENDC
        s+= bcolors.OKBLUE+"python makerun.py -p Example_classical_novae_jose.par -r Example_classical_novae_jose\n"+bcolors.ENDC
        s+= bcolors.BOLD  +"Literature:\n"+bcolors.ENDC
        s+= self.__link_lit("- Jose & Hernanz 1998","https://ui.adsabs.harvard.edu/abs/1998ApJ...494..680J/abstract")+"\n"
        s+= self.__link_lit("- Jose 2022","https://doi.org/10.5281/zenodo.6474694")+"\n"
        s+= '\n\n\n'

        s+= bcolors.UNDERLINE+"Accreting neutron stars:\n"+bcolors.ENDC
        s+= '\n'

        # Schatz et al.
        s+="Trajectory of a X-ray burst. The trajectory originates from a simulation of an accreting "\
        "neutron star. The nucleosynthesis is dominated by the "+bcolors.OKCYAN+"rp-process"+bcolors.ENDC+".\n"
        s+= bcolors.BOLD  +"Run it with:\n"+bcolors.ENDC
        s+= bcolors.OKBLUE+"python makerun.py -p Example_xrayburst_schatz.par -r Example_xrayburst_schatz\n"+bcolors.ENDC
        s+= bcolors.BOLD  +"Literature:\n"+bcolors.ENDC
        s+= self.__link_lit("- Schatz et al. 2002","https://ui.adsabs.harvard.edu/abs/2001NuPhA.688..150S/abstract")+"\n"
        s+= '\n\n\n'

        s+= bcolors.UNDERLINE+"Regular Core-collapse supernovae:\n"+bcolors.ENDC
        s+= '\n'

        # Parametrized
        s+="Simple parametric model to calculate "+bcolors.OKCYAN+"complete Si burning"+bcolors.ENDC+"."+"\n"\
           "The trajectory starts out of NSE and no initial composition is required.\n"
        s+= bcolors.BOLD  +"Run it with:\n"+bcolors.ENDC
        s+= bcolors.OKBLUE+"python makerun.py -p Example_CCSN_explosive_burning_parametrized.par -r Example_CCSN_explosive_burning_parametrized\n"+bcolors.ENDC
        s+= bcolors.BOLD  +"Literature:\n"+bcolors.ENDC
        s+= "-\n"
        s+= '\n\n'

        # Bliss et al.
        s+="Trajectory assuming a steady state model of a "+bcolors.OKCYAN+"neutrino driven wind."+bcolors.ENDC+"\n"\
           "The trajectory (CPR2) is publicly available at https://theorie.ikp.physik.tu-darmstadt.de/astro/resources.php.\n"
        s+= bcolors.BOLD  +"Run it with:\n"+bcolors.ENDC
        s+= bcolors.OKBLUE+"python makerun.py -p Example_CCSN_wind_bliss.par -r Example_CCSN_wind_bliss\n"+bcolors.ENDC
        s+= bcolors.BOLD  +"Literature:\n"+bcolors.ENDC
        s+= self.__link_lit("- Bliss et al. 2018","https://ui.adsabs.harvard.edu/abs/2018ApJ...855..135B/abstract")+"\n"
        s+= '\n\n\n'



        s+= bcolors.UNDERLINE+"Type Ia Supernova:\n"+bcolors.ENDC
        s+= '\n'
        # Bliss et al.
        s+="Parametrization of the detonation phase of a Type Ia Supernova."+" The trajectory undergoes"+bcolors.OKCYAN+" explosive burning"+\
           bcolors.ENDC+".\n"
        s+= bcolors.BOLD  +"Run it with:\n"+bcolors.ENDC
        s+= bcolors.OKBLUE+"python makerun.py -p Example_type_Ia_meakin.par -r Example_type_Ia_meakin\n"+bcolors.ENDC
        s+= bcolors.BOLD  +"Literature:\n"+bcolors.ENDC
        s+= self.__link_lit("- Meakin et al. 2009","https://ui.adsabs.harvard.edu/abs/2009ApJ...693.1188M/abstract")+"\n"
        s+= '\n\n\n'


        s+= bcolors.UNDERLINE+"AGB stars:\n"+bcolors.ENDC
        s+= '\n'

        s+="Carbon 13 and thermal pulse trajectories of a 3 Msun, Z = 0.014 (solar metallicity) star. "+\
           "Elements get synthesized within a "+bcolors.OKCYAN+"strong s-process"+bcolors.ENDC+". "
        s+="The two trajectories were calculated with the 1D stellar evolution code MESA. They were "+\
           "accessed at https://zenodo.org/record/6474686 (v.1.2.1). \n"
        s+= bcolors.BOLD  +"Run it with:\n"+bcolors.ENDC
        s+= bcolors.OKBLUE+"python makerun.py --many -p Example_AGB_cescutti.par -r Example_AGB_cescutti\n"+bcolors.ENDC
        s+= bcolors.BOLD  +"Literature:\n"+bcolors.ENDC
        s+= self.__link_lit("- Cescutti et al. 2018","https://ui.adsabs.harvard.edu/abs/2018MNRAS.478.4101C/abstract")+"\n"
        s+= self.__link_lit("- Cescutti 2022","https://doi.org/10.5281/zenodo.6474686")+"\n"
        s+= '\n\n'

        s+="Representative trajectory of a "+bcolors.OKCYAN+"weak s-process"+bcolors.ENDC+". "+\
           "The trajectory was extracted from a 25 Msun, Z = 0.014 (solar metallicity) stellar evolution model. "
        s+="It was chosen because it roughly corresponds to the average weak s-process production in massive stars "+\
           "weighted over the initial mass function. This information as well as the trajectory itself were "+\
           "accessed at https://zenodo.org/record/6474728 (v 1.1.1). \n"
        s+= bcolors.BOLD  +"Run it with:\n"+bcolors.ENDC
        s+= bcolors.OKBLUE+"python makerun.py -p Example_AGB_nishimura.par -r Example_AGB_nishimura\n"+bcolors.ENDC
        s+= bcolors.BOLD  +"Literature:\n"+bcolors.ENDC
        s+= self.__link_lit("- Hirschi et al. 2004","https://ui.adsabs.harvard.edu/abs/2004A%26A...425..649H/abstract")+"\n"
        s+= self.__link_lit("- Nishimura et al. 2017","https://ui.adsabs.harvard.edu/abs/2017MNRAS.469.1752N/abstract")+"\n"
        s+= self.__link_lit("- Pignatari & Hirschi 2022","https://zenodo.org/record/6474728")+"\n"
        s+= '\n\n\n'


        s+= bcolors.UNDERLINE+"Hydrostatic burning:\n"+bcolors.ENDC
        s+= '\n'

        # Hydrogen Burning
        s+="Hydrostatic "+bcolors.OKCYAN+"hydrogen burning"+bcolors.ENDC+" for a constant density of 100 g/ccm and "
        s+="different temperatures (between 1e7 K and 4e7 K). To be able to plot the output, the hdf5 output has to "
        s+="be enabled and configured in the Makefile.\n"
        s+= bcolors.BOLD  +"Run it with:\n"+bcolors.ENDC
        s+= bcolors.OKBLUE+"python makerun.py -p Example_hydrostatic_hydrogen_burning.par -r Example_hydrostatic_hydrogen_burning --many --val_min=1 --val_max=4 --val_it=0.1\n"+bcolors.ENDC
        s+= bcolors.BOLD  +"Literature:\n"+bcolors.ENDC
        s+= "-\n"
        s+= '\n\n'

        # Carbon-oxygen Burning
        s+="Hydrostatic "+bcolors.OKCYAN+"carbon-oxygen burning"+bcolors.ENDC+" for a constant density of 1e9 g/ccm and "
        s+="a temperature of 3 GK.\n"
        s+= bcolors.BOLD  +"Run it with:\n"+bcolors.ENDC
        s+= bcolors.OKBLUE+"python makerun.py -p Example_CO_burning.par -r Example_CO_burning\n"+bcolors.ENDC
        s+= bcolors.BOLD  +"Literature:\n"+bcolors.ENDC
        s+= "-\n"
        s+= '\n\n'

        # Dardelet et al.
        s+="A simple model, simulating an "+bcolors.OKCYAN+"i-process"+bcolors.ENDC+" within a one-zone nuclear reaction network.\n"
        s+= bcolors.BOLD  +"Run it with:\n"+bcolors.ENDC
        s+= bcolors.OKBLUE+"python makerun.py --many -p Example_i_process_dardelet.par -r Example_i_process_dardelet\n"+bcolors.ENDC
        s+= bcolors.BOLD  +"Literature:\n"+bcolors.ENDC
        s+= self.__link_lit("- Dardelet et al. 2015","https://ui.adsabs.harvard.edu/abs/2015arXiv150505500D/abstract")+"\n"
        s+= '\n\n'


        s = "'"+s+"'"
        return s







if __name__ == '__main__':
    # Test the class
    example = examplecase()
    # example.prepare_examples('Example_BigBang_many.par')
    # example.copy_scripts('Example_BigBang_many.par','test')
    print(example)
