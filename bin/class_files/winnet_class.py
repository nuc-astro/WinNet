# @Author: Moritz Reichert <mreichert>
# @Date:   2017-08-04T23:19:25+02:00
# @Last modified by:   mreichert
# @Last modified time: 2017-09-17T13:12:16+02:00



# Winnet class to analyze winnet runs
# Author: Moritz Reichert
# Date  : 08.05.17

from random                     import randint
from shutil                     import copyfile
from matplotlib                 import gridspec
from matplotlib.patches         import Rectangle, Arrow, Circle,Wedge,PathPatch
from matplotlib                 import colors,cm,colorbar
from matplotlib.textpath        import TextPath
from matplotlib.font_manager    import FontProperties
from matplotlib.collections     import PatchCollection,LineCollection
from matplotlib.offsetbox       import AnchoredOffsetbox, TextArea
from mpl_toolkits.axes_grid1    import make_axes_locatable
from mpl_toolkits.mplot3d       import Axes3D
from matplotlib                 import cm
from types                      import FunctionType
from scipy.interpolate          import interp1d
from .nucleus_class             import nucleus
import matplotlib               as mpl
import matplotlib.animation     as animation
import numpy                    as np
import pandas                   as pd
import matplotlib.pyplot        as plt
import h5py
import os
import inspect
import sys
import scipy.interpolate
import math
import multiprocessing


class winnet(object):

    def __init__(self,runpath):
        """
        Class for analyzing Winnet runs
        """

        # Make sure that a directory path is given and save the path into self.__path
        if runpath[-1]!= '/':
            runpath += '/'
        self.__path = runpath

        #Database for element names
        self.__elname      = ('neutron','h','he','li','be','b','c','n','o','f','ne','na','mg','al','si','p','s','cl','ar','k','ca','sc','ti','v','cr','mn','fe',
          'co','ni','cu','zn','ga','ge','as','se','br','kr','rb','sr','y','zr','nb','mo','tc','ru','rh','pd','ag','cd','in','sn','sb',
          'te', 'i','xe','cs','ba','la','ce','pr','nd','pm','sm','eu','gd','tb','dy','ho','er','tm','yb','lu','hf','ta','w','re','os',
          'ir','pt','au','hg','tl','pb','bi','po','at','rn','fr','ra','ac','th','pa','u','np','pu','am','cm','bk','cf','es','fm','md',
          'no','lr','rf','db','sg','bh','hs','mt','ds','rg','ub','ut','uq','up','uh','us','uo')
        # Note: ub is also known as cn

        #Storage for snapshot data
        self.__time             = []
        self.__allmassfractions = []
        self.__neutrons         = []
        self.__protons          = []


        #Storage for Finab data
        self.__finab_nuclei     = []
        self.__finab_hydrogen   = None
        self.__finab_mass       = []
        self.__finab_protons    = []
        self.__finab_neutrons   = []
        self.__finab_ab         = []



        #Storage for timescale data
        self.__timescale_time     = []
        self.__timescale_temp     = []
        self.__timescale_tau_ga   = []
        self.__timescale_tau_ag   = []
        self.__timescale_tau_ng   = []
        self.__timescale_tau_gn   = []
        self.__timescale_tau_pg   = []
        self.__timescale_tau_gp   = []
        self.__timescale_tau_np   = []
        self.__timescale_tau_pn   = []
        self.__timescale_tau_an   = []
        self.__timescale_tau_na   = []
        self.__timescale_tau_beta = []


        #Storage for mainout data
        self.__mainout_iteration    = []
        self.__mainout_time         = []
        self.__mainout_temp         = []
        self.__mainout_dens         = []
        self.__mainout_ye           = []
        self.__mainout_rad          = []
        self.__mainout_yn           = []
        self.__mainout_yp           = []
        self.__mainout_ya           = []
        self.__mainout_ylight       = []
        self.__mainout_yheavy       = []
        self.__mainout_zbar         = []
        self.__mainout_abar         = []
        self.__mainout_entropy      = []


        #Storage for tracked nuclei
        self.__track_time           = []
        self.__track_Y              = []  #[[y1],[y2],..]
        self.__track_X              = []  #[[y1],[y2],..]
        self.__track_nuclei_names   = []

        # Seed nuclei
        self.__seed_nuclei          = []

        #Storage for param.out
        self.__param_parameter = []
        self.__param_data      = []

        #Default hdf5name
        self.__hdf5_name = "WinNet_data.h5"


    def get_man(self):
        """
          Returns a string with a manual
        """
        man = '1. Create a Winnet instance. E.g. Win = winnet("path to winnet folder") \n'+\
              '2. Load everything. For example: Win.read_run(), you can also specify what you want to read with Win.read_snapshots() or Win.read_finab()\n' +\
              '3. Analyze your data. Win.get_Z_X() returns the atomic number and the mass fraction of the run'
        return man

    def get_methods(self):
        """
        Returns a list of all methods of the class
        """
        return [method for method in dir(self) if callable(getattr(self, method))]


    def __read_flow(self,filename):
        """
          Reads the flow of a given file (only one flow file)
          Reading all files would result in a too high memory consumption.

          Originally from Christian Winteler and Urs Frischknecht.
          Edited by Moritz Reichert (2.4.2017)

          INPUT: flow.dat file name including path
          OUTPUT: dictionary,
          for example ff['time'] will contain the time

          example calls:
          flux=fdic.ff('../m25z01/LOGS/flow.dat')
          then flux['time'] contains the time where the fluxes were calculated
        """
        if os.path.isfile(filename):
            f = open(filename)
            lines = f.readlines()[:4]
            f.close()

            # read the stored variable names and store them in "names"
            nvars = len(lines[0].split()+lines[2].split())

            fout={'names':lines[0].split()+lines[2].split()}
            fout[fout['names'][0]] = float(lines[1].split()[0])
            fout[fout['names'][1]] = float(lines[1].split()[1])
            fout[fout['names'][2]] = float(lines[1].split()[2])

            # read data columns
            nlines = len(lines)-3
            ncols  = len(lines[2].split())
            b= np.loadtxt(filename,skiprows=3,unpack=True)

            for i in range(ncols):
                fout[lines[2].split()[i]]= b[i]

        elif os.path.isfile(os.path.join(self.__path,self.__hdf5_name)):
            f=h5py.File(os.path.join(self.__path,self.__hdf5_name),"r")
            if "/flows" in f:
                ######################################################################################################
                num = ""
                for s in filename:
                    if s.isdigit():
                        num+=s
                number = int(num)
                print(number,filename)
                if "/flows/"+str(number) in f:
                    entry =  "/flows/"+str(number)
                    fout={'nin' :f[entry+"/n_in"][:],
                          'nout':f[entry+"/n_out"][:],
                          'zin' :f[entry+"/p_in"][:],
                          'zout':f[entry+"/p_out"][:],
                          'yin' :f[entry+"/y_in"][:],
                          'yout':f[entry+"/y_out"][:],
                          'flow':f[entry+"/flow"][:],
                          'time':f[entry+"/time"][:],
                          'temp':f[entry+"/temp"][:],
                          'dens':f[entry+"/dens"][:],
                         }
                else:
                    raise Exception("Flow out of range ("+str(number)+").")

            else:
                raise Exception("No flows in hdf5 file.")
            f.close()
        else:
            raise Exception("Could not open file "+str(filename))
        return fout





    def read_run(self):
        """
          Read a whole run with snapshots and so on.
        """
        self.read_snapshots()
        self.read_finab()
        self.read_timescales()
        self.read_mainout()
        self.read_template()

    def read_template(self):
        """
          Read the templatefile or the run (param.out)
        """
        path = os.listdir(self.__path)
        for p in path:
            if p.endswith(".par"):
                path_param = p
                break
        path = os.path.join(self.__path,path_param)

        #Read template file param.out
        self.__param_parameter,self.__param_data = np.loadtxt(path,unpack=True,dtype=str,usecols=[0,1],delimiter='=')

        # Remove blanks
        self.__param_parameter = [x.strip() for x in self.__param_parameter]


    def read_seeds(self):
        """
          Read the seed composition of the run. Stores an array of nuclei class. If the seed file was not found, it stores the value np.nan
        """
        try:
             # Read the template to get path to the seed file
             self.read_template()
             name_of_seed = self.__param_data[list(self.__param_parameter).index('seed_file')].strip().replace('"','')

             path = self.__path+name_of_seed
             # Absolute or relative path?
             if os.path.isfile(path): #Relative path
                 path = path
             elif os.path.isfile(name_of_seed): #Absolute path
                 path = name_of_seed
             else: #Path not found
                 self.__seed_nuclei = np.nan
                 return

             isotopes,massfractions = np.loadtxt(path,unpack=True,dtype=str,usecols=[0,1])

             # Cast the names to nuclei
             list_nucs = [nucleus(name=x) for x in isotopes]

             # Set the abundances
             for i in range(len(list_nucs)):
                 list_nucs[i].set_Y(float(massfractions[i])/float(list_nucs[i].get_A()))
             self.__seed_nuclei = list_nucs

        except:
            print("Error when reading seeds!")




    def get_seed_nuclei(self):
        """
          Returns the seed nuclei as a list of the instance nuclei
        """
        return self.__seed_nuclei

    def __get_at_temp(self,temp,y_array):
        """
          Returns a given quantity at a given temperature.
          If the trajectory heats up several times between the temperature
          the last one is taken.
        """
        try:
            # to make the code shorter
            mt = self.__mainout_temp
            my = y_array
            # First get the index
            for ind,t in enumerate(mt[::-1]):
                if t>temp:
                    found_ind = len(mt)-ind-1
                    break
            # Interpolate linear
            out = np.interp(temp,[mt[found_ind],mt[found_ind+1]],[my[found_ind],my[found_ind+1]])
            return out


        except:
            print("Read mainout first!")
            sys.stop(1)


    def get_ye_at_temp(self,temp):
        """
          Returns the electron fraction at a given temperature.
          If the trajectory heats up several times between the temperature
          the last one is taken
        """
        return self.__get_at_temp(temp,self.__mainout_ye)


    def get_entr_at_temp(self,temp):
        """
          Returns the entropy at a given temperature.
          If the trajectory heats up several times between the temperature
          the last one is taken
        """
        return self.__get_at_temp(temp,self.__mainout_entropy)


    def get_density_at_temp(self,temp):
        """
          Returns the entropy at a given temperature.
          If the trajectory heats up several times between the temperature
          the last one is taken
        """
        return self.__get_at_temp(temp,self.__mainout_density)


    def get_yn_at_temp(self,temp):
        """
          Returns the entropy at a given temperature.
          If the trajectory heats up several times between the temperature
          the last one is taken
        """
        return self.__get_at_temp(temp,self.__mainout_yn)


    def read_snapshots(self):
        """
          Reads the snapshots from @runpath@/snaps/
        """
        path_hdf5    = os.path.join(self.__path,self.__hdf5_name)
        path         = self.__path+'snaps/'
        amount_snaps = len(os.listdir(path))

        # Read it from the ascii files
        if amount_snaps != 0:
            for i in range(amount_snaps):
                currentfile = 'snapsh_'+str((i+1)).zfill(4)+'.dat'

                #Read time
                count = 0
                with open(path+currentfile) as f:

                    for line in f.readlines():
                        count += 1
                        if count == 2:
                            self.__time.append(float(line.split()[0]))
                        if count > 2:
                            break

                #Read mass fractions
                self.__neutrons,self.__protons,massfractions = np.loadtxt(path+currentfile,unpack=True,skiprows=3,usecols=[0,1,3])

                self.__allmassfractions.append(massfractions)

        elif os.path.isfile(path_hdf5):
            f = h5py.File(path_hdf5,"r")
            if "/snapshots" in f:
                self.__time              = f["/snapshots/time"][:]
                self.__neutrons          = f["/snapshots/N"][:]
                self.__protons           = f["/snapshots/Z"][:]
                self.__allmassfractions  = f["/snapshots/Y"][:,:]*f["/snapshots/A"][:]
            else:
                raise Exception("There are no snapshots to read!")
            f.close()
        else:
            raise Exception("There are no snapshots to read!")


    def __read_snapshot_time(self):
        """
          Read only the time from the snapshot files
        """
        path         = self.__path+'snaps/'
        amount_snaps = len(os.listdir(path))

        for i in range(amount_snaps):
            currentfile = 'snapsh_'+str((i+1)).zfill(4)+'.dat'

            #Read time
            count = 0
            with open(path+currentfile) as f:

                for line in f.readlines():
                    count += 1
                    if count == 2:
                        self.__time.append(float(line.split()[0]))
                    if count > 2:
                        break
        self.__time = np.array(self.__time)


    def read_finab(self,finab_name="finab.dat"):
        """
          Reads the file finab.dat in the run folder
        """

        # Read the file (File format : 1:A 2:Z 3:N 4:Yi 5:Xi)
        self.__finab_mass, self.__finab_protons, self.__finab_neutrons, self.__finab_ab = np.loadtxt(os.path.join(self.__path,finab_name),usecols=[0,1,2,3],unpack=True)

        #Create instance of nuclei
        self.__finab_nuclei = [nucleus(Z=int(self.__finab_protons[i]),N=int(self.__finab_neutrons[i]),Y=self.__finab_ab[i]) for i in range(len(self.__finab_mass)) ]
        self.__finab_nuclei.sort()

        #Get hydrogen out of it
        for nuc in self.__finab_nuclei:
            if (nuc.get_Z() == 1 and nuc.get_N() == 0):
                self.__finab_hydrogen = nuc
                break


    def read_timescales(self,subtract_first=True):
        """
          Read the timescales from a Winnet calculation
        """

        path = os.path.join(self.__path, 'timescales.dat')
        path_hdf5 = os.path.join(self.__path, self.__hdf5_name)

        if os.path.isfile(path):
            self.__timescale_time, self.__timescale_temp, self.__timescale_tau_ga ,self.__timescale_tau_ag, self.__timescale_tau_ng , self.__timescale_tau_gn, self.__timescale_tau_pg, self.__timescale_tau_gp, self.__timescale_tau_np, self.__timescale_tau_pn, self.__timescale_tau_an, \
            self.__timescale_tau_na, self.__timescale_tau_ap, self.__timescale_tau_pa, self.__timescale_tau_beta,self.__timescale_tau_alpha,self.__timescale_tau_nfiss,self.__timescale_tau_sfiss,self.__timescale_tau_bfiss = np.loadtxt(path,unpack=True,skiprows=1)
            # self.__timescale_time, self.__timescale_temp, self.__timescale_tau_ng , self.__timescale_tau_gn, self.__timescale_tau_pg, self.__timescale_tau_gp, self.__timescale_tau_np, self.__timescale_tau_pn, self.__timescale_tau_an, self.__timescale_tau_na, self.__timescale_tau_beta = np.loadtxt(path,unpack=True,skiprows=1)
            if subtract_first:
                self.__timescale_time = self.__timescale_time-self.__timescale_time[0]
        elif os.path.isfile(path_hdf5):
            f = h5py.File(path_hdf5,"r")
            if "/timescales" in f:
                self.__timescale_time      = f["/timescales/time"][:]
                self.__timescale_tau_ng    = f["/timescales/tau_ng"][:]
                self.__timescale_tau_gn    = f["/timescales/tau_gn"][:]
                self.__timescale_tau_pg    = f["/timescales/tau_pg"][:]
                self.__timescale_tau_gp    = f["/timescales/tau_gp"][:]
                self.__timescale_tau_np    = f["/timescales/tau_np"][:]
                self.__timescale_tau_pn    = f["/timescales/tau_pn"][:]
                self.__timescale_tau_an    = f["/timescales/tau_an"][:]
                self.__timescale_tau_na    = f["/timescales/tau_na"][:]
                self.__timescale_tau_ga    = f["/timescales/tau_ga"][:]
                self.__timescale_tau_ag    = f["/timescales/tau_ag"][:]
                self.__timescale_tau_ap    = f["/timescales/tau_ap"][:]
                self.__timescale_tau_pa    = f["/timescales/tau_pa"][:]
                self.__timescale_tau_beta  = f["/timescales/tau_beta"][:]
                self.__timescale_tau_alpha = f["/timescales/tau_alpha"][:]
                self.__timescale_tau_nfiss = f["/timescales/tau_nfiss"][:]
                self.__timescale_tau_sfiss = f["/timescales/tau_sfiss"][:]
                self.__timescale_tau_bfiss = f["/timescales/tau_bfiss"][:]
                if subtract_first:
                    self.__timescale_time = self.__timescale_time-self.__timescale_time[0]
            else:
                raise Exception('Timescales were not calculated in this run.')
        else:
            raise Exception('Timescales were not calculated in this run.')

        return


    def read_mainout_fix(self):
        """
          Read the mainout with numbers like 1.8-392
        """
        path = self.__path + 'mainout.dat'


        f = open(path)
        line_number = 0
        for line in f.readlines():

            line_number +=1

            if line_number<=4:
                continue

            string_split    = line.split()
            split_corrected = []

            for i in string_split:
                try:
                    c_float = float(i)
                except:
                    c_float = 0.
                split_corrected.append(c_float)

            self.__mainout_iteration.append(split_corrected[0 ])
            self.__mainout_time     .append(split_corrected[1 ])
            self.__mainout_temp     .append(split_corrected[2 ])
            self.__mainout_dens     .append(split_corrected[3 ])
            self.__mainout_ye       .append(split_corrected[4 ])
            self.__mainout_rad      .append(split_corrected[5 ])
            self.__mainout_yn       .append(split_corrected[6 ])
            self.__mainout_yp       .append(split_corrected[7 ])
            self.__mainout_ya       .append(split_corrected[8 ])
            self.__mainout_ylight   .append(split_corrected[9 ])
            self.__mainout_yheavy   .append(split_corrected[10])
            self.__mainout_zbar     .append(split_corrected[11])
            self.__mainout_abar     .append(split_corrected[12])
            self.__mainout_entropy  .append(split_corrected[13])


    def read_tracked_nuclei(self):
        """
        Read the nuclei data contained in the file tracked_nuclei. This file is helpfull, when you dont want to look at the full output.
        """
        path = self.__path + 'tracked_nuclei.dat'
        # Check that the file exists
        if (os.path.isfile(path)):
            # Read the first line of the file and extract the nuclei names from the header
            with open(path, 'r') as f:
                first_line = f.readline()

            first_line = first_line.replace('Y(','')
            first_line = first_line.replace(')','')
            first_line = first_line.replace('#','')
            self.__track_nuclei_names = first_line.split()[1:]      # Since the first entry is "time"

            # Read the time and the abundances of the nuclei
            alldata = np.loadtxt(path,unpack=True)
            self.__track_time = alldata[0]
            self.__track_Y    = alldata[1:]

            # Save also the mass fraction
            nuclei_helper_instance = np.array([nucleus(x).get_A() for x in self.__track_nuclei_names])
            self.__track_X         = np.array(list(map(lambda x,y:  x*y ,self.__track_Y,nuclei_helper_instance)))
        else:
            print('No nuclei are tracked!')




    def read_mainout(self,subtract_first=False):
        """
          Read the mainout from a Winnet calculation
        """
        path = self.__path + 'mainout.dat'
        self.__mainout_iteration, self.__mainout_time, self.__mainout_temp, self.__mainout_dens, self.__mainout_ye, self.__mainout_rad, self.__mainout_yn, self.__mainout_yp, self.__mainout_ya, self.__mainout_ylight, self.__mainout_yheavy, self.__mainout_zbar, self.__mainout_abar, self.__mainout_entropy = np.loadtxt(path,skiprows = 3,unpack=True,usecols=np.arange(14),ndmin=1)

        if isinstance(self.__mainout_iteration,np.float64):
            (self.__mainout_iteration, self.__mainout_time, self.__mainout_temp, self.__mainout_dens, self.__mainout_ye, self.__mainout_rad, self.__mainout_yn, self.__mainout_yp, self.__mainout_ya, self.__mainout_ylight, self.__mainout_yheavy, self.__mainout_zbar, self.__mainout_abar, self.__mainout_entropy) = ([self.__mainout_iteration], [self.__mainout_time], [self.__mainout_temp], [self.__mainout_dens], [self.__mainout_ye], [self.__mainout_rad], [self.__mainout_yn], [self.__mainout_yp], [self.__mainout_ya], [self.__mainout_ylight], [self.__mainout_yheavy], [self.__mainout_zbar], [self.__mainout_abar], [self.__mainout_entropy])
        if subtract_first:
            self.__mainout_time = self.__mainout_time-self.__mainout_time[0]


    def read_engen(self):
        """
          Read the generated energy
        """
        path = self.__path + 'generated_energy.dat'
        # Check that the file exists
        if (os.path.isfile(path)):
            with open(path, 'r') as f:
                first_line = f.readline()
                second_line = f.readline()

            self.__engen_keywords = second_line.split()[1:]
            self.__engen = np.loadtxt(path,unpack=True)
        else:
            print('No generated energy was calculated!')

    def get_engen(self):
        """
          Return the generated energy
        """
        return self.__engen

    def get_tracked_nuclei_names(self):
        return self.__track_nuclei_names

    def get_tracked_nuclei(self,name,massfraction=True):
        """
          Get the abundance of a nuclei that was tracked.
          Input:
            name         - Name of the nucleus (e.g. ti44)
            massfraction - Return mass fractions or abundances?
        """
        # Get the lowercase of the input
        name  = name.lower()
        try:
            # Index of the names and track_Y is the same, so get it and return the values
            index = self.__track_nuclei_names.index(name)
            if massfraction:

                return self.__track_X[index]
            else:
                return self.__track_Y[index]
        except:
            print('"'+name+'" not contained in tracked nuclei.')
            return None


    def get_tracked_time(self):
        """
          Get the time, contained in the file tracked_nuclei.dat
        """
        return self.__track_time


    def get_mainout_time(self):
        """
          Get the time, contained in mainout.dat
        """
        return self.__mainout_time


    def get_mainout_yn(self):
        """
          Get the neutron abundance, contained in mainout.dat
        """
        return self.__mainout_yn

    def get_mainout_yp(self):
        """
          Get the proton abundance, contained in mainout.dat
        """
        return self.__mainout_yp

    def get_mainout_ya(self):
        """
          Get the alpha abundance, contained in mainout.dat
        """
        return self.__mainout_ya

    def get_mainout_yheavy(self):
        """
          Get the yheavy abundance, contained in mainout.dat
        """
        return self.__mainout_yheavy

    def get_is_crashed(self):
        """
          Is the run crashed?
        """
        if not os.path.exists(self.__path+'finab.dat'):
            return True
        else:
            return False


    def get_trajectory(self):
        """
          Get all data from the trajectory
          Output is:
            time[s], radius[cm], dens[g/ccm], temp[GK], Ye, entr[kb/nucleon]
        """

        if len(self.__param_parameter) == 0:
            print('Call read_template first! Your script might not work now!')
            return

        # Get path to trajectory file
        ind = self.__param_parameter.index('trajectory_file')
        traj_path = self.__param_data[ind].replace('"','').strip()
        time, radius, dens, temp, Ye, entr = np.loadtxt(traj_path, unpack=True, skiprows=2)
        return time,radius,dens,temp,Ye,entr

    def get_trajectory_path(self):
        """
          Get all data from the trajectory
          Output is the path to the trajectory
        """

        if len(self.__param_parameter) == 0:
            print('Call read_template first! Your script might not work now!')
            return

        # Get path to trajectory file
        ind = self.__param_parameter.index('trajectory_file')
        traj_path = self.__param_data[ind].replace('"','').strip()
        return traj_path


    def get_mainout(self):
        """
          Returns the data from mainout.dat
          Output is:
          iteration, time, temperature, density, Ye, Radius, Yn, Yp, Y_alpha, Y_light, Y_heavy, Zbar, Abar, Entropy
        """
        return self.__mainout_iteration, self.__mainout_time, self.__mainout_temp, self.__mainout_dens, self.__mainout_ye, self.__mainout_rad, self.__mainout_yn, self.__mainout_yp, self.__mainout_ya, self.__mainout_ylight, self.__mainout_yheavy, self.__mainout_zbar, self.__mainout_abar, self.__mainout_entropy


    def get_temperature(self):
        """
        Get the temperature [GK] from mainout.dat, it also returns the time
        """
        if len(self.__mainout_time) != 0:
            return self.__mainout_time, self.__mainout_temp
        else:
            print('Call read_mainout first! Your script might not work now')


    def get_density(self):
        """
        Get the density [g/ccm] from mainout.dat, it also returns the time
        """
        if len(self.__mainout_time) != 0:
            return self.__mainout_time, self.__mainout_dens
        else:
            print('Call read_mainout first! Your script might not work now')


    def get_ye(self):
        """
        Get the electron fraction from mainout.dat, it also returns the time
        """
        if len(self.__mainout_time) != 0:
            return self.__mainout_time, self.__mainout_ye
        else:
            print('Call read_mainout first! Your script might not work now')


    def get_radius(self):
        """
        Get the radius [km] from mainout.dat, it also returns the time
        """
        if len(self.__mainout_time) != 0:
            return self.__mainout_time, self.__mainout_rad
        else:
            print('Call read_mainout first! Your script might not work now')


    def get_ylight(self):
        """
        Get lighter elements (A<=4) from mainout.dat, it also returns the time
        """
        if len(self.__mainout_time) != 0:
            return self.__mainout_time, self.__mainout_ylight
        else:
            print('Call read_mainout first! Your script might not work now')


    def get_yheavy(self):
        """
        Get heavier elements (A>4) from mainout.dat, it also returns the time
        """
        if len(self.__mainout_time) != 0:
            return self.__mainout_time, self.__mainout_yheavy
        else:
            print('Call read_mainout first! Your script might not work now')


    def get_zbar(self):
        """
        Get average proton number from mainout.dat, it also returns the time
        """
        if len(self.__mainout_time) != 0:
            return self.__mainout_time, self.__mainout_zbar
        else:
            print('Call read_mainout first! Your script might not work now')


    def get_abar(self):
        """
        Get average mass number from mainout.dat, it also returns the time
        """
        if len(self.__mainout_time) != 0:
            return self.__mainout_time, self.__mainout_abar
        else:
            print('Call read_mainout first! Your script might not work now')


    def get_abar_at_time(self,t):
        """
        Get average mass number from mainout.dat at time t
        """
        if len(self.__mainout_time) != 0:
            interfunc = interp1d(self.__mainout_time,self.__mainout_abar)
            return interfunc(t)
        else:
            print('Call read_mainout first! Your script might not work now')

    def get_entropy(self):
        """
        Get entropy [kB/baryon] from mainout.dat, it also returns the time
        """
        if len(self.__mainout_time) != 0:
            return self.__mainout_time, self.__mainout_entropy
        else:
            print('Call read_mainout first! Your script might not work now')


    def get_timescales(self):
        """
        Returns the timescales of the run
        Ouput is:
            Time, Temperature[GK], ng, gn, pg, gp, np, pn, an, na, beta (all [1/s])
        """
        if len(self.__timescale_time) == 0:
            print('Call read_timescales first !! The script might not work now!')
            return

        return self.__timescale_time, self.__timescale_temp,self.__timescale_tau_ga , self.__timescale_tau_ag, self.__timescale_tau_ng , self.__timescale_tau_gn, self.__timescale_tau_pg, self.__timescale_tau_gp, self.__timescale_tau_np, self.__timescale_tau_pn, self.__timescale_tau_an, self.__timescale_tau_na, self.__timescale_tau_beta


    def get_final_Z_Y(self):
        """
          Get abundances over proton number
          Note: read_finab has to be called in advance
        """

        if len(self.__finab_nuclei) == 0:
            print('Call read_finab first !! The script might not work now!')
            return

        #Sort it
        for nuc in self.__finab_nuclei:
            nuc.set_sortcriteria('Z')
        self.__finab_nuclei.sort()


        old_Z = -1

        fin_Z = []
        fin_Y = []

        count = -1

        for nuc in self.__finab_nuclei:
            if nuc.get_Z() == old_Z:
                fin_Y[count] += nuc.get_Y()
            else:
                fin_Z.append(nuc.get_Z())
                fin_Y.append(nuc.get_Y())
                count += 1

            old_Z = nuc.get_Z()

        return np.array(fin_Z),np.array(fin_Y)


    def get_final_Z_X(self):
        """
          Get mass fractions over proton number
          Note: read_finab has to be called in advance
        """

        if len(self.__finab_nuclei) == 0:
            print('Call read_finab first !! The script might not work now!')
            return

        #Sort it
        for nuc in self.__finab_nuclei:
            nuc.set_sortcriteria('Z')
        self.__finab_nuclei.sort()

        old_Z = -1

        fin_Z = []
        fin_X = []

        count = -1

        for nuc in self.__finab_nuclei:

            if nuc.get_Z() == old_Z:
                fin_X[count] += nuc.get_X()
            else:
                fin_Z.append(nuc.get_Z())
                fin_X.append(nuc.get_X())
                count += 1

            old_Z = nuc.get_Z()

        return fin_Z,fin_X


    def __get_x_to_seed(self,temperature,x):
        """
          Helper function to calculate Y(neutron) Y(proton) Y(alpha) to seed ratio without making too many repetitions
          Input
            Temperature - temperature in GK for which the ratio is calculated
            x           - Y(n), Y(p) or Y(alpha) of the mainout
          Note: read_mainout should get called first
        """

        mainout_length = len(self.__mainout_iteration)

        # Be sure that mainout was read in
        if mainout_length == 0:
            print('Call read_mainout first! Returning!')
            return

        reverse_temp = self.__mainout_temp[::-1]

        # Get the correct index temperature
        reversed_index = 0
        for temp in reverse_temp:
            reversed_index += 1
            if temp >= temperature: # Assume monotonic increasing of the reversed temperature array
                break

        # Transform reverse index to forward one
        index = mainout_length - reversed_index - 1


        # Interpolate (linearly)
        x_int      = interp1d(self.__mainout_temp[index:mainout_length],x[index:mainout_length])
        yheavy_int = interp1d(self.__mainout_temp[index:mainout_length],self.__mainout_yheavy[index:mainout_length])
        # Get the correct value at the input temperature
        x_at_temp      = x_int(temperature)
        yheavy_at_temp = yheavy_int(temperature)
        # Finally calculate the ratio
        x_to_seed = np.log10(x_at_temp / yheavy_at_temp)


        return x_to_seed

    def get_alpha_to_seed(self,temperature):
        """
          Get the alpha to seed ratio at given temperature [GK]
          Note: read_mainout has get called first
        """
        alpha_to_seed = self.__get_x_to_seed(temperature,self.__mainout_ya) # Call helper function
        return alpha_to_seed

    def get_n_to_seed(self,temperature):
        """
          Get the alpha to seed ratio at given temperature [GK]
          Note: read_mainout has get called first
        """
        neutron_to_seed = self.__get_x_to_seed(temperature,self.__mainout_yn) # Call helper function
        return neutron_to_seed

    def get_p_to_seed(self,temperature):
        """
          Get the alpha to seed ratio at given temperature [GK]
          Note: read_mainout has get called first
        """
        proton_to_seed = self.__get_x_to_seed(temperature,self.__mainout_yp) # Call helper function
        return proton_to_seed

    def get_final_abar(self):
        """
          Get the final abar from finab.dat. Therefore one has to call read_finab first.
        """
        #abar = sum(Xi)/sum(Yi)
        abar = sum(np.array(self.__finab_mass) * np.array(self.__finab_ab) ) / sum(np.array(self.__finab_ab))
        return abar

    def get_final_abar_heavy(self,min_mass):
        """
          Get the final abar calculated for A>min_mass from finab.dat. Therefore one has to call read_finab first.
        """
        #abar = sum(Xi)/sum(Yi)
        for nuc in self.__finab_nuclei:
            nuc.set_sortcriteria('A')

        nuc = sorted(self.__finab_nuclei)

        a_tmp = [x.get_A() for x in nuc]
        Y_tmp = [x.get_Y() for x in nuc]

        for i in range(len(a_tmp)):
            if a_tmp[i] > min_mass:
                break

        abar_heavy = sum(np.array(a_tmp[i:]) * np.array(Y_tmp[i:]) ) / sum(np.array(Y_tmp[i:]))

        return abar_heavy


    def get_finab_nuclei(self):
        """
          Returns a list containing the nuclei contained in finab (list of instance nucleus)
        """
        return self.__finab_nuclei

    def get_tracer_nr(self):
        """
          Returns the tracer number of the tracer of the current instance. If there is no digit in the foldername of the run, the return value is -1.
          e.g. runpath is abc/dfg/tracer_4523.dat -> return value is 4523
        """
        run = os.path.basename(os.path.normpath(self.__path))
        tr_nr = ''
        for s in run:
            if s.isdigit():
                tr_nr = tr_nr+s
        # Try to convert it to an integer
        try:
            tr_nr = int(tr_nr)
        except:
            tr_nr = -1

        return tr_nr


    def plot_final_isotopes(self,figure=None,axlabel=True,lower_limit=1e-10,isotopes=None,isotope_label=True,ignore_isotope_names=None,text_clipping=True,\
    ytext_pos_factor= 1.2,xtext_pos_offset=-0.5,ytext_pos_offset=0.0,nuc_input=None,fig_is_ax=False,color_dict=None,**kwargs):
        """
          Plot the massfraction of the final isotopes over mass number.
          Input:
            lower_limit             : lowest massfraction that will be plotted (default 10^-10)
            isotopes                : isotopes to be plotted. If none (default), all isotopes of finab will get plotted. Type is list that contains the Z number
            isotope_label           : Plot the labels of the isotopes?
            ignore_isotope_names    : which name label should be ignored? (list contains Z number)
            text_clipping           : Text clipping yes or no? ( Should the text shown outside the plotting area?
            ytext_pos_factor        : offset factor for the y position of the labels. Only has an effect when isotope_label is true
            xtext_pos_offset        : offset for the x position of the labels. Only has an effect when isotope_label is true
            nuc_input               : Input of nuclei. If None is given, finab nuclei are used. (instance of list of nucleus)
            fig_is_ax               : parameter figure either a figure object or an axes object?
            color_dict              : Dictionary with atomic numbers and colors.
        """


        # Make a figure if no figure was given
        if figure == None:
            final_fig = plt.figure()
        else:
            final_fig = figure

        if not fig_is_ax:
            # Get the axis instance from figure
            ax = final_fig.gca()
        else:
            ax = figure

        if nuc_input is None:
            plotted_nuclei = self.__finab_nuclei
        else:
            plotted_nuclei = nuc_input

        # Sort the nuclei (atomic number criteria)
        for nuc in plotted_nuclei:
            nuc.set_sortcriteria('Z')
        nuc = sorted(plotted_nuclei)

        # Go through nuclei and make Z-cluster
        old_Z = nuc[0].get_Z()
        A_cluster = []
        X_cluster = []
        for n in nuc:
            current_Z = n.get_Z()
            current_A = n.get_A()
            current_X = n.get_X()

            # Cycle for to small mass fractions
            if current_X <= lower_limit:
                continue


            if current_Z != old_Z: # Plot it

                # Cycle for uninteressting isotopes
                if isotopes != None:
                    if old_Z not in isotopes:
                        A_cluster = [current_A]
                        X_cluster = [current_X]
                        old_Z     = current_Z
                        continue

                if color_dict is None:
                    pl = ax.plot(A_cluster,X_cluster,marker='.',**kwargs)
                else:
                    pl = ax.plot(A_cluster,X_cluster,marker='.',color=color_dict[int(old_Z)],**kwargs)

                # Set isotope labels
                if isotope_label:
                    if X_cluster != []:
                        if (ignore_isotope_names == None) or (not old_Z in ignore_isotope_names):
                            pos_y       = max(X_cluster)
                            pos_x       = A_cluster[X_cluster.index(max(X_cluster))]
                            if color_dict is None:
                                color     = pl[-1].get_color() # Get current color
                                ax.text(pos_x+xtext_pos_offset,pos_y*ytext_pos_factor+ytext_pos_offset,nam,color=color,fontsize=15,ha='center').set_clip_on(text_clipping)
                            else:
                                ax.text(pos_x+xtext_pos_offset,pos_y*ytext_pos_factor+ytext_pos_offset,nam,color=color_dict[int(old_Z)],fontsize=15,ha='center').set_clip_on(text_clipping)


                A_cluster = [current_A]
                X_cluster = [current_X]
                old_Z     = current_Z
            else:                  # Collect data
                A_cluster.append(current_A)
                X_cluster.append(current_X)

            # Get the name and make sure that the first letter is in uppercase
            nam = n.get_elementname()
            nam = nam[0].upper()+nam[1:]




        # Set the labels
        if axlabel:
            ax.set_ylabel('Mass fraction X')
            ax.set_xlabel('Mass number A')
            ax.set_yscale('log')

        return final_fig


    def get_final_A_X(self):
        """
          Get mass fractions over proton number
          Note: read_finab has to be called in advance
        """

        if len(self.__finab_nuclei) == 0:
            print('Call read_finab first !! The script might not work now!')
            return

        #Sort it
        for nuc in self.__finab_nuclei:
            nuc.set_sortcriteria('A')
        self.__finab_nuclei.sort()

        old_A = -1

        fin_A = []
        fin_X = []

        count = -1

        for nuc in self.__finab_nuclei:

            if nuc.get_A() == old_A:
                fin_X[count] += nuc.get_X()
            else:
                fin_A.append(nuc.get_A())
                fin_X.append(nuc.get_X())
                count += 1

            old_A = nuc.get_A()

        return fin_A,fin_X

    def get_final_Z_A_X(self):
        """
          Get mass fractions over proton number
          Note: read_finab has to be called in advance
        """

        if len(self.__finab_nuclei) == 0:
            print('Call read_finab first !! The script might not work now!')
            return

        #Sort it
        for nuc in self.__finab_nuclei:
            nuc.set_sortcriteria('A')
        self.__finab_nuclei.sort()

        old_A = -1

        fin_A = []
        fin_Z = []
        fin_X = []

        count = -1

        for nuc in self.__finab_nuclei:

            fin_A.append(nuc.get_A())
            fin_Z.append(nuc.get_Z())
            fin_X.append(nuc.get_X())
            count += 1

            old_A = nuc.get_A()
        fin_A = np.array(fin_A)
        fin_Z = np.array(fin_Z)
        fin_X = np.array(fin_X)
        return fin_Z,fin_A,fin_X


    def get_final_A_Y(self):
        """
          Get Abundance over proton number
          Note: read_finab has to be called in advance
        """

        if len(self.__finab_nuclei) == 0:
            print('Call read_finab first !! The script might not work now!')
            return

        #Sort it
        for nuc in self.__finab_nuclei:
            nuc.set_sortcriteria('A')
        self.__finab_nuclei.sort()

        old_A = -1

        fin_A = []
        fin_Y = []

        count = -1

        for nuc in self.__finab_nuclei:

            if nuc.get_A() == old_A:
                fin_Y[count] += nuc.get_Y()
            else:
                fin_A.append(nuc.get_A())
                fin_Y.append(nuc.get_Y())
                count += 1

            old_A = nuc.get_A()

        return fin_A,fin_Y



    def get_final_Z_eps(self):
        """
          Get abundances over proton number
          Note: read_finab has to be called in advance
        """

        if len(self.__finab_nuclei) == 0:
            print('Call read_finab first !! The script might not work now!')
            return

        fin_Z,fin_Y = self.get_Z_Y()

        #log eps_A = log (Y_A/Y_H) + 12
        Y_eps = [np.log10(x / self.__finab_hydrogen.get_Y()) +12. for x in fin_Y]


        return fin_Z,Y_eps


    def get_time_evolution(self,nuc):
        """
          Return time and mass fraction for a given nucleus.
          Note: read_snapshots has to be called in advance!
        """

        if isinstance(nuc, str):
            nuc = nucleus(nuc)

        p = nuc.get_Z()
        n = nuc.get_N()


        for i in range(len(self.__neutrons)):
            current_n = self.__neutrons[i]
            current_p = self.__protons[i]


            if (current_n == n and current_p == p):
                index = i
                break
        try:
            massfraction_out = np.array(self.__allmassfractions)[:,index]
        except:
            raise Exception(nuc.get_name(), 'not in list!')

        return self.__time , massfraction_out



    def set_finab_nuclei(self,nucleilist):
        """
          Set the final abundances with a list of nuclei
          Input:
            nucleilist  - List of instance of nuclei
        """


        # Super stupid to save it so often, in principle self.__finab_nuclei contains all these informations!
        self.__finab_mass       = [x.get_A() for x in nucleilist]
        self.__finab_protons    = [x.get_Z() for x in nucleilist]
        self.__finab_neutrons   = [x.get_N() for x in nucleilist]
        self.__finab_ab         = [x.get_Y() for x in nucleilist]
        # Save the instance of nuclei
        self.__finab_nuclei = nucleilist
        self.__finab_nuclei.sort()

        #Get hydrogen out of it
        for nuc in self.__finab_nuclei:
            if (nuc.get_Z() == 1 and nuc.get_N() == 0):
                self.__finab_hydrogen = nuc
                break


    def plot_engen(self,figure=None,axlabel=True,fig_is_ax=False,**kwargs):
        """
        Plot the final massfraction over mass number
        Returns a figure object
        """
        # Make a figure if no figure was given
        if figure == None:
            final_fig = plt.figure()
        else:
            final_fig = figure

        if fig_is_ax:
            ax = final_fig
        else:
            # Get the axis instance from figure
            ax = final_fig.gca()

        ax.set_yscale('log')
        ax.set_xscale('log')
        ax.set_ylabel('Generated energy [erg/g/s]')
        time = self.__engen[0]
        for ind,k in enumerate(self.__engen_keywords):
            if ind==0:
                continue
            elif ind==1:
                ax.plot(time,self.__engen[ind],label=k,lw=2,ls='--')
            else:
                ax.plot(time,self.__engen[ind],label=k)
        return final_fig

    def plot_final_A_X(self,figure=None,axlabel=True,fig_is_ax=False,**kwargs):
        """
        Plot the final massfraction over mass number
        Returns a figure object
        """
        # Make a figure if no figure was given
        if figure == None:
            final_fig = plt.figure()
        else:
            final_fig = figure

        if fig_is_ax:
            ax = final_fig
        else:
            # Get the axis instance from figure
            ax = final_fig.gca()

        # Collect data
        A,X = self.get_final_A_X()

        # Plot data
        ax.plot(A,X,**kwargs)
        ax.set_yscale('log')

        # Make labels
        if axlabel:
            ax.set_ylabel('Mass fraction')
            ax.set_xlabel('Mass number')

        return final_fig


    def plot_final_Z_X(self,figure=None,axlabel=True,fig_is_ax=False,**kwargs):
        """
        Plot the final massfraction over mass number
        Returns a figure object
        """

        # Make a figure if no figure was given
        if figure == None:
            final_fig = plt.figure()
        else:
            final_fig = figure

        if fig_is_ax:
            ax = final_fig
        else:
            # Get the axis instance from figure
            ax = final_fig.gca()

        # Collect data
        Z,X = self.get_final_Z_X()

        # Plot data
        ax.plot(Z,X,**kwargs)
        ax.set_yscale('log')

        # Make labels
        if axlabel:
            ax.set_ylabel('Mass fraction')
            ax.set_xlabel('Atomic number')

        return final_fig


    def plot_final_A_Y(self,figure=None,axlabel=True,fig_is_ax=False,**kwargs):
        """
        Plot the final abundances over mass number
        Returns a figure object
        """

        # Make a figure if no figure was given
        if figure == None:
            final_fig = plt.figure()
        else:
            final_fig = figure

        if fig_is_ax:
            ax = final_fig
        else:
            # Get the axis instance from figure
            ax = final_fig.gca()

        # Collect data
        A,Y = self.get_final_A_Y()

        # Plot data
        ax.plot(A,Y,**kwargs)
        ax.set_yscale('log')

        # Make labels
        if axlabel:
            ax.set_ylabel('Abundance')
            ax.set_xlabel('Mass number')

        return final_fig


    def plot_final_Z_Y(self,figure=None,axlabel=True,fig_is_ax=False,**kwargs):
        """
        Plot the final abundances over mass number
        Returns a figure object
        """

        # Make a figure if no figure was given
        if figure == None:
            final_fig = plt.figure()
        else:
            final_fig = figure

        if fig_is_ax:
            ax = final_fig
        else:
            # Get the axis instance from figure
            ax = final_fig.gca()

        # Collect data
        Z,Y = self.get_final_Z_Y()

        # Plot data
        ax.plot(Z,Y,**kwargs)
        ax.set_yscale('log')

        # Make labels
        if axlabel:
            ax.set_ylabel('Abundance')
            ax.set_xlabel('Atomic number')

        return final_fig



    def plot_temperature(self,figure=None,axlabel=True,**kwargs):
        """
        Plot the temperature [GK]
        Returns a figure object
        """

        # Make a figure if no figure was given
        if figure == None:
            final_fig = plt.figure()
        else:
            final_fig = figure

        # Get the axis instance from figure
        ax = final_fig.gca()

        # Collect data
        time,temp = self.get_temperature()

        # Plot data
        ax.plot(time,temp,**kwargs)
        ax.set_yscale('log')

        # Make labels
        if axlabel:
            ax.set_ylabel('Temperature [GK]')
            ax.set_xlabel('Time [s]')

        return final_fig


    def plot_density(self,figure=None,axlabel=True,**kwargs):
        """
        Plot the density [g/ccm]
        Returns a figure object
        """

        # Make a figure if no figure was given
        if figure == None:
            final_fig = plt.figure()
        else:
            final_fig = figure

        # Get the axis instance from figure
        ax = final_fig.gca()

        # Collect data
        time,dens = self.get_density()

        # Plot data
        ax.plot(time,dens,**kwargs)
        ax.set_yscale('log')

        # Make labels
        if axlabel:
            ax.set_ylabel('Density [g/ccm]')
            ax.set_xlabel('Time [s]')

        return final_fig


    def plot_abar(self,figure=None,axlabel=True,**kwargs):
        """
        Plot average mass number
        Returns a figure object
        """

        # Make a figure if no figure was given
        if figure == None:
            final_fig = plt.figure()
        else:
            final_fig = figure

        # Get the axis instance from figure
        ax = final_fig.gca()

        # Collect data
        time,abar = self.get_abar()

        # Plot data
        ax.plot(time,abar,**kwargs)

        # Make labels
        if axlabel:
            ax.set_ylabel(r'<$\bar{A}$>')
            ax.set_xlabel('Time [s]')

        return final_fig


    def plot_zbar(self,figure=None,axlabel=True,**kwargs):
        """
        Plot average proton number
        Returns a figure object
        """

        # Make a figure if no figure was given
        if figure == None:
            final_fig = plt.figure()
        else:
            final_fig = figure

        # Get the axis instance from figure
        ax = final_fig.gca()

        # Collect data
        time,zbar = self.get_zbar()

        # Plot data
        ax.plot(time,zbar,**kwargs)

        # Make labels
        if axlabel:
            ax.set_ylabel(r'<$\bar{Z}$>')
            ax.set_xlabel('Time [s]')

        return final_fig


    def plot_radius(self,figure=None,axlabel=True,**kwargs):
        """
        Plot the radius
        Returns a figure object
        """

        # Make a figure if no figure was given
        if figure == None:
            final_fig = plt.figure()
        else:
            final_fig = figure

        # Get the axis instance from figure
        ax = final_fig.gca()

        # Collect data
        time,radius = self.get_radius()

        # Plot data
        ax.plot(time,radius,**kwargs)

        # Make labels
        if axlabel:
            ax.set_ylabel(r'Radius [km]')
            ax.set_xlabel('Time [s]')

        return final_fig


    def plot_trajectory(self,figure=None,**kwargs):
        """
        Plot all quantities of the trajectory from trajectory file
        """

        # Make a figure if no figure was given
        if figure == None:
            final_fig = plt.figure()
        else:
            final_fig = figure


        # Get the axis instance from figure
        gs           = gridspec.GridSpec(2,3)
        radiusplot   = plt.subplot(gs[0,0])
        densplot     = plt.subplot(gs[0,1])
        tempplot     = plt.subplot(gs[0,2])
        yeplot       = plt.subplot(gs[1,0])
        entrplot     = plt.subplot(gs[1,1])
        dentempplot  = plt.subplot(gs[1,2])

        plt.subplots_adjust(hspace=0.5,wspace=0.5)

        # Collect data
        time,radius,dens,temp,ye,entr = self.get_trajectory()

        radiusplot.plot(time,radius,**kwargs)
        radiusplot.set_ylabel('Radius [cm]')
        radiusplot.set_xlabel('Time [s]')

        densplot.plot(time,dens,**kwargs)
        densplot.set_ylabel('Density [g/ccm]')
        densplot.set_xlabel('Time [s]')

        tempplot.plot(time,temp,**kwargs)
        tempplot.set_ylabel('Temperature [GK]')
        tempplot.set_xlabel('Time [s]')

        yeplot.plot(time,ye,**kwargs)
        yeplot.set_ylabel('Ye')
        yeplot.set_xlabel('Time [s]')

        entrplot.plot(time,entr,**kwargs)
        entrplot.set_ylabel('Entropy [kB/nucleon]')
        entrplot.set_xlabel('Time [s]')

        dentempplot.plot(dens,temp,**kwargs)
        dentempplot.set_ylabel('Temperature [GK]')
        dentempplot.set_xlabel('Density [g/ccm]')

        return final_fig





    def plot_flow_range(self,start=1,end=0,outputpath='./',threads=5,silence=False,plotaxis=[0,0,0,0],ilabel=1,imlabel=0,imagic=0,prange=4,iplot=2):
        """
          Plot all flow files in a range from start to end. If end is 0 (default) all flow files will get plotted.
          The files will added to outputpath.
          This routine is also parallized, threads determines the number of cores used.
          Silence: should the routine print something?
        """
        # TODO parallize it

        if outputpath[-1] != '/':
            outputpath += '/'

        if end==0:
            end = len(os.listdir(self.__path+'flow/'))

        for i in range(start,end+1,1):
            fig = self.plot_flow(i,plotaxis=[1,20,1,20],ilabel=1,imlabel=0,imagic=0,prange=4,iplot=2)
            fig.savefig(outputpath+'flow_{n:04d}'.format(n=i))


    def plot_integrated_flow(self,figure=None,fig_is_ax=False,plotaxis=[0,0,0,0],ilabel=1,imlabel=0,imagic=0,prange=4,iplot=2):
        """
         Plot flows over all timesteps
        """
        if figure == None:
            final_fig = plt.figure()
        else:
            final_fig = figure

        if fig_is_ax:
            ax = final_fig
        else:
            # Get the axis instance from figure
            ax = final_fig.gca()

        # Count the flows in the folder
        all_flows = np.array(os.listdir(self.__path+'flow/'))
        flow_names = all_flows[np.array(['flow' in x for x in all_flows])]

        integrated_flow = None

        for flownumber in np.arange(1,len(flow_names)):
            inputfile = self.__path+'flow/flow_{n:04d}.dat'.format(n=flownumber)
            inputfile_next = self.__path+'flow/flow_{n:04d}.dat'.format(n=flownumber+1)
            # Give warning if flow file does not exist
            if (not os.path.exists(inputfile)) or (not os.path.exists(inputfile_next) ):
                print('Warning: flow_{n:04d}.dat'.format(n=flownumber)+' does not exist!'  )
                continue

            # Read the flow file
            ff      = self.__read_flow(inputfile)
            ff_next = self.__read_flow(inputfile_next)

            timestep = ff_next['time']-ff['time']
            print(timestep)
            flow = ff['flow']
            if integrated_flow is None:
                integrated_flow = flow
            else:
                integrated_flow+=flow*timestep

        integrated_flow = np.log10(integrated_flow+1e-99)
        flow = integrated_flow
        sort = flow.argsort()
        flow = flow[sort]
        nin  = ff['nin'][sort]
        zin  = ff['zin'][sort]
        yin  = ff['yin'][sort]
        nout = ff['nout'][sort]
        zout = ff['zout'][sort]
        yout = ff['yout'][sort]
        maxflow = max(flow)
        nzmax = int(max(max(zin),max(zout)))+1
        nnmax = int(max(max(nin),max(nout)))+1
        nzycheck = np.zeros([nnmax,nzmax,3])
        for i in range(len(nin)):
          ni = int(nin[i])
          zi = int(zin[i])
          nzycheck[ni,zi,0] = 1
          nzycheck[ni,zi,1] = yin[i]
          if iplot==1 or iplot==2:
            no = int(nout[i])
            zo = int(zout[i])
            nzycheck[no,zo,0] = 1
            nzycheck[no,zo,1] = yout[i]
            if iplot==1 and flow[i]>maxflow-prange:
              nzycheck[ni,zi,2] = 1
              nzycheck[no,zo,2] = 1

        #### create plot

        ## define axis and plot style (colormap, size, fontsize etc.)
        if plotaxis==[0,0,0,0]:
          xdim=10
          ydim=6
        else:
          dx = plotaxis[1]-plotaxis[0]
          dy = plotaxis[3]-plotaxis[2]
          ydim = 10
          xdim = ydim*dx/dy

        format = 'png'#'pdf'


        # color map choice for abundances
        cmapa = cm.jet
        # color map choice for arrows
        cmapr = cm.jet#autumn
        # if a value is below the lower limit its set to white
        cmapa.set_under(color='w')
        cmapr.set_under(color='w')
        # set value range for abundance colors (log10(Y))
        norma = colors.Normalize(vmin=-10,vmax=-3)
        # declare scalarmappable for abundances
        a2c = cm.ScalarMappable(norm=norma,cmap=cmapa)
        a=[]
        a2c.set_array(np.array(a))
        # set x- and y-axis scale aspect ratio to 1
        ax.set_aspect('equal')
        #print time,temp and density on top
        temp = '%8.3e' %ff['temp']
        time = '%8.3e' %ff['time']
        dens = '%8.3e' %ff['dens']


        ## only abundance plotted
        if iplot==0 or iplot==2:
          patches = []
          color = []

          for i in range(nzmax):
            for j in range(nnmax):
              if nzycheck[j,i,0]==1:
                # abundance
                yab = nzycheck[j,i,1]
                col = a2c.to_rgba(np.log10(yab))
                xy = j-0.5,i-0.5
                rect = Rectangle(xy,1,1,ec='k',color=col,picker=True)
                rect.set_zorder(1)
                ax.add_patch(rect)

          cb=plt.colorbar(a2c)

          # colorbar label
          cb.set_label('log$_{10}$(Y)')

        # Add black frames for stable isotopes
        f = open(os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))+'/data/stableiso.dat')

        #head = f.readline()
        stable = []

        for line in f.readlines():
          tmp = line.split()
          zz = int(tmp[2])
          nn = int(tmp[1])
          xy = nn,zz
          circ = Circle(xy,radius=0.1,fc='k')
          circ.set_zorder(2)
          ax.add_patch(circ)

        if iplot==1 or iplot==2:
          apatches = []
          acolor = []
          m = 0.8/prange
          vmax=np.ceil(max(flow))
          vmin=max(flow)-prange
          b=-vmin*m+0.1
          normr = colors.Normalize(vmin=vmin,vmax=vmax)
          ymax=0.
          xmax=0.

          for i in range(len(nin)):
            x = nin[i]
            y = zin[i]
            dx = nout[i]-nin[i]
            dy = zout[i]-zin[i]

            if flow[i]>=vmin:
              arrowwidth = flow[i]*m+b
              arrow = Arrow(x,y,dx,dy, width=arrowwidth)
              if xmax<x:
                xmax=x
              if ymax<y:
                ymax=y
              acol = flow[i]
              apatches.append(arrow)
              acolor.append(acol)

              if iplot==1:
                xy = x-0.5,y-0.5
                rect = Rectangle(xy,1,1,ec='k',fc='None',fill='False',lw=1.)
                rect.set_zorder(2)
                ax.add_patch(rect)
                xy = x+dx-0.5,y+dy-0.5
                rect = Rectangle(xy,1,1,ec='k',fc='None',fill='False',lw=1.)
                rect.set_zorder(2)
                ax.add_patch(rect)


          a = PatchCollection(apatches, cmap=cmapr, norm=normr)
          a.set_array(np.array(acolor))
          a.set_zorder(3)
          ax.add_collection(a)
          cb = plt.colorbar(a)

          # colorbar label
          cb.set_label('log$_{10}$(f)')

        # decide which array to take for label positions
        iarr = 0
        if iplot==1: iarr=2

        # plot element labels
        if ilabel==1:
          for z in range(nzmax):
            try:
              nmin = min(np.argwhere(nzycheck[:,z,iarr]))[0]-1
              ax.text(nmin,z,self.__elname[z],horizontalalignment='center',verticalalignment='center',fontsize='medium',clip_on=True)
            except:
              continue

        # plot mass numbers
        if imlabel==1:
          for z in range(nzmax):
             for n in range(nnmax):
                a = z+n
                if nzycheck[n,z,iarr]==1:
                  ax.text(n,z,a,horizontalalignment='center',verticalalignment='center',fontsize='small',clip_on=True)

        # plot lines at magic numbers
        if imagic==1:
          ixymagic=[2, 8, 20, 28, 50, 82, 126]
          nmagic = len(ixymagic)
          for magic in ixymagic:
            if magic<=nzmax:
              try:
                xnmax = max(np.argwhere(nzycheck[:,magic,iarr]))[0]
                line = ax.plot([-0.5,xnmax+0.5],[magic-0.5,magic-0.5],lw=3.,color='r',ls='--')
                line = ax.plot([-0.5,xnmax+0.5],[magic+0.5,magic+0.5],lw=3.,color='r',ls='--')
              except ValueError:
                dummy=0
            if magic<=nnmax:
              try:
                yzmax = max(np.argwhere(nzycheck[magic,:,iarr]))[0]
                line = ax.plot([magic-0.5,magic-0.5],[-0.5,yzmax+0.5],lw=3.,color='r',ls='--')
                line = ax.plot([magic+0.5,magic+0.5],[-0.5,yzmax+0.5],lw=3.,color='r',ls='--')
              except ValueError:
                dummy=0

        # set axis limits
        if plotaxis==[0,0,0,0]:
          if iplot==2 or iplot==0:
            xmax=max(nin)
            ymax=max(zin)
          ax.axis([-0.5,xmax+0.5,-0.5,ymax+0.5])
        else:
          ax.axis([plotaxis[0]-0.5,plotaxis[1]+0.5,plotaxis[2]-0.5,plotaxis[3]+0.5])

        # set x- and y-axis label
        ax.set_xlabel('neutron number')
        ax.set_ylabel('proton number')

        def onpick1(event):
        #  print event.artist
          if isinstance(event.artist, Rectangle):
            patch = event.artist
            pn = int(patch.get_x()+0.5)
            pz = int(patch.get_y()+0.5)
            pab = '%8.3e' %nzycheck[pn,pz,1]
            tmp = (pn+pz)*nzycheck[pn,pz,1]
            pmf = '%8.3e' %tmp
            print( self.__elname[pz] + str(pn+pz) + '  , Y = ' + pab + '  X = ' + pmf )

        return final_fig

    def plot_flow(self,flownumber,plotaxis=[0,0,0,0],ilabel=1,imlabel=0,imagic=0,prange=4,iplot=2):
        """
          Plot a flow with a given flownumber.

          Originally from Christian Winteler and Urs Frischknecht
          Edited by Moritz Reichert (2.4.2017)

          Plotting options are:
          plotaxis  :   Set axis limit: If default [0,0,0,0] the complete range in (N,Z) will be plotted, i.e. all isotopes, else specify the limits in plotaxis = [xmin,xmax,ymin,ymax]
          ilabel    :   elemental labels off/on [0/1]
          imlabel   :   label for isotopic masses off/on [0/1]
          prange    :   flow is plotted over "prange" dex. If flow < maxflow-prange it is not plotted
          iplot     :   plot handler abundance/flux/abundance+flux plot [0/1/2]
        """


        inputfile = self.__path+'flow/flow_{n:04d}.dat'.format(n=flownumber)

        ff = self.__read_flow(inputfile)
        try:
          flow = np.log10(ff['flow']+1.e-99)
          sort = flow.argsort()
          flow = flow[sort]
          nin  = ff['nin'][sort]
          zin  = ff['zin'][sort]
          yin  = ff['yin'][sort]
          nout = ff['nout'][sort]
          zout = ff['zout'][sort]
          yout = ff['yout'][sort]
          maxflow = max(flow)
          nzmax = int(max(max(zin),max(zout)))+1
          nnmax = int(max(max(nin),max(nout)))+1
        except KeyError:
          if iplot==1 or iplot==2:
            sys.exit('The read file does not contain flow data')
          nin  = ff['nin']
          zin  = ff['zin']
          yin  = ff['yin']
          nnmax = int(max(nin))+1
          nzmax = int(max(zin))+1

        nzycheck = np.zeros([nnmax,nzmax,3])
        for i in range(len(nin)):
          ni = int(nin[i])
          zi = int(zin[i])
          nzycheck[ni,zi,0] = 1
          nzycheck[ni,zi,1] = yin[i]
          if iplot==1 or iplot==2:
            no = int(nout[i])
            zo = int(zout[i])
            nzycheck[no,zo,0] = 1
            nzycheck[no,zo,1] = yout[i]
            if iplot==1 and flow[i]>maxflow-prange:
              nzycheck[ni,zi,2] = 1
              nzycheck[no,zo,2] = 1

        #### create plot

        ## define axis and plot style (colormap, size, fontsize etc.)
        if plotaxis==[0,0,0,0]:
          xdim=10
          ydim=6
        else:
          dx = plotaxis[1]-plotaxis[0]
          dy = plotaxis[3]-plotaxis[2]
          ydim = 10
          xdim = ydim*dx/dy

        format = 'png'#'pdf'
        params = {'axes.labelsize':  15,
                  'text.fontsize':   15,
                  'legend.fontsize': 15,
                  'xtick.labelsize': 15,
                  'ytick.labelsize': 15,
                  'text.usetex': True}
        plt.rcParams.update(params)
        fig=plt.figure(figsize=(xdim,ydim),dpi=100)
        axx = 0.06
        axy = 0.10
        axw = 0.90
        axh = 0.8
        ax=plt.axes([axx,axy,axw,axh])

        # color map choice for abundances
        cmapa = cm.jet
        # color map choice for arrows
        cmapr = cm.jet#autumn
        # if a value is below the lower limit its set to white
        cmapa.set_under(color='w')
        cmapr.set_under(color='w')
        # set value range for abundance colors (log10(Y))
        norma = colors.Normalize(vmin=-10,vmax=-3)
        # declare scalarmappable for abundances
        a2c = cm.ScalarMappable(norm=norma,cmap=cmapa)
        a=[]
        a2c.set_array(np.array(a))
        # set x- and y-axis scale aspect ratio to 1
        ax.set_aspect('equal')
        #print time,temp and density on top
        temp = '%8.3e' %ff['temp']
        time = '%8.3e' %ff['time']
        dens = '%8.3e' %ff['dens']

        box1 = TextArea("t : " + time + " s~~/~~T$_{9}$ : " + temp + "~~/~~$\\rho_{b}$ : " \
                      + dens + ' g/cm$^{3}$', textprops=dict(color="k"))
        anchored_box = AnchoredOffsetbox(loc=3,
                        child=box1, pad=0.,
                        frameon=False,
                        bbox_to_anchor=(0., 1.02),
                        bbox_transform=ax.transAxes,
                        borderpad=0.,
                        )
        ax.add_artist(anchored_box)

        ## only abundance plotted
        if iplot==0 or iplot==2:
          patches = []
          color = []

          for i in range(nzmax):
            for j in range(nnmax):
              if nzycheck[j,i,0]==1:
                # abundance
                yab = nzycheck[j,i,1]
                col = a2c.to_rgba(np.log10(yab))
                xy = j-0.5,i-0.5
                rect = Rectangle(xy,1,1,ec='k',color=col,picker=True)
                rect.set_zorder(1)
                ax.add_patch(rect)

          cb=plt.colorbar(a2c)

          # colorbar label
          cb.set_label('log$_{10}$(Y)')

        # Add black frames for stable isotopes
        f = open(os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))+'/data/stableiso.dat')

        #head = f.readline()
        stable = []

        for line in f.readlines():
          tmp = line.split()
          zz = int(tmp[2])
          nn = int(tmp[1])
          xy = nn,zz
          circ = Circle(xy,radius=0.1,fc='k')
          circ.set_zorder(2)
          ax.add_patch(circ)

        if iplot==1 or iplot==2:
          apatches = []
          acolor = []
          m = 0.8/prange
          vmax=np.ceil(max(flow))
          vmin=max(flow)-prange
          b=-vmin*m+0.1
          normr = colors.Normalize(vmin=vmin,vmax=vmax)
          ymax=0.
          xmax=0.

          for i in range(len(nin)):
            x = nin[i]
            y = zin[i]
            dx = nout[i]-nin[i]
            dy = zout[i]-zin[i]

            if flow[i]>=vmin:
              arrowwidth = flow[i]*m+b
              arrow = Arrow(x,y,dx,dy, width=arrowwidth)
              if xmax<x:
                xmax=x
              if ymax<y:
                ymax=y
              acol = flow[i]
              apatches.append(arrow)
              acolor.append(acol)

              if iplot==1:
                xy = x-0.5,y-0.5
                rect = Rectangle(xy,1,1,ec='k',fc='None',fill='False',lw=1.)
                rect.set_zorder(2)
                ax.add_patch(rect)
                xy = x+dx-0.5,y+dy-0.5
                rect = Rectangle(xy,1,1,ec='k',fc='None',fill='False',lw=1.)
                rect.set_zorder(2)
                ax.add_patch(rect)


          a = PatchCollection(apatches, cmap=cmapr, norm=normr)
          a.set_array(np.array(acolor))
          a.set_zorder(3)
          ax.add_collection(a)
          cb = plt.colorbar(a)

          # colorbar label
          cb.set_label('log$_{10}$(f)')

        # decide which array to take for label positions
        iarr = 0
        if iplot==1: iarr=2

        # plot element labels
        if ilabel==1:
          for z in range(nzmax):
            try:
              nmin = min(np.argwhere(nzycheck[:,z,iarr]))[0]-1
              ax.text(nmin,z,self.__elname[z],horizontalalignment='center',verticalalignment='center',fontsize='medium',clip_on=True)
            except:
              continue

        # plot mass numbers
        if imlabel==1:
          for z in range(nzmax):
             for n in range(nnmax):
                a = z+n
                if nzycheck[n,z,iarr]==1:
                  ax.text(n,z,a,horizontalalignment='center',verticalalignment='center',fontsize='small',clip_on=True)

        # plot lines at magic numbers
        if imagic==1:
          ixymagic=[2, 8, 20, 28, 50, 82, 126]
          nmagic = len(ixymagic)
          for magic in ixymagic:
            if magic<=nzmax:
              try:
                xnmax = max(np.argwhere(nzycheck[:,magic,iarr]))[0]
                line = ax.plot([-0.5,xnmax+0.5],[magic-0.5,magic-0.5],lw=3.,color='r',ls='--')
                line = ax.plot([-0.5,xnmax+0.5],[magic+0.5,magic+0.5],lw=3.,color='r',ls='--')
              except ValueError:
                dummy=0
            if magic<=nnmax:
              try:
                yzmax = max(np.argwhere(nzycheck[magic,:,iarr]))[0]
                line = ax.plot([magic-0.5,magic-0.5],[-0.5,yzmax+0.5],lw=3.,color='r',ls='--')
                line = ax.plot([magic+0.5,magic+0.5],[-0.5,yzmax+0.5],lw=3.,color='r',ls='--')
              except ValueError:
                dummy=0

        # set axis limits
        if plotaxis==[0,0,0,0]:
          if iplot==2 or iplot==0:
            xmax=max(nin)
            ymax=max(zin)
          ax.axis([-0.5,xmax+0.5,-0.5,ymax+0.5])
        else:
          ax.axis([plotaxis[0]-0.5,plotaxis[1]+0.5,plotaxis[2]-0.5,plotaxis[3]+0.5])

        # set x- and y-axis label
        ax.set_xlabel('neutron number')
        ax.set_ylabel('proton number')

        def onpick1(event):
        #  print event.artist
          if isinstance(event.artist, Rectangle):
            patch = event.artist
            pn = int(patch.get_x()+0.5)
            pz = int(patch.get_y()+0.5)
            pab = '%8.3e' %nzycheck[pn,pz,1]
            tmp = (pn+pz)*nzycheck[pn,pz,1]
            pmf = '%8.3e' %tmp
            print( self.__elname[pz] + str(pn+pz) + '  , Y = ' + pab + '  X = ' + pmf )

        fig.canvas.mpl_connect('pick_event', onpick1)
        #~ fig.savefig('test.png',dpi=400)
        return fig


    def plot_ye(self,figure=None,axlabel=True,**kwargs):
        """
        Plot average proton number
        Returns a figure object
        """

        # Make a figure if no figure was given
        if figure == None:
            final_fig = plt.figure()
        else:
            final_fig = figure

        # Get the axis instance from figure
        ax = final_fig.gca()

        # Collect data
        time,ye = self.get_ye()

        # Plot data
        ax.plot(time,ye,**kwargs)

        # Make labels
        if axlabel:
            ax.set_ylabel(r'Y$_e$')
            ax.set_xlabel('Time [s]')

        return final_fig


    def plot_entropy(self,figure=None,axlabel=True,**kwargs):
        """
        Plot the entropy [kb/nucleon]
        Returns a figure object
        """

        # Make a figure if no figure was given
        if figure == None:
            final_fig = plt.figure()
        else:
            final_fig = figure

        # Get the axis instance from figure
        ax = final_fig.gca()

        # Collect data
        time,entr = self.get_entropy()

        # Plot data
        ax.plot(time,entr,**kwargs)

        # Make labels
        if axlabel:
            ax.set_ylabel(r'Entropy [kB/nucleon]')
            ax.set_xlabel('Time [s]')

        return final_fig


    def plot_timescales(self,figure=None,axlabel=True,**kwargs):
        """
        Plot the timescales for a Winnet run (ng, gn, pg, gp, np, pn, an, na, beta (all [1/s]))
        """

        # Make a figure if no figure was given
        if figure == None:
            final_fig = plt.figure()
        else:
            final_fig = figure

        # Get the axis instance from figure
        ax = final_fig.gca()

        # Collect data
        ts = self.get_timescales()

        time        = ts[0]
        timescales  = ts[2:]
        labels      = [r'$\gamma$-a',r'a-$\gamma$',r'n-$\gamma$',r'$\gamma$-n',r'p-$\gamma$',r'$\gamma$-p',r'n-p',r'p-n',r'$\alpha$-n',r'n-$\alpha$',r'$\beta$']
        colors      = ['red','blue','green','saddlebrown','cyan','magenta',"tab:orange"]

        # Plot data
        ax.set_yscale('log')
        ax.set_xscale('log')
        for i in range(len(labels)):
            if i % 2 == 0:
                col = colors[int(i/2)]
                ax.plot(time,timescales[i],ls='-',label=labels[i],color=col,**kwargs)
            else:
                ax.plot(time,timescales[i],ls='--',label=labels[i],color=col,**kwargs)

        # Make labels
        if axlabel:
            ax.set_ylabel(r'Timescale [s]')
            ax.set_xlabel('Time [s]')

        return final_fig


    def plot_sunet(self,figure=None):
        """
          Plot the contained nuclei of the run. The information is extracted from param.out
        """
        if len(self.__param_parameter) == 0:
            print('Call read_template first! Your script might not work now!')
            return
        # Get path to sunet file
        ind = self.__param_parameter.index('net_source')
        filewithnuclei = self.__param_data[ind].replace('"','').strip()

        stableisotopes  = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))+'/data/stableiso.dat'

        #Get stable nuclei
        stable_n,stable_z  = np.loadtxt(stableisotopes,unpack = True, usecols=(1,2))

        #Get the nuclei
        nucleilist  = np.loadtxt(filewithnuclei,unpack = True,dtype=str)

        max_p = 0
        max_n = 0

        #make them to nucleiclass
        nucleiclass = []
        for nuc in nucleilist:
            nucleiclass.append(nucleus(nuc))
            if max_p < nucleus(nuc).get_Z():
                max_p = nucleus(nuc).get_Z()
            if max_n < nucleus(nuc).get_N():
                max_n = nucleus(nuc).get_N()

        amountneutrons  = max_n + 1
        amountprotons   = max_p + 1

        # Make a figure if no figure was given
        if figure == None:
            fig = plt.figure(figsize=(amountneutrons/10.,amountprotons/10.))
        else:
            fig = figure

        ax  = fig.gca()

        #~ alphanuclei  = [2,4,6,8,10,12,14,16,18,20,22,24,26,28]
        alphanuclei  = []
        magicnumbers = [8,20,28,50,82,126]

        for nuc in nucleiclass:# plot nuclei
            z       = nuc.get_Z()-0.5
            n       = nuc.get_N()-0.5
            xy      = n,z

            stable = False
            for i in range(len(stable_z)):
                if (stable_z[i] == nuc.get_Z()) and (stable_n[i] == nuc.get_N()):
                    stable = True
                    break

            # Uncomment to color the iron region
            #~ if nuc.get_Z()+nuc.get_N()>=50 and nuc.get_Z()+nuc.get_N()<=62:
                #~ rectcolor = '#6969ae'
            #~ else:
            rectcolor = '#d4e5ff'

            if not stable:
                rect    = Rectangle(xy,1,1,ec='gray',zorder=1,lw=0.5,color=rectcolor)
            else:
                rect    = Rectangle(xy,1,1,ec='k',zorder=2,lw=1,color=rectcolor)

            ax.add_patch(rect)

            if ((nuc.get_N() in alphanuclei) and (nuc.get_Z() == nuc.get_N())):
                xy = n + 0.05,z+0.05
                rect    = Rectangle(xy,0.90,0.90,ec='red',zorder=3,lw=0.7,color=rectcolor)
                ax.add_patch(rect)

        #set the elementlabels
        for z in range(amountprotons):

            min_n = amountneutrons

            for nuc in nucleiclass:
                if (z == nuc.get_Z() and min_n > nuc.get_N()):
                    min_n = nuc.get_N()
                    name  = nuc.get_elementname()
                    if name == 'neutron':
                        name = 'n'
                    if name == 'h':
                        name = 'p'

            if z==0:
                ax.text(min_n - 1 - 0.4*len(name), z-0.5, name, fontsize = 5)
            else:
                if ('l' in name) or ('i' in name):
                    ax.text(min_n - 1 - 0.3*len(name), z-0.2, name, fontsize = 5)
                elif 'm' in name:
                    ax.text(min_n - 1 - 0.5*len(name), z-0.2, name, fontsize = 5)
                else:
                    ax.text(min_n - 1 - 0.4*len(name), z-0.2, name, fontsize = 5)

        #set the magic numbers
        for i in range(len(magicnumbers)):
            magic = magicnumbers[i]

            #neutron magic
            min_z = 10000
            max_z = -1
            min_N = -0.5
            max_N = 10000
            for nuc in nucleiclass:
                if nuc.get_N() == magic:
                    if nuc.get_Z()<min_z:
                        min_z=float(nuc.get_Z())
                    if nuc.get_Z()>max_z:
                        max_z=float(nuc.get_Z())

                if nuc.get_Z() == magic:
                    if nuc.get_N()>min_N:
                        min_N=float(nuc.get_N())
                    if nuc.get_N()<max_N:
                        max_N=float(nuc.get_N())

            width = 0.7
            dash  = (3,2)
            #neutron magic
            ax.plot([magic-0.5, magic-0.5], [min_z-0.5, max_z+0.5], color='k',linestyle='--',lw=width, dashes=dash)
            ax.plot([magic+0.5, magic+0.5], [min_z-0.5, max_z+0.5], color='k',linestyle='--',lw=width, dashes=dash)
            #proton magic
            ax.plot([max_N-0.5,min_N+0.5],[magic-0.5, magic-0.5], color='k',linestyle='--',lw=width, dashes=dash)
            ax.plot([max_N-0.5,min_N+0.5],[magic+0.5, magic+0.5], color='k',linestyle='--',lw=width, dashes=dash)

        #remove the axis
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)

        for axis in ['left','right','top','bottom']:
            ax.spines[axis].set_linewidth(0)

        ax.set_xlim(-1,amountneutrons)
        ax.set_ylim(-1,amountprotons)

        return fig


    def plot_mainout(self,figure=None,**kwargs):
        """
        Plot all information that are contained in mainout.dat
            iteration, time, temperature, density, Ye, Radius, Yn, Yp, Y_alpha, Y_light, Y_heavy, Zbar, Abar, Entropy
        """
        # Make a figure if no figure was given
        if figure == None:
            final_fig = plt.figure(figsize=(15,15))
        else:
            final_fig = figure

        # Get the axis instance from figure
        gs           = gridspec.GridSpec(3,4)
        tempplot     = plt.subplot(gs[0,0])
        densplot     = plt.subplot(gs[0,1],sharex=tempplot)
        yeplot       = plt.subplot(gs[0,2],sharex=tempplot)
        radiusplot   = plt.subplot(gs[0,3],sharex=tempplot)
        ynplot       = plt.subplot(gs[1,0],sharex=tempplot)
        ypplot       = plt.subplot(gs[1,1],sharex=tempplot)
        yaplot       = plt.subplot(gs[1,2],sharex=tempplot)
        ylightplot   = plt.subplot(gs[1,3],sharex=tempplot)
        yheavyplot   = plt.subplot(gs[2,0],sharex=tempplot)
        zbarplot     = plt.subplot(gs[2,1],sharex=tempplot)
        abarplot     = plt.subplot(gs[2,2],sharex=tempplot)
        entrplot     = plt.subplot(gs[2,3],sharex=tempplot)

        plt.subplots_adjust(hspace=0.5,wspace=0.5)

        # Get data
        it,time,temp,dens,ye,rad,yn,yp,ya,ylight,yheavy,zbar,abar,entr = self.get_mainout()

        tempplot.plot(time,temp)
        tempplot.set_ylabel('Temperature [GK]')
        tempplot.set_xlabel('Time [s]')
        tempplot.set_xscale('log')

        densplot.plot(time,dens)
        densplot.set_yscale('log')
        densplot.set_ylabel('Density [g/ccm]')
        densplot.set_xlabel('Time [s]')

        yeplot.plot(time,ye)
        yeplot.set_ylabel('Ye')
        yeplot.set_xlabel('Time [s]')

        radiusplot.plot(time,rad)
        radiusplot.set_yscale('log')
        radiusplot.set_ylabel('Radius [km]')
        radiusplot.set_xlabel('Time [s]')

        ynplot.plot(time,yn)
        ynplot.set_ylabel('Yn')
        ynplot.set_xlabel('Time [s]')
        ynplot.set_yscale('log')

        ypplot.plot(time,yp)
        ypplot.set_ylabel('Yp')
        ypplot.set_xlabel('Time [s]')
        ypplot.set_yscale('log')

        yaplot.plot(time,ya)
        yaplot.set_ylabel('Yalpha')
        yaplot.set_xlabel('Time [s]')
        yaplot.set_yscale('log')

        ylightplot.plot(time,ylight)
        ylightplot.set_ylabel('Ylight')
        ylightplot.set_xlabel('Time [s]')
        ylightplot.set_yscale('log')

        yheavyplot.plot(time,yheavy)
        yheavyplot.set_ylabel('Yheavy')
        yheavyplot.set_xlabel('Time [s]')
        ylightplot.set_yscale('log')

        zbarplot.plot(time,zbar)
        zbarplot.set_ylabel('Zbar')
        zbarplot.set_xlabel('Time [s]')

        abarplot.plot(time,abar)
        abarplot.set_ylabel('Abar')
        abarplot.set_xlabel('Time [s]')

        entrplot.plot(time,entr)
        entrplot.set_ylabel('Entropy [kB/nucleon]')
        entrplot.set_xlabel('Time [s]')

        return final_fig

    def animate_nuclear_chart(self,figure=None,plot_magic=True,time_title=True,min_X=1e-8,max_X=None,cmap= cm.viridis,element_labels=True,**kwargs):
        """
          Make an animation of the mass fractions in the nuclear chart
        """
        # Make a figure if no figure was given
        if figure == None:
            final_fig = plt.figure(figsize=(10,5))
        else:
            final_fig = figure

        ax = final_fig.gca()
        ax.set_aspect(1)

        # Read time from the snapshots if not done already
        if len(self.__time)==0:
            self.read_snapshots()

        stable_n,stable_p,stable_a= np.loadtxt(\
                                    os.path.dirname(os.path.abspath(\
                                    inspect.getfile(inspect.currentframe())))+'/data/stableiso.dat',\
                                    unpack=True,usecols=[1,2,3])
        # Set a norm
        norm = mpl.colors.LogNorm(vmin=min_X, vmax=max_X,clip=False)

        def __init_anim():
            if time_title:
                snaps_time = self.__time[0]
                self.__helper_title = ax.text(0.02,0.95,"t = "+"{:.2e}".format(snaps_time)+" s",transform=ax.transAxes)
            current_n=self.__neutrons
            current_p=self.__protons
            mafra_out = self.__allmassfractions[0]
            rectangles = []
            edgecolors = []
            zorders = []
            stable_list = []
            linewidths = []
            for i in range(len(current_n)):
                n_tmp = current_n[i]
                p_tmp = current_p[i]
                xy     = n_tmp-0.5,p_tmp-0.5
                is_stable = bool(sum((stable_n==n_tmp) & (stable_p==p_tmp)))

                if mafra_out[i]>= min_X:
                    color_tmp = cmap(norm(mafra_out[i]))
                else:
                    color_tmp = "w"

                if not is_stable:
                    edgecolors.append("gray")
                    zorders.append(1)
                    rect = Rectangle(xy,1,1,zorder=1,lw=0.8)
                    linewidths.append(0.8)
                    stable_list.append(False)
                else:
                    edgecolors.append("k")
                    zorders.append(2)
                    linewidths.append(1.2)
                    stable_list.append(True)

                rect = Rectangle(xy,1,1)
                rectangles.append(rect)
            stable_list=np.array(stable_list)

            self.__helper_stab = stable_list
            # Sort it as the zorder is buggy in patch collections
            # Set the stable once as last so that they are drawn in the foreground
            rectangles =list(np.array(rectangles)[~stable_list])+list(np.array(rectangles)[stable_list])
            edgecolors=list(np.array(edgecolors)[~stable_list])+list(np.array(edgecolors)[stable_list])
            sorted_mafra=np.array(list(np.array(mafra_out)[~stable_list])+list(np.array(mafra_out)[stable_list]))
            linewidths=np.array(list(np.array(linewidths)[~stable_list])+list(np.array(linewidths)[stable_list]))
            # Create the collection
            # In principle, individual patches work as well, but with worse performance
            pc = PatchCollection(rectangles)
            self.__helper_pc = pc
            # Make it white
            cmap.set_under("w")
            pc.set(array=sorted_mafra, cmap=cmap,norm=norm,edgecolors=edgecolors,linewidths=linewidths,zorder=-30)
            ax.add_collection(pc)

            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="3%", pad=0.05)
            final_fig.colorbar(pc,shrink=0.5,label="Mass fraction", cax=cax)

            # Plot magic numbers?
            if plot_magic:
                magicnumbers = [8,20,28,50,82,126]
                #set the magic numbers
                lines = []
                for i in range(len(magicnumbers)):
                    magic = magicnumbers[i]

                    #neutron magic
                    min_z = 10000
                    max_z = -1
                    min_N = -0.5
                    max_N = 10000
                    for ind,p in enumerate(current_p):
                        if current_n[ind] == magic:
                            if current_p[ind]<min_z:
                                min_z=float(current_p[ind])
                            if current_p[ind]>max_z:
                                max_z=float(current_p[ind])

                        if current_p[ind] == magic:
                            if current_n[ind]>min_N:
                                min_N=float(current_n[ind])
                            if current_n[ind]<max_N:
                                max_N=float(current_n[ind])

                    width = 1.3
                    dash  = (3,2)
                    if min_z!=10000:
                        #neutron magic
                        # a,= ax.plot([magic-0.5, magic-0.5], [min_z-0.5, max_z+0.5], color='k',linestyle='--',lw=width, dashes=dash,zorder=100)
                        l = [(magic-0.5, min_z-0.5),( magic-0.5,max_z+0.5)]
                        lines.append(l)
                        # a,= ax.plot([magic+0.5, magic+0.5], [min_z-0.5, max_z+0.5], color='k',linestyle='--',lw=width, dashes=dash,zorder=100)
                        l = [(magic+0.5, min_z-0.5), (magic+0.5, max_z+0.5)]
                        lines.append(l)
                    if max_N!=10000:
                        #proton magic
                        # a,= ax.plot([max_N-0.5,min_N+0.5],[magic-0.5, magic-0.5], color='k',linestyle='--',lw=width, dashes=dash,zorder=100)
                        l = [(max_N-0.5,magic-0.5),(min_N+0.5, magic-0.5)]
                        lines.append(l)
                        # a,= ax.plot([max_N-0.5,min_N+0.5],[magic+0.5, magic+0.5], color='k',linestyle='--',lw=width, dashes=dash,zorder=100)
                        l = [(max_N-0.5,magic+0.5),(min_N+0.5, magic+0.5)]
                        lines.append(l)

                self.__helper_line_segments = LineCollection(lines,
                                               linewidths=1.3,
                                               linestyles='--',
                                               colors="k")
                ax.add_collection(self.__helper_line_segments)

            if element_labels:
                #set the elementlabels
                labels = []
                for z in current_p:
                    min_n = min(current_n[current_p==z])

                    n = nucleus(N=int(min_n),Z=int(z))
                    name  = n.get_elementname()
                    if name == 'neutron':
                        name = 'n'
                    elif name == 'h':
                        name = 'p'
                    else:
                        name = name[0].upper()+name[1:]

                    l = TextPath((min_n - 1.8 , z-0.25), name, size=0.9,ha="right",va='center')
                    labels.append(PathPatch(l))

                pc = PatchCollection(labels,color='k')
                ax.add_collection(pc)

            ax.set_xlabel('Neutron number')
            ax.set_ylabel('Proton number')

            if plot_magic and time_title:
                return self.__helper_pc,self.__helper_line_segments,self.__helper_title,
            if plot_magic and (not time_title):
                return self.__helper_pc,self.__helper_line_segments,
            else:
                return self.__helper_pc,

        def __animate(ind):
            if time_title:
                snaps_time = self.__time[ind]
                self.__helper_title.set_text("t = "+"{:.2e}".format(snaps_time)+" s")

            mafra_out = self.__allmassfractions[ind]
            sorted_mafra=np.array(list(np.array(mafra_out)[~self.__helper_stab])+list(np.array(mafra_out)[self.__helper_stab]))
            self.__helper_pc.set(array=sorted_mafra)

            if plot_magic and time_title:
                return self.__helper_pc,self.__helper_line_segments,self.__helper_title,
            if plot_magic and (not time_title):
                return self.__helper_pc,self.__helper_line_segments,
            else:
                return self.__helper_pc,

        anim = animation.FuncAnimation(final_fig, __animate, init_func=__init_anim,
                               frames=len(self.__allmassfractions), blit=True)

        return anim




    def plot_nuclear_chart_at(self,time,figure=None,fig_is_ax=False,plot_magic=True,
                              colorbar_inset=False,colorbar_position=[0.27, 0.8, 0.5, 0.025],colorbar=True,
                              axes_label=True,time_title=True,min_X=1e-8,max_X=None,cmap= cm.viridis,element_labels=True,
                              nuclei_linewidths=0.8,**kwargs):
        """
            Plot the abundances into the nuclear chart at a given time
        """
        # Make a figure if no figure was given
        if figure == None:
            final_fig = plt.figure(figsize=(10,5))
        else:
            final_fig = figure

        if not fig_is_ax:
            # Get the axis instance from figure
            ax = final_fig.gca()
        else:
            ax = figure
            final_fig = plt.gcf()

        ax.set_aspect(1)

        cmap = cmap.copy() # Avoid annoying deprecation warning

        stable_n,stable_p,stable_a = np.loadtxt(\
                         os.path.dirname(os.path.abspath(\
                         inspect.getfile(inspect.currentframe())))+'/data/stableiso.dat',unpack=True,usecols=[1,2,3])

        # Read time from the snapshots if not done already
        if len(self.__time)==0:
            self.__read_snapshot_time()

        # Get the difference for the time index
        ind  =  np.argmin(np.abs(self.__time - time))

        # corresponding snapshot file
        currentfile = 'snapsh_'+str((ind+1)).zfill(4)+'.dat'
        path        = self.__path+'snaps/'
        snaps_time  = self.__time[ind]
        if time_title:
            ax.set_title("t = "+"{:.2e}".format(snaps_time)+" s")
        current_n,current_p,mafra_out = np.loadtxt(path+currentfile,unpack=True,skiprows=3,usecols=[0,1,3])

        # Set the default upper value
        if max_X is None:
            max_X = np.max(mafra_out)

        # Set a norm
        norm = mpl.colors.LogNorm(vmin=min_X, vmax=max_X,clip=False)

        rectangles = []
        edgecolors = []
        zorders = []
        stable_list = []
        linewidths = []
        for i in range(len(current_n)):
            n_tmp = current_n[i]
            p_tmp = current_p[i]
            xy     = n_tmp-0.5,p_tmp-0.5
            is_stable = bool(sum((stable_n==n_tmp) & (stable_p==p_tmp)))

            if mafra_out[i]>= min_X:
                color_tmp = cmap(norm(mafra_out[i]))
            else:
                color_tmp = "w"

            if not is_stable:
                edgecolors.append("gray")
                zorders.append(1)
                rect = Rectangle(xy,1,1,zorder=1,lw=nuclei_linewidths)
                linewidths.append(0.8)
                stable_list.append(False)
            else:
                edgecolors.append("k")
                zorders.append(2)
                linewidths.append(1.2)
                stable_list.append(True)

            rect = Rectangle(xy,1,1)
            rectangles.append(rect)
        stable_list=np.array(stable_list)

        # Sort it as the zorder is buggy in patch collections
        # Set the stable once as last so that they are drawn in the foreground
        rectangles =list(np.array(rectangles)[~stable_list])+list(np.array(rectangles)[stable_list])
        edgecolors=list(np.array(edgecolors)[~stable_list])+list(np.array(edgecolors)[stable_list])
        sorted_mafra=np.array(list(np.array(mafra_out)[~stable_list])+list(np.array(mafra_out)[stable_list]))
        linewidths=np.array(list(np.array(linewidths)[~stable_list])+list(np.array(linewidths)[stable_list]))
        # Create the collection
        # In principle, individual patches work as well, but with worse performance
        pc = PatchCollection(rectangles)
        # Make it white
        cmap.set_under("w")
        pc.set(array=sorted_mafra, cmap=cmap,norm=norm,edgecolors=edgecolors,linewidths=linewidths)
        ax.add_collection(pc)

        if colorbar:
            if not colorbar_inset:
                divider = make_axes_locatable(ax)
                cax = divider.append_axes("right", size="3%", pad=0.05)
                final_fig.colorbar(pc,shrink=0.5,label="Mass fraction", cax=cax)
            else:
                cax = final_fig.add_axes(colorbar_position)
                final_fig.colorbar(pc,shrink=0.5,label="Mass fraction", cax=cax, orientation='horizontal')
                cax.xaxis.set_label_position('top')
                cax.xaxis.set_ticks_position('top')





        # Plot magic numbers?
        if plot_magic:
            magicnumbers = [8,20,28,50,82,126]
            #set the magic numbers
            for i in range(len(magicnumbers)):
                magic = magicnumbers[i]

                #neutron magic
                min_z = 10000
                max_z = -1
                min_N = -0.5
                max_N = 10000
                for ind,p in enumerate(current_p):
                    if current_n[ind] == magic:
                        if current_p[ind]<min_z:
                            min_z=float(current_p[ind])
                        if current_p[ind]>max_z:
                            max_z=float(current_p[ind])

                    if current_p[ind] == magic:
                        if current_n[ind]>min_N:
                            min_N=float(current_n[ind])
                        if current_n[ind]<max_N:
                            max_N=float(current_n[ind])

                width = 1.3
                dash  = (3,2)
                if min_z!=10000:
                    #neutron magic
                    ax.plot([magic-0.5, magic-0.5], [min_z-0.5, max_z+0.5], color='k',linestyle='--',lw=width, dashes=dash,zorder=100)
                    ax.plot([magic+0.5, magic+0.5], [min_z-0.5, max_z+0.5], color='k',linestyle='--',lw=width, dashes=dash,zorder=100)
                if max_N!=10000:
                    #proton magic
                    ax.plot([max_N-0.5,min_N+0.5],[magic-0.5, magic-0.5], color='k',linestyle='--',lw=width, dashes=dash,zorder=100)
                    ax.plot([max_N-0.5,min_N+0.5],[magic+0.5, magic+0.5], color='k',linestyle='--',lw=width, dashes=dash,zorder=100)

        if element_labels:
            font0 = FontProperties()
            font0.set_weight("ultralight")
            # font0.set_weight("light")
            #set the elementlabels
            labels = []
            for z in current_p:
                min_n = min(current_n[current_p==z])

                n = nucleus(N=int(min_n),Z=int(z))
                name  = n.get_elementname()
                if name == 'neutron':
                    name = 'n'
                elif name == 'h':
                    name = 'p'
                else:
                    name = name[0].upper()+name[1:]

                l = TextPath((min_n - 1.8 , z-0.25), name, size=0.9,ha="right",va='center',prop=font0)
                labels.append(PathPatch(l))

            pc = PatchCollection(labels,color='k')
            ax.add_collection(pc)
        if axes_label:
            ax.set_xlabel('Neutron number')
            ax.set_ylabel('Proton number')
        return final_fig




    def get_A_X_at_time(self,time,A_range=None):
        """
          Return the abundance at a given time.
          Output is:
            Mass number, Mass fraction, time of the nearest snapshot
        """
        plotting_time = time

        if len(self.__time) == 0:
            print('Call read_snapshots first !! The script might not work now!')
            return

        # Convert to numpy array and floats
        self.__time = np.array([float(x) for x in self.__time])

        # Get the difference for the time index
        diff = [abs(x - time) for x in self.__time]
        # Get the closest snapshot
        min_diff = min(diff)
        ind  = list(diff).index(min_diff)

        # Get Neutrons and protons and mass number of the isotopes
        current_n = self.__neutrons
        current_p = self.__protons
        current_a = self.__protons + self.__neutrons

        mafra_out = np.array(self.__allmassfractions)[ind,:]

        # Add up mass fraction for equal mass numbers
        out_x = np.zeros(len(current_a))
        for i in range(len(current_a)):
            out_x[int(current_a[i])] += mafra_out[i]

        if A_range == None:
            # Remove zeros from behind
            for i in range(len(current_a)):
                if out_x[len(current_a)-1-i] >= 1e-20:
                    break
            count = i
        else:
            count = len(current_a)-1-A_range

        return np.arange(len(current_a))[0:len(current_a)-1-count], out_x[0:len(current_a)-1-count], self.__time[ind]  # Return A and X and the taken time


    def plot_A_X_time(self,figure=None,axlabel=True,**kwargs):
        """
          Plot massfractions over mass number and time in a 3D plot. This plot is ugly somehow!
          read_snapshots has to be called first!
        """
        # Make a figure if no figure was given
        if figure == None:
            final_fig = plt.figure()
        else:
            final_fig = figure

        # Get the axis instance from figure
        ax = final_fig.add_subplot(111, projection='3d')

        rang = 200
        A_2d = np.arange(rang)+1.
        X_2d = []
        # Collect data from every snapshot
        for t in self.__time:
            A_2d_tmp, X_2d_tmp, dummy = self.get_A_X_at_time(float(t),A_range=rang)
            X_2d.append(np.clip(np.log([max(x,1e-20) for x in np.array(X_2d_tmp)]), -10, 1, out=np.array(X_2d_tmp)))

        time = [np.log(float(x)) for x in self.__time]

        A_2d, time = np.meshgrid(A_2d, time)

        ax.plot_surface(time, A_2d, X_2d,cstride=1,shade=True)
        #~ ax.plot_wireframe(time, A_2d, X_2d)
        ax.set_ylabel('Mass number A')
        ax.set_xlabel('Log (Time [s])')
        ax.set_zlabel('Log (Mass fraction)')
        ax.set_zlim(-10,0)
        return final_fig
