from __future__                 import print_function
from scipy.interpolate          import interp1d
import numpy                    as np
import math
import os
import sys
import inspect

class nucleus(object):
    """
    nucleus: contains basic information of a nuclei. For example the name, the name in the network, the amount of protons, neutrons and the mass number
    """

  # constructor
    def __init__(self,name='',Z=-1,N=-1,Y=0.,quiet=False,init_prop=False):
        """
        Input:
          name       - name of the nuclei. All informations are extracted from this name.
        """
        self.__input_name = name
        self.__name  = name.lower()
        self.__A     = -1
        self.__Z     = -1
        self.__N     = -1
        self.__warn  = False
        self.__quiet = quiet
        self.__ab    = Y
        self.__is_stable = False
        #Sorting criteria
        self.__sort_Z = True

        #Special cases
        if self.__name == 'p':
            self.__name = 'h1'
        if self.__name == 'n' or self.__name == 'neutrons':
            self.__name = 'neutron'
        if self.__name == 'd':
            self.__name = 'h2'
        if self.__name == 't':
            self.__name = 'h3'
        if self.__name == 'al-6':
            self.__name = 'al26'
        if self.__name == 'al*6':
            self.__name = 'al26'

        self.__elementname = ''
        self.__elname      = ('neutron','h','he','li','be','b','c','n','o','f','ne','na','mg','al','si','p','s','cl','ar','k','ca','sc','ti','v','cr','mn','fe',
          'co','ni','cu','zn','ga','ge','as','se','br','kr','rb','sr','y','zr','nb','mo','tc','ru','rh','pd','ag','cd','in','sn','sb',
          'te', 'i','xe','cs','ba','la','ce','pr','nd','pm','sm','eu','gd','tb','dy','ho','er','tm','yb','lu','hf','ta','w','re','os',
          'ir','pt','au','hg','tl','pb','bi','po','at','rn','fr','ra','ac','th','pa','u','np','pu','am','cm','bk','cf','es','fm','md',
          'no','lr','rf','db','sg','bh','hs','mt','ds','rg','cn','nh','fl','mc','lv','ts','og',"uue","ubn","ubu","ubb","ubt","ubq","ubp",
          "ubh","ubs","ubo","ube","utn","utu","utb","utt","utq","utp","uth","uts","uto","ute","uqn","uqu")


        if (Z!=-1 and N!=-1):
            if ((Z >= len(self.__elname)) and not self.__quiet):
                print('Proton number exceeds the one of database ',Z)
            self.__A = Z + N
            self.__name = self.__elname[Z] + str(self.__A)

        self.__initialize()

        if init_prop:
            self.__init_properties()

        # Check if an error occured
        self.__not_a_isotope = (self.__Z == -1)



    def __init_properties(self):
        """
          Initialize basic properties, like stableness etc.
        """
        # Check if it is stable
        stable_nuc = np.loadtxt(os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))+'/data/stableiso.dat',unpack=True,usecols=[0],dtype=str)
        if self.__name in stable_nuc:
            self.__is_stable = True
        else:
            self.__is_stable = False



    def is_stable(self):
        """
          Is it a stable nucleus?
        """
        return self.__is_stable


    def get_not_a_isotope(self):
        """
        Check if the input was a valid isotope. True if not, False if it is valid
        """
        return self.__not_a_isotope

    def __initialize(self):
        """
        Extract the information of protonnumber neutronnumber out of the name.
        """
        self.__A  = ''.join(filter(lambda x: x.isdigit(),self.__name))

        if (self.__A != ''):
            self.__A = int(self.__A)
        else:
            self.__A = 1

        self.__elementname = ''.join(filter(lambda x: not x.isdigit(),self.__name)).strip()

        if self.__elementname in self.__elname:
            self.__Z = self.__elname.index(self.__elementname)   #protonnumber
            if self.__elementname.lower()=="n" and self.__A==1:
                self.__Z = 0
        else:
            self.__Z = -1
            if not self.__quiet:
                print('Warning: The element '+self.__elementname+' is not known by nuclei class.')
            self.__warn = True

        self.__N = self.__A - self.__Z
        self.__nrnname='i'+self.__name

        self.__N = self.__A - self.__Z

  #~ Getters
    def get_A(self):
        """
        get_A : Get the mass number of the nucleus
        """
        return self.__A

    def get_elnames(self):
        """
        get_elnames : Get elementnames
        """
        return self.__elname

    def get_Z(self):
        """
        get_Z : Get the proton number of the nucleus
        """
        return self.__Z

    def get_N(self):
        """
        get_N : Get the neutron number of the nucleus
        """
        return self.__N

    def get_Y(self):
        """
        get_Y : Get the abundance of the nucleus
        """
        return self.__ab

    def get_X(self):
        """
        get_X : Get the mass fraction of the nucleus
        """
        return self.__ab * self.__A

    def get_elementname(self):
        """
        get_elementname : Get the corresponding name of the element (as a string)
        """
        return self.__elementname

    def get_name(self):
        """
        get_name : Get the full name of the nucleus. The name is the elementname plus the massnumber as a string
        """
        if self.__name=="neutron1":
            outp="n"
        elif self.__name=="h1":
            outp='p'
        elif self.__name=="h2":
            outp='d'
        elif self.__name=="h3":
            outp='t'
        else:
            outp=self.__name
        return outp

    def get_input_name(self):
        """
        get_input_name : Get the name as it was inputted
        """
        return self.__input_name




    def get_seedline(self):
        """
        get_seedline : Get the line for the nucleus as it is necessary for the seed file
        """
        mafra = '%1.6e' % float(self.get_X())
        nam = self.get_name()
        if nam.strip() == 'neutron':
            nam = 'n'
        elif nam.strip() == 'h1':
            nam = 'p'
        elif nam.strip() == 'h2':
            nam = 'd'
        elif nam.strip() == 'h3':
            nam = 't'

        line = nam.rjust(5) + ' '*7 + mafra
        return line

    def set_X(self,X):
        """
        set_X : set the massfraction of the nucleus
        """
        self.__ab = X / self.__A

    def set_Y(self,Y):
        """
        set_X : set the massfraction of the nucleus
        """
        self.__ab = Y

    def set_sortcriteria(self,criteria):
        """
        set_sortcriteria : Set the criteria for sorting (possible values are "A" and "Z")
        """
        if criteria == 'Z':
            self.__sort_Z = True
        if criteria == 'A':
            self.__sort_Z = False

    def __gt__(self,other):
        """
        A nucleus is greater if it has an higher proton number
        """
        if self.__sort_Z == True:   #Sort Z
            if (self.__Z > other.get_Z()):
                return True
            elif (self.__Z == other.get_Z()):
                if (self.__A >= other.get_A()):
                    return True
                else:
                    return False
            else:
                return False
        else:                       #Sort A
            if (self.__A > other.get_A()):
                return True
            elif (self.__A == other.get_A()):
                if (self.__Z >= other.get_Z()):
                    return True
                else:
                    return False
            else:
                return False
