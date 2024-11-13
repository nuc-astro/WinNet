import numpy as np
import pandas as pd
from winvn_class import winvn
from numba import jit


class ngamma_eq(object):

    def __init__(self, path_winvn="winvne_v2.0.dat"):
        """
        Constructor for the class
        """
        # Set the paths
        self.__path_winvn = path_winvn

        # Read the files
        self.__read_winvn()
        self.__calc_Sn_winvn()

        # Define constants
        self.__kb = 8.617343e-11  # MeV/K
        self.__mu = 1.660539e-24  # g
        self.__hix = 1.036427e-18
        self.__NA = 6.02214076e23  # 1/mol
        self.__hbar = 6.582119e-22  # MeV s



    def __read_winvn(self):
        """
            Read the file containing the partition functions
        """
        w = winvn(self.__path_winvn)
        w.read_winvn()
        df = w.get_dataframe()
        pf   = df["function"].values
        Z    = df["Z"].values
        N    = df["N"].values
        spin = df["spin"].values

        # Create 2D array with Z and N, fill non existing with nans
        # pf is of type function
        self.__pf_Z_N = np.zeros((Z.max()+1,N.max()+1),dtype=object)
        self.__spin_Z_N = np.zeros((Z.max()+1,N.max()+1))
        self.__pf_Z_N[:] = lambda x: np.nan
        self.__pf_Z_N[Z,N] = pf

        self.__spin_Z_N[Z,N] = spin


    def __calc_Sn_winvn(self):
        """
          Calculate the neutron separation energy from the mass excess in the winvn
        """
        w = winvn(self.__path_winvn)
        w.read_winvn()
        df = w.get_dataframe()

        mexc_n = df.loc["n","mass excess"]

        # Put mass excesses in 2D array
        Zwinvn = df["Z"].values
        Nwinvn = df["N"].values
        mexc   = df["mass excess"].values

        # Create 2D array with Z and N, fill non existing with nans
        self.__Sn_Z_N_winvn = np.zeros((Zwinvn.max()+1,Nwinvn.max()+1))
        self.__Sn_Z_N_winvn[:] = np.nan
        self.__Sn_Z_N_winvn[Zwinvn,Nwinvn] = mexc

        # Now calculate the neutron separation energy
        self.__Sn_Z_N_winvn[:,1:] = (self.__Sn_Z_N_winvn[:,0:-1] + mexc_n) - self.__Sn_Z_N_winvn[:,1:]

        # Put this into the Sn array
        self.__Sn_Z_N = self.__Sn_Z_N_winvn
        # ALso A Z and N
        self.__A = df["A"].values
        self.__Z = df["Z"].values
        self.__N = df["N"].values



    @staticmethod
    @jit(nopython=True)
    def helper_calc(Sn_array,density,temperature,yn,pf,spin):
        """
            Helper function to calculate the path of the r-process. Necessary to use numba.
        """

        def helper_calc_ratio(Z, N, T, ndens, Sn_Z_N, pf_Z_N, spin_Z_N):
            # Define constants
            kb = 8.617343e-11  # MeV/K
            hix = 1.036427e-18
            NA = 6.02214076e23  # 1/mol
            hbar = 6.582119e-22  # MeV s

            # kbT in MeV
            kbT = kb * T * 1e9
            A = Z + N
            Snp1 = Sn_Z_N[Z, N + 1]
            Sn = Sn_Z_N[Z, N]
            pf = pf_Z_N[Z, N]
            pf1 = pf_Z_N[Z, N + 1]

            sp = 2 * spin_Z_N[Z, N] + 1
            sp1 = 2 * spin_Z_N[Z, N + 1] + 1

            a = np.log((A + 1) / A) * 3 / 2
            b = (Snp1) / (kbT)
            c = np.log(pf * sp / (2 * pf1 * sp1))

            d = np.log(2 * np.pi * hbar**2 / (hix * kbT)) * 1.5
            e = np.log(ndens * NA)

            ratio = a + b + c + d + e
            return ratio

        ndens = density*yn

        # Make zerolike self.__S2n_Z_N_winvn
        Zmax = Sn_array.shape[0]-1
        Nmax = Sn_array.shape[1]-1
        logratios = np.zeros((Zmax+1, Nmax+1))*np.nan

        path_Z = []
        path_N = []
        for Z in range(0,Zmax):
            for N in range(0,Nmax):
                if np.isnan(Sn_array[Z,N]):
                    continue

                logratios[Z,N] = helper_calc_ratio(Z,N,temperature,ndens,Sn_array,pf,spin)

            # Find the first occurence of ratio <= 1
            if np.all(np.isnan(logratios[Z,:])):
                continue
            if Z<=2:
                continue

            log_ratios = logratios[Z, :]
            log_cumprod = np.nancumsum(log_ratios)  # Cumulative sum of logarithmic values

            maxval = np.where(~np.isnan(Sn_array[Z,:]) & (Sn_array[Z,:]>0))[0][-1]
            minval = np.where(~np.isnan(Sn_array[Z,:]))[0][0]

            idx = np.minimum(np.maximum(np.argmax(log_cumprod)+1,minval),maxval)

            if len(path_Z) >= 1:
                if path_Z[-1] != Z:
                    path_Z.append(Z)
                    path_N.append(path_N[-1]-1)

            path_Z.append(Z)
            path_N.append(idx)

        return path_Z, path_N



    def calc_r_process_path(self,density,temperature,yn):
        """
            Get the path of the r-process
        """

        # Precompute the shape once
        rows, cols = self.__pf_Z_N.shape

        # Use np.empty to pre-allocate the pf array
        pf = np.empty((rows, cols))
        def apply_callable(Z, N):
            func = self.__pf_Z_N[Z, N]
            if callable(func):
                return func(temperature)
            return np.nan  # or some default value if not callable

        # Apply the function across the array using np.vectorize
        vectorized_calc = np.vectorize(apply_callable)
        pf = vectorized_calc(np.arange(rows)[:, None], np.arange(cols))

        self.path_Z, self.path_N = self.helper_calc(self.__Sn_Z_N,density,temperature,yn,pf,self.__spin_Z_N)


    # @jit(nopython=True)
    def __calc_ratio(self,Z,N,T,ndens):
        # kbT in MeV
        kbT  = self.__kb*T*1e9
        A    = Z+N
        Snp1 = self.__Sn_Z_N[Z,N+1]
        Sn   = self.__Sn_Z_N[Z,N]
        pf   = self.__pf_Z_N[Z,N](T)
        pf1  = self.__pf_Z_N[Z,N+1](T)

        sp = 2*self.__spin_Z_N[Z,N]+1
        sp1 = 2*self.__spin_Z_N[Z,N+1]+1

        a = np.log((A+1)/A)*3/2
        b = (Snp1)/(kbT)
        # with warnings.catch_warnings():
            # warnings.simplefilter("ignore", RuntimeWarning)
        c = np.log(pf*sp/(2*pf1*sp1))

        d = np.log(2*np.pi*self.__hbar**2 / (self.__hix*kbT)) * 1.5
        e = np.log(ndens*self.__NA)

        # with warnings.catch_warnings():
            # warnings.simplefilter("ignore", RuntimeWarning)
        ratio = np.exp(a+b+c+d+e,dtype=np.float128)

        return ratio

if __name__ == "__main__":

    # Create the object
    ng = ngamma_eq("/home/mreichert/data/Networks/comparison_winNet/WinNet-dev/data/winvne_v2.0.dat")
    ng.calc_r_process_path(1e6,4,5e-1)
