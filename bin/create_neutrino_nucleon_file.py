# Author: M. Reichert
# Date  : 11.04.2023
import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as const
from scipy.integrate import quad


class nunucleon(object):

    def __init__(self):
        """
          Class to calculate the charged particle neutrino-nucleon cross sections
          and average cross sections.
        """

        # Define constants
        self.MeV2erg = 1.60218e-6     # Conversion MeV -> Erg
        self.mec2    = 0.510998910e0  # MeV
        self.mnc2    = 939.56533e0    # MeV
        self.mpc2    = 938.271998e0   # MeV
        self.Gf      = 1.436e-49      # erg/ccm
        self.c       = const.c*100    # cm/s
        self.hbar    = const.hbar*1e7 # erg*s
        self.ga      = 1.26           # Axial vector coupling constant

        # Energy difference between neutron and proton
        self.Dnp = self.mnc2-self.mpc2

        # Define constants as in Horowitz 2002
        self.cv = 1
        self.ca = self.ga
        self.F2 = 3.706

        # Calculate sigma 0. Equation 9 of Burrows 2006
        self.sigma_0 = 4*self.Gf**2*(self.mec2*self.MeV2erg)**2/((self.c*self.hbar)**4*np.pi)

        # Define the temperature Grid that is used in WinNet
        self.T_grid = np.array([2.8e0,3.5e0,4.0e0,5.0e0,6.4e0,8.0e0,10.0e0])


    def __cor(self,E,nu):
        """
          Correction for weak magnetism and recoil according to Horowitz 2002
          (https://ui.adsabs.harvard.edu/abs/2002PhRvD..65d3001H/abstract, Equation 22).
        """
        # Make the constant writing a bit shorter
        cv = self.cv
        ca = self.ca
        F2 = self.F2

        # Express energy in units of mnc2
        E =  E/self.mnc2

        # Equation 22 in Horowitz 2002
        R = cv**2*(1+4*E+16/3.*E**2)+3*ca**2*(1+4/3*E)**2
        # Difference between neutrino and antineutrino
        if nu:
            R = R + 4*(cv+F2)*ca*E*(1+4/3*E)
        else:
            R = R - 4*(cv+F2)*ca*E*(1+4/3*E)
        R = R + 8/3*cv*F2*E**2
        R = R+5/3*E**2*(1+2/5*E)*F2**2
        # Calculate denominator
        denom = ((cv**2+3*ca**2)*(1+2*E)**3)

        R = R/denom

        return R

    def sigma_nu_n(self,E):
        """
          Calculate the electron neutrino cross section according to Equation 10 of
          Burrows et al. 2006 (https://ui.adsabs.harvard.edu/abs/2006NuPhA.777..356B/abstract).
          Outputs the cross section in 10^42 cm^2.
        """
        res = self.sigma_0*1e42*((1+3*self.ga**2)/4.)*((E+self.Dnp)/self.mec2)**2
        res = res * (1-(self.mec2/(E+self.Dnp))**2.0)**(0.5)*self.WM(E,0)
        return res

    def WM(self,E,mode=0):
        """
          Weak magnetism and recoil corrections.
          Mode 0 returns a higher order of this correction according to Horowitz 2002,
          mode 1 returns a simplified value of Burrows et al. 2006.
        """
        if mode==0:
            res = self.__cor(E,True)
        elif mode==1:
            res = 1+1.1*E/self.mnc2
        else:
            raise ValueError("Mode not defined!")

        return res


    def sigma_anu_p(self,E):
        """
          Calculate the electron anti neutrino cross section according to Equation 11 of
          Burrows et al. 2006 (https://ui.adsabs.harvard.edu/abs/2006NuPhA.777..356B/abstract).
          Outputs the cross section in 10^42 cm^2.
        """
        res = self.sigma_0*1e42*((1+3*self.ga**2)/4.)*((E-self.Dnp)/self.mec2)**2
        res = res * (1-(self.mec2/(E-self.Dnp))**2.0)**(0.5)*self.WMbar(E,0)

        return res


    def WMbar(self,E,mode=0):
        """
         Weak magnetism and recoil corrections.
         Mode 0 returns a higher order of this correction according to Horowitz 2002,
         mode 1 returns a simplified value of Burrows et al. 2006.
        """
        if mode==0:
            res = self.__cor(E,False)
        elif mode==1:
            res = 1-7.1*E/self.mnc2
        else:
            raise ValueError("Mode not defined!")

        return res


    def Fermi_Dirac(self,E,T):
        """
          Calculate the Fermi-Dirac distribution.
        """
        res = E**2/(1+np.exp(E/T))
        return res


    def sigma_avg_nu_n(self,T):
        """
          Calculate the average neutrino cross section:
            nu_e + n -> e + p
          The average neutrino cross section is calculated according
          to the integral of the cross sections multiplied by the normalized
          Fermi-Dirac distribution.
        """
        integrand = lambda E: self.sigma_nu_n(E) * self.Fermi_Dirac(E, T)
        result, _ = quad(integrand, 0, np.inf)

        # Integrate the Fermi-Dirac distribution for the normalization
        integrand = lambda E: self.Fermi_Dirac(E, T)
        result2, _ = quad(integrand, 0, np.inf)

        return result/result2

    def averageE_nu_n(self,T):
        """
          Calculate the average energy of the absorped neutrino for the reaction:
            nu_e + n -> e + p
        """
        integrand = lambda E: self.sigma_nu_n(E) * self.Fermi_Dirac(E, T)*E
        result, _ = quad(integrand, 0, np.inf)

        integrand = lambda E: self.Fermi_Dirac(E, T)
        result2, _ = quad(integrand, 0, np.inf)

        result3 = self.sigma_avg_nu_n(T)
        return (result/result2)/result3


    def averageE_anu_p(self,T):
        """
          Calculate the average energy of the absorped neutrino for the reaction:
            bar(nu)_e + p -> n + positron
        """
        integrand = lambda E: self.sigma_anu_p(E) * self.Fermi_Dirac(E, T)*E
        result, _ = quad(integrand, self.Dnp+self.mec2, np.inf)

        integrand = lambda E: self.Fermi_Dirac(E, T)
        result2, _ = quad(integrand, self.Dnp+self.mec2, np.inf)

        result3 = self.sigma_avg_anu_p(T)
        return (result/result2)/result3


    def sigma_avg_anu_p(self,T):
        """
          Calculate the average neutrino cross section:
            bar(nu)_e + p -> n + positron
          The average neutrino cross section is calculated according
          to the integral of the cross sections multiplied by the normalized
          Fermi-Dirac distribution.
        """
        integrand = lambda E: self.sigma_anu_p(E) * self.Fermi_Dirac(E, T)
        result, _ = quad(integrand, self.Dnp+self.mec2, np.inf)

        # Integrate the Fermi-Dirac distribution for the normalization
        integrand = lambda E: self.Fermi_Dirac(E, T)
        result2, _ = quad(integrand, self.Dnp+self.mec2, np.inf)
        return result/result2


    def create_WinNet_file_cross_section(self,file_name):
        """
          Create the WinNet file with the calculated cross sections
          and the average energy of the absorped neutrino.
        """
        # The temperature grid that is used in WinNet
        Tgrid = self.T_grid

        # Create a line in the correct format, note that also
        # the source label is important here.
        out = "         n    p                             nen"+"\n"
        # Calculate and write the neutrino cross sections
        line = ""
        for ind,T in enumerate(Tgrid):
            line += "{:6.2f}".format(self.sigma_avg_nu_n(T))+"    "
        line = line.rstrip()
        out += line+"\n"
        # Calculate and write the average energy of the absorped neutrino
        line = ""
        for ind,T in enumerate(Tgrid):
            line += "{:6.2f}".format(self.averageE_nu_n(T))+"    "
        line = line.rstrip()
        out += line+"\n"

        out += "         p    n                            nebp"+"\n"
        line = ""
        for ind,T in enumerate(Tgrid):
            line += "{:6.2f}".format(self.sigma_avg_anu_p(T))+"    "
        line = line.rstrip()
        out += line+"\n"
        # Calculate and write the average energy of the absorped neutrino
        line = ""
        for ind,T in enumerate(Tgrid):
            line += "{:6.2f}".format(self.averageE_anu_p(T))+"    "
        line = line.rstrip()
        out += line+"\n"

        # Save the file
        with open(file_name,"w") as f:
            f.write(out)



if __name__ == "__main__":
    """
      Create the cross section file when executing this script.
      Furthermore, plot the cross sections for a range of temperatures.
    """
    # File name to save the cross sections and average energies of absorbed neutrinos
    filename_crosssection     = "neunucleons.dat"
    # Create an instance of the class
    n_class = nunucleon()
    n_class.create_WinNet_file_cross_section(filename_crosssection)

    # Prepare everything for plotting
    fig = plt.figure(figsize=(5,3))

    # Create a range of neutrino temperatures
    # and calculate the cross sections
    n_temp = 100
    T_more = np.linspace(2,13.0,n_temp)
    nu     = np.zeros(n_temp)
    anu    = np.zeros(n_temp)
    for ind,T in enumerate(T_more):
        nu[ind]  = n_class.sigma_avg_nu_n(T)
        anu[ind] = n_class.sigma_avg_anu_p(T)

    # Plot them
    plt.plot(T_more,nu,label=r"$\nu_e n$")
    plt.plot(T_more,anu,label=r"$\bar{\nu}_e p$")

    # Now the same for the Grid that is used in WinNet
    Tgrid = np.array([2.8e0,3.5e0,4.0e0,5.0e0,6.4e0,8.0e0,10.0e0])
    nu     = np.zeros(len(Tgrid))
    anu    = np.zeros(len(Tgrid))
    for ind,T in enumerate(Tgrid):
        nu[ind]  = n_class.sigma_avg_nu_n(T)
        anu[ind] = n_class.sigma_avg_anu_p(T)

    # Also plot this
    plt.scatter(Tgrid,nu,facecolors='none',edgecolors='k')
    plt.scatter(Tgrid,anu,facecolors='none',edgecolors='k')

    # Make the axis more pretty
    plt.xlabel("T [MeV]")
    plt.ylabel(r"$<\sigma> (10^{-42}$ cm$^2)$")
    plt.legend()
    plt.tight_layout()
    # Show the plot
    plt.show()
