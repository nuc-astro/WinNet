# Author: M. Reichert
import numpy             as np
import matplotlib.pyplot as plt
import h5py
import os

# This script will plot the different contributions to the energy generation
# from the (cold) CNO-cycles, pp-chain, and Ne-Na, as well as Mg-Al cycles.
# Note that given temperatures should not be smaller than 1d-2GK
# which is the lower limit of the reaclib fits.
# The energy generation of the reactions is not written to the output individually
# and, as a consequence, the plot is only an approximation to some degree.

# Get names of the runs, ignore hidden files
paths = [i for i in os.listdir(".") if i[0] != "."]

# Empty list to store various quantites

# Loop through all folders
for p in paths:
    # Skip non folders
    if not os.path.isdir(p):
        continue
    if p.strip() == "network_data":
        continue
    if not os.path.isfile(os.path.join(p,"WinNet_data.h5")):
        raise Exception("Could not find WinNet_data.h5 in '"+str(p)+"'. "\
                        "Are you sure you enabled hdf5 output in the makefile "\
                        "before running?")

    # Open the file
    f = h5py.File(os.path.join(p,"WinNet_data.h5"),"r")

    # A bit handwaving formula to get the energy plateau at the end.
    # This plateau can be estimated in a more physical motivated way rather
    # then fitting it by hand as done here. To estimate it, one can
    # consider when the pariticipating nuclei are in equilibrium and take
    # this time. A detailed derivation is given in
    # e.g., Nuclear physics of stars, C. Iliadis, Sect. 5.1.
    time_plateau = 10**(-4/2*float(p)+17.5)

    # Total energy
    time = f["energy/time"][:]
    idx_plateau = np.argmin(abs(time-time_plateau))

    # Get nuclear properties
    A = f["energy/A"][:]
    Z = f["energy/Z"][:]
    N = f["energy/N"][:]

    # In the following we separate the energy generation from
    # the different cycles. The reactions and cycles can e.g.,
    # be seen at https://cococubed.com/code_pages/burn_hydrogen.shtml .

    # CNO-cycles and pp-chain
    idx_p   = (A == 1) & (Z ==1)
    idx_d   = (A == 2) & (Z ==1)
    idx_he3 = (A == 3) & (Z ==2)
    idx_he4 = (A == 4) & (Z ==2)
    idx_be7 = (A == 7) & (Z ==4)
    idx_be8 = (A == 8) & (Z ==4)
    idx_b8  = (A == 8) & (Z ==5)
    idx_li7 = (A == 7) & (Z ==3)
    idx_c12 = (A == 12) & (Z ==6)
    idx_c13 = (A == 13) & (Z ==6)
    idx_n13 = (A == 13) & (Z ==7)
    idx_n14 = (A == 14) & (Z ==7)
    idx_n15 = (A == 15) & (Z ==7)
    idx_o15 = (A == 15) & (Z ==8)
    idx_o16 = (A == 16) & (Z ==8)
    idx_o17 = (A == 17) & (Z ==8)
    idx_o18 = (A == 18) & (Z ==8)
    idx_f17 = (A == 17) & (Z ==9)
    idx_f18 = (A == 18) & (Z ==9)
    idx_f19 = (A == 19) & (Z ==9)

    e_detailed_pg    = f["energy/detailed (p,g)"][idx_plateau,:]
    e_detailed_bet   = f["energy/detailed decay"][idx_plateau,:]
    e_detailed_ap    = f["energy/detailed (a,p)"][idx_plateau,:]
    e_detailed_other = f["energy/detailed other"][idx_plateau,:]
    e_detailed_ag    = f["energy/detailed (a,g)"][idx_plateau,:]

    #############
    # CNO cycle #
    #############

    # CNO cycle 1
    # (p,g) reactions
    c13pg = e_detailed_pg[idx_c13]
    n14pg = e_detailed_pg[idx_n14]
    c12pg = e_detailed_pg[idx_c12]
    # Decays
    n13d = e_detailed_bet[idx_n13]
    o15d = e_detailed_bet[idx_o15]
    # (p,a) reactions
    n15pa = -e_detailed_ap[idx_n15]
    # Total CNO 1
    e_cno = c13pg+n14pg+c12pg+n13d+o15d+n15pa

    # CNO cycle 2
    # (p,g) reactions
    n15pg = e_detailed_pg[idx_n15]
    o16pg = e_detailed_pg[idx_o16]
    # Decays
    f17d = e_detailed_bet[idx_f17]
    # (p,a) reactions
    o17pa = -e_detailed_ap[idx_o17]
    # Total CNO 1 + CNO 2
    e_cno = e_cno+n15pg+o16pg+f17d+o17pa

    # CNO cycle 3
    # (p,g) reactions
    o17pg = e_detailed_pg[idx_o17]
    # Decays
    f18pg = e_detailed_bet[idx_f17]
    # (p,a) reactions
    o18pa = -e_detailed_ap[idx_o18]
    # Total CNO 1 + CNO 2 + CNO 3
    e_cno = e_cno+o17pg+f18pg+o18pa

    # CNO cycle 4
    # (p,g) reactions
    o18pg = e_detailed_pg[idx_o18]
    # (p,a) reactions
    f19pa = -e_detailed_ap[idx_f19]
    # Total CNO 1 + CNO 2 + CNO 3
    e_cno = e_cno+o18pg+f19pa

    #############
    # PP-Chains #
    #############
    ppg = e_detailed_pg[idx_p]
    dpg = e_detailed_pg[idx_d]
    pd  = e_detailed_bet[idx_p]
    # Branch 1
    he3he3 = e_detailed_other[idx_he3]
    # Branch 2 and 3
    he3a   = e_detailed_ag[idx_he3]
    be7d   = e_detailed_bet[idx_be7]
    li7p   = e_detailed_other[idx_li7]
    be7pg  = e_detailed_pg[idx_be7]
    b8d    = e_detailed_bet[idx_b8]
    be8o   = e_detailed_other[idx_be8]

    # Total energy generation of pp-chain
    e_pp = ppg+dpg+he3he3+pd+he3a+li7p+be7d+b8d+be8o

    #########################
    # Ne-Na and Mg-Al chain #
    #########################
    idx_ne20 = (A == 20) & (Z ==10)
    idx_ne21 = (A == 21) & (Z ==10)
    idx_ne22 = (A == 22) & (Z ==10)
    idx_na21 = (A == 21) & (Z ==11)
    idx_na22 = (A == 22) & (Z ==11)
    idx_na23 = (A == 23) & (Z ==11)
    idx_mg22 = (A == 22) & (Z ==12)
    idx_mg23 = (A == 23) & (Z ==12)
    idx_mg24 = (A == 24) & (Z ==12)
    idx_mg25 = (A == 25) & (Z ==12)
    idx_mg26 = (A == 26) & (Z ==12)
    idx_al25 = (A == 25) & (Z ==13)
    idx_al26 = (A == 26) & (Z ==13)
    idx_al27 = (A == 27) & (Z ==13)
    idx_si26 = (A == 26) & (Z ==14)
    idx_si27 = (A == 27) & (Z ==14)

    # Ne-Na chain
    #(p,g)
    ne20pg = e_detailed_pg[idx_ne20]
    ne21pg = e_detailed_pg[idx_ne21]
    na21pg = e_detailed_pg[idx_na21]
    na22pg = e_detailed_pg[idx_na22]
    ne22pg = e_detailed_pg[idx_ne22]
    #decay
    na21d = e_detailed_bet[idx_na21]
    mg22d = e_detailed_bet[idx_mg22]
    na22d = e_detailed_bet[idx_na22]
    mg23d = e_detailed_bet[idx_mg23]
    #(p,a)
    na23pa = -e_detailed_ap[idx_na23]

    e_nena = ne20pg+ne21pg+na21pg+na22pg+ne22pg+na21d+mg22d+na22d+mg23d+na23pa

    # Mg-al chain
    #(p,g)
    na23pg = e_detailed_pg[idx_na23]
    mg24pg = e_detailed_pg[idx_mg24]
    al25pg = e_detailed_pg[idx_al25]
    mg25pg = e_detailed_pg[idx_mg25]
    al26pg = np.sum(e_detailed_pg[idx_al26])
    mg26pg = e_detailed_pg[idx_mg26]
    #decay
    al25d = e_detailed_bet[idx_al25]
    si26d = e_detailed_bet[idx_si26]
    al26d = np.sum(e_detailed_bet[idx_al26])
    si27d = e_detailed_bet[idx_si27]
    #(p,a)
    al27pa = -e_detailed_ap[idx_al27]

    e_mgal = na23pg+mg24pg+al25pg+mg25pg+al26pg+mg26pg+al25d+si26d+al26d+si27d+al27pa

    # Plot Ne-Na cycle
    plt.scatter(float(p),e_nena,color="tab:orange")
    # Plot Mg-Al cycle
    plt.scatter(float(p),e_mgal,color="tab:green")
    # Plot PP-Chains
    plt.scatter(float(p),e_pp,color="tab:blue")
    # Plot CNO cycles
    plt.scatter(float(p),e_cno,color="tab:red")


# Indicate the suns core temperature
plt.axvline(1.56,lw=2,ls="--",color="lightgrey")
plt.text(1.56,5e6,r"Temperature of the sun",rotation=90,va="top",ha="right",fontsize=14,color="lightgrey")

# Create labels for the legend
plt.scatter(np.nan,np.nan,color="tab:red"   ,label="CNO-cycle")
plt.scatter(np.nan,np.nan,color="tab:blue"  ,label="pp-chain")
plt.scatter(np.nan,np.nan,color="tab:orange",label="Ne-Na cycle")
plt.scatter(np.nan,np.nan,color="tab:green",label="Mg-Al cycle")
plt.legend()

# Set the limits of the plot
plt.ylim(3e0,1e7)
plt.xlim(1,4)
plt.yscale("log")
plt.xlabel(r"Temperature [10$^7$ K]")
plt.ylabel(r"Energy [erg g$^{-1}$ s$^{-1}$]")
plt.title(r"$\rho$ = 100 g cm$^{-3}$")
plt.savefig("energy_generation.pdf",bbox_inches="tight")
plt.show()
