# Author: M. Reichert
# Date  : 15.03.2023

# Creates a file with alpha decay half-lives. The half-lifes
# are calculated by the Viola-Seaborg formula. The fit parameters
# are partly taken from Dong&Ren 2005 and have been partly fitted
# to experimental available data of the reaclib.
# To create this file, a winvn file is needed.
import numpy as np
from tqdm import tqdm


# Path to a valid winvn file
path_winvn = "winvne_v2.0.dat"


def is_stable_nucleus(Z,N):
    """
      Check if the nucleus is a stable one.
    """
    Nstable = [ 0,   1,   1,   2,   3,   4,   5,   5,   6,   6,   7,   7,   8,   8,   9,  10,  10,  10,
               11,  12,  12,  12,  13,  14,  14,  14,  15,  16,  16,  16,  17,  18,  20,  18,  20,  18,
               20,  22,  20,  22,  20,  22,  23,  24,  26,  24,  24,  25,  26,  27,  28,  28,  26,  28,
               29,  30,  30,  28,  30,  31,  32,  32,  30,  32,  33,  34,  36,  34,  36,  34,  36,  37,
               38,  40,  38,  40,  38,  40,  41,  42,  42,  40,  42,  43,  44,  46,  44,  46,  42,  44,
               46,  47,  48,  50,  48,  46,  48,  49,  50,  50,  50,  51,  52,  54,  52,  50,  52,  53,
               54,  55,  56,  52,  54,  55,  56,  57,  58,  60,  58,  56,  58,  59,  60,  62,  64,  60,
               62,  58,  60,  62,  63,  64,  66,  64,  62,  64,  65,  66,  67,  68,  69,  70,  72,  74,
               70,  72,  68,  70,  71,  72,  73,  74,  74,  70,  72,  74,  75,  76,  77,  78,  80,  78,
               74,  76,  78,  79,  80,  81,  82,  82,  78,  80,  82,  84,  82,  82,  83,  85,  86,  88,
               82,  87,  88,  90,  92,  90,  90,  91,  92,  93,  94,  96,  94,  90,  92,  94,  95,  96,
               97,  98,  98,  94,  96,  98,  99, 100, 102, 100,  98, 100, 101, 102, 103, 104, 106, 104,
              104, 105, 106, 107, 108, 108, 108, 109, 110, 112, 110, 111, 112, 113, 114, 116, 114, 116,
              114, 116, 117, 118, 120, 118, 116, 118, 119, 120, 121, 122, 124, 122, 124, 122, 124, 125,
              126]

    Zstable = [ 1,  1,  2,  2,  3,  3,  4,  5,  5,  6,  6,  7,  7,  8,  8,  8,  9, 10, 10, 10, 11, 12, 12, 12,
               13, 14, 14, 14, 15, 16, 16, 16, 16, 17, 17, 18, 18, 18, 19, 19, 20, 20, 20, 20, 20, 21, 22, 22,
               22, 22, 22, 23, 24, 24, 24, 24, 25, 26, 26, 26, 26, 27, 28, 28, 28, 28, 28, 29, 29, 30, 30, 30,
               30, 30, 31, 31, 32, 32, 32, 32, 33, 34, 34, 34, 34, 34, 35, 35, 36, 36, 36, 36, 36, 36, 37, 38,
               38, 38, 38, 39, 40, 40, 40, 40, 41, 42, 42, 42, 42, 42, 42, 44, 44, 44, 44, 44, 44, 44, 45, 46,
               46, 46, 46, 46, 46, 47, 47, 48, 48, 48, 48, 48, 48, 49, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50,
               51, 51, 52, 52, 52, 52, 52, 52, 53, 54, 54, 54, 54, 54, 54, 54, 54, 55, 56, 56, 56, 56, 56, 56,
               56, 57, 58, 58, 58, 58, 59, 60, 60, 60, 60, 60, 62, 62, 62, 62, 62, 63, 64, 64, 64, 64, 64, 64,
               65, 66, 66, 66, 66, 66, 66, 66, 67, 68, 68, 68, 68, 68, 68, 69, 70, 70, 70, 70, 70, 70, 70, 71,
               72, 72, 72, 72, 72, 73, 74, 74, 74, 74, 75, 76, 76, 76, 76, 76, 77, 77, 78, 78, 78, 78, 78, 79,
               80, 80, 80, 80, 80, 80, 80, 81, 81, 82, 82, 82, 82 ]

    Nstable = np.array(Nstable)
    Zstable = np.array(Zstable)

    is_stable = np.any((Zstable == Z) & (Nstable == N))
    return is_stable




def get_half_life(Z,N,Qval):
    """
      Calculate the half life of an alpha decay. The fit parameters
      hlog_1, and fitpars_1 are from Dong & Ren 2005
      (https://ui.adsabs.harvard.edu/abs/2005EPJA...26...69D/abstract).

      The function returns the a0 reaclib parameter as well as the half-life
    """
    hlog_1    = [0e0, 0.8937e0, 0.5720e0, 0.9380e0]
    fitpars_1 = [1.64062e0, -8.54399e0, -0.19430e0, -33.9054e0]

    fitpars_2 = [1.71182818,  -7.50480827,  -0.25315377, -30.70284443]
    hlog_2    = [0e0,0.047594788143095035,0.12139516604487433,0.39333238949044014]

    fitpars_3 = [1.70874516,  -7.52264755,  -0.25153138, -30.82453316]
    hlog_3    = [0e0, 0.2140471879925886,0.06002092026985373,0.49989717651371746]

    fitpars_4 = [1.71370933,  -7.34225699,  -0.2497752 , -30.6826417]
    hlog_4    = [0e0, -0.12418256402907601,1.179881795937427,0.7165978784278064]

    # Get correct fit parameters for the individual regions
    if Z>=82 and N>126:
        fitpars = fitpars_1
        hlog    = hlog_1
    elif 82<N<=126 and Z>=82:
        fitpars = fitpars_2
        hlog    = hlog_2
    elif 82<N<=126 and Z<82:
        fitpars = fitpars_3
        hlog    = hlog_3
    elif 50<N<=82 and 50<=Z<82:
        fitpars = fitpars_4
        hlog    = hlog_4
    else:
        return np.inf,np.inf

    if Z%2 == 0:
        if N%2 == 0:
            helper_idx = 0
        else:
            helper_idx = 1
    else:
        if N%2 == 0:
            helper_idx = 2
        else:
            helper_idx = 3

    log10Ta = (fitpars[0]*Z + fitpars[1])*Qval**(-0.5) + (fitpars[2]*Z + fitpars[3]) + hlog[helper_idx]
    try:
        Ta = 10**log10Ta
        a0 = np.log(np.log(2)/Ta)
    except OverflowError:
        a0 = np.inf
        Ta = np.inf

    return a0,Ta


def get_name(Z,N,Zdict,Ndict):
    """
      Get the name of a nucleus with Z and N from the
      dictionaries Zdict and Ndict.
    """
    names  = np.array(list(Zdict.keys()))
    Zarray = np.array([x[1] for x in np.array(list(Zdict.items()))]).astype(float)
    Narray = np.array([x[1] for x in np.array(list(Ndict.items()))]).astype(float)
    try:
        idx = np.where((Zarray==Z) & (Narray==N))[0]
        name = names[idx][0]
    except IndexError:
        name = None

    return name



# Get mass excesses and other data from winvn
mexc = {}
Z_number = {}
N_number = {}
with open(path_winvn,'r') as f:
    lines = f.readlines()
    for line in lines:
        if len(line.split())==7:
            mexc[line.split()[0].strip()]     = float(line.split()[5].strip())
            Z_number[line.split()[0].strip()] = float(line.split()[2].strip())
            N_number[line.split()[0].strip()] = float(line.split()[3].strip())

# Prepare output string
out = ""

# Get mass excess of He4
mexc_he4 = mexc['he4']

for nuc_name in tqdm(mexc.keys()):
    if Z_number[nuc_name]<=50:
        continue

    if is_stable_nucleus(Z_number[nuc_name],N_number[nuc_name]):
        continue

    # Get mass excesses of parent
    mexc_parent = mexc[nuc_name]
    # Get mass excess of daughter
    daughter_name = get_name(Z_number[nuc_name]-2,N_number[nuc_name]-2,Z_number,N_number)
    if not (daughter_name is None):
        mexc_daughter = mexc[daughter_name]
        # Get the Q-value
        Qval = mexc_parent - mexc_daughter - mexc_he4
        # Ignore nuclei with negative Q-values
        if Qval<0:
            continue
        # Get the half life and a0
        a0, Thalf = get_half_life(Z_number[nuc_name],N_number[nuc_name],Qval)

        # Save the result to a string
        # Ignore very long lived rates
        if Thalf>1e20:
            continue
        out   += nuc_name.rjust(5)+"   "+"{:.5e}".format(Thalf)+"\n"


# Save the rates
with open("alpha_decays.dat","w") as f:
    f.write(out)