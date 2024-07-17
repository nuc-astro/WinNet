# Author: M. Reichert
# Date: 04.07.2024
import numpy as np
from tqdm import tqdm
import h5py
import os



class wreader(object):
    """
       Minimalistic class to lazily read WinNet data.
    """

    def __init__(self, path):
        """
        Initialize the class
           - path: Path to the WinNet data
        """
        # The path to the WinNet run
        self.path = path

        # The path to the hdf5 file
        self.filename = os.path.join(path, "WinNet_data.h5")

        # Check whether hdf5 mode was activated or not
        if not os.path.exists(self.filename):
            self.__mode = 'ascii'
            self.__snapshot_path = os.path.join(path, "snaps")
        else:
            self.__mode = 'hdf5'


    @property
    def A(self):
        """
        Mass number
        """
        if not hasattr(self,"_wreader__A"):
            self.__read_A_Z_N()
        return self.__A

    def __read_A_Z_N(self):
        """
        Read the mass number A
        """
        # Check if Ascii (mode = 2), hdf5 (mode = 1), or not present (mode = 0)
        mode = self.check_existence('snapshot')

        if mode == 2:
            self.__N, self.__Z = np.loadtxt(os.path.join(self.__snapshot_path, "snapsh_0001.dat"),unpack=True, usecols=(0,1),
                                            skiprows=3,dtype=int)
            self.__A = self.__N + self.__Z
        elif mode == 1:
            with h5py.File(self.filename, 'r') as hf:
                self.__A = hf['snapshots/A'][:]
                self.__Z = hf['snapshots/Z'][:]
                self.__N = hf['snapshots/N'][:]
        elif mode == 0:
            error_msg = 'Failed to read A, Z, and N in "'+str(self.path)+'". '+\
                        'Snapshots not present as Ascii nor as Hdf5!'
            raise ValueError(error_msg)


    @property
    def Z(self):
        """
        Atomic number
        """
        if not hasattr(self,"_wreader__Z"):
            self.__read_A_Z_N()
        return self.__Z

    @property
    def N(self):
        """
        Neutron number
        """
        if not hasattr(self,"_wreader__N"):
            self.__read_A_Z_N()
        return self.__N


    def __read_snapshots(self):
        """
        Read the snapshots
        """
        # Check if Ascii (mode = 2), hdf5 (mode = 1), or not present (mode = 0)
        mode = self.check_existence('snapshot')

        if mode == 2:
            # Get list of files (no directories)
            snapshot_files = [f for f in os.listdir(self.__snapshot_path) if os.path.isfile(os.path.join(self.__snapshot_path, f))]
            # Now only take files that start with a 'snapsh_'
            snapshot_files = [f for f in snapshot_files if f.startswith('snapsh_')]
            self.__snapshots_time = np.zeros(len(snapshot_files))
            self.__snapshots_Y = np.zeros((len(snapshot_files), len(self.A)))
            self.__snapshots_X = np.zeros((len(snapshot_files), len(self.A)))
            for i, f in enumerate(tqdm(snapshot_files, desc='Reading snapshots')):
                fname = 'snapsh_' + str(i+1).zfill(4) + '.dat'
                with open(os.path.join(self.__snapshot_path, fname), 'r') as file:
                    lines = file.readlines()
                    self.__snapshots_time[i] = float(lines[1].split()[0])
                    for j, line in enumerate(lines[3:]):
                        self.__snapshots_Y[i, j] = float(line.split()[2])
                        self.__snapshots_X[i, j] = float(line.split()[3])
        elif mode == 1:
            with h5py.File(self.filename, 'r') as hf:
                self.__snapshots_time = hf['snapshots/time'][:]
                self.__snapshots_Y = hf['snapshots/Y'][:]
                self.__snapshots_X = self.__snapshots_Y*self.A
        elif mode == 0:
            error_msg = 'Failed to read timescales in "'+str(self.path)+'". '+\
                        'Not present as Ascii nor as Hdf5!'
            raise ValueError(error_msg)



    def check_existence(self, entry):
        """
        Check whether an entry is in the hdf5 format (return 1),
        or in the ascii format (return 2), or does not exist (return 0)
        """
        if   (entry == 'mainout'):
            value = self.__check_files('mainout.dat','mainout')
        elif (entry == 'timescales'):
            value = self.__check_files('timescales.dat','timescales')
        elif (entry == 'energy'):
            value = self.__check_files('generated_energy.dat','energy')
        elif (entry == 'snapshot'):
            value = self.__check_files('snaps/snapsh_0001.dat','snapshots')
        elif (entry == 'flows'):
            value = self.__check_files('flow/flow_0001.dat','flows')
        elif (entry == 'finab'):
            value = self.__check_files('finab.dat','finab')
        elif (entry == 'finabsum'):
            value = self.__check_files('finabsum.dat','finab')
        elif (entry == 'finabelem'):
            value = self.__check_files('finabelem.dat','finab')
        else:
            error_msg = 'Checked for unknown entry "'+str(entry)+'". '
            raise ValueError(error_msg)

        return value


    def __check_files(self, ascii_file_path, hdf5_key):
        '''
          Check if something exists in hdf5 or ascii
        '''
        value = 0
        if not os.path.exists(self.filename):
            if not os.path.exists(os.path.join(self.path, ascii_file_path)):
                value = 0
            else:
                value = 2
        else:
            with h5py.File(self.filename, 'r') as hf:
                if hdf5_key in hf.keys():
                    value = 1
                elif os.path.exists(os.path.join(self.path, ascii_file_path)):
                    value = 2
                else:
                    value = 0

        return value

    @property
    def nr_of_snaps(self):
        """
        Number of snapshots
        """
        if not hasattr(self,"_wreader__nr_of_snaps"):
            self.__read_nr_of_snaps()
        return self.__nr_of_snaps

    def __read_nr_of_snaps(self):
        """
           Read the number of snapshots
        """
        # Check if Ascii (mode = 2), hdf5 (mode = 1), or not present (mode = 0)
        mode = self.check_existence('snapshot')

        if mode == 2:
            snapshot_files = [f for f in os.listdir(self.__snapshot_path) if os.path.isfile(os.path.join(self.__snapshot_path, f))]
            snapshot_files = [f for f in snapshot_files if f.startswith('snapsh_')]
            self.__nr_of_snaps = len(snapshot_files)
        elif mode == 1:
            with h5py.File(self.filename, 'r') as hf:
                self.__nr_of_snaps = len(hf['snapshots/time'][:])
        elif mode == 0:
            error_msg = 'Failed to read nr of snapshots in "'+str(self.path)+'". '+\
                        'Not present as Ascii nor as Hdf5!'
            raise ValueError(error_msg)


    @property
    def Y(self):
        """
        Get abundance at snapshot idx
        """
        if not hasattr(self,"_wreader__snapshots_time"):
            self.__read_snapshots()
        return self.__snapshots_Y

    @property
    def X(self):
        """
        Get mass fraction at snapshot idx
        """
        if not hasattr(self,"_wreader__snapshots_time"):
            self.__read_snapshots()
        return self.__snapshots_X

    @property
    def snapshot_time(self):
        """
        Get time at snapshot idx
        """
        if not hasattr(self,"_wreader__snapshots_time"):
            self.__read_snapshots()
        return self.__snapshots_time

    @property
    def tau(self):
        """
        Get the timescale of "tau", e.g., "tau_ag"
        """
        if not hasattr(self,"_wreader__timescales"):
            self.__read_timescales()
        return self.__timescales


    def __read_timescales(self):
        """
        Read the timescales
        """
        # Check if Ascii (mode = 2), hdf5 (mode = 1), or not present (mode = 0)
        mode = self.check_existence('timescales')

        if (mode == 2):
            self.__timescales = {}
            with open(os.path.join(self.path, "timescales.dat"), 'r') as f:
                lines = f.readlines()
                header = lines[0]
                header = header.replace('#', '').replace('[s]', '').replace('[GK]', '')
                key = header.split()
            data = np.loadtxt(os.path.join(self.path, "timescales.dat"), unpack=True)
            for i, k in enumerate(key):
                self.__timescales[k] = data[i]

        elif (mode == 1):
            with h5py.File(self.filename, 'r') as hf:
                self.__timescales = {}
                for key in hf['timescales'].keys():
                    self.__timescales[key] = hf['timescales'][key][:]
        elif (mode == 0):
            error_msg = 'Failed to read timescales in "'+str(self.path)+'". '+\
                        'Not present as Ascii nor as Hdf5!'
            raise ValueError(error_msg)


    @property
    def mainout(self):
        """
        Get an entry from the mainout
        """
        if not hasattr(self,"_wreader__mainout"):
            self.__read_mainout()
        return self.__mainout

    def __read_mainout(self):
        """
        Read the mainout
        """
        # Check if Ascii (mode = 2), hdf5 (mode = 1), or not present (mode = 0)
        mode = self.check_existence('mainout')

        if mode == 2:
            columns = ['iteration', 'time', 'temp', 'dens', 'ye',
                       'rad', 'yn', 'yp', 'ya', 'ylight', 'yheavy',
                       'zbar', 'abar', 'entr']
            data = np.loadtxt(os.path.join(self.path, "mainout.dat"), unpack=True)
            self.__mainout = {}
            for i, col in enumerate(columns):
                self.__mainout[col] = data[i]
        elif mode == 1:
            with h5py.File(self.filename, 'r') as hf:
                self.__mainout = {}
                for key in hf['mainout'].keys():
                    self.__mainout[key] = hf['mainout'][key][:]
        elif mode == 0:
            error_msg = 'Failed to read mainout in "'+str(self.path)+'". '+\
                        'Not present as Ascii nor as Hdf5!'
            raise ValueError(error_msg)


    @property
    def energy(self):
        """
        Get the energy
        """
        if not hasattr(self,"_wreader__energy"):
            self.__read_energy()
        return self.__energy

    def __read_energy(self):
        """
        Read the energy
        """
        # Check if Ascii (mode = 2), hdf5 (mode = 1), or not present (mode = 0)
        mode = self.check_existence('energy')

        if mode == 2:
            columns = ['time', 'engen_tot', 'S_src', 'engen_beta', 'engen_ng_gn', 'engen_pg_gp',
                       'engen_ag_ga', 'engen_np_pn', 'engen_na_an', 'engen_ap_pa', 'engen_fiss']
            data = np.loadtxt(os.path.join(self.path, "generated_energy.dat"), unpack=True)
            self.__energy = {}
            for i, col in enumerate(columns):
                self.__energy[col] = data[i]
        elif mode == 1:
            with h5py.File(self.filename, 'r') as hf:
                self.__energy = {}
                for key in hf['energy'].keys():
                    self.__energy[key] = hf['energy'][key][:]
        elif mode == 0:
            error_msg = 'Failed to read energy in "'+str(self.path)+'". '+\
                        'Not present as Ascii nor as Hdf5!'
            raise ValueError(error_msg)

    @property
    def finab(self):
        """
        Get the final abundances from the finab.dat file
        """
        if not hasattr(self,"_wreader__finab"):
            self.__read_finab()
        return self.__finab

    def __read_finab(self):
        """
        Reader of the finab
        """
        # Check if Ascii (mode = 2), hdf5 (mode = 1), or not present (mode = 0)
        mode = self.check_existence('finab')
        if mode == 2:
            A,Z,N,Y,X = np.loadtxt(os.path.join(self.path, "finab.dat"), unpack=True)
        elif mode == 1:
            with h5py.File(self.filename, 'r') as hf:
                A = hf['finab/finab/A'][:]
                Z = hf['finab/finab/Z'][:]
                N = A-Z
                Y = hf['finab/finab/Y'][:]
                X = hf['finab/finab/X'][:]
        elif mode == 0:
            error_msg = 'Failed to read finab in "'+str(self.path)+'". '+\
                        'Not present as Ascii nor as Hdf5!'
            raise ValueError(error_msg)
        self.__finab = {"A": A, "Z": Z, "N": N, "Y": Y, "X": X}


    @property
    def finabsum(self):
        """
        Get the final abundances from the finabsum.dat file
        """
        if not hasattr(self,"_wreader__finabsum"):
            self.__read_finabsum()
        return self.__finabsum

    def __read_finabsum(self):
        """
        Reader of the finabsum
        """
        # Check if Ascii (mode = 2), hdf5 (mode = 1), or not present (mode = 0)
        mode = self.check_existence('finabsum')
        if mode == 2:
            A,Y,X = np.loadtxt(os.path.join(self.path, "finabsum.dat"), unpack=True)
        elif mode == 1:
            with h5py.File(self.filename, 'r') as hf:
                A = hf['finab/finabsum/A'][:]
                Y = hf['finab/finabsum/Y'][:]
                X = hf['finab/finabsum/X'][:]
        elif mode == 0:
            error_msg = 'Failed to read finabsum in "'+str(self.path)+'". '+\
                        'Not present as Ascii nor as Hdf5!'
            raise ValueError(error_msg)
        self.__finabsum = {"A": A, "Y": Y, "X": X}


    @property
    def finabelem(self):
        """
        Get the final abundances from the finabelem.dat file
        """
        if not hasattr(self,"_wreader__finabelem"):
            self.__read_finabelem()
        return self.__finabelem

    def __read_finabelem(self):
        """
        Reader of the finabsum
        """
        # Check if Ascii (mode = 2), hdf5 (mode = 1), or not present (mode = 0)
        mode = self.check_existence('finabelem')
        if mode == 2:
            Z, Y = np.loadtxt(os.path.join(self.path, "finabelem.dat"), unpack=True)
        elif mode == 1:
            with h5py.File(self.filename, 'r') as hf:
                Z = hf['finab/finabelem/Z'][:]
                Y = hf['finab/finabelem/Y'][:]
        elif mode == 0:
            error_msg = 'Failed to read finabelem in "'+str(self.path)+'". '+\
                        'Not present as Ascii nor as Hdf5!'
            raise ValueError(error_msg)
        self.__finabelem = {"Z": Z, "Y": Y}

    def flow_entry(self, iteration, flow_group='flows'):
        """
        Get the flow entry
        """
        mode = self.check_existence('flows')

        if mode == 2:
            fname = os.path.join(self.path, "flow/flow_" + str(iteration).zfill(4) + ".dat")
            with open(fname, 'r') as f:
                lines = f.readlines()
                header = lines[1]
                time, temp, dens = np.array(header.split()).astype(float)

            # Read the whole flow
            nin,zin,yin,nout,zout,yout,flows = np.loadtxt(fname, unpack=True, skiprows=3)
            nin = nin.astype(int)
            zin = zin.astype(int)
            nout = nout.astype(int)
            zout = zout.astype(int)
            flow = {}
            flow['time']  = time
            flow['temp']  = temp
            flow['dens']  = dens
            flow['n_in']  = nin
            flow['p_in']  = zin
            flow['y_in']  = yin
            flow['n_out'] = nout
            flow['p_out'] = zout
            flow['y_out'] = yout
            flow['flow']  = flows
        elif mode == 1:
            with h5py.File(self.filename, 'r') as hf:
                flow = {}
                for key in hf[flow_group][str(iteration)].keys():
                    flow[key] = hf[flow_group][str(iteration)][key][:]
        elif mode == 0:
            error_msg = 'Failed to read flow in "'+str(self.path)+'". '+\
                        'Not present as Ascii nor as Hdf5!'
            raise ValueError(error_msg)

        return flow

if __name__ == "__main__":
    pass