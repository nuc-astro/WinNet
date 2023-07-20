# Author: M. Reichert
import pandas           as pd
import numpy            as np
from  .nucleus_class    import nucleus
import os


class nucleus_multiple(object):
    """
      Class to deal with lists of nucleus objects
    """
    def __init__(self,names=None,A=None,Z=None,N=None,Y=None,X=None):
        """
          Initialize the class
          Input:
            nuclist     - list of nucleus names
        """
        self.names   = names
        self.A       =  A
        self.Z       =  Z
        self.N       =  N
        # Get sophisticated names, A, Z, N
        if not (names is None):
            self.__init_with_list()
        elif not(A is None) or not(Z is None) or not (N is None):
            self.__init_with_data()

        if not(Y is None):
            self.Y = np.array(Y)
        elif not(X is None):
            self.Y = np.array(X)/self.A
        else:
            self.Y = np.zeros(len(self.A))


    def __repr__(self):
        return self.df.__str__()

    def __init_with_list(self):
        """
          Initialize the class using the list of nucleinames
        """
        self.A = np.array([nucleus(n).get_A() for n in self.names])
        self.Z = np.array([nucleus(n).get_Z() for n in self.names])
        self.N = self.A-self.Z
        self.names = np.array(self.names)

    def __init_with_data(self):
        """
          Initialize the class using A,Z,N
        """
        # Sanity checks
        if np.sum([not(self.A is None), not(self.Z is None), not (self.N is None)]) <2:
            raise Exception("You need to provide at least two of A, Z, N!")

        if not (self.A is None):
            if not (self.Z is None):
                self.A = np.array(self.A)
                self.Z = np.array(self.Z)
                # Sanity check
                if len(self.A)!=len(self.Z):
                    raise Exception("A and Z had different lengths (A: "+str(len(self.A))+", Z:"+str(len(self.Z))+")!")
                self.N = self.A-self.Z
            elif not (self.N is None):
                self.A = np.array(self.A)
                self.N = np.array(self.N)
                # Sanity check
                if len(self.A)!=len(self.N):
                    raise Exception("A and N had different lengths (A: "+str(len(self.A))+", N:"+str(len(self.N))+")!")
                self.Z = self.A-self.N
            else:
                raise Exception("Something strange happened!")
        elif not (self.Z is None):
            self.N = np.array(self.N)
            self.Z = np.array(self.Z)
            if len(self.Z)!=len(self.N):
                raise Exception("N and Z had different lengths (N: "+str(len(self.N))+", Z:"+str(len(self.Z))+")!")
            self.A = self.Z+self.N
        else:
            raise Exception("Something strange happened!")

        # Sanity check
        if np.min(self.N)<0:
            raise Exception("N was smaller zero!")
        elif np.min(self.Z)<0 :
            raise Exception("Z was smaller zero!")
        elif np.min(self.A)<0:
            raise Exception("A was smaller zero!")

        # Create the names
        self.names = np.array([nucleus(N=self.N[ind],Z=self.Z[ind]).get_name() for ind in range(len(self.N))])



    def __sum_over(self,A,Y):
        max_A = max(A)
        # second_dimension = len(Y[0,:])
        Y_new = np.zeros([max_A+1,])
        for i in range(len(A)):
            Y_new[int(A[i])] += Y[i]
        return np.array(range(max_A+1)),Y_new



    @property
    def X(self):
        """
          Mass fraction
        """
        return self.Y*self.A

    @X.setter
    def X(self,value):
        """
          Mass fraction
        """
        self.Y=value/self.A


    @property
    def A_X(self):
        """
          Get sum over equal A's
        """
        atmp,xtmp = self.__sum_over(self.A,self.X)
        return np.array(atmp),np.array(xtmp)

    @property
    def Z_Y(self):
        """
          Get sum over equal Z's
        """
        ztmp,ytmp = self.__sum_over(self.Z,self.Y)
        return ztmp,ytmp

    @property
    def Z_X(self):
        """
          Get sum over equal Z's
        """
        ztmp,ytmp = self.__sum_over(self.Z,self.X)
        return ztmp,ytmp


    @property
    def Yprot(self):
        """
          Abundance of protons
        """
        mask = (self.A==1) & (self.Z==1)
        if sum(mask)!=1:
            raise Exception("Problem when getting hydrogen abundances!")
        else:
            return self.Y[mask][0]
    @property
    def Yneut(self):
        """
          Abundance of protons
        """
        mask = (self.A==1) & (self.N==1)
        if sum(mask)!=1:
            raise Exception("Problem when getting neutron abundances!")
        else:
            return self.Y[mask][0]
    @property
    def Yhe4(self):
        """
          Abundance of protons
        """
        mask = (self.A==4) & (self.Z==2)
        if sum(mask)!=1:
            raise Exception("Problem when getting alpha abundances!")
        else:
            return self.Y[mask][0]


    @property
    def df(self):
        """
          Get a dataframe out
        """
        df = pd.DataFrame()
        df["nucleus"] = self.names
        df["A"]       = self.A
        df["Z"]       = self.Z
        df["N"]       = self.N
        df["Y"]       = self.Y
        df["X"]       = self.Y*self.A
        return df


    def __merge_nuclei(self,A,Z):
        """
          Merge another list of nuclei into the own one.
        """
        A_merged = np.append(self.A,A)
        Z_merged = np.append(self.Z,Z)
        c_list = np.array([A_merged,Z_merged])
        # print(A_merged)
        # print(Z_merged)
        unique=np.unique(c_list,axis=1)
        return unique[0,:],unique[1,:]


    def __mul__(self,other):
        """
         Multiply either two instances of nuclei lists or a float with this instance.
        """
        if isinstance(other,float) or isinstance(other,int):
            self.Y=self.Y*other
        elif isinstance(other,nucleus_multiple):
            merge_A,merge_Z = self.__merge_nuclei(other.A,other.Z)
            merge_N         = merge_A-merge_Z
            Y_list          = np.zeros(len(merge_N))
            for ind,atmp in enumerate(merge_A):
                mask_self  = (self.A ==merge_A[ind]) & (self.Z ==merge_Z[ind])
                mask_other = (other.A==merge_A[ind]) & (other.Z==merge_Z[ind])
                if (np.any(mask_self)) and (np.any(mask_other)):
                    Y_list[ind] = self.Y[mask_self][0]*other.Y[mask_other][0]
            self.A = merge_A
            self.Z = merge_Z
            self.N = merge_N
            self.Y = Y_list
            self.__init_with_data()
        return self

    def __truediv__(self,other):
        """
         Multiply either two instances of nuclei lists or a float with this instance.
        """
        if isinstance(other,float) or isinstance(other,int):
            self.Y=self.Y/other
        elif isinstance(other,nucleus_multiple):
            merge_A,merge_Z = self.__merge_nuclei(other.A,other.Z)
            merge_N         = merge_A-merge_Z
            Y_list          = np.zeros(len(merge_N))
            for ind,atmp in enumerate(merge_A):
                mask_self  = (self.A ==merge_A[ind]) & (self.Z ==merge_Z[ind])
                mask_other = (other.A==merge_A[ind]) & (other.Z==merge_Z[ind])
                if (np.any(mask_self)) and (np.any(mask_other)):
                    Y_list[ind] = self.Y[mask_self][0]/other.Y[mask_other][0]
            self.A = merge_A
            self.Z = merge_Z
            self.N = merge_N
            self.Y = Y_list
            self.__init_with_data()
        return self



    def __rmul__(self,other):
        return self.__mul__(other)


    def __add__(self,other):
        """
          Add to instances of this class. The abundances are added for each nucleus.
        """
        if isinstance(other,float) or isinstance(other,int):
            self.Y=self.Y+other
        elif isinstance(other,nucleus_multiple):
            merge_A,merge_Z = self.__merge_nuclei(other.A,other.Z)
            merge_N         = merge_A-merge_Z
            Y_list          = np.zeros(len(merge_N))
            for ind,atmp in enumerate(merge_A):
                mask_self  = (self.A ==merge_A[ind]) & (self.Z ==merge_Z[ind])
                mask_other = (other.A==merge_A[ind]) & (other.Z==merge_Z[ind])
                if (np.any(mask_self)):
                    Y_list[ind] += self.Y[mask_self][0]
                if (np.any(mask_other)):
                    Y_list[ind] += other.Y[mask_other][0]
            self.A = merge_A
            self.Z = merge_Z
            self.N = merge_N
            self.Y = Y_list
            self.__init_with_data()
        return self

    def write_finab(self,path="finab.dat"):
        """
        Write the result (contained in self.__abundances) into a final abundance file. The format is:
        |A   Z   N   Yi   Xi|
        Input:
          path - Path to output file
        """
        Y = self.Y
        X = self.Y*self.A
        A = self.A
        Z = self.Z
        N = self.N

        output = np.array([A,Z,N,Y,X]).T
        np.savetxt(path,output,fmt=['%7i','%7i','%7i','    %1.9e','    %1.9e'],header='A    Z    N    Y    X')


    def write_seed(self,path=None):
        """
          Write a seed file
        """
        line=""
        for i in range(len(self.A)):
            X = self.Y[i]*self.A[i]
            mafra = '%1.6e' % float(X)
            if np.isnan(X) or X<=1e-15:
                continue
            nam = self.names[i].lower()
            # Take care of special names
            if (nam.strip() == 'neutron') or (nam.strip() == 'neutron1'):
                nam = 'n'
            elif nam.strip() == 'h1':
                nam = 'p'
            elif nam.strip() == 'h2':
                nam = 'd'
            elif nam.strip() == 'h3':
                nam = 't'

            line+= nam.rjust(5) + ' '*7 + mafra
            if i != len(self.A)-1:
                line+="\n"
        # Make a default path
        if path is None:
            path = "seed.dat"

        with open(path,"w") as f:
            f.write(line)


if __name__ == "__main__":
    """    Test the class    """

    nlist = ["o16","ni56"]
    A     = [1,16,17,56]
    Z     = [1,8,8,28]
    Y = [0.1,0.1,0.1,0.1]
    nucl = nucleus_multiple(A=A,Z=Z,Y=Y)
    A     = [2,16,10,56,56]
    Z     = [1,8,8,28,26]
    Y = [0.2,0.2,0.3,0.4,0.1]

    nucl2 = nucleus_multiple(A=A,Z=Z,Y=Y)
    print(nucl+nucl2)
    A,X=nucl2.A_X
