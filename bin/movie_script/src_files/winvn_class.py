import pandas as pd
import numpy  as np
from scipy.interpolate import interp1d



class winvn(object):

    def __init__(self,path):
        """
          Class to read and manage a winvn
        """
        self.__path = path


    def reset_winvn(self):
        self.__df = self.__df_reset

    def read_winvn(self):
        """
          Read the winvn. (stolen from Carlos)
        """
        with open(self.__path,'r') as winvn:
            """Open file with the context manager and extract data."""
            # read file content
            winvn.readline()

            # read the ugly log temperature string
            self.T_string = winvn.readline().strip()
            str_size = 3
            logT = np.array([int(self.T_string[i:i+str_size])
                    for i in range(0,len(self.T_string),str_size)])
            logT[-1]*=10 # I think there is a factor 10 missing!!!
            logT = logT*1e-2

            # skip the hacky name part, i.e. look for the last name,
            # it appears twice
            name1,name2 = ("","-")
            nameList = []
            while (name1!=name2):
                name2 = name1
                name1 = winvn.readline().strip()
                nameList.append(name1)

            # get rid of the last name which appears twice...
            #...why don't they just give the number of elements!!!
            nameList.pop()


            # initialize element data container
            eleList = []


            # now, get the data for each of these elements
            for element in nameList:
                # extract element data
                #..what is mass excess again??!! I hate c12 related quantities!!
                wv_l = winvn.readline().split()
                if len(wv_l)==7:
                    name, A, Z, N, sp, mass_excess, no_idea = wv_l
                else:
                    name, A, Z, N, sp, mass_excess = wv_l
                    no_idea="unknown"
                # I don't trust that guys! Let's check the order!
                if (name != element):
                    sys.exit("OMG!!!")
                # extract partition function
                pf_list = []
                # read three lines...
                for i in range(3):
                    pf_list.extend(winvn.readline().split())
                pf_values = np.array(pf_list,dtype=float)

                func = interp1d(logT,pf_values,bounds_error=False,kind="cubic",fill_value="extrapolate")
                # save element data in a list
                eleList.append([int(Z),int(N),name,float(A),float(sp),
                                    float(mass_excess),str(no_idea), pf_values,func] )

        # give the element data list entries names
        column_names = ['Z','N','name','A','spin','mass excess','mass model',
                        'partition function','function']
        # build table
        self.__df = pd.DataFrame(eleList)
        # set column names
        self.__df.columns=column_names
        self.__df = self.__df.set_index("name")
        self.__df["name"] = self.__df.index

        self.__df_reset = self.__df.copy()

        # Calculate binding energies
        mexc_prot = self.__df.loc["p","mass excess"]
        mexc_neut = self.__df.loc["n","mass excess"]
        self.__df["binding energy"] = mexc_prot*self.__df["Z"]+mexc_neut*self.__df["N"]-self.__df["mass excess"]

    def filter_with_sunet(self,nuclei):
        self.__df = self.__df.loc[nuclei,:]

    def write_winvn(self,path="winvn.dat"):
        """
          Save the winvn to a file.
        """
        out=""
        out+="\n"
        out+=self.T_string+"\n"

        # Create the long list of nuclei
        for ind,row in self.__df.iterrows():
            out+=row["name"].rjust(5)+"\n"
        out+=row["name"].rjust(5)+"\n"

        for ind,row in self.__df.iterrows():
            out+=row["name"].rjust(5)+5*" "+"{:.3f}".format(row["A"]).rjust(7)+\
                 " "+str(int(row["Z"])).rjust(3)+" "+str(int(row["N"])).rjust(3)+\
                 "{:.1f}".format(row["spin"]).rjust(6)+"{:.3f}".format(row["mass excess"]).rjust(10)+\
                 row["mass model"].rjust(6)+\
                 "\n"
            out+=" "
            # Now the partition functions
            for i in range(len(row['partition function'])):
                out+="{:.5E}".format(row['partition function'][i]).rjust(12)
                if (i+1)%8 == 0:
                    out+="\n "
            out = out[:-1]

        with open(path,"w") as f:
            f.write(out)


    # def get_pf_Z_N(self):


    def rate_factor(self,reactants,products,temperature):
        """
            Calculate the rate factor for a given reaction
            Reactants: List of nuclei of the reactants of the inverse reaction
            Products: List of nuclei of the products of the inverse reaction
            temperature: temperature in GK
        """
        # Nominator
        nom = 1e0
        for reac in reactants:
            nom= nom*self.__df.loc[reac,"function"](temperature)
        # Denominator
        den = 1e0
        for prod in products:
            den= den*self.__df.loc[prod,"function"](temperature)
        factor = den/nom
        return factor


    def get_dataframe(self):
        """
          Get the dataframe
        """
        return self.__df

    def set_dataframe(self,value):
        """
          Get the dataframe
        """
        self.__df = value


if __name__=="__main__":
    w = winvn("data/winvne_v2.0.dat")
    w.read_winvn()
    fac = w.rate_factor(["mg24"],["he4","ne20"],np.linspace(1,5))
    print(fac)
