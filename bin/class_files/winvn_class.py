import pandas as pd
import numpy  as np




class winvn(object):
    
    def __init__(self,path):
        """
          Class to read and manage a winvn
        """
        self.__path = path
        
        
        
    def read_winvn(self):
        """
          Read the winvn. (stolen from Carlos)
        """
        with open(self.__path,'r') as winvn:
            """Open file with the context manager and extract data."""
            # read file content
            winvn.readline()
            
            # read the ugly log temperature string
            T_string = winvn.readline().strip()
            str_size = 3
            logT = np.array([int(T_string[i:i+str_size])
                    for i in range(0,len(T_string),str_size)])
            logT[-1]*=10 # I think there is a factor 10 missing!!!
            
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
                name, A, Z, N, sp, mass_excess, no_idea = winvn.readline().split()
                # I don't trust that guys! Let's check the order!
                if (name != element): 
                    sys.exit("OMG!!!")
                # extract partition function
                pf_list = []
                # read three lines...
                for i in range(3):
                    pf_list.extend(winvn.readline().split())
                pf_values = np.array(pf_list,dtype=float)
                
                # save element data in a list
                eleList.append([int(Z),int(N),name,float(A),float(sp),
                                    float(mass_excess),str(no_idea), pf_values] )
        
        
        # give the element data list entries names
        column_names = ['Z','N','name','A','spin','mass excess','mass model',
                        'partition function'] 
        # build table
        self.__df = pd.DataFrame(eleList)
        # set column names
        self.__df.columns=column_names



    def get_dataframe(self):
        """
          Get the dataframe
        """
        return self.__df
