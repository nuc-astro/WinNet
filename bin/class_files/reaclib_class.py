# Reaclib class to analyze and change reaclib files
# Author: Moritz Reichert
# Date  : 08.05.17

from .nucleus_class             import nucleus
import pandas                   as pd
import numpy                    as np
import matplotlib.pyplot        as plt
import os
import sys
import warnings



class reaclib(object):

    def __init__(self,path,quiet=False):
        """
          Analyze a rate file for a given path.
          Input:
            path    - path to reaclib file
            quiet   - Do you want to have status informations of the progress?
        """
        # Turn off panda warning when changing copied frames
        pd.options.mode.chained_assignment = None
        # Make warnings to exceptions to be able to catch them
        warnings.filterwarnings('error')

        self.__path = path      # Path to reaclib file (independent from format)
        self.__format = 1       # 1 original v1, without chapter 9-11

        self.__quiet = quiet    # Verbosity of the class

        # give some column names
        self.__column_names = ['chapter','reactants','products','label','type','reverse','q value','error','reaction type'] + ['a'+str(i) for i in range(7)]
        self.__pd_reaclib = None
        # Have also a backup of the frame to be able to reset it
        self.__reset_reaclib = None
        # Default settings of reaclib file format
        self.__cpn = 5          # characters per nucleus
        self.__lpo = 5+6*5+8    # position of the label
        self.__lc = 4           # label characters
        self.__cpe = 13         # characters per entry in the set

        # Fission rates are not matching their chapter
        self.__ignore_fission_errors = True
        self.__fission_flags = ['ms99','fiss','mp01','sfis']


        self.__temperature_range = np.logspace(-2,1,num=20)    # Range for the temperature grid for plots and tests

        self.__test_rate_upper_limit = 1.e100

    def read_reaclib(self):
        """
          read the whole reaclib file.
        """
        # Use a different functions for different chapters
        read1 = lambda x : ([x[0]] , [x[1]])
        read2 = lambda x : ([x[0]] , [x[1],x[2]])
        read3 = lambda x : ([x[0]] , [x[1],x[2],x[3]]) if (len([a for a in x if a != '']) == 4 or (not self.__format == 1)) else ([x[0]] , [x[1],x[2],x[3],x[4]])
        read4 = lambda x : ([x[0],x[1]] , [x[2]])
        read5 = lambda x : ([x[0],x[1]] , [x[2],x[3]])
        read6 = lambda x : ([x[0],x[1]] , [x[2],x[3],x[4]])
        read7 = lambda x : ([x[0],x[1]] , [x[2],x[3],x[4],x[5]])
        read8 = lambda x : ([x[0],x[1],x[2]] , [x[3]]) if (len([a for a in x if a != '']) == 4 or (not self.__format == 1)) else ([x[0],x[1],x[2]] , [x[3],x[4]])
        read9 = lambda x : ([x[0],x[1],x[2]] , [x[3],x[4]])
        read10 = lambda x : ([x[0],x[1],x[2],x[3]] , [x[4],x[5]])
        read11 = lambda x : ([x[0]] , [x[1],x[2],x[3],x[4]])

        # create a parser which picks the 'readX' function depending on the chapter.
        r_parser = {1 : read1, 2 : read2, 3 : read3, 4 : read4, 5 : read5,
                    6 : read6, 7 : read7, 8 : read8, 9 : read9, 10 : read10,
                    11 : read11}

        check = 0
        reaclib_data = []
        # Open the reaclib file
        with open(self.__path,'r') as reaclib_file:
            count = 0
            # read file content
            for line in reaclib_file.readlines():
                count += 1
                if not self.__quiet and (count % 1000 == 0):
                    print('Reading line ' +str(count) +"         ",end='\r')

                # Check for chapter line
                if line.strip().isdigit():
                    chapter = int(line.strip())
                    continue

                # Check for empty lines after the chapter and ignore them
                if line.strip() == '':
                    continue

                check +=1
                # Line with nuclei and Q-value and so on
                if check == 1:
                    involved = [line[i:i+self.__cpn].strip() for i in range(1*self.__cpn, 7*self.__cpn, self.__cpn)]
                    reactants, products = r_parser[chapter](involved)
                    # Read all the rest
                    loc = self.__lpo
                    label = line[loc:loc+self.__lc] # reaction label
                    loc += self.__lc
                    r_type = line[loc:loc+1] # reaction type
                    loc += 1
                    r_reverse = (line[loc:loc+1] == 'v') # reverse rate flag
                    loc += 4
                    r_qval = float(line[loc:loc+12].replace('D','e').replace('d','e')) # q value
                    # Make first error tests
                    if (len([a for a in involved if a != '']) != (len(reactants) + len(products))) and (not (self.__ignore_fission_errors and (label.strip() in self.__fission_flags))):
                        error = 'Not matching its chapter!'
                    else:
                        error = ''
                elif check == 2:
                    a_factors = [line[i:i+self.__cpe] for i in range(0, 4*self.__cpe, self.__cpe)]
                elif check == 3:
                    # three factors in the second line
                    a_factors.extend([line[i:i+self.__cpe] for i in range(0, 3*self.__cpe, self.__cpe)])
                    try:
                        a_factors = [float(x.replace('D','e').replace('d','e')) for x in a_factors]
                    except:
                        error = 'Fitting value is not a float!'
                    # gather information about the rate
                    reaction_type = self.__classify_reaction(chapter,reactants,products,r_type,label,a_factors)

                    # Winnet requires neutrons first and heavier nuclei after that
                    if label not in self.__fission_flags:
                        reactants = sorted(reactants)
                        products = sorted(products)
                    # Save the rate
                    rate_data = [chapter,reactants,products,label,r_type,r_reverse,r_qval,error,reaction_type]
                    rate_data.extend(a_factors)

                    reaclib_data.append(rate_data)
                    check = 0

            # Make the last output
            if not self.__quiet:
                print('Reading line ' +str(count) +"         ")

            # create data frame
            self.__pd_reaclib = pd.DataFrame(reaclib_data)
            # set column names
            self.__pd_reaclib.columns=self.__column_names
            # Create a copy of it to be able to reset the dataframe
            self.__reset_reaclib = self.__pd_reaclib.copy()


    def reset_reaclib(self):
        """
          Reset the reaclib to the loaded one.
        """
        self.__pd_reaclib = self.__reset_reaclib.copy()



    def __classify_reaction(self,chapter,reactants,products,weak,label,afactors):
        """
          Classify a reaction if it is gamma-n, n-gamma, ...
          Input:
            chapter   - reaclib chapter of the reaction (int)
            reactants - all reactants of the reaction (list of strings)
            products  - all products of the reaction (list of strings)
            weak      - 'w' if weak, else ''
            label     - reaction label to identify fission rates
            afactors  - check whether a rate is constant in dependence of the temperature
          Output:
            Returns a string, containing the type of reaction
        """

        # is fission reaction?
        if label in self.__fission_flags:
            reaction_type = 'fission'
        # is alpha decay?
        elif (chapter == 2) and ('he4' in products) and (sum(afactors[1:])==0):
            reaction_type = 'a-decay'
        # is weak reaction?
        elif weak.strip() =='w':
            reaction_type = 'weak'
        # is n-gamma?
        elif (chapter == 4) and ('n' in reactants):
            reaction_type = 'n-gamma'
        # is gamma-n?
        elif (chapter == 2) and ('n' in products):
            reaction_type = 'gamma-n'
        # is p-gamma?
        elif (chapter == 4) and ('p' in reactants):
            reaction_type = 'p-gamma'
        # is gamma-p?
        elif (chapter == 2) and ('p' in products):
            reaction_type = 'gamma-p'
        # is a-gamma?
        elif (chapter == 4) and ('he4' in reactants):
            reaction_type = 'a-gamma'
        # is gamma-a?
        elif (chapter == 2) and ('he4' in products) and (sum(afactors[1:])!=0):
            reaction_type = 'gamma-a'
        # is n-p?
        elif (chapter == 5) and ('n' in reactants) and ('p' in products):
            reaction_type = 'n-p'
        # is p-n?
        elif (chapter == 5) and ('p' in reactants) and ('n' in products):
            reaction_type = 'p-n'
        # is n-a?
        elif (chapter == 5) and ('n' in reactants) and ('he4' in products):
            reaction_type = 'n-a'
        # is a-n?
        elif (chapter == 5) and ('he4' in reactants) and ('n' in products):
            reaction_type = 'a-n'
        # is p-a?
        elif (chapter == 5) and ('p' in reactants) and ('he4' in products):
            reaction_type = 'p-a'
        # is a-p?
        elif (chapter == 5) and ('he4' in reactants) and ('p' in products):
            reaction_type = 'a-p'
        else:
            reaction_type = 'other'
        return reaction_type


    def get_statistics(self):
        """
          Get the amount of different rate types (n-gamma, gamma-n, weak ,alpha-n,...)
          Output:
            Returns a pandas series, containing the counts of the different reaction types and the total number of reactions
        """
        # Count the occurance of the types and append the amount of reactions
        return (pd.value_counts(self.__pd_reaclib['reaction type'].values).append(pd.Series({'Total reactions' :self.__pd_reaclib.count()[0]})))


    def set_rate(self,dataframe):
        """
          Set the dataframe (containing all reactions) of the class with the given one
        """
        self.__pd_reaclib = dataframe.copy()


    def __reaclib_fit_function(self,temperature,coefficients):
        """
          Returns the value of the rate for a given temperature
          Input:
            temperature  - Temperature [GK]
            coefficients - The 7 fit coefficients from the reaclib format
        """
        # Calculate the exponent
        exponent = coefficients[0] + coefficients[6] * np.log(temperature)

        for i in range(5):
            exponent += coefficients[i+1]*temperature**((2.*(i+1)-5.)/3.)
        # Calculate the rate
        try:
            rate = np.exp(exponent)
        except RuntimeWarning:
            rate = np.inf
        return rate

    def test_reaclib(self, test_weak=True, test_mass=True, test_overflow=True):
        """
          Test the reaclib for failures
        """
        # Amount of reactions
        total_reactions = self.__pd_reaclib['chapter'].count()
        # Loop through the entries
        for index, row in self.__pd_reaclib.iterrows():

            if not self.__quiet and (index % 100 == 0):
                percentage = str(int((index+1.)/total_reactions * 10000.) / 100.).ljust(4)
                print('Running tests: '+percentage.ljust(4) + " done!          ",end='\r')

            # Has errors anyways
            if row['error'] != '':
                continue

            # Ignore fission for test or not
            if self.__ignore_fission_errors and (self.__pd_reaclib['label'].iloc[index] in self.__fission_flags):
                continue

            # Test for consistency of weak rates
            if test_weak:
                if ((row['type'].strip() == 'w') and self.__test_weakconsistency(row)):
                    self.__pd_reaclib['error'].iloc[index] = 'Weak rate is not consistent!'
                    continue

            # Test for mass consistency
            if test_mass:
                if self.__test_massconsistency(row) and test_mass:
                    self.__pd_reaclib['error'].iloc[index] = 'Mass inconsistency of the rate!'
                    continue

            # Test for overflows of the rate
            if test_overflow:
                if self.__test_overflows(row):
                    self.__pd_reaclib['error'].iloc[index] = 'Rate exceeds value of ' + str(self.__test_rate_upper_limit)
                    continue

        if not self.__quiet:
            print('Running tests: 100 done!          ')


    def get_contained_nuc(self):
        """
          Get all nuclei contained in reaclib file.
          Output:
            list of strings, containing the names of all nuclei contained in the reaclib
        """
        flatten = lambda l: [item for sublist in l for item in sublist] # From https://stackoverflow.com/questions/952914/making-a-flat-list-out-of-list-of-lists-in-python
        prods  = self.__pd_reaclib["products"].tolist()
        reacts = self.__pd_reaclib["reactants"].tolist()
        prods.extend(reacts)
        return list(set(flatten(prods)))

    def test_net_consistency(self,path):
        """
          Test if nuclei are included in the net file, but not in the reaclib. Be careful, reaclib contains nuclei where you are trapped in.
          There are reactions producing that nuclei, but no reaction to destroy it. (e.g. "f16" in the 2017 version of reaclib)
          Input:
            path - path to net file of WinNet
          Output:
            list of strings, containing problematic nuclei. If all nuclei have a corresponding reaction, an empty list is returned
        """
        # Nuclei from network file
        sunet_nucs    = np.loadtxt(path,unpack=True,dtype=str)
        # Nuclei from reaclib
        reaction_nucs = self.get_contained_nuc()

        # Get problematic nuclei
        problem_nuc = []
        for nuc in sunet_nucs:
            if nuc not in reaction_nucs:
                problem_nuc.append(nuc)

        return problem_nuc

    def __test_massconsistency(self,row):
        """
          Test the reaction rate for overflows!
        """
        error = False
        reactants = list(row['reactants'])
        products  = list(row['products'])

        # Transform to Z values
        reactants_A = [nucleus(x).get_A() for x in reactants]
        products_A  = [nucleus(x).get_A() for x in products]

        if sum(reactants_A) != sum(products_A):
            error = True
        return error


    def __test_overflows(self,row,temp_range=None,upper_lim=None):
        """
          Test the reaction rate for overflows!
        """
        error = False
        coefficients = list(row['a0':'a6'])

        if upper_lim is None:
            upper_lim = self.__test_rate_upper_limit
        if temp_range is None:
            temp_range = self.__temperature_range

        # loop trough temperatures
        for t in temp_range:
            if self.__reaclib_fit_function(t,coefficients) >= upper_lim:
                error = True
                break
        return error


    def __test_weakconsistency(self,row):
        """
          Test the reaction for consistency
        """
        error = False
        reactants = list(row['reactants'])
        products  = list(row['products'])

        # Transform to Z values
        reactants_Z = [nucleus(x).get_Z() for x in reactants]
        products_Z  = [nucleus(x).get_Z() for x in products]

        if (row['reaction type'].strip() == 'a-decay'):
            return error

        if sum(reactants_Z) == sum(products_Z):
            error = True
        return error


    def get_rate_at_temp(self,temperature,reactants,products):
        """
          Get a reaction rate at given temperature
          Input:
            temperature - Temperature [GK]
            reactants   - list with reactants
            products    - list with products
          Output:
            Reaction rate at a given temperature
        """
        # Get matching reaction, got damn this was complicated
        matching_reactants = self.__pd_reaclib[[x == sorted(reactants) for x in self.__pd_reaclib.reactants.values]]
        reactions = matching_reactants[[x == sorted(products) for x in matching_reactants.products.values]]
        # print(sum((self.__pd_reaclib.reactants.values == ['n']).astype(int)))
        # Get the reaction rate
        rat = 0.
        for index, row in reactions.iterrows():
            rat += self.__reaclib_fit_function(temperature,list(row['a0':'a6']))
        return rat




    def plot_rate(self,reactants,products,figure=None,axlabel=True,temp_grid=None,**kwargs):
        """
          Plot a rate.
          Input:
            reactants - list of reactants (list of strings)
            products  - list of products (list of strings)
        """
        # Make a figure if no figure was given
        if figure == None:
            final_fig = plt.figure()
        else:
            final_fig = figure

        ax = final_fig.gca()

        rat_list = []

        if temp_grid is None:
            temp_grid = self.__temperature_range

        for t in temp_grid:
            rat_list.append(self.get_rate_at_temp(t,reactants,products))
        ax.plot(temp_grid,rat_list,**kwargs)

        if axlabel:
            ax.set_xlabel('Temperature [GK]')
            ax.set_ylabel('Rate')
            ax.set_yscale('log')
            ax.set_xscale('log')

        return final_fig


    def get_rate_error(self):
        """
          Get a list with all rates that produce errors
        """
        return self.__pd_reaclib[self.__pd_reaclib['error'] != '']


    def get_rate_error_html(self,path):
        """
          Get an html file with all rates that contain an error
          Input:
            path - Path to the html file that will be created
        """
        f = open(path,'w')
        out = self.get_rate_error().to_html()
        f.write(out)
        f.close()

    def __convert_series_to_string(self,series):
        """
          Convert a pandas series to a reaclib string.
          Input:
            series - pandas series of one reaction
        """
        # Get one list with both, reactants and products
        clean_copy = series.copy()
        merged_nuclei = list(clean_copy['reactants'])
        merged_nuclei.extend(list(clean_copy['products']))
        # Make it out of 5 spaces again
        merged_nuclei = [x.rjust(5) for x in merged_nuclei]
        # Because reverse was saved as bool
        if series['reverse']:
            reverse = 'v'
        else:
            reverse = ' '

        five_spaces = 5 * ' '
        # First line:
        # Reaclib format: 5 spaces + 6 times 5spaces for nuclei + 8 spaces + 4 characters for the label + 1 character for a flag + 1 for reverse +
        # 3 spaces + 12 characters for the Q-value + 10 spaces
        out = five_spaces + ''.join(merged_nuclei) +(6-len(merged_nuclei)) * five_spaces + 8*' '+series['label'].ljust(4) + series['type'].ljust(1) + reverse +  \
              3*' ' + '%12.5e' % (series['q value']) + 10*' '
        out += '\n'
        # Second line:
        # a0 - a3 followed by 22 spaces
        out += ''.join(['%13.6e' % x for x in series['a0':'a3']])+ 22*' '
        out += '\n'
        # Third line:
        # a4 - a6 followed by 35 spaces
        out +=''.join(['%13.6e' % x for x in series['a4':'a6']])+ 35*' '
        out += '\n'
        return out

    def __sort_prods_reacs(self,in_list):
        """
          Sort the products and reactants of a given dataframe.
          A nucleus is greater if it has an higher proton number
        """
        helplist = [nucleus(x) for x in in_list]
        helplist = [x.get_input_name() for x in sorted(helplist)]
        return helplist



    def save_reaclib(self,filename,dataframe=None,sort=False):
        """
          Save the pandas dataframe to a reaclib file with the correct format.
          Input:
            filename  - Filename of the output.
            dataframe - Pandas dataframe that contains reaction rates. If none, the reaction rates of the class are used
        """
        # String to cache the output, before writing it to file
        reaclib_out = ''

        if dataframe is None:
            input_reaclib = self.__pd_reaclib.copy()
        else:
            input_reaclib = dataframe.copy()

        # Get amount of reactions for status
        amount_reactions = input_reaclib.count()[0]
        # Show more content of the string (in principle not necessary, but nice to have)
        pd.set_option('display.max_colwidth', -1)

        # Sort the reactions to correct order again
        if sort:
            if not self.__quiet:
                print("Sorting database. This may take a while!")
                print("Sorting products.")

            input_reaclib['products']  = input_reaclib['products'].apply(self.__sort_prods_reacs)
            if not self.__quiet:
                print("Sorting reactants.")
            input_reaclib['reactants'] = input_reaclib['reactants'].apply(self.__sort_prods_reacs)

            input_reaclib['Z'] = input_reaclib['reactants'].apply(lambda x: nucleus(x[0]).get_Z() )
            input_reaclib['A'] = input_reaclib['reactants'].apply(lambda x: nucleus(x[0]).get_A() )
            input_reaclib.sort_values(['Z','A'],inplace=True)




        if not self.__quiet:
            print('Saving reaclib to file.')
            print('Creating correct format, this may take a while.')

        # Create reaclib format and save it temporary into dataframe
        input_reaclib['tmp_out'] = input_reaclib.apply(self.__convert_series_to_string,axis=1)
        if not self.__quiet:
            print('Saving chapterwise.')

        # Determine the maximal chapter (necessary for different formats)
        max_chapter = input_reaclib['chapter'].max()

        # Cycle through chapters
        if not isinstance(max_chapter,int):
            if not self.__quiet:
                print('Warning: Chapter was not an integer, it was '+str(max_chapter)+' instead!')


        for i in range(1,int(max_chapter)+1):
            if not self.__quiet:
                print('Saving chapter '+str(i)+' of '+str(max_chapter)+'.      ',end='\r')

            # Set the chapter before the rates
            reaclib_out += str(i) + '\n\n\n'
            # Get the entries of the chapters
            tmp = (input_reaclib[input_reaclib['chapter']==i])
            reaclib_out +=(''.join(tmp['tmp_out'].values))

        if not self.__quiet:
            print('Saving to file.         ')

        # Save reaclib to file
        with open(filename,'w') as f_out:
            f_out.write(reaclib_out)

        # Remove the format again
        del input_reaclib['tmp_out']


    def get_critical_low_temperature_rates(self,min_temperature=1e-4,amount_points=20,max_rate=1.e100):
        """
          Get rates that are critical for low temperatures (< 1e-2 GK)
        """
        temp_points = np.logspace(np.log10(min_temperature),-2,num=amount_points)

        # Amount of reactions
        total_reactions = self.__pd_reaclib['chapter'].count()
        # Store the index of problematic reactions in list
        problematic_index = []
        # Order the index
        self.__pd_reaclib = self.__pd_reaclib.reset_index()
        # Loop through the entries
        for index, row in self.__pd_reaclib.iterrows():

            if not self.__quiet and (index % 100 == 0):
                percentage = str(int((index+1.)/total_reactions * 10000.) / 100.).ljust(4)
                print('Get critical rates: '+percentage.ljust(4) + " done!          ",end='\r')

            # Test for overflows of the rate
            if self.__test_overflows(row,temp_range=temp_points,upper_lim=max_rate):
                problematic_index.append(index)
                continue

        if not self.__quiet:
            print('Running tests: 100 done!          ')
        problematic_index = np.array(problematic_index)

        return self.__pd_reaclib.iloc[problematic_index]


    def drop_errors(self):
        """
          Drop errors from dataframe. Need to run test_reaclib() first
        """
        # Count the errors
        if not self.__quiet:
            amount_errors = self.__pd_reaclib[self.__pd_reaclib['error'] != ''].count()[0]
            print('Dropped '+str(amount_errors)+' rates.')
        # Only keep the rates without errors
        self.__pd_reaclib = self.__pd_reaclib[self.__pd_reaclib['error'] == '']


    def get_dataframe(self,reaction_type=None,chapter=None):
        """
          Returns the pandas dataframe, containing all reactions (after applying a filter) of the reaclib
          Input:
            reaction_type - only these reaction types (n-gamma,..), none for dont filter the type    (list of strings or string)
            chapter       - only these chapters                   , none for dont filter the chapter (list of integers or integer)
        """
        return self.__filter_rates(self.__pd_reaclib.copy(),reaction_type,chapter)


    def __analyze_update(self,dataframe):
        """
          Analyze a dataframe that contains the column "version".
          Output:
           Returns the amount of new inserted rates
        """
        # Count the values of new and old rates in the frame
        count = pd.value_counts(dataframe['version'].values)
        if 'new' in count.index:
            new_count = count['new']
        else:
            new_count = 0

        if 'old' in count.index:
            old_count = count['old']
        else:
            old_count = 0
        return new_count

    def get_all_reactiontypes(self):
        """
          Get all types of reactions, necessary to write in "reaction_type" when filtering rates
          Output:
            List of strings, containing the possible reaction types
        """
        types = self.get_statistics()
        return (list(types.index)[:-1])


    def is_included(self,series):
        """
          Is the rate included in the dataframe? Compared are only products and reactants.
          Returns true if yes, false if not
        """
        # Filter for the chapter
        chapter_reactions = self.__pd_reaclib[self.__pd_reaclib['chapter']==series['chapter']]
        reactants = series['reactants']
        products  = series['products']
        # Filter for reactants
        for r in reactants:
            rates = chapter_reactions[[r in x for x in chapter_reactions['reactants']]]
        # Filter for products
        if not rates.empty:
            for p in products:
                rates = rates[[p in x for x in rates['products']]]

        if rates.empty:
            return False
        else:
            return True




    def __filter_rates(self,dataframe,reaction_type,chapter):
        """
          Filter the reaction rates by type and chapter.
          Input:
            dataframe     - Pandas dataframe that contains all reactions
            reaction_type - only these reaction types (n-gamma,..), none for dont filter the type    (list of strings or string)
            chapter       - only these chapters                   , none for dont filter the chapter (list of integers or integer)
          Output:
            Filtered dataframe, that contains only the reactions with given type and chapter
        """
        # First convert everything to a list if the input is no list
        if isinstance(reaction_type,str):
            reaction_type=[reaction_type]
        if isinstance(chapter,int):
            chapter=[chapter]
        # Filter for the correct reactions
        # Distinguish three cases
        # First case: Filter for both
        if reaction_type != None and chapter != None:
            dataframe = dataframe[dataframe['chapter'].isin(chapter) & dataframe['reaction type'].isin(reaction_type)]
        # Second case: Filter for the chapter only
        elif reaction_type == None and chapter != None:
            dataframe = dataframe[dataframe['chapter'].isin(chapter)]
        # Third case: Filter for the reaction_type only
        elif reaction_type != None and chapter == None:
            dataframe = dataframe[dataframe['reaction type'].isin(reaction_type)]
        # (Fourth case: Don't do anything)
        return dataframe


    def update(self,reaclib_path=None,dataframe=None,reaction_type=None,chapter=None,ignore_label=False,ignore_reverse=False,ignore_type=False,replace_specific_label=None):
        """
          Update the reaclib with rates from other reaclib. Dont include new rates in the file. Both reaclibs have to be in the "old version" (without chapter 9-10) or in the
          "new version", with these chapters.
          Input:
            reaclib_path          - path of the reaclib that contains new rates
            dataframe             - pandas dataframe that contains reactions. If this is given, you don't have to give a reaclib_path.
            reaction_type         - Reaction type that will be updated (string or list of strings), for possible values see ".get_statistics()"
            chapter               - Chapter that will be updated (list of integer or integer)
            ignore_label          - ignore the label for updating and only look at the reactions itself
            ignore_reverse        - ignore if the reaction is reverse
            ignore_type           - ignore the type
            replace_specific_label- Should one specific label be replaced. E.g. Moeller 93 with Moeller 03 rates. For a list of the Jina labels see https://groups.nscl.msu.edu/jina/reaclib/db/labels.php
                                    if yes, the input should be a list with [old label, new label]
        """
        # Function to categorize 3 cases
        def categorize(x):
            # First case: There are new and old rates for the label (-> take the new ones)
            if ('new' in x) and ('old' in x):
                return 1
            # Second case: There are only new rates (-> don't take these rates)
            if ('new' in x) and not ('old' in x):
                return 2
            # Third case: There are only old rates (-> take all old rates)
            if not ('new' in x) and ('old' in x):
                return 3

        if not self.__quiet:
            print('Updating reaclib.')

        # Get the input values
        # First case: reaclib path is given, so use this as new dataframe
        if reaclib_path is not None and dataframe is None:
            # Reading other reaclib to dataframe and bring it to correct shape
            new_rates = reaclib(reaclib_path,quiet=True)
            new_rates.read_reaclib()
            new_rates_df = new_rates.get_dataframe()
        # Second case: dataframe is given, so use these reaction rates
        elif reaclib_path is None and dataframe is not None:
            new_rates_df = dataframe
        # Third case: None of them are given, raise an error
        elif reaclib_path is None and dataframe is None:
            raise Exception('Pass either a path of a reaclib or a pandas dataframe.')
        # Fourth case: Both of them are given, raise an error
        elif reaclib_path is not None and dataframe is not None:
            raise Exception('Too much input. Pass either reaclib_path or pandas dataframe.')

        # Filter for the correct reactions that will be merged
        new_rates_df= self.__filter_rates(new_rates_df,reaction_type,chapter)

        # Replace a specific label?
        if not replace_specific_label is None:
            oldlabel = replace_specific_label[0]
            newlabel = replace_specific_label[1]
            # Do not ignore labels for this
            ignore_label = False
            # self.__pd_reaclib.loc[self.__pd_reaclib['label'] == oldlabel,'label'] = 'totaldummylabel' # Something stupid that will never be in the reaclib
            # new_rates_df.loc[new_rates_df['label'] == newlabel,'label']           = 'totaldummylabel' # Same for the other dataframe
            # print(new_rates_df[new_rates_df['label'] == 'totaldummylabel'])

        # Pandas is really bad in dealing with lists as entry, so create a unique string for categorizing. (he4+c12 -> "chapterhe4c12label")
        reactant_string = new_rates_df['reactants'].apply(lambda x: ''.join(x))
        product_string  = new_rates_df['products'].apply(lambda x: ''.join(x))
        # Create a unique string for the reaction
        new_rates_df['sort_string']  = new_rates_df['chapter'].astype(str) + reactant_string + product_string

        if not ignore_label:
            if not replace_specific_label is None:
                new_rates_df.loc[new_rates_df['label']==newlabel,'sort_string'] +='totaldummylabel'
                new_rates_df.loc[new_rates_df['label']!=newlabel,'sort_string'] +='dontreplacethis'
            else:
                new_rates_df['sort_string'] += new_rates_df['label']
        if not ignore_reverse:
            new_rates_df['sort_string'] += new_rates_df['reverse'].astype(str)
        if not ignore_type:
            new_rates_df['sort_string'] += new_rates_df['type'].astype(str)
        # "Flag" the origin of the rates
        new_rates_df['version'] = 'new'


        # Pandas is really bad in dealing with lists as entry, so create a unique string for categorizing. (he4+c12 -> "chapterhe4c12label")
        reactant_string = self.__pd_reaclib['reactants'].apply(lambda x: ''.join(x))
        product_string  = self.__pd_reaclib['products'].apply(lambda x: ''.join(x))


        # Create also a unique string for the the old reactions
        self.__pd_reaclib['sort_string']  = self.__pd_reaclib['chapter'].astype(str) + reactant_string + product_string
        if not ignore_label:
            if not replace_specific_label is None:
                self.__pd_reaclib.loc[self.__pd_reaclib['label']==oldlabel,'sort_string']+= 'totaldummylabel'
            else:
                self.__pd_reaclib['sort_string'] += self.__pd_reaclib['label']
        if not ignore_reverse:
            self.__pd_reaclib['sort_string'] += self.__pd_reaclib['reverse'].astype(str)
        if not ignore_type:
            self.__pd_reaclib['sort_string'] += self.__pd_reaclib['type'].astype(str)

        self.__pd_reaclib['version'] = 'old'
        # Merge the dataframes
        frames = [self.__pd_reaclib,new_rates_df]
        merged_frame = pd.concat(frames)

        # Group them by the created string, create a new string like "oldoldnewold" and categorize that string
        category =  merged_frame.groupby(merged_frame['sort_string'])['version'].agg(lambda col: ''.join(col)).apply(categorize)
        # Join this result to the original frame
        merged_frame = merged_frame.join(category, on='sort_string',rsuffix='_cat')
        # Filter for the correct reactions
        self.__pd_reaclib = merged_frame[((merged_frame['version'] == 'old') & (merged_frame['version_cat'] == 3)) | ((merged_frame['version'] == 'new') & (merged_frame['version_cat'] == 1))]

        if not self.__quiet:
            new_rates_count = self.__analyze_update(self.__pd_reaclib)
            print('Updated '+str(new_rates_count)+' rates.')
            print('Updating done!')


        if not replace_specific_label is None:
            self.__pd_reaclib.loc[self.__pd_reaclib['label'] == 'totaldummylabel','label'] = newlabel

        # Make sure to remove everything again. Remove the version and sort_string again
        del self.__pd_reaclib['version']
        del self.__pd_reaclib['version_cat']
        del self.__pd_reaclib['sort_string']



    def add_new_rates(self,reaclib_path=None,dataframe=None,reaction_type=None,chapter=None,ignore_label=True,ignore_reverse=True,ignore_type=False):
        """
          Add new rates to existing ones
          Input:
            reaclib_path          - path of the reaclib that contains new rates
            dataframe             - pandas dataframe that contains reactions. If this is given, you don't have to give a reaclib_path.
            reaction_type         - Reaction type that will be updated (string or list of strings), for possible values see ".get_statistics()"
            chapter               - Chapter that will be updated (list of integer or integer)
            ignore_label          - ignore the label for updating and only look at the reactions itself
            ignore_reverse        - ignore if the reaction is reverse
            ignore_type           - ignore the type
            replace_specific_label- Should one specific label be replaced. E.g. Moeller 93 with Moeller 03 rates. For a list of the Jina labels see https://groups.nscl.msu.edu/jina/reaclib/db/labels.php
                                    if yes, the input should be a list with [old label, new label]
        """
        # Function to categorize 3 cases
        def categorize(x):
            # First case: There are new and old rates for the label (-> take the new ones)
            if ('new' in x) and ('old' in x):
                return 1
            # Second case: There are only new rates (-> don't take these rates)
            if ('new' in x) and not ('old' in x):
                return 2
            # Third case: There are only old rates (-> take all old rates)
            if not ('new' in x) and ('old' in x):
                return 3

        if not self.__quiet:
            print('Updating reaclib.')

        # Get the input values
        # First case: reaclib path is given, so use this as new dataframe
        if reaclib_path is not None and dataframe is None:
            # Reading other reaclib to dataframe and bring it to correct shape
            new_rates = reaclib(reaclib_path,quiet=True)
            new_rates.read_reaclib()
            new_rates_df = new_rates.get_dataframe()
        # Second case: dataframe is given, so use these reaction rates
        elif reaclib_path is None and dataframe is not None:
            new_rates_df = dataframe
        # Third case: None of them are given, raise an error
        elif reaclib_path is None and dataframe is None:
            raise Exception('Pass either a path of a reaclib or a pandas dataframe.')
        # Fourth case: Both of them are given, raise an error
        elif reaclib_path is not None and dataframe is not None:
            raise Exception('Too much input. Pass either reaclib_path or pandas dataframe.')

        # Filter for the correct reactions that will be merged
        new_rates_df= self.__filter_rates(new_rates_df,reaction_type,chapter)

        # Pandas is really bad in dealing with lists as entry, so create a unique string for categorizing. (he4+c12 -> "chapterhe4c12label")
        reactant_string = new_rates_df['reactants'].apply(lambda x: ''.join(x))
        product_string  = new_rates_df['products'].apply(lambda x: ''.join(x))
        # Create a unique string for the reaction
        new_rates_df['sort_string']  = new_rates_df['chapter'].astype(str) + reactant_string + product_string

        if not ignore_label:
            new_rates_df['sort_string'] += new_rates_df['label']
        if not ignore_reverse:
            new_rates_df['sort_string'] += new_rates_df['reverse'].astype(str)
        if not ignore_type:
            new_rates_df['sort_string'] += new_rates_df['type'].astype(str)
        # "Flag" the origin of the rates
        new_rates_df['version'] = 'new'

        # Pandas is really bad in dealing with lists as entry, so create a unique string for categorizing. (he4+c12 -> "chapterhe4c12label")
        reactant_string = self.__pd_reaclib['reactants'].apply(lambda x: ''.join(x))
        product_string  = self.__pd_reaclib['products'].apply(lambda x: ''.join(x))

        # Create also a unique string for the the old reactions
        self.__pd_reaclib['sort_string']  = self.__pd_reaclib['chapter'].astype(str) + reactant_string + product_string
        if not ignore_label:
            self.__pd_reaclib['sort_string'] += self.__pd_reaclib['label']
        if not ignore_reverse:
            self.__pd_reaclib['sort_string'] += self.__pd_reaclib['reverse'].astype(str)
        if not ignore_type:
            self.__pd_reaclib['sort_string'] += self.__pd_reaclib['type'].astype(str)

        self.__pd_reaclib['version'] = 'old'
        # Merge the dataframes
        frames = [self.__pd_reaclib,new_rates_df]
        merged_frame = pd.concat(frames)

        # Group them by the created string, create a new string like "oldoldnewold" and categorize that string
        category =  merged_frame.groupby(merged_frame['sort_string'])['version'].agg(lambda col: ''.join(col)).apply(categorize)
        # Join this result to the original frame
        merged_frame = merged_frame.join(category, on='sort_string',rsuffix='_cat')
        # Filter for the correct reactions
        self.__pd_reaclib = merged_frame[((merged_frame['version'] == 'old') & (merged_frame['version_cat'] == 3)) | ((merged_frame['version'] == 'new') & (merged_frame['version_cat'] == 2))]

        if not self.__quiet:
            new_rates_count = self.__analyze_update(self.__pd_reaclib)
            print('Added '+str(new_rates_count)+' rates.')
            print('Adding done!')


        # Make sure to remove everything again. Remove the version and sort_string again
        del self.__pd_reaclib['version']
        del self.__pd_reaclib['version_cat']
        del self.__pd_reaclib['sort_string']





    def add_rate(self,dataframe):
        """
          Add a given dataframe to the reaclib
        """
        self.__pd_reaclib = pd.concat([self.__pd_reaclib, dataframe], axis=0)


    def filter_with_sunet(self,nucleus_names):
        """
          Filter the rates with contained nuclei.
          Input:
            Contained nuclei of the network
        """
        def contained(a_list,ref_list=None):
            """
            Check if all elements of a list(a_list) are contained in another list(ref_list)
            """
            if ref_list is None:
                ref_list = nucleus_names

            for element in a_list:
                if element not in ref_list:
                    return False
            return True


        # Get the dataframe locally
        pd_dataframe = self.__pd_reaclib

        # Check that the rates are contained in sunet
        pd_dataframe['filter_reactants'] = pd_dataframe['reactants'].apply(contained)
        pd_dataframe['filter_products']  = pd_dataframe['products'].apply(contained)
        # Check if the reaction should be contained for our given sunet
        pd_dataframe['keep'] = list(map(lambda f_reac, f_prods: f_reac and f_prods, pd_dataframe['filter_reactants'], pd_dataframe['filter_products'] ))
        pd_dataframe = pd_dataframe[pd_dataframe['keep'] == True]

        # Delete the additional columns again
        del pd_dataframe['filter_reactants']
        del pd_dataframe['filter_products']
        del pd_dataframe['keep']

        # Set the rates
        self.__pd_reaclib = pd_dataframe





    def get_label_statistic(self):
        """
          Get a statistic which labels are involved in the database. An overview of all labels an their meaning can be found
          at https://groups.nscl.msu.edu/jina/reaclib/db/labels.php
        """
        # Count the occurance of the types and append the amount of reactions
        return (pd.value_counts(self.__pd_reaclib['label'].values).append(pd.Series({'Total reactions' :self.__pd_reaclib.count()[0]})))


if __name__ == '__main__':

    print('Test reaclib class')
    r = reaclib('../testfiles/reaclib_testfiles/reaclib.Christian.121204_beta_sf',quiet=False)
    r.read_reaclib()
    print(r.get_dataframe())
    print(r.get_statistics())
    # print(r.get_dataframe())
    # r.save_reaclib('unchanged2.dat',r.get_critical_low_temperature_rates())
    # r.get_all_reactiontypes()
    # print(r.get_critical_low_temperature_rates())
    # r.update('../testfiles/reaclib_testfiles/20161205ReaclibV2.0no910',reaction_type=None,ignore_label=False,ignore_reverse=False,ignore_type=False)
    # r.save_reaclib('all_ignore_none.dat')
    # r.test_reaclib()
    # r.drop_errors()
    # r.get_rate_error_html('reaclib_errors.html')
