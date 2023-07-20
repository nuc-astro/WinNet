#!/usr/bin/env python

import os, sys, subprocess, time
from scipy.interpolate import interp1d
import numpy as np

# helper function to run shell commands and grab output
def shell_command(cmdlist):
   p = subprocess.Popen(cmdlist, stdout=subprocess.PIPE,\
                                 stderr=subprocess.PIPE)
   out, err = p.communicate()
   return out


def compare_default (stest, strial, tol,l_limit, flog):
   """
    compare text data from two files
    - ignore all strings which start with '#'
    - use relative tolerance 'tol'
    - compare test file 'stest' against trial file 'strial'
    - flog is a log file
   """
   i1 = 0
   i2 = 0
   lines1 = 0
   lines2 = 0
   errcode = 0
   for i1 in range(len(strial)):
      s1 = strial[i1].strip()
      # skip all empty lines and lines which start with '#'
      if (len(s1)==0) or (s1[0]=='#'):
         continue

      # same for the test file
      i2min = i2
      for i2 in range(i2min,len(stest)):
         s2 = stest[i2].strip()
         if (len(s2)==0) or (s2[0]=='#'):
            continue
         break

      # if reached the end of second file
      if i2==len(stest):
         break

      # compare strings field-by-field
      fields1 = s1.split()
      fields2 = s2.split()
      if len(fields1)!=len(fields2):
         errcode += 1
      else:
         for i in range(len(fields1)):
            try:
               f1 = float(fields1[i])
               f2 = float(fields2[i])

                #Check for lower limit
               if ((abs(f1)<=l_limit) and (abs(f2)<=l_limit)):
                   continue
                #Check for tolerance
               if abs(f2-f1)/(abs(f1)+1e-20)>tol:
                  errcode += 1
                  break
            except ValueError:
               if fields1[i]!=fields2[i]:
                  errcode += 1
                  break

      lines1 += 1
      lines2 += 1
      i2 += 1

   if errcode!=0:
      flog.write("FAILED, files differ significantly in %d lines\n" % errcode)
   elif lines1!=lines2:
      flog.write("FAILED, files have different number of lines\n")
      errcode = 1000

   return errcode




def compare_lists( stest, strial, x_col, y_col, tol, l_limit, flog):
   """
    compare x- and y- data of two lists from two files via linear interpolation
    - ignore all strings which start with '#'
    - use relative tolerance 'tol'
    - compare test file 'stest' against trial file 'strial'
    - x_col is the index of the x-column, this column should be monotonic increasing
    - y_col is the index of the y-column
    - flog is a log file
   """
   i = 0
   errcode = 0
   # Convert them to lists
   if isinstance(y_col,int):
       y_col = [y_col]
   if isinstance(tol,int):
       tol = [tol]

   # Get the amount of y-cols to compare with. There are more than 1 y-column, but only one x-column
   amount_ycols = len(y_col)

   x_trial = []
   y_trial = [[] for i in range(amount_ycols)]

   #Create x- and y-list of trial file
   for i in range(len(strial)):
      s1 = strial[i].strip()
      # skip all empty lines and lines which start with '#'
      if (len(s1)==0) or (s1[0]=='#'):
         continue

      try:
          # Only one x-value
          x_trial.append(float(s1.split()[x_col]))
          # Multiple y-values
          for j in range(amount_ycols):
              y_trial[j].append(float(s1.split()[y_col[j]]))
      except:
          errcode = 1


   x_test = []
   y_test = [[] for i in range(amount_ycols)]
   # Create x- and y-list of test file
   for i in range(len(stest)):
      s1 = stest[i].strip()
      # skip all empty lines and lines which start with '#'
      if (len(s1)==0) or (s1[0]=='#'):
         continue

      try:
          # Only one x-value
          x_test.append(float(s1.split()[x_col]))
          # Multiple y-values
          for j in range(amount_ycols):
              y_test[j].append(float(s1.split()[y_col[j]]))
      except:
          errcode = 1


    # Get the first and the last index of the x-value list
   start_x = max(x_trial[0],x_test[0])
   end_x   = max(start_x,min(x_trial[-1],x_test[-1]))
   comp_ind = []
   for i in range(len(x_trial)):
       if (x_trial[i] >= start_x) and (x_trial[i] <= end_x):
           comp_ind.append(i)

   # Get start and end comparison index
   comp_ind_start = min(comp_ind)
   comp_ind_end   = max(comp_ind)

   # Cycle through all lists
   for i in range(amount_ycols):
      # Now we can compare the two lists
      # Take the x-values of trial as interpolation points
      try:
          test_function  = interp1d(x_test,y_test[i]  ,kind='linear')
          trial = np.array(y_trial[i][comp_ind_start:comp_ind_end])
          test  = np.array(test_function(x_trial[comp_ind_start:comp_ind_end]))
          diff_list = list(map(lambda x,y: abs(1.-abs(x)/(abs(y)+1e-20)),trial,test))
          # Lower limit
          diff_list =np.array(diff_list)
          diff_list[(test<l_limit) & (trial<l_limit)] = 0
      except:
          errcode=4

      # Check for empty list
      try:
          max_diff  = max(diff_list)
      except:
          max_diff = 0.
          errcode=3


      # Is it above the tolerance?
      if max_diff > tol[i]:
          errcode = 2

      # Write errors to logfile
      if errcode==1:
          flog.write("FAILED, error in reading columns\n")
      elif errcode==2:
          flog.write("FAILED, maximum difference %4.2e" % max_diff + " was above the tolerance of %4.2e"  % tol[i]+" in column "+str(y_col[i])+"\n")
      elif errcode==3:
          flog.write("FAILED, could not find a maximum difference"+"\n")
      elif errcode==4:
          flog.write("FAILED, interpolation failed"+"\n")

      if errcode != 0:
         # Set the errorcode to some undefined value
         errcode = 100

   return errcode


def compare_analytic( stest, strial, x_col, y_col, function, tol, flog):
   """
    compare x- and y- data of two lists from two files via linear interpolation
    - ignore all strings which start with '#'
    - use relative tolerance 'tol'
    - compare test file 'stest' against an analytic solution
    - x_col is the index of the x-column
    - y_col is the index of the y-column
    - function is a list of functions, representing the behavior of the specific columns
    - flog is a log file
   """
   i = 0
   errcode = 0
   # Convert them to lists
   if isinstance(y_col,int):
       y_col = [y_col]
   if isinstance(tol,int):
       tol = [tol]

   # Get the amount of y-cols to compare with. There are more than 1 y-column, but only one x-column
   amount_ycols = len(y_col)

   x_trial = []
   y_trial = [[] for i in range(amount_ycols)]
   max_diff = 0
   #Create x- and y-list of trial file
   for i in range(len(strial)):
      s1 = strial[i].strip()
      # skip all empty lines and lines which start with '#'
      if (len(s1)==0) or (s1[0]=='#'):
         continue

      try:
          # Only one x-value
          x_trial.append(float(s1.split()[x_col]))
          # Multiple y-values
          for j in range(amount_ycols):
              y_trial[j].append(float(s1.split()[y_col[j]]))
      except:
          errcode = 1


   x_test = []
   y_test = [[] for i in range(amount_ycols)]
   # Create x- and y-list of test file
   for i in range(len(stest)):
      s1 = stest[i].strip()
      # skip all empty lines and lines which start with '#'
      if (len(s1)==0) or (s1[0]=='#'):
         continue

      try:
          # Only one x-value
          x_test.append(float(s1.split()[x_col]))
          # Multiple y-values
          for j in range(amount_ycols):
              y_test[j].append(float(s1.split()[y_col[j]]))
      except:
          errcode = 1



   # Cycle through all lists
   for i in range(amount_ycols):
      # Now we can compare the analytic solution to all given columns
      # Get the actual analytic function
      current_f = function[i]
      # Calculate the analytic solution
      try:
          analytic_y = current_f(np.array(x_test))
      except:
          # Something wrong with the x-input or the input function
          errcode=2

      # Calculate the difference in percent
      diff_list = list(map(lambda x,y: abs(1.-abs(x)/(abs(y)+1e-20)),y_test[i],analytic_y))
      # Check for empty list
      try:
          max_diff  = max(diff_list)
      except:
          errcode=3

      # Is it above the tolerance?
      if max_diff > tol[i]:
          errcode = 4

      # Write errors to logfile
      if errcode==1:
          flog.write("FAILED, error in reading columns\n")
      elif errcode==2:
          flog.write("FAILED, could not calculate analytic solution for column "+str(y_col[i])+"\n")
      elif errcode==3:
          flog.write("FAILED, could not find a maximum difference"+"\n")
      elif errcode==4:
          flog.write("FAILED, maximum difference %4.2e" % max_diff + " was above the tolerance of %4.2e"  % tol[i]+" in column "+str(y_col[i])+"\n")

      if errcode != 0:
         # Set the errorcode to some undefined value
         errcode = 100

   return errcode


class testcase:
   """ Class for winnet testcases

       Attributes:
       - testname: the name of the testcase
       - basedir:  winnet base directory (where the Makefile is located)
       - program:  path to the program that is being tested
       - testdir:  test directory (normally $WINNET/test/<TESTNAME>
       - parfile:  name of the parameter file that is used
       - logfile:  where the test log will be written

       Methods:
       - deploy()
       - launch()
       - monitor()
       - analyze()
       - cleanup()
   """
   def __init__(self, tc, bdir):
      self.testname = tc
      self.basedir = bdir
      self.program = bdir + "/bin/winnet"

   def deploy(self):
      #--- unpack testcase ---------------------------------------------------
      if not hasattr(self,'testdir'):
         self.testdir = self.basedir + "/test/" + self.testname
      if os.path.exists(self.testdir): # delete if test directory already exists
         subprocess.call(["rm","-rf", self.testdir])
      if not hasattr(self,'logfile'):
         self.logfile = self.basedir+"/test/"+self.testname+".log"
      if os.path.exists(self.logfile):
         subprocess.call(["rm","-rf",self.logfile])
      if not hasattr(self,'arcfile'):
         self.arcfile = self.basedir + "/test/" + self.testname + ".tgz"
      if not os.path.exists(self.arcfile):
         print("ERROR: archive file " + self.arcfile + " is missing")
         return 404

      os.chdir(self.basedir + "/test")
      subprocess.call(["tar","axf",self.arcfile])
      if not os.path.isdir(self.testdir):
         print("ERROR: archive "+self.arcfile+" is corrupt. Exiting...")
         return 101

      #--- parameter file ---------------------------------------------------
      if not hasattr(self,'parfile'):
         # the unpacked testcase directory must contain one (and only one!)
         # parameter file, matching the template '*.par'
         os.chdir(self.testdir)
         out = os.listdir('./')
         parfile_count = 0
         for fname in out:

            if fname.endswith('.par'):
               self.parfile = fname
               parfile_count = parfile_count + 1
         if parfile_count != 1:
            print("ERROR: cannot determine parameter file for the test")
            return 103
      else: # check if the parfile exists
         if not os.path.exists(self.parfile):
            print("ERROR: specified paramter file "+self.parfile+" is missing")
            return 404
      self.parfile_name = os.path.basename(self.parfile)

      #--- prepare temporary simulation directory ----------------------------
      os.chdir(self.testdir)
      if os.path.exists("testrun"):
         subprocess.call(["rm","-rf","testrun"])
      os.mkdir("testrun")
      # substitute @WINNET@  -> basedir, @TESTDIR@ -> testdir (with UNIX sed)
      command = "sed 's%@WINNET@%" +self.basedir+"%g;" \
              +      "s%@TESTDIR@%"+self.testdir+"%g' "\
              + self.parfile+" > " \
              + self.testdir + "/testrun/" + self.parfile_name
      subprocess.call(command,shell=True)
      # create directoris for snapshots and flow files
      os.chdir("testrun")
      os.mkdir("snaps")
      os.mkdir("flow")

      #--- compile code if necessary -----------------------------------------
      if not os.path.exists(self.program):
         os.chdir(self.basedir)
         subprocess.call(["make","clean"])
         subprocess.call("make")

      #--- return successful termination code --------------------------------
      return 0


   def launch(self):
      #--- run the program, don't forget to save process ID  -----------------
      # (the following assumes that the shell is bash; might not work in more
      #  primitive shells; will not work with csh or tcsh)
      os.chdir(self.testdir+"/testrun")
      command= self.program+" "+self.parfile_name+" >OUT 2>ERR &\n"
      subprocess.call("ulimit -s unlimited\n"  + command + \
      """echo $! > PID
      if [ $? -ne 0 ]; then
        echo "Test simulation with PID #`cat PID` probably failed to launch."
        echo "Check """+self.testdir+"""/testrun/ERR for errors"
      fi
      """, shell=True)
      return 0

   def monitor(self):
      #--- monitor the progress ----------------------------------------------
      # must have OUT file in the test suite
      os.chdir(self.testdir+"/testrun")
      if not os.path.exists("../OUT"):
         print("\nWARNING: test suite doesn't have example OUT file!")
         return 0

      out = shell_command(["wc","-l","../OUT"])
      expected_lines = int(out.split()[0])
      current_lines = 0

      # obtain process ID
      out = shell_command(["cat","PID"])
      PID = int(out.split()[0])
      still_running = True
      i = 0
      while still_running:
         time.sleep(0.1)
         if i%2==1:
            out = shell_command(["ps","h",str(PID)])
            still_running = (len(out.split())>0)
            #
            out = shell_command(["wc","-l","OUT"])
            current_lines = int(out.split()[0])
         #
         n1 = self.len_spaces * current_lines / expected_lines
         n1 = int(min(n1, self.len_spaces))
         n2 = self.len_spaces - n1
         q  = 100.0 * current_lines / expected_lines
         q  = min(100.0, q)
         i = i + 1
         sys.stdout.write("\r [" + self.whichtest + "] " + self.testname \
                                 + ": .."+"."*n1  + "-/|\\"[i%4] +" "*n2   \
                                 + (" [%5.1f" % q) + "%]")
         sys.stdout.flush()

      return 0


   def analyze(self):
      os.chdir(self.testdir+"/testrun")
      errcode = 0
      flog = open(self.logfile,'w')
      if not hasattr(self, 'checklist'):
         flog.write("WARNING: no files defined for checking\n")
         return 1
      flog.write("Comparing files:\n")
      for outfile in list(self.checklist.keys()):
         flog.write(" - %-20s: " % outfile)


         # call specified method with given tolerance
         meth     = self.checklist[outfile]['method']
         # read trial file
         if os.path.exists('../'+outfile) and (meth != 'exists'):
            f_trial = open('../'+outfile)
            s_trial = f_trial.read().split('\n')
            f_trial.close()
         elif (meth == 'exists'):
            pass # No trial file needed
         else:
            errcode += 1
            flog.write("FAILED, missing trial file in the test archive\n")
            continue

         # read test file
         if os.path.exists(outfile):
            if (meth != 'exists'):
                f_test = open(outfile)
                s_test = f_test.read().split('\n')
                f_test.close()
         else:
            errcode += 1
            flog.write("FAILED, missing file\n")
            continue

        # Set the tolerance, not necessary for exists method
         if (meth != 'exists'):
             tol      = self.checklist[outfile]['tolerance']
         #Check for lower limit of the comparison
         if 'lowerlimit' in self.checklist[outfile]:
             l_limit  = self.checklist[outfile]['lowerlimit']
         else:
             l_limit = 0.   #Default value
         try:
             if (meth == 'default'):
                errc = compare_default (s_test, s_trial, tol,l_limit, flog)
                errcode += errc
             elif(meth == 'listcompare'):
                x_column  = self.checklist[outfile]['x_column']
                y_column  = self.checklist[outfile]['y_column']
                errc      = compare_lists (s_test, s_trial, x_column, y_column, tol, l_limit, flog)
                errcode  += errc
             elif(meth == 'analyticcompare'):
                x_column  = self.checklist[outfile]['x_column']
                y_column  = self.checklist[outfile]['y_column']
                functions = self.checklist[outfile]['function']
                errc      = compare_analytic (s_test, s_trial, x_column, y_column, functions, tol, flog)
                errcode  += errc
             elif(meth == 'exists'):
                errc = 0
             else: #Insert new methods here
                errcode += 1
                flog.write("FAILED, unknown file comparison method '%s'\n" % meth)
                continue
         except:
             errcode += 1
             errc = 10
             flog.write("FAILED, test method failed.\n")

         if errc==0:
            flog.write("OK\n")

      flog.close()
      return errcode

   def cleanup(self):
      if os.path.exists(self.logfile):
         subprocess.call(["rm","-rf",self.logfile])
      if os.path.exists(self.testdir):
         subprocess.call(["rm","-rf",self.testdir])

