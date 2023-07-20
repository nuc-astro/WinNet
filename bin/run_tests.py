#!/usr/bin/env python
########################################################################
# A Python script to run one or more tests from the test/ directory
#
# 2014 O. Korobkin
########################################################################

import os, sys, optparse
import time, subprocess
from testcase_class import *

#--- check arguments and print help message if no arguments given ------
parse = optparse.OptionParser()
parse.set_usage("""
   ./run_tests.py [--keep|-k] [--help|-h] t1 [t2 [...]]

 - unpack and launch test[s] in the test/ directory.

Arguments:
  t1 [t2 [..]]: tests to run """)

#--- parse options and arguments ---------------------------------------
parse.add_option("-k","--keep",action="store_false", \
                 dest="cleanup_after_test", default=True, \
                 help="do not erase test directory after it succeeded")

(options, args) = parse.parse_args()

if len(args) < 1:
   print("ERROR: Wrong argument string")
   parse.print_help()
   sys.exit(1)

#--- source and base directories  --------------------------------------
maindir = os.getcwd()
program = maindir+"/bin/winnet"
testdir = maindir+"/test"
bindir  = maindir+"/bin"
if not os.path.isdir(testdir): # check if the test directory exists
    print("ERROR: test directory does not exist. Perhaps you are not running")
    print("       this script from the WINNET root directory with the ")
    print("       `make tests` command.")
    sys.exit(1)

#--- main loop      ----------------------------------------------------
print("Running tests:")
n_tests = len(args)
n_failed = 0
for i in range(len(args)):
   # skip options
   tc= args[i]
   if tc[0]=='-':
      if not (tc=="-k" or tc=="--keep"):
         print ("ERROR: unknown option.")
         parse.print_help()
      continue

   t = testcase(tc, maindir)
   nd = 3 if len(args)>99 else 2 if len(args)>9 else 1
   wt_format = "%" + str(nd) + "d/%" + str(nd) + "d"
   t.whichtest = wt_format % (i+1,len(args))
   t.len_spaces = 45 - len(tc)
   if not os.path.exists(testdir + "/" + tc + ".py"):
      print("ERROR: testcase "+tc+" does not exist!")
      n_failed += 1
      continue
   exec(compile(open(testdir + "/" + tc + ".py", "rb").read(), testdir + "/" + tc + ".py", 'exec'))

   res= t.deploy()
   if res!=0:
      print("ERROR: while unpacking test "+tc)
      n_failed += 1
      continue

   res= t.launch()
   if res!=0:
      print("ERROR: when launching test "+tc)
      n_failed += 1
      continue

   t.monitor()
   res= t.analyze()
   if res!=0:
      print("\r [" + t.whichtest + "] " + tc \
                  + ": .."+"."*t.len_spaces + "... [FAIL]    ")
      n_failed += 1
   else:
      print("\r [" + t.whichtest + "] " + tc \
                  + ": .."+"."*t.len_spaces + "... [ OK ]    ")
      if options.cleanup_after_test:
         t.cleanup() # only if passed

if n_failed>0:
   print(" - %d out of %d tests have failed." % (n_failed, n_tests))
else:
   print(" - all tests have passed successfully!")
