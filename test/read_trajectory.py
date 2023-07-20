#!/usr/bin/env python

# Choose the tolerance relatively large, because it the comparison is a linear interpolation!
t.checklist = { \
   'mainout.dat'  : { 'method':'listcompare', 'tolerance':[1e-10,1e-10] ,'x_column':1, 'y_column':[2,3]}, \
}
t.program = t.basedir + "/bin/winnet"
t.testdir = t.basedir + "/test/" + t.testname
t.logfile = t.testdir + ".log"
t.arcfile = t.testdir + ".tar.gz"
