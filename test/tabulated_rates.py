#!/usr/bin/env python

t.checklist = { \
   'tracked_nuclei.dat'    : \
   { 'method':'listcompare', 'tolerance':[5e-4 for i in range(5)] ,\
   'x_column':0, 'y_column':[i+1 for i in range(5)] ,"lowerlimit":1e-6} \
}
t.program = t.basedir + "/bin/winnet"
t.testdir = t.basedir + "/test/" + t.testname
t.logfile = t.testdir + ".log"
t.arcfile = t.testdir + ".tar.xz"
