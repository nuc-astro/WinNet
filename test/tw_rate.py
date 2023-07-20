#!/usr/bin/env python

t.checklist = { \
   'mainout.dat'  : { 'method':'listcompare', 'tolerance':[5e-4,5e-4] ,'x_column':1, 'y_column':[6,7]}, \
   'finab.dat'    : { 'method':'default', 'tolerance':5e-6 }, \
   'finabsum.dat' : { 'method':'default', 'tolerance':5e-6 }, \
}
t.program = t.basedir + "/bin/winnet"
t.testdir = t.basedir + "/test/" + t.testname
t.logfile = t.testdir + ".log"
t.arcfile = t.testdir + ".tar.xz"
