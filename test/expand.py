#!/usr/bin/env python

t.checklist = {  \
   'mainout.dat'  : { 'method':'listcompare', 'tolerance':[1e-2,1e-2,1e-7] ,'x_column':1, 'y_column':[2,3,5]}, \
   'finab.dat'    : { 'method':'exists'}}
t.program = t.basedir + "/bin/winnet"
t.testdir = t.basedir + "/test/" + t.testname
t.logfile = t.testdir + ".log"
t.arcfile = t.testdir + ".tar.xz"
