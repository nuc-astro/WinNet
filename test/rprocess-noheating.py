#!/usr/bin/env python

t.checklist = { \
   'mainout.dat'  : { 'method':'default', 'tolerance':1e-5 }, \
   'finab.dat'    : { 'method':'default', 'tolerance':1e-4 }, \
   'finabsum.dat' : { 'method':'default', 'tolerance':1e-4 }, \
}
t.program = t.basedir + "/bin/winnet"
t.testdir = t.basedir + "/test/" + t.testname
t.logfile = t.testdir + ".log"
t.arcfile = t.testdir + ".tgz"


