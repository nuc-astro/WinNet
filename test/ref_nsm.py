#!/usr/bin/env python

t.checklist = { \
   'mainout.dat'  : { 'method':'default', 'tolerance':1e-15 }, \
   'finab.dat'    : { 'method':'default', 'tolerance':1e-15 }, \
   'finabsum.dat' : { 'method':'default', 'tolerance':1e-15 }, \
}
t.program = t.basedir + "/bin/winnet"
t.testdir = t.basedir + "/test/" + t.testname
t.logfile = t.testdir + ".log"
t.arcfile = t.testdir + ".tgz"


