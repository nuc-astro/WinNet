#!/usr/bin/env python

t.checklist = { \
   'finab.dat'    : { 'method':'default', 'tolerance':1e-4 }, \
}
t.program = t.basedir + "/bin/winnet"
t.testdir = t.basedir + "/test/" + t.testname
t.logfile = t.testdir + ".log"
t.arcfile = t.testdir + ".tar.gz"
