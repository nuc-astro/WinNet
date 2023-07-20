#!/usr/bin/env python

# Choose the tolerance relatively large, because it the comparison is a linear interpolation!
t.checklist = { \
   'finab.dat'  : {  'method':'default', 'tolerance':1.0e-5 }, \
}
t.program = t.basedir + "/bin/winnet"
t.testdir = t.basedir + "/test/" + t.testname
t.logfile = t.testdir + ".log"
t.arcfile = t.testdir + ".tar.gz"
