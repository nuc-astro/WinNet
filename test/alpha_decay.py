#!/usr/bin/env python

t.checklist = { \
   'tracked_nuclei.dat'    : { 'method':'default', 'tolerance':2.0e-2,"lowerlimit":1e-10 }, \
}
t.program = t.basedir + "/bin/winnet"
t.testdir = t.basedir + "/test/" + t.testname
t.logfile = t.testdir + ".log"
t.arcfile = t.testdir + ".tar.xz"
