#!/usr/bin/env python

t.checklist = { \
   'finab.dat'    : { 'method':'exists' }, \
   'WinNet_data.h5' : { 'method':'exists' }, \
}
t.program = t.basedir + "/bin/winnet"
t.testdir = t.basedir + "/test/" + t.testname
t.logfile = t.testdir + ".log"
t.arcfile = t.testdir + ".tar.gz"
