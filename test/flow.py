#!/usr/bin/env python

t.checklist = { \
   'flow/flow_0010.dat'    : { 'method':'default', 'tolerance':2.0e-2 }, \
   'flow/flow_0020.dat'    : { 'method':'default', 'tolerance':2.0e-2 }, \
   'flow/flow_0030.dat'    : { 'method':'default', 'tolerance':2.0e-2 }, \
   'flow/flow_0040.dat'    : { 'method':'default', 'tolerance':2.0e-2 }, \
   'flow/flow_0050.dat'    : { 'method':'default', 'tolerance':2.0e-2 }, \
   'flow/flow_0060.dat'    : { 'method':'default', 'tolerance':2.0e-2 }, \
}
t.program = t.basedir + "/bin/winnet"
t.testdir = t.basedir + "/test/" + t.testname
t.logfile = t.testdir + ".log"
t.arcfile = t.testdir + ".tar.gz"
