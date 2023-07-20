#!/usr/bin/env python

def neutron_decay(x):
    """
      Function to calculate the decay of a neutron
    """
    X_n0  = 0.01
    alpha = np.exp(-6.786182)
    n_theory = X_n0 * np.exp(-alpha*x)
    return n_theory




t.checklist = { \
   'mainout.dat'  : { 'method':'analyticcompare', 'tolerance':[5e-4,5e-4] ,'x_column':1, 'y_column':[6,7],'function':[neutron_decay, lambda x: 1. - neutron_decay(x)]}, \
   'finab.dat'    : { 'method':'default', 'tolerance':5e-4 }, \
   'finabsum.dat' : { 'method':'default', 'tolerance':5e-4 }, \
}
t.program = t.basedir + "/bin/winnet"
t.testdir = t.basedir + "/test/" + t.testname
t.logfile = t.testdir + ".log"
t.arcfile = t.testdir + ".tgz"
