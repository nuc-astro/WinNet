#!/usr/bin/env python
t.checklist = { \
   'finab.dat'    : { 'method':'default', 'tolerance':2.0e-2,"lowerlimit":1e-9 }, \
}
t.program = t.basedir + "/bin/winnet"
t.testdir = t.basedir + "/test/" + t.testname
t.logfile = t.testdir + ".log"
t.arcfile = t.testdir + ".tar.gz"


def prepare_function():
    # Create a prepared network folder
    prepared_network_path = "../network_data"
    os.chdir(t.testdir+"/testrun")
    command= t.program+" "+t.parfile_name+" "+prepared_network_path+" >Prepare_OUT 2>Prepare_ERR \n"
    subprocess.call("ulimit -s unlimited\n"  + command + \
    """echo $! > PID
    if [ $? -ne 0 ]; then
      echo "Test simulation with PID #`cat PID` probably failed to prepare."
      echo "Check """+t.testdir+"""/testrun/ERR for errors"
    fi
    """, shell=True)

    # Delete the network data from the parameter file to see if it works only with the prepared network folder
    command = "sed -i '/^net_source/d' "+t.parfile_name+ '&& sed -i "/^reaclib_file/d" '+t.parfile_name +\
              '&& sed -i "/^isotopes_file/d" '+t.parfile_name
    subprocess.call(command, shell=True)
    return 0


t.prepare_function = prepare_function