import sys
import os
here = os.path.realpath("."); path = sys.path.insert(0,here+"/scripts")
try:
    import test as T
except ImportError:
    print "\nError!!! test module can not be imported. Be sure that you're executing the test from the main directory, using make for that.\n"

## NOTE: For this test to succeed, it is necessary:
## 1. Either to set the EPS value (accuracy) to a very low value in RT_SolveMatrix.c
## 2. Or to deactivate the radiative energy incremental guess (multiply d[l] by 0.0
## in the loop of RT_GuessErad.c)

description1 = """Testing a restart with the fargoRT setup
Initial run is on 6 processors.\n"""
description2 = """Restarting the simulation with 4 processors.\n"""
RestartTest = T.GenericTest(testname = "RESTART_TESTRT",
                            flags1 = "SETUP=fargoRT UNITS=0 PARALLEL=1 FARGO_DISPLAY=NONE GPU=0",
                            launch1 = "mpirun -np 6 ./fargo3d -m",
                            description1 = description1,
                            flags2 = "SETUP=fargoRT PARALLEL=1 FARGO_DISPLAY=NONE GPU=0",
                            launch2 = "mpirun -quiet -np 4 ./fargo3d -S 1 -m",
                            description2 = description2,
                            parfile = "setups/fargoRT/fargoRT.par",
                            verbose = False,
                            plot=False,
                            field = "gasdens",
                            compact = True,
                            parameters = {'dt':1.0, 'ninterm':1 ,'ntot':3,
                                          'nx':30, 'ny':30, 'nz':15},
                            clean = True,
                            restore = True,
                            n = 2)

RestartTest.set_commands(command1 = "mkdir RESTART_TESTRT/test2; cp -r RESTART_TESTRT/test1/* RESTART_TESTRT/test2/",
                         command2 = None)
RestartTest.run()
