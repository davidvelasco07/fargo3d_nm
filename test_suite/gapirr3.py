import sys
import os
here = os.path.realpath("."); path = sys.path.insert(0,here+"/scripts")
try:
    import test as T
except ImportError:
    print "\nError!!! test module can not be imported. Be sure that you're executing the test from the main directory, using make for that.\n"


description1 = """Testing a setup with irradiation and density floor 1: CGS vs scalefree
Initial run is on 2 CPU processes in CGS.\n"""
description2 = """Repeating the simulation with 8 CPU processes (scale free simulation)\n"""
RestartTest = T.GenericTest(testname = "GAPIRR3",
                            flags1 = "SETUP=gapirr UNITS=CGS  PARALLEL=1 FARGO_DISPLAY=NONE GPU=0 MPICUDA=0 RESCALE=1",
                            launch1 = "mpirun -np 2 ./fargo3d -m",
                            description1 = description1,
                            flags2 = "SETUP=gapirr UNITS=0 PARALLEL=1 FARGO_DISPLAY=NONE GPU=0 MPICUDA=0 RESCALE=0",
                            launch2 = "mpirun -np 8 ./fargo3d -m",
                            description2 = description2,
                            parfile = "setups/gapirr/gapirr.par",
                            verbose = True,
                            plot = False,
                            field = "gasdens",
                            compact = True,
                            parameters = {'dt':0.5, 'ninterm':1 ,'ntot':3,
                                          'nx':5, 'ny':30, 'nz':15},
                            clean = True,
                            restore = False,
                            n = 3)

RestartTest.set_commands(command1 = "mkdir GAPIRR3/test2; cp -r GAPIRR3/test1/* GAPIRR3/test2/; rm -f GAPIRR3/test2/gas*3.dat; make mrproper",
                         command2 = None)
RestartTest.run()
