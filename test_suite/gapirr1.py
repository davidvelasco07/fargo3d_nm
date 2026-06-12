# Test used for backward compatibility nightly built
import sys
import os
here = os.path.realpath("."); path = sys.path.insert(0,here+"/scripts")
try:
    import test as T
except ImportError:
    print "\nError!!! test module can not be imported. Be sure that you're executing the test from the main directory, using make for that.\n"


description1 = """Testing a setup with irradiation and density floor 1: GPU vs CPU
Initial run is on 2 GPU processes.\n"""
description2 = """Repeating the simulation with 5 processors.\n"""
RestartTest = T.GenericTest(testname = "GAPIRR1",
                            flags1 = "SETUP=gapirr UNITS=0 PARALLEL=1 FARGO_DISPLAY=NONE GPU=1 MPICUDA=0",
                            launch1 = "mpirun -np 2 ./fargo3d -D 0 -m",
                            description1 = description1,
                            flags2 = "SETUP=gapirr PARALLEL=1 FARGO_DISPLAY=NONE GPU=0 UNITS=MKS RESCALE=1",
                            launch2 = "mpirun -np 5 ./fargo3d -m",
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

RestartTest.set_commands(command1 = "mkdir GAPIRR1/test2; cp -r GAPIRR1/test1/* GAPIRR1/test2/; rm -f GAPIRR1/test2/gas*3.dat; make mrproper",
                         command2 = None)
RestartTest.run()
