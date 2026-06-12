import sys
import os
here = os.path.realpath("."); path = sys.path.insert(0,here+"/scripts")
try:
    import test as T
except ImportError:
    print "\nError!!! test module cannot be imported." + \
        "Be sure that you are executing the test from the main directory\n"

description1 = """
Building the scale free version of the setup planetesimalsRT
and running it on one processor. 
"""

description2 = """
Building the cgs version of the setup planetesimalsRT
and running it on one processor. 
"""

flags1 = "SETUP=planetesimalsRT PARALLEL=0 FARGO_DISPLAY=NONE GPU=0 RESCALE=0 UNITS=0 MPICUDA=0"
flags2 = "SETUP=planetesimalsRT PARALLEL=0 FARGO_DISPLAY=NONE GPU=0 RESCALE=1 UNITS=CGS MPICUDA=0"
DimTest = T.GenericTest(testname = "DIM_TEST_PRT",
                        flags1 = flags1,
                        launch1 = "./fargo3d -m",
                        description1 = description1,
                        flags2 = flags2,
                        launch2 = "./fargo3d -m",
                        description2 = description2,
                        parfile = "setups/planetesimalsRT/planetesimalsRT.par",
                        verbose = False,
                        clean=True,
                        plot=False,
                        restore=False,
                        field = "gasdens",
                        compact = True,
                        parameters = {'dt':1.0, 'ninterm':1, 'ntot':1,
                                      'nx':30, 'ny':30, 'nz':15})
DimTest.run()
