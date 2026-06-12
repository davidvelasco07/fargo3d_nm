import sys
import os
here = os.path.realpath("."); path = sys.path.insert(0,here+"/scripts")
try:
    import test as T
except ImportError:
    print "\nError!!! test module can not be imported." + \
        "Be sure that you are executing the test from the main directory\n"

description1 = """
Building the parallel version of the setup fargoRT,
and running it with one processor. 
"""

description2 = """
Running the setup fargoRT with four processors.
"""

flags = "SETUP=fargoRT PARALLEL=1 FARGO_DISPLAY=NONE GPU=0"
MpiTest = T.GenericTest(testname = "MPI_TEST_RT",
                        flags1 = flags,
                        launch1 = "mpirun -np 1 ./fargo3d -m",
                        description1 = description1,
                        flags2 = flags,
                        launch2 = "mpirun -np 4 ./fargo3d -m",
                        description2 = description2,
                        parfile = "setups/fargoRT/fargoRT.par",
                        verbose = False,
                        clean=True,
                        plot=False,
                        field = "gasdens",
                        compact = True,
                        parameters = {'dt':1.0, 'ninterm':1, 'ntot':1,
                                      'ny':20, 'nz':20, 'nx':20})
MpiTest.run()
