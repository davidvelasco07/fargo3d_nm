import sys                                                                                               
import os                                                                                                
here = os.path.realpath("."); path = sys.path.insert(0,here+"/scripts")                                  
try:                                                                                                     
    import test as T                                                                                     
except ImportError:                                                                                      
    print "\nError!!! test module can not be imported. Be sure that you're executing the test from the main directory, using make for that.\n"

description1 = """Running p3dadi setup without ghostsX on the CPU.\n"""                                   
description2 = """Running p3dadi setup with ghostsX on the GPU.\n"""                                      
Ghostsp3dadiTest = T.GenericTest(testname = "p3dadi_Ghosts_TEST", 
                           flags1 = "SETUP=p3dadi GHOSTSX=0 PARALLEL=0 FARGO_DISPLAY=NONE GPU=0",         
                           launch1 = "./fargo3d -m",
                           description1 = description1,                                                  
                           flags2 = "SETUP=p3dadi GHOSTSX=1 PARALLEL=0 FARGO_DISPLAY=NONE GPU=1",                   
                           launch2 = "./fargo3d -D 0 -m",
                           description2 = description2,
                           parfile = "setups/p3dadi/p3dadi.par",                                           
                           verbose = True,
                           plot=False,
                           field = "gasdens",                                                            
                           compact = True,
                           parameters = {'dt':0.2, 'ninterm':1, 'ntot':1,                                
                                         'ny':64, 'nx':64, 'nz':10},
                           clean = True) 
Ghostsp3dadiTest.run()  

