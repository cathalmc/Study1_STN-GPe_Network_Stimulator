import os
import sys
import numpy as np
import time as time
import random
from SimulationControl import SimControl

SimControllerMaster = SimControl()
t0 = time.time()


SimControllerMaster.params["name"] = "SG_July_NormSingle" ##### MAKE SURE NAME IS A SINGLE CONTINUOUS STRING SO SUBMITTER DOESNT GET CONFUSED
SimControllerMaster.params["Network_type"] = sys.argv[3]
SimControllerMaster.params["simtime"] = 3000

SimControllerMaster.params["recip"] = 1
SimControllerMaster.params["k"] = 10
SimControllerMaster.params["n"] = 500
SimControllerMaster.params["h"] = 0.05

if sys.argv[3] == "Spatial":
    pparams = np.geomspace(0.01,1,50)
elif sys.argv[3] == "Small_world":
    raise Exception
elif sys.argv[3] == "Scale_free":
    raise Exception
elif sys.argv[3] == "SBlock":
    raise Exception
elif sys.argv[3] == "Regular":
    raise Exception
else:
    raise Exception


stride = int(sys.argv[1]) #number of cores being used to run the simulations
core = int(sys.argv[2]) #the core id (0 to n-1)
pparams2 = pparams[np.arange(int(core),len(pparams),stride,dtype=int)]


for p in pparams2:
    SimControllerActual = SimControl(SimControllerMaster) #copy constructor

    SimControllerActual.params["p"] = p
    SimControllerActual.params["weight"] = 0
    SimControllerActual.params["GSweight"] = 0
    SimControllerActual.params["GGweight"] = 0
    
    os.system(SimControllerActual.Generate_Command())



t1 = time.time()

print("Total Time taken {:.3f}".format(t1-t0))


