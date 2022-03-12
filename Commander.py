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
SimControllerMaster.params["simtime"] = 5000

SimControllerMaster.params["recip"] = 1
SimControllerMaster.params["k"] = 10
SimControllerMaster.params["n"] = 1000
SimControllerMaster.params["h"] = 0.01

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


max_k = 500
kvals= np.array(list(set([int(i+0.5) for i in np.geomspace(2,max_k,max_k-1)])),dtype=int)
stride = int(sys.argv[1])
core = int(sys.argv[2]) +0.5 #+0.5 so that it rounds instead of truncates
kvals = kvals[np.arange(int(core),len(kvals),stride,dtype=int)]

for p in pparams2:
    SimControllerActual = SimControl(SimControllerMaster) #copy constructor

    SimControllerActual.params["p"] = p
    SimControllerActual.params["weight"] = 0
    SimControllerActual.params["GSweight"] = 0
    SimControllerActual.params["GGweight"] = 0
    
    os.system(SimControllerActual.Generate_Command())



t1 = time.time()

print("Total Time taken {:.3f}".format(t1-t0))


