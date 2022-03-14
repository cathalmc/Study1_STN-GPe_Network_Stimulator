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
nt= SimControllerMaster.params["Network_type"]
print(f"Using network: {nt}")

SimControllerMaster.params["recip"] = 1
SimControllerMaster.params["n"] = 500
SimControllerMaster.params["h"] = 0.02

if sys.argv[3] == "Spatial":
    SimControllerMaster.params["p"] = 0.01
elif sys.argv[3] == "Small_world":
    SimControllerMaster.params["p"] = 0.01
elif sys.argv[3] == "Scale_free":
    SimControllerMaster.params["p"] =  1
elif sys.argv[3] == "SBlock":
    SimControllerMaster.params["p"] = 1-0.01
elif sys.argv[3] == "Regular":
    SimControllerMaster.params["p"] = 0
else:
    raise Exception


max_k = 400
kvals= np.array(list(set([int(i+0.5) for i in np.geomspace(2,max_k,max_k-1)])),dtype=int)
stride = int(sys.argv[1])
core = int(sys.argv[2]) +0.5 #+0.5 so that it rounds instead of truncates
kvals = kvals[np.arange(int(core),len(kvals),stride,dtype=int)]

for k in kvals:
    SimControllerActual = SimControl(SimControllerMaster) #copy constructor

    SimControllerActual.params["k"] = k
    SimControllerActual.params["weight"] = 0
    SimControllerActual.params["GSweight"] = 0
    SimControllerActual.params["GGweight"] = 0
    runcommand = SimControllerActual.Generate_Command()
    print(runcommand)
    os.system(runcommand)



t1 = time.time()

print("Total Time taken {:.3f}".format(t1-t0))


