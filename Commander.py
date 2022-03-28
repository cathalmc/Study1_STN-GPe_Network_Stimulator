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
#nt= SimControllerMaster.params["Network_type"]
#print(f"Using network: {nt}")

SimControllerMaster.params["recip"] = 1
SimControllerMaster.params["n"] = 500
SimControllerMaster.params["k"] = 10
SimControllerMaster.params["h"] = 0.02

ps = {"Spatial": 0.01,
      "Small_world":0.01 ,
      "Scale_free": 1, 
      "SBlock": 1-0.01,
      "Regular": 0,}
      
#SimControllerMaster.params["p"] = ps[nt]    

run_list=np.geomspace(0.01,2,100)

stride = int(sys.argv[1])
core = int(sys.argv[2])

torun = np.array(run_list)[np.arange(int(core),len(run_list),stride,dtype=int)]

for d in torun:
    SimControllerActual = SimControl(SimControllerMaster) #copy constructor
    
    SimControllerActual.params["p"] = d
    
    SimControllerActual.params["weight"] = 0
    SimControllerActual.params["GSweight"] = 0
    SimControllerActual.params["GGweight"] = 0
    runcommand = SimControllerActual.Generate_Command()
    #print(runcommand)
    os.system(runcommand)



t1 = time.time()

print("Total Time taken {:.3f}".format(t1-t0))

#max_k = 90
# replicates = 5
# stride = int(sys.argv[1])
# core = int(sys.argv[2])

# kvals= list(set([int(i+0.5) for i in np.geomspace(2,max_k,max_k-1)])) #geometrically spaced values
# kvals=np.array([kv for _ in range(replicates) for kv in kvals ],dtype=int) #add in replicates
# kvals = sorted(kvals[np.arange(int(core),len(kvals),stride,dtype=int)],reverse=core%2) #assign values to core

