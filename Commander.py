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
#nt= SimControllerMaster.params["Network_type"]
#print(f"Using network: {nt}")

SimControllerMaster.params["recip"] = 1
SimControllerMaster.params["n"] = 200
SimControllerMaster.params["k"] = 10
SimControllerMaster.params["h"] = 0.01

SimControllerMaster.params["p"] = 1

      
#SimControllerMaster.params["p"] = ps[nt]    
ps = {"ImprovedSpatial": np.linspace(0,6,1000),
    "Small_world":np.linspace(1e-3,1,1000) ,
        "Scale_free": np.linspace(1e-4,4,1000), 
        "SBlock": 1-np.linspace(1e-3,1,1000),
        "Regular": 0,}


run_list=ps[SimControllerMaster.params["Network_type"]]

torun=[]
for i in [5,10,15,20,25,30,40,60,100,120,140,150]:
    for j in [0.1,0.5,1,2,3,5]: #length of the negative current input in biphasic DBS
        for k in [0.02,0.1,0.2,0.5,1.0,2.0]:
            torun.append({"StimFrequency":i,"StimAmplitude":j,"StimSites":k})

#replicates = 20
#torun=np.array([kv for _ in range(replicates) for kv in run_list ],dtype=int) 

stride = int(sys.argv[1])
core = int(sys.argv[2])

torun = torun[np.arange(int(core),len(torun),stride,dtype=int)]


for d in torun:
    SimControllerActual = SimControl(SimControllerMaster) #copy constructor
    
    SimControllerActual.params["StimFrequency"] = d["StimFrequency"]
    SimControllerActual.params["StimAmplitude"] = d["StimAmplitude"]
    SimControllerActual.params["StimSites"] = d["StimSites"]
    
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

#ps = {"Spatial": 0.01,"Small_world":0.01 ,"Scale_free": 1, "SBlock": 1-0.01,"Regular": 0,}

