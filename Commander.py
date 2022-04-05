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

try:
    tester = sys.argv[4]
except:
    tester = None
    
testrun=False
if tester=="tester":
    testrun=True
    print("Validating")

SimControllerMaster.params["recip"] = 1
SimControllerMaster.params["n"] = 500
SimControllerMaster.params["k"] = 10
SimControllerMaster.params["h"] = 0.02

#SimControllerMaster.params["p"] = ps[nt]    
ps = {"ImprovedSpatial": np.linspace(0,10,1000),
    "Small_world":np.geomspace(1e-4,1,1000) ,
        "Scale_free": np.linspace(1e-4,4,1000), 
        "SBlock": 1-np.geomspace(1e-4,1,1000),
        "Regular": 0,}
        
psC = {"ImprovedSpatial": 3,
    "Small_world":2e-3,
        "Scale_free": 1, 
        "SBlock": 1-(2e-2),
        "Regular": 0,}
        
SimControllerMaster.params["p"] = psC[sys.argv[3]]

#torun=ps[SimControllerMaster.params["Network_type"]]
torun=np.arange(0,501,5)

replicates = 10
torun2=np.array([kv for _ in range(replicates) for kv in torun ],dtype=int) 

torun = []
for v in torun2:
    for net in [sys.argv[3]]:#,"Small_world","Scale_free", "SBlock","Regular"]:
        torun.append( {"net":net,"p":psC[net],"tostim":v}) 

stride = int(sys.argv[1])
core = int(sys.argv[2])
torun = np.array(torun)
torun = torun[np.arange(int(core),len(torun),stride,dtype=int)]


for d in torun:
    SimControllerActual = SimControl(SimControllerMaster) #copy constructor
    SimControllerActual.params["Network_type"] = d["net"]
    SimControllerActual.params["p"] = d["p"]
    SimControllerActual.params["StimSites"] = d["tostim"]
    
    SimControllerActual.params["StimAmplitude"] = 15
    SimControllerActual.params["StimFrequency"] = 140

    SimControllerActual.params["weight"] = 0
    SimControllerActual.params["GSweight"] = 0
    SimControllerActual.params["GGweight"] = 0
    runcommand = SimControllerActual.Generate_Command()
    
    if testrun:
        print(d)
        #print(runcommand)
    else:
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

#for i in [5,10,15,20,25,30,40,60,100,120,140,150]:
#    for j in [0.1,0.5,1,2,3,5]: #length of the negative current input in biphasic DBS
#        for k in [0.02,0.1,0.2,0.5,1.0,2.0]:
#            torun.append({"StimFrequency":i,"StimAmplitude":j,"StimSites":k})
#torun=np.array(torun)

