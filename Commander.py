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
SimControllerMaster.params["simtime"] = 4000

try:
    tester = sys.argv[4]
except:
    tester = None
    
testrun=False
if tester=="tester":
    testrun=True
    print("Validating")

SimControllerMaster.params["recip"] = 1
SimControllerMaster.params["n"] = 100
#SimControllerMaster.params["k"] = 20
SimControllerMaster.params["h"] = 0.03

nitit=200  
ps = {"ImprovedSpatial": np.linspace(0.5,15,nitit),#np.linspace(0,8.5,1000),
    "Small_world":np.geomspace(1e-3,1,nitit),#np.geomspace(1e-4,1,1000) ,
        "Scale_free": np.linspace(1e-4,4,nitit), 
        "SBlock": 1-np.geomspace(1e-3,1,nitit),
        "Regular": 0,}
        
psC = {"ImprovedSpatial": 7,
    "Small_world":1e-3,
        "Scale_free": 1, 
        "SBlock": 1-(1e-2),
        "Regular": 0,}
        
replicates = 1

max_k = 70

kvals= list(set([int(i+0.5) for i in np.geomspace(2,max_k,max_k-1)])) #geometrically spaced values
torun2 = np.array([kv for _ in range(replicates) for kv in kvals ],dtype=int) #add in replicates
torun = []

for v in torun2:
    for net in [sys.argv[3]]:#,"Small_world","Scale_free", "SBlock","Regular"]:
        torun.append( {"net":net,"p":psC[net],"k":v}) 


stride = int(sys.argv[1])
core = int(sys.argv[2])
torun = np.array(torun)
torun = torun[np.arange(int(core),len(torun),stride,dtype=int)]


for d in torun:
    SimControllerActual = SimControl(SimControllerMaster) #copy constructor
    SimControllerActual.params["Network_type"] = d["net"]
    SimControllerActual.params["p"] = d["p"]
    SimControllerActual.params["k"] = d["k"]
    SimControllerActual.params["StimSites"] = 0
    
    SimControllerActual.params["StimAmplitude"] = 1 #80
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

#max_k=25
#kvals = list(set([int(i+0.5) for i in np.geomspace(2,max_k,max_k-1)])) #geometrically spaced values
#replicates = 4
#torun2=np.array([kv for _ in range(replicates) for kv in kvals]) 

