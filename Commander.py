import os
import sys
import numpy as np
import time as time
import random
from SimulationControl import SimControl

SimControllerMaster = SimControl()
t0 = time.time()


SimControllerMaster.params["name"] = "SG_July_NormSingle" ##### MAKE SURE NAME IS A SINGLE CONTINUOUS STRING SO SUBMITTER DOESNT GET CONFUSED
#SimControllerMaster.params["Network_type"] = sys.argv[3]
SimControllerMaster.params["simtime"] = 5000
#nt= SimControllerMaster.params["Network_type"]
#print(f"Using network: {nt}")

SimControllerMaster.params["recip"] = 1
#SimControllerMaster.params["n"] = 100
SimControllerMaster.params["h"] = 0.02

ps = {"Spatial": 0.01,
      "Small_world":0.01 ,
      "Scale_free": 1, 
      "SBlock": 1-0.01,
      "Regular": 0,}
      
#SimControllerMaster.params["p"] = ps[nt]      

d3 = [{'Network_type': 'Scale_free', 'n': '1000', 'k': 10},
 {'Network_type': 'Scale_free', 'n': '1000', 'k': 11},
 {'Network_type': 'Scale_free', 'n': '1000', 'k': 13},
 {'Network_type': 'Scale_free', 'n': '1000', 'k': 14},
 {'Network_type': 'Scale_free', 'n': '1000', 'k': 15},
 {'Network_type': 'Scale_free', 'n': '1000', 'k': 16},
 {'Network_type': 'Scale_free', 'n': '1000', 'k': 18},
 {'Network_type': 'Scale_free', 'n': '1000', 'k': 19},
 {'Network_type': 'Scale_free', 'n': '1000', 'k': 20},
 {'Network_type': 'Scale_free', 'n': '1000', 'k': 21},
 {'Network_type': 'Scale_free', 'n': '1000', 'k': 3},
 {'Network_type': 'Scale_free', 'n': '1000', 'k': 23},
 {'Network_type': 'Scale_free', 'n': '1000', 'k': 24},
 {'Network_type': 'Scale_free', 'n': '1000', 'k': 25},
 {'Network_type': 'Scale_free', 'n': '1000', 'k': 26},
 {'Network_type': 'Scale_free', 'n': '1000', 'k': 28},
 {'Network_type': 'Scale_free', 'n': '1000', 'k': 29},
 {'Network_type': 'Scale_free', 'n': '1000', 'k': 30},
 {'Network_type': 'Scale_free', 'n': '1000', 'k': 31},
 {'Network_type': 'Scale_free', 'n': '1000', 'k': 33},
 {'Network_type': 'Scale_free', 'n': '1000', 'k': 34},
 {'Network_type': 'Scale_free', 'n': '1000', 'k': 35},
 {'Network_type': 'Scale_free', 'n': '1000', 'k': 36},
 {'Network_type': 'Scale_free', 'n': '1000', 'k': 38},
 {'Network_type': 'Scale_free', 'n': '1000', 'k': 39},
 {'Network_type': 'Scale_free', 'n': '1000', 'k': 40},
 {'Network_type': 'Scale_free', 'n': '1000', 'k': 4},
 {'Network_type': 'Scale_free', 'n': '1000', 'k': 41},
 {'Network_type': 'Scale_free', 'n': '1000', 'k': 43},
 {'Network_type': 'Scale_free', 'n': '1000', 'k': 5},
 {'Network_type': 'Scale_free', 'n': '1000', 'k': 438},
 {'Network_type': 'Scale_free', 'n': '1000', 'k': 44},
 {'Network_type': 'Scale_free', 'n': '1000', 'k': 443},
 {'Network_type': 'Scale_free', 'n': '1000', 'k': 448},
 {'Network_type': 'Scale_free', 'n': '1000', 'k': 45},
 {'Network_type': 'Scale_free', 'n': '1000', 'k': 453},
 {'Network_type': 'Scale_free', 'n': '1000', 'k': 46},
 {'Network_type': 'Scale_free', 'n': '1000', 'k': 463},
 {'Network_type': 'Scale_free', 'n': '1000', 'k': 468},
 {'Network_type': 'Scale_free', 'n': '1000', 'k': 473},
 {'Network_type': 'Scale_free', 'n': '1000', 'k': 478},
 {'Network_type': 'Scale_free', 'n': '1000', 'k': 48},
 {'Network_type': 'Scale_free', 'n': '1000', 'k': 489},
 {'Network_type': 'Scale_free', 'n': '1000', 'k': 49},
 {'Network_type': 'Scale_free', 'n': '1000', 'k': 494},
 {'Network_type': 'Scale_free', 'n': '1000', 'k': 500},
 {'Network_type': 'Scale_free', 'n': '1000', 'k': 50},
 {'Network_type': 'Scale_free', 'n': '1000', 'k': 51},
 {'Network_type': 'Scale_free', 'n': '1000', 'k': 53},
 {'Network_type': 'Scale_free', 'n': '1000', 'k': 54},
 {'Network_type': 'Scale_free', 'n': '1000', 'k': 55},
 {'Network_type': 'Scale_free', 'n': '1000', 'k': 56},
 {'Network_type': 'Scale_free', 'n': '1000', 'k': 58},
 {'Network_type': 'Scale_free', 'n': '1000', 'k': 59},
 {'Network_type': 'Scale_free', 'n': '1000', 'k': 60},
 {'Network_type': 'Scale_free', 'n': '1000', 'k': 6},
 {'Network_type': 'Scale_free', 'n': '1000', 'k': 61},
 {'Network_type': 'Scale_free', 'n': '1000', 'k': 63},
 {'Network_type': 'Scale_free', 'n': '1000', 'k': 64},
 {'Network_type': 'Scale_free', 'n': '1000', 'k': 65},
 {'Network_type': 'Scale_free', 'n': '1000', 'k': 66},
 {'Network_type': 'Scale_free', 'n': '1000', 'k': 68},
 {'Network_type': 'Scale_free', 'n': '1000', 'k': 69},
 {'Network_type': 'Scale_free', 'n': '1000', 'k': 70},
 {'Network_type': 'Scale_free', 'n': '1000', 'k': 73},
 {'Network_type': 'Scale_free', 'n': '1000', 'k': 8},
 {'Network_type': 'Scale_free', 'n': '1000', 'k': 9}]
 
stride = int(sys.argv[1])
core = int(sys.argv[2])

torun = d3[np.arange(int(core),len(d3),stride,dtype=int)]

for dct in torun:
    SimControllerActual = SimControl(SimControllerMaster) #copy constructor
    
    SimControllerActual.params["Network_type"] = dct["Network_type"]
    SimControllerActual.params["n"] = dct["n"]
    SimControllerActual.params["k"] = dct["k"]
    SimControllerActual.params["p"] = ps[dct["Network_type"]]
    
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

