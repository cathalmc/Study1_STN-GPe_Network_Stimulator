#now including noise and realistic connections

import neuron
h = neuron.h

from Cortical_Basal_Ganglia_Cell_Classes import GP_Neuron_Type, STN_Neuron_Type

import pyNN.neuron as sim
from pyNN.utility import get_simulator, init_logging, normalized_filename
from pyNN.parameters import Sequence
from pyNN.random import RandomDistribution as rnd

import scipy.signal as sg
import numpy as np

import time as time

from datetime import datetime

import matplotlib as mpl
import matplotlib.pyplot as plt

from useful_functions_sonic import *
import networkx as nx
from scipy.stats import kurtosis
import random
from additional_functions import *

now = datetime.now().strftime("%d-%m-%Y %H:%M:%S")
t0 = time.time()
### Set up parameters
params = {  "k":float(sys.argv[1]),
                "p": float(sys.argv[2]),
                "simtime": float(sys.argv[3]), 
                "delay":float(sys.argv[4]),
                "recip":float(sys.argv[5]),
                "weight":float(sys.argv[6]),
                "name":str(sys.argv[7]),
                "STNbias":float(sys.argv[8]),
                "Dummy":int(float(sys.argv[9])),
                "GPebias":float(sys.argv[10]),
                "Network_type":sys.argv[11],
                "GSweight": float(sys.argv[12]),
                "GGweight": float(sys.argv[13]),
                "Strweight": float(sys.argv[14]),
                "CTXweight": float(sys.argv[15]),
                "convergence": int(sys.argv[16]),
                "n":int(sys.argv[17]),
                "h":float(sys.argv[18]),
                "Notebook":int(sys.argv[19]),
                "StimSites":int(sys.argv[20]),
                "StimAmplitude":float(sys.argv[21]),
                "StimFrequency":float(sys.argv[22])
                }
                
                
def cdb(x):
    if x%3==0:
        return 1
    elif x%3==1:
        return -1/cbf
    else:
        return 0

amp=params["StimAmplitude"]
fs=params["StimFrequency"]
t = 1000/fs

stimstrt = 100   
pw = 0.2 #pulse width
cbf = 10 #charge balance factor
num_to_stim = params["StimSites"]
simtime = params["simtime"]
stimstop=simtime


recip = params['recip']
k = int(2*params['k']/(1+recip))            
params["weight"] = SGw(k)
params["GSweight"] = GSw(k)
params["GGweight"] = GGw(k)

np.random.seed(params["Dummy"])

converg = params["convergence"]
p = params['p']

simtime = params['simtime']    # (ms)
delays = params['delay']*np.array([2,2,1]) #GPe delay is shorter

n= params['n'] 

dt         = params['h']           # (ms) #############################################!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! best 0.03, use 0.1 for quick
max_syn_delay  = 2*params['delay']           # (ms)
sim.setup(timestep=dt, max_delay=max_syn_delay)

if params["StimSites"]>0:
    print("Performing external stimulation")
    x = sorted(np.arange(stimstrt,stimstop,t).tolist()+np.arange(stimstrt+pw,stimstop,t).tolist()+np.arange(stimstrt+(pw*(1+cbf)),stimstop,t).tolist())
    y = [amp*cdb(i) for i in range(len(x))]
    cs3 = sim.StepCurrentSource(times=x, amplitudes = y)


if not(params["Network_type"]=="Stochastic_block"):
    k = int(k)

gen_functions = {   "Scale_free":SG_ScaleFree,
                    "Small_world": SG_SmallWorld,
                    "Spatial":SG_ExponentialSpatial,
                    "Regular":SG_Regular,
                    "SBlock":SG_SBlock,
                    "ImprovedSpatial":SG_SpatialImproved}


network_gen_function = gen_functions[params["Network_type"]]
STG_list,GTS_list,GTG_list,graph_measures = network_gen_function(n,k,p,params["recip"])


### Create the neurons
v_init = rnd('uniform', (-66.0, -56.0))
a = STN_Neuron_Type(bias_current=params["STNbias"])
b = GP_Neuron_Type(bias_current=params["GPebias"])

STN_cells = sim.Population(n, a ,initial_values={'v': v_init}, label="STN_cells")
GPe_cells = sim.Population(n, b ,initial_values={'v': v_init}, label="GPe_cells")
Striatal_Pop = sim.Population(converg*n,sim.SpikeSourcePoisson(rate=3) , label='Striatal Neuron Spike Source')
CTX_Pop = sim.Population(converg*n,sim.SpikeSourcePoisson(rate=10) , label='Cortical Neuron Spike Source')

### Establish connections
SG = sim.FromListConnector(STG_list, column_names=["weight"])
GS = sim.FromListConnector(GTS_list, column_names=["weight"])
GG = sim.FromListConnector(GTG_list, column_names=["weight"])
syn_STNGPe = sim.StaticSynapse(delay= delays[0])
syn_GPeSTN = sim.StaticSynapse( delay= delays[1])
syn_GPeGPe = sim.StaticSynapse( delay= delays[2])

    
syn_StriatalGPe = sim.StaticSynapse(weight=params['Strweight'] , delay=0)
syn_CTXSTN = sim.StaticSynapse(weight=params['CTXweight'] , delay=0)

prj_STNGPe = sim.Projection(STN_cells,GPe_cells, SG, syn_STNGPe, source='soma(0.5)', receptor_type='AMPA') #AMPA
prj_GPeSTN = sim.Projection(GPe_cells, STN_cells, GS, syn_GPeSTN, source='soma(0.5)', receptor_type='GABAa')
prj_GPeGPe = sim.Projection(GPe_cells,GPe_cells, GG, syn_GPeGPe, source='soma(0.5)', receptor_type='GABAa') #GABAa
prj_StriatalGPe = sim.Projection(Striatal_Pop, GPe_cells, sim.FixedNumberPreConnector(converg), 
                             syn_StriatalGPe, source='soma(0.5)', receptor_type='GABAa')
prj_CTXSTN = sim.Projection(CTX_Pop, STN_cells, sim.FixedNumberPreConnector(converg), 
                             syn_CTXSTN, source='soma(0.5)', receptor_type='AMPA') ##!!!!!!!!!!

STNNoise = [sim.NoisyCurrentSource(mean=0, stdev = 0.05, start=0,stop =simtime,dt=1.0) for i in range(n)] #was 20*dt which was usually .6ms
GPeNoise = [sim.NoisyCurrentSource(mean=0, stdev = 0.05, start=0,stop =simtime,dt=1.0) for i in range(n)]


for i,cell in enumerate(STN_cells):
    cell.inject(STNNoise[i])
    

for i,cell in enumerate(GPe_cells):
    cell.inject(GPeNoise[i])

if num_to_stim>0:
    #all_edges = update_index(STG_list,0,n,0) + update_index(GTS_list,n,0,0) + update_index(GTG_list,n,n,0)
    #to_stim = [i for i in get_best_nodes(all_edges,2*n,ind=1) if i<n]
    #to_stim = [i for i in get_central_nodes(all_edges) if i<n]
    #for i in range(num_to_stim):
    #    if i<len(to_stim)-1:
    #        STN_cells[to_stim[i]].inject(cs3)
    
    for cell in STN_cells[np.random.choice(n,num_to_stim,replace=False)]:
        cell.inject(cs3)
        
### Recording
STN_cells.record('spikes')
GPe_cells.record('spikes')
recordables = ['soma(0.5).v', 'AMPA.i','GABAa.i']
STN_cells.record(recordables)
GPe_cells.record(recordables)

### Run
sim.run(simtime)
t1 = time.time()

### Extract data
STN = STN_cells.get_data().segments[0]
GPe = GPe_cells.get_data().segments[0]

### Process data
data_dict = fill_dict(STN,GPe,dt,simtime)
for key in params:
    data_dict[key] = params[key] 

data_dict["When"] = now
data_dict["Time_taken"] = t1-t0

#update output dictionary with graph measures
for key in graph_measures:
    data_dict[key] = graph_measures[key]
    
dat = [data_dict]

def AccuSave(data_dict):
    eigs=sorted(data_dict["eigs"])
    lam2=eigs[1]
    lamN=eigs[-1]
    smallsave = {key:val for key,val in data_dict.items() if not(key=='eigs')}
    smallsave["lam2"]=lam2
    smallsave["lamN"]=lamN
    smallsave["R"]=lamN/lam2
    olddat = np.load('Acumulated_data.npy',allow_pickle=True)
    newdat = np.append(olddat,smallsave)
    np.save("Acumulated_data.npy",newdat)



if params["Notebook"] == 1:
    np.save(f"B:\SimOutput\\NotebookData.npy",np.array({"STN":STN,"GPe":GPe,"data":data_dict}))
    AccuSave(data_dict)

elif params["Notebook"] == 2:
    np.save(f"B:\SimOutput\\NotebookData.npy",np.array({"STN":STN,"GPe":GPe,"data":data_dict}))
    np.save(f"B:\SimOutput\\NotebookNetworkData.npy",np.array({"STG_list":STG_list,"GTS_list":GTS_list,"GTG_list":GTG_list,"graph_measures":graph_measures}))
    AccuSave(data_dict)
else:
    initname =  ''.join('_{:.3f}'.format(params[x]) if type(params[x])==float else '_{:}'.format(params[x]) for x in params)
    filename = 'Run_data_{:}_t{:.5f}.npy'.format(initname,float(t1-t0))
    np.save(filename, dat) 