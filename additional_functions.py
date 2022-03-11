import neuron
h = neuron.h
import random
import numpy.matlib
from scipy import linalg as LA
from scipy.io import loadmat
import copy

#import os
import sys

import pyNN.neuron as sim
from pyNN.utility import get_simulator, init_logging, normalized_filename
from pyNN.parameters import Sequence
from pyNN.random import RandomDistribution as rnd

import scipy.signal as sg
import scipy.signal as signal

import numpy as np

import time as time

from datetime import datetime

from Cortical_Basal_Ganglia_Cell_Classes import GP_Neuron_Type, STN_Neuron_Type

from useful_functions_sonic import *

import networkx as nx
from scipy.stats import kurtosis
from scipy.stats import poisson
from scipy.stats import entropy
from collections import Counter

from random import shuffle


def get_peak_freq(LFP, dt):
    X = abs(np.fft.fft(LFP-np.mean(LFP)))
    freq = np.linspace(0,1000/dt,len(X))

    Y = sg.savgol_filter(X, 2*int(1/dt)+1, 3, deriv = 0, mode = 'nearest')
    Y = Y[0:int(1000/dt)]

    S_Peak_freq = freq[np.nonzero([Y==np.max(Y)])[1][0]]
    S_Peak_val = np.max(Y)
    return S_Peak_freq,S_Peak_val


def KuramotoPhase(pop,simtime=5000):
    main_lis= [] # list of each spike train, converted to a numpy array 
    #(its a 2d array because later on the distance matrix function requires this)
    for sp in pop.spiketrains:
        s_lis = [ float(i) for i in sp if i>0] #add the condition that we only measure spikes after things settle
        main_lis.append(np.array(s_lis))
    phases = []
    #fast and clever way of interpolating a phase within the interspike interval
    #each spike is +2*pi larger than the previous
    #interpolate this at whatever resolution you want
    #the interpolated signal mod 2*pi gives the phase in radians as a function of time
    for x in main_lis:
        factor = np.arange(2*np.pi,2*np.pi*len(x)+1,2*np.pi) 
        inx = np.arange(0,simtime)
        try:
            y = np.interp(inx,x,factor)
        except:
            y = np.ones(len(inx))
        phases.append(y%(2*np.pi))
    complex_phases = np.array([np.exp(1j*np.array(i)) for i in phases])
    R = abs(complex_phases.mean(axis=0))
    
    #ignore first 500 and last 50 ms to give the system time to settle, and to remove edge artefacts
    return np.mean(R[500:-50])


    

def burst_indicator(trains):
    '''
    Counts the proportion of interspike intervals less than 10ms
    If this proportion is not == 0 or ==1 then there may be some bursting
    Not useful on its own
    '''
    inter_spike_intervals = [ np.diff([j for j in i if j>1000]) for i in trains]
    maximum_ISI = 15
    proportion = [np.mean([1 if i<maximum_ISI else 0 for i in j]) for j in inter_spike_intervals if len(j)>2]
    if len(proportion)==0:
        return 0
    else:
        return np.mean(proportion)

def Spike_CV(trains):
    inter_spike_intervals = [ np.diff([j for j in i if j>500]) for i in trains]
    mn = [np.mean(i) for i in inter_spike_intervals]
    sd = [np.std(i) for i in inter_spike_intervals]
    CVs = [i/j for i,j in zip(sd,mn)]
    
    return np.mean(CVs)

def YaoFuglevand(pop,dt):
    simtime = int(len(np.array(pop.filter(name='soma(0.5).v')[0]).mean(axis=1))*dt)
    
    #Add up all the binned spike counts and compare them to whats expected under asynchronous (Poisson random) conditions
    
    main_lis= [] # list of each spike train, converted to a numpy array 
    for sp in pop.spiketrains:
        if len(sp)>5:
            s_lis = [ float(i) for i in sp if (i>500)] #add the condition that we only measure spikes after things settle
            main_lis.append(np.array(s_lis))

    trains = []
    if len(main_lis) == 0:
        return 0
    factor = 2 #bin width in ms
    
    #bin the data
    for sp in main_lis:
        x = np.zeros(int((simtime-500)/factor))
        ind = sp/factor

        x[ind.astype(int)-int(500/factor)]=1
        trains.append(x)

    atrains = np.array(trains)
    
    #theres no need to sort it, but it makes for nice plotting
    y = sorted(atrains.sum(axis=0).tolist())
    
    NB = len(y)
    m = np.mean(y)

    x = Counter(y)
    dat = []
    for c in x:
        if c>1: #2 or more synchronous spikes
            cx = c*(c-1)/2 #account for permutations
            dat.append(x[c]*cx) #weighted count of spikes
            #print(c,x[c]*cx)
    
    #compare with the expected number of spikes under asynchronous conditions by
    odat=[]
    for c in range(int(max(x))):
        if c>1:
            cx = c*(c-1)/2
            Nx = NB*poisson.pmf(c,m)#*((m**c)/(np.math.factorial(c)))*np.exp(-m) #poisson distribution
            odat.append(Nx*cx)
    
    val = (np.sum(dat)-np.sum(odat))/np.sum(odat)
    
    if 0<val<1e4:
        return val
    else:
        return 0




def get_hist(de,bounds = None):
    #bins the data into size 4ms bins and returns the location of the bins and their heights
    if bounds == None:
        #a = np.histogram(de, bins =np.arange(np.percentile(de,5),np.percentile(de,95),2))
        a = np.histogram(de, bins =np.linspace(np.percentile(de,5),np.percentile(de,95),21))
    else:
        #a = np.histogram(de, bins =np.arange(bounds[0],bounds[1]+1,4))
        a = np.histogram(de, bins =np.linspace(bounds[0],bounds[1],21))
    # np.linspace(np.percentile(de,5),np.percentile(de,95),11)
    
    
    b = a[1][:-1] + np.diff(a[1])/2
    c = a[0]
    return b,c


def Chi_Synchrony(popdata):
    '''
    Takes in data from a PyNN population segment 

    Returns synchrony measure from the scholarpedia article: 
        i.e the ratio of population variance to individual neuron variance

    '''
    sig_vi = np.var(popdata,axis = 0).mean()
    sig_V = np.mean(popdata,axis=1).var()
    
    return sig_V/sig_vi

def get_peak_freq(LFP,dt):
    X = abs(np.fft.fft(LFP-np.mean(LFP))) #subtract minimum value to remove DC component
    freq = np.linspace(0,1000/dt,len(X)) #frequency vector in Hz (the 1000 is because dt is given in ms)
    
    #CHANGE, no more filtering
    #Y = sg.savgol_filter(X, 2*int(1/dt)+1, 3, deriv = 0, mode = 'nearest') #filter window is roughly 2ms long
    #Y = Y[0:int(1000/dt)]
    
    
    Y = X[0:int(1000/dt)]

    S_Peak_freq = freq[np.nonzero([Y==np.max(Y)])[1][0]] #find the frequency with largest amplitude
    S_Peak_val = np.max(Y) 
    return S_Peak_freq,S_Peak_val

def recreate_list(prj):
    #Gets the edge list from the neuron projection object
    pynnlist = prj.get(["weight"], format="list",gather =True)
    plist = [[x[0],x[1]] for x in pynnlist]  
    return plist
    

def PoweratFreq(sig,freq,fs=1000/0.03):
    #Calculate the power within +/- 4Hz of the freq for a given sig
    wid = 1.5
    try:
        assert freq-wid>0
    except:
        return 0
    #calculate the PSD with welches method
    f,Pxx = signal.welch(sig-np.min(sig),fs,nperseg=int(len(sig)/4)) 

    #get indices corresponding to frequencies in range freq-4<x<freq+4
    indices = np.nonzero(np.logical_and(f>freq-wid,f<freq+wid)) 
    
    #get the change in frequencies
    df = f[1]-f[0]
    
    #integrate (quadrature)
    sum_squares = sum([df*y for y in Pxx[indices]])
    
    return sum_squares


def fill_reciprocal(el):
    #efficient way of producing the edge list, that when combined with el, gives a reciprocally connected graph

    el = tuple(map(tuple, el))
    perm = [1,0] #swap first and second column
    idx = np.empty_like(perm)
    idx[perm] = np.arange(len(perm))
    el2 = np.array(el)[:, idx]

    et= tuple(map(tuple, el2))

    properedges = set(et)-set(el)
    return list(properedges)

    
def standalone_get_reciprocal(elist):
    el1 = np.array(elist) 
    #this is an efficient way to reorder columns in a numpy array
    perm = [1,0] #swap first and second column
    idx = np.empty_like(perm)
    idx[perm] = np.arange(len(perm))
    el2 = el1[:, idx]
    return el2.tolist()
    

    
def update_index(elist,offset0,offset1,inplace = True):
    if inplace:
        for e in elist:
            e[0] = e[0]+offset0
            e[1] = e[1]+offset1
    else:
        newlist = []
        for i,e in enumerate(elist):
            newlist.append( list(e) ) 
            newlist[i][0] = e[0]+offset0
            newlist[i][1] = e[1]+offset1
            
        return newlist  





def fast_ER3(n,frac,reciprocal = False):
    #will return either a directed or undirected edge list for an Erdos Reyni graph
    flag = 0
    itercount = 0
    no_edges = int(frac*(n*(n-1)/2))
    
    while (flag==0):
        if itercount>30:
            return None
        G = nx.gnm_random_graph(n, no_edges)
        flag = (list(nx.selfloop_edges(G)) == []) #dont want self loops, but dont care if G is disconnected
        #print("Iterating "+str(itercount))
        itercount+=1
        
    #el = list(G.edges())
    el = [list(i) for i in G.edges]
    
    if reciprocal:
        #the edgelist of an undirected network needs to be changed for a directed net by adding the reciprocal connection
        el1 = np.array(el) 
        #this is an efficient way to reorder columns in a numpy array
        perm = [1,0] #swap first and second column
        idx = np.empty_like(perm)
        idx[perm] = np.arange(len(perm))
        el2 = el1[:, idx]
        el3 = np.concatenate((el1,el2)).tolist()
        return el3
    else:
        return el

def fast_ER4(n,k):
    #fast_ER3 returns a edgelist designed for a undirected network, and as a result the degree distribution is wrong
    #here we fix that, by taking a random sample of edges from the fully connected graph
    #I could probably make this more efficient but its already very fast so I cant be arsed
    x = fast_ER3(n,1,True)
    z = random.sample(x,k*n)
    return z



def split(n):
    if n%2==1:
        print("dont use this function on its own")
        raise Exception
    return [int(i/2) if i%2==0 else int((n+i)/2) for i in range(n)]

def get_degree_sequeunce(n,m,init_size,reciprocal):
    #Get the degree sequence which gives the exact number of edges required
    #If you want to generate a net with a specific number of connections ,not every edge will have same degree
    
    #number of nodes (2n) times number of connections m, minus connections that are already there
    
    factor = 1 if reciprocal else 2
    if reciprocal:
        factor = 1
        existing_edges = init_size*(init_size-1)/2
    else:
        factor = 2
        existing_edges = (init_size)**2 
    
    num_edges = n*m - existing_edges #factor*
    valid_neurons = factor*(n-init_size)
    node_degrees = {i:int(num_edges/valid_neurons) for i in range(factor*n) if not((i in range(n,n+init_size)) or (i in range(0,init_size)) )}

    h = 0
    for d in node_degrees:
        h+=node_degrees[d]
        
    #difference between desired number of edges and actual 
    offset = num_edges-h 

    count = 0
    flag = 0
    while not(flag):
        for d, v in node_degrees.items():
            if random.randint(0,1):
                node_degrees[d] = v+1
                count +=1
                if count>=offset:
                    flag = 1
                    break
    h=0
    for d in node_degrees:
        h+=node_degrees[d]
    

    
    return node_degrees
    
def scale_free_BA(n,k,alpha=1, reciprocal = False):
    #new and improved, fixes all the problems with the old method, and adds a bunch new ones (maybe)
    #m = m/2
    
    factor = 1 if reciprocal else 2
    initial_cluster_size = int(k)
    m = int(k/2) if reciprocal else k
    
    
    #Barabasi Albert non linear scale free network
    #for setting up the bipartite STN - GPe network
    #all nodes have equal out-degree, variable in-degree (change return statement to get opposite)




    #create the initial hub as a fully connected graph
    #el1 = fast_ER3(factor*initial_cluster_size,1) #,reciprocal = True)
    
    if reciprocal:
        edge_list = fast_ER3(initial_cluster_size,1,True) #el1 + standalone_get_reciprocal(el1)
        
    else:
        el2 = [(i,j+n) for i in range(initial_cluster_size) for j in range(initial_cluster_size)]
        S1, S2 = np.split(np.random.permutation(el2), [int(len(el2)/2)] )
        S1 = S1.tolist()
        S2 = standalone_get_reciprocal(S2)
        edge_list = S1+S2

    probs = {i:0 for i in range(factor*n)}


    #to satisfy requirment that there are n*m edges going from S to G and G to S, we need to adjust the degree of each node
    #this is necessary because having that initial fully connected block messes with the count
    per_node_deg_list = get_degree_sequeunce(n,m,initial_cluster_size,reciprocal)
    
    if reciprocal:
        noderng = range(factor*n)
    else:
        noderng = split(factor*n)

    for i in noderng:
        if i in range(0,initial_cluster_size):
            continue
            #skip over nodes whos index is in the initial cluster 
        if i in range(n,n+initial_cluster_size):
            continue

        degs = {q:0 for q in range(factor*n)}
        for e in edge_list: #this is expensive, might be a more efficient way
            degs[e[0]]+=1
            degs[e[1]] +=1 #in degrees

        total_deg = 0
        for d in degs:
            total_deg+= (degs[d])**alpha

        for p in probs:
            try:
                probs[p]= ((degs[p])**alpha)/total_deg
            except:
                print(i)

        #deg_list = [degs[d] for d in degs]
        prob_list = [probs[d] for d in probs]

        tempnpnodes = [q for q,val in enumerate(prob_list)]

        assert len(prob_list)==len(tempnpnodes)

        ##GPe consideration
        if reciprocal:
                prob_list[i] = 0
        #set it so that connections cannot be established within the same population
        else:
            if i<n:
                for b in range(0,n):
                    prob_list[b] = 0
            else:
                for b in range(n,2*n):
                    prob_list[b] = 0

        q=0
        mx = n-1
        for j in range(per_node_deg_list[i]):
            flag = 0

            target = random.choices(tempnpnodes, weights=prob_list, k=1)[0]
            #target = np.random.choice(tempnpnodes, 1, p=prob_list/np.sum(prob_list))[0]
            proposition = [i,target]

            while  (proposition in edge_list):
                prob_list[tempnpnodes.index(target)] = 0
                if sum(prob_list)==0:
                    #if probabilities sum to zero we move on to the next node, should only happen in first few nodes
                    print(f"Probabilities sum to zero: i = {i}, q = {q}, j = {j}")
                    flag = 1
                    break

                target = random.choices(tempnpnodes, weights=prob_list, k=1)[0]
                proposition = [i,target]

                q+=1
                if q>mx:
                    print('never satisfied')
                    raise Exception
            if flag == 1:
                break

            edge_list.append(proposition)
            if reciprocal:
                rev = [proposition[1],proposition[0]]
                edge_list.append(rev)

    #reverse
    return [[e[1],e[0]] for e in  edge_list]  




def get_degree(lst,ind,n,rtn=False):
    degree = {i:0 for i in range(n)}
    for e in lst:
        degree[e[ind]] +=1

    deglist = np.array([[i,degree[i]] for i in degree])

    if rtn:
        return deglist.tolist()
    else:
        plt.figure(figsize=(12,8))
        plt.bar(deglist[:,0],deglist[:,1] )
        




from scipy.sparse import dok_matrix
import scipy.sparse as spr



def get_partial_reciprocal(init_edges,n,recip=1,SF=False):
    #First make sure the indices are correct, we start with the STN to GPe edges
    #STN nodes are numbers less than n, GPe nodes are greater than or equal n
    #init_edges = update_index(init_edges,0,n,0) 
    
    if not(SF):
        init_edges = update_index(init_edges,0,n,0)
    
    #Split the STN to GPe edges randomly into two groups
    S1, S2 = np.split(np.random.permutation(init_edges), [int(len(init_edges)/2)] )
    #keep the first group as is
    S1 = S1.tolist()
    #invert the edges in the second group, these are thus GPe to STN edges
    S2 = standalone_get_reciprocal(S2)
    
    #we now have an STN-GPe network with non over lapping edges going from STN to GPe and GPe to STN 
    #Reciprocity is 0 and the mean degree is halved 
    unique_edges = S1+S2
    
    #calculate the total number of reciprocal edges to add
    recip_number = int(recip*len(unique_edges))
    #randomly select the appropriate number edges to be made reciprocal
    recip_edges = random.sample(unique_edges,k=recip_number)
    #add these reciprocal to the edge list and return
    new_el = unique_edges +  ([] if len(recip_edges)==0 else standalone_get_reciprocal(recip_edges))
    
    
    SG = [[e[0],e[1]-n] for e in new_el if e[0]<n]
    GS = [[e[0]-n,e[1]] for e in new_el if e[0]>=n]
    
    return SG, GS

def set_reciprocal(el,n,recip):
    '''
    Gets the unique edges of a list by looking at upper triangual form of adjacence matrix
    Then creates a new adjacency matrix by reversing half of the edges
    The makes each edege reciprocal with probability recip, by sampling
    
    '''
    R=[]
    for e in el:
        if e[1]>=e[0]:
            R.append((e[0],e[1]))
    unique_edges = list(set(R))

    S1, S2 = np.split(np.random.permutation(unique_edges), [int(len(unique_edges)/2)] )

    S1 = S1.tolist()
    S2 = standalone_get_reciprocal(S2)

    unique_edges = S1+S2
    recip_number = int(recip*len(unique_edges))
    recip_edges = random.sample(unique_edges,k=recip_number)
    new_el = unique_edges +  ([] if len(recip_edges)==0 else standalone_get_reciprocal(recip_edges))
    
    return new_el
    
def fastSW(n,k,p):
    #create an undirected SW network then add in the extra edges to make it directed, then return the edge list
    x = nx.connected_watts_strogatz_graph(n,k,p)
    el = []
    for ee in x.edges:
        e = list(ee)
        el.append([e[0],e[1]])
        el.append([e[1],e[0]])

    return el
   


def prob(x,cd):
    return np.exp(-x/cd)

def ExponentialSpatial(n,k,p,distxy):
    edges = set()
    for selected in range(n):
        row = [[i,val] for i,val in enumerate(distxy[selected])]
        indices = [i[0] for i in row ]
        distances  = [prob(i[1],p) for i in row ]
        targets = np.random.choice(indices, k, replace=False, p = distances/sum(distances))
        for t in targets:
            edges.add((selected,t))

    return list(edges)

def ExponentialSpatialReciprocal(n,k,p,distyy):
    edges = set()
    for selected in range(n):
        row = [[i,val] for i,val in enumerate(distyy[selected]) if not(i==selected)]
        indices = [i[0] for i in row ]
        distances  = [prob(i[1],p) for i in row ]
        targets = np.random.choice(indices, int(1.2*k), replace=False, p = distances/sum(distances))

        start = len(edges)
        for t in targets:
            edges.add((selected,t))
            edges.add((t,selected))
            end = len(edges)
            if (end-start)>=k:
                break
    return list(edges)



def get_spks(dat):
    spk_lis= [] # list of each spike train, converted to a numpy array 
    #(its a 2d array because later on the distance matrix function requires this)
    for sp in dat.spiketrains:
        s_lis = [ float(i) for i in sp] #add the condition that we only measure spikes after things settle
        spk_lis.append(np.array(s_lis))
    return spk_lis

def SBlock(n,k,p):
    b=4
    s = int(n/b)

    pi = min((p*k)/(s-1), 1)

    internal_edges = pi*b*s*(s-1)

    external_edges = n*k - internal_edges
    
    pe = max(0,external_edges/(s*s*b*(b-1)))

    sizes = [s,s,s,s]
    probs = [[pi, pe, pe, pe], 
             [pe, pi, pe, pe], 
             [pe, pe, pi, pe],
             [pe, pe, pe, pi]
            ]
    
    x =  nx.stochastic_block_model(sizes, probs,directed=True,selfloops=False)
    return [list(i) for i in x.edges]

def RP(x,vec):
    x = max(x,2)
    x=np.log10(x)
    a0,a1,a2,b1,b2 = vec
    return (a0+a1*x+a2*x*x)/(1+b1*x+b2*x*x)

def SGw(k):
    return RP(k,[ 0.12529176, -0.07874646,  0.01889676,  1.3292782 ,  0.155241  ])

def GSw(k):
    return RP(k,[ 2.4018365 ,  0.44080655,  0.05658033, -1.4868729 ,  1.68979823])

def GGw(k):
    return RP(k,[ 0.00738542, -0.00900335,  0.00296832, -1.47771821,  0.57017723])


def weight_list_pynn_direct(el,wf,n):
    weight = get_degree(el,1,n,rtn=True)
    return [ [e[0],e[1],wf(weight[e[1]][1])]  for e in el ]


def calculate_weigenvalues(elist,n):
    #got it right this time i swear
    #Motter 2005 calculation
    #In degree normalised laplacian
    #Calculate other graph stuff here for efficiency
    
    A = np.zeros((n,n))
    D = np.zeros((n,n))
    for c in elist:
        A[c[0],c[1]] = -1
        
    indegrees = [-np.sum(A[:,i]) for i in range(n)]

    for i in range(n):
        A[i,i]=indegrees[i]
        D[i,i]=(1/indegrees[i]) if indegrees[i]>0 else 0

    L = np.matmul(D,A)
    try:
        Lapeigenvalues = np.linalg.eigvals(L)
        Fr = (np.linalg.norm(L,'fro'))**2
        LapDeparture = np.sqrt(max(0,Fr - sum([abs(i)**2 for i in Lapeigenvalues])))
        condition_number = np.linalg.cond(L)
    except:
        Lapeigenvalues = "Fail"
        Fr = "Fail"
        LapDeparture = "Fail"
        condition_number = "Fail"
        print("Eig calculation failed")
    

    
    F = nx.DiGraph()
    F.add_edges_from( [(e[0],e[1]) for e in elist])
    
    graph_measures = {"mean_indegree":np.mean(indegrees),
                        "var_indegree":np.var(indegrees)}
    
    try:
        graph_measures["ASPL"] = nx.average_shortest_path_length(F)
    except:
        graph_measures["ASPL"] = 100
        print("ASPL Failed")
    try:
        graph_measures["diameter"] = nx.diameter(F)
    except:
        graph_measures["diameter"] = 100
        print("Diameter Failed")

    graph_measures["eigs"] = Lapeigenvalues
    graph_measures["LapDeparture"] = LapDeparture
    graph_measures["Condition"]=condition_number
    
    
    return graph_measures


def calc_network_measures(SG,GS,GG,n):
    #Make the edgelists weighted
    #These edge lists will be used in the pynn simulation
    STG_list = weight_list_pynn_direct(SG,SGw,n) 
    GTS_list = weight_list_pynn_direct(GS,GSw,n)
    GTG_list = weight_list_pynn_direct(GG,GGw,n)
    
    #copy the edge lists so they can be normalsied and used for network analysis
    lsg = copy.deepcopy(STG_list)
    lgs = copy.deepcopy(GTS_list)
    lgg = copy.deepcopy(GTG_list)

    all_normed_edges = update_index(lsg,0,n,0) + update_index(lgs,n,0,0) + update_index(lgg,n,n,0)
    
    graph_measures = calculate_weigenvalues(all_normed_edges,2*n)
    
    return STG_list,GTS_list,GTG_list, graph_measures
    

def SG_SBlock(n,k,p,r=1):

    STG_list,GTS_list = get_partial_reciprocal(SBlock(n,k,p),n,recip=r)
 
    GTG_list = set_reciprocal(SBlock(n,k,p),n,r)

    return calc_network_measures(STG_list,GTS_list,GTG_list,n)

def SG_ExponentialSpatial(n,k,p=0.1,r=1):
    ##ER random connections
    dim = 2
    a = np.random.uniform(0,1,size = [n,dim])
    b = np.random.uniform(0,1,size = [n,dim])
    
    distxy = distance_matrix(a,b,p=2) #distances from points x to points y
    distyy = distance_matrix(b,b,p=2) #distances between points y
    
    STG_list,GTS_list = get_partial_reciprocal(ExponentialSpatial(n,k,p,distxy),n,recip=r)
  
    GTG_list = set_reciprocal(ExponentialSpatialReciprocal(n,k,p,distyy),n,r)

    return calc_network_measures(STG_list,GTS_list,GTG_list,n)

def SG_Regular(n,k,p=None,r=None):
    r=1
    reg = nx.random_regular_graph(k,n)

    el  = []
    for e in reg.edges:
        el.append(tuple(e))
        el.append((e[1],e[0])  )
        
    STG_list,GTS_list = get_partial_reciprocal(el,n,recip=r)

    GTG_list = set_reciprocal(el,n,r)

    return calc_network_measures(STG_list,GTS_list,GTG_list,n)
    
def SG_ScaleFree(n,k,p,r=1):
    SGGS = scale_free_BA(n,k,alpha=p,reciprocal = False)
    STG_list,GTS_list = get_partial_reciprocal(SGGS,n,recip=r, SF=True)

    GTG_list = set_reciprocal(scale_free_BA(n,k,alpha=p, reciprocal = True),n,r)
    
    return calc_network_measures(STG_list,GTS_list,GTG_list,n)

   
def SG_SmallWorld(n,k,p,r=1):
    STG_list,GTS_list = get_partial_reciprocal(fastSW(n,k,p),n,recip=r)
    GTG_list = set_reciprocal(fastSW(n,k,p),n,r)

    return calc_network_measures(STG_list,GTS_list,GTG_list,n) 
    

def fill_dict(STNdata,GPedata,dt,simtime,currents = True):
    low_cutoff = int(500/dt) #number of segments corresponding to first 500ms

    SLFP = np.array(STNdata.filter(name='soma(0.5).v')[0]).mean(axis=1)#(axis=1)
    GLFP = np.array(GPedata.filter(name='soma(0.5).v')[0]).mean(axis=1)
    cut_SLFP = SLFP[low_cutoff::]
    cut_GLFP = GLFP[low_cutoff::]
    
    #Largest frequencies in the power spectrum
    S_peak_freq, S_peak_val = get_peak_freq(cut_SLFP,dt)
    G_peak_freq, G_peak_val = get_peak_freq(cut_GLFP,dt)

    #Spike frequencies: mean and coefficient of variation
    #STNSpikes = np.array(medspkfreqs(STNdata.spiketrains))
    Smean = meanCVspkfreq(STNdata.spiketrains,simtime)

    #GPeSpikes = np.array(medspkfreqs(GPedata.spiketrains))
    Gmean = meanCVspkfreq(GPedata.spiketrains,simtime)
    
    nps = int((len(cut_GLFP)-low_cutoff)/4)
    freq, csd = sg.csd(P_SLFP,P_GLFP, fs =1000/dt,nperseg = nps)
    
    coh = (abs(csd)**2)
    ang_csd = np.angle(csd)
    max_coh_ind = np.argmax(coh)
    freq_at_max = freq[max_coh_ind]
    max_coh_val = coh[max_coh_ind]
    phase_at_max = ang_csd[max_coh_ind]
    phase_ms = 1000*phase_at_max/(2*freq_at_max*np.pi)


    Stats = {
             "Peak STN Freq": (S_peak_freq, S_peak_val),
             "Peak GPe Freq": (G_peak_freq, G_peak_val),
             "STN synchrony": Chi_Synchrony(STNdata.filter(name='soma(0.5).v')[0]),
             "GPe synchrony": Chi_Synchrony(GPedata.filter(name='soma(0.5).v')[0]),
             "SMean": np.mean(Smean),
             "GMean": np.mean(Gmean),
             "STNk":KuramotoPhase(STNdata,simtime),
            "GPek":KuramotoPhase(GPedata, simtime),
            "SISIstd":Spike_CV(STNdata.spiketrains),
            "GISIstd":Spike_CV(GPedata.spiketrains),
             "Max CSD": max(abs(csd)),
             "Max COH": max_coh_val,
             "Max COH Freq":freq_at_max,
             "Max COH Phase": phase_at_max,
             "STNyf":YaoFuglevand(STNdata,dt),
            "GPeyf":YaoFuglevand(GPedata,dt),
             

            }
    if currents:
        #measures calculated from synpatic potentials
        SAMPA = np.array(STNdata.filter(name='AMPA.i')[0]).mean(axis=1).mean()
        GAMPA = np.array(GPedata.filter(name='AMPA.i')[0]).mean(axis=1).mean() 
        
        SGABA = np.array(STNdata.filter(name='GABAa.i')[0]).mean(axis=1).mean()
        GGABA = np.array(GPedata.filter(name='GABAa.i')[0]).mean(axis=1).mean()
         
        current_stats = {
            "SAMPA":SAMPA,
            "GAMPA":GAMPA,
            "SGABA":SGABA,
            "GGABA":GGABA,
            }
            
        for key in current_stats:
            Stats[key] = current_stats[key]
            
            
    
    return Stats
