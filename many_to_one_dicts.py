import numpy as np
import pickle
import os
import sys

folder = str(sys.argv[1])
directory = os.fsencode(folder) 

all_dicts = []

for file in os.listdir(directory):
        filename = os.fsdecode(file)
        if filename[-4::] ==".npy":
                try:
                	dict_array = np.load(folder+'/'+filename, allow_pickle = True)
                	all_dicts.extend(dict_array)
                except:
                	print("Unable to open \n", filename)
tosave = np.array(all_dicts)
np.save(folder+'.npy',tosave)