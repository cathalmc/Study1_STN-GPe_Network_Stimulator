import numpy as np

class SimControl_lite():
    '''Class for generating the smaller string command to supply to the 
       STN-GPe simulator'''
       
    def __init__(self, orig=None):
        if orig is None:
            self.non_copy_constructor()
        else:
            self.copy_constructor(orig)
            
    def non_copy_constructor(self):
        '''Create command with default values'''
        self.params = {  
                "FileName":"TestNetwork.npy",
                "simtime":1500, 
                "Dummy":252332,
                "n":20,
                "h":0.02,
                "Notebook":0,
                "StimSites":0,
                "StimAmplitude":0.00001,
                "StimFrequency":1
                }
                
    def copy_constructor(self, orig):
        '''Create command by copying values of another SimControl object'''
        self.params={}
        for key,val in orig.params.items():
            self.params[key] = val
            

    def Generate_Command(self):
        sub_params=self.params
        '''Create the command string which will be ran by the os
        First map the dictionary values to an index
        This part seems redundent but it ensures that the dictionary items are always given in the correct order '''
        list_map = {  
                "FileName":0,
                "simtime":1, 
                "Dummy":2,
                "n":3,
                "h":4,
                "Notebook":5,
                "StimSites":6,
                "StimAmplitude":7,
                "StimFrequency":8
                }
        
        dict_list = [np.nan for i in range(len(list_map))]

        for key in list_map:
            dict_list[list_map[key]] = sub_params[key]
        strparams = ''.join(' {:}'.format(sub_params[x]) for x in sub_params)

        progname = "SimulationFromEL.py"
        ip = "python " + progname + strparams
        return ip



class SimControl():
    '''Class for generating the string command to supply to the 
       STN-GPe simulator'''
       
    def __init__(self, orig=None):
        if orig is None:
            self.non_copy_constructor()
        else:
            self.copy_constructor(orig)
            
    def non_copy_constructor(self):
        '''Create command with default values'''
        self.params = {  "k":10,
                "p": 1.0,
                "simtime": 5000, 
                "delay":2,
                "recip":1,
                "weight":0.02,
                "name":"testrun",
                "STNbias":-1.0,
                "Dummy":26278342,
                "GPebias":-0.2,
                "Network_type": "Small_world",
                
                "GSweight": 3.0,
                "GGweight": 0.02,

                "Strweight": 0.25,
                "CTXweight": 0.25,
                "convergence": 20,
                
                "n":100,
                "h":0.01,
                "Notebook":0,
                "StimSites":0,
                "StimAmplitude":15,
                "StimFrequency":140,
                
                }
                
    def copy_constructor(self, orig):
        '''Create command by copying values of another SimControl object'''
        self.params={}
        for key,val in orig.params.items():
            self.params[key] = val
            

    def Generate_Command(self):
        sub_params=self.params
        '''Create the command string which will be ran by the os
        First map the dictionary values to an index
        This part seems redundent but it ensures that the dictionary items are always given in the correct order '''
        list_map = {  "k":0,
                    "p": 1,
                    "simtime":2, 
                    "delay":3,
                    "recip":4,
                    "weight":5,
                    "name":6,
                    "STNbias":7,
                    "Dummy":8,
                    "GPebias":9,
                    "Network_type":10,
                    "GSweight": 11,
                    "GGweight": 12,
                    "Strweight": 13,
                    "CTXweight": 14,
                    "convergence": 15,
                    "n":16,
                    "h":17,
                    "Notebook":18,
                    "StimSites":19,
                "StimAmplitude":20,
                "StimFrequency":21,
                    }
        dict_list = [np.nan for i in range(len(list_map))]

        for key in list_map:
            dict_list[list_map[key]] = sub_params[key]
        strparams = ''.join(' {:}'.format(sub_params[x]) for x in sub_params)

        progname = "Simulation.py"
        ip = "python " + progname + strparams
        return ip


