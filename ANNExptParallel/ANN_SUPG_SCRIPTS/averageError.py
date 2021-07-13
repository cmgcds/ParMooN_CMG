import numpy as np
from numpy import genfromtxt
import matplotlib.pyplot as plt
import sys

## import epsnumber list
from SetParameters import *




# %%

if __name__ == "__main__":
    
    # N_iterations = 2
    
    for name in LISTFOLDERNAMES:
        
        filename = name + "_" + str(0) + ".list"
        temp  = genfromtxt(filename,delimiter=',')
        e = np.zeros(temp.shape)
        normArray = []  ## to store the norm of Data
        for i in range(N_iterations):
            extension = i
            filename = name + "_" + str(i) + ".list"
            array = genfromtxt(filename,delimiter=',')
            normArray.append(np.linalg.norm(array))   ## Compute norm of the Array
            
        ## Compute z Scores 
        mean = np.mean(normArray)
        std_dev = np.std(normArray)
        
        divisor = N_iterations
        for i in range(N_iterations):
            extension = i
            filename = name + "_" + str(i) + ".list"
            ## Get the Z index score
            absZScore = abs((normArray[i] - mean)/std_dev)
            if( (absZScore < 1.8) and normArray[i] < 200 ):
                array = genfromtxt(filename,delimiter=',')
                e += array
            else:
                print (" File name : ", filename)
                divisor -= 1

            
        
        e = e/divisor
        ListName      = name
        listExtension = ".list" 
        savename = ListName + listExtension
        np.savetxt(savename, 
            e,
            delimiter =",", 
            fmt ='%10.9f')

# %%
# %%
# import numpy as np
# from numpy import genfromtxt
# LISTFOLDERNAMES = ['BestTDS100','BestTDS200','BestTDS400','WorstTDS100','WorstTDS200','WorstTDS400']
# for name in LISTFOLDERNAMES:
#     print (" ------------------------------ LIST FOLDER NAMES ------------------------------- " )
#     print (" ####### FOLDER " , name , "   ###############")
#     normArray = []
#     for i in range(20):
#         fileName = name + "_" + str(i) + ".list"
#         arr = np.genfromtxt(fileName,delimiter=',')
#         normArray.append(np.linalg.norm(arr))
        
    
#     ## Compute z Scores 
#     mean = np.mean(normArray)
#     std_dev = np.std(normArray)
    
#     for i in range(20):
#         print (fileName , " : " , normArray[i] , " Z : " , (normArray[i] - mean)/std_dev )



