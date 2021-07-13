### Script to setup all the execution for the automated scripts
import os
import numpy as np
import subprocess
import inspect

## define the parameters

# Number of Iterations
N_iterations = os.environ["N_iterations"];
N_iterations = int(N_iterations)
# Pe_nr values
epsnumbersList = [100000,1000000,10000000,100000000,1000000000,10000000000,100000000000,1000000000000]

## Model Values for Folder names
LISTFOLDERNAMES = ['BestTDS100','BestTDS200','BestTDS400','WorstTDS100','WorstTDS200','WorstTDS400']

# It has to be less than or equal to N_iterations
iterValue = 1

# Quadrature rule integration 
QUAD_RULE_NO = 3

## Define The Exact Solution fucntion here
def ExactSolution(x,eps):
    term1 = np.exp(- ((1-x)/eps ) )
    term2 = np.exp( -1.0/eps  )
    return  x - ( (term1 - term2)/ (1 - term2 ))



## Checking
if(iterValue >= N_iterations):
    iterValue = 0


def printExactSolution():
    lines = inspect.getsource(ExactSolution)
    print(lines)


# arrayDeclare = "declare -a epsnumbersList=("
# for i in epsnumbersList:
#     arrayDeclare = arrayDeclare + str(i)
#     if(epsnumbersList.index(i) != len(epsnumbersList) -1 ):
#         arrayDeclare = arrayDeclare + ","
# arrayDeclare = arrayDeclare + ")"



# ## Create a shell script 

# file1 = open('epsList.sh', 'w')
# file1.write(arrayDeclare)
# file1.close()

# os.system(". epsList.sh")