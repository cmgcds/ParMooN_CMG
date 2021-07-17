import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc


def plotValid(runNumber, sampleNumber, exptNumber):
    #_______________________________________________________
    # Set paths and variables
    #_______________________________________________________
    currDir = os.getcwd();
    # Output path for sensitivity analysis (./output)
    outputPath = os.getcwd()+'/output';

    #_______________________________________________________
    # Plot for Expt1
    projectName = "Expt"+str(exptNumber);
    projectOutputDir = outputPath+'/'+projectName; 
    runDir = projectOutputDir+'/'+str(runNumber); 
    sampleDir = runDir+'/'+str(sampleNumber);
    data1 = np.loadtxt(sampleDir+"/testResults.csv",delimiter=',');



    s1color = 'red';
    s2color = 'green';
    s3color = 'darkblue';
    s4color = 'black';

    s1size = '5';
    s2size = '5';
    s3size = '5';
    s4size = '2';

    corr1 = np.corrcoef(data1[:,0], data1[:,1])[0,1];
    label1 = "Expt 1 (Corr: %0.2f"%corr1 + ")";

    print(corr1);

    pass;


if __name__ == "__main__":
    from cycler import cycler

    plt.rcParams["text.usetex"] = True
    plt.rcParams["font.family"] = "serif"
    plt.rcParams["font.serif"] = ["Computer Modern"]

    plt.rcParams['xtick.labelsize'] = 20
    plt.rcParams['ytick.labelsize'] = 20
    plt.rcParams['axes.prop_cycle'] = cycler(color=['darkblue', '#d62728', '#2ca02c', '#ff7f0e', '#bcbd22', '#8c564b', '#17becf', '#9467bd', '#e377c2', '#7f7f7f'])


    ''' 
    # BEST and WORST networks for all dataset sizes:
    # accoring to [L1, MS, Min, Max] errors
    __________________________________
    # TDS: 6
    Best Network for the Avg dataset:
    [4732, 4110, 1585, 1409]
    Worst Network for the Avg dataset:
    [2605, 17, 6597, 17]
    __________________________________
    # TDS: 12
    Best Network for the Avg dataset:
    [1565, 1491, 1455, 2037]
    Worst Network for the Avg dataset:
    [2605, 2605, 4673, 2605]
    __________________________________
    # TDS: 25
    Best Network for the Avg dataset:
    [1553, 4272, 7476, 2325]
    Worst Network for the Avg dataset:
    [2605, 2605, 4289, 3021]
    __________________________________
    # TDS: 50
    Best Network for the Avg dataset:
    [4450, 1045, 2281, 2325]
    Worst Network for the Avg dataset:
    [5821, 5821, 4093, 5821]
    __________________________________
    # TDS: 100
    Best Network for the Avg dataset:
    [1438, 1053, 3626, 2037]
    Worst Network for the Avg dataset:
    [6597, 6597, 6597, 6597]
    __________________________________
    # TDS: 200
    Best Network for the Avg dataset:
    [1702, 1041, 3391, 2037]
    Worst Network for the Avg dataset:
    [2365, 6597, 4601, 6597]
    __________________________________
    # TDS: 400
    Best Network for the Avg dataset:
    [1438, 1438, 6847, 2037]
    Worst Network for the Avg dataset:
    [2365, 4301, 4290, 3749]
    __________________________________
    # TDS: 800
    Best Network for the Avg dataset:
    [1438, 1557, 6418, 2037]
    Worst Network for the Avg dataset:
    [1500, 2365, 4290, 3749]
    __________________________________

    '''
    # Name of the project
    runNumber = 6; # TDS 400

    # Best MSE for TDS 400
    BestMSEID = 1438;
    # Worst MSE for TDS 400
    WorstMSEID = 4301;

    sampleNumber = BestMSEID;
    sampleNumber = WorstMSEID;

    for exptNumber in range(50):
        plotValid(runNumber, sampleNumber, exptNumber);
        pass;
