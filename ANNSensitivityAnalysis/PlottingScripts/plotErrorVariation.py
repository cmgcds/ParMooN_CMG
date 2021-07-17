import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import dataFunctions as DF

from plotError import plotError

def getError(projectName, runNumber):
    curDir = os.getcwd();
    #_______________________________________________________
    # Set paths and variables
    #_______________________________________________________
    # Output path for sensitivity analysis (./output)
    outputPath = os.getcwd()+'/output';
    # Location for storing test case (name is specified in inputData.py)
    projectOutputDir = outputPath+'/'+projectName; 
    runDir = projectOutputDir+'/'+str(runNumber); 

    # Change into the run directory inside output/caseName
    os.chdir(runDir);

    #_______________________________________________________
    # Prepare Sample data
    #_______________________________________________________

    #DF.createOutputSpace(runDir);

    # Read metadata
    # Load input space for sensitivity analysis
    inputData = np.loadtxt("inputSpace.dat");
    # Total number of samples
    numberOfSamples = inputData.shape[0];

    # Load the output space for sensitivity analysis
    outputData = np.loadtxt("outputSpace.dat");

    #_______________________________________________________
    # Process the data 
    #_______________________________________________________
    p5 = np.percentile(outputData,5, axis=0);
    p5={'L1Error':p5[0], 'L2Error':p5[1], 'MinError':p5[2], 'MaxError':p5[3], 'MSError':p5[4]};

    p95 = np.percentile(outputData,95, axis=0);
    p95={'L1Error':p95[0], 'L2Error':p95[1], 'MinError':p95[2], 'MaxError':p95[3], 'MSError':p95[4]};

    pMean = np.mean(outputData, axis=0);
    pMean={'L1Error':pMean[0], 'L2Error':pMean[1], 'MinError':pMean[2], 'MaxError':pMean[3], 'MSError':pMean[4]};

    pStd = np.std(outputData, axis=0);
    pStd={'L1Error':pStd[0], 'L2Error':pStd[1], 'MinError':pStd[2], 'MaxError':pStd[3], 'MSError':pStd[4]};

    os.chdir(curDir);
    return (pMean, p5, p95, pStd);


def plotErrorVariation(projectName):

    FontSize = 9;

    curDir = os.getcwd();

    TotalRuns = 8;

    # 4 metrics: mean value, 5 percentile, 95 percentile, standard deviation
    L1Error = np.zeros(shape=(TotalRuns, 4));
    L2Error = np.zeros(shape=(TotalRuns, 4));
    MSError = np.zeros(shape=(TotalRuns, 4));
    MinError = np.zeros(shape=(TotalRuns, 4));
    MaxError = np.zeros(shape=(TotalRuns, 4));

    projectDir = os.getcwd()+"/output/"+projectName+"/";

    for runNumber in range(TotalRuns):
        print (runNumber);
        (pMean, p5, p95, pStd) = getError(projectName , runNumber);

        L1Error[runNumber, 0] = pMean['L1Error'];
        L2Error[runNumber, 0] = pMean['L2Error'];
        MinError[runNumber, 0] = pMean['MinError'];
        MaxError[runNumber, 0] = pMean['MaxError'];
        MSError[runNumber, 0] = pMean['MSError'];

        L1Error[runNumber,  1] = p5['L1Error'];
        L2Error[runNumber,  1] = p5['L2Error'];
        MinError[runNumber, 1] = p5['MinError'];
        MaxError[runNumber, 1] = p5['MaxError'];
        MSError[runNumber,  1] = p5['MSError'];

        L1Error[runNumber,  2] = p95['L1Error'];
        L2Error[runNumber,  2] = p95['L2Error'];
        MinError[runNumber, 2] = p95['MinError'];
        MaxError[runNumber, 2] = p95['MaxError'];
        MSError[runNumber,  2] = p95['MSError'];

        L1Error[runNumber,  3] = pStd['L1Error'];
        L2Error[runNumber,  3] = pStd['L2Error'];
        MinError[runNumber, 3] = pStd['MinError'];
        MaxError[runNumber, 3] = pStd['MaxError'];
        MSError[runNumber,  3] = pStd['MSError'];

        pass;

    #_______________________________________________________
    # plot the data: 
    #_______________________________________________________

    x = np.arange(TotalRuns);

    fig, axs = plt.subplots(2,2, figsize=(6.4,4.8), dpi=300, constrained_layout=True);

    # Title:
    plt.suptitle("Variation of error vs. size of training dataset", fontsize=10);

    # Legend location
    legendLocation = "upper right";

    #_______________________________________________________
    # L1 error plot
    #_______________________________________________________
    location = (0,0);
    axs[location].semilogy(x, L1Error[:,0], '-', color ='black', linewidth='1.5', label = 'Mean error');
    axs[location].semilogy(x, L1Error[:,1], '--', color ='green', linewidth='2', label = 'p5 error');
    axs[location].semilogy(x, L1Error[:,2], '-.', color ='red', linewidth='2', label = 'p95 error');
    axs[location].semilogy(x, L1Error[:,3], '.', color ='blue', linewidth='3', label = 'St. Dev.');

    axs[location].set_xlabel(r"TDS");
    axs[location].set_xticks([0,1,2,3,4,5,6,7]);
    axs[location].set_xticklabels(['6','12','25','50','100','200','400','800']);

    axs[location].set_ylabel(r"Aggr. $L_1$ Error");
    axs[location].set_ylim(10**-1, 10**3);
    axs[location].legend(loc=legendLocation, fontsize=FontSize);

    #_______________________________________________________
    # MS error plot
    #_______________________________________________________
    location = (1,0);
    axs[location].semilogy(x, MSError[:,0], '-', color ='black', linewidth='1.5', label = 'Mean error');
    axs[location].semilogy(x, MSError[:,1], '--', color ='green', linewidth='2', label = 'p5 error');
    axs[location].semilogy(x, MSError[:,2], '-.', color ='red', linewidth='2', label = 'p95 error');
    axs[location].semilogy(x, MSError[:,3], '.', color ='blue', linewidth='3', label = 'St. Dev.');

    axs[location].set_xlabel(r"TDS");
    axs[location].set_xticks([0,1,2,3,4,5,6,7]);
    axs[location].set_xticklabels(['6','12','25','50','100','200','400','800']);

    axs[location].set_ylabel(r"Aggr. MS Error");
    axs[location].set_ylim(10**-4, 10**4);
    axs[location].legend(loc=legendLocation, fontsize=FontSize);

    #_______________________________________________________
    # Rel. Min error plot
    #_______________________________________________________
    location = (0,1);
    axs[location].semilogy(x, MinError[:,0], '-', color ='black', linewidth='1.5', label = 'Mean error');
    axs[location].semilogy(x, MinError[:,1], '--', color ='green', linewidth='2', label = 'p5 error');
    axs[location].semilogy(x, MinError[:,2], '-.', color ='red', linewidth='2', label = 'p95 error');
    axs[location].semilogy(x, MinError[:,3], '.', color ='blue', linewidth='3', label = 'St. Dev.');

    axs[location].set_xlabel(r"TDS");
    axs[location].set_xticks([0,1,2,3,4,5,6,7]);
    axs[location].set_xticklabels(['6','12','25','50','100','200','400','800']);

    axs[location].set_ylabel(r"Aggr. Min Error");
    axs[location].set_ylim(10**-2, 10**2);
    axs[location].legend(loc=legendLocation, fontsize=FontSize);

    #_______________________________________________________
    # Max error plot
    #_______________________________________________________
    location = (1,1);
    axs[location].semilogy(x, MaxError[:,0], '-', color ='black', linewidth='1.5', label = 'Mean error');
    axs[location].semilogy(x, MaxError[:,1], '--', color ='green', linewidth='2', label = 'p5 error');
    axs[location].semilogy(x, MaxError[:,2], '-.', color ='red', linewidth='2', label = 'p95 error');
    axs[location].semilogy(x, MaxError[:,3], '.', color ='blue', linewidth='3', label = 'St. Dev.');

    axs[location].set_xlabel(r"TDS");
    axs[location].set_xticks([0,1,2,3,4,5,6,7]);
    axs[location].set_xticklabels(['6','12','25','50','100','200','400','800']);

    axs[location].set_ylabel(r"Aggr. Max Error");
    axs[location].set_ylim(10**4, 10**10);
    axs[location].legend(loc=legendLocation, fontsize=FontSize);

    #plt.ticklabel_format(axis="y", style="sci", scilimits=(0,00))
   
    os.chdir(projectDir);

    plt.savefig("ErrorVariation.pdf");

    os.chdir(curDir);

    pass;



if __name__ == "__main__":
    from cycler import cycler
    plt.rcParams["text.usetex"] = True
    plt.rcParams["font.family"] = "serif"
    plt.rcParams["font.serif"] = ["Computer Modern"]
    plt.rcParams['xtick.labelsize'] = 10
    plt.rcParams['ytick.labelsize'] = 10
    plt.rcParams['axes.labelsize'] = 10
    plt.rcParams['axes.prop_cycle'] = cycler(color=['darkblue', '#d62728', '#2ca02c', '#ff7f0e', '#bcbd22', '#8c564b', '#17becf', '#9467bd', '#e377c2', '#7f7f7f'])
    plotErrorVariation("Avg");
