# NOTE: This file generates a frequency plot w.r.t. hyperparameters of the best 
# performing networks on three metrics (simultaneously)

import os
import time 

from datetime import datetime
from decimal import Decimal

import numpy as np

from matplotlib import pylab as plt
import pandas as pd


def plotFrequencies(projectName, runNumber):
    #_______________________________________________________
    # Set paths and variables
    #_______________________________________________________
    currDir = os.getcwd();
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

    # Load input space for sensitivity analysis
    inputData = np.loadtxt("inputSpace.dat");
    # Total number of samples
    numberOfSamples = inputData.shape[0];
    # Total number of input parameters i.e. dimension of Input space for sensitivity analysis
    numberOfInputParam = inputData.shape[1]-1;
    # Total number of output results (4) 
    # i.e. min error, max error, L1 error, L2 error
    numberOfOutputParam = 4;


    # Load the output space for sensitivity analysis
    outputData = np.loadtxt("outputSpace.dat");
    L1Error = outputData[:,0];
    L2Error = outputData[:,1];
    MinError = outputData[:,2];
    MaxError = outputData[:,3];
    MSError = outputData[:,4];


    #_______________________________________________________
    # Process the data 
    #_______________________________________________________
    p5 = np.percentile(outputData,5, axis=0);
    p5={'L1Error':p5[0], 'L2Error':p5[1], 'MinError':p5[2], 'MaxError':p5[3], 'MSError':p5[4]};

    p95 = np.percentile(outputData,95, axis=0);
    p95={'L1Error':p95[0], 'L2Error':p95[1], 'MinError':p95[2], 'MaxError':p95[3], 'MSError':p95[4]};
    #pMean = np.mean(outputData, axis=0);


    # Find out the NHL1, NHL2 and NHL3 data
    flag1 = 0;
    flag2 = 0;
    for i in range(numberOfSamples):
        if (inputData[i,1] == 1):
            flag1 = i+1;
        elif (inputData[i,1] == 2):
            flag2 = i+1;


    
    # Extract information for the NHL3
    inputData1 = inputData[0:];
    outputData1 = outputData[0:];

    print(p95['L1Error'], p95['MSError'], p95['MaxError']);
    list1 = inputData1[(outputData1[:,0])>p95['L1Error']];
    l1set = set([tuple(x) for x in list1])

    list2 = inputData1[(outputData1[:,4])>p95['MSError']];
    msset = set([tuple(x) for x in list2])

    list3 = inputData1[(outputData1[:,3])>p95['MaxError']];
    maxset = set([tuple(x) for x in list3])

    list4 = inputData1[(outputData1[:,2])>p95['MinError']];
    minset = set([tuple(x) for x in list4])


    # Select the networks with maximum L1, MS and Min error (i.e. min error bounded below)
    supersetIP = (np.array([x for x in l1set & msset & minset])).astype(int);
    print(np.shape(supersetIP));
    supersetOP = (outputData[supersetIP[:,0],:]);

    worstNetworks = (np.argmax(supersetOP,axis=0));
    print("L1Worst: ", supersetIP[worstNetworks[0]], "\nMSWorst: ", supersetIP[worstNetworks[4]], "\nMaxWorst: ", supersetIP[worstNetworks[3]]);

    print(np.shape(supersetIP));



    # Process the data
    NHL     = np.zeros(3);
    OPLTYPE = np.zeros(4);
    HL1TYPE = np.zeros(4);
    HL2TYPE = np.zeros(4);
    HL3TYPE = np.zeros(4);
    HL1DIM  = np.zeros(3);
    HL2DIM  = np.zeros(3);
    HL3DIM  = np.zeros(3);

    # NHL
    # corresponding to inputData[:,1] 
    NHL[0] = np.sum(supersetIP[:,1] == 1); 
    NHL[1] = np.sum(supersetIP[:,1] == 2);
    NHL[2] = np.sum(supersetIP[:,1] == 3);

    # Output layer type
    # corresponding to inputData[:,2] 
    OPLTYPE[0] = np.sum(supersetIP[:,2] == 0); 
    OPLTYPE[1] = np.sum(supersetIP[:,2] == 1);
    OPLTYPE[2] = np.sum(supersetIP[:,2] == 3);
    OPLTYPE[3] = np.sum(supersetIP[:,2] == 4);

    # HL1  layer type
    # corresponding to inputData[:,4] 
    HL1TYPE[0] = np.sum(supersetIP[:,4] == 0); 
    HL1TYPE[1] = np.sum(supersetIP[:,4] == 1);
    HL1TYPE[2] = np.sum(supersetIP[:,4] == 3);
    HL1TYPE[3] = np.sum(supersetIP[:,4] == 4);

    # HL2  layer type
    # corresponding to inputData[:,6] 
    HL2TYPE[0] = np.sum(supersetIP[:,6] == 0); 
    HL2TYPE[1] = np.sum(supersetIP[:,6] == 1);
    HL2TYPE[2] = np.sum(supersetIP[:,6] == 3);
    HL2TYPE[3] = np.sum(supersetIP[:,6] == 4);

    # HL3  layer type
    # corresponding to inputData[:,8] 
    HL3TYPE[0] = np.sum(supersetIP[:,8] == 0); 
    HL3TYPE[1] = np.sum(supersetIP[:,8] == 1);
    HL3TYPE[2] = np.sum(supersetIP[:,8] == 3);
    HL3TYPE[3] = np.sum(supersetIP[:,8] == 4);

    # HL1  layer dim
    # corresponding to inputData[:,3] 
    HL1DIM[0] = np.sum(supersetIP[:,3] == 5); 
    HL1DIM[1] = np.sum(supersetIP[:,3] == 10);
    HL1DIM[2] = np.sum(supersetIP[:,3] == 15);

    # HL2  layer dim
    # corresponding to inputData[:,5] 
    HL2DIM[0] = np.sum(supersetIP[:,5] == 5); 
    HL2DIM[1] = np.sum(supersetIP[:,5] == 10);
    HL2DIM[2] = np.sum(supersetIP[:,5] == 15);

    # HL3  layer dim
    # corresponding to inputData[:,7] 
    HL3DIM[0] = np.sum(supersetIP[:,7] == 5); 
    HL3DIM[1] = np.sum(supersetIP[:,7] == 10);
    HL3DIM[2] = np.sum(supersetIP[:,7] == 15);


    NHLX = np.array([1,2,3]);
    DIMX = np.array([5,10,15]);
    TYPEX = np.array([1,2,3,4]);
    YLABEL = r"No. of ANNs"
    #YTICKS = [0,10,20,30,40,50];
    YTICKS = [0,2,4,6];

    

    fig, axs = plt.subplots(3,3, figsize=(6.4,6.4), dpi=300, constrained_layout=True);
    Size = 14;
    Rot = 45;

    # NHL frequencies
    location = (0,0);
    axs[location].bar(NHLX, NHL, width=0.5);
    axs[location].set_yticks(YTICKS);
    axs[location].set_ylabel(YLABEL,size=Size);
    axs[location].set_xlabel(r"NHL",size=Size);

    # HL1-D frequencies
    location = (0,1);
    axs[location].bar(DIMX, HL1DIM, width=2.5);
    axs[location].set_yticks(YTICKS);
    axs[location].set_ylabel(YLABEL,size=Size);
    axs[location].set_xlabel(r"HL1-D" ,size=Size);

    # HL2-D frequencies
    location = (0,2);
    axs[location].bar(DIMX, HL2DIM, width=2.5);
    axs[location].set_yticks(YTICKS);
    axs[location].set_ylabel(YLABEL,size=Size);
    axs[location].set_xlabel(r"HL2-D" ,size=Size);

    # HL3-D frequencies
    location = (1,0);
    axs[location].bar(DIMX, HL3DIM, width=2.5);
    axs[location].set_yticks(YTICKS);
    axs[location].set_ylabel(YLABEL,size=Size);
    axs[location].set_xlabel(r"HL3-D" ,size=Size);

    #plt.savefig("BestNetworkDim.pdf");


    #fig, axs = plt.subplots(2,2, figsize=(6.4,4.8), dpi=300, constrained_layout=True);
    #Size = 14;
    #Rot = 60;
    # HL1-A frequencies
    location = (2,0);
    axs[location].bar(TYPEX, HL1TYPE, width=0.75);
    axs[location].set_yticks(YTICKS);
    axs[location].set_ylabel(YLABEL,size=Size);
    axs[location].set_xlabel(r"HL1-A" ,size=Size);
    axs[location].set_xticks([1,2,3,4]);
    axs[location].set_xticklabels([r"AF-1",r"AF-2",r"AF-3",r"AF-4"]);
    axs[location].tick_params(axis='x', rotation=Rot)

    # HL2-A frequencies
    location = (2,1);
    axs[location].bar(TYPEX, HL2TYPE, width=0.75);
    axs[location].set_yticks(YTICKS);
    axs[location].set_ylabel(YLABEL,size=Size);
    axs[location].set_xlabel(r"HL2-A" ,size=Size);
    axs[location].set_xticks([1,2,3,4]);
    axs[location].set_xticklabels([r"AF-1",r"AF-2",r"AF-3",r"AF-4"]);
    axs[location].tick_params(axis='x', rotation=Rot)

    # HL3-A frequencies
    location = (2,2);
    axs[location].bar(TYPEX, HL3TYPE, width=0.75);
    axs[location].set_yticks(YTICKS);
    axs[location].set_ylabel(YLABEL,size=Size);
    axs[location].set_xlabel(r"HL3-A" ,size=Size);
    axs[location].set_xticks([1,2,3,4]);
    axs[location].set_xticklabels([r"AF-1",r"AF-2",r"AF-3",r"AF-4"]);
    axs[location].tick_params(axis='x', rotation=Rot)

    # OP-A frequencies
    location = (1,2);
    axs[location].bar(TYPEX, OPLTYPE, width=0.75);
    axs[location].set_yticks(YTICKS);
    axs[location].set_ylabel(YLABEL,size=Size);
    axs[location].set_xlabel(r"OP-A"  ,size=Size);
    axs[location].set_xticks([1,2,3,4]);
    axs[location].set_xticklabels([r"AF-1",r"AF-2",r"AF-3",r"AF-4"]);
    axs[location].tick_params(axis='x', rotation=Rot)


    # Glossary
    location = (1,1);
    axs[location].set_axis_off();
    axs[location].text(0.25,0.8,r"\underline{Key}",fontsize=16);
    axs[location].text(-0.05,0.55,r"AF-1: Sigmoid",fontsize=14);
    axs[location].text(-0.05,0.35,r"AF-2: Identity",fontsize=14);
    axs[location].text(-0.05, 0.15,r"AF-3: L-ReLU",fontsize=14);
    axs[location].text(-0.05, -0.05,r"AF-4: TanH",fontsize=14);
    axs[location].set_yticks([]);
    axs[location].set_yticklabels([]);
    axs[location].set_xticks([]);
    axs[location].set_xticklabels([]);



    plt.savefig("WorstNetwork.pdf");

    os.chdir(currDir);
    plt.savefig("WorstNetwork.pdf");
    pass;


if __name__=="__main__":
    from cycler import cycler

    plt.rcParams["text.usetex"] = True
    plt.rcParams["font.family"] = "serif"
    plt.rcParams["font.serif"] = ["Computer Modern"]

    plt.rcParams['xtick.labelsize'] = 14
    plt.rcParams['ytick.labelsize'] = 14
    plt.rcParams['axes.prop_cycle'] = cycler(color=['darkblue', '#d62728', '#2ca02c', '#ff7f0e', '#bcbd22', '#8c564b', '#17becf', '#9467bd', '#e377c2', '#7f7f7f'])

    # Name of the project
    projectName = "Avg";
    runNumber = 4;
    plotFrequencies(projectName, runNumber);
