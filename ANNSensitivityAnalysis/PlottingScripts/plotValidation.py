import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc


def plotValid(runNumber, sampleNumber):
    #_______________________________________________________
    # Set paths and variables
    #_______________________________________________________
    currDir = os.getcwd();
    # Output path for sensitivity analysis (./output)
    outputPath = os.getcwd()+'/output';

    #_______________________________________________________
    # Plot for Expt1
    projectName = "Expt1";
    projectOutputDir = outputPath+'/'+projectName; 
    runDir = projectOutputDir+'/'+str(runNumber); 
    sampleDir = runDir+'/'+str(sampleNumber);
    data1 = np.loadtxt(sampleDir+"/testResults.csv",delimiter=',');

    #_______________________________________________________
    # Plot for Expt2
    projectName = "Expt2";
    projectOutputDir = outputPath+'/'+projectName; 
    runDir = projectOutputDir+'/'+str(runNumber); 
    sampleDir = runDir+'/'+str(sampleNumber);
    data2 = np.loadtxt(sampleDir+"/testResults.csv",delimiter=',');

    #_______________________________________________________
    # Plot for Expt3
    projectName = "Expt3";
    projectOutputDir = outputPath+'/'+projectName; 
    runDir = projectOutputDir+'/'+str(runNumber); 
    sampleDir = runDir+'/'+str(sampleNumber);
    data3 = np.loadtxt(sampleDir+"/testResults.csv",delimiter=',');


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
    corr2 = np.corrcoef(data2[:,0], data2[:,1])[0,1];
    label2 = "Expt 2 (Corr: %0.2f"%corr2 + ")";
    corr3 = np.corrcoef(data3[:,0], data3[:,1])[0,1];
    label3 = "Expt 3 (Corr: %0.2f"%corr3 + ")";

    plt.figure(figsize=(6.4,4.8), dpi=300, constrained_layout = True);

    plt.loglog(data1[:,0], data1[:,1],"o",color=s1color,markersize=s1size,label= label1);
    plt.loglog(data2[:,0], data2[:,1],"s",color=s2color,markersize=s2size,label=label2);
    plt.loglog(data3[:,0], data3[:,1],"^",color=s3color,markersize=s3size,label=label3);

    lineData = np.linspace(10**-8, 1,1000);
    plt.plot(lineData, lineData,linewidth='1',color=s4color);

    plt.legend(loc=4, fontsize=16);

    plt.xlabel(r"Reference $\tau$", fontsize=16);
    plt.xlim(10**-8, 1);
    plt.ylim(10**-8, 1);

    plt.ylabel(r"ANN Predicted $\tau$", fontsize=16);
    #print(np.corrcoef(data[:,0], data[:,1]));
    plt.grid();


    plt.savefig("ValidBestMSE.pdf");
    #plt.savefig("ValidBestL1E.pdf");
    #plt.savefig("ValidWorstL1E.pdf");
    #plt.savefig("ValidWorstMSE.pdf");
    pass;


if __name__ == "__main__":
    from cycler import cycler

    plt.rcParams["text.usetex"] = True
    plt.rcParams["font.family"] = "serif"
    plt.rcParams["font.serif"] = ["Computer Modern"]

    plt.rcParams['xtick.labelsize'] = 16
    plt.rcParams['ytick.labelsize'] = 16
    plt.rcParams['axes.prop_cycle'] = cycler(color=['darkblue', '#d62728', '#2ca02c', '#ff7f0e', '#bcbd22', '#8c564b', '#17becf', '#9467bd', '#e377c2', '#7f7f7f'])

    # Name of the project
    runNumber = 5;
    BestMSEID = 1041;
    BestL1EID = 1702;
    WorstMSEID = 6597;
    WorstL1EID = 2365;
    sampleNumber = BestMSEID;
    plotValid(runNumber, sampleNumber);
