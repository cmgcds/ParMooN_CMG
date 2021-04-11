import numpy as np
import matplotlib.pyplot as plt
from pylab import figure, cm
from matplotlib.colors import LogNorm
import os

def getIndices(projectName, runNumber, NHL):
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

    if (NHL == 0):
        dataL1 = np.load("sobolIndices_L1Error.npz");
        dataMS = np.load("sobolIndices_MSError.npz");
        dataMin = np.load("sobolIndices_MinError.npz");
        dataMax = np.load("sobolIndices_MaxError.npz");
    elif (NHL == 1):
        dataL1 = np.load("sobolIndices_NH1_L1Error.npz");
        dataMS = np.load("sobolIndices_NH1_MSError.npz");
        dataMin =np.load("sobolIndices_NH1_MinError.npz");
        dataMax =np.load("sobolIndices_NH1_MaxError.npz");
    elif (NHL == 2):
        dataL1 = np.load("sobolIndices_NH2_L1Error.npz");
        dataMS = np.load("sobolIndices_NH2_MSError.npz");
        dataMin =np.load("sobolIndices_NH2_MinError.npz");
        dataMax =np.load("sobolIndices_NH2_MaxError.npz");
    elif (NHL == 3):
        dataL1 = np.load("sobolIndices_NH3_L1Error.npz");
        dataMS = np.load("sobolIndices_NH3_MSError.npz");
        dataMin =np.load("sobolIndices_NH3_MinError.npz");
        dataMax =np.load("sobolIndices_NH3_MaxError.npz");

        pass;

    os.chdir(curDir);
    return (dataL1, dataMS, dataMin, dataMax);








def getPlotParameters(NHL):

    train = ['6','12','25','50','100','200','400','800'];

    # X axis for hyperparameters
    x = []
    if (NHL == 0):
        x = ['NHL', 'OP-A', 'HL1-D', 'HL1-A', 'HL2-D', 'HL2-A', 'HL3-D', 'HL3-A'];
        title = "Sobol' indices";
    elif (NHL == 1):
        x = ['OP-A', 'HL1-D', 'HL1-A'];
        title = "Sobol' indices for ANNs with one hidden layer";
    elif (NHL == 2):
        x = ['OP-A', 'HL1-D', 'HL1-A', 'HL2-D', 'HL2-A'];
        title = "Sobol' indices for ANNs with two hidden layers";
    elif (NHL == 3):
        x = ['OP-A', 'HL1-D', 'HL1-A', 'HL2-D', 'HL2-A', 'HL3-D', 'HL3-A'];
        title = "Sobol' indices for ANNs with three hidden layers";
        pass;

    linestyles = ['-o','-^','-s','-P','-X','-D','-*','-p'];

    markerstyles = ['o','^','s','P','X','D','*','p'];

    return x, train, title, linestyles, markerstyles;









def plotIndices(projectName, NHL):
    curDir = os.getcwd();

    projectDir = os.getcwd()+"/output/"+projectName+"/";

    TotalRuns = 8;

    indicesCount = [8,3,5,7];

    # List of the First and the Total Sobol indices
    indicesListFirst = np.zeros(shape=(TotalRuns,4,indicesCount[NHL]));
    indicesListTotal = np.zeros(shape=(TotalRuns,4,indicesCount[NHL]));

    for runNumber in range(TotalRuns):
        print(runNumber);
        (dataL1, dataMS, dataMin, dataMax) = getIndices(projectName, runNumber,NHL);
        indicesListFirst[runNumber,0] = dataL1['sobol1'];
        indicesListFirst[runNumber,1] = dataMS['sobol1'];
        indicesListFirst[runNumber,2] = dataMin['sobol1'];
        indicesListFirst[runNumber,3] = dataMax['sobol1'];

        indicesListTotal[runNumber,0] = dataL1['soboltotal'];
        indicesListTotal[runNumber,1] = dataMS['soboltotal'];
        indicesListTotal[runNumber,2] = dataMin['soboltotal'];
        indicesListTotal[runNumber,3] = dataMax['soboltotal'];

        pass;


    #__________________________________________
    # Plot Sobol indices with x axis as hyperP
    #__________________________________________
    x, train, title, linestyles, markerstyles = getPlotParameters(NHL);

    FigSize = (6.4,4.8); 
    #__________________________________________
    # Plot FO Sobol indices
    #__________________________________________

    fig, axs = plt.subplots(2,2, figsize=FigSize, dpi=300);

    # Title:
    #plt.suptitle(title);

    # Legend location
    legendLocation = "upper center";
    
    for i in range(np.size(train)):
        #__________________________________________
        # L1 error plot
        #__________________________________________
        location = (0,0);

        axs[location].plot(x,indicesListFirst[i,0], markerstyles[i], label= 'TDS: '+train[i]);
        #axs[location].legend(loc=legendLocation, fontsize=8, ncol=2);
        axs[location].set_ylim(-0.1,1.1);
        axs[location].set_ylabel(r"FO Indices");
        axs[location].tick_params(axis='x', rotation=60)
        axs[location].set_title(r"$L_1$ Error");

        #__________________________________________
        # MS error plot
        #__________________________________________
        location = (0,1);

        axs[location].plot(x,indicesListFirst[i,1], markerstyles[i], label= 'TDS: '+train[i]);
        #axs[location].legend(loc=legendLocation, fontsize=8, ncol=2);
        axs[location].set_ylim(-0.1,1.1);
        axs[location].set_ylabel(r"FO Indices");
        axs[location].tick_params(axis='x', rotation=60)
        axs[location].set_title("MS Error");

        #__________________________________________
        # Min error plot
        #__________________________________________
        location = (1,0);

        axs[location].plot(x,indicesListFirst[i,2], markerstyles[i], label= 'TDS: '+train[i]);
        #axs[location].legend(loc=legendLocation, fontsize=8, ncol=2);
        axs[location].set_ylim(-0.1,1.1);
        axs[location].set_ylabel(r"FO Indices");
        axs[location].tick_params(axis='x', rotation=60)
        axs[location].set_title("Min Error");

        #__________________________________________
        # Max error plot
        #__________________________________________
        location = (1,1);

        axs[location].plot(x,indicesListFirst[i,3], markerstyles[i], label= 'TDS: '+train[i]);
        #axs[location].legend(loc=legendLocation, fontsize=8, ncol=2);
        axs[location].set_ylim(-0.1,1.1);
        axs[location].set_ylabel(r"FO Indices");
        axs[location].set_title("Max Error");
        axs[location].tick_params(axis='x', rotation=60)
        handles, labels = axs[location].get_legend_handles_labels()
        pass;

    fig.tight_layout(rect=[0,0,1,0.9])
    fig.legend(handles, labels, loc='upper center',ncol=4)


    os.chdir(projectDir);
    plt.savefig("FOSobolIndices"+str(NHL)+".pdf");
    os.chdir(curDir);
    
    #__________________________________________
    # Plot TO Sobol indices
    #__________________________________________

    fig, axs = plt.subplots(2,2, figsize=FigSize, dpi=300);


    # Title:
    #plt.suptitle(title);

    # Legend location
    legendLocation = "upper center";
    
    for i in range(np.size(train)):
        #__________________________________________
        # L1 error plot
        #__________________________________________
        location = (0,0);

        axs[location].plot(x,indicesListTotal[i,0], markerstyles[i], label= 'TDS: '+train[i]);
        axs[location].set_ylim(-0.1,1.1);
        axs[location].set_ylabel(r"TO Indices");
        axs[location].tick_params(axis='x', rotation=60)
        axs[location].set_title(r"$L_1$ Error");

        #__________________________________________
        # MS error plot
        #__________________________________________
        location = (0,1);

        axs[location].plot(x,indicesListTotal[i,1], markerstyles[i], label= 'TDS: '+train[i]);
        axs[location].set_ylim(-0.1,1.1);
        axs[location].set_ylabel(r"TO Indices");
        axs[location].tick_params(axis='x', rotation=60)
        axs[location].set_title("MS Error");

        #__________________________________________
        # Min error plot
        #__________________________________________
        location = (1,0);

        axs[location].plot(x,indicesListTotal[i,2], markerstyles[i], label= 'TDS: '+train[i]);
        axs[location].set_ylim(-0.1,1.1);
        axs[location].set_ylabel(r"TO Indices");
        axs[location].tick_params(axis='x', rotation=60)
        axs[location].set_title("Min Error");

        #__________________________________________
        # Max error plot
        #__________________________________________
        location = (1,1);

        axs[location].plot(x,indicesListTotal[i,3], markerstyles[i], label= 'TDS: '+train[i]);
        axs[location].set_ylim(-0.1,1.1);
        axs[location].set_ylabel(r"TO Indices");
        axs[location].set_title("Max Error");
        axs[location].tick_params(axis='x', rotation=60)
        handles, labels = axs[location].get_legend_handles_labels()
        pass;

    fig.tight_layout(rect=[0,0,1,0.9])
    fig.legend(handles, labels, loc='upper center',ncol=4)


    os.chdir(projectDir);
    plt.savefig("TOSobolIndices"+str(NHL)+".pdf");
    os.chdir(curDir);



    #__________________________________________
    # Plot FO Sobol indices variation
    #__________________________________________

    fig, axs = plt.subplots(2,2, figsize=FigSize, dpi=300);


    # Legend location
    legendLocation = "upper center";

    for i in range(np.size(x)):
        #__________________________________________
        # L1 error plot
        #__________________________________________
        location = (0,0);
        axs[location].set_ylim(-0.1,1.1);
        axs[location].set_ylabel(r"FO Indices");
        axs[location].set_xlabel("Training dataset size");
        axs[location].plot(train,indicesListFirst[:,0,i], linestyles[i], label= ': '+x[i]);
        axs[location].tick_params(axis='x', rotation=60)
        axs[location].set_title(r"$L_1$ Error");
        handles, labels = axs[location].get_legend_handles_labels()

        #__________________________________________
        # MS error plot
        #__________________________________________
        location = (0,1);
        axs[location].set_ylim(-0.1,1.1);
        axs[location].set_ylabel(r"FO Indices");
        axs[location].set_xlabel("Training dataset size");
        axs[location].plot(train,indicesListFirst[:,1,i], linestyles[i], label= ': '+x[i]);
        axs[location].tick_params(axis='x', rotation=60)
        axs[location].set_title("MS Error");

        #__________________________________________
        # Min error plot
        #__________________________________________
        location = (1,0);
        axs[location].set_ylim(-0.1,1.1);
        axs[location].set_ylabel(r"FO Indices");
        axs[location].set_xlabel("Training dataset size");
        axs[location].plot(train,indicesListFirst[:,2,i], linestyles[i], label= ': '+x[i]);
        axs[location].tick_params(axis='x', rotation=60)
        axs[location].set_title("Min Error");

        #__________________________________________
        # Max error plot
        #__________________________________________
        location = (1,1);
        axs[location].set_ylim(-0.1,1.1);
        axs[location].set_ylabel(r"FO Indices");
        axs[location].set_xlabel("Training dataset size");
        axs[location].plot(train,indicesListFirst[:,3,i], linestyles[i], label= ': '+x[i]);
        axs[location].tick_params(axis='x', rotation=60)
        axs[location].set_title("Max Error");


        pass;

    fig.tight_layout(rect=[0,0,1,0.9])
    fig.legend(handles, labels, loc='upper center',ncol=4)

    os.chdir(projectDir);
    plt.savefig("FOSobolIndices"+str(NHL)+"Variation.pdf");
    os.chdir(curDir);

    

    #__________________________________________
    # Plot TO Sobol indices variation
    #__________________________________________

    fig, axs = plt.subplots(2,2, figsize=FigSize, dpi=300);


    # Legend location
    legendLocation = "upper center";

    for i in range(np.size(x)):
        #__________________________________________
        # L1 error plot
        #__________________________________________
        location = (0,0);
        axs[location].set_ylim(-0.1,1.1);
        axs[location].set_ylabel(r"TO Indices");
        axs[location].set_xlabel("Training dataset size");
        axs[location].plot(train,indicesListTotal[:,0,i], linestyles[i], label= ': '+x[i]);
        axs[location].tick_params(axis='x', rotation=60)
        axs[location].set_title(r"$L_1$ Error");
        handles, labels = axs[location].get_legend_handles_labels()

        #__________________________________________
        # MS error plot
        #__________________________________________
        location = (0,1);
        axs[location].set_ylim(-0.1,1.1);
        axs[location].set_ylabel(r"TO Indices");
        axs[location].set_xlabel("Training dataset size");
        axs[location].plot(train,indicesListTotal[:,1,i], linestyles[i], label= ': '+x[i]);
        axs[location].tick_params(axis='x', rotation=60)
        axs[location].set_title("MS Error");

        #__________________________________________
        # Min error plot
        #__________________________________________
        location = (1,0);
        axs[location].set_ylim(-0.1,1.1);
        axs[location].set_ylabel(r"TO Indices");
        axs[location].set_xlabel("Training dataset size");
        axs[location].plot(train,indicesListTotal[:,2,i], linestyles[i], label= ': '+x[i]);
        axs[location].tick_params(axis='x', rotation=60)
        axs[location].set_title("Min Error");

        #__________________________________________
        # Max error plot
        #__________________________________________
        location = (1,1);
        axs[location].set_ylim(-0.1,1.1);
        axs[location].set_ylabel(r"TO Indices");
        axs[location].set_xlabel("Training dataset size");
        axs[location].plot(train,indicesListTotal[:,3,i], linestyles[i], label= ': '+x[i]);
        axs[location].tick_params(axis='x', rotation=60)
        axs[location].set_title("Max Error");


        pass;

    fig.tight_layout(rect=[0,0,1,0.9])
    fig.legend(handles, labels, loc='upper center',ncol=4)

    os.chdir(projectDir);
    plt.savefig("TOSobolIndices"+str(NHL)+"Variation.pdf");
    os.chdir(curDir);


    if (NHL == 0):
        VMIN = 1e-3;
        VMAX = 1e0
        #__________________________________________
        # Plot SO Sobol indices
        #__________________________________________
        # Get indices for training size of 200 (i.e. run number 5, with 0 NHL i.e. all 7536 expts)
        i = 5;
        (dataL1, dataMS, dataMin, dataMax) = getIndices(projectName, i,0);
        indicesListSecond = np.zeros(shape=(4,8,8));
        indicesListSecond[0] = dataL1['sobol2'];
        indicesListSecond[1] = dataMS['sobol2'];
        indicesListSecond[2] = dataMin['sobol2'];
        indicesListSecond[3] = dataMax['sobol2'];


        for i in range(4):
            # Set the autocorrelation to 1
            np.fill_diagonal(indicesListSecond[i], 1.0);
            pass;



        FigSize = (6.4,4.8);
        fig, axs = plt.subplots(2,2, figsize=FigSize, dpi=300);
        locations = [(0,0),(0,1),(1,0),(1,1)];
        # Legend location
        legendLocation = "upper center";
        
        #__________________________________________
        # L1 error plot
        #__________________________________________
        location = locations[0];
        im = axs[location].matshow(indicesListSecond[0],cmap=plt.cm.get_cmap('Blues', 6), norm=LogNorm(vmin=VMIN, vmax=VMAX));
        axs[location].set_title(r"$L_1$ Error");
        axs[location].set_yticks(np.arange(8));
        axs[location].set_yticklabels(x);
        axs[location].tick_params(axis='x', bottom=False, top=False, labeltop = False)
        axs[location].grid();

        #__________________________________________
        # MS error plot
        #__________________________________________
        location = locations[1];
        axs[location].matshow(indicesListSecond[1],cmap=plt.cm.get_cmap('Blues'), norm=LogNorm(vmin=VMIN, vmax=VMAX));
        axs[location].set_title(r"MS Error");
        axs[location].tick_params(axis='x', bottom=False, top=False, labeltop = False, labelbottom = False)
        axs[location].tick_params(axis='y', left=False, right=False, labelleft = False)
        axs[location].grid();

        #__________________________________________
        # Min error plot
        #__________________________________________
        location = locations[2];
        axs[location].matshow(indicesListSecond[2],cmap=plt.cm.get_cmap('Blues'), norm=LogNorm(vmin=VMIN, vmax=VMAX));
        axs[location].set_title(r"Min Error");
        axs[location].set_xticks(np.arange(8));
        axs[location].set_xticklabels(x);
        axs[location].xaxis.tick_bottom();
        axs[location].set_yticks(np.arange(8));
        axs[location].set_yticklabels(x);
        axs[location].tick_params(axis='x', rotation=75)
        axs[location].grid();

        #__________________________________________
        # Max error plot
        #__________________________________________
        location = locations[3];
        axs[location].matshow(indicesListSecond[3],cmap=plt.cm.get_cmap('Blues'), norm=LogNorm(vmin=VMIN, vmax=VMAX));
        axs[location].set_title(r"Max Error");
        axs[location].set_xticks(np.arange(8));
        axs[location].set_xticklabels(x);
        axs[location].xaxis.tick_bottom();
        axs[location].tick_params(axis='x', rotation=75)
        axs[location].tick_params(axis='y', left=False, right=False, labelleft = False)
        axs[location].grid();




        fig.colorbar(im, ax=axs[0, :2], shrink=0.5, location='top');
        fig.tight_layout(rect=[0.2,0,0.8,0.8])

        os.chdir(projectDir);
        plt.savefig("SOSobolIndices.pdf");
        os.system("pdfcrop SOSobolIndices.pdf SOSobolIndices.pdf");
        os.chdir(curDir);

        

    



if __name__== "__main__":
    from cycler import cycler
    plt.rcParams["font.family"] = "monospace"
    plt.rcParams['axes.grid'] = True
    plt.rcParams['axes.prop_cycle'] = cycler(color=['darkblue', '#d62728', '#2ca02c', '#ff7f0e', '#bcbd22', '#8c564b', '#17becf', '#9467bd', '#e377c2', '#7f7f7f'])
    for i in range(4):
        plotIndices("Avg",i);
