import os
from datetime import datetime
from decimal import Decimal
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from dataFunctions import *

import Run1HL as H1
import Run2HL as H2
import Run3HL as H3

def runExperiment(projectName, trainingSize):

    #______________________________________
    # Create paths and manage directories
    #______________________________________
    # Output directory
    outputPath = os.getcwd()+'/output';

    # Location for storing test case. Change 'ANN' to another name for different project
    projectOutputDir = outputPath+'/'+projectName; 

    # If the project of the given project-name doesn't exist, then create it
    if not os.path.exists(projectOutputDir):
        os.makedirs(projectOutputDir)

    # Location of the current run of the project (note: there could be multiple runs for the same project) created on different occasions
    # Each run in-turn will have multiple folders each corresponding to a sample point from input space
    # If no runs exist so far, a folder '1' will be automatically created for a run
    # Else a new folder will be created having number as 'last maximum number + 1'
    listDir = os.listdir(projectOutputDir);
    if (listDir == []):
        thisRun = '0';
    else:
        listDir.sort(key=lambda f: int(''.join(filter(str.isdigit, f))))
        lastRun = int(listDir[-1]);
        thisRun = str(lastRun + 1);

    thisRunDir = projectOutputDir+'/'+thisRun;
    # Create a folder for this run
    os.mkdir(thisRunDir);

    #_______________________________________________________
    # Create Metadata file
    #_______________________________________________________
    startDateTime = datetime.now() 
    date = startDateTime.strftime("%d-%b-%Y")
    time = startDateTime.strftime("%H:%M:%S")

    #_______________________________________________________
    # Start writing metadata in a file inside the current project run
    #_______________________________________________________
    metadataFile = open(thisRunDir+'/metadata.dat', 'a');
    metadataFile.write('Case name = '+ projectName);
    metadataFile.write('\nTest run number = '+ thisRun);
    metadataFile.write('\nAddress = '+ thisRunDir);
    metadataFile.write('\nDate = '+ date);
    metadataFile.write('\nStarting time (h:m:s) = '+ time);
    metadataFile.close();


    #_______________________________________________________
    # Start writing input space in a file inside the current project run
    #_______________________________________________________
    inputSpaceFile = open(thisRunDir+'/inputSpace.dat', 'a');
    inputSpaceFile.write('#sampleNo. NHL OPLTYPE HL_0_DIM HL_0_TYPE HL_1_DIM HL_1_TYPE HL_2_DIM HL_2_TYPE EPOCHS\n');

    inputSpaceFile.close();

    # Copy the training and the testing dataset
    os.system('cp trainingData.csv '+thisRunDir+'/.');
    os.system('cp testingData.csv '+thisRunDir+'/.');
    # Write the input and test dataset information into metadata file
    metadataFile = open(thisRunDir+'/metadata.dat', 'a');
    metadataFile.write("\nTraining data size :"+str(trainingSize));
    metadataFile.close();

    # Run simulations
    H1.runSimulations(thisRunDir);
    H2.runSimulations(thisRunDir);
    H3.runSimulations(thisRunDir);


    #_______________________________________________________
    # Create output space file from the simulation
    #_______________________________________________________
    createOutputSpace(thisRunDir);

    #_______________________________________________________
    # Finish writing metadata 
    #_______________________________________________________
    metadataFile = open(thisRunDir+'/metadata.dat', 'a');
    endDateTime = datetime.now() 
    time = endDateTime.strftime("%H:%M:%S")
    metadataFile.write('\nFinal time (h:m:s) = ' + time); 
    metadataFile.close();
    metadataFile.close();


if __name__ == "__main__":
    createTestingDataset("allData", 10);
    createTrainingAndValidationDataset("allData",10,10);
    runExperiment("ANN",10);
