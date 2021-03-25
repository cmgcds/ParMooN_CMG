import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def getAveragePerformance(dirList, avgDirName):
    # dirList: list of directories containing experiment results e.g. ["Expt1", "Expt2", ...]
    # Note: dirList will contain the directories from ./output/
    # avgDirName: name of the directory to store the average results
    # Note: the avgDirName will be created inside ./output/.

    # Number of training datasets used i.e. size [6, 12 ,25, 50, 100, 200, 400 ,800] total 8
    trainingDataSize = 8;

    # NOTE: dirList is a list containing names of the runs e.g. ['Expt1', 'Expt2', ...]
    # NOTE: the relative or the abs address is NOT given

    # Get the location of the current directory
    curDir = os.getcwd();

    # This gives the total number of neural networks 
    inputDataSize = np.shape(np.loadtxt(curDir+'/output/'+dirList[0]+'/0/inputSpace.dat'))[0];

    # empty list
    outputDataList = []


    for dirName in dirList:
        data = np.zeros(shape=(trainingDataSize,inputDataSize,5));
        for i in range(trainingDataSize):
            data[i] = np.loadtxt(curDir+'/output/'+dirName+'/'+str(i)+'/outputSpace.dat');
            pass;
        outputDataList.append(data);
        pass;

    # outputDataList now contains the outputSpace from all the runs of all the directories 

    # Find the average of all the directories 
    averageData = np.mean(np.array(outputDataList), axis=0);

    # If average directory doesn't exist, create it
    avgDirPath = curDir + '/output/'+ avgDirName;

    if not os.path.exists(avgDirPath):
        os.makedirs(avgDirPath);
        for i in range(trainingDataSize):
            os.makedirs(avgDirPath+'/'+str(i));
            pass;
        pass;

    # Save the average output space for each training data size (i.e. in all 8 dirs)

    print("Saving average data...");
    for i in range(trainingDataSize):
        np.savetxt(avgDirPath+'/'+str(i)+'/outputSpace.dat', averageData[i]);
        pass;

    # Save the input space as well
    for i in range(trainingDataSize):
        os.system("cp "+curDir+'/output/'+dirList[0]+'/'+str(i)+'/inputSpace.dat '+ avgDirPath+'/'+str(i)+'/.');
        # Copy a metadata file from the first directory of the argument list. This is used only to know the size of the training data for plotting later
        os.system("cp "+curDir+'/output/'+dirList[0]+'/'+str(i)+'/metadata.dat '+ avgDirPath+'/'+str(i)+'/.');
        pass;


    print("Done");



def createOutputSpace(dirName):
    # NOTE: dirName is the name of the specific run full address e.g. ./output/Expt1/0 
    curDir = os.getcwd();

    # Change the directory to argument
    os.chdir(dirName);
    # Load input space for sensitivity analysis
    inputData = np.loadtxt("inputSpace.dat");
    # Total number of samples
    numberOfSamples = inputData.shape[0];
    # Prepare output of model for samples
    L1Error = np.zeros(numberOfSamples);
    L2Error = np.zeros(numberOfSamples);
    MinError = np.zeros(numberOfSamples);
    MaxError = np.zeros(numberOfSamples);
    MSError = np.zeros(numberOfSamples);

    # Create the output data arrays
    for n in range(numberOfSamples):
        sampleDir = dirName+'/'+str(n);
        h_data = np.loadtxt(sampleDir+'/testResults.csv', delimiter=',')
        # Relative maximum error
        MaxError[n] = max(h_data[:,2]/h_data[:,0]);
        # Relative minimum error
        MinError[n] = min(h_data[:,2]/h_data[:,0]);
        # Relative L1 norm
        L1Error[n] = np.sum(abs(h_data[:,0] - h_data[:,1]))/ np.sum(abs(h_data[:,0]));
        # Relative L2 norm
        L2Error[n] = np.sqrt(np.sum((h_data[:,0] - h_data[:,1])**2))/ np.sqrt(np.sum(h_data[:,0]**2));
        # Mean squared error
        MSError[n] = np.square(h_data[:,0] - h_data[:,1]).mean();
        pass;

    # Save the output space:
    outputData = np.vstack((L1Error, L2Error, MinError, MaxError, MSError)).T;

    #print("Saving outputSpace...");
    np.savetxt("outputSpace.dat", outputData);
    #print("Done.");

    # Change the directory back to previous
    os.chdir(curDir);
    pass;


def createTestingDataset(datasetName, size):
    # datasetName: name of the parent dataset without the .csv extension i.e. "allData" for allData.csv
    # size: size of the required testing dataset. i.e. number of testing points

    # create a dataset of given 'size' from 'datasetName' by random selection
    shuffleDataset(datasetName);

    df = pd.read_csv(datasetName+"1.csv", delimiter=',');

    dfNew = df.head(size);

    dfNew.to_csv("testingData.csv", index=False, header=False);

    os.system("rm "+datasetName+"1.csv");

    pass;


def createTrainingDataset(datasetName, size):
    # datasetName: name of the parent dataset without the .csv extension i.e. "allData" for allData.csv
    # size: size of the required testing dataset. i.e. number of testing points

    # create a dataset of given 'size' from 'datasetName' by random selection
    shuffleDataset(datasetName);

    df = pd.read_csv(datasetName+"1.csv", delimiter=',');

    dfNew = df.head(size);

    dfNew.to_csv("trainingData.csv", index=False, header=False);

    os.system("rm "+datasetName+"1.csv");

    pass;

def createTrainingAndValidationDataset(datasetName, trainingSize, validationPercentage):
    # Creates a dataset of required training examples given by trainingSize
    # Also accounts for the validation examples based on the validation percentage 
    # Randomly selected from 'datasetName'
    trainingPercentage = 100. - validationPercentage;

    N = int(trainingSize * 100. / trainingPercentage);

    shuffleDataset(datasetName);
    
    df = pd.read_csv(datasetName+"1.csv", delimiter=',');

    dfNew = df.head(N);

    dfNew.to_csv("trainingData.csv", index=False, header=False);

    os.system("rm "+datasetName+"1.csv");

def shuffleDataset(fileName):
    # Read file
    df = pd.read_csv(fileName+".csv", delimiter=',');

    # Shuffle file
    df = df.reindex(np.random.permutation(df.index));

    # Write the file without the additional counter index (i.e. 1, 2, 3 on the left size)
    # NOTE: the data is written to the same file twice so that the total count is > 1600
    # This is required for the training data of size 800 with 50% validation percentage
    # It won't have any impact on the training as anything beyond 800 points will be used 
    # for validation anyway. 
    df.to_csv(fileName+"1.csv", index = False);
    df.to_csv(fileName+"1.csv", mode='a',index = False, header=False);

    pass;


def createScript(data):
    # Note: data is a dictionary
    # create a new script for the simulation

    file1 = open("runScript.dat", "w");
    
    #______________________________________________
    file1.write("# Script for running the ANN code \n");
    #______________________________________________

    file1.write("\nANN_NHL: "+str(data["NHL"]));

    file1.write("\nANN_IPDATADIM: "+str(data["IPDATADIM"]));
    file1.write("\nANN_IPLDIM: "+str(data["IPLDIM"]));
    file1.write("\nANN_OPLDIM: "+str(data["OPLDIM"]));
    file1.write("\nANN_OPLTYPE: "+str(data["OPLTYPE"]));

    for i in range(data["NHL"]):
        file1.write("\nANN_HL_"+str(i)+"_DIM: " + str(data["HL_"+str(i)+"_DIM"]));
        file1.write("\nANN_HL_"+str(i)+"_TYPE: " + str(data["HL_"+str(i)+"_TYPE"]));
        pass;

    file1.write("\nANN_TRAINING_DATASET_NAME: "+str(data["TRAINING_DATASET_NAME"]));

    file1.write("\nANN_TESTING_DATASET_NAME: "+str(data["TESTING_DATASET_NAME"]));

    file1.write("\nANN_VALIDATION_DATA_PERCENTAGE: "+str(data["VALIDATION_DATA_PERCENTAGE"]));

    file1.write("\nANN_OPTIMIZER_CODE: "+ str(data["OPTIMIZER_CODE"]));

    file1.write("\nANN_OPTIMIZER_STEP_SIZE: "+ str(data["OPTIMIZER_STEP_SIZE"]));

    file1.write("\nANN_SGD_BATCH_SIZE: "+ str(data["SGD_BATCH_SIZE"]));

    file1.write("\nANN_MAX_ITERATIONS: "+ str(data["MAX_ITERATIONS"]));

    file1.write("\nANN_TOLERANCE: "+ str(data["TOLERANCE"]));

    file1.write("\nANN_DROPOUT_RATIO: " + str(data["DROPOUT_RATIO"]));

    file1.write("\nANN_EPOCHS: " + str(data["EPOCHS"]));

    file1.write("\nANN_SAVE_DATA_FILE: " + str(data["SAVE_DATA_FILE"]));

    file1.close();

    pass;
    


