import os
import time 
import random

from datetime import datetime
from decimal import Decimal

import numpy as np

import openturns as ot
from openturns.viewer import View
import openturns.viewer as viewer
from matplotlib import pylab as plt
import pandas as pd


def getMetamodel(coordinates, observations, outputNumber):
    print("Creating Kriging metamodel...");
    # prepare training data from coordinates and observations for Openturns
    input_train = ot.Sample(coordinates)
    output_train = ot.Sample([[ui] for ui in observations[:,outputNumber]])

    # Fit
    inputDimension = np.shape(coordinates)[1];
    print("INPUT DIM: ",inputDimension);
    # Basis functions (linear)
    #basis = ot.LinearBasisFactory(inputDimension).build()
    # Basis functions (constant)
    basis = ot.LinearBasisFactory(inputDimension).build()
    # Statistical noise
    covarianceModel = ot.SquaredExponential([1.]*inputDimension, [1.0])
    covarianceModel.setActiveParameter([]);
    # Run Kriging algorithm
    algo = ot.KrigingAlgorithm(input_train, output_train, covarianceModel, basis, )
    algo.run()
    # extract metamodel
    result = algo.getResult()
    krigingMetamodel = result.getMetaModel()
    print("Done.");
    return result, krigingMetamodel;


def plotKrigingMetamodel(projectName, runNumber, size):
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

    # Load the output space for sensitivity analysis
    outputData = np.loadtxt("outputSpace.dat");

    # Find out the NHL1, NHL2 and NHL3 data
    flag1 = 0;
    flag2 = 0;
    for i in range(inputData.shape[0]):
        if (inputData[i,1] == 1):
            flag1 = i+1;
        elif (inputData[i,1] == 2):
            flag2 = i+1;

    # Select the NH1 data alone
    inputData = inputData[:flag1];
    outputData = outputData[:flag1];

    # Total number of samples
    numberOfSamples = inputData.shape[0];
    print("\nNumber of samples: " + str(numberOfSamples));
    # Total number of input parameters i.e. dimension of Input space for sensitivity analysis
    numberOfInputParam = 7;
    # Total number of output results (4) 
    # i.e. min error, max error, L1 error, L2 error
    numberOfOutputParam = 4;

    # Extract useful information about arrays etc. Refer Run.py to find out the sequence of data
    sampleNumber = inputData[:,0]; # not used
    NHL = inputData[:,1]; # not used
    OPLTYPE = inputData[:,2];
    HL_0_DIM = inputData[:,3];
    HL_0_TYPE = inputData[:,4];
    HL_1_DIM = inputData[:,5];  # Will NOT be used
    HL_1_TYPE = inputData[:,6];  # Will NOT be used
    HL_2_DIM = inputData[:,7];  # Will NOT be used
    HL_2_TYPE = inputData[:,8];  # Will NOT be used


    L1Error = outputData[:,0];
    L2Error = outputData[:,1];  # Not used, MSError used instead
    MinError = outputData[:,2];
    MaxError = outputData[:,3];
    MSError = outputData[:,4];

    #_______________________________________________________
    # Metamodel using OpenTurns
    #_______________________________________________________

    # Coordinates in input space for samples
    coordinates = np.zeros(shape=(numberOfSamples, numberOfInputParam));
    coordinates[:,0] = OPLTYPE;
    coordinates[:,1] = HL_0_DIM;
    coordinates[:,2] = HL_0_TYPE;
    coordinates[:,3] = HL_1_DIM;
    coordinates[:,4] = HL_1_TYPE;
    coordinates[:,5] = HL_2_DIM;
    coordinates[:,6] = HL_2_TYPE;

    # Observations for these coordinates
    observations = np.zeros(shape=(numberOfSamples, numberOfOutputParam));

    observations[:,0] = L1Error;
    observations[:,1] = MSError;
    observations[:,2] = MinError;
    observations[:,3] = MaxError;





    input_names = ['OPLTYPE', 'HL_0_DIM', 'HL_0_TYPE', 'HL_1_DIM', 'HL_1_TYPE', 'HL_2_DIM', 'HL_2_TYPE']

    output_names = ['L1Error', 'MSError','MinError','MaxError'];
    #_______________________________________________________
    # Find sobol indices
    #_______________________________________________________

    
    ot.RandomGenerator.SetSeed(int(1000*time.time()))
    distributionNHL = ot.Uniform(1,3);
    distributionOPLTYPE = ot.Uniform(0,4);
    distributionHL0DIM = ot.Uniform(5,15);
    distributionHL0TYPE = ot.Uniform(0,4);
    distributionHL1DIM = ot.Uniform(5,15);
    distributionHL1TYPE = ot.Uniform(0,4);
    distributionHL2DIM = ot.Uniform(5,15);
    distributionHL2TYPE = ot.Uniform(0,4);

    distributionList = [distributionOPLTYPE, distributionHL0DIM, distributionHL0TYPE, distributionHL1DIM, distributionHL1TYPE, distributionHL2DIM, distributionHL2TYPE];
    distribution = ot.ComposedDistribution(distributionList)


    for i in range(numberOfOutputParam):
        # Get a fresh meta model for ith output
        print("Solving for ",output_names[i]," \n");

        # Create test arrays for input and output
        OP = np.zeros(shape=(20,1));
        indices = random.sample(range(1, numberOfSamples), 20)
        OP[:,0] = np.array(observations[indices,i]);
        X_test = ot.Sample(np.array(coordinates[indices,:])); 
        Y_test = ot.Sample(OP); 

        coordinates_new = np.delete(coordinates, indices, 0);
        observations_new = np.delete(observations, indices, 0);
        

        print("1. Getting metamodel...");
        krigingResult, metamodel = getMetamodel(coordinates, observations,i);
        #krigingResult, metamodel = getMetamodel(coordinates_new, observations_new,i);
        #krigingResult, metamodel = getMetamodel(coordinates, observations,i);


        val = ot.MetaModelValidation(X_test, Y_test,metamodel)
        Q2 = val.computePredictivityFactor()[0]
        r = val.getResidualSample()

        graph = val.drawValidation()
        graph.setTitle("Q2 = %.2f%%" % (100*Q2))
        view = viewer.View(graph)

        plt.show()

        pass;
    os.chdir(currDir);
    pass;


if __name__=="__main__":
    #os.sched_setaffinity(0,{i for i in range(28)})

    # Name of the project
    projectName = "Avg";

    # Data size for running Sobol analysis. i.e. these many samples will be generated using the metamodel
    size = 100000;

    runNumber = 3;
    plotKrigingMetamodel(projectName, runNumber, size); 
