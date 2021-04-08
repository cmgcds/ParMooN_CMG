import os
from datetime import datetime
from decimal import Decimal

import numpy as np

import openturns as ot
from openturns.viewer import View
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
    return krigingMetamodel;


def getSobolIndices(projectName, runNumber, size):
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
    numberOfInputParam = 3;
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

    # Observations for these coordinates
    observations = np.zeros(shape=(numberOfSamples, numberOfOutputParam));

    observations[:,0] = L1Error;
    observations[:,1] = MSError;
    observations[:,2] = MinError;
    observations[:,3] = MaxError;





    input_names = ['OPLTYPE', 'HL_0_DIM', 'HL_0_TYPE']

    output_names = ['L1Error', 'MSError','MinError','MaxError'];
    #_______________________________________________________
    # Find sobol indices
    #_______________________________________________________

    distributionList = [ot.Uniform(0.0,1.0)] * numberOfInputParam;
    distribution = ot.ComposedDistribution(distributionList)


    for i in range(numberOfOutputParam):
        # Get a fresh meta model for ith output
        print("Solving for ",output_names[i]," \n");

        print("1. Getting metamodel...");
        metamodel = getMetamodel(coordinates, observations,i);
        # create new experiement
        sie = ot.SobolIndicesExperiment(distribution, size)
        print("2. Generating input and output data...");
        # generate fresh input data
        inputDesign = sie.generate()
        inputDesign.setDescription(input_names)
        # generate corresponding output data
        outputDesign = metamodel(inputDesign)

        # perform Sobol analysis using Saltelli algorithm
        print("3. Running Saltelli algorithm...");
        sensitivityAnalysis = ot.SaltelliSensitivityAlgorithm(inputDesign, outputDesign, size)

        # Get First order indices
        sobol1=sensitivityAnalysis.getFirstOrderIndices();
        # Get Total order indices
        soboltotal=sensitivityAnalysis.getTotalOrderIndices()

        print(sobol1); 
        print(soboltotal);
        print("Saving Sobol indices...");
        np.savez("sobolIndices_NH1_"+output_names[i], sobol1=sobol1, soboltotal=soboltotal);
        print("Done.");

        # Draw the results
        '''
        graph = sensitivityAnalysis.draw()
        graph.setLegendFontSize(18);
        graph.setTitle(output_names[i]);
        ot.Show(graph);
        '''
        pass;
    os.chdir(currDir);
    pass;


if __name__=="__main__":
    #os.sched_setaffinity(0,{i for i in range(28)})

    # Name of the project
    projectName = "Avg";

    # Data size for running Sobol analysis. i.e. these many samples will be generated using the metamodel
    size = 100000;

    # Save Sobol' indices
    for runNumber in range(8):
        print("\n\n Training set number: ", runNumber, "  \n\n");
        getSobolIndices(projectName, runNumber, size); 
