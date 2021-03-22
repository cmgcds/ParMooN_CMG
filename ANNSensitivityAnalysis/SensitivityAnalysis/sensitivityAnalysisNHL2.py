# Sensitivity analysis of the NHL2 network
import os
from datetime import datetime
from decimal import Decimal

import numpy as np
import matplotlib.pyplot as plt

import openturns as ot
from openturns.viewer import View
import pandas as pd

# Name of the project
projectName = "ANN";

# Enter the run number i.e. folder output/ANN/runNUmber
runNumber = 0;

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

# Load input space for sensitivity analysis
inputData = np.loadtxt("inputSpace.dat");
inputData = inputData[72:720];
# Total number of samples
numberOfSamples = inputData.shape[0];
# Total number of input parameters i.e. dimension of Input space for sensitivity analysis
# Hidden layer not included in this study, hence -2
numberOfInputParam = 6;
# Total number of output results (4) 
# i.e. min error, max error, L1 error, L2 error
numberOfOutputParam = 4;

# Extract useful information about arrays etc. Refer Run.py to find out the sequence of data
sampleNumber = inputData[:,0]; # Will NOT be used
NHL = inputData[:,1]; # Will NOT be used
OPLTYPE = inputData[:,2];
HL_0_DIM = inputData[:,3];
HL_0_TYPE = inputData[:,4];
HL_1_DIM = inputData[:,5];
HL_1_TYPE = inputData[:,6];
HL_2_DIM = inputData[:,7]; # Will NOT be used
HL_2_TYPE = inputData[:,8]; # Will NOT be used
EPOCHS = inputData[:,9];


# Prepare output of model for samples
L1Error = np.zeros(numberOfSamples);
L2Error = np.zeros(numberOfSamples);
MinError = np.zeros(numberOfSamples);
MaxError = np.zeros(numberOfSamples);

for n in range(numberOfSamples):
    sampleDir = runDir+'/'+str(n);
    h_data = np.loadtxt(sampleDir+'/testResults.csv', delimiter=',')
    MaxError[n] = max(h_data[:,2]);
    MinError[n] = min(h_data[:,2]);
    L1Error[n] = np.mean(abs(h_data[:,0] - h_data[:,1]));
    L2Error[n] = np.sqrt(np.sum((h_data[:,0] - h_data[:,1])**2));
    pass;

# Save the output space:
outputData = np.vstack((L1Error, L2Error, MinError, MaxError)).T;

print("Saving outputSpace...");
np.savetxt("outputSpace_NHL2.dat", outputData);
print("Done.");

# make sure that we are in the run dir
os.chdir(runDir);

#_______________________________________________________
# Metamodel using OpenTurns
#_______________________________________________________

# Coordinates in input space for samples
coordinates = np.zeros(shape=(numberOfSamples, numberOfInputParam));
#coordinates[:,0] = NHL;
coordinates[:,0] = OPLTYPE;
coordinates[:,1] = HL_0_DIM;
coordinates[:,2] = HL_0_TYPE;
coordinates[:,3] = HL_1_DIM;
coordinates[:,4] = HL_1_TYPE;
#coordinates[:,5] = HL_2_DIM;
#coordinates[:,6] = HL_2_TYPE;
coordinates[:,5] = EPOCHS;

# Observations for these coordinates
observations = np.zeros(shape=(numberOfSamples, numberOfOutputParam));

observations[:,0] = L1Error;
observations[:,1] = L2Error;
observations[:,2] = MinError;
observations[:,3] = MaxError;


def getMetamodel(outputNumber):
    # prepare training data from coordinates and observations for Openturns
    input_train = ot.Sample(coordinates)
    output_train = ot.Sample([[ui] for ui in observations[:,outputNumber]])

    # Fit
    inputDimension = numberOfInputParam;
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
    return krigingMetamodel;




#input_names = ['NHL', 'OPLTYPE', 'HL_0_DIM', 'HL_0_TYPE', 'HL_1_DIM', 'HL_1_TYPE', 'HL_2_DIM', 'HL_2_TYPE','EPOCHS']
input_names = ['OPLTYPE', 'HL_0_DIM', 'HL_0_TYPE', 'HL_1_DIM', 'HL_1_TYPE','EPOCHS']

output_names = ['L1Error', 'L2Error','MinError','MaxError'];
#_______________________________________________________
# Find sobol indices
#_______________________________________________________

distributionList = [ot.Uniform(0.0,1.0)] * numberOfInputParam;
distribution = ot.ComposedDistribution(distributionList)

# Size of created data for metamodel
size = 500000

for i in range(numberOfOutputParam):
    # Get a fresh meta model for ith output
    print("Solving for ",output_names[i]," \n");

    print("1. Getting metamodel...");
    metamodel = getMetamodel(i);
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
    np.savez("sobolIndices_NHL2_"+output_names[i], sobol1=sobol1, soboltotal=soboltotal);
    print("Done.");


    # Draw the results
    '''
    graph = sensitivityAnalysis.draw()
    graph.setLegendFontSize(18);
    graph.setTitle(output_names[i]);
    ot.Show(graph);
    '''
    pass;

