import os
import time 

from datetime import datetime
from decimal import Decimal

import numpy as np

import openturns as ot
from openturns.viewer import View
import openturns.viewer as viewer
from matplotlib import pylab as plt
import pandas as pd

def linearSample(xmin,xmax,npoints):
    '''Returns a sample created from a regular grid
    from xmin to xmax with npoints points.'''
    step = (xmax-xmin)/(npoints-1)
    rg = ot.RegularGrid(xmin, step, npoints)
    vertices = rg.getVertices()
    return vertices

def plot_kriging_bounds(vLow,vUp,n_test):
    '''
    From two lists containing the lower and upper bounds of the region,
    create a PolygonArray.
    '''
    palette = ot.Drawable.BuildDefaultPalette(10)
    myPaletteColor = palette[6]
    polyData = [[vLow[i], vLow[i+1], vUp[i+1], vUp[i]] for i in range(n_test-1)]
    polygonList = [ot.Polygon(polyData[i], myPaletteColor, myPaletteColor) for i in range(n_test-1)]
    boundsPoly = ot.PolygonArray(polygonList)
    boundsPoly.setLegend("95% bounds")
    return boundsPoly


def plotBasicKriging(krigResult, meta, xMin, xMax, X, Y, level = 0.95):
    '''
    Given a kriging result, plot the data, the kriging metamodel
    and a confidence interval.
    '''
    samplesize = np.size(X)
    graphKriging = meta.draw(xMin, xMax)
    graphKriging.setLegends(["Kriging"])
    # Create a grid of points and evaluate the function and the kriging
    nbpoints = 100
    Index = 2;
    xGrid = np.ones(shape=(nbpoints,3));
    xGrid[:,Index] =  np.linspace(xMin, xMax, nbpoints);
    xGrid[:,0] =  0; 
    xGrid[:,1] =  15;  

    xGrid = ot.Sample(xGrid);
    #yFunction = g(xGrid)
    meta = krigResult.getMetaModel();
    yKrig = meta(xGrid)
    # Compute the conditional covariance
    epsilon = ot.Point(nbpoints,1.e-8)
    conditionalVariance = krigResult.getConditionalMarginalVariance(xGrid)+epsilon
    conditionalVarianceSample = ot.Sample([[cv] for cv in conditionalVariance])
    conditionalSigma = np.sqrt(conditionalVarianceSample)
    # Compute the quantile of the Normal distribution
    alpha = 1-(1-level)/2
    quantileAlpha = ot.DistFunc.qNormal(alpha)
    # Graphics of the bounds
    epsilon = 1.e-8
    dataLower = [yKrig[i,0] - quantileAlpha * conditionalSigma[i,0] for i in range(nbpoints)]
    dataUpper = [yKrig[i,0] + quantileAlpha * conditionalSigma[i,0] for i in range(nbpoints)]
    # Coordinates of the vertices of the Polygons
    # Index
    vLow = [[xGrid[i,Index],dataLower[i]] for i in range(nbpoints)]
    vUp = [[xGrid[i,Index],dataUpper[i]] for i in range(nbpoints)]
    # Compute the Polygon graphics
    boundsPoly = plot_kriging_bounds(vLow,vUp,nbpoints)
    boundsPoly.setLegend("95% bounds")
    # Validate the kriging metamodel
    #mmv = ot.MetaModelValidation(xGrid, yFunction, meta)
    #Q2 = mmv.computePredictivityFactor()[0]
    # Plot the function
    #graphFonction = ot.Curve(xGrid,yFunction)
    #graphFonction.setLineStyle("dashed")
    #graphFonction.setColor("magenta")
    #graphFonction.setLineWidth(2)
    #graphFonction.setLegend("Function")
    # Draw the X and Y observed
    cloudDOE = ot.Cloud(X, Y)
    cloudDOE.setPointStyle("circle")
    cloudDOE.setColor("black")
    cloudDOE.setLegend("Data")
    # Assemble the graphics
    graph = ot.Graph()
    graph.add(boundsPoly)
    #graph.add(graphFonction)
    graph.add(cloudDOE)
    graph.add(graphKriging)
    #graph.setLegendPosition("top")
    graph.setLegendFontSize(1.5)
    graph.setAxes(True)
    graph.setGrid(True)
    return graph

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
    numberOfInputParam = 3;
    # Total number of output results (4) 
    # i.e. min error, max error, L1 error, L2 error
    numberOfOutputParam = 4;

    # Extract useful information about arrays etc. Refer Run.py to find out the sequence of data
    sampleNumber = inputData[:,0]; # not used
    NHL = inputData[:,1]; # not used

    OPLTYPE = inputData[:,2];
    # Convert the OPLTYPE into 0,1,2,3 from OPLTYPE = [0,1,3,4]
    OPLTYPE[OPLTYPE == 3] = 2;
    OPLTYPE[OPLTYPE == 4] = 3;

    HL_0_DIM = inputData[:,3];
    HL_0_TYPE = inputData[:,4];
    # Convert the HL_0_TYPE into 0,1,2,3 from HL_0_TYPE = [0,1,3,4]
    HL_0_TYPE[HL_0_TYPE == 3] = 2;
    HL_0_TYPE[HL_0_TYPE == 4] = 3;

    HL_1_DIM = inputData[:,5];  # Will NOT be used
    HL_1_TYPE = inputData[:,6];  # Will NOT be used
    # Convert the HL_1_TYPE into 0,1,2,3 from HL_1_TYPE = [0,1,3,4]
    HL_1_TYPE[HL_1_TYPE == 3] = 2;
    HL_1_TYPE[HL_1_TYPE == 4] = 3;

    HL_2_DIM = inputData[:,7];  # Will NOT be used
    HL_2_TYPE = inputData[:,8];  # Will NOT be used
    # Convert the HL_2_TYPE into 0,1,2,3 from HL_2_TYPE = [0,1,3,4]
    HL_2_TYPE[HL_2_TYPE == 3] = 2;
    HL_2_TYPE[HL_2_TYPE == 4] = 3;


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

    observations[:,0] = L1Error/np.max(L1Error);
    observations[:,1] = MSError/np.max(MSError);
    observations[:,2] = MinError/np.max(MinError);
    observations[:,3] = MaxError/np.max(MaxError);





    input_names = ['OPLTYPE', 'HL0-D', 'HL0-A']

    output_names = ['L1Error', 'MSError','MinError','MaxError'];
    #_______________________________________________________
    # Find sobol indices
    #_______________________________________________________

    
    ot.RandomGenerator.SetSeed(int(1000*time.time()))
    distributionNHL = ot.Uniform(1,3);
    distributionOPLTYPE = ot.Uniform(0,3);
    distributionHL0DIM = ot.Uniform(5,15);
    distributionHL0TYPE = ot.Uniform(0,3);
    distributionHL1DIM = ot.Uniform(5,15);
    distributionHL1TYPE = ot.Uniform(0,3);
    distributionHL2DIM = ot.Uniform(5,15);
    distributionHL2TYPE = ot.Uniform(0,3);

    distributionList = [distributionOPLTYPE, distributionHL0DIM, distributionHL0TYPE];
    distribution = ot.ComposedDistribution(distributionList)


    for i in range(0,4):

        inputIndex = 2;
        # Get a fresh meta model for ith output
        print("Solving for ",output_names[i]," \n");

        print("1. Getting metamodel...");
        krigingResult, metamodel = getMetamodel(coordinates, observations,i);

        # Fixing number of neurons and activation functions on neurons
        x1Ref = 0; 
        x2Ref = 15;

        # Get another metamodel with two of the inputs fixed
        metamodelAtXref = ot.ParametricFunction(metamodel, [0, 1], [x1Ref, x2Ref])


        IPInd1 = np.where(coordinates[:,0] == x1Ref)[0];
        IPInd2 = np.where(coordinates[:,1] == x2Ref)[0];

        IPInd = np.intersect1d(IPInd1, IPInd2);

        IP = coordinates[IPInd, inputIndex];
        OP = observations[IPInd, i];

        # Fix a range for the plotting data 
        x2min = -0.25;
        x2max = 3.25;
        

        graph = plotBasicKriging(krigingResult,metamodelAtXref,x2min,x2max, IP, OP)
        view = viewer.View(graph);

        axes = view.getAxes();
        location = 0;
        axes[location].figure.set_size_inches(6.4,6.4)
        axes[location].figure.tight_layout(rect=[0.1,0.1,0.95,0.95])
        axes[location].legend(loc = 9, fontsize=18);
        axes[location].set_ylabel(output_names[i], fontsize=20);
        axes[location].set_xlabel(input_names[inputIndex], fontsize=20);
        axes[location].set_ylim(-0.5,2);
        axes[location].set_xlim(-0.25,3.25);
        axes[location].set_xticks([0,1,2,3]);
        axes[location].set_xticklabels(["AF1","AF2","AF3","AF4"], rotation=60);

        #plt.show();
        plt.savefig("KrigingExpt"+str(runNumber)+input_names[inputIndex]+output_names[i]+".pdf");

        pass;
    os.chdir(currDir);
    pass;


if __name__=="__main__":
    plt.rcParams["text.usetex"] = True
    plt.rcParams["font.family"] = "serif"
    plt.rcParams["font.serif"] = ["Computer Modern"]
    plt.rcParams['xtick.labelsize'] = 18
    plt.rcParams['ytick.labelsize'] = 18

    # Name of the project
    projectName = "Avg";

    # Data size for running Sobol analysis. i.e. these many samples will be generated using the metamodel
    size = 250000;

    #runNumber = 5;
    for runNumber in range(8):
        plotKrigingMetamodel(projectName, runNumber, size); 
