from dataFunctions import *

import RunExpt as exp


def run(exptName, exeName):
    # exptName: name of the experiment
    # exeName: name of the executable 

    TrainingSizes = [6,12,25,50,100,200,400,800];
    TestingSize = 50;
    ValidationPercentage = 50;


    # Create testing data
    createTestingDataset("allData", TestingSize);

    for trainingSize in TrainingSizes:
        # Create training Dataset
        createTrainingAndValidationDataset("allData", trainingSize, ValidationPercentage);

        exp.runExperiment(exptName,trainingSize, ValidationPercentage, exeName);
        pass;


if __name__=="__main__":
    run("TrialExpt", "parmoon_2D_SEQUENTIAL.exe");



