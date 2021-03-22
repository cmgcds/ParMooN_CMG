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

import RunExpt as exp


TrainingSizes = [12,25,50,100,200,400,800];
TestingSize = 50;
ValidationPercentage = 10;

name = "Expt2"


# Create testing data
createTestingDataset("allData", TestingSize);

for trainingSize in TrainingSizes:
    # Create training Dataset
    createTrainingAndValidationDataset("allData", trainingSize, ValidationPercentage);

    exp.runExperiment(name,trainingSize);
    pass;




