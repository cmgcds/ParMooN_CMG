/** ************************************************************************ 
* @brief     source file for TANNDatasetHandler
* @author    Subodh M. Joshi
* @date      06.01.21
* @History   Implementation started 24.12.20 (Subodh M. Joshi) 
 ************************************************************************  */

/** Class declaration is done in include/ANN/ANNDatasetHandler.h */

#include <ANNDatasetHandler.h>

TANNDatasetHandler::TANNDatasetHandler(TANNParamReader *paramReader){

  // Load the training and validation data
  // NOTE: allData contains the training as well as the validation data
  mlpack::data::Load(paramReader->trainingDatasetName.c_str(), this->allData, true);

  // Find the total number of samples in the dataset
  this->totalNumberOfSamples = allData.n_cols;

  // Number of training samples based on the user input
  this->numberOfValidationSamples = this->totalNumberOfSamples * paramReader->validationDataPercentage / 100.;

  // Number of validation samples based on the user input
  this->numberOfTrainingSamples = this->totalNumberOfSamples - this->numberOfValidationSamples;


  // Find the testing dataset
  mlpack::data::Load(paramReader->testingDatasetName.c_str(), this->testDataset, true);

  this->numberOfTestingSamples =  testDataset.n_cols;

  std::cout << "No. of training samples: " << this->numberOfTrainingSamples << "  No. of validation samples: " << this->numberOfValidationSamples << "   No. of testing samples: " << numberOfTestingSamples << std::endl;


  //int numberOfTrainingAndValidationSamples = this->numberOfTrainingSamples + this->numberOfValidationSamples;

  // Set the number of epochs for training
  this->epochs = paramReader->epochs;

  // Set the file name for storing the results for the tests
  this->saveDataFile = paramReader->saveDataFile;

  // Dimension of the input data (number of features at the input layer per perceptron)
  this->ipDataDim = paramReader->ipDataDim;

  int totalInputLayerDim = this->ipDataDim * paramReader->layerDim[0];

  //__________________________________________
  //1. Create datasets 
  //__________________________________________

  // Create the training dataset
  this->trainData = this->allData.submat(0,0, totalInputLayerDim-1, numberOfTrainingSamples-1);

  // Create the training labels
  this->trainLabels = this->allData.submat(this->allData.n_rows-1,0, this->allData.n_rows-1,  numberOfTrainingSamples-1);

  // Create validation dataset
  this->validationData = this->allData.submat(0,numberOfTrainingSamples, totalInputLayerDim-1, totalNumberOfSamples-1);

  // Create the validation labels
  this->validationLabels = this->allData.submat(this->allData.n_rows-1,numberOfTrainingSamples, this->allData.n_rows-1, totalNumberOfSamples-1);

  // Create the testing dataset
  this->testData = this->testDataset.submat(0,0, totalInputLayerDim-1, numberOfTestingSamples-1);

  // Create the testing labels
  this->testLabels = this->testDataset.submat(this->allData.n_rows-1,0, this->allData.n_rows-1,  numberOfTestingSamples-1);

  //__________________________________________
  //2. Train the scaler to scale the data
  //__________________________________________
  // Train the scaler based on the training data
  this->dataScaler.Fit(this->trainData);

  // Train the labelScaler 
  this->labelScaler.Fit(this->trainLabels);

  //__________________________________________
  //3. Scale the dataset and the labels
  //__________________________________________

  // Scale the training data
  this->dataScaler.Transform(this->trainData, this->trainData);

  // Scale the training Labels
  this->labelScaler.Transform(this->trainLabels, this->trainLabels);

  // Scale the validation dataset
  dataScaler.Transform(this->validationData, this->validationData);

  // Scale the validation labels
  this->labelScaler.Transform(this->validationLabels, this->validationLabels);

  // Scale the test data
  this->dataScaler.Transform(this->testData, this->testData);

  // Scale the testing labels
  this->labelScaler.Transform(this->testLabels, this->testLabels);

};


// Export the Data 
void TANNDatasetHandler::saveModel(){

  // save the data
  mlpack::data::Save("model.txt", "model",  this->allData, false);

}

TANNDatasetHandler::~TANNDatasetHandler(){};

void TANNDatasetHandler::postProcessResults(){
  if (predictionTemp.n_rows == 1){
    // Regression problem
    prediction = predictionTemp;

    // First scale back the result values
    // NOTE: The values were scaled when the arrays were created in the constructor.
    this->labelScaler.InverseTransform(this->prediction, this->prediction);
    this->labelScaler.InverseTransform(this->testLabels, this->testLabels);

    errorL1Absolute = this->computeError(testLabels, prediction, "L1","ABS");
    errorL1Relative = this->computeError(testLabels, prediction, "L1","REL");
    errorL2Absolute = this->computeError(testLabels, prediction, "L2","ABS");
    errorL2Relative = this->computeError(testLabels, prediction, "L2","REL");
    errorLInfAbsolute = this->computeError(testLabels, prediction, "LInf","ABS");
    errorL0Absolute = this->computeError(testLabels, prediction, "L0","ABS");
    errorMSE = this->computeError(testLabels, prediction, "MSE", "ABS");

    std::cout << " Test results: " << std::endl;

    std::cout << " errorL1Absolute : " << errorL1Absolute << std::endl;
    std::cout << " errorL1Relative : " << errorL1Relative << std::endl;
    std::cout << " errorL2Absolute : " << errorL2Absolute << std::endl;
    std::cout << " errorL2Relative : " << errorL2Relative << std::endl;
    std::cout << " errorLInfAbsolute : " << errorLInfAbsolute << std::endl;
    std::cout << " Min Error Absolute : " << errorL0Absolute << std::endl;
    
    writeErrorFile(testLabels, prediction);
  }
  else{
    // Classification problem 
    // std::cout << " PRED TEMP : " << predictionTemp.n_rows <<std::endl;
    std::cout << "Classification problem " << std::endl;
    prediction = arma::zeros<arma::mat>(1, predictionTemp.n_cols);
      // Find index of max prediction for each data point and store in "prediction"
    for (int i = 0; i < predictionTemp.n_cols; ++i)
    {
      // we add 1 to the max index, so that it matches the actual test labels.
      prediction(i) = arma::as_scalar(arma::find(
          arma::max(predictionTemp.col(i)) == predictionTemp.col(i), 1)) + 1;
    };

    // Find the error
      /*
      Compute the error between predictions and testLabels,
      now that we have the desired predictions.
    */
    
    int correct = arma::accu(prediction == testLabels);

    this->errorL1Relative = 1 - double(correct) / testData.n_cols;
  };

};

double TANNDatasetHandler::computeError(arma::mat referenceValue, arma::mat numericalValue, std::string norm = "L1", std::string type = "ABS"){

  int arraySize = referenceValue.n_elem;
  double numerator, denominator, temp;
  double error;

  if (norm == "MSE"){
    return mlpack::metric::SquaredEuclideanDistance::Evaluate(referenceValue, numericalValue) / (numericalValue.n_elem);
  }

  else if (norm == "L1"){
    numerator = 0.0;
    denominator = 0.0;
    for (int i=0; i < arraySize; i++){
      numerator += abs(referenceValue(i) - numericalValue(i));
      if (type == "ABS"){
        denominator = 1;
      }
      else if (type == "REL"){
        denominator += abs(referenceValue(i));
      }
      else{
        std::cout << "Wrong argument in computeError for type" << std::endl;
      };
    };

    error = numerator / denominator;
  }

  else if (norm == "L2"){
    numerator = 0.0;
    denominator = 0.0;
    for (int i=0; i < arraySize; i++){
      numerator += pow(referenceValue(i) - numericalValue(i), 2);
      if (type == "ABS"){
        denominator = 1;
      }
      else if (type == "REL"){
        denominator += pow(referenceValue(i), 2);
      }
      else{
        std::cout << "Wrong argument in computeError for type" << std::endl;
      };
    };

    error = sqrt(numerator / denominator);
  }
  else if (norm == "L0"){
    numerator = 1000;
    denominator = 0.0;
    for (int i=0; i < arraySize; i++){
      if ( abs(referenceValue(i) - numericalValue(i))  <= numerator){
        numerator = abs(referenceValue(i) - numericalValue(i));
      };
    };

    error = numerator;
  }
  else if (norm == "LInf"){
    std::ofstream file;
    numerator = 0;
    denominator = 0.0;
    for (int i=0; i < arraySize; i++){
      if ( abs(referenceValue(i) - numericalValue(i))  > numerator){
        numerator = abs(referenceValue(i) - numericalValue(i));
      };
    };

    error = numerator;
    file.close();
  };

  return error;
};
    

void TANNDatasetHandler::writeErrorFile(arma::mat referenceValue, arma::mat numericalValue){
    // Print results into a file 
    std::ofstream file;
    file.open(this->saveDataFile);
    file << "#Test results \n#Reference" << "," << "Actual" << "," << "Abs Error" <<std::endl;
    int arraySize = referenceValue.n_elem;
    for (int i=0; i < arraySize; i++){
      file << referenceValue(i) <<","<< numericalValue(i) << "," << abs(referenceValue(i) - numericalValue(i)) <<std::endl;
      };

    file.close();
    };


