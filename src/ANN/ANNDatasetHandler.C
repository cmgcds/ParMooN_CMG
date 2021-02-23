/** ************************************************************************ 
* @brief     source file for TANNDatasetHandler
* @author    Subodh M. Joshi
* @date      06.01.21
* @History   Implementation started 24.12.20 (Subodh M. Joshi) 
 ************************************************************************  */

/** Class declaration is done in include/ANN/ANNDatasetHandler.h */

#include <ANNDatasetHandler.h>

TANNDatasetHandler::TANNDatasetHandler(TANNParamReader *paramReader){

  // Load the data
  mlpack::data::Load(paramReader->datasetName.c_str(), this->allData, true);

  

  // Find the total number of samples in the dataset
  this->totalNumberOfSamples = allData.n_cols;

  // Number of training samples based on the user input
  this->numberOfTrainingSamples = this->totalNumberOfSamples * paramReader->trainingDataPercentage / 100.;

  // Number of validation samples based on the user input
  this->numberOfValidationSamples = this->totalNumberOfSamples * paramReader->validationDataPercentage / 100.;


  // Find the testing dataset
  this->numberOfTestingSamples = this->totalNumberOfSamples - (this->numberOfTrainingSamples + this->numberOfValidationSamples);

  int numberOfTrainingAndValidationSamples = this->numberOfTrainingSamples + this->numberOfValidationSamples;

  // Set the number of epochs for training
  this->epochs = paramReader->epochs;

  // Dimension of the input data (number of features at the input layer per perceptron)
  this->ipDataDim = paramReader->ipDataDim;

  int totalInputLayerDim = this->ipDataDim * paramReader->layerDim[0];

  // Create the training dataset
  this->trainData = this->allData.submat(0,0, totalInputLayerDim-1, numberOfTrainingSamples-1)/paramReader->featureScalingConstant;

  // Create the training labels
  this->trainLabels = this->allData.submat(this->allData.n_rows-1,0, this->allData.n_rows-1,  numberOfTrainingSamples-1);

  // Create validation dataset
  this->validationData = this->allData.submat(0,numberOfTrainingSamples, totalInputLayerDim-1, numberOfTrainingAndValidationSamples-1)/paramReader->featureScalingConstant;

  // Create the validation labels
  this->validationLabels = this->allData.submat(this->allData.n_rows-1,numberOfTrainingSamples, this->allData.n_rows-1,  numberOfTrainingAndValidationSamples-1);

  // Create the testing dataset
  this->testData = this->allData.submat(0,numberOfTrainingAndValidationSamples, totalInputLayerDim-1, totalNumberOfSamples-1)/paramReader->featureScalingConstant;

  // Create the testing labels
  this->testLabels = this->allData.submat(this->allData.n_rows-1,numberOfTrainingAndValidationSamples, this->allData.n_rows-1,  totalNumberOfSamples-1);

};


// Export the Data 
// void TANNDatasetHandler::saveModel(){

//   // Load the data
//   mlpack::data::Save("model.xml", "model",  this->allData, false);

// }

TANNDatasetHandler::~TANNDatasetHandler(){};

void TANNDatasetHandler::postProcessResults(){
  if (predictionTemp.n_rows == 1){
    // Regression problem
    prediction = predictionTemp;
    errorL1Absolute = this->computeError(testLabels, prediction, "L1","ABS");
    errorL1Relative = this->computeError(testLabels, prediction, "L1","REL");
    errorL2Absolute = this->computeError(testLabels, prediction, "L2","ABS");
    errorL2Relative = this->computeError(testLabels, prediction, "L2","REL");
    errorLInfAbsolute = this->computeError(testLabels, prediction, "LInf","ABS");
    errorLInfRelative = this->computeError(testLabels, prediction, "LInf","REL");

    std::cout << " errorL1Absolute : " << errorL1Absolute << std::endl;
    std::cout << " errorL1Relative : " << errorL1Relative << std::endl;
    std::cout << " errorL2Absolute : " << errorL2Absolute << std::endl;
    std::cout << " errorL2Relative : " << errorL2Relative << std::endl;
    std::cout << " errorLInfAbsolute : " << errorLInfAbsolute << std::endl;
    std::cout << " errorLInfRelative : " << errorLInfRelative << std::endl;
    
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

  double numerator, denominator, temp;
  double error;
  // L1 norm of the error 
  if (norm == "L1"){
    numerator = 0.0;
    denominator = 0.0;
    for (int i=0; i < numberOfTestingSamples; i++){
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
    for (int i=0; i < numberOfTestingSamples; i++){
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
  else if (norm == "LInf"){
    std::ofstream file;
    file.open("result.txt");
    file << "Reference" << "," << "Actual" << "," << "Error" <<std::endl;
    numerator = 1000;
    denominator = 0.0;
    for (int i=0; i < numberOfTestingSamples; i++){
      file << referenceValue(i) <<","<< numericalValue(i) << "," << abs(referenceValue(i) - numericalValue(i)) <<std::endl;
      if ( abs(referenceValue(i) - numericalValue(i))  <= numerator){
        numerator = abs(referenceValue(i) - numericalValue(i));
      };
    };

    error = numerator;
    file.close();
  };

  return error;
};
    


