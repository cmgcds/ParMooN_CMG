// =======================================================================
// Purpose:     main program for ANN model for Regression problem
// Author:      Subodh Joshi
//
// History:     Implementation started on 25.02.2020
// =======================================================================

#include <ANN.h>
#include <ANNDatasetHandler.h>
// main program

int main(int argc, char* argv[])
{

  // Create a new parameter reader (takes argument .dat file)
  TANNParamReader paramReader(argv[1]);

  paramReader.print();
  // Create a new dataset handler (to create train and test datasets and labels)
  TANNDatasetHandler datasetHandler(&paramReader);

  // Create a linear regression model
  mlpack::regression::LinearRegression lr(datasetHandler.trainData, datasetHandler.trainLabels);

  arma::rowvec predOutTrain;

  lr.Predict(datasetHandler.trainData, (arma::rowvec&)predOutTrain);

  double errorTrain = datasetHandler.computeError(datasetHandler.trainLabels, predOutTrain, "MSE", "ABS");


  arma::mat predOutValid;

  lr.Predict(datasetHandler.validationData, (arma::rowvec&)predOutValid);

  double errorValid = datasetHandler.computeError(datasetHandler.validationLabels, predOutValid, "MSE", "ABS");

  std::cout << " MSE for training data: " << errorValid << std::endl;
  std::cout << " MSE for validation data: " << errorValid << std::endl;

  lr.Predict(datasetHandler.testData, (arma::rowvec&)datasetHandler.predictionTemp);
  datasetHandler.postProcessResults();

  return 0;
} // end main
