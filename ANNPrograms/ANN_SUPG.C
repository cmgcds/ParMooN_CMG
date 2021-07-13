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

  // Create a new ANN model (template arguments can be found in ./include/ANN/ANNIncludes.h)
  TANN<MEAN_SQUARED_ERROR,RANDOM_INITIALIZATION> ann(&paramReader);

  // Train
  ann.trainNetwork(&datasetHandler);

  // // Test
  ann.testNetwork(&datasetHandler);
  // std::cout << " TESTING COMPLETED " <<std::endl;
  // // Save Model  
  ann.saveModel();

  // // Load External Model
  // ann.loadExternalModel();

  // ann.testNetwork(&datasetHandler);

  // double b,eps,h;
  // b = 1;
  // eps = 1e-8;
  // h = 0.01;
  // // Predict the Valuees
  // double predictedTau = ann.predictTau(b,eps,h,&datasetHandler);

  // double Pe_nr  = (fabs(b)*h)/(2*eps);

  // double analyticalTau = h/(2*fabs(b)) * ( 1.0/tanh(Pe_nr) - 1.0/Pe_nr );


  // std::cout << " PREDICTED : " << predictedTau << "\n ANALYTICAL : " << analyticalTau <<"\n";
  // std::cout << " ABS DIFF : " << fabs(predictedTau - analyticalTau)<< " \n";
  return 0;
} // end main
