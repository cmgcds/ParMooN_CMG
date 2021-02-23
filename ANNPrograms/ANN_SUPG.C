// =======================================================================
// Purpose:     main program for ANN model for Thyroid dataset (classification)
// Author:      Subodh Joshi
//
// History:     Implementation started on 07.01.2020
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

  // Test
  ann.testNetwork(&datasetHandler);

  return 0;
} // end main
