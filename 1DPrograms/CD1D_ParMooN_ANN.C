// =======================================================================
// Purpose:     main program for 1D with new kernels of ParMooN
// Author:      Sashikumaar Ganesan
//
// History:     Implementation started on 12.12.2020
// =======================================================================

#include <ANN.h>
#include <ANNDatasetHandler.h>
// main program

int main(int argc, char* argv[])
{
  // scout << " HELOOOOOOOOOO "<<endl;
  TANNParamReader paramReader(argv[1]);

  TANNDatasetHandler datasetHandler(&paramReader);

  TANN<MEAN_SQUARED_ERROR,RANDOM_INITIALIZATION> ann(&paramReader);
  //TANN<MEAN_SQUARED_ERROR, RANDOM_INITIALIZATION> ann(&paramReader);

  ann.trainNetwork(&datasetHandler);

  // Save the model 
  // datasetHandler.saveModel();
  ann.testNetwork(&datasetHandler);

  return 0;
} // end main
