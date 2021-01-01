/** ************************************************************************ 
* @brief     source file for TANN
* @author    Subodh M. Joshi
* @date      24.12.20
* @History   Implementation started 24.12.20 (Subodh M. Joshi) 
 ************************************************************************  */

/** Class declaration is done in include/ANN/ANN.h */

#include <ANN.h>

/** Constructor 1 **/
TANN::TANN(){};

/** Constructor 2 **/
TANN::TANN(TANNParamReader *paramReader){
  /** Set all the flags false by default */
  layerFlag = false;
  modelFlag = false;

  /** Set data members */
  // Set up the number of hidden layers
  nHL = paramReader->nHL;

  // Set up the total number of layers
  nLayers = nHL + 2;

  /** Set up layers array */
  // Assign space for the layers array
  layers = new TANNLayer[nLayers];

  // Call constructor for each layer
  for (int i=0; i<nLayers; i++){
    layers[i] = TANNLayer(i, paramReader->layerDim[i], paramReader->layerType[i]);
  };

  layerFlag = true;


  /** Set up the model */

};

/** Destructor **/
TANN::~TANN(){
  assert(layerFlag == true and "layerFlag not activated in ANN creation");
  delete[] layers;
};


