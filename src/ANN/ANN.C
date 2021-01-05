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

  // Data size at the input layer. i.e. number of inputs per neuron at the input layer.
  // This could be considered as the r,g,b,intensity values per neuron in image proc.
  ipDataDim = paramReader->ipDataDim;

  /** Set up layers array */
  // Assign space for the layers array
  layers = new TANNLayer[nLayers];

  // Call constructor for each layer
  for (int i=0; i<nLayers; i++){
    layers[i] = TANNLayer(i, paramReader->layerDim[i], paramReader->layerType[i], paramReader->layerTypeInt[i]);
  };

  layerFlag = true;


  /** Set up the model */
  this->setupModel();

};

/** Destructor **/
TANN::~TANN(){
  assert(layerFlag == true and "layerFlag not activated in ANN creation");
  delete[] layers;
};


void TANN::setupModel(){

  int inSize, outSize;
  for (int i=0; i<nLayers; i++){

    if (i == 0){
      inSize = ipDataDim;
    }
    else {
      inSize = this->layers[i-1].dim;
    };

    outSize = this->layers[i].dim;

    this->model.Add<mlpack::ann::Linear<> >(inSize, outSize);
    this->applyActivation(this->layers[i].typeInt);
  };

};

void TANN::applyActivation(int activation){
  switch (activation){
    case 0:
      this->model.Add<mlpack::ann::SigmoidLayer<> >();
      break;

    case 1:
      this->model.Add<mlpack::ann::IdentityLayer<> >();
      break;

    case 2:
      this->model.Add<mlpack::ann::ReLULayer<> >();
      break;

    case 3:
      this->model.Add<mlpack::ann::TanHLayer<> >();
      break;

    case 4:
      this->model.Add<mlpack::ann::SoftPlusLayer<> >();
      break;

    default:
      this->model.Add<mlpack::ann::SigmoidLayer<> >();
      break;


  };
};
