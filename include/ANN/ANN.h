// =======================================================================
// @(#)TANN.h        1.0 24.12.20
// 
// Class:       TANN
// Purpose:     Class for Artificial Neural Networks (ANN) - Base class
//
// Author:      Subodh M. Joshi (24.12.20)
//
// History:     start of implementation 24.12.20 (Subodh M. Joshi)
//
// =======================================================================

#ifndef __TANN__
#define __TANN__

#include <ANNIncludes.h>
#include <ANNLayer.h>
#include <ANNDatasetHandler.h>

/** Class for artificial neural networks (ANN) **/
template<
  typename OutputLayerType = mlpack::ann::NegativeLogLikelihood<>,
  typename InitializationRuleType = mlpack::ann::RandomInitialization,
  typename... CustomLayers
>
class TANN
{
  public:
  /** constructor 1 */  
  TANN();

  /** constructor 2 */
  TANN(TANNParamReader *paramReader);

  /** destrcutor */
  ~TANN();

  public:
  /** Data members */

  TANNLayer* layers; // array storing individual layers

  int nHL;  // number of hidden layers

  int nLayers; // total number of layers

  int ipDataDim; // Number of inputs per neuron at the input layer

  std::string lossFunction; // Loss function in string form

  int lossFunctionInt; // int. encoding the loss function

  mlpack::ann::FFN<OutputLayerType, InitializationRuleType> model; // mlpack model for feed forward network, initialized with proper arguments in the constructor of this class

  public:
  /** Methods */
  void trainNetwork(TANNDatasetHandler *datasetHandler); 

  void testNetwork(TANNDatasetHandler *datasetHandler);

  private:
  /** Methods */

  // Function to set up the mlpack model
  void setupModel();

  // Function to apply the activation function
  void applyActivation(int activation);


  private:
  /** Flags */
  bool layerFlag;
  bool modelFlag;

};


/** Constructor 1 */
template<
  typename OutputLayerType , 
  typename InitializationRuleType, 
  typename... CustomLayers
>
TANN<OutputLayerType, InitializationRuleType, CustomLayers...>::TANN(){};

/** Constructor 2 **/
template<
  typename OutputLayerType ,
  typename InitializationRuleType, 
  typename... CustomLayers
>
TANN<OutputLayerType, InitializationRuleType, CustomLayers...>::TANN(TANNParamReader *paramReader){
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


  // Set the loss functions
  this->lossFunction = paramReader->lossFunction;
  this->lossFunctionInt = paramReader->lossFunctionInt;

  /** Set up the model */
  this->setupModel();

};

/** Destructor */
template<
  typename OutputLayerType ,
  typename InitializationRuleType, 
  typename... CustomLayers
>
TANN<OutputLayerType, InitializationRuleType, CustomLayers...>::~TANN(){
  assert(layerFlag == true and "layerFlag not activated in ANN creation");
  delete[] layers;
};

/** Setup the model */
template<
  typename OutputLayerType ,
  typename InitializationRuleType, 
  typename... CustomLayers
>
void TANN<OutputLayerType, InitializationRuleType, CustomLayers...>::setupModel(){

  int inSize, outSize;
  for (int i=0; i<nLayers; i++){

    if (i == 0){
      inSize = ipDataDim;
    }
    else {
      inSize = this->layers[i-1].dim;
    };

    outSize = this->layers[i].dim;

    this->model.template Add<mlpack::ann::Linear<> >(inSize, outSize);
    this->applyActivation(this->layers[i].typeInt);

  };

};


template<
  typename OutputLayerType ,
  typename InitializationRuleType, 
  typename... CustomLayers
>
void TANN<OutputLayerType, InitializationRuleType, CustomLayers...>::applyActivation(int activation){
  switch (activation){
    case 0:
      this->model.template Add<mlpack::ann::SigmoidLayer<> >();
      break;

    case 1:
      this->model.template Add<mlpack::ann::IdentityLayer<> >();
      break;

    case 2:
      this->model.template Add<mlpack::ann::ReLULayer<> >();
      break;

    case 3:
      this->model.template Add<mlpack::ann::TanHLayer<> >();
      break;

    case 4:
      this->model.template Add<mlpack::ann::SoftPlusLayer<> >();
      break;

    default:
      this->model.template Add<mlpack::ann::SigmoidLayer<> >();
      break;


  };
};

template<
  typename OutputLayerType ,
  typename InitializationRuleType, 
  typename... CustomLayers
>
void TANN<OutputLayerType, InitializationRuleType, CustomLayers...>::trainNetwork(TANNDatasetHandler *datasetHandler){
  
  for (int i=0; i<datasetHandler->epochs; i++){
    this->model.template Train(datasetHandler->trainData, datasetHandler->trainLabels);
  };
};

template<
  typename OutputLayerType ,
  typename InitializationRuleType, 
  typename... CustomLayers
>
void TANN<OutputLayerType, InitializationRuleType, CustomLayers...>::testNetwork(TANNDatasetHandler *datasetHandler){
  this->model.template Predict(datasetHandler->testData, datasetHandler->predictionTemp);
  datasetHandler->postProcessResults();
  std::cout << "Error : " << datasetHandler->errorL1Relative << std::endl;
};
#endif
