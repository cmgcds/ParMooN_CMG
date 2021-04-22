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

  /** Optimizer code */
  /** 0: Default (RMSProp)
   *  1: SGD
   *  2: Adam
   *  **/
  int optimizerCode;

  /** Step size for the optimizer */
  double optimizerStepSize;

  /** Stochastic GD batch size */
  int sgdBatchSize;

  /** Max number of iterations */
  int maxIterations;

  int epochs;

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

  // Verify results with training and validation datasets
  void verifyModel(TANNDatasetHandler *datasetHandler);



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
    int effectiveLayerDim;
    if (i == 0){
      // INPUT Layer
      effectiveLayerDim = ipDataDim * paramReader->layerDim[0];
    }
    else{
      // Next layers
      effectiveLayerDim =  paramReader->layerDim[i];
    };

    // Note: for i=0, the activation 'Dummy' will be automatically set
    layers[i] = TANNLayer(i, effectiveLayerDim, paramReader->layerType[i], paramReader->layerTypeInt[i]);
  };

  layerFlag = true;


  // Set the loss functions
  this->lossFunction = paramReader->lossFunction;
  this->lossFunctionInt = paramReader->lossFunctionInt;

  /** Set up the model */
  this->setupModel();

  this->optimizerCode = paramReader->optimizerCode;

  this->optimizerStepSize = paramReader->optimizerStepSize;

  this->sgdBatchSize = paramReader->sgdBatchSize;

  this->maxIterations = paramReader->maxIterations;

  this->epochs = paramReader->epochs;

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
  // Create the activation layers
  /**
   *              
   *   i = 0       1           2         3
   *
   *          |    O     |          | 
   *     O    |    O     |     O    | 
   *     O    |    O     |          |    O
   *     O    |    O     |     O    | 
   *          |    O     |          | 
   *    ______________________________________
   *    IP    |    H1    |     H2   |    OP  
   *          |          |          | 
   *          |          |          | 
   *       First       Second     Third 
   *       Activation            
   *
   * Assume that all the activation functions are ReLU except for the last one.
   * For the above example, the arrays woule look like:
   * layers[:].dim = [3, 5, 2, 1]
   * layers[:].type = [Dummy, ReLU, ReLU, LogSoftMax]
   * lyaers[:].typeInt = [Dummy, 2, 2, 5]
   *
   * **/
                                      
  for (int i=0; i<nLayers-1; i++){

    // inSize is the dim of the Neural layer on the LHS of the activation layer
    inSize = this->layers[i].dim;

    // outSize is the dim of the Neural Layer on the RHS of the activation layer 
    outSize = this->layers[i+1].dim;


    // Add regularization before the output layer
    if (i == nLayers-2){
      this->model.template Add<mlpack::ann::Dropout<>>(0.2);
    };

    // Add a linear Layer (i.e. Sum(x_i * omega_i) )
    this->model.template Add<mlpack::ann::Linear<> >(inSize, outSize);

    // Apply the nonlinear activation function
    // Note: 
    // For the first layer (Input), the activation Type and TypeInt were set to
    // a dummy value (refer the constructor as well as the diagram above)
    // Hence we skip the first entry from this->layers[0].typeInt and start reading 
    // from i+1
    this->applyActivation(this->layers[i+1].typeInt);
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

    case 5:
      this->model.template Add<mlpack::ann::LogSoftMax<> >();
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
void TANN<OutputLayerType, InitializationRuleType, CustomLayers...>::verifyModel(TANNDatasetHandler *datasetHandler){
      // Test the predictions on the training data
      arma::mat predOutTrain;
      this->model.template Predict(datasetHandler->trainData, predOutTrain);

      if (predOutTrain.n_rows == 1){
          // Regression model

          double errorTrain = datasetHandler->computeError(datasetHandler->trainLabels, predOutTrain, "L2", "REL");
          std::cout << " Train error (L2, Rel) : " << errorTrain << std::endl;

          arma::mat predOutValid;

          this->model.template Predict(datasetHandler->validationData, predOutValid);

          double errorValid = datasetHandler->computeError(datasetHandler->validationLabels, predOutValid, "L2", "REL");
          std::cout << " Validation error (L2, Rel) : " << errorValid << std::endl;

      }
      else{
          // Classification model
            
          // Get the prediction results in a single column vector
          arma::mat predictionTrain = arma::zeros<arma::mat>(1, predOutTrain.n_cols);
          for (int i = 0; i < predOutTrain.n_cols; ++i){
            // +1 because the labels start from 1 and not 0
            predictionTrain(i) = predOutTrain.col(i).index_max() + 1;
          };
          // Calculating accuracy on training data points.
          double trainAccuracy =
              arma::accu(predictionTrain == datasetHandler->trainLabels) / (double) datasetHandler->trainLabels.n_elem * 100;

          // Test the predictions on the validation data
          arma::mat predOutValid;
          this->model.template Predict(datasetHandler->validationData, predOutValid);
          // Get the prediction results in a single column vector
          arma::mat predictionValid = arma::zeros<arma::mat>(1, predOutValid.n_cols);
          for (int i = 0; i < predOutValid.n_cols; ++i){
            // +1 because the labels start from 1 and not 0
            predictionValid(i) = predOutValid.col(i).index_max() + 1;
          };
          // Calculating accuracy on validation data points.
          double validAccuracy =
              arma::accu(predictionValid == datasetHandler->validationLabels) / (double) datasetHandler->validationLabels.n_elem * 100;

          std::cout << "Accuracy: train = " << trainAccuracy << "%," << "\t valid = " << validAccuracy << "%" << std::endl;
      };
};
  

template<
  typename OutputLayerType ,
  typename InitializationRuleType, 
  typename... CustomLayers
>
void TANN<OutputLayerType, InitializationRuleType, CustomLayers...>::trainNetwork(TANNDatasetHandler *datasetHandler){
  
  switch(this->optimizerCode){
    case 0:
    {
      // Default optimizer (i.e. RMSProp)
      this->model.template Train(datasetHandler->trainData, datasetHandler->trainLabels, ens::PrintLoss(), ens::ProgressBar());
      this->verifyModel(datasetHandler);

      break;
    };

    case 1:
    {
      // SGD optimizer 
      // Set parameters for the SGD optimizer.
      ens::StandardSGD optimizer(this->optimizerStepSize, this->sgdBatchSize, this->maxIterations);
      // Declare callback to store best training weights.
      ens::StoreBestCoordinates<arma::mat> bestCoordinates;

      // Train neural network. If this is the first iteration, weights are
      // random, using current values as starting point otherwise.
      model.Train(datasetHandler->trainData,
                  datasetHandler->trainLabels,
                  optimizer,
                  //ens::PrintLoss(),
                  ens::ProgressBar(),
                  ens::EarlyStopAtMinLoss(),
                  bestCoordinates);

      // Save the best training weights into the model.
      model.Parameters() = bestCoordinates.BestCoordinates();

      this->verifyModel(datasetHandler);

      break;
    };

    case 2:
    {
      // Adam optimizer
      // Set parameters for the Adam optimizer.
      ens::Adam optimizer(this->optimizerStepSize, this->sgdBatchSize, this->maxIterations);
      // Declare callback to store best training weights.
      ens::StoreBestCoordinates<arma::mat> bestCoordinates;

      // Train neural network. If this is the first iteration, weights are
      // random, using current values as starting point otherwise.
      model.Train(datasetHandler->trainData,
                  datasetHandler->trainLabels,
                  optimizer,
                  //ens::PrintLoss(),
                  ens::ProgressBar(),
                  ens::EarlyStopAtMinLoss(),
                  bestCoordinates);

      // Save the best training weights into the model.
      model.Parameters() = bestCoordinates.BestCoordinates();

      this->verifyModel(datasetHandler);

      break;
    };

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
