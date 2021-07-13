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
   *  1: Gradient Descent
   *  2: SGD
   *  3: Adam
   *  4: L-BFGS
   *  **/
  int optimizerCode;

  /** Step size for the optimizer */
  double optimizerStepSize;

  /** Stochastic GD batch size */
  int sgdBatchSize;

  /** Max number of iterations */
  int maxIterations;

  /** Number of training epochs */
  int epochs;

  /** Dropout ratio for the regularization before the output layer */
  double dropoutRatio;

  /** Tolerance (absolute) */
  double tolerance;

  /** Save model to this file after training */
  std::string saveModelFile;

  /** Load model from this file if the loadModelFlag is 1 */
  std::string loadModelFile;

  /** Load model flag. If 1, load the model from a file instead of training */
  std::string loadModelFlag;

  public:
  /** Methods */
  void trainNetwork(TANNDatasetHandler *datasetHandler); 

  void testNetwork(TANNDatasetHandler *datasetHandler);

  // Function to save the model to the external file
  void saveModel();

  // FUnction to load the model from external source
  void loadExternalModel(char* filename);

  // Predict the Tau Value
  double predictTau(double b, double eps, double h,TANNDatasetHandler* dataHandler);

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

  this->dropoutRatio = paramReader->dropoutRatio;

  this->tolerance = paramReader->tolerance;

  this->saveModelFile = paramReader->saveModelFile;

  this->loadModelFile = paramReader->loadModelFile;

  this->loadModelFlag = paramReader->loadModelFlag;
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
      this->model.template Add<mlpack::ann::Dropout<>>(this->dropoutRatio);
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
      this->model.template Add<mlpack::ann::LeakyReLU<> >();

    case 4:
      this->model.template Add<mlpack::ann::TanHLayer<> >();
      break;

    case 5:
      this->model.template Add<mlpack::ann::SoftPlusLayer<> >();
      break;

    case 6:
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

          double errorTrain = datasetHandler->computeError(datasetHandler->trainLabels, predOutTrain, "MSE", "ABS");
          std::cout << " Train error (MSE) : " << errorTrain << std::endl;

          arma::mat predOutValid;

          this->model.template Predict(datasetHandler->validationData, predOutValid);

          double errorValid = datasetHandler->computeError(datasetHandler->validationLabels, predOutValid, "MSE", "ABS");
          std::cout << " Validation error (MSE) : " << errorValid << std::endl;

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
  // Set the max iterations based on epochs

  if (this->maxIterations != 0){
    if (this->maxIterations < datasetHandler->trainData.n_cols * this->epochs){
      this->maxIterations = datasetHandler->trainData.n_cols * this->epochs;
    };
  };
  
  switch(this->optimizerCode){
    case 0:
    {
      // Default optimizer (i.e. RMSProp)
      this->model.template Train(datasetHandler->trainData, datasetHandler->trainLabels);
      this->verifyModel(datasetHandler);

      break;
    };

    case 1:
    {
      // Gradient Descent
      ens::GradientDescent optimizer(this->optimizerStepSize, this->maxIterations, this->tolerance);

      // Train neural network. If this is the first iteration, weights are
      // random, using current values as starting point otherwise.
      model.Train(datasetHandler->trainData,
                  datasetHandler->trainLabels,
                  optimizer
                  );

      this->verifyModel(datasetHandler);

      break;
    };


    case 2:
    {
      // SGD optimizer 
      // Set parameters for the SGD optimizer.
      ens::StandardSGD optimizer(this->optimizerStepSize, this->sgdBatchSize, this->maxIterations, this->tolerance);
      // Declare callback to store best training weights.
      ens::StoreBestCoordinates<arma::mat> bestCoordinates;

      // Train neural network. If this is the first iteration, weights are
      // random, using current values as starting point otherwise.
      model.Train(datasetHandler->trainData,
                  datasetHandler->trainLabels,
                  optimizer,
                  ens::EarlyStopAtMinLoss(),
                  bestCoordinates);

      // Save the best training weights into the model.
      model.Parameters() = bestCoordinates.BestCoordinates();

      this->verifyModel(datasetHandler);

      break;
    };

    case 3:
    {
      // Adam optimizer
      // Set parameters for the Adam optimizer.
      ens::Adam optimizer(this->optimizerStepSize, this->sgdBatchSize, 0.9, 0.999, 1e-8, this->maxIterations,this->tolerance,true);
      // Declare callback to store best training weights.
      ens::StoreBestCoordinates<arma::mat> bestCoordinates;

      // Train neural network. If this is the first iteration, weights are
      // random, using current values as starting point otherwise.
      model.Train(datasetHandler->trainData,
                  datasetHandler->trainLabels,
                  optimizer,
                  ens::EarlyStopAtMinLoss(20),
                  bestCoordinates);

      // Save the best training weights into the model.
      model.Parameters() = bestCoordinates.BestCoordinates();

      this->verifyModel(datasetHandler);

      break;
    };

    case 4:
    {
      // L-BFGS algorithm
      // Number of basis dimensions 
      int numBasis = 10;
      ens::L_BFGS optimizer(numBasis, this->maxIterations);

      // Train neural network. If this is the first iteration, weights are
      // random, using current values as starting point otherwise.
      model.Train(datasetHandler->trainData,
                  datasetHandler->trainLabels,
                  optimizer
                  );

      this->verifyModel(datasetHandler);

      break;
    };

    //std::cout << "\nNetwork training done. Saving the model in " << this->saveModelFile << std::endl;
    //mlpack::data::Save(saveModelFile, "Feed Forward Neural Network",this->model);
  };
};

template<
  typename OutputLayerType ,
  typename InitializationRuleType, 
  typename... CustomLayers
>
void TANN<OutputLayerType, InitializationRuleType, CustomLayers...>::testNetwork(TANNDatasetHandler *datasetHandler){
  datasetHandler->testData.print(" TEST DATA : ");
  this->model.template Predict(datasetHandler->testData, datasetHandler->predictionTemp);
  datasetHandler->postProcessResults();
  std::cout << "Test error (MSE): " << datasetHandler->errorMSE << std::endl;
};


//Load the model from the external Saved 
template<
  typename OutputLayerType ,
  typename InitializationRuleType, 
  typename... CustomLayers
>
void TANN<OutputLayerType, InitializationRuleType, CustomLayers...>::saveModel(){
  
  bool save = mlpack::data::Save("model.txt", "model", model, false);

  if(save)
    std::cout << "Model saved successfully " << std::endl;
  else
    std::cout <<"[ERROR] :  Model not saved "<< std::endl;
};


template<
  typename OutputLayerType ,
  typename InitializationRuleType, 
  typename... CustomLayers
>
void TANN<OutputLayerType, InitializationRuleType, CustomLayers...>::loadExternalModel(char* filename){
  
  bool successState = mlpack::data::Load(filename, "model", model);
  
  if(successState)
    std::cout << "Model Loaded successfully into the ANN class " << std::endl;
  else
    std::cout <<"[ERROR] :  Model not loaded "<< std::endl;
};


// Predict the Tau Value
template<
  typename OutputLayerType ,
  typename InitializationRuleType, 
  typename... CustomLayers
>
double TANN<OutputLayerType, InitializationRuleType, CustomLayers...>::predictTau(double b, double eps, double h,TANNDatasetHandler* datasetHandler){
  
  // Create a Arma Matrix
  arma::mat inputTestData(3,1);
  inputTestData(0) = b;
  inputTestData(1) = eps;
  inputTestData(2) = h;
  // std::cout << "INPUT TEST DATA :  roes : " << inputTestData.n_rows << "  col: " << inputTestData.n_cols  << " " << inputTestData(0) << "   " << inputTestData(1) << "  " << inputTestData(2)<<std::endl;
  // inputTestData.print(" INPUT TEST DATA : ");


  // std::cout << " DATASET HANDLER : " << datasetHandler->testData.n_rows << " col : "<< datasetHandler->testData.n_cols <<std::endl;
  arma::mat outputResult;


  // TRANSFORM THE PREDICT DATA based on scaling of training
  datasetHandler->dataScaler.Transform(inputTestData, inputTestData);

  //  inputTestData.print(" INPUT TEST DATA : ");

  // std::cout << "[SUCCESS] :  Prediction Before  \n"; 

  // Predict the data using given Values
  this->model.template Predict(inputTestData, outputResult);

  // std::cout << " Prediction Scalled before  : " << outputResult <<std::endl;
  // outputResult.print("outputResult:");

  // std::cout << "[SUCCESS] :  Prediction complted  \n";
  // Scale the values back
  datasetHandler->labelScaler.InverseTransform(outputResult, outputResult);

  double PredictedTau = outputResult(0,0);
  // std::cout<< " Predicted Tau Value : " << PredictedTau << "\n";

  
  return PredictedTau;

};



#endif
