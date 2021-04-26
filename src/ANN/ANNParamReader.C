/** ************************************************************************ 
* @brief     source file for TANNParamReader
* @author    Subodh M. Joshi
* @date      28.12.20
* @History   Implementation started 28.12.20 (Subodh M. Joshi) 
 ************************************************************************  */

/** Class declaration is done in include/ANN/ANNParamReader.h */

#include <ANNParamReader.h>

TANNParamReader::TANNParamReader(char *paramFile){
  /** Initiate all flags with zero */
  this->layerDimFlag = false;
  this->layerTypeFlag = false;
  this->layerFlag = false;

  /** Read the param file */
  this->readNumberOfLayers(paramFile);
  assert(layerFlag == true);

  // Create an array to store the type of the layer (i.e. activation function)
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
   **/
  // Create an array to store the dimension of each layer (nHL + 2 for IP and OP resp.)
  this->layerDim = new int[nHL + 2];
  layerDimFlag = true;

  // Create an array to store activation functions for each layer 
  // NOTE: a dummy activation function is assigned at the input layer
  this->layerType = new std::string[nHL+2];

  // Create and array to store integer code for the layer type
  // NOTE: a dummy activation function type is assigned at the input layer
  this->layerTypeInt = new int[nHL+2];
  layerTypeFlag = true;


  // ___________________________________________
  // Set the default values 
  // ___________________________________________
    
  // Set the default value of the validation set to 7
  this->validationDataPercentage = 10;

  // Set the default value for optimizerCode to 0
  this->optimizerCode = 0;

  // Set the default value for the optimizer step size to 0.01
  this->optimizerStepSize = 0.001;

  // Set the default value of the batch size to 32
  this->sgdBatchSize = 32;

  // Set the default value for epochs to 1
  this->epochs = 1;

  // Set the default value for the maximum iterations of optimizer to 100000
  this->maxIterations = 100000;

  // Pre-set the feature Scaling constant to 1.0
  this->featureScalingConstant = 1.0;

  // Set the default value for the dropout ratio to 0.2
  this->dropoutRatio = 0.2;

  // Set the default value for the tolerance to 1e-5
  this->tolerance = 1e-5;

  // Set the default value for the save model file
  this->saveModelFile = "FFN.bin";

  // Set the default value for the load model file
  this->loadModelFile = "FFN.bin";

  // Set the default value for the save data file
  // The results for the testNetwork are stored in  this
  this->saveDataFile = "Results.csv";

  // Set the default flag for the load model file to zero
  this->loadModelFlag = 0;

  //_________________________________________
  //Read the param file
  this->readParamFile(paramFile);
};

TANNParamReader::~TANNParamReader(){
  assert(layerFlag == true and "Error while attempting to call destructor of TANNParamReader defined in src/ANN/ANNParamReader.C");

  delete[] layerDim;
  delete[] layerType;
  delete[] layerTypeInt;

};


void TANNParamReader::readNumberOfLayers(char *paramFile){

  std::ifstream file(paramFile);

  std::string line;
  std::string word1;
  std::string delimiter;
  int numberInt;
  double numberDouble;

  delimiter = ":";


  if (file.is_open()){
    while(std::getline(file, line)){
      word1 = line.substr(0,line.find(delimiter));
      line = line.erase(0,word1.length()+1);
      word1.erase(word1.find_last_not_of(" \r\t")+1);
      if(word1 == "ANN_NHL"){
        numberInt = stringToInt(line);
        this->nHL = numberInt;
      }
    };
  };

  // Set the layerFlag to true
  this->layerFlag = true;
  file.close();
};


void TANNParamReader::readParamFile(char *paramFile){
  assert (layerFlag == true);

  std::ifstream file(paramFile);

  std::string line;
  std::string word1;
  std::string delimiter;
  int numberInt;
  double numberDouble;

  delimiter = ":";


  if (file.is_open()){
    while(std::getline(file, line)){
      word1 = line.substr(0,line.find(delimiter));
      line = line.erase(0,word1.length()+1);
      word1.erase(word1.find_last_not_of(" \r\t")+1);

      if (word1 == "ANN_IPDATADIM"){
        assert(layerDimFlag == true and layerTypeFlag == true and "Error in reading ANN_NHL: from the .dat file");
        // Set input layer dim
        numberInt = stringToInt(line);
        ipDataDim = numberInt;
      }
      

      if (word1 == "ANN_IPLDIM"){
        assert(layerDimFlag == true and layerTypeFlag == true and "Error in reading ANN_NHL: from the .dat file");
        // Set input layer dim
        numberInt = stringToInt(line);
        layerDim[0] = numberInt;
        // Set a dummy activation function for the input layer
        layerTypeInt[0] = 100;
        // Set a dummy activation function for the input layer
        layerType[0] = getLayerType(100);
      }


      if (word1 == "ANN_OPLDIM"){
        assert(layerDimFlag == true and layerTypeFlag == true and "Error in reading ANN_NHL: from the .dat file");
        // Set input layer dim
        numberInt = stringToInt(line);
        layerDim[nHL+2-1] = numberInt;
      }


      if (word1 == "ANN_OPLTYPE"){
        assert(layerDimFlag == true and layerTypeFlag == true and "Error in reading ANN_NHL: from the .dat file");
        // Set input layer dim
        numberInt = stringToInt(line);
        layerTypeInt[nHL+2-1] = numberInt;
        layerType[nHL+2-1] = getLayerType(numberInt);
      }


      // Set the dim of the hidden layers
      for (int i=0; i<nHL; i++){
        std::string name = "ANN_HL_"+varToString(i)+"_DIM";
        if (word1 == name){
          assert(layerDimFlag == true and layerTypeFlag == true and "Error in reading ANN_NHL: from the .dat file");
          numberInt = stringToInt(line);
          // Remember, the zeroth layer is IP layer
          layerDim[1+i] = numberInt;
        };
      };


      // Set the type of the hidden layers
      for (int i=0; i<nHL; i++){
        std::string name = "ANN_HL_"+varToString(i)+"_TYPE";
        if (word1 == name){
          assert(layerDimFlag == true and layerTypeFlag == true and "Error in reading ANN_NHL: from the .dat file");
          numberInt = stringToInt(line);
          // Remember, the zeroth layer is IP layer
          layerTypeInt[1+i] = numberInt;
          layerType[1+i] = getLayerType(numberInt);
        };
      };


      if (word1 == "ANN_TRAINING_DATASET_NAME"){
        assert(layerDimFlag == true and layerTypeFlag == true and "Error in reading ANN_NHL: from the .dat file");
        // Set input layer dim
        line.erase(0,line.find_first_of(" \r\t")+1);
        this->trainingDatasetName = line;
      }

      if (word1 == "ANN_TESTING_DATASET_NAME"){
        assert(layerDimFlag == true and layerTypeFlag == true and "Error in reading ANN_NHL: from the .dat file");
        // Set input layer dim
        line.erase(0,line.find_first_of(" \r\t")+1);
        this->testingDatasetName = line;
      }

      // Read the percentage of the training date out of the total data
      if (word1 == "ANN_VALIDATION_DATA_PERCENTAGE")
      {
        assert(layerDimFlag == true and layerTypeFlag == true and "Error in reading ANN_NHL: from the .dat file");
        numberDouble = stringToDouble(line);
        validationDataPercentage = numberDouble;
      }

      if (word1 == "ANN_OPTIMIZER_CODE"){
        assert(layerDimFlag == true and layerTypeFlag == true and "Error in reading ANN_NHL: from the .dat file");
        // Optimizer codes:
        // 0: default (RMSProp)
        // 1: Gradient Descent
        // 2: SGD
        // 3: Adam
        // 4: L-BFGS
        numberInt = stringToInt(line);
        this->optimizerCode = numberInt;
      }

      // Read the optimizer step size
      if (word1 == "ANN_OPTIMIZER_STEP_SIZE")
      {
        assert(layerDimFlag == true and layerTypeFlag == true and "Error in reading ANN_NHL: from the .dat file");
        numberDouble = stringToDouble(line);
        optimizerStepSize = numberDouble;
      }

      if (word1 == "ANN_SGD_BATCH_SIZE"){
        assert(layerDimFlag == true and layerTypeFlag == true and "Error in reading ANN_NHL: from the .dat file");
        // Set the batch size for the stochastic gradient descent method
        numberInt = stringToInt(line);
        this->sgdBatchSize = numberInt;
      }

      if (word1 == "ANN_MAX_ITERATIONS"){
        assert(layerDimFlag == true and layerTypeFlag == true and "Error in reading ANN_NHL: from the .dat file");
        numberInt = stringToInt(line);
        this->maxIterations = numberInt;
      }

      if (word1 == "ANN_EPOCHS"){
        assert(layerDimFlag == true and layerTypeFlag == true and "Error in reading ANN_NHL: from the .dat file");
        // Set input layer dim
        numberInt = stringToInt(line);
        this->epochs = numberInt;
      }

      // Read the parameter used for feature scaling
      if (word1 == "ANN_DROPOUT_RATIO")
      {
        assert(layerDimFlag == true and layerTypeFlag == true and "Error in reading ANN_NHL: from the .dat file");
        numberDouble = stringToDouble(line);
        dropoutRatio = numberDouble;
      }

      // Read the parameter used for feature scaling
      if (word1 == "ANN_TOLERANCE")
      {
        assert(layerDimFlag == true and layerTypeFlag == true and "Error in reading ANN_NHL: from the .dat file");
        numberDouble = stringToDouble(line);
        tolerance = numberDouble;
      }

      // Read the parameter used for feature scaling
      if (word1 == "ANN_FEATURE_SCALING_CONSTANT")
      {
        assert(layerDimFlag == true and layerTypeFlag == true and "Error in reading ANN_NHL: from the .dat file");
        numberDouble = stringToDouble(line);
        featureScalingConstant = numberDouble;
      }

      // Read the file name for storing the trained model after the training is over
      if (word1 == "ANN_SAVE_MODEL_FILE"){
        assert(layerDimFlag == true and layerTypeFlag == true and "Error in reading ANN_NHL: from the .dat file");
        this->saveModelFile = line;
      }

      // Read the file name for storing the results data for the tests
      if (word1 == "ANN_SAVE_DATA_FILE"){
        assert(layerDimFlag == true and layerTypeFlag == true and "Error in reading ANN_NHL: from the .dat file");
        this->saveDataFile = line;
      }

      // Read the file name for reading the trained model from that file
      if (word1 == "ANN_LOAD_MODEL_FILE"){
        assert(layerDimFlag == true and layerTypeFlag == true and "Error in reading ANN_NHL: from the .dat file");
        this->loadModelFile = line;
      }

      // Flag for training new vs. loading existing model
      if (word1 == "ANN_LOAD_FLAG"){
        assert(layerDimFlag == true and layerTypeFlag == true and "Error in reading ANN_NHL: from the .dat file");
        // Set input layer dim
        numberInt = stringToInt(line);
        this->loadModelFlag = numberInt;
      }

    };
  }
  else{
    std::cerr << " File " + varToString(paramFile) + " can not be opened" << std::endl;
  };
  file.close();
};


void TANNParamReader::print(){
  std::cout << "___________________________________________________\n";
  std::cout << "___________________________________________________\n";
  
  std:: cout << " Param Data: " << std::endl;

  std::cout << " Training Dataset Name: " << this->trainingDatasetName << std::endl;

  std::cout << " Testing Dataset Name: " << this->testingDatasetName << std::endl;

  std::cout << " Layers info: " << std::endl;
  std::cout << "___________________________\n";
  for (int i=0; i< this->nHL + 2 ; i++){
    std::cout << " Layer Number: " << i << "   Layer Dim: " << this->layerDim[i] << "    Layer Type: " << this->layerType[i] << "   Layer Type code: " << this->layerTypeInt[i] << std::endl;
  };
  std::cout << "___________________________\n";

  std::cout << " Training Data Percentage   " << 100 - this->validationDataPercentage << "    Default value: 90 " << std::endl;
  std::cout << " Validation Data Percentage   " << this->validationDataPercentage  << "    Default value: 10" << std::endl;
  std::cout << " Optimizer Code:  " <<  this->optimizerCode << "   NOTE 0: default (RMSProp), 1:GD, 2:SGD, 3: Adam, 4:L-BFGS " << std::endl;
  std::cout << " Optimizer step size " << this->optimizerStepSize  << "    Default value: 0.001" <<  std::endl;
  std::cout << " SGD batch size " << this->sgdBatchSize  << "    Default value: 32" << std::endl;
  std::cout << " Epochs " << this->epochs  << "    Default value: 1" << std::endl;
  std::cout << " Max Iterations " << this->maxIterations  << "    Default value: 100000" << std::endl;
  std::cout << " Feature Scaling Parameter " << this->featureScalingConstant  << "    Default value: 1.0 " << std::endl;
  std::cout << " Dropout Ratio " << this->dropoutRatio  << "    Default value: 0.2 " << std::endl;
  std::cout << " Tolerance " << this->tolerance << "    Default value: 1e-5 " << std::endl;
  std::cout << "___________________________________________________\n";
};

std::string TANNParamReader::getLayerType(int number){
  switch (number){
    case 0:
      return "SigmoidLayer";
      break;
    
    case 1:
      return "IdentityLayer";
      break;
    
    case 2:
      return "ReLULayer";
      break;

    case 3:
      return "LeakyReLU";
      break;

    case 4:
      return "TanHLayer";
      break;

    case 5:
      return "SoftplusLayer";
      break;

    case 6:
      return "LogSoftMax";

    case 100:
      return "Dummy";
      break;

    default:
      return "SigmoidLayer";
      break;

  };
};

std::string TANNParamReader::getLossFunction(int number){
  switch(number){
    case 0:
      // Mean squared error
      return "MeanSquaredError";
      break;

    case 1:
      // Cross entropy error
      return "CrossEntropyError";
      break;

    default:
      return "MeanSquaredError";
      break;
  };
};
