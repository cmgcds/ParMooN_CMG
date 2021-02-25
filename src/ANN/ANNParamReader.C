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


      if (word1 == "ANN_DATASET_NAME"){
        assert(layerDimFlag == true and layerTypeFlag == true and "Error in reading ANN_NHL: from the .dat file");
        // Set input layer dim
        line.erase(0,line.find_first_of(" \r\t")+1);
        this->datasetName = line;
      }


      // Read the percentage of the training date out of the total data
      if (word1 == "ANN_TRAINING_DATA_PERCENTAGE")
      {
        assert(layerDimFlag == true and layerTypeFlag == true and "Error in reading ANN_NHL: from the .dat file");
        numberDouble = stringToDouble(line);
        trainingDataPercentage = numberDouble;
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
        // 1: SGD
        // 2: Adam
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
      if (word1 == "ANN_FEATURE_SCALING_CONSTANT")
      {
        assert(layerDimFlag == true and layerTypeFlag == true and "Error in reading ANN_NHL: from the .dat file");
        numberDouble = stringToDouble(line);
        featureScalingConstant = numberDouble;
      }

    };
  }
  else{
    std::cout << " File " + varToString(paramFile) + " can not be opened" << std::endl;
  };
  file.close();
};


void TANNParamReader::print(){
  std::cout << "___________________________________________________\n";
  std::cout << "___________________________________________________\n";
  
  std:: cout << " Param Data: " << std::endl;

  std::cout << " Dataset Name: " << this->datasetName << std::endl;

  std::cout << " Layers info: " << std::endl;
  std::cout << "___________________________\n";
  for (int i=0; i< this->nHL + 2 ; i++){
    std::cout << " Layer Number: " << i << "   Layer Dim: " << this->layerDim[i] << "    Layer Type: " << this->layerType[i] << "   Layer Type code: " << this->layerTypeInt[i] << std::endl;
  };
  std::cout << "___________________________\n";

  std::cout << " Training Data Percentage   " << this->trainingDataPercentage  << std::endl;
  std::cout << " Validation Data Percentage   " << this->validationDataPercentage  << std::endl;
  std::cout << " Optimizer Code:  " <<  this->optimizerCode << "   NOTE 0: default (RMSProp), 1:SGD, 2: Adam " << std::endl;
  std::cout << " Optimizer step size " << this->optimizerStepSize  << std::endl;
  std::cout << " SGD batch size " << this->sgdBatchSize  << std::endl;
  std::cout << " Epochs " << this->epochs  << std::endl;
  std::cout << " Max Iterations " << this->maxIterations  << std::endl;
  std::cout << " Feature Scaling Parameter " << this->featureScalingConstant  << std::endl;
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
      return "TanHLayer";
      break;

    case 4:
      return "SoftplusLayer";
      break;

    case 5:
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
