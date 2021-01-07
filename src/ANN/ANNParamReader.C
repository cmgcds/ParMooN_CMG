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

  // Create an array to store the dimension of each layer (nHL + 2 for IP and OP resp.)
  this->layerDim = new int[nHL + 2];
  layerDimFlag = true;

  // Create an array to store the type of the layer (string , "IP", "OP" or "HIDDEN")
  this->layerType = new std::string[nHL+2];

  // Create and array to store integer code for the layer type
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
      }


      if (word1 == "ANN_IPLTYPE"){
        assert(layerDimFlag == true and layerTypeFlag == true and "Error in reading ANN_NHL: from the .dat file");
        // Set input layer dim
        numberInt = stringToInt(line);
        layerTypeInt[0] = numberInt;
        layerType[0] = getLayerType(numberInt);
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
          layerDim[1+i] = numberInt;
        };
      };


      // Set the type of the hidden layers
      for (int i=0; i<nHL; i++){
        std::string name = "ANN_HL_"+varToString(i)+"_TYPE";
        if (word1 == name){
          assert(layerDimFlag == true and layerTypeFlag == true and "Error in reading ANN_NHL: from the .dat file");
          numberInt = stringToInt(line);
          layerTypeInt[1+i] = numberInt;
          layerType[i+1] = getLayerType(numberInt);
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


      if (word1 == "ANN_EPOCHS"){
        assert(layerDimFlag == true and layerTypeFlag == true and "Error in reading ANN_NHL: from the .dat file");
        // Set input layer dim
        numberInt = stringToInt(line);
        this->epochs = numberInt;
      }
    };
  }
  else{
    std::cout << " File " + varToString(paramFile) + " can not be opened" << std::endl;
  };
  file.close();
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
