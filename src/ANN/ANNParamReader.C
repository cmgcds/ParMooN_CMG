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
  this->readParamFile(paramFile);
};

TANNParamReader::~TANNParamReader(){
  assert(layerFlag == true and "Error while attempting to call destructor of TANNParamReader defined in src/ANN/ANNParamReader.C");

  delete[] layerDim;
  delete[] layerType;
  delete[] layerTypeInt;

};


void TANNParamReader::readParamFile(char *paramFile){
#ifdef _MPI
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  char line[100];
  int N_Param = 0;
  std::ifstream dat(paramFile);

  if (!dat)
  {
#ifdef _MPI
  if(rank==0)
#endif
    std::cout << "cannot open '" << paramFile << "' for input" << std::endl;
    exit(-1);
  }

  while (!dat.eof())
  {
    dat >> line;

    if (!strcmp(line, "ANN_NHL:"))
    {
      // Store value of the number of hidden layers (nHL)
      dat >> nHL;
      N_Param++;

      // Create an array to store the dimension of each layer (nHL + 2 for IP and OP resp.)
      this->layerDim = new int[nHL+2];
      layerDimFlag = true;

      // Create an array to store the type of the layer (string , "IP", "OP" or "HIDDEN")
      this->layerType = new std::string[nHL+2];

      // Create and array to store integer code for the layer type
      this->layerTypeInt = new int[nHL+2];
      layerTypeFlag = true;



    }

    // Set the dim of the IP data to input layer
    // This could be considered as number of inputs to each input neuron 
    if (!strcmp(line, "ANN_IPDATADIM:"))
    {
      assert(layerDimFlag == true and layerTypeFlag == true and "Error in reading ANN_NHL: from the .dat file");
      // Set input layer dim
      dat >> ipDataDim;
      N_Param++;
    }


    // Set the dim of the IP layer
    if (!strcmp(line, "ANN_IPLDIM:"))
    {
      assert(layerDimFlag == true and layerTypeFlag == true and "Error in reading ANN_NHL: from the .dat file");
      // Set input layer dim
      dat >> layerDim[0];
      N_Param++;
    }

    // Set the type of the IP layer
    if (!strcmp(line, "ANN_IPLTYPE:"))
    {
      assert(layerDimFlag == true and layerTypeFlag == true and "Error in reading ANN_NHL: from the .dat file");
      // Set input layer type
      int layerTypeNumber;
      dat >> layerTypeNumber;
      layerType[0] = getLayerType(layerTypeNumber);
      layerTypeInt[0] = layerTypeNumber;
      N_Param++;
    }

    // Set the dim of the OP layer
    if (!strcmp(line, "ANN_OPLDIM:"))
    {
      assert(layerDimFlag == true and layerTypeFlag == true and "Error in reading ANN_NHL: from the .dat file");
      dat >> layerDim[nHL+2-1];
      N_Param++;
    }

    // Set the dim of the OP layer
    if (!strcmp(line, "ANN_OPLTYPE:"))
    {
      assert(layerDimFlag == true and layerTypeFlag == true and "Error in reading ANN_NHL: from the .dat file");
      int layerTypeNumber;
      dat >> layerTypeNumber;
      layerType[nHL+2-1] = getLayerType(layerTypeNumber);
      layerTypeInt[nHL+2-1] = layerTypeNumber;
      N_Param++;
    }

    // Set the dim of the hidden layers
    for (int i=0; i<nHL; i++){
      std::string name = "ANN_HL_"+varToString(i)+"_DIM:";
      const char *iChar = name.c_str();

      if (!strcmp(line,iChar)){
        assert(layerDimFlag == true and layerTypeFlag == true and "Error in reading ANN_NHL: from the .dat file");
        dat >> layerDim[1+i];
        N_Param++;
      };
    };
 
    // Set the type of the hidden layers
    for (int i=0; i<nHL; i++){
      std::string name = "ANN_HL_"+varToString(i)+"_TYPE:";
      const char *iChar = name.c_str();

      if (!strcmp(line,iChar)){
        assert(layerDimFlag == true and layerTypeFlag == true and "Error in reading ANN_NHL: from the .dat file");
        int layerTypeNumber;
        dat >> layerTypeNumber;
        layerType[i+1] = getLayerType(layerTypeNumber);
        layerTypeInt[i+1] = layerTypeNumber;
        N_Param++;
      };
    };
 
   // read until end of line
    dat.getline (line, 99);

   // Set the layerFlag to true
    this->layerFlag = true;
  }
  dat.close();

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
