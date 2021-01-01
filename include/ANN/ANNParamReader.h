// =======================================================================
// @(#)TANNParamReader.h        1.0 24.12.20
// 
// Class:       TANNParamReader
// Purpose:     Class for reading the parameter file 
//
// Author:      Subodh M. Joshi (28.12.20)
//
// History:     start of implementation 28.12.20 (Subodh M. Joshi)
//
// =======================================================================

#ifndef __TANNParamReader__
#define __TANNParamReader__

#include <ANNIncludes.h>

class TANNParamReader
{
  public:
    /** constructor */
    TANNParamReader(char *paramFile);

    /** destructor */
    ~TANNParamReader();

    /** number of hidden layers in FFN */
    int nHL;

    /** array storing dim of each layer (from IP to OP) */
    int *layerDim;

    /** array storing type of each layer (from IP to OP) */
    std::string *layerType;

  private:
    /** flags */
    bool layerDimFlag;
    bool layerTypeFlag;
    bool layerFlag;
    
  public:
    /** Class methods **/

    /** read param file */
    void readParamFile(char *paramFile);

  private:
    /** Class methods private */
    std::string getLayerType(int number);
};

#endif
