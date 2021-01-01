/** ************************************************************************ 
* @brief     source file for TANNLayer
* @author    Subodh M. Joshi
* @date      24.12.20
* @History   Implementation started 24.12.20 (Subodh M. Joshi) 
 ************************************************************************  */

/** Class declaration is done in include/ANN/ANNLayer.h */

#include <ANNLayer.h>

/** Constructor 1 */
TANNLayer::TANNLayer(){};

/** Constructor 2 */
TANNLayer::TANNLayer(int rankArg, int dimArg, std::string typeArg){

  // Set rank
  this->rank = rankArg;

  // Set dimension (i.e. size of the layer)
  this->dim = dimArg;

  // Set the type of the layer
  this->type = typeArg;

};

TANNLayer::~TANNLayer(){};

