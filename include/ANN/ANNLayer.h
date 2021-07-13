// =======================================================================
// @(#)TANN.h        1.0 24.12.20
// 
// Class:       TANN
// Purpose:     Class for storing data for a single (isolated) layer of the ANN
//
// Author:      Subodh M. Joshi (24.12.20)
//
// History:     start of implementation 24.12.20 (Subodh M. Joshi)
//
// =======================================================================

#ifndef __TANN_LAYER__
#define __TANN_LAYER__

#include <ANNIncludes.h>

/** Class for artificial neural networks (ANN) **/
class TANNLayer
{
  public:
  /** constructor 1 */  
  TANNLayer();

  /** constructor 2 */
  TANNLayer(int rankArg, int dimArg, std::string typeArg, int typeIntArg);

  /** destrcutor */
  ~TANNLayer();

  public:
  /** Data members */
  /** Rank in the neural network. 
   *  e.g. 0: Input, 1: First hidden layer ..., N-1: Output 
   *  for a N layered network */
  int rank;    

  /** Size of the layer */
  int dim;

  /** Type of the network 
   * Supported types are: IP, OP, HIDDEN  */
  std::string type;

  /** Int code for the type of the network */
  int typeInt;
};

#endif
