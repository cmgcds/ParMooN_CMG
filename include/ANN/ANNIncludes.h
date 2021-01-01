// =======================================================================
// @(#)ANNIncludes.h        1.0 24.12.20
// 
// Class:       Not a class.
//
// Purpose:     A header file to collect all the #include for ANN classes 
//
// Author:      Subodh M. Joshi (24.12.20)
//
// History:     start of implementation 24.12.20 (Subodh M. Joshi)
//
// =======================================================================

#ifndef __ANNInclude__
#define __ANNInclude__

// Standard routines
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <assert.h>
#include <iomanip>
#include <string>
#include <time.h>

// MLPACK routines
#include <mlpack/core.hpp>
#include <mlpack/methods/ann/layer/layer.hpp>
#include <mlpack/methods/ann/ffn.hpp>

// Parmoon ANN file includes
#include <ANNParamReader.h>
#include <ANNFunctions.h>

#endif
