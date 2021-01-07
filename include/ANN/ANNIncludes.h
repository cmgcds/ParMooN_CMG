// =======================================================================
// @(#)ANNIncludes.h        1.0 24.12.20
// 
// Class:       Not a class.
//
// Purpose:     A header file to collect all the #include for ANN classes 
//
// Author:      Subodh M. Joshi (07.01.21)
//
// History:     start of implementation 24.12.20 (Subodh M. Joshi)
//
// =======================================================================

#ifndef __ANNINCLUDES__
#define __ANNINCLUDES__

// Standard routines
#include <iostream>
#include <sstream>
#include <fstream>
#include <assert.h>
#include <string>
#include <unistd.h>
#include <cstring>

// MLPACK routines
#include <mlpack/core.hpp>
#include <mlpack/methods/ann/layer/layer.hpp>

// MLPACK Loss functions
#include <mlpack/methods/ann/loss_functions/mean_squared_error.hpp>
#include <mlpack/methods/ann/loss_functions/cross_entropy_error.hpp>
#include <mlpack/methods/ann/loss_functions/mean_squared_logarithmic_error.hpp>
#include <mlpack/methods/ann/loss_functions/sigmoid_cross_entropy_error.hpp>
#include <mlpack/methods/ann/loss_functions/mean_absolute_percentage_error.hpp>

// MLPACK Initialization rules
#include <mlpack/methods/ann/init_rules/gaussian_init.hpp>
#include <mlpack/methods/ann/init_rules/random_init.hpp>
#include <mlpack/methods/ann/init_rules/const_init.hpp>

// MLPACK main routines
// Feed Forward Network
#include <mlpack/methods/ann/ffn.hpp>
// Recurrent Neural Network
#include <mlpack/methods/ann/rnn.hpp>


// Parmoon ANN file includes
#include <ANNFunctions.h>
#include <ANNParamReader.h>


/**************************/
// System definitions
/**************************/

// Loss Functions
#define NEGATIVE_LOG_LIKELIHOOD mlpack::ann::NegativeLogLikelihood<>
#define MEAN_SQUARED_ERROR mlpack::ann::MeanSquaredError<>
#define MEAN_SQUARED_LOGARITHMIC_ERROR mlpack::ann::MeanSquaredLogarithmicError<>
#define CROSS_ENTROPY_ERROR mlpack::ann::CrossEntropyError<>
#define SIGMOID_CROSS_ENTROPY_ERROR mlpack::ann::SigmoidCrossEntropyError<>
#define MEAN_ABSOLUTE_PERCENTAGE_ERROR mlpack::ann::MeanAbsolutePercentageError<>

// Initialization functions
#define RANDOM_INITIALIZATION mlpack::ann::RandomInitialization
#define GAUSSIAN_INITIALIZATION mlpack::ann::GaussianInitialization
#define CONST_INITIALIZATION mlpack::ann::ConstInitialization

#endif
