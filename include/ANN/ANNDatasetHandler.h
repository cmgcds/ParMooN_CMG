// =======================================================================
// @(#)TANNDatasetHandler.h        1.0 24.12.20
// 
// Class:       TANNDatasetHandler
// Purpose:     Class for handling datasets for Artificial Neural Networks (ANN)
//
// Author:      Subodh M. Joshi (25.12.20)
//
// History:     start of implementation 25.12.20 (Subodh M. Joshi)
//
// =======================================================================

#ifndef __TANNDatasetHandler__
#define __TANNDatasetHandler__

#include <ANNIncludes.h>

class ANNDadasetHandler
{
  public:
    /** Constructor */
    ANNDatasetHandler(string nameArg);

    /** Destructor */
    ~ANNDatasetHandler();

    /** Armadillo matrix storing training data */
    arma::mat trainData;

    /** Armadillo matrix storing testing data */
    arma::mat testData;

    /** Armadillo matrix storing the training lables */
    arma::mat trainLables;

    /** Armadillo matrix storing the testing lables */
    arma::mat testLables;

    /** Armadillo matrix storing the temperory (unprocessed)
     * prediction values */
    arma::mat predictionTemp;

    /** Armadillo matrix storing the processed (final)
     * prediction values */
    arma::mat prediction;

};

#endif
