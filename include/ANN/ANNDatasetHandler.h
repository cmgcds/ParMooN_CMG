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

#ifndef __TANNDATASETHANDLER__
#define __TANNDATASETHANDLER__

#include <ANNIncludes.h>


class TANNDatasetHandler
{
  public:
    /** Constructor */
    TANNDatasetHandler(TANNParamReader *paramReader);

    /** Destructor */
    ~TANNDatasetHandler();

  public:

    int totalNumberOfSamples;

    int numberOfTrainingSamples;

    int numberOfValidationSamples;

    int numberOfTestingSamples;

    int epochs;

    int ipDataDim;

    std::string saveDataFile;

    /** Armadillo matrix storing raw data */
    arma::mat allData;

    /** Armadillo matrix storing training data */
    arma::mat trainData;

    /** Armadillo matrix storing validation data */
    arma::mat validationData;

    /** Armadillo matrix storing testing data */
    arma::mat testData;

    /** Armadillo matrix storing the training lables */
    arma::mat trainLabels;

    /** Armadillo matrix storing the validation lables */
    arma::mat validationLabels;

    /** Armadillo matrix storing the testing lables */
    arma::mat testLabels;

    /** Armadillo matrix storing the temperory (unprocessed)
     * prediction values */
    arma::mat predictionTemp;

    /** Armadillo matrix storing the processed (final)
     * prediction values */
    arma::mat prediction;

    /** training data scaler,i.e. betn [0,1] */
    mlpack::data::MinMaxScaler dataScaler;

    /** Training label scaler ,i.e. betn [0,1]*/
    mlpack::data::MinMaxScaler labelScaler;


    /** error */
    double errorL1Absolute; // L1 norm
    double errorL2Absolute; // L2 norm
    double errorLInfAbsolute; // max error
    double errorL0Absolute; // min error
    double errorMSE; // Mean Sq. error based on mlpack routines

    double errorL1Relative;
    double errorL2Relative;
    double errorLInfRelative;

  public:
    /** Methods of the class */
    void postProcessResults();

    double computeError(arma::mat referenceValue, arma::mat numericalValue, std::string norm, std::string type);

    void writeErrorFile(arma::mat referenceValue, arma::mat numericalValue);

    // save the model 
    // void saveModel();
};

#endif
