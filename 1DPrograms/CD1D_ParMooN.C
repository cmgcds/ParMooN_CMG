// =======================================================================
// Purpose:     main program for 1D with new kernels of ParMooN
// Author:      Sashikumaar Ganesan
//
// History:     Implementation started on 12.12.2020
// =======================================================================

#include <Domain.h>
#include <Database.h>
#include <SystemCD1D.h>
#include <FEDatabase2D.h>

// ANN INLUDES
#include <ANN.h>
#include <ANNDatasetHandler.h>

//

// =======================================================================
// include current example
// =======================================================================

#include "../Examples/CD_1D/cd1d.h"

// main program

int main(int argc, char *argv[])
{
	TDatabase *Database = new TDatabase();
	TSystemCD1D *SystemCD;
	TFEDatabase2D *FEDatabase = new TFEDatabase2D();

	std::ostringstream os;
	os << " ";

	// Read the Param Values with Empty constructor
	// The actual TDomain Constructor is constructed inside the System1D
	TDomain *Domain = new TDomain();
	Domain->ReadParam(argv[1]);

	// Obtain FE Mesh Details from the Param File
	int N_Elements = TDatabase::ParamDB->N_ELEMENTS_1D;
	double Start_X = TDatabase::ParamDB->START_X;
	double End_X = TDatabase::ParamDB->END_X;

	int FE_Order = TDatabase::ParamDB->ANSATZ_ORDER;
	int N_DOF = N_Elements * (FE_Order) + 1;

	// cout << " NELEM : " << N_Elements<<endl;
	// cout << " Start : " << Start_X<<endl;
	// cout << " End : " << End_X<<endl;
	// cout << " FE ORDER : " << FE_Order <<endl;
	// cout << " N_DOF : " << N_DOF <<endl;

	TDatabase::ParamDB->USE_ANN_PREDICTED_TAU = 0;

	//======================================================================
	// CD1D System construction
	//======================================================================
	SystemCD = new TSystemCD1D(N_Elements, Start_X, End_X, BoundCondition_LminLMax, BoundVales, argv[1]);


	// initiaize system
	SystemCD->Init(BilinearCoeffs);

	//   SystemCD->Interpolate(Exact);

	TDatabase::ParamDB->USE_ANN_PREDICTED_TAU = 0;
	// interpolate (if needed)
	SystemCD->Solve();

	// Get the Mean Value of Derivatives
	SystemCD->getMeanValueDerivatives();

	// // Print the Solution to terminal
	// SystemCD->printSolution();

	// // Generate VTK File for the problem
	// SystemCD->generateVTK();

	// Generate plot.dat file to be used for using in Tecplot or GNUplot.
	// SystemCD->plotGNU();
	double *Solution = SystemCD->getSolutionVector();
	std::ofstream ofile;
	std::string name = std::to_string(long (TDatabase::ParamDB->PE_NR));
	std::string file = "supgSol_" + name + ".csv";
	ofile.open(file.c_str());
	for (int i = 0; i < N_DOF; i++)
		ofile << Solution[i] << ",";
	ofile << "\n";

	CloseFiles();

	// ------------------------------------------- CREATE A NEW SYSTEM FOR THE ANN MATRIX: ----------------------------------------------- //
	TDatabase::ParamDB->USE_ANN_PREDICTED_TAU =1;
	//-- ANN PARAMETER INCLUSION --- //
	//   // Create a new parameter reader (takes argument .dat file)
	TANNParamReader paramReader(argv[1]);
	// paramReader.print();

	// Create a new dataset handler (to create train and test datasets and labels)
	TANNDatasetHandler datasetHandler(&paramReader);

	// Create a new ANN model (template arguments can be found in ./include/ANN/ANNIncludes.h)
	TANN<MEAN_SQUARED_ERROR, RANDOM_INITIALIZATION> ann(&paramReader);

	char* filenameModel = argv[2];
	// Load External Model
	ann.loadExternalModel(filenameModel);
	// TDatabase::ParamDB->PE_NR = std::stod(argv[2]);z
	double eps = 1.0 / TDatabase::ParamDB->PE_NR;
	cout << " EPS VALUE : " << eps << endl;
	double b = 1.0;
	double h = (End_X - Start_X) / N_Elements;
	double predictedTau = ann.predictTau(b, eps, h, &datasetHandler);
	// TDatabase::ParamDB->GLOBAL_TAU = predictedTau;
	
	TSystemCD1D *SystemCD_ann = new TSystemCD1D(N_Elements, Start_X, End_X, BoundCondition_LminLMax, BoundVales, argv[1]);

	// initiaize system
	SystemCD_ann->Init(BilinearCoeffs);

	SystemCD_ann->InitialiseANNParameters(&datasetHandler,&paramReader,&ann);
	
	// interpolate (if needed)
	SystemCD_ann->Solve();

	// Get the Mean Value of Derivatives
	// SystemCD_ann->getMeanValueDerivatives();

	double Pe_nr = (fabs(b) * h) / (2 * eps);

	double analyticalTau = h / (2 * fabs(b)) * (1.0 / tanh(Pe_nr) - 1.0 / Pe_nr);
	// Generate VTK File for the problem
	SystemCD_ann->generateVTK();

	cout << " =======================================================================" << endl;
	cout << "   PREDCITED TAU  : " << predictedTau << endl;
	cout << "   ANALYTICAL TAU : " << analyticalTau << endl;
	cout << " =======================================================================" << endl;

	double *Solution1 = SystemCD_ann->getSolutionVector();
	for (int i = 0; i < N_DOF; i++)
		ofile << Solution1[i] << ",";
	ofile << "\n";
	


	//------ GAlerkin Solution ----- //
	TSystemCD1D* SystemCD_Galerkin = new TSystemCD1D(N_Elements, Start_X, End_X, BoundCondition_LminLMax, BoundVales, argv[1]);


	// initiaize system
	SystemCD_Galerkin->Init(BilinearCoeffs);

	//   SystemCD_Galerkin->Interpolate(Exact);

	TDatabase::ParamDB->USE_ANN_PREDICTED_TAU = 0;
	TDatabase::ParamDB->DISCTYPE = 1;
	// interpolate (if needed)
	SystemCD_Galerkin->Solve();

	// Get the Mean Value of Derivatives
	SystemCD_Galerkin->getMeanValueDerivatives();

	// // Print the Solution to terminal
	// SystemCD_Galerkin->printSolution();


	// Generate plot.dat file to be used for using in Tecplot or GNUplot.
	// SystemCD_Galerkin->plotGNU();
	double *Solution_Galerkin = SystemCD_Galerkin->getSolutionVector();

	for (int i = 0; i < N_DOF; i++)
		ofile << Solution_Galerkin[i] << ",";
	ofile<<"\n";
	ofile.close();


	// Write out the tau parameters 
	std::ofstream ofile_tau;
	std::string name_tau = std::to_string(long (TDatabase::ParamDB->PE_NR));
	std::string file_tau = "supgTau_" + name_tau + ".csv";
	
	ofile_tau.open(file_tau.c_str());

	ofile_tau << predictedTau << ","<<analyticalTau<<endl;
	
	ofile_tau.close();

	return 0;
} // end main
