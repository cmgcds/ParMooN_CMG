// =======================================================================
//
// Purpose:     main program for scalar equations with new kernels of ParMooN
//
// Author:      Sashikumaar Ganesan
//
// History:     Implementation started on 08.08.2014

// =======================================================================

#include <Domain.h>
#include <Database.h>
#include <SystemTCD2D.h>
#include <SystemTCD2D_ALE.h>
#include <FEDatabase2D.h>
#include <FESpace2D.h>
#include <SquareStructure2D.h>
#include <Structure2D.h>
#include <QuadAffin.h>
#include <DirectSolver.h>
#include <Assemble2D.h>
#include <Output2D.h>
#include <LinAlg.h>
#include <CD2DErrorEstimator.h>
#include <MainUtilities.h>
#include <TimeDiscRout.h>

#include <string.h>
#include <sstream>
#include <MooNMD_Io.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>

// =======================================================================
// include current example
// =======================================================================
#include "../Main_Users/Thivin/SPADE/advection.h"
// #include "../Examples/TCD_2D/SinCos1.h"
// #include "../Examples_All/TCD_2D/Time3.h"
// #include "../Examples/TCD_2D/exp_0.h"
//    #include "../Examples/TCD_2D/exp_2.h"
// #include "../Examples_All/TCD_2D/exp_1.h"
// #include "../Main_Users/Sashi/TCD_2D/Hemker.h"

int main(int argc, char *argv[])
{
	int i, j, l, m, N_SubSteps, ORDER, N_Cells, N_DOF, img = 1, N_G;
	int N_Active;

	double *sol, *rhs, *oldrhs, t1, t2, errors[5], Linfty;
	double tau, end_time, *defect, olderror, olderror1, hmin, hmax;

	bool UpdateStiffnessMat, UpdateRhs, ConvectionFirstTime;
	char *VtkBaseName;
	const char vtkdir[] = "VTK";

	TDomain *Domain;
	TDatabase *Database = new TDatabase();
	TFEDatabase2D *FEDatabase = new TFEDatabase2D();
	TCollection *coll;
	TFESpace2D *Scalar_FeSpace, *fesp[1];
	TFEFunction2D *Scalar_FeFunction_Mean;
	TOutput2D *Output;
	TSystemTCD2D *SystemMatrix_Mean;
    TSystemTCD2D *SystemMatrix_Mode;
    TSystemTCD2D *SystemMatrix_Coefficient;
	TAuxParam2D *aux;
	MultiIndex2D AllDerivatives[3] = {D00, D10, D01};

	std::ostringstream os;
	os << " ";

	// ======================================================================
	// set the database values and generate mesh
	// ======================================================================
	// set variables' value in TDatabase using argv[1] (*.dat file), and generate the MESH based
	Domain = new TDomain(argv[1]);

	if (TDatabase::ParamDB->PROBLEM_TYPE == 0)
		TDatabase::ParamDB->PROBLEM_TYPE = 2;
	OpenFiles();

	Database->WriteParamDB(argv[0]);
	Database->WriteTimeDB();
	ExampleFile();

	/* include the mesh from a mesh generator, for a standard mesh use the 
   * build-in function. The GEOFILE describes the boundary of the domain. */
	if (TDatabase::ParamDB->MESH_TYPE == 0)
	{
		Domain->ReadGeo(TDatabase::ParamDB->GEOFILE);
		OutPut("PRM-GEO used for meshing !!!" << endl);
	} // ParMooN  build-in Geo mesh
	else if (TDatabase::ParamDB->MESH_TYPE == 1)
	{
		Domain->GmshGen(TDatabase::ParamDB->GEOFILE);
		OutPut("GMSH used for meshing !!!" << endl);
	}											 //gmsh mesh
	else if (TDatabase::ParamDB->MESH_TYPE == 2) //triangle mesh
	{
		OutPut("Triangle.h used for meshing !!!" << endl);
		// TriaReMeshGen(Domain);
	}
	else
	{
		OutPut("Mesh Type not known, set MESH_TYPE correctly!!!" << endl);
		exit(0);
	}

	// refine grid up to the coarsest level
	for (i = 0; i < TDatabase::ParamDB->UNIFORM_STEPS; i++)
		Domain->RegRefineAll();

	// write grid into an Postscript file
	if (TDatabase::ParamDB->WRITE_PS)
		Domain->PS("Domain.ps", It_Finest, 0);

	// create output directory, if not already existing
	if (TDatabase::ParamDB->WRITE_VTK)
		mkdir(vtkdir, 0777);

	//=========================================================================
	// construct all finite element spaces
	//=========================================================================
	ORDER = TDatabase::ParamDB->ANSATZ_ORDER;

	coll = Domain->GetCollection(It_Finest, 0);
	N_Cells = coll->GetN_Cells();
	OutPut("N_Cells (space) : " << N_Cells << endl);

	// fespaces for scalar equation
	Scalar_FeSpace = new TFESpace2D(coll, (char *)"fe space", (char *)"solution space",
									BoundCondition, 1, NULL);


	N_DOF = Scalar_FeSpace->GetN_DegreesOfFreedom();
	N_Active = Scalar_FeSpace->GetActiveBound();
	OutPut("dof all      : " << setw(10) << N_DOF << endl);
	OutPut("dof active   : " << setw(10) << N_Active << endl);

	//======================================================================
	// construct all finite element functions
	//======================================================================
	sol = new double[N_DOF];
	rhs = new double[N_DOF];
	oldrhs = new double[N_DOF];

	memset(sol, 0, N_DOF * SizeOfDouble);
	memset(rhs, 0, N_DOF * SizeOfDouble);

	Scalar_FeFunction_Mean = new TFEFunction2D(Scalar_FeSpace, (char *)"sol", (char *)"sol", sol, N_DOF);

	//interpolate the initial value
	Scalar_FeFunction_Mean->Interpolate(InitialCondition);

	//======================================================================
	// SystemMatrix construction and solution
	//======================================================================
	// Disc type: GALERKIN (or) SDFEM  (or) UPWIND (or) SUPG (or) LOCAL_PROJECTION
	// Solver: AMG_SOLVE (or) GMG  (or) DIRECT
	SystemMatrix_Mean = new TSystemTCD2D(Scalar_FeSpace, GALERKIN, DIRECT);

	// initilize the system matrix with the functions defined in Example file
	SystemMatrix_Mean->Init(DO_Mean_Equation_Coefficients, BoundCondition, BoundValue);

    // -0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0---0-0--0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0--00-0-0-0-0-0-0-0-0-0-0-0--0-0-0-0-//
    //------------------------------------------ MEAN EQUATION SETUP -----------------------------------------------------//
    // -0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0---0-0--0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0--00-0-0-0-0-0-0-0-0-0-0-0--0-0-0-0-//

    int N_terms_Mean = 3;  // Number of Shape function derivatives required ( in this case 3, N, NX, NY )
    MultiIndex2D Derivatives_MatrixARhs_Mean[3] = { D00, D10, D01};
    int SpacesNumbers_MatrixARhs_Mean[3] = { 0, 0, 0  };
    int N_Matrices_MatrixARhs_Mean = 1;
    int RowSpace_MatrixARhs_Mean[1] = { 0 };
    int ColumnSpace_MatrixARhs_Mean[1] = { 0 };
    int N_Rhs_MatrixARhs_Mean = 1;
    int RhsSpace_MatrixARhs_Mean[1] = { 0 };


    SystemMatrix_Mean->Init_WithDiscreteform(DO_Mean_Equation_Coefficients,BoundCondition, BoundValue,"DO_LINEAR_MEAN", "DO_LINEAR_MEAN",
                                                N_terms_Mean, Derivatives_MatrixARhs_Mean, SpacesNumbers_MatrixARhs_Mean,
                                                N_Matrices_MatrixARhs_Mean, N_Rhs_MatrixARhs_Mean,
                                                RowSpace_MatrixARhs_Mean,ColumnSpace_MatrixARhs_Mean, RhsSpace_MatrixARhs_Mean,
                                                DO_Mean_Equation_Assembly, DO_Mean_Equation_Coefficients,
                                                NULL);

    // Aux Setup for the RHS -- There is no Aux for the Mean equation, So set the values as NULL
    fesp[0] = Scalar_FeSpace;
    aux = new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);

    /* -0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0---0-0--0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0--00-0-0-0-0-0-0-0-0-0-0-0--0-0-0-0-
    --------------------------------------[[[ END  ]]] MEAN EQUATION SETUP -----------------------------------------------------
     -0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0---0-0--0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0--00-0-0-0-0-0-0-0-0-0-0-0--0-0-0-0- */



    // -0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0---0-0--0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0--00-0-0-0-0-0-0-0-0-0-0-0--0-0-0-0-//
    //-------------------------------------- MODE EQUATION SETUP -----------------------------------------------------//
    // -0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0---0-0--0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0--00-0-0-0-0-0-0-0-0-0-0-0--0-0-0-0-//


    int N_terms_Mode = 3;  // Number of Shape function derivatives required ( in this case 3, N, NX, NY )
    MultiIndex2D Derivatives_MatrixARhs_Mode[3] = { D00, D10, D01};
    int SpacesNumbers_MatrixARhs_Mode[3] = { 0, 0, 0  };
    int N_Matrices_MatrixARhs_Mode = 1;
    int RowSpace_MatrixARhs_Mode[1] = { 0 };
    int ColumnSpace_MatrixARhs_Mode[1] = { 0 };
    int N_Rhs_MatrixARhs_Mode = 1;
    int RhsSpace_MatrixARhs_Mode[1] = { 0 };


    SystemMatrix_Mode->Init_WithDiscreteform(DO_Mode_Equation_Coefficients,BoundCondition, BoundValue,"DO_LINEAR_Mode", "DO_LINEAR_Mode",
                                                N_terms_Mode, Derivatives_MatrixARhs_Mode, SpacesNumbers_MatrixARhs_Mode,
                                                N_Matrices_MatrixARhs_Mode, N_Rhs_MatrixARhs_Mode,
                                                RowSpace_MatrixARhs_Mode,ColumnSpace_MatrixARhs_Mode, RhsSpace_MatrixARhs_Mode,
                                                DO_Mode_Equation_Assembly, DO_Mode_Equation_Coefficients,
                                                NULL);

    // Set up a FE VECT FUNCTION TO STORE ALL THE Components of CTilde  

 

    int TimeLinear_FESpaces_DO = 1;
	int TimeLinear_Fct_DO = 1; // \tilde(C)
	int TimeLinear_ParamFct_DO = 1;
	int TimeLinear_FEValues_DO = 3;
	int TimeLinear_Params_DO = 3;
	int TimeNSFEFctIndex_DO[3] = {0, 0, 0};  
	MultiIndex2D TimeNSFEMultiIndex_DO[3] = {D00, D01, D10};
	ParamFct *TimeNSFct_DO[1] 		= { DO_Mode_RHS_Aux_Param };
	int TimeNSBeginParam_DO[1] 		= { 0 };

    TFEFunction2D* fefct_RHS[4];
	TFESpace2D* fesp_RHS[2];


    // Set Up Aux  Param for the given MOde 
    // The mode equation needs the entire values of the C_tilde matrix array 
    TAuxParam2D* aux_RHS_DO = new TAuxParam2D (	TimeLinear_FESpaces_DO, <FE VECT FUNCTION>, TimeLinear_ParamFct_DO,
												TimeLinear_FEValues_DO,
												fesp_RHS,
												TimeNSFct_DO,
												TimeNSFEMultiIndex_DO,
												TimeLinear_Params_DO, TimeNSBeginParam_DO);




    // -0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0---0-0--0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0--00-0-0-0-0-0-0-0-0-0-0-0--0-0-0-0-//
    //--------------------------------------[[[ END  ]]] MODE EQUATION SETUP -----------------------------------------------------//
    // -0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0---0-0--0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0--00-0-0-0-0-0-0-0-0-0-0-0--0-0-0-0-//


	// assemble the system matrix with given aux, sol and rhs
	// aux is used to pass  addition fe functions (eg. mesh velocity) that is nedded for assembling,
	// otherwise, just pass with NULL
	SystemMatrix_Mean->AssembleMRhs(NULL, sol, rhs);

	//======================================================================
	// produce outout at t=0
	//======================================================================
	VtkBaseName = TDatabase::ParamDB->VTKBASENAME;
	Output = new TOutput2D(2, 2, 1, 1, Domain);

	Output->AddFEFunction(Scalar_FeFunction_Mean);



	//     Scalar_FeFunction->Interpolate(Exact);
	if (TDatabase::ParamDB->WRITE_VTK)
	{
		os.seekp(std::ios::beg);
		if (img < 10)
			os << "VTK/" << VtkBaseName << ".0000" << img << ".vtk" << ends;
		else if (img < 100)
			os << "VTK/" << VtkBaseName << ".000" << img << ".vtk" << ends;
		else if (img < 1000)
			os << "VTK/" << VtkBaseName << ".00" << img << ".vtk" << ends;
		else if (img < 10000)
			os << "VTK/" << VtkBaseName << ".0" << img << ".vtk" << ends;
		else
			os << "VTK/" << VtkBaseName << "." << img << ".vtk" << ends;
		Output->WriteVtk(os.str().c_str());
		img++;
	}


	coll->GetHminHmax(&hmin, &hmax);
	OutPut("h_min : " << hmin << " h_max : " << hmax << endl);

	//    TDatabase::TimeDB->TIMESTEPLENGTH =  hmax;

	// TDatabase::TimeDB->TIMESTEPLENGTH =  hmax;
	//======================================================================
	// time disc loop
	//======================================================================
	// parameters for time stepping scheme
	m = 0;
	N_SubSteps = GetN_SubSteps();
	end_time = TDatabase::TimeDB->ENDTIME;

	UpdateStiffnessMat = FALSE; //check BilinearCoeffs in example file
	UpdateRhs = TRUE;			//check BilinearCoeffs in example file
	ConvectionFirstTime = TRUE;

	// time loop starts
	while (TDatabase::TimeDB->CURRENTTIME < end_time)
	{
		m++;
		TDatabase::TimeDB->INTERNAL_STARTTIME = TDatabase::TimeDB->CURRENTTIME;

		for (l = 0; l < N_SubSteps; l++) // sub steps of fractional step theta
		{
			SetTimeDiscParameters(1);

			if (m == 1)
			{
				OutPut("Theta1: " << TDatabase::TimeDB->THETA1 << endl);
				OutPut("Theta2: " << TDatabase::TimeDB->THETA2 << endl);
				OutPut("Theta3: " << TDatabase::TimeDB->THETA3 << endl);
				OutPut("Theta4: " << TDatabase::TimeDB->THETA4 << endl);
			}

			tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
			TDatabase::TimeDB->CURRENTTIME += tau;

			OutPut(endl
				   << "CURRENT TIME: ");
			OutPut(TDatabase::TimeDB->CURRENTTIME << endl);

			//copy rhs to oldrhs
			memcpy(oldrhs, rhs, N_DOF * SizeOfDouble);

			// unless the stiffness matrix or rhs change in time, it is enough to
			// assemble only once at the begning
			if (UpdateStiffnessMat || UpdateRhs || ConvectionFirstTime)
			{
				SystemMatrix_Mean->AssembleARhs(NULL, sol, rhs);

				// M:= M + (tau*THETA1)*A
				// rhs: =(tau*THETA4)*rhs +(tau*THETA3)*oldrhs +[M-(tau*THETA2)A]*oldsol
				// note! sol contains only the previous time step value, so just pass
				// sol for oldsol
				SystemMatrix_Mean->AssembleSystMat(oldrhs, sol, rhs, sol);
				ConvectionFirstTime = FALSE;
			}

			// solve the system matrix
			SystemMatrix_Mean->Solve(sol, rhs);

			// restore the mass matrix for the next time step
			// unless the stiffness matrix or rhs change in time, it is not necessary to assemble the system matrix in every time step
			if (UpdateStiffnessMat || UpdateRhs)
			{
				SystemMatrix_Mean->RestoreMassMat();
			}
		} // for(l=0;l<N_SubSteps;l++)

		//======================================================================
		// produce outout
		//======================================================================
		if (m == 1 || m % TDatabase::TimeDB->STEPS_PER_IMAGE == 0)
			if (TDatabase::ParamDB->WRITE_VTK)
			{
				os.seekp(std::ios::beg);
				if (img < 10)
					os << "VTK/" << VtkBaseName << ".0000" << img << ".vtk" << ends;
				else if (img < 100)
					os << "VTK/" << VtkBaseName << ".000" << img << ".vtk" << ends;
				else if (img < 1000)
					os << "VTK/" << VtkBaseName << ".00" << img << ".vtk" << ends;
				else if (img < 10000)
					os << "VTK/" << VtkBaseName << ".0" << img << ".vtk" << ends;
				else
					os << "VTK/" << VtkBaseName << "." << img << ".vtk" << ends;
				Output->WriteVtk(os.str().c_str());
				img++;
			}

            //======================================================================
		// measure errors to known solution
		//======================================================================
		if (TDatabase::ParamDB->MEASURE_ERRORS)
		{
            fesp[0] = Scalar_FeSpace;
            aux = new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);
			Scalar_FeFunction_Mean->GetErrors(Exact, 3, AllDerivatives, 2, L2H1Errors, DO_Mean_Equation_Coefficients, aux, 1, fesp, errors);

			OutPut("time: " << TDatabase::TimeDB->CURRENTTIME);
			OutPut(" L2: " << errors[0]);
			OutPut(" H1-semi: " << errors[1] << endl);

			errors[3] += (errors[0] * errors[0] + olderror * olderror) * TDatabase::TimeDB->TIMESTEPLENGTH / 2.0;
			olderror = errors[0];
			OutPut(TDatabase::TimeDB->CURRENTTIME << " L2(0,T;L2) " << sqrt(errors[3]) << " ");

			errors[4] += (errors[1] * errors[1] + olderror1 * olderror1) * TDatabase::TimeDB->TIMESTEPLENGTH / 2.0;
			OutPut("L2(0,T;H1) " << sqrt(errors[4]) << endl);
			olderror1 = errors[1];

			if (Linfty < errors[0])
				Linfty = errors[0];

			OutPut("Linfty " << Linfty << endl);
		} //  if(TDatabase::ParamDB->MEASURE_ERRORS)



	} // while(TDatabase::TimeDB->CURRENTTIME< end_time)

	//======================================================================
	// produce final outout
	//======================================================================

	if (TDatabase::ParamDB->WRITE_VTK)
	{
		os.seekp(std::ios::beg);
		if (img < 10)
			os << "VTK/" << VtkBaseName << ".0000" << img << ".vtk" << ends;
		else if (img < 100)
			os << "VTK/" << VtkBaseName << ".000" << img << ".vtk" << ends;
		else if (img < 1000)
			os << "VTK/" << VtkBaseName << ".00" << img << ".vtk" << ends;
		else if (img < 10000)
			os << "VTK/" << VtkBaseName << ".0" << img << ".vtk" << ends;
		else
			os << "VTK/" << VtkBaseName << "." << img << ".vtk" << ends;
		Output->WriteVtk(os.str().c_str());
		img++;
	}

	CloseFiles();

	return 0;
} // end main
