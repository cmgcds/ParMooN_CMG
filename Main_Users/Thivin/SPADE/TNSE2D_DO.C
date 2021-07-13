// =======================================================================
//
// Purpose:     main program for solving a time-dependent NSE equation in ParMooN
//
// Author:      Sashikumaar Ganesan
//
// History:     Implementation started on 03.09.2014

// =======================================================================
#include <Domain.h>
#include <Database.h>
#include <FEDatabase2D.h>
#include <LinAlg.h>
#include <FESpace2D.h>
#include <SystemTNSE2D.h>
#include <SquareStructure2D.h>
#include <Structure2D.h>
#include <Output2D.h>
#include <CD2DErrorEstimator.h>
#include <MainUtilities.h>
// #include <TimeUtilities.h>
#include <TNSE2D_ParamRout.h>
#include <TimeDiscRout.h>

#include <string.h>
#include <sstream>
#include <MooNMD_Io.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>

#ifdef INTELMKLBLAS

#endif
#include<mkl.h>
#include<omp.h>
#include "../Examples/TNSE_2D/DrivenCavity.h" // smooth sol in unit square


// FOR DO
#include<random>
// Call the RK Routines for Solving the Co-efficient equation
#include"Do_RK.h"

// PAram Routing file for DO 
void TimeNSParams_DO(double *in, double *out)
{
	out = in;
}




// double* RK_4(double* phi_0, double t0)
// {
//   double* y0 = phi_0;
//   double  x0 = t0; 
//   double* k1 = h * (CoefficientEquation(x0, y0));
//   double* k2 = h * (CoefficientEquation((x0+h/2), (y0+k1/2)));
//   double* k3 = h * (CoefficientEquation((x0+h/2), (y0+k2/2)));
//   double* k4 = h * (CoefficientEquation((x0+h), (y0+k3)));
//   double* k = (k1+2*k2+2*k3+k4)/6;

//   double* yn = y0 + k;
// }





int main(int argc, char *argv[])
{
    //------- Initialise all the Essential Variables  ----- //
    TDomain *Domain;
    TDatabase *Database = new TDatabase();
    TFEDatabase2D *FEDatabase = new TFEDatabase2D();

    // Collection of Cells
    TCollection *coll, *mortarcoll = NULL;

    // -- finite Element Space fo Velocity and pressure
    TFESpace2D *Velocity_FeSpace, *Pressure_FeSpace, *fesp[2];

    // -- Fe Functions for Velocity and pressure --- //
    TFEVectFunct2D *Velocity_mean;
    TFEFunction2D *u1_mean, *u2_mean, *Pressure_mean, *fefct[2];


    // Class to Write the System output to VTK File
    TOutput2D *Output;

    // -- Parameter to pass older variables into the Assembly function --
    TAuxParam2D *aux, *NSEaux_error;

    // -- Multi Index to be used while defining the Assembly function
    MultiIndex2D AllDerivatives[3] = {D00, D10, D01};

    // Variables initialisation for output and Definitions
    const char vtkdir[] = "VTK";
    char *PsBaseName, *VtkBaseName, *GEO;
    char UString[] = "u";
    char PString[] = "p";
    char NameString[] = "VMS";

    // Set the output stream to write VTK files
    std::ostringstream os;
    os << " ";

    mkdir(vtkdir, 0777);

    // Write the accessories output file and check parameter consistencies for NSE
    OpenFiles();
    OutFile.setf(std::ios::scientific);

    Database->CheckParameterConsistencyNSE();
    Database->WriteParamDB(argv[0]);
    Database->WriteTimeDB();
    ExampleFile();

    // ------ DOMAIN CREATION ------------------ //
    // Create a Domain based on the input mesh File provided.
    Domain = new TDomain(argv[1]);
    GEO = TDatabase::ParamDB->GEOFILE;
    Domain->Init(NULL, GEO);

    // refine grid up to the coarsest level
    for (int i = 0; i < TDatabase::ParamDB->UNIFORM_STEPS; i++)
        Domain->RegRefineAll();

    int ORDER = TDatabase::ParamDB->ANSATZ_ORDER;
    int NSEType = TDatabase::ParamDB->NSTYPE;
    int Disctype = TDatabase::ParamDB->DISCTYPE;
    int pressure_space_code, velocity_space_code;
    coll = Domain->GetCollection(It_Finest, 0);
    int N_Cells = coll->GetN_Cells();
    OutPut("N_Cells : " << N_Cells << endl);

    // fespaces for velocity and pressure
    GetVelocityAndPressureSpace(coll, BoundCondition, mortarcoll, Velocity_FeSpace,
                                Pressure_FeSpace, &pressure_space_code,
                                TDatabase::ParamDB->VELOCITY_SPACE,
                                TDatabase::ParamDB->PRESSURE_SPACE);

    // defaulty inf-sup pressure space will be selected based on the velocity space, so update it in database
    TDatabase::ParamDB->INTERNAL_PRESSURE_SPACE = pressure_space_code;
    velocity_space_code = TDatabase::ParamDB->VELOCITY_SPACE;

    int N_U = Velocity_FeSpace->GetN_DegreesOfFreedom();
    int N_P = Pressure_FeSpace->GetN_DegreesOfFreedom();
    int N_TotalDOF = 2 * N_U + N_P;

    OutPut("Dof Velocity : " << setw(10) << 2 * N_U << endl);
    OutPut("Dof Pressure : " << setw(10) << N_P << endl);
    OutPut("Total Dof all: " << setw(10) << N_TotalDOF << endl);



    //--- Setup arrays and Functions for System Matrices and Arrays ( AX=b )
    double* sol     = new double[N_TotalDOF]; 
    double* rhs     = new double[N_TotalDOF];
    double* oldrhs  = new double[N_TotalDOF];
     
    memset(sol, 0, N_TotalDOF*SizeOfDouble);
    memset(rhs, 0, N_TotalDOF*SizeOfDouble);

    Velocity_mean = new TFEVectFunct2D(Velocity_FeSpace, UString,  UString,  sol, N_U, 2);
    u1_mean = Velocity_mean->GetComponent(0);
    u2_mean = Velocity_mean->GetComponent(1);
    Pressure_mean = new TFEFunction2D(Pressure_FeSpace, PString,  PString,  sol+2*N_U, N_P);
    
    //  interpolate the initial solution based on the initial condition of the problem
    u1_mean->Interpolate(InitialU1);
    u2_mean->Interpolate(InitialU2);
    Pressure_mean->Interpolate(InitialP);


    // -- Error checking ---- //

    if ( TDatabase::ParamDB->DISCTYPE != 1)
    {
        cerr << "Only Disctype 1 can be selected " <<endl;
        exit(0);
    }

    // Note : The below is very bad way to write the code  
    // `

    #include "DO_implementation.h"

    cout << " --- DO SETUP PART COMPLETED --------------" <<endl;


    // ------ Set up Discrete Forms for RHS Assembly of DO Part of NSE -------------- //
	TFEFunction2D* fefct_RHS[4];
	TFESpace2D* fesp_RHS[2];

	fesp_RHS[0] = Velocity_FeSpace;

	double* fluctVelocity = new double[N_U*2]();

	// Set up a vectfunction for Instentaneous Velocity
	TFEVectFunct2D* FluctuatingVelocity = new TFEVectFunct2D(Velocity_FeSpace, UString,  UString,  FluctuatingVelocityArray, N_U, 2*N_S);

	int TimeNSN_FESpaces_DO = 1;
	int TimeNSN_Fct_DO = 2; // \tilde(u,v), \bar(u,v)
	int TimeNSN_ParamFct_DO = 1;
	int TimeNSN_FEValues_DO = 3;
	int TimeNSN_Params_DO = 3;
	int TimeNSFEFctIndex_DO[3] = {0, 0, 0};  
	MultiIndex2D TimeNSFEMultiIndex_DO[3] = {D00, D01, D10};
	ParamFct *TimeNSFct_DO[1] 		= { TimeNSParams_DO };
	int TimeNSBeginParam_DO[1] 		= { 0 };


	TAuxParam2D* aux_RHS_DO = new TAuxParam2D (	TimeNSN_FESpaces_DO, FluctuatingVelocity, TimeNSN_ParamFct_DO,
												TimeNSN_FEValues_DO,
												fesp_RHS,
												TimeNSFct_DO,
												TimeNSFEMultiIndex_DO,
												TimeNSN_Params_DO, TimeNSBeginParam_DO);
	
    // Define the Aux for the TNSE Function
    fesp[0] = Velocity_FeSpace;
    fefct[0] = u1_mean;
    fefct[1] = u2_mean;  

    aux =  new TAuxParam2D(TimeNSN_FESpaces2, TimeNSN_Fct2, TimeNSN_ParamFct2,
                           TimeNSN_FEValues2,
                           fesp, fefct,
                           TimeNSFct2,
                           TimeNSFEFctIndex2, TimeNSFEMultiIndex2,
                           TimeNSN_Params2, TimeNSBeginParam2);  
 

    
	TSystemTNSE2D* SystemMatrix = new TSystemTNSE2D(Velocity_FeSpace, Pressure_FeSpace, Velocity_mean, Pressure_mean, sol, rhs, Disctype, NSEType, DIRECT
	#ifdef __PRIVATE__
										,
										Projection_space, NULL, NULL
	#endif
	);

	// initilize the system matrix with the functions defined in Example file
	// last argument is aux that is used to pass additional fe functions (eg. mesh velocity)
	SystemMatrix->Init(LinCoeffs, BoundCondition, U1BoundValue, U2BoundValue, aux, NSEaux_error,aux_RHS_DO);

	// assemble M, A matrices and rhs
	SystemMatrix->Assemble(sol, rhs);

    cout << " Assembly Done " <<endl;

	//======================================================================
	// produce outout
	//======================================================================
	VtkBaseName = TDatabase::ParamDB->VTKBASENAME;
	Output = new TOutput2D(2, 2, 1, 1, Domain);

	Output->AddFEVectFunct(Velocity_mean);
	Output->AddFEFunction(Pressure_mean);

	int img = 0;
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
	// time disc loop
	//======================================================================
	// parameters for time stepping scheme
	int m = 0;
	int N_SubSteps = GetN_SubSteps();
	double oldtau = 1.;
	double end_time = TDatabase::TimeDB->ENDTIME;
	int limit = TDatabase::ParamDB->SC_NONLIN_RES_NORM_MIN_SADDLE;
	int Max_It = TDatabase::ParamDB->SC_NONLIN_MAXIT_SADDLE;
	double tau;
	double* defect;
	double residual;
	double impuls_residual;

	// time loop starts
	while (TDatabase::TimeDB->CURRENTTIME < end_time)
	{ // time cycle
		m++;
		TDatabase::TimeDB->INTERNAL_STARTTIME = TDatabase::TimeDB->CURRENTTIME;
		for (int l = 0; l < N_SubSteps; l++) // sub steps of fractional step theta
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

			OutPut(endl<< "CURRENT TIME: ");
			OutPut(TDatabase::TimeDB->CURRENTTIME << endl);

			//copy sol, rhs to olssol, oldrhs
			memcpy(oldrhs, rhs, N_TotalDOF * SizeOfDouble);

			// assemble only rhs, nonlinear matrix for NSE will be assemble in fixed point iteration
			// not needed if rhs is not time-dependent
			if (m != 1)
			{
				SystemMatrix->AssembleRhs(sol, rhs);
			}
			else
			{
				SystemMatrix->Assemble(sol, rhs);
			}

			//scale B matices and assemble NSE-rhs based on the \theta time stepping scheme
			SystemMatrix->AssembleSystMat(tau / oldtau, oldrhs, rhs, sol);
			oldtau = tau;

			// calculate the residual
			defect = new double[N_TotalDOF];
			memset(defect, 0, N_TotalDOF * SizeOfDouble);

			SystemMatrix->GetTNSEResidual(sol, defect);

			//correction due to L^2_O Pressure space
			if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
				IntoL20Vector2D(defect + 2 * N_U, N_P, pressure_space_code);

			residual = Ddot(N_TotalDOF, defect, defect);
			impuls_residual = Ddot(2 * N_U, defect, defect);
			OutPut("Nonlinear iteration step   0");
			OutPut(setw(14) << impuls_residual);
			OutPut(setw(14) << residual - impuls_residual);
			OutPut(setw(14) << sqrt(residual) << endl);

			//======================================================================
			//Solve the system
			//Nonlinear iteration of fixed point type
			//======================================================================
			for (int j = 1; j <= Max_It; j++)
			{
				// Solve the NSE system
				SystemMatrix->Solve(sol);

				if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
					IntoL20FEFunction(sol + 2 * N_U, N_P, Pressure_FeSpace, velocity_space_code, pressure_space_code);

				//no nonlinear iteration for Stokes problem
				if (TDatabase::ParamDB->FLOW_PROBLEM_TYPE == STOKES)
					break;

				// restore the mass matrix for the next nonlinear iteration
				SystemMatrix->RestoreMassMat();

				// assemble the system matrix with given aux, sol and rhs
				SystemMatrix->AssembleANonLinear(sol, rhs);

				// assemble system mat, S = M + dt\theta_1*A
				SystemMatrix->AssembleSystMatNonLinear();

				// get the residual
				memset(defect, 0, N_TotalDOF * SizeOfDouble);
				SystemMatrix->GetTNSEResidual(sol, defect);

				//correction due to L^2_O Pressure space
				if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
					IntoL20Vector2D(defect + 2 * N_U, N_P, pressure_space_code);

				residual = Ddot(N_TotalDOF, defect, defect);
				impuls_residual = Ddot(2 * N_U, defect, defect);
				OutPut("nonlinear iteration step " << setw(3) << j);
				OutPut(setw(14) << impuls_residual);
				OutPut(setw(14) << residual - impuls_residual);
				OutPut(setw(14) << sqrt(residual) << endl);

				if (sqrt(residual) <= limit)
					break;

			} // for(j=1;j<=Max_It;j++)
			  /*           cout << " test VHM main " << endl;
                    exit(0);      */
			// restore the mass matrix for the next time step
			SystemMatrix->RestoreMassMat();

		} // for(l=0;l<N_SubSteps;




        

        // get the Co-ordinates of the Current Grid Point
        double* gridCoordArray = new double(N_U*2);
        TFEVectFunct2D* gridCoordFEVect = new TFEVectFunct2D(Velocity_FeSpace, UString,  UString,  gridCoordArray, N_U, 2);
        // Obtains the Co-ordinate values of the Current Grid
        gridCoordFEVect->GridToData();

        double* X_grid = gridCoordFEVect->GetComponent(0)->GetValues();
        double* Y_grid = gridCoordFEVect->GetComponent(1)->GetValues();


        // List of variables needed for RK-4
        //Dx_p_titlde,Dy_p_titlde
        //u_titlde, Dx_u_titlde, DDx_u_titlde, Dy_u_titlde, DDy_u_titlde
        //v_titlde , Dx_v_titlde, DDx_v_titlde, Dy_v_titlde, DDy_v_titlde,
        //u_bar,v_bar,Dx_u_bar,Dy_u_bar,Dx_v_bar,Dy_v_bar

		// 
		
		// 
		

		double* p_tilde = new double[N_P*N_S];
		TFEVectFunct2D* p_tilde_vectFunction = new TFEVectFunct2D(Pressure_FeSpace,PString,PString,p_tilde,N_P,N_S);

		// p_tilde 
		double* p_tilde_extended = new double[N_U * N_S];

		//Dx_p_tilde 
		double* Dx_p_tilde = new double[N_U * N_S];

		//Dy_p_tilde 
		double* Dy_p_tilde = new double[N_U * N_S];

		// u_tilde      -- is a FEVectFunction with N_S Components 
		double* u_tilde  = new double[N_U * N_S];
		u_tilde = ModeVector; 
		TFEVectFunct2D* u_tildeVectfunction = new TFEVectFunct2D(Velocity_FeSpace,UString,UString,u_tilde,N_U,N_S);
		
		// v_tilde      -- is a vectfunction with N_S Components
		double* v_tilde  = new double[N_U * N_S];
		v_tilde = ModeVector_2;
		TFEVectFunct2D* v_tildeVectfunction = new TFEVectFunct2D(Velocity_FeSpace,UString,UString,v_tilde,N_U,N_S);
		
		// Dx_u_tilde   -- is a Array with N_S Componetns
		double* Dx_u_tilde  = new double[N_U * N_S];

		// Dy_u_tilde   -- is a Array with N_S Componetns
		double* Dy_u_tilde  = new double[N_U * N_S];

		// Dx_v_tilde   -- is a Array with N_S Componetns
		double* Dx_v_tilde  = new double[N_U * N_S];

		// Dy_v_tilde   -- is a Array with N_S Componetns
		double* Dy_v_tilde  = new double[N_U * N_S];

		// DDx_u_tilde   -- is a Array with N_S Componetns
		double* DDx_u_tilde  = new double[N_U * N_S];

		// DDy_u_tilde   -- is a Array with N_S Componetns
		double* DDy_u_tilde  = new double[N_U * N_S];

		// DDx_v_tilde   -- is a Array with N_S Componetns
		double* DDx_v_tilde  = new double[N_U * N_S];

		// DDy_v_tilde   -- is a Array with N_S Componetns
		double* DDy_v_tilde  = new double[N_U * N_S];

		// Dx_u_bar   -- is a Array with N_S Componetns
		double* Dx_u_bar  = new double[N_U * N_S];

		// Dx_v_bar   -- is a Array with N_S Componetns
		double* Dx_v_bar  = new double[N_U * N_S];

		// Dy_u_bar   -- is a Array with N_S Componetns
		double* Dy_u_bar  = new double[N_U * N_S];

		// Dy_v_bar   -- is a Array with N_S Componetns
		double* Dy_v_bar  = new double[N_U * N_S];

		double* u_bar = u_tildeVectfunction->GetComponent(0)->GetValues();
		double* v_bar = u_tildeVectfunction->GetComponent(1)->GetValues();


        int N_Coord = gridCoordFEVect->GetComponent(0)->GetLength();


		// Obtain the values at all the co-ordinate points 
		// Need to find Derivatives 
		for ( int j = 0 ; j < N_S; j++)
		{
			for( int i = 0 ; i < N_Coord; i++)
			{
				double x = X_grid[i];
				double y = Y_grid[i];
				double* values = new double[3];

				u_tildeVectfunction->GetComponent(j)->FindGradient(x,y,values);
				Dx_u_tilde[i] = values[1];
				Dy_u_tilde[i] = values[2];

				v_tildeVectfunction->GetComponent(j)->FindGradient(x,y,values);
				Dx_u_tilde[i] = values[1];
				Dy_u_tilde[i] = values[2];

				u1_mean->FindGradient(x,y,values);
				Dx_u_bar[i] = values[1];
				Dy_u_bar[i] = values[2];


				u2_mean->FindGradient(x,y,values);
				Dx_v_bar[i] = values[1];
				Dy_v_bar[i] = values[2];

				//Pressure
				p_tilde_vectFunction->GetComponent(j)->FindGradient(x,y,values);
				Dx_p_tilde[i] = values[1];
				Dy_p_tilde[i] = values[2];

			}

		}

		double h = 0.1;
		double t0;

		//CoeffVector    -- make sure to send the ptilde extended version for p_tilde
		Solve_RK4( h, CoeffVector, t0, u_tilde, v_tilde, 
				 Dx_u_tilde, Dy_u_tilde, DDx_u_tilde, DDy_u_tilde,
				 DDx_v_tilde, DDy_v_tilde,
				 Dx_v_tilde, Dy_v_tilde ,u_bar , v_bar, Dx_u_bar,Dy_u_bar,
				 Dx_v_bar,Dy_v_bar, p_tilde_extended, Dx_p_tilde, Dy_p_tilde , Cov);


	


        // Producing Output
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

	} // while(TDatabase::TimeDB->CURRENTTIME< e

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
}