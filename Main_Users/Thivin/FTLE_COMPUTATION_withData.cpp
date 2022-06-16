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
#include <MainUtilities.h> // #include <TimeUtilities.h>
#include <TNSE2D_ParamRout.h>
#include <TimeDiscRout.h>
#include<omp.h>

#include <string.h>
#include <sstream>
#include <MooNMD_Io.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <FTLE.h>
#include<algorithm>

// =======================================================================
// include current example
// =======================================================================
#include "../Main_Users/Thivin/DoubleGyre.h" //   in unit square
// #include "../Examples/TNSE_2D/Bsp3.h" // smooth sol in unit square
// #include "../Examples_All/TNSE_2D/Benchmark2.h"
// #include "../Examples/TNSE_2D/SinCos.h" // smooth sol in unit square
// #include "../Examples/TNSE_2D/channel.h" // smooth sol in unit square
// =======================================================================
// main program
// =======================================================================
int main(int argc, char *argv[])
{
	// ======================================================================
	//  declaration of variables
	// ======================================================================
	int i, j, l, m, N_Cells, ORDER, N_U, N_P, N_L, N_TotalDOF, img = 1, pressure_space_code;
	int Max_It, NSEType, velocity_space_code, N_SubSteps, Disctype;

	double *sol, *rhs, *oldrhs, *defect, t1, t2, residual, impuls_residual;
	double limit, AllErrors[7], end_time, oldtau, tau;

	TDomain *Domain;
	TDatabase *Database = new TDatabase();
	TFEDatabase2D *FEDatabase = new TFEDatabase2D();
	TCollection *coll, *mortarcoll = NULL;
	TFESpace2D *Velocity_FeSpace, *Pressure_FeSpace, *fesp[2];
	TFEVectFunct2D *Velocity;
	TFEFunction2D *u1, *u2, *Pressure, *fefct[2];
	TOutput2D *Output;
	TSystemTNSE2D *SystemMatrix;
	TAuxParam2D *aux, *NSEaux_error;
	MultiIndex2D AllDerivatives[3] = {D00, D10, D01};
	std::vector<double *> U_Velocity_Array;
	std::vector<double *> V_Velocity_Array;
#ifdef __PRIVATE__
	TFESpace2D *Projection_space;
#endif

	const char vtkdir[] = "VTK";
	char *PsBaseName, *VtkBaseName, *GEO;
	char UString[] = "u";
	char PString[] = "p";
	char NameString[] = "VMS";

	std::ostringstream os;
	os << " ";

	mkdir(vtkdir, 0777);

	omp_set_num_threads(24);

	// ======================================================================
	// set the database values and generate mesh
	// ======================================================================
	/** set variables' value in TDatabase using argv[1] (*.dat file), and generate the MESH based */
	Domain = new TDomain(argv[1]);

	OpenFiles();
	OutFile.setf(std::ios::scientific);

	Database->CheckParameterConsistencyNSE();
	Database->WriteParamDB(argv[0]);
	Database->WriteTimeDB();
	ExampleFile();

	/* include the mesh from a meshgenerator, for a standard mesh use the build-in function */
	// standard mesh
	GEO = TDatabase::ParamDB->GEOFILE;
	Domain->Init(NULL, GEO);

	// refine grid up to the coarsest level
	for (i = 0; i < TDatabase::ParamDB->UNIFORM_STEPS; i++)
		Domain->RegRefineAll();

	if (TDatabase::ParamDB->WRITE_PS)
	{
		// write grid into an Postscript file
		os.seekp(std::ios::beg);
		os << "Domain"
		   << ".ps" << ends;
		Domain->PS(os.str().c_str(), It_Finest, 0);
	}

	//=========================================================================
	// construct all finite element spaces
	//=========================================================================
	ORDER = TDatabase::ParamDB->ANSATZ_ORDER;
	NSEType = TDatabase::ParamDB->NSTYPE;
	Disctype = TDatabase::ParamDB->DISCTYPE;

	coll = Domain->GetCollection(It_Finest, 0);
	N_Cells = coll->GetN_Cells();
	OutPut("N_Cells : " << N_Cells << endl);

	// fespaces for velocity and pressure
	GetVelocityAndPressureSpace(coll, BoundCondition, mortarcoll, Velocity_FeSpace,
								Pressure_FeSpace, &pressure_space_code,
								TDatabase::ParamDB->VELOCITY_SPACE,
								TDatabase::ParamDB->PRESSURE_SPACE);

	// defaulty inf-sup pressure space will be selected based on the velocity space, so update it in database
	TDatabase::ParamDB->INTERNAL_PRESSURE_SPACE = pressure_space_code;
	velocity_space_code = TDatabase::ParamDB->VELOCITY_SPACE;

	N_U = Velocity_FeSpace->GetN_DegreesOfFreedom();
	N_P = Pressure_FeSpace->GetN_DegreesOfFreedom();
	N_TotalDOF = 2 * N_U + N_P;

	OutPut("Dof Velocity : " << setw(10) << 2 * N_U << endl);
	OutPut("Dof Pressure : " << setw(10) << N_P << endl);
	OutPut("Total Dof all: " << setw(10) << N_TotalDOF << endl);

#ifdef __PRIVATE__
	if (Disctype == VMS_PROJECTION)
	{
		if (TDatabase::ParamDB->VMS_LARGE_VELOCITY_SPACE == 0)
			Projection_space = new TFESpace2D(coll, NameString, UString, BoundCondition, DiscP_PSpace, 0, mortarcoll);
		else
			Projection_space = new TFESpace2D(coll, NameString, UString, BoundCondition, DiscP_PSpace, 1, mortarcoll);

		N_L = Projection_space->GetN_DegreesOfFreedom();
		OutPut("Dof Projection : " << setw(10) << N_L << endl);
	}
#endif

	//======================================================================
	// construct all finite element functions
	//======================================================================
	sol = new double[N_TotalDOF];
	rhs = new double[N_TotalDOF];
	oldrhs = new double[N_TotalDOF];

	memset(sol, 0, N_TotalDOF * SizeOfDouble);
	memset(rhs, 0, N_TotalDOF * SizeOfDouble);

	Velocity = new TFEVectFunct2D(Velocity_FeSpace, UString, UString, sol, N_U, 2);
	u1 = Velocity->GetComponent(0);
	u2 = Velocity->GetComponent(1);
	Pressure = new TFEFunction2D(Pressure_FeSpace, PString, PString, sol + 2 * N_U, N_P);

	//  interpolate the initial solution
	u1->Interpolate(InitialU1);
	u2->Interpolate(InitialU2);
	Pressure->Interpolate(InitialP);

	U_Velocity_Array.emplace_back(u1->GetValues());
	V_Velocity_Array.emplace_back(u2->GetValues());

	// I
	int numFiles = (TDatabase::TimeDB->ENDTIME - TDatabase::TimeDB->STARTTIME) / (TDatabase::TimeDB->TIMESTEPLENGTH) + 1;
	double *solVector_U = new double[numFiles*N_U]();
	double *solVector_V = new double[numFiles*N_U]();
	
	
	//======================================================================
	// SystemMatrix construction and solution
	//======================================================================
	// Disc type: GALERKIN
	// Solver: AMG_SOLVE (or) GMG  (or) DIRECT
	
	double hmin, hmax;
	coll->GetHminHmax(&hmin, &hmax);
	OutPut("h_min : " << hmin << " h_max : " << hmax << endl);
	//      TDatabase::TimeDB->TIMESTEPLENGTH = hmin;
	//       cout<<TDatabase::TimeDB->TIMESTEPLENGTH<<"\n";

	//======================================================================
	// produce outout
	//======================================================================
	VtkBaseName = TDatabase::ParamDB->VTKBASENAME;
	Output = new TOutput2D(2, 2, 1, 1, Domain);

	Output->AddFEVectFunct(Velocity);
	Output->AddFEFunction(Pressure);

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
	m = 0;
	N_SubSteps = GetN_SubSteps();
	oldtau = 1.;
	end_time = TDatabase::TimeDB->ENDTIME;
	limit = TDatabase::ParamDB->SC_NONLIN_RES_NORM_MIN_SADDLE;
	Max_It = TDatabase::ParamDB->SC_NONLIN_MAXIT_SADDLE;
	memset(AllErrors, 0, 7. * SizeOfDouble);

	// time loop starts
	while (TDatabase::TimeDB->CURRENTTIME < end_time)
	{ // time cycle
		m++;
        SetTimeDiscParameters(1);
        tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
        TDatabase::TimeDB->CURRENTTIME += tau;

        //Update the solution at everytimestep using the interpolation Initial conditoin
        u1->Interpolate(InitialU1);
        u2->Interpolate(InitialU2);
        Pressure->Interpolate(InitialP);
		
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
		
		//Copy the current Solution to the Solution Array 
		cout << " M Val : " << m << endl;
		for ( int i = 0 ; i < N_U ; i++)
		{
			solVector_U[(m-1)*N_U + i] = sol[i];
			solVector_V[(m-1)*N_U + i] = sol[N_U + i];
		}
		
		// cout << " SolActual_u : " << Ddot(N_U, sol, sol) <<endl;
		// cout << " SolActual_v : " << Ddot(N_U, sol+N_U, sol+N_U) <<endl;
		// cout << " SolActual_u : " << Ddot(N_U, solVector_U + (m-1)*N_U, solVector_U + (m-1)*N_U) <<endl;
		// cout << " SolActual_u : " << Ddot(N_U, solVector_V + (m-1)*N_U, solVector_V + (m-1)*N_U) <<endl;


	} // while(TDatabase::TimeDB->CURRENTTIME< endl
	
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

	// Save all Solutions into a file for Reuse
	std::ofstream file_uSol("U_Solution.txt");
	std::ofstream file_vSol("V_Solution.txt");

	for ( int i = 0 ; i < numFiles ; i++)
	{
		for( int dof =0 ; dof < N_U; dof++)
		{
			file_uSol << solVector_U[i*N_U + dof]<<",";
		}
		file_uSol<<endl;
	}
	
	for ( int i = 0 ; i < numFiles ; i++)
	{
		for( int dof =0 ; dof < N_U; dof++)
		{
			file_vSol << solVector_V[i*N_U + dof]<<",";
		}
		file_vSol<<endl;
	}
	

	//FEVectfunction for all components. 
	
	TFEVectFunct2D* VelocityAll_U = new TFEVectFunct2D(Velocity_FeSpace,NameString,NameString,solVector_U,N_U,numFiles);
	TFEVectFunct2D* VelocityAll_V = new TFEVectFunct2D(Velocity_FeSpace,NameString,NameString,solVector_V,N_U,numFiles);

	// for ( int i = 0 ; i < numFiles;i++)
	// {
	// 	double* a = VelocityAll_U->GetComponent(i)->GetValues();
	// 	double* b = VelocityAll_V->GetComponent(i)->GetValues();
	// 	int len = N_U;

	// 	cout << " NOrm U : " << Ddot(N_U,a,a)  << " Norm V : " << Ddot(N_U,b,b)<<endl;
	// }
	// exit(0);

	// Create more particles, using higher order FEspace
	TFESpace2D* ftleFespace = new TFESpace2D(coll,"ftle","ftle",BoundCondition,4,NULL);
	
	FTLE *ftle = new FTLE(ftleFespace,Velocity, 5,VelocityAll_U,VelocityAll_V);

	cout << " NUMFILES : " << numFiles<<endl;

    //VTK PARAMETERS
    TOutput2D* OutputFTLE;
    double* FTLEValues = ftle->FTLEValues.data();
    int N_Particles = ftle->N_Particles;
    TFEFunction2D* FTLE_FeFunction = new TFEFunction2D(ftleFespace, (char *)"FTLE", (char *)"FTLE", FTLEValues, N_Particles);
    OutputFTLE = new TOutput2D(2, 2, 1, 1, Domain);
    
    OutputFTLE->AddFEFunction(FTLE_FeFunction);

	double T = int(1 / (TDatabase::TimeDB->CURRENTTIMESTEPLENGTH))*15;
    img = 0;
    mkdir("FTLE", 0777);

	for ( int i = 0 ; i < 1; i++)
	{
		cout << " Started for " << i <<endl;
		double time = omp_get_wtime();
		ftle->computeFTLE(TDatabase::TimeDB->CURRENTTIMESTEPLENGTH,T,i);
		cout << " Time taken : " << omp_get_wtime() - time << " sec"<<endl;
		
        os.seekp(std::ios::beg);
		if (img < 10)
			os << "FTLE/" << "FTLE" << ".0000" << img << ".vtk" << ends;
		else if (img < 100)
			os << "FTLE/" << "FTLE" << ".000" << img << ".vtk" << ends;
		else if (img < 1000)
			os << "FTLE/" << "FTLE" << ".00" << img << ".vtk" << ends;
		else if (img < 10000)
			os << "FTLE/" << "FTLE" << ".0" << img << ".vtk" << ends;
		else
			os << "FTLE/" << "FTLE" << "." << img << ".vtk" << ends;
		OutputFTLE->WriteVtk(os.str().c_str());
		img++;

        cout << " Completed for " << i <<endl;
	}
    // Output for VTK 
	CloseFiles();

	return 0;
} // end main
