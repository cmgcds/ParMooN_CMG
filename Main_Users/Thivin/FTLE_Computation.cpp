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

#include <string.h>
#include <sstream>
#include <MooNMD_Io.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <FTLE.h>
#include<algorithm>
#include<cstring>

// =======================================================================
// include current example
// =======================================================================
#include "../Examples/TNSE_2D/DrivenCavity.h" //   in unit square
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


	//======================================================================
	// construct all finite element functions
	//======================================================================
	sol = new double[N_TotalDOF];



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


    std::ifstream file_uSol("U_Solution.txt");
    std::ifstream file_vSol("V_Solution.txt");

	std::string mainLine;
	int lines = 0;
	int valCount = 0;
    while (std::getline(file_uSol,mainLine))
    {
		lines++;
		
        std::istringstream linestream(mainLine);

        std::string value;
        int k = 0;
    
        while(std::getline(linestream,value,',') )
        {
			if(k >= N_U) break;
            double val = std::stod(value);
            solVector_U[(lines-1)*N_U + k] = val;
            k+=1;
        }
		valCount=k;
	//   cout << " NormSolU : " << Ddot(N_U,solVector_U + (lines-1)*N_U,solVector_U +  (lines-1)*N_U) <<endl;
    }
    file_uSol.close();
	// cout << " Lines Read : " << lines << " k : " << valCount <<endl;

	lines = 0;
    while (std::getline(file_vSol,mainLine,'\n'))
    {
		lines++;
        std::istringstream linestream(mainLine);

        std::string value;
        int k = 0;
        while(std::getline(linestream,value,',') && k < N_U)
        {
			if(k >= N_U) break;
            double val = std::stod(value);			
            solVector_V[(lines-1)*N_U + k] = val;
            k+=1;
        }
		valCount=k;
		// cout << " NormSolV : " << Ddot(N_U,solVector_V+ (lines-1)*N_U,solVector_V + (lines-1)*N_U) <<endl;
    }

    file_vSol.close();
	// cout << " Lines Read : " << lines <<  " k : " << valCount<<endl;

	// cout << " NormSolU : " << Ddot(N_U,solVector_U,solVector_U) <<endl;
	// cout << " NormSolV : " << Ddot(N_U,solVector_V,solVector_V) <<endl;
    
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


	// // Save all Solutions into a file for Reuse
	// std::ofstream file_uSol_2("U_Solution_2.txt");
	// std::ofstream file_vSol_2("V_Solution_2.txt");

	// for ( int i = 0 ; i < numFiles ; i++)
	// {
	// 	for( int dof =0 ; dof < N_U; dof++)
	// 	{
	// 		file_uSol_2 << solVector_U[dof]<<",";
	// 	}
	// 	file_uSol_2<<endl;
	// }
	
	// for ( int i = 0 ; i < numFiles ; i++)
	// {
	// 	for( int dof =0 ; dof < N_U; dof++)
	// 	{
	// 		file_vSol_2 << solVector_V[dof]<<",";
	// 	}
	// 	file_vSol_2<<endl;
	// }
	

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
	TFESpace2D* ftleFespace = new TFESpace2D(coll,"ftle","ftle",BoundCondition,3,NULL);
	
	FTLE *ftle = new FTLE(ftleFespace,Velocity, 4,VelocityAll_U,VelocityAll_V);

	cout << " NUMFILES : " << numFiles<<endl;
	mkdir("FTLE", 0777);
	double T = 5;
	for ( int i = 0 ; i < numFiles-T; i++)
	{
		ftle->computeFTLE(TDatabase::TimeDB->CURRENTTIMESTEPLENGTH,T,i);
		cout << " Completed for " << i <<endl;
	}
	
	CloseFiles();

	return 0;
} // end main
