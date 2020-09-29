#include <Domain.h>
#include <Database.h>
#include <FEDatabase3D.h>
#include <SquareMatrix3D.h>
#include <SquareStructure3D.h>
#include <DiscreteForm3D.h>
#include <AssembleMat3D.h>
#include <Assemble3D.h>
#include <DirectSolver.h>
#include <FESpace3D.h>
#include <Output3D.h>
#include <DirectSolver.h>
#include <string.h>
#include <math.h>
#include <ctime>
#include <LinAlg.h>
#include <FEVectFunct3D.h>
#include <sys/stat.h>
#include <vector>
#include <exception>
#include <fstream>
#include <LinAlg.h>
#include "DeformMesh3D.h"
#include <FE3D_ALE.h>


using namespace std ;


void FE_BoundCondition(int BdID, double x, double y, double z, BoundCond &cond)
{
	if(BdID == 3 || BdID == 0)
	{
		cond = DIRICHLET;
		//cout<<" BD ID for Dirichlet : "<<BdID<<endl;
	}
	else 
		cond = NEUMANN;
}



void ALE_BoundCondition(int BdID, double x, double y, double z, BoundCond &cond)
{
	if(BdID == 3 || BdID == 0)
	{
		cond = DIRICHLET;
		//cout<<" BD ID for Dirichlet : "<<BdID<<endl;
	}
	else 
		cond = NEUMANN;
}

// Boundary Values for all the dimensions of the boundary domain. 
void ALE_BoundValue_X(int BdComp, double x, double y, double z, double &value)
{
    value = 0;
}

void ALE_BoundValue_Y(int BdComp, double x, double y, double z, double &value)
{
  		value = 0;
}

void ALE_BoundValue_Z(int BdComp, double x, double y, double z, double &value)
{
  value = 0;
}

void Grid_BoundCondition(int BdID, double x, double y, double z, BoundCond &cond)
{
	if(BdID == 3 || BdID == 0)
	{
		cond = DIRICHLET;
		//cout<<" BD ID for Dirichlet : "<<BdID<<endl;
	}
	else 
		cond = NEUMANN;
}

// Boundary Values for all the dimensions of the boundary domain. 
void Grid_BoundValue_X(int BdComp, double x, double y, double z, double &value)
{
    value = 0;
}

void Grid_BoundValue_Y(int BdComp, double x, double y, double z, double &value)
{
	if(BdComp == 3)
	{
		value = 0.4*x ;
	}
	else
	{
		value = 0;
	}
}

void Grid_BoundValue_Z(int BdComp, double x, double y, double z, double &value)
{
  value = 0;
}

void GridCoeffs(int n_points, double *x, double *y,double *z,
        double **parameters, double **coeffs)
{
  int i;
  double *coeff;
  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    coeff[0] = 1;    coeff[1] = 0;  coeff[2] = 0; coeff[3] = 0;  coeff[4] = 0;
  }
}


void getFreeSurfaceBoundaryIds(std::vector<int>& Boundids)
{
	int size = 1;    // NOTE : Enter the Size of the Boundary ID in the given 
	// Boundids->resize(size);  // number of bd'ids to be considered as Free Surface
	Boundids.push_back(3);
}

void getFreeSlipBoundaryIds(std::vector<int>& FreeSlipBoundids)
{
	int size = 4;    // NOTE : Enter the Size of the Boundary ID in the given 
	FreeSlipBoundids.push_back(1);
	FreeSlipBoundids.push_back(2);
	FreeSlipBoundids.push_back(4);
	FreeSlipBoundids.push_back(5);
}







int main (int argc, char* argv[])
{
    TDatabase *Database = new TDatabase();
	TDomain *Domain = new TDomain(argv[1]);
	TFEDatabase3D *FEDatabase = new TFEDatabase3D();
	TCollection *coll;
	TFESpace3D *fesp[1],*ferhs[3];
	int ORDER;

	// Vector For storing BdId's
	

	// --- Parameters required for VTK   --- //
	const char vtkdir[] = "VTK";
	std::ostringstream os;
	os << " ";
	mkdir(vtkdir, 0777);
	char UString[] = "T";
	char NameString[] = "name";
	char CString[] = "C";
	// -- End of parameters for VTK -- //
	cout << " THIVIN -- Selected mesh type : "<< TDatabase::ParamDB->MESH_TYPE<< endl;
	// -- Get the mesh datastructure into ParMooN's Data structure
	if(TDatabase::ParamDB->MESH_TYPE==0)
		Domain->Init(TDatabase::ParamDB->BNDFILE, TDatabase::ParamDB->GEOFILE);

	else if(TDatabase::ParamDB->MESH_TYPE==1){
		// cout << "Hey Bro" << endl;
		Domain->ReadGeo(TDatabase::ParamDB->GEOFILE);
	}
	
	
	for(int i=0;i<TDatabase::ParamDB->UNIFORM_STEPS;i++)
		Domain->RegRefineAll();
	
    coll=Domain->GetCollection(It_Finest, 0);
    
	ORDER  = 1;

	FE3D *fe3d = new FE3D[coll->GetN_Cells()];
	for ( int i = 0  ; i <coll->GetN_Cells() ; i++ )    
		fe3d[i] = C_Q2_3D_H_M;


    TFESpace3D *fespace = new TFESpace3D(coll, NameString, UString, FE_BoundCondition, fe3d);    // FOR VELOCITY FE SPACE
	cout << "FESPACE Velocity - Declared " << endl;
	cout << " VELOCITY -- fespace->GetN_DegreesOfFreedom();&&&&&&&&&&&&&&&&&&  : "<< fespace->GetN_DegreesOfFreedom()  <<endl;

	FE3D *fe3d_mesh = new FE3D[coll->GetN_Cells()];
	for ( int i = 0  ; i <coll->GetN_Cells() ; i++ )    
		fe3d_mesh[i] = C_Q1_3D_H_M;

	TFESpace3D *fespace_mesh = new TFESpace3D(coll, NameString, UString, FE_BoundCondition, fe3d_mesh); 
	cout << " MESH -- fespace->GetN_DegreesOfFreedom();&&&&&&&&&&&&&&&&&&  : "<< fespace_mesh->GetN_DegreesOfFreedom()  <<endl;

	int N_DOF = fespace->GetN_DegreesOfFreedom();
	int N_DOF_mesh = fespace_mesh->GetN_DegreesOfFreedom();
	int N_Active = fespace->GetActiveBound();
	cout << "N_DOF : "<<N_DOF<<endl;
	cout << " NCells : "<<fespace->GetN_Cells();



	double *meshvelo = new double[3*N_DOF_mesh]();
	double *meshvelo_Old = new double[3*N_DOF_mesh]();
	double *Velocity = new double[3*N_DOF];




	TFEVectFunct3D *MeshVelo_FEvect = new TFEVectFunct3D(fespace_mesh, (char*)"C", (char*)"C", meshvelo, N_DOF_mesh, 3);
	TFEVectFunct3D *Velocity_FEvect = new TFEVectFunct3D(fespace, (char*)"C", (char*)"C", Velocity, N_DOF, 3);

	char *VtkBaseName;
	VtkBaseName = TDatabase::ParamDB->VTKBASENAME;
  	TOutput3D* Output = new TOutput3D(2, 2, 1, 1, Domain);
	Output->AddFEVectFunct(Velocity_FEvect);
	os.seekp(std::ios::beg);
	int output_write_counter = 0;
    os <<  "VTK/"<<VtkBaseName<<".0000" <<output_write_counter <<".vtk" << ends;
	cout<<" hi "<<  endl;
	Output->WriteVtk(os.str().c_str());
	output_write_counter =  output_write_counter + 1;


	///////////////////// --- TEST CODE BLOCK ---------------- ////////////////////////////////////////////////
	BoundCondFunct3D *GridBoundaryConditions[3];
	BoundValueFunct3D *GridBoundValues[3];

	GridBoundaryConditions[0] = Grid_BoundCondition;
	GridBoundaryConditions[1] = Grid_BoundCondition;
	GridBoundaryConditions[2] = Grid_BoundCondition;
	GridBoundValues[0] = Grid_BoundValue_X;
	GridBoundValues[1] = Grid_BoundValue_Y;
	GridBoundValues[2] = Grid_BoundValue_Z;
	
	deformMesh3D* defMesh = new deformMesh3D(coll,GridBoundaryConditions,GridBoundValues);
	defMesh->moveMesh();


	os.seekp(std::ios::beg);
    os <<  "VTK/"<<VtkBaseName<<".0000" <<output_write_counter <<".vtk" << ends;
	Output->WriteVtk(os.str().c_str());
	output_write_counter =  output_write_counter + 1;

	cout << " -------------------------- MOVED MESH IN MAIN FILE --------------------------------     " <<endl;



	///////////////////// ---END---  TEST CODE BLOCK ---------------- ////////////////////////////////////////////////

	// ------------------------- ALE Boundary Condition ------------------------------------ //

    BoundCondFunct3D *ALEBoundaryConditions[3];
	BoundValueFunct3D *ALEBoundValues[3];

	ALEBoundaryConditions[0] = FE_BoundCondition;
	ALEBoundaryConditions[1] = FE_BoundCondition;
	ALEBoundaryConditions[2] = FE_BoundCondition;
	ALEBoundValues[0] = ALE_BoundValue_X;
	ALEBoundValues[1] = ALE_BoundValue_Y;
	ALEBoundValues[2] = ALE_BoundValue_Z;
	



    // Call the class for mesh movement
    FE3D_ALE* mesh_deform = new FE3D_ALE(coll,ALEBoundaryConditions,ALEBoundValues);
    // Pick the Free surface Bdid from example file
	std::vector<int> BoundIds;
	getFreeSurfaceBoundaryIds(BoundIds);

	std::vector<int> freeSlipBoundIds;
	getFreeSlipBoundaryIds(freeSlipBoundIds);

	///////////////// ---- CODE BLOCK ( ALE MESH MOVEMENT ) --------------- /////////////////
	
	
	//Function to Pick Free Surface DOF
	mesh_deform->Pick_free_surface_DOFs(fespace_mesh, BoundIds,coll);

	// Function to pick up FreeSli Boundary DOF's
	mesh_deform->pickDOFsOfFreeSlipBoundaries(fespace_mesh,freeSlipBoundIds,BoundIds);

	// Function to Calculate the Mesh Velocity at the Surface based on the current Velocity(n+1)
	mesh_deform->get_surface_normals(MeshVelo_FEvect,Velocity_FEvect);    

	mesh_deform->get_surface_normals_slipBoundary(MeshVelo_FEvect);

	mesh_deform->get_surface_normals_slipBoundary_EdgeNodes(MeshVelo_FEvect);

	// Copy the New Values of The Mesh Velocity to Mesh Velcity Old before Calculating the Interior Nodes Solution
	memcpy(meshvelo_Old,meshvelo,sizeof(double)*3*N_DOF_mesh);

	// Function to Interpolate mesh velcoty to the Interior Nodes
	mesh_deform->Solve_mesh_velocity_interior(fespace_mesh,MeshVelo_FEvect);

	cout << " ------------ MESH VELOCITY INTERIOR SOLVED ----------------------" <<endl;

	double time_step = 0.00001;
	mesh_deform->move_mesh_ale(MeshVelo_FEvect,meshvelo_Old,time_step);
	///////////////// ---------- CODE BLOCK ( ALE MESH MOVEMENT ) --------------- /////////////////

	// Function to pick 

	//Delete Function Variables

	//mesh_deform->moveMesh();
	cout << " After  Calling class --------------- Success" << endl;
	os.seekp(std::ios::beg);
    os <<  "VTK/"<<VtkBaseName<<".0000" <<output_write_counter <<".vtk" << ends;
	cout<<" Before VTK "<<  endl;
	Output->WriteVtk(os.str().c_str());
	cout<<" After VTK "<<  endl;
	output_write_counter =  output_write_counter + 1;
	cout << " Program Finised"<<endl;
}
