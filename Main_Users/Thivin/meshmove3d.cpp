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


using namespace std;

// Boundary Condition for all the dimensions of the boundary domain
void BoundCondition(int BdID, double x, double y, double z, BoundCond &cond)
{
	if((fabs(y)<1e-6) || (fabs(y-1)<1e-6) )
	{
		cond = DIRICHLET;
	}
	else 
		cond = NEUMANN;
}


// Boundary Values for all the dimensions of the boundary domain. 
void BoundValue_X(int BdComp, double x, double y, double z, double &value)
{
    value = 0;
}

void BoundValue_Y(int BdComp, double x, double y, double z, double &value)
{
	if((fabs(y-1)<1e-6) )
		value = 10;
	else
  		value = 0;
}

void BoundValue_Z(int BdComp, double x, double y, double z, double &value)
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
    coeff[0] = 1;

    coeff[1] = 0;
    coeff[2] = 0;
    coeff[3] = 0;

    coeff[4] = 0;
  }
}

void Assembly_linear_elasticity_3D(double quad_wt, double *coeff, double *param, double hK, 
                  double **derivatives, int *N_BaseFuncts, double ***LocMatrices, double **LocRhs)
{
  double *N = derivatives[0], *Nx = derivatives[1], *Ny = derivatives[2] , *Nz = derivatives[3];
  double **K11, **K12, **K13, **K21, **K22, **K23 , **K31, **K32, **K33, *F1, *F2 , *F3;
  K11 = LocMatrices[0];
  K12 = LocMatrices[1];
  K13 = LocMatrices[2];
  K21 = LocMatrices[3];
  K22 = LocMatrices[4];
  K23 = LocMatrices[5];
  K31 = LocMatrices[6];
  K32 = LocMatrices[7];
  K33 = LocMatrices[8];
  
  F1 = LocRhs[0];
  F2 = LocRhs[1];
  F3 = LocRhs[1];

  for (int i = 0; i < N_BaseFuncts[0]; i++){
    for (int j = 0; j < N_BaseFuncts[0]; j++){
		K11[i][j] += quad_wt * (Nx[i]*Nx[j] + Ny[i]*Ny[j] + Nz[i]*Nz[j]);
		K22[i][j] += quad_wt * (Nx[i]*Nx[j] + Ny[i]*Ny[j] + Nz[i]*Nz[j]);
		K33[i][j] += quad_wt * (Nx[i]*Nx[j] + Ny[i]*Ny[j] + Nz[i]*Nz[j]);  

      //cout << "Nx[ " << i << "] = " << Nx[i] << endl;
      	K12[i][j] += 0.;
		K13[i][j] += 0.;
      	K21[j][i] += 0.;
		K23[j][i] += 0.;
		K31[j][i] += 0.;
		K32[j][i] += 0.;
		
      /* NON LINEAR PART */
    }
    /* RHS */

    F1[i] = 0.;
    F2[i] = 0.;
	F3[i] = 0.;
  }


}







int main (int argc, char* argv[])
{
	TDatabase *Database = new TDatabase();
	TDomain *Domain = new TDomain(argv[1]);
	TFEDatabase3D *FEDatabase = new TFEDatabase3D();
	TCollection *coll;
	TFESpace3D *fesp[1],*ferhs[3];
	int ORDER;

	// -- Parameters required for VTK   -- //
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

	
    TFESpace3D *fespace = new TFESpace3D(coll, NameString, UString, BoundCondition, ORDER);
	cout << "FESPACE - Declared " << endl;

	int N_DOF = fespace->GetN_DegreesOfFreedom();
	int N_Active = fespace->GetActiveBound();

	cout << endl << "N_Cells = " << coll -> GetN_Cells() << endl << endl;
	cout << "Degrees of Freedom = " << N_DOF << endl << endl << "N_Active = " << N_Active << endl << endl;

	double *sol = new double[3*N_DOF]();
	double *rhs = new double[3*N_DOF]();
	double *gridpos = new double[3*N_DOF]();

	TFEVectFunct3D *vecfunc = new TFEVectFunct3D(fespace, (char*)"C", (char*)"C", gridpos, N_DOF, 3);
	vecfunc -> GridToData();
    cout << " Grid to data - Done "<< endl;

	// Aux Param
	TAuxParam3D *Meshaux = new TAuxParam3D(1, 0, 0, 0, &fespace, NULL, NULL, NULL, NULL, 0, NULL);
	cout << " Thivin " << endl;

	// -------------------------- start of Discrete Form - Equation  ------------------ //

	int N_Terms = 4;   // Number of terms in the All derivatives index  	 
  	int *SpacesNumbers = new int[N_Terms](); // Spacenumbers -- array of values which denotes the SPACE numbers that eac terms mentioned above belongs to
  	int N_Matrices = 9;   // Number of block matrices needed to create our final stiffness matrix , For NSE 2d it will be 4 A matrices and 2 B matrices 	
	int N_RHS = 3;   // The number of components in the RHS of the matrix
	// The row space and column space of the each matrix being constructed
  	//Note:  This is used for mixed finite elements like NSE , where the B matrix needs to compute <grad(u),p>
	//NOte :  Bt and B will have alternate row and column spaces , so these neeeds to be filled accordingly
	int *rowspace = new int[N_Matrices]();
	int *columnspace = new int[N_Matrices]();
	int *rhsspace = new int[N_RHS]();

	MultiIndex3D AllDerivatives[4] = {D000, D100, D010,D001};

	TDiscreteForm3D* discreteform = new TDiscreteForm3D(UString, UString, N_Terms, AllDerivatives,
                                        SpacesNumbers, N_Matrices, N_RHS, rowspace, columnspace, rhsspace,
										                    Assembly_linear_elasticity_3D, GridCoeffs, NULL); 

	cout << " Thivin " << endl;
	// ------------------------- END OF DISCRETE FORM EQUATION --------------------------- //

	// --------------------- START OF MATRIX STRUCTURE DECLARATION -------------------//

	TSquareStructure3D *sqstructure = new TSquareStructure3D(fespace);
	sqstructure -> Sort();
	TSquareMatrix3D *SQMATRICES_GRID[9];
	double *RHS[3];
	double *Entries[9];
	int *GridKCol, *GridRowPtr;
	
	

	TSquareMatrix3D *SqmatrixG11 = new TSquareMatrix3D(sqstructure);
	TSquareMatrix3D *SqmatrixG12 = new TSquareMatrix3D(sqstructure);
	TSquareMatrix3D *SqmatrixG13 = new TSquareMatrix3D(sqstructure);
	TSquareMatrix3D *SqmatrixG21 = new TSquareMatrix3D(sqstructure);
	TSquareMatrix3D *SqmatrixG22 = new TSquareMatrix3D(sqstructure);
	TSquareMatrix3D *SqmatrixG23 = new TSquareMatrix3D(sqstructure);
	TSquareMatrix3D *SqmatrixG31 = new TSquareMatrix3D(sqstructure);
	TSquareMatrix3D *SqmatrixG32 = new TSquareMatrix3D(sqstructure);
	TSquareMatrix3D *SqmatrixG33 = new TSquareMatrix3D(sqstructure);

	SQMATRICES_GRID[0] = SqmatrixG11;
	SQMATRICES_GRID[1] = SqmatrixG12;
	SQMATRICES_GRID[2] = SqmatrixG13;
	SQMATRICES_GRID[3] = SqmatrixG21;
	SQMATRICES_GRID[4] = SqmatrixG22;
	SQMATRICES_GRID[5] = SqmatrixG23;
	SQMATRICES_GRID[6] = SqmatrixG31;
	SQMATRICES_GRID[7] = SqmatrixG32;
	SQMATRICES_GRID[8] = SqmatrixG33;


	RHS[0] = rhs;
    RHS[1] = rhs + N_DOF;
	RHS[2] = rhs + 2*N_DOF;

	Entries[0] = SqmatrixG11->GetEntries();
	Entries[1] = SqmatrixG12->GetEntries();
	Entries[2] = SqmatrixG13->GetEntries();
	Entries[3] = SqmatrixG21->GetEntries();
	Entries[4] = SqmatrixG22->GetEntries();
	Entries[5] = SqmatrixG23->GetEntries();
	Entries[6] = SqmatrixG31->GetEntries();
	Entries[7] = SqmatrixG32->GetEntries();
	Entries[8] = SqmatrixG33->GetEntries();

	GridKCol = sqstructure->GetKCol();
	GridRowPtr = sqstructure->GetRowPtr();


	// ------------------ END OF MATRIX STRUCURE DECLARATIONS ----------------//

	// --------------- START OF BOUNDARY FUNCTIONS -------------------//
	//Boundary Conditions
	BoundCondFunct3D *GridBoundaryConditions[3];
	BoundValueFunct3D *GridBoundValues[3];

	GridBoundaryConditions[0] = BoundCondition;
	GridBoundaryConditions[1] = BoundCondition;
	GridBoundaryConditions[2] = BoundCondition;
	GridBoundValues[0] = BoundValue_X;
	GridBoundValues[1] = BoundValue_Y;
	GridBoundValues[2] = BoundValue_Z;
	
	// --------------- END OF BOUNDARY FUNCTIONS -------------------//

	
	
	// ---------------- START OF ASSEMBLY 3D FUNCTION -----------//
	fesp[0] = fespace;
    ferhs[0] = fespace;
    ferhs[1] = fespace;
	ferhs[2] = fespace;
	
	TDiscreteForm3D *DiscreteFormGrid;
	InitializeDiscreteFormGrid(DiscreteFormGrid, GridCoeffs);
	TAssembleMat3D *MeshMatAssemble = new TAssembleMat3D(1, &fespace, 9, SQMATRICES_GRID, 0, NULL, 0, NULL, NULL, DiscreteFormGrid, GridBoundaryConditions, GridBoundValues, Meshaux);
	MeshMatAssemble->Init();
	MeshMatAssemble->Reset();
	MeshMatAssemble->Assemble3D();


	

	DirectSolver(SQMATRICES_GRID, 3, 3, sol, rhs);


	

	TFEVectFunct3D *displacements = new TFEVectFunct3D(fespace, (char*)"C", (char*)"C", sol, N_DOF, 3);
	// ------ SOLVER ----------------------------- //
	//DirectSolver(SQMATRICES_GRID, 3, 3, sol, rhs);
	// Copy solution of the system to the grid pos array , to be added to the existing co-ordinates
	for ( int i=0; i< 3*N_DOF ; i++){
    	gridpos[i] += sol[i];
		//cout << "Sol[" << i << "] : " << sol[i] << endl;
		}
	// Transfer the new grid position as current grid co ordinates in the system
	vecfunc -> DataToGrid();
	
	for ( int i = 0 ; i < N_Matrices ; i++)
		SQMATRICES_GRID[i]->Reset();

    for (int i_rhs = 0; i_rhs < 3*N_DOF; i_rhs++)
      rhs[i_rhs] = 0;



	// const char vtkdir[] = "VTK";
	char *VtkBaseName;
	// char UString[] = "u";
	// char PString[] = "p";
	// // std::ostringstream os;
	// os << " ";

	VtkBaseName = TDatabase::ParamDB->VTKBASENAME;
  	TOutput3D* Output = new TOutput3D(2, 2, 1, 1, Domain);
	Output->AddFEVectFunct(vecfunc);
	os.seekp(std::ios::beg);
	int output_write_counter = 0;
    os <<  "VTK/"<<VtkBaseName<<".0000" <<output_write_counter <<".vtk" << ends;
	cout<<" hi "<<  endl;
	Output->WriteVtk(os.str().c_str());
	output_write_counter =  output_write_counter + 1;
		
	return 0;

}
