/** ************************************************************************ 
 @brief     source file for deformMesh3D 
             Moves the computational domain based on the boundary conditions provided by solving 2D linear Elasticity equation

             Parent Class : <same> deformMesh3D

             Parameters Required for Constructors :
              fespace ( pointer of TFESpace2D)
              BoundCondition_x ( function pointer of boundary condition of X)
              BoundCondition_y ( function pointer of boundary condition of Y)
              BoundValue_x     ( function pointer of boundary value of X)
              BoundValue_y     ( function pointer of boundary value of Y)


* @author    Pon Suganth Elangovan , Thivin Anandh
* @date      4-Sep-2019
* @History   
 ************************************************************************  */

#include <Domain.h>
#include <Database.h>
#include <FEDatabase3D.h>
#include <SquareMatrix3D.h>
#include <SquareStructure3D.h>
#include <DiscreteForm3D.h>
#include <BoundFace.h>
#include <BoundComp3D.h>
#include <Assemble3D.h>
#include <AssembleMat3D.h>
#include <DirectSolver.h>
#include <FE3D.h>
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
#include <algorithm>
#include <LinAlg.h>
#include <list>
#include <map>
#include <tuple>

#include <FEDatabase3D.h>
#include <BoundFace.h>
#include <IsoJointEqN.h>
#include <HexaAffin.h>
#include <HexaTrilinear.h>
#include <HexaIsoparametric.h>
#include <TetraAffin.h>
#include <TetraIsoparametric.h>
#include <FEFunction3D.h>
#include <InterfaceJoint3D.h>
#include <NodalFunctional3D.h>



#include "DeformMesh3D.h"


void deformMesh3D::GridCoeffs(int n_points, double *x, double *y,double *z, double **parameters, double **coeffs)
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

void deformMesh3D::Assembly_poisson_3D(double quad_wt, double *coeff, double *param, double hK, 
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

void deformMesh3D::moveMesh()
{
	char UString[] = "T";
	char NameString[] = "name";
	char CString[] = "C";
    int ORDER  = 1;

	

	// Initialise all the FE3D elements to Tetra

    TFESpace3D* fespace = new TFESpace3D(coll, NameString, UString, GridBoundaryConditions[0], ORDER);

	int N_DOF = fespace->GetN_DegreesOfFreedom();
	int N_Active = fespace->GetActiveBound();

	cout  << "N_Cells = " << coll -> GetN_Cells() << endl;
	cout << "Degrees of Freedom = " << N_DOF  << "    N_Active = " << N_Active << endl;

	double *sol = new double[3*N_DOF]();
	double *rhs = new double[3*N_DOF]();
	double *gridpos = new double[3*N_DOF]();

	TFEVectFunct3D *vecfunc = new TFEVectFunct3D(fespace, (char*)"C", (char*)"C", gridpos, N_DOF, 3);
	vecfunc -> GridToData();

	// Aux Param
	TAuxParam3D *Meshaux = new TAuxParam3D(1, 0, 0, 0, &fespace, NULL, NULL, NULL, NULL, 0, NULL);

	// -------------------------- start of Discrete Form - Equation  ------------------ //

	int N_Terms = 4;   // Number of terms in the All derivatives index  	 
  	int *SpacesNumbers = new int[N_Terms](); // Spacenumbers -- array of values which denotes the SPACE numbers that eac terms mentioned above belongs to
  	int N_Matrices = 9;   // Number of block matrices needed to create our final stiffness matrix , For NSE 2d it will be 4 A matrices and 2 B matrices 	
	int N_RHS = 3;   // The number of components in the RHS of the matrix
	// The row space and column space of the each matrix being constructed
  	//Note:  This is used for mixed finite elements like NSE , where the B matrix needs to compute <grad(u),p>
	//NOte :  Bt and B will have alternate row and column spaces , so these neeeds to be filled accordingly
	int *rowspace = new int[N_Matrices](); // we are initialising it as zero because all blocks belong to the same FE space
	int *columnspace = new int[N_Matrices](); // we are initialising it as zero because all blocks belong to the same FE space
	int *rhsspace = new int[N_RHS](); // we are initialising it as zero because all blocks belong to the same FE space

	MultiIndex3D AllDerivatives[4] = {D000, D100, D010,D001};

	TDiscreteForm3D* discreteform = new TDiscreteForm3D(UString, UString, N_Terms, AllDerivatives,
                                        SpacesNumbers, N_Matrices, N_RHS, rowspace, columnspace, rhsspace,
										                    Assembly_poisson_3D, GridCoeffs, NULL); 
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


	// ---------------- START OF ASSEMBLY 3D FUNCTION -----------//
	fesp[0] = fespace;   // Type of FE Space to be used for Blocks in A Matrix
    ferhs[0] = fespace;  // Type of FE Space to be used for Blocks in Rhs Matrix
    ferhs[1] = fespace;  // Type of FE Space to be used for Blocks in Rhs Matrix
	ferhs[2] = fespace;  // Type of FE Space to be used for Blocks in Rhs Matrix
	
	// Get an Instance for 3D Matrix Assembly
	TAssembleMat3D *MeshMatAssemble = new TAssembleMat3D(1, &fespace, 9, SQMATRICES_GRID, 0, NULL, 3, RHS, ferhs, discreteform, GridBoundaryConditions, GridBoundValues, Meshaux);
	MeshMatAssemble->Init();
	MeshMatAssemble->Reset();
	MeshMatAssemble->Assemble3D();

	memset(Entries[1] + GridRowPtr[N_Active], 0, (GridRowPtr[N_DOF] - GridRowPtr[N_Active])*SizeOfDouble);
	memset(Entries[2] + GridRowPtr[N_Active], 0, (GridRowPtr[N_DOF] - GridRowPtr[N_Active])*SizeOfDouble);
	memset(Entries[3] + GridRowPtr[N_Active], 0, (GridRowPtr[N_DOF] - GridRowPtr[N_Active])*SizeOfDouble);
	memset(Entries[5] + GridRowPtr[N_Active], 0, (GridRowPtr[N_DOF] - GridRowPtr[N_Active])*SizeOfDouble);
	memset(Entries[6] + GridRowPtr[N_Active], 0, (GridRowPtr[N_DOF] - GridRowPtr[N_Active])*SizeOfDouble);
	memset(Entries[7] + GridRowPtr[N_Active], 0, (GridRowPtr[N_DOF] - GridRowPtr[N_Active])*SizeOfDouble);



	// N_Active = Non Dirichlet DOF's
	int N_BDDof = N_DOF - N_Active;

	// ------ Start of SOLVER ----------------------------- // 
    PardisoDirectSolver(SQMATRICES_GRID, 3, 3, sol, rhs);

	// ------ SOLVER ----------------------------- //
	for ( int i=0; i< 3*N_DOF ; i++)
    	gridpos[i] += sol[i];
	
	cout << " SOL NORM : " << Ddot(3*N_DOF,sol,sol)<<endl;
	
	// Transfer the new grid position as current grid co ordinates in the system
	vecfunc -> DataToGrid();
	

	for ( int i = 0 ; i < N_Matrices ; i++)
		SQMATRICES_GRID[i]->Reset();

	cout << " CLASS -  DeformMesh 3D : --------  Mesh has been Moved ---------------- " <<endl;

}


// Constructor
deformMesh3D::deformMesh3D(TCollection* _coll,	BoundCondFunct3D **_GridBoundaryConditions, BoundValueFunct3D **_GridBoundValues)
{
    GridBoundaryConditions  = _GridBoundaryConditions;
    GridBoundValues  = _GridBoundValues;
	coll              	= _coll;
    
}

void deformMesh3D::get_surface_normals(TFEVectFunct3D* MeshVelo_FEvect,TFEVectFunct3D* Velocity_FEvect )
{

	//variable Declarations
	cout << " REACHED - Get surface Normals " <<endl;
	TBaseCell* currentCell;  TVertex* currentVertex;	
	FE3D FEId, FEId_velocity; TFE3D *ele , *ele_Velocity;
	BF3DRefElements RefElement , RefElement_Velocity;	
	RefTrans3D RefTrans , RefTransVelocity;TRefTrans3D *F_K , *F_K_V;
	TBaseFunct3D *bf;
	TFESpace3D *Velocity_FESpace,*Mesh_FESpace;
	TFEFunction3D*  Velocity_FEFunction[3];


	//  Get the Compunents of Velocity_FEvect and MeshVelo_fevect
	int meshVelComp = MeshVelo_FEvect->GetN_Components();
	int velComp = Velocity_FEvect->GetN_Components();

	//get FE Space for Velocity & Mesh
	Velocity_FESpace = Velocity_FEvect->GetFESpace3D();
	Mesh_FESpace = MeshVelo_FEvect->GetFESpace3D();	

	
	double* VelocityArray = Velocity_FEvect->GetComponent(0)->GetValues();
	double* MeshVelocityArray = MeshVelo_FEvect->GetComponent(0)->GetValues();

	// Create an FEFunction array for all the the FEVectfunction3d objects and to Store Value Arrays from the respective FEFunctions
	for ( int i = 0 ; i < velComp ; i++)
	{   
		Velocity_FEFunction[i] = Velocity_FEvect->GetComponent(i);
	}
	


	// Get the Lengths of the Value Arrays of Respective FEFunction Values
	int MeshVeloLength = MeshVelo_FEvect->GetComponent(0)->GetLength();
	int VelocityLength = MeshVelo_FEvect->GetComponent(0)->GetLength();

	cout << " &&&&&&&&&&&&7Velocity Length : " << VelocityLength <<endl;
	cout << " TEST SIZE : " << 3*VelocityLength - 1<<endl;
	cout << " TEST : " << VelocityArray[3*VelocityLength - 1]  << endl;
	for ( int i = 0; i < VelocityLength ; i++)
	{
		VelocityArray[i                     ] = i;
		VelocityArray[i + VelocityLength    ] = i + 100;
		VelocityArray[i + (VelocityLength*2)] = i + 200;
		
	}

	cout << " Velocity Array Assigned " <<endl;

	double norm1 = 0.,norm2 = 0.,norm3 = 0.,normLen = 0.;
	
	xi_1.resize(N_freeSurfaceVertex,0); // Arrays to store temp XI valu on 2d Face
	xi_2.resize(N_freeSurfaceVertex,0); // Arrays to store temp XI valu on 2d Face 

	xi_freeSurfaceNodes.resize(N_freeSurfaceVertex,0); // Arrays to store temp XI valu on 2d Face
	eta_freeSurfaceNodes.resize(N_freeSurfaceVertex,0); // Arrays to store temp XI valu on 2d Face
	zeta_freeSurfaceNodes.resize(N_freeSurfaceVertex,0); // Arrays to store temp XI valu on 2d Face

	// Resize the normal Vectors based on the number of free surface vertex available
	meshVelocityatFreeSurface_1.resize(N_freeSurfaceVertex,0);
	meshVelocityatFreeSurface_2.resize(N_freeSurfaceVertex,0);
	meshVelocityatFreeSurface_3.resize(N_freeSurfaceVertex,0);

	//Calculate normal for Each Node in the Joint. 
	// Hard Code the Values of Joint id's for the given HEXAHEADRAL - TRILINEAR
	for ( int i = 0 ; i < N_freeSurfaceVertex ; i ++)
	{
		int vertex_number = i;
		cout << " --------------------------- vertex number " <<i <<"---------------------- " <<endl;
		cout << " Local DOF of the Vertex : " << freeSurfaceVertexLocal[i] ;
		cout  << " Global DOF of the Vertex : " << freeSurfaceVertex[i] <<endl;
		int cellNr = freeSurfaceCells[i];
		currentCell =  coll->GetCell(cellNr);

		// Parameters for the Mesh Velocity FESPACE	---- /////
		FEId = Mesh_FESpace->GetFE3D(cellNr, currentCell);
		ele = TFEDatabase3D::GetFE3D(FEId);
		RefElement = TFEDatabase3D::GetRefElementFromFE3D(FEId);
		RefTrans3D referenceTransformation =  TFEDatabase3D::GetRefTrans3D_IDFromFE3D(FEId);
		// --END-- Parameters for the Mesh Velocity FESPACE	---- /////


		// Parameters for the  Velocity FESPACE	---- /////
		FEId_velocity = Velocity_FESpace->GetFE3D(cellNr, currentCell);
		ele_Velocity = TFEDatabase3D::GetFE3D(FEId_velocity);
		RefElement_Velocity = TFEDatabase3D::GetRefElementFromFE3D(FEId_velocity);
		RefTrans3D RefTransVelocity =  TFEDatabase3D::GetRefTrans3D_IDFromFE3D(FEId_velocity);
		// --END-- Parameters for the Mesh Velocity FESPACE	---- /////

		//----------- Get the Local to Global Numbering of the Velocity FE Space and MESH Fe Space  ----  //
		int* GlobalNumbers_Velocity = Velocity_FESpace->GetGlobalNumbers();
		int* BeginIndex_velocity = Velocity_FESpace->GetBeginIndex();
		int* Numbers_Velocity  = GlobalNumbers_Velocity + BeginIndex_velocity[cellNr];

		int* GlobalNumbers_Mesh = Mesh_FESpace->GetGlobalNumbers();
		int* BeginIndex_Mesh = Velocity_FESpace->GetBeginIndex();
		int* Numbers_Mesh = GlobalNumbers_Mesh + BeginIndex_Mesh[cellNr];
		//----END---- Get the Local to Global Numbering of the Velocity FE Space and MESH Fe Space  ----  //

		// SETUP the Values of s and T that needs to be sent to "GetOuterNormal"/ getTangentVectors to get the normal at the point
		double x_coord ,y_coord,z_coord;
		double xi_temp , eta_temp , zeta_temp;
		Mesh_FESpace->GetDOFPosition(freeSurfaceVertex[i], x_coord,y_coord, z_coord);
		// cout << "X : "<< x_coord <<"  Y : "<< y_coord <<"   Z : "<< z_coord <<endl;
		int JointNumber = freeSurfaceJoints[i];
		
		

		////////////// CODE BLOCK - A - One Time Assignment for finding the xi eta zeta of the Free surface Nodes - ///////////////
		// NOte : Move this Code block to a newer Function so that it is executed only once ////////////////////
		// *****IMPORTANT*** Note 2 : This Code Block is Executed only for mesh Velocity Space , not on Velocity FE Space.

		switch(referenceTransformation)     // Reftrans of Velocity
    	{
			case HexaTrilinear: // HEXATRILINEAR 
			{	
				//RefTrans = HexaTrilinear;
				F_K = TFEDatabase3D::GetRefTrans3D(referenceTransformation);
				((THexaTrilinear *)F_K)->SetCell(currentCell);

				int localNodeNUmber = freeSurfaceVertexLocal[vertex_number];

				switch ( localNodeNUmber)
				{
					case 0:
					{
						xi_freeSurfaceNodes[i] = -1; eta_freeSurfaceNodes[i] = -1; zeta_freeSurfaceNodes[i] = -1;
						break;
					}
					case 1:
					{
						xi_freeSurfaceNodes[i] = 1; eta_freeSurfaceNodes[i] = -1; zeta_freeSurfaceNodes[i] = -1;
						break;
					}
					case 2:
					{
						xi_freeSurfaceNodes[i] = 1; eta_freeSurfaceNodes[i] = 1; zeta_freeSurfaceNodes[i] = -1;
						break;
					}
					case 3:
					{
						xi_freeSurfaceNodes[i] = -1; eta_freeSurfaceNodes[i] = 1; zeta_freeSurfaceNodes[i] = -1;
						break;
					}
					case 4:
					{
						xi_freeSurfaceNodes[i] = -1; eta_freeSurfaceNodes[i] = -1; zeta_freeSurfaceNodes[i] = 1;
						break;
					}
					case 5:
					{
						xi_freeSurfaceNodes[i] = 1; eta_freeSurfaceNodes[i] = -1; zeta_freeSurfaceNodes[i] = 1;
						break;
					}
					case 6:
					{
						xi_freeSurfaceNodes[i] = 1; eta_freeSurfaceNodes[i] = 1; zeta_freeSurfaceNodes[i] = 1;
						break;
					}
					case 7:
					{
						xi_freeSurfaceNodes[i] = -1; eta_freeSurfaceNodes[i] = 1; zeta_freeSurfaceNodes[i] = 1;
						break;
					}
					default:
						cout << " ERROR IN CLASS - DEFORMMESH3D - Function - Pick Surface Displacements " <<endl;
						cout << " The Local node NUmber " << localNodeNUmber << " Exceeds the Max Value of 8 for Hexaheadron " <<endl;
						exit(0); /// EXIT SEQUENCE

						
				}


				switch (JointNumber)
				{
					case 0:
					{
						xi_1[i] = xi_freeSurfaceNodes[i]; xi_2[i] = eta_freeSurfaceNodes[i]; break;
					}
					case 1:
					{
						xi_1[i] = zeta_freeSurfaceNodes[i]; xi_2[i] = xi_freeSurfaceNodes[i]; break;
					} 
					case 2:
					{
						xi_1[i] = zeta_freeSurfaceNodes[i]; xi_2[i] = eta_freeSurfaceNodes[i]; break;
					}
					case 3:
					{
						xi_1[i] = zeta_freeSurfaceNodes[i]; xi_2[i] = -xi_freeSurfaceNodes[i]; break;
					}
					case 4:
					{
						xi_1[i] = eta_freeSurfaceNodes[i]; xi_2[i] = zeta_freeSurfaceNodes[i];break;
					}
					case 5:
					{
						xi_1[i] = eta_freeSurfaceNodes[i]; xi_2[i] = xi_freeSurfaceNodes[i];break;
					}
					default:
					{
						cout << " Error in Class : DeformMesh3d - Function : get_surface_Displacements"<<endl;
						cout << " The Given Joint Id for the HEXAHEADRAL cell "<< freeSurfaceJoints[i] << " Does not exist "<<endl;
					}
				}

				break;

			}
			default:
			{
				cout << " ERROR TYPE : Method Not yet Implemented "<<endl;
				cout << " ERROR DESC : Method to Calculate the Normal has not been Implemented for REFRTANS3D type : "<< referenceTransformation <<endl;
				cout << " ERROR LOCATION : Class DeformMesh3D , Function get_surface_normals " <<endl;
				exit(0);
			}
		}
		////////////// END - CODE BLOCK - A - One Time Assignment for finding the xi eta zeta of the Free surface Nodes - ///////////////
		
		
		
		//////////// CODE BLOCK A2 - Calculate Joint Normal using "MESH" velocity FE Space  ////////////////////
		////// Note : This Code Block is Repeated Block of Code to calculate Normal of a Face at the mesh Face //////
		////// Note : This Code is for "MESH" Velocity FE Space only /////////////////
						
		switch(referenceTransformation)     // Reftrans of MESH Velocity
    	{
			case HexaTrilinear: // HEXATRILINEAR 
			{	
				//RefTrans = HexaTrilinear;
				F_K = TFEDatabase3D::GetRefTrans3D(referenceTransformation);
				((THexaTrilinear *)F_K)->SetCell(currentCell);
				( (THexaTrilinear *) F_K )->GetOuterNormal(JointNumber,xi_1[i],xi_2[i],
											norm1, norm2, norm3);
				//cout << "Xi : "<< xi_temp <<"  eta : "<< eta_temp <<" Zeta : "<< zeta_temp <<endl;
				break;
			}
			default:
			{
				cout << " ERROR TYPE : Method Not yet Implemented "<<endl;
				cout << " ERROR DESC : Method to Calculate the Normal has not been Implemented for REFRTANS3D type : "<< referenceTransformation <<endl;
				cout << " ERROR LOCATION : Class DeformMesh3D , Function get_surface_normals " <<endl;
				exit(0);
			}
		}

		normLen = sqrt(norm1*norm1 + norm2*norm2 + norm3*norm3);
		double norm_1  = norm1/normLen; 
		double norm_2 = norm2/normLen; 
		double norm_3 = norm3/normLen; 

		// cout << " Norm Vector  : [ " << norm_1 << ", " << norm_2 << ", " << norm_3 << " ]  Norm len : "  << normLen <<endl;

	//////////// -END- CODE BLOCK A2 - Calculate Joint Normal using "MESH" velocity FE Space  ////////////////////






	//////////// - CODE BLOCK B - To get the Velocity Values at the nodal Points ///////////////////////
	// Note - This Code Block will be repeated Every Time when the meshmovement is called /////
	// This Part is Called for "Velocity FE Space"  to Obtain the values of Velocity at nodal Points

		double* u_values = new double[3]();
		//Get the Global Co-ordinates of the DOF in Velocity FESpace
		TBaseCell* VelocityCell =  MovCells[i];
		int freeSurface_globalID = freeSurfaceVertex[i];
		switch(RefTransVelocity)
    	{
			case HexaTrilinear: 
			{
				F_K_V = TFEDatabase3D::GetRefTrans3D(RefTransVelocity);		
				((THexaTrilinear *)F_K)->SetCell(currentCell);
				// ------ Calculate the Values of the Velocity Function at the Nodal Points   ------- //
				bf = ele_Velocity->GetBaseFunct3D();
  				int N_BaseFunct = bf->GetDimension();
				int BaseVectDim = bf->GetBaseVectDim();
				// Declare Variables and Spaces for Arrays
				double* uorig = new double[N_BaseFunct*BaseVectDim];
				
				bf->GetDerivatives(D000, xi_freeSurfaceNodes[i], eta_freeSurfaceNodes[i], zeta_freeSurfaceNodes[i], uorig);

				// For D000 , the u_ref Value is the u Orig Value ( Refer GetOrigValues function in THexaTrilinear)
				double val_1 = 0. ;
				// cout<< "(xi,eta,zeta) :  ( "<< xi_freeSurfaceNodes[i] <<", "<<eta_freeSurfaceNodes[i]<<", " <<zeta_freeSurfaceNodes[i] <<" ) " <<endl;
				double* VelocityArray_Component;
				for (int i = 0 ; i < 3 ; i++){
					VelocityArray_Component = VelocityArray + i*VelocityLength;
					for(int j=0;j<N_BaseFunct;j++)
					{
						val_1 = VelocityArray_Component[Numbers_Velocity[j]];
						// cout << " VEL VAL : (" << val_1 << ", " << uorig[j] << ")  ";
						u_values[i] +=  uorig[j]*val_1;				
					}
					// cout<< " uVal["<<i<<"] : "<<u_values[i]<< "   " <<endl;
				}
				cout<<endl;
				delete[] uorig;
				break;
			}
			default:
			{
				cout << " ERROR TYPE : Method Not yet Implemented "<<endl;
				cout << " ERROR DESC : Method to Calculate the Normal has not been Implemented for REFRTANS3D type : "<< referenceTransformation <<endl;
				cout << " ERROR LOCATION : Class DeformMesh3D , Function: get_surface_normals " <<endl;
				exit(0);
			}
		}


	
		//////////// - END -- CODE BLOCK B - To get the Normal Vectors and Velocity Values at the nodal Points ///////////////////////
		//cout << " Velocity Vector  : [ " << u_values[0] << ", " << u_values[1] << ", " << u_values[2] << " ] " <<endl;
		double alpha = 0;
		alpha = u_values[0]*norm_1 + u_values[1]*norm_2 + u_values[2]*norm_3;	
		// cout << " Alpha : " << alpha <<endl;
		double* MeshVelocity_ValuesArray = MeshVelocityArray;
		
	
		// ASsign the ith component of Mesh Velocity as ( u1.n1 + u2.n2 + u3.n3 ) * Normal
        // HARDCODED - For 
		meshVelocityatFreeSurface_1[vertex_number] = alpha*norm_1;   // Value of w1 of Global DOF freeSurfaceVertex[i]
		meshVelocityatFreeSurface_2[vertex_number] = alpha*norm_2;   // Value of w2 of Global DOF freeSurfaceVertex[i]
		meshVelocityatFreeSurface_3[vertex_number] = alpha*norm_3;   // Value of w3 of Global DOF freeSurfaceVertex[i]


		delete[] u_values;
	}


	for ( int i = 0 ; i < N_freeSurfaceVertex ; i++)
			cout << " Global No : " << 	freeSurfaceVertex[i] << "  Val 1 : "<< meshVelocityatFreeSurface_1[i] << 
			" val 2 : "<< meshVelocityatFreeSurface_2[i] 
			 << " Val 3 : "<< meshVelocityatFreeSurface_3[i] << endl;
 	cout << " Size of freeSurfaceCells : " << freeSurfaceCells.size()<<endl;
	cout << " Size of freeSurfaceJoints: " << freeSurfaceJoints.size()<<endl;
	cout << " Size of freeSurfaceVertex : " << freeSurfaceVertex.size()<<endl;
}


//
void deformMesh3D::Pick_free_surface_DOFs(TFESpace3D* _fespace, std::vector<int> BoundIds,TCollection* coll)
{

	N_freeSurfaceJoints = 0;
	N_freeSurfaceCells = 0;
	N_freeSurfaceVertex = 0 ;
	fespace_alemeshMovement = _fespace;
	
	// -------------------- Variable Declaraions -------------------------//
		
    int* N_Movfaces; 					// Local Face id's of the boundary Domain 
	TBaseCell* currentCell;   			// Cell pointer for the current cell
	TJoint* Joint;						// Pointer to Joint in 3D Cell
	TBoundFace* Bdface;					// Pointer to Boundary Face in 3D Cell
	TBoundComp *BoundComp;
	TVertex* currentVertex;


	int N_Cells = coll->GetN_Cells();   // N_cells in given FE system
	int N_Vertices; 
	int N_Joints;
	int TMP_DOF = 0;
	const int *TmpFV,*TmpLen;
	int MaxLen;
	int *GlobalDOF;
	int *BeginIndex;
	int* GlobalNumbers;

	GlobalNumbers = fespace_alemeshMovement->GetGlobalNumbers();
	BeginIndex = fespace_alemeshMovement->GetBeginIndex();	

	// for ( int k = 0 ; k < N_Cells*8 ; k++)
	// {
	// 	cout << "loc : " << k <<  "  Glob  " <<   GlobalNumbers[k] << "  " <<endl;
 	// }

	// Loop over all the cells to collect the DATA VALUES
	for( int cellId = 0 ; cellId < N_Cells ; cellId++)
	{
		currentCell = coll->GetCell(cellId); 	    // Obatin the pointer to current cell from Collection
		currentCell->SetGlobalCellNo(cellId);
		GlobalDOF = GlobalNumbers + BeginIndex[cellId];
								
		FE3D elementId = fespace_alemeshMovement->GetFE3D(cellId, currentCell);
        TFE3D *element = TFEDatabase3D::GetFE3D(elementId);
        TFEDesc3D *fedesc = element->GetFEDesc3D();
	
		N_Joints = currentCell->GetN_Joints();
		bool cell_setflag = FALSE;

		currentCell->GetShapeDesc()->GetFaceVertex(TmpFV, TmpLen, MaxLen);

		// Obtain the Vertex and face information of the current cell
		// Loop Over the Joints in the cell 
		for ( int jointId = 0 ; jointId < N_Joints ; jointId++)  // Joints refer to Faces in 3D
		{
			Joint = currentCell->GetJoint(jointId);
			if(Joint->GetType() == BoundaryFace)
			{
				// cout << "  Boundary Face " <<endl;
				Bdface = (TBoundFace*)Joint;
	  			BoundComp = Bdface->GetBoundComp();
	  			int bdid = BoundComp->GetID();
				std::vector<int>::iterator it = std::find(BoundIds.begin(), BoundIds.end(), bdid);
				if(it != BoundIds.end())                     // CHANGE THIS HARDCODED VALUE
				{
					N_freeSurfaceJoints++;
					if(cell_setflag == FALSE){
						N_freeSurfaceCells++;					
						cell_setflag = TRUE ;
					}   
					int *JointDOF = fedesc->GetJointDOF(jointId);
					N_Vertices = TmpLen[jointId];
					N_Movfaces++;
					//freeSurfaceJoints.emplace_back(jointId);
					for ( int vert = 0 ; vert < fedesc->GetN_JointDOF() ; vert++)
					{
						int local_vertex_no = TmpFV[jointId*MaxLen + vert];
						currentVertex = currentCell->GetVertex(vert) ;
						int glob_vertex_no = GlobalDOF[JointDOF[vert]];
						if (std::find(freeSurfaceVertex.begin(), freeSurfaceVertex.end(), glob_vertex_no) == freeSurfaceVertex.end()) {
							freeSurfaceVertex.push_back(glob_vertex_no);
							freeSurfaceVertexLocal.emplace_back(local_vertex_no);
							freeSurfaceCells.emplace_back(cellId);
							freeSurfaceJoints.emplace_back(jointId);
							N_freeSurfaceVertex++;
						}
						if(currentVertex->GetClipBoard() != -5 || N_freeSurfaceVertex == 0)
						{						
							currentVertex->SetClipBoard(-5);
							
						}
					}
				}
			}
		}
	}


	// NOW Remove Duplicate GLOBAL DOF's in the Given System 
	
	// cout << " ----- GLOBAL DOF's After remove " <<endl;
	// for ( int k = 0 ; k < freeSurfaceVertex.size() ; k++)
	// 	cout << " " << freeSurfaceVertex[k]  << "  Joint :  "<<freeSurfaceJoints[k] << "  Cells : "  << freeSurfaceCells[k] << endl;
	
	// cout << " ----------------------------" <<endl;

 

	// Collect all the Cells and Save them in Tvertex** POinter
	MovBoundVert = new TVertex*[N_freeSurfaceVertex];
	MovCells = new TBaseCell*[N_freeSurfaceVertex];
	MovJoints = new TJoint*[N_freeSurfaceJoints];

	for( int i = 0  ; i < N_freeSurfaceVertex ; i++){
		MovCells[i] = coll->GetCell(freeSurfaceCells[i]);
	}

	cout << " No of Free Surface DOF after DOF picking: " << N_freeSurfaceVertex<<endl;
}




void deformMesh3D::Solve_mesh_velocity_interior(TFESpace3D* Meshfespace, TFEVectFunct3D* MeshVelocityVectFunction3D)
{
	char UString[] = "T";
	char NameString[] = "name";
	char CString[] = "C";

    int N_DOF = Meshfespace->GetN_DegreesOfFreedom();
	int N_Active = Meshfespace->GetActiveBound();
    
    cout  << "N_Cells = " << coll -> GetN_Cells() << endl;
	cout << "Degrees of Freedom = " << N_DOF  << "    N_Active = " << N_Active << endl;
	double *sol = new double[3*N_DOF]();
	double *rhs =  new double[3*N_DOF]();
	
	double* meshVelocityValues;
	meshVelocityValues = MeshVelocityVectFunction3D->GetValues();
    
    TAuxParam3D *Meshaux = new TAuxParam3D(1, 0, 0, 0, &Meshfespace, NULL, NULL, NULL, NULL, 0, NULL);
    
    int N_Terms = 4;    
  	int *SpacesNumbers = new int[N_Terms](); 
  	int N_Matrices = 9;    	
	int N_RHS = 3;   
	int *rowspace = new int[N_Matrices]();
	int *columnspace = new int[N_Matrices](); 
	int *rhsspace = new int[N_RHS]();
    
    
    MultiIndex3D AllDerivatives[4] = {D000, D100, D010,D001};
     
    TDiscreteForm3D* discreteform = new TDiscreteForm3D(UString, UString, N_Terms, AllDerivatives,
                                        SpacesNumbers, N_Matrices, N_RHS, rowspace, columnspace, rhsspace,
										                    Assembly_poisson_3D, GridCoeffs, NULL); 
    
    // --------------------- START OF MATRIX STRUCTURE DECLARATION -------------------//

	TSquareStructure3D *sqstructure = new TSquareStructure3D(Meshfespace);
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

	SQMATRICES_GRID[0] = SqmatrixG11;  SQMATRICES_GRID[1] = SqmatrixG12; SQMATRICES_GRID[2] = SqmatrixG13;
	SQMATRICES_GRID[3] = SqmatrixG21;  SQMATRICES_GRID[4] = SqmatrixG22; SQMATRICES_GRID[5] = SqmatrixG23;
	SQMATRICES_GRID[6] = SqmatrixG31;  SQMATRICES_GRID[7] = SqmatrixG32; SQMATRICES_GRID[8] = SqmatrixG33;

	RHS[0] = rhs; RHS[1] = rhs + N_DOF; RHS[2] = rhs + 2*N_DOF;

	Entries[0] = SqmatrixG11->GetEntries(); Entries[1] = SqmatrixG12->GetEntries(); Entries[2] = SqmatrixG13->GetEntries();
	Entries[3] = SqmatrixG21->GetEntries(); Entries[4] = SqmatrixG22->GetEntries(); Entries[5] = SqmatrixG23->GetEntries();
	Entries[6] = SqmatrixG31->GetEntries(); Entries[7] = SqmatrixG32->GetEntries(); Entries[8] = SqmatrixG33->GetEntries();
	
	GridKCol = sqstructure->GetKCol(); GridRowPtr = sqstructure->GetRowPtr();
	
	// ---------------- START OF ASSEMBLY 3D FUNCTION -----------//
	fesp[0] = Meshfespace;   // Type of FE Space to be used for Blocks in A Matrix
    ferhs[0] = Meshfespace;  // Type of FE Space to be used for Blocks in Rhs Matrix
    ferhs[1] = Meshfespace;  // Type of FE Space to be used for Blocks in Rhs Matrix
	ferhs[2] = Meshfespace;  // Type of FE Space to be used for Blocks in Rhs Matrix
	
    TAssembleMat3D *MeshMatAssemble = new TAssembleMat3D(1, &Meshfespace, 9, SQMATRICES_GRID, 0, NULL, 3, RHS, ferhs, discreteform, GridBoundaryConditions, GridBoundValues, Meshaux);
	MeshMatAssemble->Init();
	MeshMatAssemble->Reset();
	MeshMatAssemble->Assemble3D();

	// N_Active = Non Dirichlet DOF's
	int N_BDDof = N_DOF - N_Active;

	// Transfer the Solution Calculated from "U" as a Dirichlet boundary to the Mesh Velocity Solution
	for ( int i = 0 ; i < N_freeSurfaceVertex ; i++)
	{
		int globDOF = freeSurfaceVertex[i];
		rhs[globDOF           ] =  meshVelocityatFreeSurface_1[i];
		rhs[globDOF +  1*N_DOF] =  meshVelocityatFreeSurface_2[i];
		rhs[globDOF +  2*N_DOF] =  meshVelocityatFreeSurface_3[i];
	}

	auto start =  clock();
	DirectSolver(SQMATRICES_GRID, 3, 3, sol, rhs);
	auto end =  clock();
    auto duration = (end -start)/double(CLOCKS_PER_SEC)*1000; 
	
	memcpy( meshVelocityValues, sol , SizeOfDouble*3*N_DOF);

	// Release the Matrix Storage Parameters
	for ( int i = 0 ; i < N_Matrices ; i++)
		SQMATRICES_GRID[i]->Reset();

    for (int i_rhs = 0; i_rhs < 3*N_DOF; i_rhs++)
      rhs[i_rhs] = 0;


	std::cout << " Time Taken to Solve : "<<duration << " ms"<<std::endl;
    
}

// Move the Mesh as per the displacement calculated by Avergaes of MESH Velocities ( old and new )


    void deformMesh3D::move_mesh_ale(TFEVectFunct3D* MeshVelocityVectFunction3D,double* meshVelocity_old , double time_step)
{
	// Create a Vect Fucntion 2D to get the Mesh Co ordinates 
	TFESpace3D* Gridfespace = MeshVelocityVectFunction3D->GetFESpace3D();
	int N_DOF = Gridfespace->GetN_DegreesOfFreedom();
	cout << " N_ DOF in move mesh class  : " << N_DOF<<endl;
	double* gridpos = new double[3*N_DOF]();
	TFEVectFunct3D* gridFEVectfunc3D = new TFEVectFunct3D(Gridfespace, (char*)"C", (char*)"C", gridpos, N_DOF, 3);

	// value array of current Mesh velocity
	double* meshVelocityvalues = MeshVelocityVectFunction3D->GetValues();

	// get the Co-ordinates as data array using grid to data
	gridFEVectfunc3D->GridToData();

	// Now the increment in New position of each Vertes is given by ( (W_new - W_old * 0.5) * time_step )
	for ( int i = 0 ; i < 3*N_DOF ; i ++ ){
		gridpos[i] += (meshVelocityvalues[i] - meshVelocity_old[i])*0.5 * time_step;
	}
	// for ( int k = 0 ; k < 3*N_DOF ; k++)
	// 	cout << "  Grid pos ["<<k<<"] : "  << gridpos[k]  <<endl;

	// Now Move the Value from the array to Grid ( vertex Co ordinates )
	gridFEVectfunc3D->DataToGrid();
	delete[] gridpos;

	cout << " CLASS : DEFORM MESH 3D --  ' Mesh Has been Moved ' " <<endl;

}

