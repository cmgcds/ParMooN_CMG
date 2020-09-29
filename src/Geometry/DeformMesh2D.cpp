/** ************************************************************************ 
* @brief     source file for deformMesh2D 
             Moves the computational domain based on the boundary conditions provided by solving 2D linear Elasticity equation

             Parent Class : <same> deformMesh2D

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
#include <DiscreteForm2D.h>
#include <FEDatabase2D.h>
#include <LinAlg.h>
#include <FESpace2D.h>
#include <SquareStructure2D.h>
#include <Structure2D.h>
#include <Output2D.h>
#include <LinAlg.h>
#include <Assemble2D.h>
#include <DirectSolver.h>
#include <MainUtilities.h>
#include <NSE2D_ParamRout.h>
#include <LinAlg.h> 
#include <string.h>
#include <sstream>
#include <stdlib.h>
#include <numeric>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "DeformMesh2D.h"

// Test 
#include <vector>
#include<fstream>

using namespace std;

deformMesh2D::deformMesh2D(TCollection* coll1, BoundCondFunct2D* BoundCondition_X1, 
BoundCondFunct2D* BoundCondition_Y1, BoundValueFunct2D* BoundValue_X1,  BoundValueFunct2D* BoundValue_Y1)
{
    BoundCondition_X  = BoundCondition_X1;
    BoundCondition_Y  = BoundCondition_Y1;
    BoundValue_X      = BoundValue_X1;
    BoundValue_Y      = BoundValue_Y1;
    coll              = coll1;

    // Call the initialising function 
    deformMesh2D::moveMesh();
}



void deformMesh2D::Assembly_poisson_2D(double quad_wt, double *coeff, double *param,
                                      double hK, double **derivatives, int *N_BaseFuncts,double ***LocMatrices, double **LocRhs)
{

  double *N = derivatives[0], *Nx = derivatives[1], *Ny = derivatives[2];
  double **K11, **K12, **K21, **K22, *F1, *F2;
  K11 = LocMatrices[0];
  K12 = LocMatrices[1];
  K21 = LocMatrices[2];
  K22 = LocMatrices[3];

  F1 = LocRhs[0];
  F2 = LocRhs[1];
  double mu = 1;

  for (int i = 0; i < N_BaseFuncts[0]; i++){
    for (int j = 0; j < N_BaseFuncts[0]; j++){
     

       K11[i][j] += quad_wt * (Nx[i]*Nx[j] + Ny[i]*Ny[j]); 

      //cout << "Nx[ " << i << "] = " << Nx[i] << endl;

      K12[i][j] += 0.0;
      K21[j][i] += 0.0;

      K22[i][j] += quad_wt * (Nx[i]*Nx[j] + Ny[i]*Ny[j]);

      /* NON LINEAR PART */
    }
    /* RHS */

    F1[i] = 0.;

    F2[i] = 0.;
  }

}


void deformMesh2D::BilinearCoeffs(int n_points, double *x, double *y,double **parameters, double **coeffs)
{
  double angle = 0, v1, v2;
  int i;
  double *coeff;

  //v1 = cos(angle);
  //v2 = sin(angle);
}

void deformMesh2D::hello()
{
  cout <<  " " << endl;
}

void deformMesh2D::SolveGridEquation(double **Entries, double *sol, double *rhs,
                       int *KCol, int *RowPtr, int N_DOF)
{
  int i,j,k, col, Diognal;
  double *Entries11, *Entries12, *Entries21, *Entries22;
  double sum1, sum2, max_sum1, max_sum2;
  int start, end;

  double max_error, error=1.e-4;
  int iter;

  Entries11 = Entries[0];
  Entries12 = Entries[1];
  Entries21 = Entries[2];
  Entries22 = Entries[3];
   
  max_error = 1.; iter = 0;
  while(max_error>error)
  {
    max_error = 0.0; iter++;
    for(i=0;i<N_DOF;i++)
    {
      start = RowPtr[i];
      end = RowPtr[i+1];
      sum1 = rhs[i];
      sum2 = rhs[i+N_DOF];
      for(k=start;k<end;k++)
      {
        col = KCol[k];
        if (col==i)  Diognal = k;
	sum1 -= Entries11[k] * sol[col] 
	      +Entries12[k] * sol[col+N_DOF];
	sum2 -= Entries21[k] * sol[col]
              +Entries22[k] * sol[col+N_DOF];
      } // endfor k
        sol[i] += sum1/Entries11[Diognal];
        sol[i+N_DOF] += sum2/Entries22[Diognal];
      if(max_error<fabs(sum1/Entries11[Diognal])) max_error = fabs(sum1/Entries11[Diognal]);
      if(max_error<fabs(sum2/Entries22[Diognal])) max_error = fabs(sum2/Entries22[Diognal]);
    } // endfor i
    if(iter == 10) break;
  } // end while
// OutPut("Grid Solver: Number iteration "<<iter<<endl);
// exit(0);

}




void deformMesh2D::moveMesh()
{   
    
    //FEVect funciton 2D to get the co-ordinates of the grid from FESPACE
    TFEVectFunct2D *displacements;
   
   	  /* ----------- Finite element Space Used ----------- /
	   *  Order of the Element ,Shape functions , Derivative matices  */
    TFESpace2D *fesp[1],*ferhs[2];
    
    // Aux parameters - generates values like u,x u,y,  which are used for the construction of local matrices
    // NOTE : None of the values are generated from AUX in this routine , since the weak form of elasticity equation 
    // does not need any special terms . So we will pass the AUX pointer as "NULL"  to the discrete form.
    TAuxParam2D *aux;

    TDiscreteForm2D *discreteform;

    // The derivatives required for solving equations
    MultiIndex2D AllDerivatives[3] = {D00, D10, D01};
    TFESpace2D *fespace;
    
   
    /* For VTK */
    const char vtkdir[] = "VTK";
	  char *PsBaseName, *VtkBaseName, *GEO;
	  char UString[] = "u";
	  char PString[] = "p";
	  std::ostringstream os;
	  os << " ";
    /* End for VTK */



    // FE SPACE 
    // Define FEspace using the obtained collections
    // We are interested in only findig the displacements of the mesh at the corner nodes, so it is enough to pass the Finite element order as "1"
    fespace = new TFESpace2D(coll, (char *)"name", (char *)"description", BoundCondition_X, 1, NULL);


    // Get the number of Unknowns and DOF and allocate memory for the respective arrays
    int N_U = fespace->GetN_DegreesOfFreedom();
    int N_Active = fespace->GetActiveBound();
    int N_DOF = 2*N_U; // for x and Y 
    
    double* gridpos = new double[N_DOF]();
    double* rhs = new double[N_DOF](); 
    double* solution = new double[N_DOF]();   

    // define a Vect FUNC 2D for the storing the grid co-ordinates
    TFEVectFunct2D* gridCoordinates = new TFEVectFunct2D(fespace, (char*)"C", (char*)"C", gridpos, N_U, 2);
    
    // get the co-ordinates(X,Y) of the grid to the datastructure of FEVECT2d into the array "gridpos"
    gridCoordinates -> GridToData();
     
    // Declare a FEVECTFunc2d for the solution array which will hold all values of solution of the linear elastic system in 2D
    displacements = new TFEVectFunct2D(fespace, (char *)"disp", (char *)"description", solution, N_U, 2);

     // Update the number of FESPACES that needs to be used as an array in the FESPACE2D object
     fesp[0] = fespace;

    // Aux - default Null type for AUX
    aux = new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);

    
  /* ------------------------------- start of DISCRETE FORM ---------------------------- */
  // N_Terms  -  Number of terms and derivatives needed , in this case N , Nx, Ny
  int N_Terms = 3;
  // Spacenumbers -- array of values which denotes the SPACE numbers that eac terms mentioned above belongs to
  int *SpacesNumbers = new int[N_Terms]();
	// Number of block matrices needed to create our final stiffness matrix
 	// For NSE 2d it will be 4 A matrices and 2 B matrices
  int N_Matrices = 4;
  // The number of components in the RHS of the matrix
	int N_RHS = 2;
	// The row space and column space of the each matrix being constructed
	int *rowspace = new int[N_Matrices]();
	int *columnspace = new int[N_Matrices]();
	int *rhsspace = new int[N_RHS]();


    discreteform = new TDiscreteForm2D(UString,UString,N_Terms,AllDerivatives,SpacesNumbers,N_Matrices,N_RHS,
                                        rowspace,columnspace,rhsspace,Assembly_poisson_2D,BilinearCoeffs,NULL);

    /** constructor with assembling using parameters */


  /* ------------------------------- end of DISCRETE FORM ---------------------------- */

   
   /*------------------ Start of MATRICES ----------------- */ 
    // Matrix - structure declaration for square matrices
    // This is a module which converts gives an interface between normal matrices and the matrices
    // stored in the CSR format //
    TSquareStructure2D *sqstructure = new TSquareStructure2D(fespace);
    sqstructure->Sort();

    // Once the strucute of the matrix is finalised for a particular FE space
    // assign a matrix for that particular structure
    double *RHS[2];
    TSquareMatrix2D *sqmatrixA11, *sqmatrixA12, *sqmatrixA21, *sqmatrixA22, *SQMATRICES[4];
    sqmatrixA11 = new TSquareMatrix2D(sqstructure);
    sqmatrixA12 = new TSquareMatrix2D(sqstructure);
    sqmatrixA21 = new TSquareMatrix2D(sqstructure);
    sqmatrixA22 = new TSquareMatrix2D(sqstructure);

    SQMATRICES[0] = sqmatrixA11;
    SQMATRICES[1] = sqmatrixA12;
    SQMATRICES[2] = sqmatrixA21;
    SQMATRICES[3] = sqmatrixA22;


    RHS[0] = rhs;
    RHS[1] = rhs + N_U;

  /*------------------ End  of MATRICES ----------------- */ 


  /*------------------ Start of BOUNDARY CONDITIONS and values----------------- */ 
    BoundCondFunct2D *BoundaryConditions_Neumann[2];
    BoundValueFunct2D *BoundaryValues_Neumann[2], *BoundaryValues_Iter2[2];

    BoundaryConditions_Neumann[0] = BoundCondition_X;
    BoundaryConditions_Neumann[1] = BoundCondition_Y;

    // Assigning the boundary  value function pointers to the above declared pointers
    BoundaryValues_Neumann[0] = BoundValue_X;
    BoundaryValues_Neumann[1] = BoundValue_Y;

/*------------------ Start of BOUNDARY CONDITIONS and values----------------- */ 

  // Array of FEspaces that needs to be passed to the Assembly2D system 

    fesp[0] = fespace;
    ferhs[0] = fespace;
    ferhs[1] = fespace;
    

    
    // Assemble 2D - Functions
          Assemble2D(1, fesp,
                 4, SQMATRICES,
                 0, NULL,
                 2, RHS, ferhs,
                 discreteform,
                 BoundaryConditions_Neumann,
				 BoundaryValues_Neumann,
                 aux);
    
   
    double* solution_iterative = new double[N_DOF]();
    double** ENTRIES = new double*[4];

    ENTRIES[0] = SQMATRICES[0]->GetEntries();
    ENTRIES[1] = SQMATRICES[1]->GetEntries();
    ENTRIES[2] = SQMATRICES[2]->GetEntries();
    ENTRIES[3] = SQMATRICES[3]->GetEntries();

    
    // Solve the System using DIRECT SOLVER part
    DirectSolver(sqmatrixA11, sqmatrixA12, sqmatrixA21, sqmatrixA22, rhs, rhs + N_U, solution_iterative, solution_iterative + N_U);


    cout << " --- Before Iterative Solver ------" << endl;
    // SolveGridEquation(ENTRIES, solution_iterative, rhs,
    //                    SQMATRICES[0]->GetKCol(), SQMATRICES[0]->GetRowPtr(), N_U);
    cout << " --- After Iterative Solver ------" << endl;
            
    // ------------------- PRINT THE MATRIX TO FILE --------------------------------- //
    
   
    // std::vector<vector < double>> a(N_DOF);

    // for (int i =0 ; i<N_DOF ; i++)
    //    a[i] = vector<double>(N_DOF);

    // for ( int i = 0 ; i < N_U ; i++)
    // {
    //   for (int j=0 ; j<N_U ; j++)
    //   {
    //     a[i][j] = sqmatrixA11->get(i,j);
    //     a[i][j+N_U] = sqmatrixA12->get(i,j);
    //     a[i+N_U][j] = sqmatrixA21->get(i,j);
    //     a[i+N_U][j+N_U] = sqmatrixA22->get(i,j);
    //   }
    // }

    // string filename_matrix = "step_mat" + to_string(TDatabase::ParamDB->UNIFORM_STEPS) + ".mat";
    // string filename_rhs = "step_rhs" + to_string(TDatabase::ParamDB->UNIFORM_STEPS) + ".mat";
    // string filename_sol = "step_sol" + to_string(TDatabase::ParamDB->UNIFORM_STEPS) + ".mat";
    // ofstream file(filename_matrix);

    // for(int i=0;i<N_DOF; i++)
    // {
    //   for(int j=0;j<N_DOF; j++)
    //   {
    //     file<<a[i][j]<<"\t";
    //   }
    //     file<<endl;
    // }

    // file.close();

    // ofstream file2(filename_rhs);
    // for(int i=0;i<N_DOF; i++)
    //   file2<<rhs[i]<<"\t";
    
    // file2.close();

    // ofstream file3(filename_sol);
    // for(int i=0;i<N_DOF; i++)
    //   file3<<solution_iterative[i]<<"\t";
    
    // file3.close();


    // cout  << " Grid Co-ordinates inside class - Before addition" << endl;
    // for( int i = 0 ; i < N_DOF ; i++)
    //   cout << gridpos[i] << "\t";
    // cout << endl;

    // -------------------End of PRINT THE MATRIX TO FILE --------------------------------- //

    // Add the Currnt Solution of the system to the existing grid points to move the grids to new position
    for ( int i=0; i< N_DOF ; i++)
        gridpos[i] += solution_iterative[i];

    // Reset matrices to free up the memory requirement
    sqmatrixA11->Reset();
    sqmatrixA12->Reset();
    sqmatrixA21->Reset();
    sqmatrixA22->Reset();
    for (int i_rhs = 0; i_rhs < N_DOF; i_rhs++)
      rhs[i_rhs] = 0;

  // Push the new computed values of the X,Y array to the actual coputational grid data structure
  gridCoordinates -> DataToGrid();


  cout << " Class : DeformMeshElasticity -- >  Mesh is successfully moved " << endl;

}