/** ************************************************************************ 
* @brief     source file for DeformMeshElasticity 
             Moves the computational domain based on the boundary conditions provided by solving 2D linear Elasticity equation

             Parent Class : <same> deformMeshElasticity

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

#ifndef _DEFORMMESH3D_
#define _DEFORMMESH3D_

class deformMesh3D{

private:
    TCollection* coll;
    BoundCondFunct3D **GridBoundaryConditions;
	  BoundValueFunct3D **GridBoundValues;
    

    //TFESpace2D *fespace;
    int N_DOF;

  // Local assembly function , which has the definition for linear elasticity assembly for local matrices in 2D
    static void Assembly_poisson_3D(double quad_wt, double *coeff, double *param, double hK, 
                  double **derivatives, int *N_BaseFuncts, double ***LocMatrices, double **LocRhs);

    // Function which is required for local matrix assembly in linear elasticity equation
    // NOTE : None of the values mentioned in this function is used in current calculation 
    //        it is just added for the dependency purposes
    static void GridCoeffs(int n_points, double *x, double *y,double *z, double **parameters, double **coeffs);

    //Function to Move the mesh 

    static void SolveGridEquation(double **Entries, double *sol, double *rhs,
                        int *KCol, int *RowPtr, int N_DOF);

public:

  TFESpace3D *fesp[1],*ferhs[3];


  // Initialising the Variables for Class - ALE MESH MOVEMENT
  TFESpace3D* fespace_alemeshMovement;
  TBaseCell** MovCells; 				
  TVertex **MovBoundVert;
  TJoint **MovJoints;
  int N_freeSurfaceJoints;
	int N_freeSurfaceCells;
	int N_freeSurfaceVertex;

  std::vector<int> freeSurfaceVertex;
  std::vector<int> freeSurfaceCells;
  std::vector<int> freeSurfaceJoints;
  std::vector<int> freeSurfaceVertexLocal;

  std::vector<double> meshVelocityatFreeSurface_1;
  std::vector<double> meshVelocityatFreeSurface_2;
  std::vector<double> meshVelocityatFreeSurface_3;

  std::vector<double>  xi_1; // Arrays to store temp XI valu on 2d Face
	std::vector<double>  xi_2; // Arrays to store temp XI valu on 2d Face 


  std::vector<double> xi_freeSurfaceNodes;
  std::vector<double> eta_freeSurfaceNodes;
  std::vector<double> zeta_freeSurfaceNodes;





/* CONSTRUCTOR /
 */
deformMesh3D(TCollection* _coll,BoundCondFunct3D** _GridBoundaryConditions, BoundValueFunct3D** _GridBoundValues);


// Main function , which solves the linear elasticity equation in 2d to move the computational domain.
void moveMesh();

// Pick up the Nodes based on the Boundary ID provided by the user
void Pick_free_surface_DOFs(TFESpace3D* fespace, std::vector<int> boundIds ,TCollection* coll );

// function to Calculate the Free surface Displacements of the Free Surface
  void get_surface_normals(TFEVectFunct3D* Velocity_FEvect,TFEVectFunct3D* MeshVelo_FEvect );

  // Function to Interpolate Mesh Velocity to the Interior Points
  void Solve_mesh_velocity_interior(TFESpace3D* fespace, TFEVectFunct3D* MeshVelocityVectFunction3D);

  // Function to Move the Mesh for ALE
  void move_mesh_ale(TFEVectFunct3D* MeshVelocityVectFunction3D,double* meshVelocity_old ,  double time_step);



/* DESTRUCTOR */
~deformMesh3D();


};

#endif