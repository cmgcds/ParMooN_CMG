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


#ifndef _DEFORMMESH2D_
#define _DEFORMMESH2D_

class deformMesh2D{

private:
  TCollection* coll;
  BoundCondFunct2D* BoundCondition_X;
  BoundCondFunct2D* BoundCondition_Y;
  BoundValueFunct2D* BoundValue_X;
  BoundValueFunct2D* BoundValue_Y;
  //TFESpace2D *fespace;
  int N_DOF;

  // Local assembly function , which has the definition for linear elasticity assembly for local matrices in 2D
  static void  Assembly_poisson_2D(double quad_wt, double *coeff, double *param,
					double hK, double **derivatives, int *N_BaseFuncts,
					double ***LocMatrices, double **LocRhs);

  // Function which is required for local matrix assembly in linear elasticity equation
  // NOTE : None of the values mentioned in this function is used in current calculation 
  //        it is just added for the dependency purposes
    static void  BilinearCoeffs(int n_points, double *x, double *y,
        double **parameters, double **coeffs);

   static void SolveGridEquation(double **Entries, double *sol, double *rhs,
                       int *KCol, int *RowPtr, int N_DOF);

public:
/* CONSTRUCTOR /
Param 1  --  the fespace pointer 
param 2 -- Boundary condition function in 2D - for X direction
param 3 -- Boundary condition function in 2D - for Y direction
param 4 -- Boundary Value function in 2D for X direction
param 5 -- Boundary value function in 2D for Y direction  */
deformMesh2D(TCollection* coll, BoundCondFunct2D BoundCondition_X, 
BoundCondFunct2D BoundCondition_Y, BoundValueFunct2D BoundValue_X,  BoundValueFunct2D BoundValue_Y);

/* DESTRUCTOR */
~deformMesh2D();

// Main function , which solves the linear elasticity equation in 2d to move the computational domain.
void moveMesh();

void hello();



};

#endif