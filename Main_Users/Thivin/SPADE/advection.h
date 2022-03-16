// ======================================================================
// instationary problem
// ======================================================================

#include <MacroCell.h>
#include <IsoBoundEdge.h>
#include <IsoInterfaceJoint.h>

extern "C"
{
 #include <gridgen.h>    
  
  void triangulate(char*, struct triangulateio*,
                   struct triangulateio*, struct triangulateio*);
}



/// ========================================================================
// example file
// ========================================================================

#define __SIN3__

void ExampleFile()
{
  OutPut("Example: advection.h" << endl); 
}

// exact solution
void Exact(double x, double y, double *values)
{
  double t;

  t = TDatabase::TimeDB->CURRENTTIME;

  values[0] = exp(t)*(sin(2*Pi*x)*sin(2*Pi*y));
  values[1] = exp(t)*2*Pi*cos(2*Pi*x)*sin(2*Pi*y);
  values[2] = exp(t)*2*Pi*sin(2*Pi*x)*cos(2*Pi*y);
  values[3] = 0;
}

// kind of boundary condition (for FE space needed)
void BoundCondition(int BdComp, double t, BoundCond &cond)
{
//  cond = NEUMANN;

//  if(BdComp == 1 || BdComp == 3 || BdComp == 2 )
//  {
   cond = DIRICHLET;
//  }
}

// value of boundary condition
void BoundValue(int BdComp, double Param, double &value)
{
  value = 0;
}

// initial conditon
void InitialCondition(double x, double y, double *values)
{
  double t;

  t = TDatabase::TimeDB->CURRENTTIME;
  // values[0] = exp(t)*(sin(2*Pi*x)*sin(2*Pi*y));
}


void DO_Mean_Equation_Coefficients(int n_points, double *X, double *Y,
        double **parameters, double **coeffs)
{
  double eps=1/TDatabase::ParamDB->PE_NR;
  double b1=1, b2=-1, c=1;
  int i;
  double *coeff;
  double x, y;
  double t = TDatabase::TimeDB->CURRENTTIME;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    
    x = X[i];
    y = Y[i];

    coeff[0] = 0;
    coeff[1] = -sin(t)*0.2;
    coeff[2] = cos(t)*0.2;
    coeff[3] = 0;

    coeff[4] = 0.0;
  }
}

void DO_Mode_Equation_Coefficients(int n_points, double *X, double *Y,
        double **parameters, double **coeffs)
{
  double eps=1/TDatabase::ParamDB->PE_NR;
  double b1=1, b2=-1, c=1;
  int i;
  double *coeff;
  double x, y;
  double t = TDatabase::TimeDB->CURRENTTIME;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    
    x = X[i];
    y = Y[i];

    coeff[0] = 0;
    coeff[1] = -sin(t)*0.2;
    coeff[2] = cos(t)*0.2;
    coeff[3] = 0;

    coeff[4] = 0.0;
  }
}

// ======================================================================
// ASSEMBLY FUNCTION
// This fucntion will be called for Every Quadrature Point inside a Cell for Local Assembly
// ======================================================================
 void DO_Mean_Equation_Assembly(double quad_wt, double *coeff, double *param,
                     double hK, double **derivatives, int *N_BaseFuncts, double ***LocMatrices, double **LocRhs)
 {

	// The below values N_x, N_y , etc are ARRAY of Values 
    // Which provides Value of a Particular Shape function or its derivative ( Based on a DOF )
    // at the given Quadrature POint. 
    // Here the Size of the Array is equal to the NUmber of DOF in the cell 
    
    // For Eg : Orig0[1]  gives the derivative of Shape Function number 1 ( Shape function of DOF 1 ) w.r.t x
    // at the given Quadrature Point.
    // Take these values based on the MULTIINDEX DERIVATIVE
    double *N = derivatives[0];
    double *Nx = derivatives[1];
    double *Ny = derivatives[2];
    
	
	double **A11, *F1;
    double val = 0.;
    
	
	A11 = LocMatrices[0];  // Local Stiffenss matrix

    F1 = LocRhs[0];


    double c0 = coeff[0]; // nu
    double b1 = coeff[1]; // b1
    double b2 = coeff[2]; // b2
    double c  = coeff[3]; // c
    double f  = coeff[4]; //f

    int N_DOF_perCell = N_BaseFuncts[0];

    for (int i = 0; i < N_DOF_perCell; i++)  // Test
    {
        // Assemble RHS 
        F1[i] =  0.0; // TO DO
        double val = 0;
        for (int j = 0; j < N_DOF_perCell; j++) // Ansatz
        {

            val  +=  (b1*Nx[j]  + b2*Ny[j])* N[i];//   TO DO 
            // val  +=  c0*((Nx[j] * Nx[i])   + (Ny[j] * Ny[i]) );
            // val  +=  c*N[i]*N[j];

            A11[i][j] += val*quad_wt;
        } // endfor j
    }     // endfor i

 }


 // ======================================================================
// ASSEMBLY FUNCTION
// This fucntion will be called for Every Quadrature Point inside a Cell for Local Assembly
// ======================================================================
 void DO_Mode_Equation_Assembly(double quad_wt, double *coeff, double *param,
                     double hK, double **derivatives, int *N_BaseFuncts, double ***LocMatrices, double **LocRhs)
 {

	// The below values N_x, N_y , etc are ARRAY of Values 
    // Which provides Value of a Particular Shape function or its derivative ( Based on a DOF )
    // at the given Quadrature POint. 
    // Here the Size of the Array is equal to the NUmber of DOF in the cell 
    
    // For Eg : Orig0[1]  gives the derivative of Shape Function number 1 ( Shape function of DOF 1 ) w.r.t x
    // at the given Quadrature Point.
    // Take these values based on the MULTIINDEX DERIVATIVE
    double *N = derivatives[0];
    double *Nx = derivatives[1];
    double *Ny = derivatives[2];
    
	
	double **A11, *F1;
    double val = 0.;
    
	
	A11 = LocMatrices[0];  // Local Stiffenss matrix

    F1 = LocRhs[0];


    double c0 = coeff[0]; // nu
    double b1 = coeff[1]; // b1
    double b2 = coeff[2]; // b2
    double c  = coeff[3]; // c
    double f  = coeff[4]; //f

    int N_DOF_perCell = N_BaseFuncts[0];

    for (int i = 0; i < N_DOF_perCell; i++)  // Test
    {
        // Assemble RHS 
        F1[i] =  0.0; // TO DO
        double val = 0;
        for (int j = 0; j < N_DOF_perCell; j++) // Ansatz
        {

            val  +=  (b1*Nx[j]  + b2*Ny[j])* N[i];//   TO DO 
            val  +=  c0*((Nx[j] * Nx[i])   + (Ny[j] * Ny[i]) );
            val  +=  c*N[i]*N[j];

            A11[i][j] += val*quad_wt;
        } // endfor j
    }     // endfor i

 }


void DO_Mode_RHS_Aux_Param (double *in, double *out)
{
	out = in;
}
