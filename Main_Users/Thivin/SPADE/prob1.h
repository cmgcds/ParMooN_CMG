// Navier-Stokes problem, Driven cavity
//
// u(x,y) = unknown
// p(x,y) = unknown

void ExampleFile()
{
    OutPut("Example: DrivenCavity.h" << endl);
}

// ========================================================================
// exact solution   -- NO need to take into account for this Assignment
// ========================================================================
void Exact(double x, double y, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

///////////// BOUNDARY COMPONENT ID  /////////////////////////////////
/* 
    FOr a Square Domain, The edge y= 0    ID  = 0
                         The Edge x = 1   ID  = 1
                         The edge y = 1   ID  = 2
                         The edge x = 0   ID  = 3   
                         Inner Circle     ID  = 4
*/
//////


/////////// HEMKER PROBLEM /////////////////////////
// Simple perturbed Eliptic Equation 
//                       Eps( grad)^u + Ux  = 0 
// Similar to Convection diffusion advection equation ( CD2D) with b1 = 1 , b2 = 0, c = 0 
// Boundary Conditions, at the unit circle, the Value of unknown wil be "1" ( Diriclet )
// Its zero on inlet and walls, No flux on outlet



// ========================================================================
// boundary conditions
// Provide The Boundary Condition as "DIRICHLET" or "NEUMANN"
// ========================================================================
// TO DO

void BoundCondition(int BdComp, double t, BoundCond &cond)
{
    switch(BdComp)
    {
	case 0:
	case 1:
	case 2:
	    cond = NEUMANN;
	    break;
	default:
	    cond = DIRICHLET;
    }
}


// value of boundary condition
void BoundValue(int BdComp, double Param, double &value)
{
  switch(BdComp)
  {
    case 1:
      value = 0;
      break;
    case 4:
      value = 1;
      break;
    default:
      value = 0;
  }
}

// initial conditon
void InitialCondition(double x,  double y, double *values)
{
  values[0] = 0;
}

void BoundConditionAdjoint(int BdComp, double t, BoundCond &cond)
{
    switch(BdComp)
    {
	case 0:
	case 2:
	case 3:
	    cond = NEUMANN;
	    break;
	default:
	    cond = DIRICHLET;
    }
}


// value of boundary condition
void BoundValueAdjoint(int BdComp, double Param, double &value)
{
    value = 0;
}

void BilinearCoeffs(int n_points, double *x, double *y,
        double **parameters, double **coeffs)
{
  double eps=1/TDatabase::ParamDB->RE_NR;
  double angle = 0, v1, v2;
  int i;
  double *coeff;

  v1 = cos(angle);
  v2 = sin(angle);

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    coeff[0] = eps;
    coeff[1] = v1;
    coeff[2] = v2;
    coeff[3] = 0;

    coeff[4] = 0;
  }
}

// ======================================================================
// ASSEMBLY FUNCTION
// This fucntion will be called for Every Quadrature Point inside a Cell for Local Assembly
// ======================================================================
 void HemkerAssembly(double quad_wt, double *coeff, double *param,
                     double hK, double **derivatives, int *N_BaseFuncts, double ***LocMatrices, double **LocRhs)
 {

	// The below values N_x, N_y , etc are ARRAY of Values 
    // Which provides Value of a Particular Shape function or its derivative ( Based on a DOF )
    // at the given Quadrature POint. 
    //
    // Here the Size of the Array is equal to the NUmber of DOF in the cell 
    
    // For Eg : Orig0[1]  gives the derivative of Shape Function number 1 ( Shape function of DOF 1 ) w.r.t x
    // at the given Quadrature Point.

    double *N = derivatives[0], *Nx = derivatives[1], *Ny = derivatives[2];
    
	
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
        for (int j = 0; j < N_DOF_perCell; j++) // Ansatz
        {

            A11[i][j]  +=  (b1*Nx[i][j]  + b2*Ny[i][j])* N[i][j] ;//   TO DO 
        } // endfor j
    }     // endfor i

 }


