// ======================================================================
// instationary problem
// ======================================================================
#include <TimeConvDiff2D.h>

void ExampleFile()
{
  
#define __ROBINBC__ 
#define __TERAHERTZ__ 

  TDatabase::ParamDB->INTERNAL_CONVECTION_EQ_VELOFIELD=0;
    
  OutPut("Example: Time1.h" << endl);
}

// exact solution
void Exact(double x, double y, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

// exact solution
void Initial(double x, double y, double *values)
{
 values[0] = 300;
}

// kind of boundary condition (for FE space needed)
void BoundCondition(int BdComp, double t, BoundCond &cond)
{
    switch(BdComp)
     {
      case 0 : 
         cond = DIRICHLET;  
       break;     
       
      case  1:
         cond = ROBIN;
      break;      
      case 2 : 
         cond = NEUMANN;  
       break;
            
      case  3:
         cond = NEUMANN;
      break;  
      
       default:
            OutPut("Unknown BoundCondition" << endl);
            exit(4711);;
     }
}

// value of boundary condition
void BoundValue(int BdComp, double Param, double &value)
{
  value = 0;
}

void BilinearCoeffs_T(int n_points, double *X, double *Y,
        double **parameters, double **coeffs)
{
  double eps=1./TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff, *param;
  double x, y;
  double k0, alpha,  char_L;
  
  alpha = TDatabase::ParamDB->P1;
  char_L = TDatabase::ParamDB->P4;
  k0 = TDatabase::ParamDB->P5;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    param = parameters[i];

    x = X[i];
    y = Y[i];

    coeff[0] = eps;
    coeff[1] = 0;
    coeff[2] = 0;
    coeff[3] = 0;

    if(x<=1)
     {
      coeff[4] =  (alpha*exp(alpha*y*char_L))/(1e3*Pi*k0);
     }
    else
     {coeff[4] = 0; }
  }
}




void GetExampleFileData(BoundCondFunct2D **BoundaryConditions, BoundValueFunct2D **BoundValues, 
                        DoubleFunct2D **InitiaValues, CoeffFct2D **BilinearCoeffs, 
                        int &N_PBEqns, int &N_IndepntScalarEqns, int *Disctypes)
{

  TDatabase::ParamDB->INTERNAL_CONVECTION_EQ_VELOFIELD=0;


   N_IndepntScalarEqns = 1;
   N_PBEqns = 0;

   BilinearCoeffs[0] =  BilinearCoeffs_T;
   BoundaryConditions[0] = BoundCondition;
   BoundValues[0] = BoundValue;
   InitiaValues[0] = Initial;
   Disctypes[0] = GALERKIN;
}





