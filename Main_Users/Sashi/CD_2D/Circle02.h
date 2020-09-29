// ======================================================================
//
// ======================================================================
#include <ConvDiff2D.h>

void ExampleFile()
{
  OutPut("Example: Circle02.h" << endl) ;
}

// exact solution
void Exact(double x, double y, double *values)
{
  values[0]= 4*(1-x*x-y*y)*x*y;
  values[1]= 4*y-12*x*x*y-4*y*y*y;
  values[2]= 4*x-4*x*x*x-12*x*y*y;
  values[3]= -48*x*y;
}

// kind of boundary condition (for FE space needed)
void BoundCondition(int BdComp, double t, BoundCond &cond)
{
 cond = DIRICHLET;
}

// value of boundary condition
void BoundValue(int BdComp, double Param, double &value)
{
  value = 0;
}

void BilinearCoeffs(int n_points, double *X, double *Y,
        double **parameters, double **coeffs)
{
  static double eps=1./TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff,*param; 
  double x, y ,c;

 c=1.;
  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
   param = parameters[i];

    x = X[i];
    y = Y[i];

    coeff[0] = eps; // diffusion coeffient
    coeff[1] = 0; // b1
    coeff[2] = 0; // b2
    coeff[3] = 1; // c

   
   

   // right-hand side
    coeff[4]= 48*x*y*eps+c*(4*x*y-4*x*x*x*y-4*x*y*y*y);

  }
}
