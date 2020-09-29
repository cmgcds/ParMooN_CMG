// ======================================================================
//
// ======================================================================
#include <ConvDiff2D.h>

void ExampleFile()
{
  OutPut("Example: Circle03.h" << endl) ;
}

// exact solution
void Exact(double x, double y,  double *values)
{
  double r, theeta;

  r = sqrt(x*x + y*y);
  theeta = atan2(y,x);

 if(theeta<0)
  theeta= 2.*Pi + theeta;

  values[0]= pow(r,2./3.)* sin(2./3.*theeta);
  values[1]= 2./3.*pow(r,-1./3.)*sin(2./3.*theeta-theeta);
  values[2]= 2./3.*pow(r,-1./3.)*cos(2./3.*theeta-theeta);
  values[3]= 0;
}

// kind of boundary condition (for FE space needed)
void BoundCondition(int BdComp, double t, BoundCond &cond)
{
 cond = DIRICHLET;
}

// value of boundary condition
void BoundValue(int BdComp, double Param, double &value)
{
double x, y, theeta;

if(BdComp==2)
{
  theeta = Param*3.*180./2.;
 // cout << " theeta " << theeta << endl;
  value=sin((2./3.)*theeta*Pi/180.);
}
 else
  value = 0;

}

void BilinearCoeffs(int n_points, double *X, double *Y,
        double **parameters, double **coeffs)
{
  static double eps=1./TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff,*param; 
 
  for(i=0;i<n_points;i++)
  {
   coeff = coeffs[i];
   param = parameters[i];

    coeff[0] = eps; // diffusion coeffient
    coeff[1] = 0; // b1
    coeff[2] = 0; // b2
    coeff[3] = 0; // c
  
      // right-hand side
    coeff[4]= 0;
  }
}
