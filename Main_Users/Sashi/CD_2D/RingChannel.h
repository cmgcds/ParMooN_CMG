// ======================================================================
// magnetic fluid in a ring channel
// ======================================================================
#include <ConvDiff2D.h>

void ExampleFile()
{
  OutPut("Example: TwoInteriorLayers.h" << endl) ;
}

// exact solution
void Exact(double x, double y, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

// kind of boundary condition (for FE space needed)
void BoundCondition(int i, double t, BoundCond &cond)
{
  cond = DIRICHLET;
}

// value of boundary condition
void BoundValue(int BdComp, double Param, double &value)
{
  switch(BdComp)
  {
    case 0: value=0; // DIRICHLET
            break;
    case 1: value=Param; // DIRICHLET
            break;
    case 2: value=1; // DIRICHLET
            break;
    case 3: value=1-Param; // DIRICHLET
            break;
    default: cout << "wrong boundary part number" << endl;
            break;
  }
}

void BilinearCoeffs(int n_points, double *X, double *Y,
        double **parameters, double **coeffs)
{
  static double mu=TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff, x, y;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    x = X[i];
    y = Y[i];

    if((fabs(x-0.5) < 0.01666666666666666667) && 
       (fabs(y-0.5) < .020833333333333333333) )
    {
      coeff[0] = mu*x;
    }
    else
    {
      coeff[0] = x;
    }

    coeff[1] = 0;
    coeff[2] = 0;
    coeff[3] = 0;

    coeff[4] = 0;
  }
}

