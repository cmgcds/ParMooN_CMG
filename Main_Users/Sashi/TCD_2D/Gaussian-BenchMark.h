
void ExampleFile()
{
  OutPut("Example: Gaussian-BenchMark.h" << endl);
}


// exact solution
void Exact(double x, double y, double *values)
{
  double t = TDatabase::TimeDB->CURRENTTIME;
  double t1 = cos(t+Pi/4);
  double t2 = sin(t+Pi/4);
  
  double tt = (x-0.3*t1)*(x-0.3*t1) + (y-0.3*t2)*(y-0.3*t2);
  double val = exp(-10*tt);
 
  double valx = -20*val*(x-0.3*t1);
  double valy = -20*val*(y-0.3*t2);

  values[0] = val;
    
  values[1] = valx;
   
  values[2] = valy;
}


// kind of boundary condition (for FE space needed)

void BoundCondition(int BdComp, double t, BoundCond &cond)
{

    cond = DIRICHLET;
}


// value of boundary condition
void BoundValue(int BdComp, double Param, double &value)
{
   double t = TDatabase::TimeDB->CURRENTTIME;
   
   double x = cos(2*Pi*Param);
   double y = sin(2*Pi*Param);
   double t1 = cos(t+Pi/4);
   double t2 = sin(t+Pi/4);
   double val = (x-0.3*t1)*(x-0.3*t1)+(y-0.3*t2)*(y-0.3*t2);

   value = exp(-10*val);
   //value = 0;

}


// initial conditon
void InitialCondition(double x,  double y, double *values)
{
  double t1 = cos(Pi/4);
  double t2 = sin(Pi/4);
  double val = (x-0.3*t1)*(x-0.3*t1)+(y-0.3*t2)*(y-0.3*t2);

  values[0] = exp(-10.0*val);
}


void BilinearCoeffs(int n_points, double *X, double *Y,
double **parameters, double **coeffs)
{
  double eps = 1.0/TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff;                                  // *param;
  double x, y;
  
  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    // param = parameters[i];

    x = X[i];
    y = Y[i];
      
    // diffusion
    coeff[0] = 0.0;
    // convection in x direction
    coeff[1] = -y;
    // convection in y direction
    coeff[2] = x;
    // reaction
    coeff[3] = 0.0;
    // rhs
    coeff[4] = 0.0;
  }
}

