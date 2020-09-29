// ======================================================================
// two boundary layer problem, from [JMT97]
// ======================================================================
// #include <ConvDiff2D.h>

void ExampleFile()
{
  OutPut("Example: TwoBoundaryLayersLaplace.h" << endl) ;
}

// exact solution
void Exact(double x, double y, double *values)
{
  static double Reynolds = TDatabase::ParamDB->RE_NR;
  double t29,t33,t37,t41,t43,t44,t45,t46,t48;

  values[0]= x*y*y -y*y*exp(-2.0*(1.0-x)*Reynolds)
       -x*exp(-3.0*(1.0-y)*Reynolds)
       +exp(-(5.0-2.0*x-3.0*y)*Reynolds);

  t29 = y*y;
  t33 = exp(-2.0*(1.0-x)*Reynolds);
  t37 = exp(-3.0*(1.0-y)*Reynolds);
  t41 = exp(-(5.0-2.0*x-3.0*y)*Reynolds);
  t43 = t29*Reynolds*t33;
  t44 = Reynolds*t41;
  values[1]  = t29-2.0*t43-t37+2.0*t44;
  t45 = x*y;
  t46 = y*t33;
  t48 = x*Reynolds*t37;
  values[2]  = 2.0*t45-2.0*t46-3.0*t48+3.0*t44;  

  values[3] = 0; // laplace
}

// kind of boundary condition (for FE space needed)
void BoundCondition(int i, double t, BoundCond &cond)
{
  cond = DIRICHLET;
}

// value of boundary condition
void BoundValue(int BdComp, double Param, double &value)
{
  static double Reynolds = TDatabase::ParamDB->RE_NR;
  double x;

  switch(BdComp)
  {
    case 0: value = -Param*exp(-3*Reynolds) + exp((-5+2*Param)*Reynolds);
            break;
    case 1: value=0;
            break;
    case 2: value=0;
            break;
    case 3: x=1.0-Param;
            value = -exp(-2.0*Reynolds)+2.0*exp(-2.0*Reynolds)*Param
                    -exp(-2.0*Reynolds)*Param*Param
                    +exp(-1.0*(2.0+3.0*Param)*Reynolds);
            break;
    default: cout << "wrong boundary part number" << endl;
            break;
  }
}

void BilinearCoeffs(int n_points, double *X, double *Y,
        double **parameters, double **coeffs)
{
  static double eps=1/TDatabase::ParamDB->RE_NR;
  static double Reynolds=TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff, *param;
  double x, y;
  double u,t29,t33,t37,t41,t43,t44,t45,t46,t48,t49,t51,t52,t54;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    param = parameters[i];

    coeff[0] = eps; 
    coeff[1] = 0; // b1
    coeff[2] = 0; // b2
    coeff[3] = 0; // c

    x = X[i];
    y = Y[i];

    if (x<0) x=0;
    if (y<0) y=0;
    if (x>1) x=1;
    if (y>1) y=1;

    t29 = y*y;
    t33 = exp(-2.0*(1.0-x)*Reynolds);
    t37 = exp(-3.0*(1.0-y)*Reynolds);
    t41 = exp(-(5.0-2.0*x-3.0*y)*Reynolds);
//    cout << x << t33 << " " << t37 << " " << t41 << endl;
    u = x*t29-t29*t33-x*t37+t41;
    t43 = t29*Reynolds*t33;
    t44 = Reynolds*t41;
    t45 = x*y;
    t46 = y*t33;
    t48 = x*Reynolds*t37;
    t49 = Reynolds*Reynolds;
    t51 = t29*t49*t33;
    t52 = t49*t41;
    t54 = x*t49*t37;

    // right-hand side
    coeff[4]= -(-4.0*t51+13.0*t52+2.0*x-2.0*t33-9.0*t54)/Reynolds;
  }
}

