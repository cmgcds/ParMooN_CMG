// ======================================================================
// Magnetic field in a sandwich hexagon
// ======================================================================
#include <ConvDiff3D.h>

void ExampleFile()
{
  OutPut("Example: HexaMF1Layer.h" << endl) ;
}

// exact solution
void Exact(double x, double y, double z, double *values)
{
  static double height = TDatabase::ParamDB->DRIFT_Z;
  static double Middle = 0.5*TDatabase::ParamDB->DRIFT_Z;
  static double H1 = TDatabase::ParamDB->P3;
  static double H2 = TDatabase::ParamDB->P4;
  static double LH = TDatabase::ParamDB->P1;
  static double b = (0.5 - LH/height)>0 ? (0.5 - LH/height) : 0;

  if(z>Middle)
  {
    // above the layer
    values[0] = (z-Middle);
    values[1] = 0;
    values[2] = 0;
    values[3] = 1;
    values[4] = 0;
  }
  else
  {
    if(z>b*height)
    {
      // in the layer
      values[0] = (z-Middle)*H1/H2;
      values[1] = 0;
      values[2] = 0;
      values[3] = H1/H2;
      values[4] = 0;
    }
    else
    {
      // below the layer
      values[0] = height*(b-0.5)*H1/H2 - b*height + z ;
      values[1] = 0;
      values[2] = 0;
      values[3] = 1;
      values[4] = 0;
    }
  }
}

// kind of boundary condition (for FE space needed)
void BoundCondition(double x, double y, double z, BoundCond &cond)
{
  static double height = TDatabase::ParamDB->DRIFT_Z;

  if(fabs(z)<1e-8 || fabs(z-height)<1e-8)
    cond = DIRICHLET;
  else
    cond = NEUMANN;
}

// value of boundary condition
void BoundValue(double x, double y, double z, double &value)
{
  static double height = TDatabase::ParamDB->DRIFT_Z;
  static double Middle = 0.5*height;
  static double H1 = TDatabase::ParamDB->P3;
  static double H2 = TDatabase::ParamDB->P4;
  static double LH = TDatabase::ParamDB->P1;
  static double b = (0.5 - LH/height)>0 ? (0.5 - LH/height) : 0;

  value = 0;

  if(fabs(z)<1e-8)
    value = height*( (H1/H2)*(b-0.5)-b );

  if(fabs(z-height)<1e-8)
    value = Middle;
}

void BilinearCoeffs(int n_points, double *x, double *y, double *z,
        double **parameters, double **coeffs)
{
  static double chi0 = TDatabase::ParamDB->P2-1;
  static double H0 = TDatabase::ParamDB->P4;
  static double Ms = TDatabase::ParamDB->P5;
  static double gamma = 3*chi0*H0/Ms;
  int i;
  double *coeff, *param, H, arg;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    param = parameters[i];

    H = param[0];

    coeff[0] = 1;
    coeff[1] = 0;
    coeff[2] = 0;
    coeff[3] = 0;
    coeff[4] = 0;
    coeff[5] = 0;

    arg = gamma*H;
    if(arg > 1e-2)
      coeff[6] = Ms/(H*H0) * ( 1/tanh(arg) - 1/arg ) + 1;
    else
      coeff[6] = chi0 + 1;

    coeff[7] = 1; // mu_vacuum
  }
}

void SetMuR(int N_Points, double **Coeffs, double **Params, TBaseCell *cell)
{
  int i, ID;
  double *Coeff;

  ID = cell->GetSubGridID();

  for(i=0;i<N_Points;i++)
  {
    Coeff = Coeffs[i];
    if(ID == 1) // Fluid
      Coeff[0] = Coeff[6];
    else // air/vacuum
      Coeff[0] = Coeff[7];
  } // endfor i
}

void ParametersMF(double *in, double *out)
{
  out[0] = sqrt(in[3]*in[3] + in[4]*in[4] + in[5]*in[5]);
}

static int N_FESpacesMF = 1;
static int N_FEFunctionsMF = 1;
static int N_ParamFctMF = 1;
static int N_FEValuesMF = 3;
static int FEValuesFctIndexMF[3] = { 0, 0, 0 };
static MultiIndex3D FEValuesMultiIndexMF[3] = { D100, D010, D001 };
static int N_ParametersMF = 1;
static int BeginParametersMF[1] = { 0 };
static ParamFct *ParameterFctMF[1] = { ParametersMF };

