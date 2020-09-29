// ======================================================================
// magnetic potential in the drop
// ======================================================================
#include <ConvDiff2D.h>

void ExampleFile()
{
  OutPut("Example: MFDrop_axBEM+FEM.h" << endl) ;
}

// exact solution for linear case
// ------------------------------
void Exact(double x, double y, double *values)
{

}

// kind of boundary condition (for FE space needed)
void BoundCondition(int BdComp, double t, BoundCond &cond)
{
  if(BdComp==0)
    cond = DIRICHLET;
  if(BdComp==1)
    //cond = DIRICHLET; // if we use Phi from BEM
  cond = NEUMANN; //if we use Q from BEM
  if(BdComp==2)
    cond = NEUMANN;
}

// Dirichlet boundary condition in any case
void DirichletBoundCondition(int BdComp, double t, BoundCond &cond)
{
  cond = DIRICHLET;
}

// zero boundary condition
void ZeroBoundValue(int BdComp, double t, double &value)
{
  value = 0;
}

// value of boundary condition
void BoundValue(int BdComp, double t, double &value)
{
  switch(BdComp)
  {
    case 0: value = 0;
            break;
    case 1: value = 0; 
            break;
    case 2: value = 0;
            break;
    default: cout << "Wrong boundary part number" << endl;
            break;
  }
}

void SimpleBilinearCoeffs(int n_points, double *x, double *y,
        double **parameters, double **coeffs)
{
  int i;
  double *coeff, *param;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    param = parameters[i];

    coeff[0] = 1;
    coeff[1] = 0;
    coeff[2] = 0;
    coeff[3] = 0;
    coeff[4] = 0;
  }
}

void BilinearCoeffs(int n_points, double *x, double *y,
        double **parameters, double **coeffs)
{
  static double mu = TDatabase::ParamDB->P1;

  int i;
  double *coeff, *param;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    param = parameters[i];

    coeff[0] = mu*x[i];
    coeff[1] = 0;
    coeff[2] = 0;
    coeff[3] = 0;
    coeff[4] = 0;
  }
}

void ParametersMF(double *in, double *out)
{
  out[0] = sqrt(in[2]*in[2] + in[3]*in[3]);
}

// This information is needed for the instances of the class TAuxParam
static int N_FESpacesMF = 1;
static int N_FEFunctionsMF = 1;
static int N_ParamFctMF = 1;
static int N_FEValuesMF = 2;
static ParamFct *ParameterFctMF[1] = { ParametersMF };
static int FEValuesFctIndexMF[2] = { 0, 0};
static MultiIndex2D FEValuesMultiIndexMF[2] = { D10, D01 };
static int N_ParametersMF = 1;
static int BeginParametersMF[1] = { 0 };

int GridN_Terms = 2;
MultiIndex2D GridDerivatives[2] = { D10, D01 };
int GridSpacesNumbers[2] = { 0, 0 };
int GridN_Matrices = 4;
int GridRowSpace[4] = { 0, 0, 0, 0};
int GridColumnSpace[4] = { 0, 0, 0, 0 };
int GridN_Rhs = 2;
int GridRhsSpace[2] = {0,0};
