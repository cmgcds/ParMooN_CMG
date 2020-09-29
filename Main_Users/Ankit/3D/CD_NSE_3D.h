// channel with circular cross section



#include <InterfaceJoint3D.h>
#include <IsoInterfaceJoint3D.h>
#include <IsoBoundFace.h>
#include <MacroCell.h>
#include <BdSphere.h>
//#include <tetgen.h>


void ExampleFile()
{
  OutPut("Example: MHD3D.h" << endl) ;
  
  #define __Cylinder__   
}

// ========================================================================
// Boundary Condition For Navier Stokes
// ========================================================================
void InitialU1(double x, double y, double z, double *values)
{
  values[0] = 0;
}

void InitialU2(double x, double y, double z, double *values)
{
  values[0] = 0;
}

void InitialU3(double x, double y, double z, double *values)
{
  values[0] = 0;
}

void InitialP(double x, double y,  double z, double *values)
{
  values[0] = 0;
}

// exact solution
void ExactU1(double x, double y,  double z, double *values)
{

}

void ExactU2(double x, double y,  double z, double *values)
{

}

void ExactU3(double x, double y,  double z, double *values)
{
 

}

void ExactP(double x, double y,  double z, double *values)
{

}

// kind of boundary condition (for FE space needed)
void BoundCondition_NSE(int CompID, double x, double y, double z, BoundCond &cond)
{
  TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 0;
    if(CompID == 4)
        cond = NEUMANN;
    else
        cond  = DIRICHLET;
}

// value of boundary condition
void U1BoundValue(int CompID, double x, double y, double z, double &value)
{
   if (CompID == 2)
   {
      value = z*(4-z)*y*(4-y)/16;
   }
   else
   {
      value = 0;
   }
}

// value of boundary condition
void U2BoundValue(int CompID, double x, double y, double z, double &value)
{
  value = 0;
}

// value of boundary condition
void U3BoundValue(int CompID, double x, double y, double z, double &value)
{
  value = 0;
}

// coefficients for Stokes form: A, B1, B2, f1, f2
void LinCoeffs(int n_points, double *X, double *Y, double *Z,
               double **parameters, double **coeffs)
{
  static double eps = 1./TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff, x, y, z;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
      
    coeff[0] = eps;
    coeff[1] = 0;
    coeff[2] = 0;
    coeff[3] = 0;
  }
}

// ========================================================================
// Boundary Condition For Convection Diffusion
// ========================================================================
void InitialCondition(double x, double y, double z, double *values)
{ 
 double t = 0;
 double k = 0.1;
  
  values[0] = 0;
}

void BoundCondition_CD(int CompID,double x, double y, double z, BoundCond &cond)
{
    if(CompID == 2 || CompID == 3 || CompID == 5)
        cond = DIRICHLET;
    else
        cond = NEUMANN;
}

// value of boundary condition
void BoundValue(int BdComp,double x, double y, double z, double &value)
{ 
    if(BdComp == 3)
        value = 1;
    else if(BdComp > 9)
        cout << "CD wrong boundary part number" << endl;
    else
        value = 0;
}

// BilinearCoeffs for Heat 
void BilinearCoeffs(int n_points, double *X, double *Y, double *Z,
               double **parameters, double **coeffs)
{
//   double eps= 1.0/TDatabase::ParamDB->RE_NR;
  double Pe = TDatabase::ParamDB->RE_NR*Prandtl_Number;
  double eps=1/Pe;
  int i;
  double *coeff;                                  // *param;
  double *param;
  double c = 0;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    param = parameters[i];

    // diffusion
    coeff[0] = eps;
    if(TDatabase::ParamDB->INTERNAL_CONVECTION_EQ_VELOFIELD)
    {
        // convection in x direction
        coeff[1] = param[0];
        // convection in y direction
        coeff[2] = param[1];
        // convection in z direction
        coeff[3] = param[2];
    }

    // reaction
    coeff[4] = c;
     // rhs
    coeff[5] = 0;
    coeff[6] = 0;

  }
}
