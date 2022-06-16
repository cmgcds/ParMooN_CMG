// Navier-Stokes, channel flow


void ExampleFile()
{
  OutPut("Example: Channel with slip BC" << endl) ;
}

#define __BENCH__


#include <Joint.h>
#include <BoundEdge.h>
#include <BoundComp.h>
#include <FE2D.h>
#include <FEDesc2D.h>

#include <MainUtilities.h>


#include <MacroCell.h>
#include <BoundEdge.h>
#include <IsoBoundEdge.h>
#include <gridgen.h>
#include <IsoInterfaceJoint.h>
#include <BdLine.h>
#include <BdCircle.h>
#include <GridCell.h>

#include <QuadAffin.h>
#include <QuadBilinear.h>
#include <QuadIsoparametric.h>
#include <TriaAffin.h>
#include <TriaIsoparametric.h>


// ========================================================================
// initial solution
// ========================================================================
void InitialU1(double x, double y, double *values)
{
  values[0] = 0;
}

void InitialU2(double x, double y, double *values)
{
  values[0] = 0;
}

void InitialP(double x, double y, double *values)
{
  values[0] = 0;
}
// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

void ExactU2(double x, double y, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

void ExactP(double x, double y, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

// ========================================================================
// boundary conditions
// ========================================================================
void BoundCondition(int i, double t, BoundCond &cond)
{
  switch(i)
  {
      case 2:  case 0:  //case 2: Walls
      cond = NEUMANN;
      
	  break;

   case 3:  //case 2: // case 3 : // case 4: 
  // TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION = 1;
  cond = DIRICHLET;
  break;

    case 1:  
	    cond = NEUMANN;
    break;
    
  }
  TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 0;
}

void U1BoundValue(int BdComp, double Param, double &value)
{
  // Param = Param + 0.5;
  
  switch(BdComp)
  {
    case 0: value=0.0; //bottom
            break;
   case 3: value=  01*(((Param)*(1-Param))/0.25); //1*(((Param)*(1-Param))/0.25);//1.5*(1.0 - (pow(0.5*Param,2))/0.25); // inlet
  //  value = 1;
// case 3:
//     if(Param > 0.005 and Param < 0.995){
//       value =5;
//     }
//     else
//       value = 0;
// value = 1;
            break;
    case 2: value = 0;  //upper wall
            break;
    case 1: value =0;  // outlet
            break;

    default: cout << "wrong boundary part number: " << BdComp << endl;
  }
}

void U2BoundValue(int BdComp, double Param, double &value)
{
  value = 0;
  if(BdComp>5) cout << "wrong boundary part number: " << BdComp << endl;

}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *x, double *y,
               double **parameters, double **coeffs)
{
  double eps = 1/TDatabase::ParamDB->RE_NR;
  double t = TDatabase::TimeDB->CURRENTTIME;
  int i;
  double *coeff;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    coeff[0] = eps;
    coeff[1] = 0; // 2*sin(x[i]); // f1
    coeff[2] = 0; // f2
  }
}


