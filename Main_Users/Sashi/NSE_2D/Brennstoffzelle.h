// Navier-Stokes problem 
// backward facing step flow

#define __CHANNEL30__

#include <math.h>
#include <Constants.h>
#include <Joint.h>
#include <BoundEdge.h>
#include <BoundComp.h>
#include <FE2D.h>
#include <FEDesc2D.h>

// ========================================================================
// exact solution
// ========================================================================
// first component of the velocity u1
// values[0] = u1
// values[1] = \partial u1/ \partial x
// values[2] = \partial u1/ \partial y
// values[3] = Laplacian of u1
// if the functions are not known or if they are not needed, set them to 0

void ExampleFile()
{
  OutPut("Example: Brennstoffzelle.h" << endl) ;
}

void ExactU1(double x, double y, double *values)
{
  
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

// second component of the velocity u2
// values[0] = u2
// values[1] = \partial u2/ \partial x
// values[2] = \partial u2/ \partial y
// values[3] = Laplacian of u2
// if the functions are not known or if they are not needed, set them to 0

void ExactU2(double x, double y, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

// pressure p
// values[0] = p
// values[1] = \partial p/ \partial x
// values[2] = \partial p/ \partial y
// values[3] = Laplacian of p
// if the functions are not known or if they are not needed, set them to 0

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

// type of the boundary condition
// possible types:
// DIRICHLET, NEUMANN, SLIP_FRICTION_PENETRATION_RESISTANCE
void BoundCondition(int BdComp, double t, BoundCond &cond)
{
  // no-slip on the whole boundary
  switch(BdComp)
  {
    case 0: 
       cond = SLIP_FRICTION_PENETRATION_RESISTANCE;
      break;
    case 1: 
      cond = NEUMANN;
      TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 0;
      break;
    case 2: 
       cond = SLIP_FRICTION_PENETRATION_RESISTANCE;
     break;
    case 3: 
      if (t>0.5)
        cond = SLIP_FRICTION_PENETRATION_RESISTANCE;
      else
         cond = DIRICHLET;
      break;
  }
  TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION = 1;
}

// boundary values of u1
// counting the boundary parts starts at 0
void U1BoundValue(int BdComp, double Param, double &value)
{
  switch(BdComp)
  {
    case 0: 
      value=0;
      break;
    case 1: 
      value=0;
      break;
    case 2: 
      value=0;
      break;
    case 3: 
      // inflow boundary
      if (Param>=0.5)
        value=0;
      else
        value = 24*Param*(0.5-Param);
      break;
    default: cout << "wrong boundary part number" << endl;
            break;
  }
}

// boundary values of u2
void U2BoundValue(int BdComp, double Param, double &value)
{
  switch(BdComp)
  {
    case 0: value=0;
            break;
    case 1: value=0;
            break;
    case 2: value=0;
            break;
    case 3: value=0;
            break;
    default: cout << "wrong boundary part number" << endl;
            break;
  }
}

// ========================================================================
// coefficients for the Navier--Stokes equations
// viscosity \nu and right hand side f 
// ========================================================================
void LinCoeffs(int n_points, double *X, double *Y,
               double **parameters, double **coeffs)
{
  // set nu
  static double nu=1/TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff;

  // the coefficients are needed for a set of quadrature points
  // loop over all points
  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    
    // coeff[0] is the viscosity
    coeff[0] = nu;
    
    // coeff[1] is the first component of the right hand side 
    coeff[1] = 0;
    coeff[2]=  0;           
  }
}



/** calculate reattachment point */
void GetReattachmentPoint(TFEFunction2D *u1fct,
                          double *reattachment_point)
{
  int i,j,k, N_Cells;
  double xi, eta;
  TBaseCell *cell;
  TCollection *Coll;
  FE2D FE_ID;
  TFE2D *FE_Obj;
  RefTrans2D RefTrans;
  TBaseFunct2D *bf;
  int N_BaseFunct;
  TFESpace2D *FESpace2D;
  double *uorig, *uxorig, *uyorig, *uref, *uxiref, *uetaref, *u1;
  TJoint *joint;
  TBoundEdge *boundedge;
  TBoundComp *BoundComp;
  
  int *Numbers, comp, found, N_Edges;
  double u, ux, uy, x, y, *Values;
  double val, x_min, val_min,  x_max, val_max, eps=1e-12;
  double xv[4], valv[4];
  int *GlobalNumbers, *BeginIndex;

  FESpace2D = u1fct->GetFESpace2D();
  BeginIndex = FESpace2D->GetBeginIndex();
  GlobalNumbers = FESpace2D->GetGlobalNumbers();
  Values = u1fct->GetValues();  

  reattachment_point[0] =0;
  reattachment_point[1] =0;
  reattachment_point[2] =0;

  Coll = FESpace2D->GetCollection();
  N_Cells = Coll->GetN_Cells();

  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    N_Edges=cell->GetN_Edges();
    found = 0;
    for(j=0;j<N_Edges;j++)              // loop over all edges of cell
    {
      joint=cell->GetJoint(j);
      if ((joint->GetType() == BoundaryEdge)||
          (joint->GetType() == IsoBoundEdge)) // boundary edge 
      {
        
        boundedge=(TBoundEdge *)joint;  
        BoundComp=boundedge->GetBoundComp();  // get boundary component
        comp=BoundComp->GetID();              // boundary id 
        if ((comp==0)||(comp==2))
        {
          found = 1;
          break;
        }
      }
    }

    if (!found) continue;

    FE_ID = FESpace2D->GetFE2D(i, cell);
    FE_Obj = TFEDatabase2D::GetFE2D(FE_ID);
    RefTrans = FE_Obj->GetRefTransID();

    // get base function object
    bf = FE_Obj->GetBaseFunct2D();
    N_BaseFunct = bf->GetDimension();
    
    uorig = new double[N_BaseFunct];
    uxorig = new double[N_BaseFunct];
    uyorig = new double[N_BaseFunct];
    
    uref = new double[N_BaseFunct];
    uxiref = new double[N_BaseFunct];
    uetaref = new double[N_BaseFunct];
    
    // set cell for reference transformation
    TFEDatabase2D::SetCellForRefTrans(cell, RefTrans);
    for (j=0;j<N_Edges;j++)
      xv[j] = -1;
    for (j=0;j<N_Edges;j++)
    {
      y = cell->GetVertex(j)->GetY();
      // point not on boundary
      // fo`r slip b.c.
      if (comp==0)
      {
        if  (TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION)
        {
          if (fabs(y+0.5)>eps)
            continue;
        }
        else
        {      
          // for Dirichlet b.c.
          if (fabs(y+0.5)<eps)
            continue;
        }
      }
      else
      {
        if  (TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION)
        {
          if (fabs(y-0.5)>eps)
            continue;
        }
        else
        {      
          // for Dirichlet b.c.
          if (fabs(y-0.5)<eps)
            continue;
        }
      }

      x = cell->GetVertex(j)->GetX();
      xv[j] = x;
      // find local coordinates of the given point
      //cout << " x: " << x << endl;
      //cout << " y: " << y << endl;
      TFEDatabase2D::GetRefFromOrig(RefTrans, x, y, xi, eta);
      //cout << " xi: " << xi << endl;
      //cout << "eta: " << eta << endl;

      bf->GetDerivatives(D00, xi, eta, uref);
      bf->GetDerivatives(D10, xi, eta, uxiref);
      bf->GetDerivatives(D01, xi, eta, uetaref);
      
      TFEDatabase2D::GetOrigValues(RefTrans, xi, eta, N_BaseFunct,
                uref, uxiref, uetaref, uorig, uxorig, uyorig);

      u = 0;
      Numbers = GlobalNumbers + BeginIndex[i];
      for(k=0;k<N_BaseFunct;k++)
      {
        val = Values[Numbers[k]];
        u += uorig[k]*val;
      }
      OutPut(x << " " << y << " " << u << endl);
      valv[j] = u;
     } // endfor 
    for (j=0;j<N_Edges;j++)
    {
      // vertex not on boudary
      if ((xv[j] < 0)||(xv[j]>25))
        continue;
      for (k=j+1;k<N_Edges;k++)
      {
        // vertex not on boudary
        if ((xv[k] < 0)||(xv[k]>25))
          continue;
        if (comp==0)
        {
          if ((xv[j] < xv[k])&&(valv[j]<=0)&&(valv[k]>=0))
          { 
            if (reattachment_point[0]>0)
              OutPut("FIRST lower reattachment point a: " << reattachment_point[0] << endl);
            reattachment_point[0] = xv[j] - valv[j]*(xv[k]-xv[j])/(valv[k]-valv[j]);
          }
          if ((xv[j] > xv[k])&&(valv[j]>=0)&&(valv[k]<=0))
          { 
            if (reattachment_point[0]>0)
              OutPut("FIRST lower reattachment point b: " << reattachment_point[0] << endl);
            reattachment_point[0] = xv[k] - valv[k]*(xv[j]-xv[k])/(valv[j]-valv[k]);
          }
        }
        else
        {
          if ((xv[j] < xv[k])&&(valv[j]>=0)&&(valv[k]<=0))
          { 
            if (reattachment_point[1]>0)
              OutPut("FIRST upper reattachment point : " << reattachment_point[1] << endl);
            reattachment_point[1] = xv[j] - valv[j]*(xv[k]-xv[j])/(valv[k]-valv[j]);
          }
          if ((xv[j] > xv[k])&&(valv[j]<=0)&&(valv[k]>=0))
          { 
            if (reattachment_point[1]>0)
              OutPut("FIRST upper reattachment point : " << reattachment_point[1] << endl);
            reattachment_point[1] = xv[k] - valv[k]*(xv[j]-xv[k])/(valv[j]-valv[k]);
          }
            if ((xv[j] < xv[k])&&(valv[j]<=0)&&(valv[k]>=0))
          { 
            if (reattachment_point[2]>0)
              OutPut("SECOND upper reattachment point : " << reattachment_point[2] << endl);
            reattachment_point[2] = xv[j] - valv[j]*(xv[k]-xv[j])/(valv[k]-valv[j]);
          }
          if ((xv[j] > xv[k])&&(valv[j]>=0)&&(valv[k]<=0))
          { 
            if (reattachment_point[2]>0)
              OutPut("SECOND upper reattachment point : " << reattachment_point[2] << endl);
            reattachment_point[2] = xv[k] - valv[k]*(xv[j]-xv[k])/(valv[j]-valv[k]);
          }
        }
      }
    }
    delete uorig;
    delete uxorig;
    delete uyorig;
    delete uref;
    delete uxiref;
    delete uetaref;

  } // endfor

}
