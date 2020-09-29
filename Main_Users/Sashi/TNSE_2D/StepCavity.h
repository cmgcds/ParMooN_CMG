// time dependent Navier-Stokes problem 
// flow in a channel over a cavity and across a step
// 
// u(x,y) = ???
// p(x,y) = ???

//#define __CHANNELSTEP__

#include <Joint.h>
#include <BoundEdge.h>
#include <BoundComp.h>
#include <FE2D.h>
#include <FEDesc2D.h>

void ExampleFile()
{
  OutPut("Example: StepCavity.h" << endl) ;
}

// ========================================================================
// initial solution
// ========================================================================
void InitialU1(double x, double y, double *values)
{
  //values[0] = 1 ;
  // values[0] = y*(10-y)/25.0; 
  values[0] = 0 ;
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
  case 8:
     cond = NEUMANN;
     TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 0;
     break;
  case 10:
    cond = DIRICHLET;
    break;
  case 0:
  case 9:
    cond = SLIP_FRICTION_PENETRATION_RESISTANCE;
    TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION = 1;
    break; 
  default :
    //cond = SLIP_FRICTION_PENETRATION_RESISTANCE;
    //TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION = 1;
    cond = DIRICHLET;
    break;
  }
}

void U1BoundValue(int BdComp, double Param, double &value)
{
  switch(BdComp)
  {
  case 0: 
  case 1: 
  case 2: 
  case 3: 
  case 4: 
  case 5: 
  case 6: 
  case 7: 
  case 8: 
  case 9: 
    value=0;
    break;
  case 10:  
//    value = Param*(1-Param)*4.0;
    value = 1.0;
    break;
  default: cout << "wrong boundary part number" << endl;
    exit(4711);
  }
}

void U2BoundValue(int BdComp, double Param, double &value)
{
  value=0;
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *x, double *y,
               double **parameters, double **coeffs)
{
  static double eps = 1/TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    coeff[0] = eps;
    coeff[1] = 0; // f1
    coeff[2] = 0; // f2
  }
}

/** calculate reattachment point */
void GetReattachmentPoint(TFEFunction2D *u1fct,
                          double &reattachment_point)
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

  reattachment_point =0;

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
        if (comp==4)
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
      xv[j] = 0;
    for (j=0;j<N_Edges;j++)
    {
      y = cell->GetVertex(j)->GetY();
     // point not on boundary
      //if (fabs(y)>eps)
      //  continue;
      // for Dirichlet b.c.
      if (fabs(y)<eps)
        continue;
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
      
      TFEDatabase2D::GetOrigValues(RefTrans, xi, eta, bf, (TGridCell *)cell,
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
      if (xv[j] < 6.0) continue;
      for (k=j+1;k<N_Edges;k++)
      {
        if (xv[k] < 6.0) continue;
        if ((xv[j] < xv[k])&&(valv[j]<=0)&&(valv[k]>=0))
        { 
          if (reattachment_point>0)
            OutPut("FIRST reattachment point : " << reattachment_point << endl);
          reattachment_point = xv[j] - valv[j]*(xv[k]-xv[j])/(valv[k]-valv[j]);
       }
        if ((xv[j] > xv[k])&&(valv[j]>=0)&&(valv[k]<=0))
        { 
          if (reattachment_point>0)
            OutPut("FIRST reattachment point : " << reattachment_point << endl);
          reattachment_point = xv[k] - valv[k]*(xv[j]-xv[k])/(valv[j]-valv[k]);
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
