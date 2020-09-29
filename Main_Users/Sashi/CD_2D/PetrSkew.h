// ======================================================================
// smooth solution problem
// ======================================================================
#define __TWO_INTERIOR_LAYERS__
#include <ConvDiff2D.h>
#include <Joint.h>
#include <BoundEdge.h>
#include <BoundComp.h>
#include <FE2D.h>
#include <FEDesc2D.h>

void ExampleFile()
{
  OutPut("Example: PetrSkew.h" << endl) ;
}
// exact solution
void Exact(double x, double y, double *values)
{

  values[0]= 0;
  values[1]= 0;
  values[2]= 0;
  values[3]= 0;
}

// kind of boundary condition (for FE space needed)
void BoundCondition(int i, double t, BoundCond &cond)
{
   cond = DIRICHLET;
}

// value of boundary condition
void BoundValue(int BdComp, double Param, double &value)
{
   switch (BdComp)
   {
      case 0:
         value = 1;
         break;
      case 1:
      case 2:   
         value = 0;
         break;
      case 3:
         if (Param<=0.8)
            value = 0;
         else
            value = 1;
         break;
   }
}

void BilinearCoeffs(int n_points, double *X, double *Y,
        double **parameters, double **coeffs)
{
  static double eps=1/TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff, *param;
  double t1,t3,t4,t5,t7,t9,t10,t12,t15,t30;
  double x, y;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    param = parameters[i];
    x = X[i];
    y = Y[i];

    coeff[0] = eps;
    coeff[1] = 2.0/sqrt(13.0);
    coeff[2] = 3.0/sqrt(13.0);
    coeff[3] = 0;
    coeff[4] = 0;
  }
}

void ComputeExtremalValues(int N, double *sol, double  *values)
{
   int i;
   double max, min;

   min = 1e10;
   max = -1e10;
   
   for(i=0;i<N;i++)
   {
      if(sol[i] > max)
         max = sol[i];
      if(sol[i] < min)
         min = sol[i];
   }

   values[0] = min;
   values[1] = max;
}
   
/** compute curve of the outflow boundary */
void ComputeOutflowBoundary(int level, TFEFunction2D *ufct)
{
  const int max_bound_points = 10000;
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
  double val, x_min, val_min,  x_max, val_max, eps=1e-10;
  double y_coord[max_bound_points], uval[max_bound_points], min;
  int *GlobalNumbers, *BeginIndex, bound_points, index;

  FESpace2D = ufct->GetFESpace2D();
  BeginIndex = FESpace2D->GetBeginIndex();
  GlobalNumbers = FESpace2D->GetGlobalNumbers();
  Values = ufct->GetValues();  

  Coll = FESpace2D->GetCollection();
  N_Cells = Coll->GetN_Cells();

  bound_points = 0;
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    N_Edges=cell->GetN_Edges();
    found = 0;
    for(j=0;j<N_Edges;j++)              // loop over all edges of cell
    {                                   // find edges on boundary part 3  
      joint=cell->GetJoint(j);          // this is x=0   
      if ((joint->GetType() == BoundaryEdge)||
          (joint->GetType() == IsoBoundEdge)) // boundary edge 
      {
        
        boundedge=(TBoundEdge *)joint;  
        BoundComp=boundedge->GetBoundComp();  // get boundary component
        comp=BoundComp->GetID();              // boundary id 
        if (comp==3)
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
    {
      // check all vertices if they are on the boundary 
      x = cell->GetVertex(j)->GetX();
      // point not on boundary
      if (fabs(x)>eps)
         continue;    
      y = cell->GetVertex(j)->GetY();
      //check if this node is already in the array of boundpoints
      found = 0;
      for (k=bound_points-1;k>=0; k--)
      {
         if (fabs(y- y_coord[k]) < eps)
         {
            found = 1;
            break;
         }
      }
      if (found)
         continue;
      // new node
      y_coord[bound_points] = y;
      bound_points++;      
      if ( bound_points > max_bound_points)
      {
         OutPut("PetrSkew.h: maximal number of boundary points reached !!!" << endl);
         exit(4711);
      }
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
      //OutPut(x << " " << y << " " << u << endl);
      uval[bound_points-1] = u;
    } // endfor 
    delete uorig;
    delete uxorig;
    delete uyorig;
    delete uref;
    delete uxiref;
    delete uetaref;
  } // endfor
  // order the arrays

  for (i=0;i<bound_points; i++)
  {
     min = 1e6;
     for(j=i;j<bound_points; j++)
     {
        if (y_coord[j]< min)
        {
           min = y_coord[j];
           index = j;
        }
     }
     // change the entries
     y_coord[index] = y_coord[i];
     y_coord[i] = min;
     val = uval[i];
     uval[i] = uval[index];
     uval[index] = val;
  }
  for (i=0;i<bound_points; i++)
  {
     OutPut("outflow " << level << " " <<  y_coord[i] << 
            " " <<  uval[i] << endl);
  }
}
/** compute curve of the outflow boundary */
void ComputeDiagonal(int level, TFEFunction2D *ufct)
{
  const int max_bound_points = 10000;
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
  
  int *Numbers, comp, found, N_Edges, N_add, diag =0;
  double u, ux, uy, x, y, *Values, x_add;
  double val, x_min, val_min,  x_max, val_max, eps=1e-10;
  double x_coord[max_bound_points], uval[max_bound_points], min;
  int *GlobalNumbers, *BeginIndex, bound_points, index;

  FESpace2D = ufct->GetFESpace2D();
  BeginIndex = FESpace2D->GetBeginIndex();
  GlobalNumbers = FESpace2D->GetGlobalNumbers();
  Values = ufct->GetValues();  

  Coll = FESpace2D->GetCollection();
  N_Cells = Coll->GetN_Cells();

  bound_points = 0;
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    N_Edges=cell->GetN_Edges();
    found = 0;

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
    x_add = 0;
    N_add = 0;
    diag = 0;
    for (j=0;j<=N_Edges;j++)
    {
      // check all vertices if they are on y = 0.5
      if (j <  N_Edges)
      {
         x = cell->GetVertex(j)->GetX();
         y = cell->GetVertex(j)->GetY();
         //check if this node is already in the array of boundpoints
         if (fabs(y-0.5)>eps)
            continue;
         x_add += x;
         N_add++;
      }
      else
      {
         // no mesh cell on y=0.5
         if (!diag)
            continue;
         // dof in the centre of the mesh cell 
         else
         {
            x = x_add/N_add;
            y = 0.5;
         }
      }
      diag++;
      found = 0;
      for (k=bound_points-1;k>=0; k--)
      {
         if (fabs(x - x_coord[k]) < eps)
         {
            found = 1;
            break;
         }
      }
      if (found)
         continue;
      // new node
      x_coord[bound_points] = x;
      bound_points++;      
      if ( bound_points > max_bound_points)
      {
         OutPut("PetrSkew.h: maximal number of boundary points reached !!!" << endl);
         exit(4711);
      }
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
      //OutPut(x << " " << y << " " << u << endl);
      uval[bound_points-1] = u;
    } // endfor 
    delete uorig;
    delete uxorig;
    delete uyorig;
    delete uref;
    delete uxiref;
    delete uetaref;
  } // endfor
  // order the arrays

  for (i=0;i<bound_points; i++)
  {
     min = 1e6;
     for(j=i;j<bound_points; j++)
     {
        if (x_coord[j]< min)
        {
           min = x_coord[j];
           index = j;
        }
     }
     // change the entries
     x_coord[index] = x_coord[i];
     x_coord[i] = min;
     val = uval[i];
     uval[i] = uval[index];
     uval[index] = val;
  }
  for (i=0;i<bound_points; i++)
  {
     OutPut("centerline " << level << " " <<  x_coord[i] << 
            " " <<  uval[i] << endl);
  }
}
