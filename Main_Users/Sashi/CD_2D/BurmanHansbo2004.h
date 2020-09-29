// ======================================================================
// nonhomogeneous problem
// Burman, Hansbo, CMAME 2004
// ======================================================================
#define __BURMAN_HANSBO_2004__

#include <ConvDiff2D.h>
#include <Joint.h>
#include <BoundEdge.h>
#include <BoundComp.h>
#include <FE2D.h>
#include <FEDesc2D.h>

void ExampleFile()
{
  OutPut("Example: BurmanHansbo2004.h" << endl) ;
}
// exact solution (this is the solution for eps = 0)
void Exact(double x, double y, double *values)
{

  if (y<1-1e-7)
  {
      values[0] = exp(-5*x/(1-y*y));
      values[1] = -5*exp(-5*x/(1-y*y))/(1-y*y);
      values[2] = -10*x*y*exp(-5*x/(1-y*y))/((1-y*y)*(1-y*y));
  }
  else
  {
      values[0] = 0;
      values[1] = 0;
      values[2] = 0;
  }
  values[3]= 0;
}

// kind of boundary condition
void BoundCondition(int i, double t, BoundCond &cond)
{
      cond = DIRICHLET;
}

// value of boundary condition
void BoundValue(int BdComp, double Param, double &value)
{
    switch(BdComp)
    {
	case 0:
	    value = 0;
	    break;
	case 1:
	    if (Param<0.8)
		value = 0;
	    else
		value = 1;
	    break;
	case 2:
	case 3:
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
  double x, y, arg;

  arg = -125*Pi/180;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    param = parameters[i];
    x = X[i];
    y = Y[i];

    coeff[0] = eps;
    coeff[1] = cos(arg);
    coeff[2] = sin(arg);
    //OutPut(coeff[1] << " " << coeff[2] << " ");
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
   
/** compute solution at x=0.5 and y=0.5 */
void ComputeCutLines(int level, TFEFunction2D *ufct)
{
  double h, x=1,values[3],y, max = 0;
  int i, bound_points = 201;
  h = 1.0/(bound_points-1);
  for (i=0;i<bound_points; i++)
  {
      x = 0.5;
      y = i*h;
      ufct->FindGradient(x,y,values);
      if (values[0] > max)
	  max = values[0];
     OutPut("cut-x " << level << " " <<  y << 
            " " <<  values[0] << endl);
   }

  for (i=0;i<=bound_points; i++)
  {
      y = 0.5;
      x = i*h;
      ufct->FindGradient(x,y,values);
      if (values[0] > max)
	  max = values[0];
      OutPut("cut-y " << level << " " <<  x << 
            " " <<  values[0] << endl);
   }

     
}
/** compute cuts only in points of the mesh */
void ComputeCutsFE(int level, TFEFunction2D *ufct)
{
  const int max_bound_points = 10001;
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
  return;
  int *Numbers, comp, found, N_Edges, N_add, diag =0;
  double u, ux, uy, x, y, *Values, x_add;
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
      // check all vertices if they are on the boundary 
      if (j <  N_Edges)
      {
         x = cell->GetVertex(j)->GetX();
         y = cell->GetVertex(j)->GetY();
         //check if this node is already in the array of boundpoints
         if (fabs(x-0.5)>eps)
            continue;
         x_add += x;
         N_add++;
      }
      else
      {
         // no mesh cell on x=0.5
         if (!diag)
            continue;
         else
         {
            x = x_add/N_add;
            y = x;
         }
      }
      diag++;
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
         OutPut("TwoInteriorLeyers.h: maximal number of boundary points reached !!!" << endl);
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
     OutPut("diagonal " << level << " " <<  y_coord[i]*sqrt(2.0) << 
            " " <<  uval[i] << endl);
  }
}
