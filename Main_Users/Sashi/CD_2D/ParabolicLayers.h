// ======================================================================
// convection skew to the domain
// Hughes, Mallet, Mizukami 1986
// ======================================================================
#define __PARABOLIC_LAYERS__

#include <ConvDiff2D.h>
#include <Joint.h>
#include <BoundEdge.h>
#include <BoundComp.h>
#include <FE2D.h>
#include <FEDesc2D.h>
#include <FEDatabase2D.h>
#include <MainRoutines2D.h>

void ExampleFile()
{
  OutPut("Example: ParabolicLayers.h" << endl) ;
}
// exact solution (this is the solution for eps = 0,
// solution of reduced problem)
void Exact(double x, double y, double *values)
{
    values[0] = x;
    values[1] = 1;
    values[2] = 0;
    values[3] = 0;
}

// kind of boundary condition
void BoundCondition(int i, double t, BoundCond &cond)
{
    cond = DIRICHLET;
}

// value of boundary condition
void BoundValue(int BdComp, double Param, double &value)
{
    value = 0;
}

void BilinearCoeffs(int n_points, double *X, double *Y,
        double **parameters, double **coeffs)
{
  static double eps=1/TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff, *param;
  double x, y;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    param = parameters[i];
  
    coeff[0] = eps;
    coeff[1] = 1;
    coeff[2] = 0;
    coeff[3] = 0;
    coeff[4] = 1;
  }
}
// kind of boundary condition
void BoundConditionAdjoint(int i, double t, BoundCond &cond)
{
    cond = DIRICHLET;
}

// value of boundary condition
void BoundValueAdjoint(int BdComp, double Param, double &value)
{
    value = 0;
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
   values[1] = max - 1;
   values[2] = max - min;
}
   
/** compute solution at x=0.5 and y=0.5 */
void ComputeCutLines(int level, TFEFunction2D *ufct, double *values,
		     double meshsize)
{
  double h, x=1,val[3],y, max = 0, min = 0, x1 = 0, x2 = 0;
  double valref, yv[4], eps=1e-12;
  int i, j, N_Edges, bound_points, first = 1;
  int order = TDatabase::ParamDB->ANSATZ_ORDER;
  TCollection *Coll;
  TFESpace2D *FESpace2D;
  TBaseCell *cell;
  
  // compute distance of grid lines parallel to x--axis
  FESpace2D = ufct->GetFESpace2D();
  Coll = FESpace2D->GetCollection();
  cell = Coll->GetCell(0);
  N_Edges=cell->GetN_Edges();
  for(j=0;j<N_Edges;j++)              // loop over all edges of cell
  {
      yv[j] = cell->GetVertex(j)->GetY();
  }
  h = fabs(yv[0]-yv[1]);
  if (fabs(yv[1]-yv[2]) > h)
      h = fabs(yv[1]-yv[2]);
  if (N_Edges == 3)
  {
      if (fabs(yv[2]-yv[0]) > h)
	  h = fabs(yv[2]-yv[0]);
  }
  else
  {
      if (fabs(yv[2]-yv[3]) > h)
	  h = fabs(yv[2]-yv[3]);
      if (fabs(yv[3]-yv[0]) > h)
	  h = fabs(yv[3]-yv[0]);
  }

  OutPut(endl<<"ONLY FOR MESHES PARALLEL TO X-AXES !!!"<<endl);
  //bound_points =  (int)(sqrt(2.0)/meshsize+1e-4) + 1;
  //bound_points = 33;
  //h = 1.0/(bound_points-1);
  // for DG
  switch(order)
    {
    case -11:
    case -110:
      order = 1;
      break;
    case -12:
    case -120:
      order = 2;
      break;
    case -13:
    case -130:
      order = 3;
      break;
    case -14:
    case -140:
      order = 4;
      break;
    }
  h /= order;
  bound_points = (int) (1/h+1e-6);
  bound_points++;
  OutPut(order << " bound_points "<<  bound_points <<endl);
  OutPut(endl);

  // cut at x = 0.5;
  x = 0.5+eps;
  y = 0.5;
  ufct->FindGradient(x,y,val);
  valref = val[0];
  
  for (i=1;i<bound_points-1; i++)
  {
      y = i*h;
      ufct->FindGradient(x,y,val);
      OutPut("cut-x " << level << " " <<  y << " " <<  val[0] << endl);
      // max{u_h(0.5,y)-u_h(0.5,0.5)}
      if (val[0] - valref > max)
	  max = val[0] - valref;
      // min{u_h(0.5,y) - u_h(0.5,0.5)}
      if (val[0] - valref < min)
	  min = val[0] - valref;
  }
  values[0] = max;
  values[1] = min;
  // u_h(0.5,0.5) - 0.5
  values[2] = valref-0.5;

  // cut at y = 0.5; u_h(\bar x,0.5-\bar x)
  x = 1-h;
  y = 0.5;
  OutPut("last point " << 1.0/h << " " << x << " " << y);
  ufct->FindGradient(x,y,val); 
  values[3] = val[0] - x;
  OutPut(" " << x << " " << val[0] << endl);
}

// computation of some global errors, only for P1 or Q1 !!!
void ComputeGlobalErrors(TFEFunction2D *ufct, double *values)
{
  TBaseCell *cell;
  TCollection *Coll;
  TFESpace2D *FESpace2D;
  TJoint *joint;
  RefTrans2D RefTrans;
  TBaseFunct2D *bf;
  FE2D FE_ID;
  TFE2D *FE_Obj;
  TBoundEdge *boundedge;
  TBoundComp *BoundComp;
  int N_BaseFunct, comp;
  double xi, eta, eps = 1e-6, edgeval[4];
  double *uorig, *uxorig, *uyorig, *uref, *uxiref, *uetaref, u;
  double linfty = 0, l2 = 0, l2interior = 0;
  double *Values, fevalue[4], area, local, x, y, val;
  int *GlobalNumbers, *BeginIndex, index, N_Cells, N_Edges;
  int i, j, k, boundary_cell, *Numbers;

  FESpace2D = ufct->GetFESpace2D();
  BeginIndex = FESpace2D->GetBeginIndex();
  GlobalNumbers = FESpace2D->GetGlobalNumbers();
  Values = ufct->GetValues();  

  Coll = FESpace2D->GetCollection();
  N_Cells = Coll->GetN_Cells();

  // loop over all edges
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    N_Edges=cell->GetN_Edges();
    area = cell->GetMeasure();
    boundary_cell = 0;
    Numbers = GlobalNumbers + BeginIndex[i];
    for(j=0;j<N_Edges;j++)              // loop over all edges of cell
    {                                  
	/* joint=cell->GetJoint(j);         
      if ((joint->GetType() == BoundaryEdge)||
          (joint->GetType() == IsoBoundEdge)) // boundary edge 
      {
        boundedge=(TBoundEdge *)joint;  
        BoundComp=boundedge->GetBoundComp();  // get boundary component
        comp=BoundComp->GetID();              // boundary id 
	if (comp!=3)
	    boundary_cell = 1;
	    }*/
      fevalue[j] = 0;
    }

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
      // compute coordinates
      x = cell->GetVertex(j)->GetX();
      y = cell->GetVertex(j)->GetY();
      // cell with node at boundary
      if (fabs(y)<eps || fabs(1-y) < eps || fabs(1-x) < eps)
	  boundary_cell = 1;
	  
      // find local coordinates of the given point
      //cout << " x: " << x << endl;
      //cout << " y: " << y << endl;
      TFEDatabase2D::GetRefFromOrig(RefTrans, x, y, xi, eta);
      //cout << " xi: " << xi << endl;
      //cout << "eta: " << eta << endl;

      bf->GetDerivatives(D00, xi, eta, uref);
      bf->GetDerivatives(D10, xi, eta, uxiref);
      bf->GetDerivatives(D01, xi, eta, uetaref);
  
      TFEDatabase2D::GetOrigValues(RefTrans, xi, eta, bf, Coll, (TGridCell *) cell,
                uref, uxiref, uetaref, uorig, uxorig, uyorig);
      // compute value of fe function at (x,y)
      u = 0;
      for(k=0;k<N_BaseFunct;k++)
      {
        val = Values[Numbers[k]];
        u += uorig[k]*val;
      }
      //OutPut(x << " " << y << " " << u << endl);
      // update global errors, do not consider boundary points
      if ((fabs(u-x) > linfty)&& fabs(y)>eps && fabs(1-y) > eps && fabs(1-x) > eps)
      {
	  linfty = fabs(u-x);OutPut(x << " " << y << " " <<  fabs(u-x) << endl);
	  //OutPut(x << " " << y << endl);
      }
      if (fabs(y)>eps && fabs(1-y) > eps && fabs(1-x) > eps)
      {
	  if (N_Edges==3)
	      l2 += (u-x)*(u-x)/6.0; // six mesh cells at every interior node
	  else
	      l2 += (u-x)*(u-x)/4.0; // four mesh cells at every interior node
      }
      
      fevalue[j] = u-x;
    } // endfor (j) N_Edges

    // l2 error in interior mesh cells
    if (!boundary_cell)
    {
	local = 0;
	if (N_Edges == 3)
	{
	    edgeval[0] = (fevalue[0] + fevalue[1])/2.0;
	    edgeval[1] = (fevalue[1] + fevalue[2])/2.0;
	    edgeval[2] = (fevalue[2] + fevalue[0])/2.0;
	}
	else
	{
	    edgeval[0] = (fevalue[0] + fevalue[1])/2.0;
	    edgeval[1] = (fevalue[1] + fevalue[2])/2.0;
	    edgeval[2] = (fevalue[2] + fevalue[3])/2.0;
	    edgeval[3] = (fevalue[3] + fevalue[0])/2.0;
	}

	for (j=0;j<N_Edges;j++)
	    local+= edgeval[j]*edgeval[j];
	l2interior += area*local/N_Edges;
    }

    delete uorig;
    delete uxorig;
    delete uyorig;
    delete uref;
    delete uxiref;
    delete uetaref;
  } // endfor
    
  values[0] = linfty;
  values[1] = sqrt(l2);
  values[2] = sqrt(l2interior);
 }
// computation of some global errors, only for P1 or Q1 !!!
void ComputeGlobalErrorsIrr(TFEFunction2D *ufct, double *values)
{
  TBaseCell *cell;
  TCollection *Coll;
  TFESpace2D *FESpace2D;
  TJoint *joint;
  RefTrans2D RefTrans;
  TBaseFunct2D *bf;
  FE2D FE_ID;
  TFE2D *FE_Obj;
  TBoundEdge *boundedge;
  TBoundComp *BoundComp;
  int N_BaseFunct, comp;
  double xi, eta, eps = 1e-6, edgeval[4];
  double *uorig, *uxorig, *uyorig, *uref, *uxiref, *uetaref, u, ux, uy;
  double para_val = 0, para_der = 0, exp_der = 0, l2interior = 0;
  double *Values, fevalue[4], area, diam, local, x, y, val, *diffs;
  int *GlobalNumbers, *BeginIndex, index, N_Cells, N_Edges;
  int i, j, k, boundary_cell, *Numbers, N_dof, *dofs, current_dof;

  FESpace2D = ufct->GetFESpace2D();
  BeginIndex = FESpace2D->GetBeginIndex();
  GlobalNumbers = FESpace2D->GetGlobalNumbers();
  Values = ufct->GetValues();  
  N_dof = FESpace2D->GetActiveBound();
 
  diffs = new double[N_dof];
  memset(diffs,0,N_dof*SizeOfDouble);
  dofs = new int[N_dof];
  memset(dofs,0,N_dof*SizeOfInt);

  Coll = FESpace2D->GetCollection();
  N_Cells = Coll->GetN_Cells();

  // loop over all edges
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    N_Edges=cell->GetN_Edges();
    diam = cell->GetDiameter();
    if (N_Edges==3)
	area = diam*diam / 4.0;
    else
	area = diam*diam/2.0;
    // compute index of global d.o.f.
    Numbers = GlobalNumbers + BeginIndex[i];    
    boundary_cell = 0;
   
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
      // compute coordinates
      x = cell->GetVertex(j)->GetX();
      y = cell->GetVertex(j)->GetY();
      // cell with node at boundary
      if (fabs(y)<eps || fabs(1-y) < eps || fabs(1-x) < eps)
	  boundary_cell = 1;
	  
      // find local coordinates of the given point
      //cout << " x: " << x << endl;
      //cout << " y: " << y << endl;
      TFEDatabase2D::GetRefFromOrig(RefTrans, x, y, xi, eta);
      //cout << " xi: " << xi << endl;
      //cout << "eta: " << eta << endl;

      bf->GetDerivatives(D00, xi, eta, uref);
      bf->GetDerivatives(D10, xi, eta, uxiref);
      bf->GetDerivatives(D01, xi, eta, uetaref);
      
      TFEDatabase2D::GetOrigValues(RefTrans, xi, eta, bf, Coll, (TGridCell *)cell,
                uref, uxiref, uetaref, uorig, uxorig, uyorig);
      // compute value of fe function at (x,y)
      u = ux = uy = 0;

      for(k=0;k<N_BaseFunct;k++)
      {
        val = Values[Numbers[k]];
        u += uorig[k]*val;
        ux+= uxorig[k]*val;
        uy+= uyorig[k]*val;
     }
      //OutPut(x << " " << y << " " << u << endl);
      // update errors in parabolic layer, do not consider boundary points
      if ((x > 1e-6)&&(x<0.9)&&(y>1e-6)&&(y <= 0.1))
      {
        if (u-x > para_val)	
	   para_val = u-x;
        if (-uy > para_der)
            para_der = -uy; 
      }
      if ((x > 1e-6)&&(x<0.9)&&(y>=0.9)&&(y < 1-1e-6))
      {
         if (u-x > para_val)	
	   para_val = u-x;
        if (uy > para_der)
            para_der = uy; 
      }
      if ((x >= 0.9)&&(x<1-1e-6)&&(y>0.1)&&(y < 0.9))
      {
        if (ux > exp_der)
            exp_der = ux; 
      }
      // correspondence of local vertex to local dof
      if (N_Edges == 3)
         current_dof = Numbers[j];
      else
      {
         switch(j)
         {
            case 1: current_dof = Numbers[3];
                    break;
            case 3: current_dof = Numbers[1];
                    break;
            default:   current_dof = Numbers[j];
                    break;
          }
      } 
      if ((u<x)&&(current_dof<N_dof))
      {
         diffs[current_dof] += (u-x)*(u-x);
         dofs[current_dof]++;
     }
    } // endfor (j) N_Edges

    delete uorig;
    delete uxorig;
    delete uyorig;
    delete uref;
    delete uxiref;
    delete uetaref;
  } // endfor
  
  for (i=0;i<N_dof;i++)
  { 
    if (dofs[i] > 0){
      l2interior += diffs[i]/dofs[i];
      //OutPut(i << " " << dofs[i] << " " << diffs[i]  << " : ");
      }
  }

  delete diffs;
  delete dofs;

  values[0] = para_val;
  values[1] = para_der;
  values[2] = exp_der;
  values[3] = sqrt(l2interior);
 }

void CheckWrongNeumannNodes(TCollection *Coll, TFESpace2D *fespace,
			    int &N_neum_to_diri, int* &neum_to_diri,
			    int* &neum_to_diri_bdry, 
			    double* &neum_to_diri_param)
{
  CheckWrongNeumannNodesUnitSquareDiri(Coll, fespace,  N_neum_to_diri,
				       neum_to_diri, neum_to_diri_bdry,
				       neum_to_diri_param);
}

