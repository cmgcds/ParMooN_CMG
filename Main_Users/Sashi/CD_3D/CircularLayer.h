// ======================================================================
// three boundary layer problem
// ======================================================================
#define __CIRCULAR_LAYER__

#include <ConvDiff3D.h>
#include <Joint.h>
#include <BoundFace.h>
#include <BoundComp.h>
#include <FE3D.h>
#include <FEDesc3D.h>

void ExampleFile()
{
  OutPut("Example: CircularLayer.h" << endl) ;
}

// exact solution
void Exact(double x, double y, double z, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}

// kind of boundary condition (for FE space needed)
void BoundCondition(double x, double y, double z, BoundCond &cond)
{
   if( (fabs(x)<1e-6) || (fabs(y)<1e-6) || (fabs(z-1)<1e-6))
      cond = DIRICHLET;
   else
      cond = NEUMANN;
}

// value of boundary condition
void BoundValue(double x, double y, double z, double &value)
{
   if ((x-0.25)*(x-0.25) + y*y + (z-0.75)*(z-0.75) < 0.01)
      value = 1;
   else
      value = 0;
}

void BilinearCoeffs(int n_points, double *X, double *Y, double *Z,
        double **parameters, double **coeffs)
{
  static double eps=1/TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff, *param;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    coeff[0] = eps;
    coeff[1] = 0.5;
    coeff[2] = 0.75;
    coeff[3] = -0.75;
    coeff[4] = 0;
    coeff[5] = 0;
  }
}

/** calculate characteristic values */
void ComputeOutflowBoundary(int level, TFEFunction3D *ufct, double *errors)
{
  const int max_bound_points = 10000;
  int i,j,k, N_Cells;
  double xi, eta, zeta;
  TBaseCell *cell;
  TCollection *Coll;
  FE3D FE_ID;
  TFE3D *FE_Obj;
  RefTrans3D RefTrans;
  TBaseFunct3D *bf;
  int N_BaseFunct;
  TFESpace3D *FESpace3D;
  double *uorig, *uxorig, *uyorig, *uzorig;
  double *uref, *uxiref, *uetaref, *uzetaref, u1, x1, y1, z1;
  TJoint *joint;
  TBoundFace *boundface;
  TBoundComp3D *BoundComp;
/*  
  int *Numbers, comp, found, N_Faces;
  double u, ux, uy, uz, x, y, z, *Values;
  double val, x_min, val_min,  x_max, val_max, eps=1e-10;
  double y_coord[max_bound_points], uval[max_bound_points], min;
  int *GlobalNumbers, *BeginIndex, bound_points, index, first, first0;
  double vertx[4], verty[4];
  int vertices;

  FESpace3D = ufct->GetFESpace3D();
  BeginIndex = FESpace3D->GetBeginIndex();
  GlobalNumbers = FESpace3D->GetGlobalNumbers();
  Values = ufct->GetValues();  

  Coll = FESpace3D->GetCollection();
  N_Cells = Coll->GetN_Cells();
  val_min = 1e10;
  val_max = 1e-10;
  first0 =1 ;
  bound_points = 0;
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    N_Faces=cell->GetN_Faces();
    found = 0;
    for(j=0;j<N_Faces;j++)              // loop over all edges of cell
    {                                   // find edges on boundary part 3  
      joint=cell->GetJoint(j);          // this is x=0   
      if (joint->GetType() == BoundaryFace) // boundary face
      {        
        boundface = (TBoundFace *)joint;  
        BoundComp = boundface->GetBoundComp();  // get boundary component
        comp=BoundComp->GetID();              // boundary id 
        if (comp==0)
        {
          found = 1;
          break;
        }
      }
    }

    if (!found) continue;

    FE_ID = FESpace3D->GetFE3D(i, cell);
    FE_Obj = TFEDatabase3D::GetFE3D(FE_ID);
    RefTrans = FE_Obj->GetRefTransID();

    // get base function object
    bf = FE_Obj->GetBaseFunct3D();
    N_BaseFunct = bf->GetDimension();
    
    uorig = new double[N_BaseFunct];
    uxorig = new double[N_BaseFunct];
    uyorig = new double[N_BaseFunct];
    uzorig = new double[N_BaseFunct];
    
    uref = new double[N_BaseFunct];
    uxiref = new double[N_BaseFunct];
    uetaref = new double[N_BaseFunct];
    uzetaref = new double[N_BaseFunct];
    
    // set cell for reference transformation
    TFEDatabase3D::SetCellForRefTrans(cell, RefTrans);
    first = 1;
    vertices = 0;
    //OutPut("outflow " << level << endl);
    for (j=0;j<N_Faces;j++)
    {
      // check all vertices if they are on the boundary 
      z = cell->GetVertex(j)->GetZ();
      // point not on boundary
      if (fabs(z)>eps)
         continue;    
      y = cell->GetVertex(j)->GetY();
      x = cell->GetVertex(j)->GetX();
      vertx[vertices] = x;
      verty[vertices] = y;
      vertices++;
      // find local coordinates of the given point
      //cout << " x: " << x << endl;
      //cout << " y: " << y << endl;
      TFEDatabase3D::GetRefFromOrig(RefTrans, x, y, z, xi, eta, zeta);
      //cout << " xi: " << xi << endl;
      //cout << "eta: " << eta << endl;

      bf->GetDerivatives(D000, xi, eta, zeta, uref);
      bf->GetDerivatives(D100, xi, eta, zeta, uxiref);
      bf->GetDerivatives(D010, xi, eta, zeta, uetaref);
      bf->GetDerivatives(D001, xi, eta, zeta, uzetaref);
      
      TFEDatabase3D::GetOrigValues(RefTrans, xi, eta, zeta, N_BaseFunct,
                uref, uxiref, uetaref, uzetaref, uorig, uxorig, uyorig, uzorig);

      u = 0;
      Numbers = GlobalNumbers + BeginIndex[i];
      for(k=0;k<N_BaseFunct;k++)
      {
        val = Values[Numbers[k]];
        u += uorig[k]*val;
      }
      if (u<val_min)
	  val_min = u;
      if (u> val_max)
	  val_max = u;
      //OutPut("outflow " << level << " " << x << " " << y << " " << z << " " << u  << endl);
      if (first)
      {
         first = 0;
         x1 = x;
         y1 = y;
         z1 = z;
         u1 = u;
      }
      if (first0)
      {
         first0 = 0;
         x1 = x;
         y1 = y;
         z1 = z;
         u1 = u;
      }
    } // endfor faces;
    z = 0;
    // compute values for barycenter of edges and face
    // for (j=0;j<0;j++)
    for (j=0;j<=vertices;j++)
    {
	switch(vertices)
	{
	    case 3:
		switch(j)
		{
		    case 0:
			x = (vertx[0] + vertx[1])/2;
			y = (verty[0] + verty[1])/2;
			break;
		    case 1:
			x = (vertx[0] + vertx[2])/2;
			y = (verty[0] + verty[2])/2;
			break;
		    case 2:
			x = (vertx[2] + vertx[1])/2;
			y = (verty[2] + verty[1])/2;
			break;
		    case 3:
			x = (vertx[0] + vertx[1] + vertx[2])/3;
			y = (verty[0] + verty[1] + verty[2])/3;
			break;
		}
		break;
	    case 4:
		switch(j)
		{
		    case 0:
			x = (vertx[0] + vertx[1])/2;
			y = (verty[0] + verty[1])/2;
			break;
		    case 1:
			x = (vertx[1] + vertx[2])/2;
			y = (verty[1] + verty[2])/2;
			break;
		    case 2:
			x = (vertx[2] + vertx[3])/2;
			y = (verty[2] + verty[3])/2;
			break;
		    case 3:
			x = (vertx[3] + vertx[0])/2;
			y = (verty[3] + verty[0])/2;
			break;
		    case 4:
			x = (vertx[0] + vertx[1] + vertx[2]+ vertx[3])/4;
			y = (verty[0] + verty[1] + verty[2]+ verty[3])/4;
			break;
		}
		break;
	    default:
		OutPut("CircularLayer.h, wrong number of vertices" << endl);
		exit(4711);
	}
      // find local coordinates of the given point
      //cout << " x: " << x << endl;
      //cout << " y: " << y << endl;
      TFEDatabase3D::GetRefFromOrig(RefTrans, x, y, z, xi, eta, zeta);
      //cout << " xi: " << xi << endl;
      //cout << "eta: " << eta << endl;

      bf->GetDerivatives(D000, xi, eta, zeta, uref);
      bf->GetDerivatives(D100, xi, eta, zeta, uxiref);
      bf->GetDerivatives(D010, xi, eta, zeta, uetaref);
      bf->GetDerivatives(D001, xi, eta, zeta, uzetaref);
      
      TFEDatabase3D::GetOrigValues(RefTrans, xi, eta, zeta, N_BaseFunct,
                uref, uxiref, uetaref, uzetaref, uorig, uxorig, uyorig, uzorig);

      u = 0;
      Numbers = GlobalNumbers + BeginIndex[i];
      for(k=0;k<N_BaseFunct;k++)
      {
        val = Values[Numbers[k]];
        u += uorig[k]*val;
      }
      if (u<val_min)
	  val_min = u;
      if (u> val_max)
	  val_max = u;
    }
    //OutPut("outflow " << level << " " << x1 << " " << y1 << " " << z1 << " " << u1  << endl);
    delete uorig;
    delete uxorig;
    delete uyorig;
    delete uzorig;
    delete uref;
    delete uxiref;
    delete uetaref;
    delete uzetaref;
  } // endfor
  errors[0] = val_min;
  errors[1] = val_max;
  errors[2] = 0;
}

void ComputeOutflowBoundary_GetGrad(int level, TFEFunction3D *ufct, double *val)
{
  double h, z=0,values[4],y,x,min,max;    
  int i, j, bound_points = 128;
  h = 1.0/bound_points;

  min = 1e10;
  max = -1e10;

  for (i=0;i<=bound_points; i++)
  {
      y = i*h;
      OutPut(i<<endl);
      for (j=0; j <=bound_points; j++)
      {
	  x = j*h;
	  ufct->FindGradient(x,y,z,values);
	  if (values[0] > max)
	      max = values[0];
	  if (values[0] < min)
	      min = values[0];
      }
  }
  val[0] = min;
  val[1] = max;
  val[2] = 0.0;
*/
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

