// ======================================================================
// convection skew to the domain
// Hughes, Mallet, Mizukami 1986
// ======================================================================
#define __HMM_1986__

#include <ConvDiff2D.h>
#include <Joint.h>
#include <BoundEdge.h>
#include <BoundComp.h>
#include <FE2D.h>
#include <FEDesc2D.h>

void ExampleFile()
{
  OutPut("Example: HMM1986.h" << endl) ;
}
// exact solution (this is the solution for eps = 0)
void Exact(double x, double y, double *values)
{
    values[0] = 0;
    values[1] = 0;
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
    switch(BdComp)
    {
	case 0:
	case 1:
	    value = 0;
	    break;
	case 2:
	    if (Param < 1e-6)
		value = 0;
	    else
		value = 1;
	    break;
	case 3:
	    if (Param<0.3)
		value = 1;
	    else
		value = 0;
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

  arg = -Pi/3;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    param = parameters[i];
    x = X[i];
    y = Y[i];

    coeff[0] = eps;
    coeff[1] = cos(arg);
    coeff[2] = sin(arg);
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
   values[1] = max - 1;
   values[2] = max - min;
}
// kind of boundary condition
void BoundConditionAdjoint(int i, double t, BoundCond &cond)
{
    cond = DIRICHLET;
}
void BoundValueAdjoint(int BdComp, double Param, double &value)
{
  value = 0;
}

/** compute solution at x=0.5 and y=0.5 */
void ComputeCutLines(int level, TFEFunction2D *ufct, double *values,
		     double meshsize)
{
  double h, x=1,val[3],y, max = 0, x1 = 0, x2 = 0, size;
  double intlow = 0, intupp = 0, c, x3;
  int i, bound_points = 100001, first = 1, start, end;
  
  h = 1.0/(bound_points-1);
  size = meshsize/sqrt(2);
  //OutPut("size "<< 1.0/size <<endl);
  //OutPut(endl);

  // smearing of the boundary layer at y = 0.25
  y = 0.25;
  start = (int) (0.175/h);
  end = (int) (0.4/h);
  for (i=start;i<=end; i++)
  {
      x = i*h;
      ufct->FindGradient(x,y,val);
      // OutPut("cut-y " << level << " " <<  x << " " <<  val[0] << endl);
      if ((val[0] >= 0.1)&&(x1==0))
	 x1 = x; 
      if ((val[0] >= 0.9)&&(x2==0))
      {
	 x2 = x;
	 break;
      }
  }

  start = (int) (0.74/h);
  end = (int) (1/h);
  end = start-1;
  for (i=start;i<=end; i++)
  {
      x = i*h;
      ufct->FindGradient(x,y,val);
       if ((x>=0.75)&&(x<=1-size))
      {
	  c = 1;
	  // first integration point
	  if (first)
	  {
	      c = 0.5;
	      first = 0;
	  }
	  // last integration point
	  if (x+h>1-size)
	      c = 0.5;
	  if (val[0] - 1 > 0)
	      intupp += (val[0]-1)*c;
	  else
	      intlow += (1-val[0])*c;
      }
  }
  
  start = (int) (0.9/h);
  end = (int) (1/h);
  start = end+1;
  for (i=end;i>=start; i--)
  {
      x = i*h;
      ufct->FindGradient(x,y,val);
      if (val[0] >= 0.9)
      {
	  x3 = x;
	  break;
      }
  }

  x = 0.8;
  y = 1.0/64.0;
  ufct->FindGradient(x,y,val);
  x3 = val[0];
  x = 63.0/64.0;
  y = 0.5;
  ufct->FindGradient(x,y,val);
  if (val[0] < x3)
      x3 = 1 - val[0];
  else
      x3 = 1 - x3;
  
  //OutPut(x1 << " " << x2 << " " << x3 << endl);
  values[0] = x2 - x1;
  //values[1] = intupp * h;
  //values[2] = intlow * h;
  //values[3] = x3;  
  //values[1] = x3;
  //OutPut("intsmear " << x2 - x1 << " ");
  
}

// computation of some global errors, only for P1 or Q1 !!!
void ComputeGlobalErrors(TFEFunction2D *ufct, double *values,
			 int N_neum_to_diri,
			 int *neum_to_diri)
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
  double u;
  double linftyexp = 0, l2exp = 0, l2int = 0, linftyint = 0, linftyintm = 0;
  double linftysmear = 0, l2smear = 0, linftyint_sold2 = 0;
  double *Values, fevalue[4], x, y, val, x0, y0;
  double *right_layer, *right_smear, *int_layer;
  int *GlobalNumbers, *BeginIndex, index, N_Cells, N_Edges, found;
  int i, j, k, *Numbers, nvert, N_U, dof, N_Active, vert_out;
  
  FESpace2D = ufct->GetFESpace2D();
  BeginIndex = FESpace2D->GetBeginIndex();
  GlobalNumbers = FESpace2D->GetGlobalNumbers();
  N_Active = FESpace2D->GetActiveBound();
  N_Active = FESpace2D->GetHangingBound();
  Values = ufct->GetValues();  
  N_U = ufct->GetLength();

  Coll = FESpace2D->GetCollection();
  N_Cells = Coll->GetN_Cells();

  right_layer = new double[3*N_U];
  right_smear = right_layer + N_U;
  int_layer = right_smear + N_U;
  memset(right_layer,0,3*N_U*SizeOfDouble);

  // loop over all edges
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    N_Edges=cell->GetN_Edges();
    
    x = 0;
    y = 0;
    // compute bary centre of mesh cell
    vert_out = 0;
    for(j=0;j<N_Edges;j++)              // loop over all edges of cell
    {                                  
      x0 = cell->GetVertex(j)->GetX();
      y0 = cell->GetVertex(j)->GetY();
      if ((x0>0.5)&&(x0<0.7))
	  vert_out++;
      if ((x0<=0.5)&&(y0<0.1))
	  vert_out++;
      x += x0;
      y += y0;
    }
    if (vert_out)
	continue;
    x /= N_Edges;
    y /= N_Edges;

    FE_ID = FESpace2D->GetFE2D(i, cell);
    FE_Obj = TFEDatabase2D::GetFE2D(FE_ID);

    // get base function object
    bf = FE_Obj->GetBaseFunct2D();
    N_BaseFunct = bf->GetDimension();

    // compute index of global d.o.f.
    Numbers = GlobalNumbers + BeginIndex[i];    

    for(j=0;j<N_BaseFunct;j++)              // loop over all base functions
    {
	dof = Numbers[j];
	// no Dirichlet nodes
	if (dof>=N_Active)
	    continue;
	found = 0;
	for (k=0;k<N_neum_to_diri;k++)
	{
	    if (dof == neum_to_diri[k])
	    {
		found++;
		break;
	    }
	}
	if (found)
	    continue;
	u = Values[dof];
	// interior layer
	if ((x<=0.5)&&(y>=0.1))
	{
	    if (u>1)
	    {
		//OutPut("o " << x << " " << y << " " << u-1 << endl); 
		int_layer[dof] = u-1;
		if (u-1>linftyint)
		  linftyint = u-1;
	    }
	    if (u<0)
	    {
		//OutPut("u " << x << " " << y << " " << u << endl); 
		int_layer[dof] = u;
		if (-u > linftyintm)
		    linftyintm = -u;
	    }
	}
	if ((x<=0.5)&&(y>=0.25))
	{
	    if (u>1)
	    {
		//OutPut("o " << x << " " << y << " " << u-1 << endl); 
		if (u-1>linftyint_sold2)
		  linftyint_sold2 = u-1;
	    }
	    if (u<0)
	    {
		//OutPut("u " << x << " " << y << " " << u << endl); 
		if (-u > linftyint_sold2)
		    linftyint_sold2 = -u;
	    }
	}
	// oscillations at exponential layer
	if (x>=0.7)
	{
	    if (u>1)
	    {
                //if (u-1>1e-8)
		    //OutPut("exp " << x << " " << y << " " << u-1 << endl); 
		right_layer[dof] = (u-1);
		if (u-1>linftyexp)
		    linftyexp = u-1;
	    }
	}
	// smearing at exponential layers
	if (x>=0.7)
	{
	    if (u<1)
	    {
	      right_smear[dof] = (1-u);
	      //OutPut(x << " " << y << " " << u << " "<< dof << endl);
	      if (1-u >  linftysmear)
	      {
		  linftysmear = 1-u;
	      }
	    }
      }
    } // endfor (j) N_BaseFunct
  } // endfor
  
/*  OutPut("l2int " << sqrt(l2int) << " liint " << linftyint+linftyintm 
	 << " l2exp " << sqrt(l2exp) << " liexp " << linftyexp
	 << " l2smear " << sqrt(l2smear) << " lismear " << linftysmear << endl);

  OutPut(" linftyint " << linftyint << " linftyintm " << linftyintm <<
	 " linftyexp " << linftyexp << " " << " linftyint_sold2 "
      << linftyint_sold2 << " " << endl);
*/
  values[1] = sqrt(Ddot(N_Active,int_layer,int_layer));
  values[2] = sqrt(Ddot(N_Active,right_layer,right_layer));
  values[3] =  sqrt(Ddot(N_Active,right_smear,right_smear));
  values[4] = linftyint_sold2;
  values[5] = linftyexp;
  delete right_layer;
 }

// computation of some global errors, only for P1 or Q1 !!!
void ComputeScores(double *value, double *score)
{
    if (value[1]<=1e-4)
	score[1] = 4;
    else
    {
	if (value[1] <= 2e-4)
	    score[1] = 3;
	else
	{
	    if (value[1] <= 1e-2)
		score[1] = 2;
	    else
	    {
		if  (value[1] <= 2e-2)
		    score[1] = 1;
		else
		{
		    if  (value[1] <= 1e-1)
			score[1] = 0;
		    else
		    {
			if  (value[1] <= 2e-1)
			    score[1] = -2;
			else
			    score[1] = -4;
		    }
		}
	    }
	}
    }

    if (value[2]<=1e-5)
	score[2] = 4;
    else
    {
	if (value[2] <= 2e-5)
	    score[2] = 3;
	else
	{
	    if (value[2] <= 1e-3)
		score[2] = 2;
	    else
	    {
		if  (value[2] <= 2e-3)
		    score[2] = 1;
		else
		{
		    if  (value[2] <= 1e-1)
			score[2] = 0;
		    else
		    {
			if  (value[2] <= 2e-1)
			    score[2] = -2;
			else
			    score[2] = -4;
		    }
		}
	    }
	}
    }

    if (value[0]<=4e-2)
	score[0] = 2;
    else
    {
	if (value[0] <= 5e-2)
	    score[0] = 1.5;
	else
	{
	    if (value[0] <= 6e-2)
		score[0] = 1;
	    else
	    {
		if  (value[0] <= 6.2e-2)
		    score[0] = 0.5;
		else
		{
		    if  (value[0] <= 8e-2)
			score[0] = 0;
		    else
		    {
			if  (value[0] <= 8.2e-2)
			    score[0] = -1;
			else
			    score[0] = -2;
		    }
		}
	    }
	}
    }
    if (value[3]<=1e-4)
	score[3] = 2;
    else
    {
	if (value[3] <= 1e-3)
	    score[3] = 1.5;
	else
	{
	    if (value[3] <= 1e-2)
		score[3] = 1;
	    else
	    {
		if  (value[3] <= 3e-2)
		    score[3] = 0.5;
		else
		{
		    if  (value[3] <= 5e-1)
			score[3] = 0;
		    else
		    {
			if  (value[3] <= 7e-1)
			    score[3] = -1;
			else
			    score[3] = -2;
		    }
		}
	    }
	}
    }
    return;
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

