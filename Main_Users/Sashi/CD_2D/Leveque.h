// ======================================================================
// Leveque problem
// ======================================================================
#define __LEVEQUE__

#include <ConvDiff2D.h>
#include <Joint.h>
#include <BoundEdge.h>
#include <BoundComp.h>
#include <FE2D.h>
#include <FEDesc2D.h>

void ExampleFile()
{
  OutPut("Example: Leveque.h" << endl) ;
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
void BoundCondition(int i, double Param, BoundCond &cond)
{
  switch(i)
    {
    case 0:
      if ((Param >= 0.1) && (Param <= 0.9))
	cond = DIRICHLET;
      else
	cond = NEUMANN;
      break;
    case 1:
      cond = NEUMANN;
      break;
    case 2:
      cond = NEUMANN;
      break;
    case 3:
      cond = DIRICHLET;
      break;
    }
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
	    value = 0;
	    break;
	case 2:
	  value = 0;
	  break;
	case 3:
	  value = 1;
	  break;
    }
}

void BilinearCoeffs(int n_points, double *X, double *Y,
        double **parameters, double **coeffs)
{
  double eps=1/TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff, *param;
  double x, y;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    param = parameters[i];
    y = Y[i];

    coeff[0] = 1;
    coeff[1] = 4 * y*(1-y)/eps;
    coeff[2] = 0;
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
   


// computation of some global errors, only for P1 or Q1 !!!
double ComputeSherwoodNumber(TFEFunction2D *ufct)
{
  int i, j, k, N_Cells, N_Edges, comp;
  double t0, t1, val[3], integral = 0, length, x0, x1, xp, value;
  TBaseCell *cell;
  TCollection *Coll;
  TFESpace2D *FESpace2D;
  TJoint *joint;
  TBoundEdge *boundedge;
  TBoundComp *BoundComp;
  
  // data for Gauss3 formula
  double w[]={0.2777777777777778,0.4444444444444444,0.2777777777777778};
  double x[]={0.11270166537925831149, 0.5, 0.88729833462074168851};
  for(int i=0;i<3;i++)
  {
    w[i] *= 2;
    x[i] = 2*x[i]-1;
  }

  FESpace2D = ufct->GetFESpace2D();
  Coll = FESpace2D->GetCollection();
  N_Cells = Coll->GetN_Cells();

  // loop over all edges
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    N_Edges=cell->GetN_Edges();
    // find boundary edge
    for(j=0;j<N_Edges;j++)
      {
        joint = cell->GetJoint(j);
	// boundary edge found
        if (joint->GetType() == BoundaryEdge)
        {
	  // get parameters
	  boundedge = (TBoundEdge *)joint;
	  BoundComp = boundedge->GetBoundComp();
	  boundedge->GetParameters(t0, t1);
          comp=BoundComp->GetID();
	  // not on the electrode
	  if ((comp > 0)|| (t0 < 0.1) || (t1 > 0.9))
	    continue;
	  // length of the edge
	  x0 = t0 * 10;
	  x1 = t1 * 10;
	  length = x1 - x0;
	  // compute integral on the edge
	  value = 0;
	  // apply quadrature rule
	  for (k=0;k<3;k++)
	    {
	      // compute point on edge
	      xp = x0 + (x[k]+1)*length/2.0;
	      // find gradient
	      ufct->FindGradientLocal(cell, i, xp, 0, val);
	      // update value
	      value += -val[2]*w[k];
	    }
	  integral += value * length/2.0;
	}  
      }
  }
  
  return(integral);
}

void CheckWrongNeumannNodes(TCollection *Coll, TFESpace2D *fespace,
			    int &N_neum_to_diri, int* &neum_to_diri,
			    int* &neum_to_diri_bdry, 
			    double* &neum_to_diri_param)
{
  const int max_entries = 8193;  
  int i, j, N_, min_val;
  int N_Cells, N_V, diri_counter = 0, found, diri_counter_1 = 0;
  int *global_numbers, *begin_index, *dof;
  int boundary_vertices[4], tmp_diri[max_entries], tmp_bdry[max_entries];
  double x[4], y[4], eps = 1e-6, tmp_param[max_entries];
  TBaseCell *cell;
  TVertex *vertex;
  FE2D CurrentElement;

  // number of mesh cells
  N_Cells = Coll->GetN_Cells();
  // array with global numbers of d.o.f.
  global_numbers = fespace->GetGlobalNumbers();
  // array which points to the beginning of the global numbers in
  // global_numbers for each mesh cell
  begin_index = fespace->GetBeginIndex();

  diri_counter = 0;
 
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    N_V = cell->GetN_Vertices();
    found = 0;
    for (j=0;j<N_V;j++)
    {
      // read coordinates of the mesh cell
      boundary_vertices[j] = 0;
      vertex = cell->GetVertex(j);
      vertex->GetCoords(x[j], y[j]);
      // vertex on the inlet
      if (fabs(x[j])<eps)
      {
	   boundary_vertices[j] = 1;
	   found++;
      }
      // vertex on the electrode
      if ((fabs(y[j])<eps) && (fabs(x[j]+eps)>=1) &&  (fabs(x[j]-eps)<=9))
      {
	boundary_vertices[j] = 1;
	found++;
      }
    }
    // no cell with edge with vertex on the boundary
    if (found<2) 
	continue;
    // finite element on the mesh cell
    CurrentElement = fespace->GetFE2D(i, cell);
    // number of basis functions (= number of d.o.f.)
    N_ = TFEDatabase2D::GetN_BaseFunctFromFE2D(CurrentElement);
    // the array which gives the mapping of the local to the global d.o.f.
    dof = global_numbers+begin_index[i];
    switch(CurrentElement)
    {
	// P_1, Q_1
	case C_P1_2D_T_A:
	case C_Q1_2D_Q_A:
	case C_Q1_2D_Q_M:
	    for (j=0;j<N_V;j++)
	    {
		// vertex on the boundary
		if (boundary_vertices[j])
		{
		    if (CurrentElement==C_P1_2D_T_A)
			tmp_diri[diri_counter] = dof[j];
		    else
		    {
			if (j<2)
			    tmp_diri[diri_counter] = dof[j];
			else
			{
			    if (j==2)
				tmp_diri[diri_counter] = dof[3];
			    else
				tmp_diri[diri_counter] = dof[2];
			}
		    }
		    if (diri_counter > max_entries)
		    {
			OutPut("tmp_diri too short !!!"<<endl);
			exit(4711);
		    }
		    if (fabs(y[j])<eps)
		    {
			tmp_bdry[diri_counter] = 0;
			tmp_param[diri_counter] = x[j];
		    }
		    if (fabs(1-y[j])<eps)
		    {
			tmp_bdry[diri_counter] = 2;
			tmp_param[diri_counter] = 1-x[j];
		    }
		    if (fabs(x[j])<eps)
		    {
			tmp_bdry[diri_counter] = 3;
			tmp_param[diri_counter] = 1-y[j];
		    }
		    if (fabs(1-x[j])<eps)
		    {
			tmp_bdry[diri_counter] = 1;
			tmp_param[diri_counter] = y[j];
		    }
		    diri_counter++;
		}
	    }
	    break;
	default:
	    OutPut("CheckNeumannNodesForVelocity not implemented for element "
		   << CurrentElement << endl);
	    OutPut("code can be run without CheckNeumannNodesForVelocity, just delete the exit" << endl);
	    exit(4711);
    }	    
  }
 
  // condense
  for (i=0;i<diri_counter;i++)
  {
      if (tmp_diri[i] == -1)
	  continue;
      diri_counter_1++;
      for (j=i+1;j<diri_counter;j++)
      {
	  if (tmp_diri[i] == tmp_diri[j])
	  {
	      tmp_diri[j] = -1;
	  }
      }
  }

  OutPut("CheckNeumannNodesForVelocity: N_neum_to_diri " << diri_counter_1 << endl);
  N_neum_to_diri = diri_counter_1;
  // allocate array for the indices
  neum_to_diri = new int[diri_counter_1];
  // allocate array for the corresponding boundary numbers
  neum_to_diri_bdry = new int[diri_counter_1];
  // allocate array for the corresponding boundary parameters
  neum_to_diri_param = new double[diri_counter_1];
  // fill array and sort
  for (i=0;i<diri_counter_1;i++)
  {
      min_val = tmp_diri[0];
      found = 0;
      for (j=1;j<diri_counter;j++)
      {
	  if ((tmp_diri[j]>0) && ((tmp_diri[j] < min_val) || 
				  (min_val == -1)))
	  {
	       min_val =  tmp_diri[j];
	       found = j;
	  }
      }
      neum_to_diri[i] = tmp_diri[found];
      neum_to_diri_bdry[i] = tmp_bdry[found];
      neum_to_diri_param[i] = tmp_param[found];
      tmp_diri[found] = -1;
  }

/*  for (i=0;i<diri_counter_1;i++)
  {
      OutPut(i << " " << neum_to_diri[i] << " " << neum_to_diri_bdry[i] <<
	     " " << neum_to_diri_param[i] <<  endl);
	     }*/
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
