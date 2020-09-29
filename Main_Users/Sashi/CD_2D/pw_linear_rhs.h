// ======================================================================
// Example with piecewise linear right hand side
// ======================================================================
#define __PW_LINEAR_RHS__
#include <ConvDiff2D.h>
void ExampleFile()
{
  OutPut("Example: pw_linear_rhs.h" << endl) ;
  TDatabase::ParamDB->INTERNAL_QUAD_RULE = 98;
  TDatabase::ParamDB->INTERNAL_PROBLEM_IDENTITY = 100;
}

// reduced solution
void Exact(double x, double y, double *values)
{
    if (fabs(x-0.5)>0.25 || fabs(y-0.5)>0.25)
        values[0] = 0;
    else
        values[0] = -16.*(x-0.25)*(x-0.75);
    if (fabs(x-0.5)>0.25 || fabs(y-0.5)>0.25)
        values[1] = 0;
    else
        values[1] = -32.*(x-0.5);
  values[2] = 0;
  values[3] = 0;
}

void BoundCondition(int BdComp, double t, BoundCond &cond)
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
  double eps=1/TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff, *param, x,y;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    x = X[i];
    y = Y[i];

    coeff[0] = eps;
    coeff[1] = 1;
    coeff[2] = 0;
    coeff[3] = 0;
    if (fabs(x-0.5)>0.25 || fabs(y-0.5)>0.25)
        coeff[4] = 0;
    else
        coeff[4] = -32.*(x-0.5);
  }
}

// computation of some global errors, only for P1 or Q1 !!!
void ComputeLocalExtrema(TFEFunction2D *ufct, double *values)
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
  double *Values, fevalue[4], area, diam, local, x, y, val;
  int *GlobalNumbers, *BeginIndex, index, N_Cells, N_Edges;
  int i, j, k, boundary_cell, *Numbers;
  double extr[4];

  extr[0] = -1;
  extr[1] = -1;
  extr[2] = -1;

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
    diam = cell->GetDiameter();
    if (N_Edges==3)
	area = diam*diam / 4.0;
    else
	area = diam*diam/2.0;
    for(j=0;j<N_Edges;j++)              // loop over all edges of cell
    {                                  
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
      if (x<-1.5)
	  continue;
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
      // compute value of fe function at (x,y)
      u = 0;
      Numbers = GlobalNumbers + BeginIndex[i];
      for(k=0;k<N_BaseFunct;k++)
      {
        val = Values[Numbers[k]];
        u += uorig[k]*val;
      }
      //OutPut(x << " " << y << " " << u << endl);
      // strip with the circle
      if ((fabs(x)>=0.4)&&(fabs(x)<=0.6))
      {
	  if ((u<=0)&&(fabs(u)>extr[0]))
	      extr[0] = fabs(u);
     }
      if (fabs(x)>=0.8)
      {
	  if ((u<=0)&&(fabs(u)>extr[1]))
	      extr[1] = fabs(u);
	  if ((u>=0)&&(fabs(u)>extr[2]))
	      extr[2] = fabs(u);
      }
    } // endfor (j) N_Edges
  
    delete uorig;
    delete uxorig;
    delete uyorig;
    delete uref;
    delete uxiref;
    delete uetaref;
  } // endfor
    
  //OutPut( extr[1] << " " <<  extr[2] << endl);
  values[0] = extr[0];
  values[1] = extr[1]+extr[2];
 }
/** compute solution at x=0.5 */
void ComputeCutLines(int level, TFEFunction2D *ufct, double *values,
		     double meshsize)
{
  double h, x=1,val[3],y, max = 0, x1 = 0, x2 = 0, size;
  double intlow = 0, intupp = 0, c, x3;
  int i, bound_points = 10001, first = 1, start, end;
  
  h = 1.0/(bound_points-1);
  size = meshsize/sqrt(2);
  OutPut("size "<< 1.0/size <<endl);
  //OutPut(endl);

  x = 0.5;
  x1 = 0;
  start = (int) (0.0/size);
  end = (int) (1.0/size);
  for (i=start;i<=end; i++)
  {
      y = i*size;
      ufct->FindGradient(x,y,val);
      //OutPut("cut-y " << level << " " <<  y << " " <<  val[0] << endl);
      if (val[0] < x1)
	 x1 = val[0]; 
  }
  OutPut(setprecision(7) << " " << x1 <<endl);
  values[0] = x1;
  return;
}

void CheckWrongNeumannNodes(TCollection *Coll, TFESpace2D *fespace,
			    int &N_neum_to_diri, int* &neum_to_diri,
			    int* &neum_to_diri_bdry, 
			    double* &neum_to_diri_param)
{
    const int max_entries = 4096;  
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
      // vertex on the upper lid
      if ((fabs(x[j])<eps)||(fabs(y[j])<eps)||(fabs(1-x[j])<eps)||(fabs(1-y[j])<eps))
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

  for (i=0;i<diri_counter_1;i++)
  {
      OutPut(i << " " << neum_to_diri[i] << " " << neum_to_diri_bdry[i] <<
	     " " << neum_to_diri_param[i] <<  endl);
  }
}

