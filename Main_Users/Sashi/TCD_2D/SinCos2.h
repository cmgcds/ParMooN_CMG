// ======================================================================
// instationary problem
// ======================================================================

#define __SINCOS1__

/// ========================================================================
// example file
// ========================================================================

void ExampleFile()
{
  OutPut("Example: SinCos2.h" << endl); 
  TDatabase::ParamDB->INTERNAL_STEADY_STATE_MATRICES_OR_RHS = 0;
}
// exact solution
void Exact(double x, double y, double *values)
{
  double t;

  t = TDatabase::TimeDB->CURRENTTIME;

  values[0] = sin(t*t)*cos(Pi*(x-0.5))*y*(1-y);
  values[1] = -sin(t*t)*Pi*sin(Pi*(x-0.5))*y*(1-y);
  values[2] = sin(t*t)*cos(Pi*(x-0.5))*(1-2*y);
  values[3] =  -sin(t*t)*Pi*Pi*cos(Pi*(x-0.5))*y*(1-y)-sin(t*t)*cos(Pi*(x-0.5))*2;
}

// kind of boundary condition (for FE space needed)
void BoundCondition(int BdComp, double t, BoundCond &cond)
{
    if (TDatabase::ParamDB->SOLD_PARAMETER_TYPE == FEM_FCT)
	cond = NEUMANN;
    else
	cond = DIRICHLET;
}

// value of boundary condition
void BoundValue(int BdComp, double Param, double &value)
{
  double t;
  
  value = 0.0;
}

// initial conditon
void InitialCondition(double x, double y, double *values)
{
  double t;

  t = TDatabase::TimeDB->CURRENTTIME;

  values[0] = sin(t*t)*cos(Pi*(x-0.5))*y*(1-y);
}

// kind of boundary condition (for FE space needed)
/*void BoundConditionNSE(int BdComp, double t, BoundCond &cond)
{
  cond = DIRICHLET;
}
*/

void BilinearCoeffs(int n_points, double *X, double *Y,
        double **parameters, double **coeffs)
{
  double eps=1/TDatabase::ParamDB->RE_NR;
  double a=2, b=-1, c=1;
  int i;
  double *coeff, *param;
  double x, y;
  double t = TDatabase::TimeDB->CURRENTTIME;
  double tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;

  // previous discrete time
  tau = t-tau;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    param = parameters[i];

    x = X[i];
    y = Y[i];

    coeff[0] = eps;
    coeff[1] = a;
    coeff[2] = b;
    coeff[3] = c;

    coeff[4] = cos(t*t)*2*t*cos(Pi*(x-0.5))*y*(1-y)
	- eps* ( -sin(t*t)*Pi*Pi*cos(Pi*(x-0.5))*y*(1-y)-sin(t*t)*cos(Pi*(x-0.5))*2) 
	- a*sin(t*t)*Pi*sin(Pi*(x-0.5))*y*(1-y) + b*sin(t*t)*cos(Pi*(x-0.5))*(1-2*y)
	+ c*sin(t*t)*cos(Pi*(x-0.5))*y*(1-y);
    // rhs from previous time step
    coeff[5] = cos(tau*tau)*2*tau*cos(Pi*(x-0.5))*y*(1-y)
	- eps* ( -sin(tau*tau)*Pi*Pi*cos(Pi*(x-0.5))*y*(1-y)-sin(tau*tau)*cos(Pi*(x-0.5))*2) 
	- a*sin(tau*tau)*Pi*sin(Pi*(x-0.5))*y*(1-y) + b*sin(tau*tau)*cos(Pi*(x-0.5))*(1-2*y)
	+ c*sin(tau*tau)*cos(Pi*(x-0.5))*y*(1-y);
  }
}

/****************************************************************/
/* finds the nodes which are Neumann and should be Dirichlet    */
/* for FEM_FCT schemes                                          */
/****************************************************************/

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
	    OutPut( tmp_diri[diri_counter] << " " <<
		    tmp_bdry[diri_counter] << " " << tmp_param[diri_counter] << endl);
            diri_counter++;
          }
        }
	OutPut(endl);
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
      if ((tmp_diri[j]>-1) && ((tmp_diri[j] < min_val) ||
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

