// ==========================================================================
// instationary problem
// ==========================================================================

//===========================================================================
// example file
// =========================================================================

//==========================================================================
//three solid bodies (a slotted cylider, a hump and a cone) on a unit square,
//velocy field v=(0.5-y,x-0.5), counterclockwise rotation about the center
//of the unit square (0,1)X(0,1)
//==========================================================================
#define __ROTATING_BODIES__

void ExampleFile()
{
  OutPut("Example: RotatingBodies for level set" << endl);
  OutPut("Example: 2PhaseLevelSet.h" << endl);
}


// ========================================================================
// initial solution
// ========================================================================
void InitialU1(double x, double y, double *values)
{
  values[0] =  0.5 - y;
}

void InitialU2(double x, double y, double *values)
{
  values[0] = x - 0.5;
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
  cond = DIRICHLET;
}

void U1BoundValue(int BdComp, double Param, double &value)
{
//   double t = TDatabase::TimeDB->CURRENTTIME;
// 
//   switch(BdComp)
//   {
//     case 0: value=0;
//             break;
//     case 1: value=0;
//             break;
//     case 2: if(fabs(Param-0.5)<0.4999) 
//               value=1;
//             else
//               value=0;
//             break;
//     case 3: value=0;
//             break;
//     default: cout << "wrong boundary part number" << endl;
//             break;
//   }
  value = 0;
}

void U2BoundValue(int BdComp, double Param, double &value)
{
  value = 0;
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *X, double *Y,
               double **parameters, double **coeffs)
{
  int i;
  double *coeff, x, y; 
  static double nu=1./TDatabase::ParamDB->RE_NR;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    x = X[i];
    y = Y[i];

    coeff[0] = nu;

    coeff[1] = 0;  // f1
    coeff[2] = 0;  // f2
  }
}


// ========================================================================
// Routines for Level Set
// ========================================================================

// exact solution
void Exact_LS(double x, double y, double *values)
{
  double hump, conical, x0, y0;
  double t = TDatabase::TimeDB->CURRENTTIME;

  //rotation
  x0 = cos(t)*x+sin(t)*y-0.5*cos(t)-0.5*sin(t)+0.5;
  y0 = -sin(t)*x+cos(t)*y+0.5*sin(t)-0.5*cos(t)+0.5;

  values[0] = 0;
  //cylinder
  if((x0-0.5)*(x0-0.5)+(y0-0.75)*(y0-0.75)<=0.15*0.15)
  {
    if(fabs(x0-0.5) >= 0.025 || y0 >= 0.85)
      values[0] = 1;
  }
  //conical body
  conical = 1.0/0.15*sqrt((x0-0.5)*(x0-0.5)+(y0-0.25)*(y0-0.25));
  if(conical <= 1.0)
    values[0] = 1.0 - conical;
  //hump
  hump = 1.0/0.15*sqrt((x0-0.25)*(x0-0.25)+(y0-0.5)*(y0-0.5));
  if(hump <= 1.0)
    values[0] = 0.25*(1.0+(cos(Pi*hump)));

  values[1]=0;
  values[2]=0;
  values[3]=0;
}


// kind of boundary condition (for FE space needed)

void BoundCondition_LS(int BdComp, double t, BoundCond &cond)
{
  if (TDatabase::ParamDB->SOLD_PARAMETER_TYPE== FEM_FCT)
    cond = NEUMANN;
  else
    cond = DIRICHLET;
}


// value of boundary condition
void BoundValue_LS(int BdComp, double Param, double &value)
{
  value = 0;
}


// initial conditon
void InitialCondition_LS(double x,  double y, double *values)
{
  double hump, conical;

//   values[0] = 0;
//   //cylinder
//   if((x-0.5)*(x-0.5)+(y-0.75)*(y-0.75)<= 0.15*0.15)
//   {
//     if(fabs(x-0.5) >= 0.025 || y >= 0.85)
//       values[0] = 1.0;
//   }
//   //conical body
//   conical = 1.0/0.15*sqrt((x-0.5)*(x-0.5)+(y-0.25)*(y-0.25));
//   if(conical <= 1.0)
//     values[0] = 1.0 - conical;
//   //hump
//   hump = 1.0/0.15*sqrt((x-0.25)*(x-0.25)+(y-0.5)*(y-0.5));
//   if(hump <= 1.0)
//     values[0] = 0.25*(1.0+(cos(Pi*hump)));
//   
  
  
   //zero Levelset
   values[0] = 0.15*0.15 -  ((x-0.5)*(x-0.5)+(y-0.75)*(y-0.75)) ;  
  
}


void BilinearCoeffs_LS(int n_points, double *X, double *Y,
double **parameters, double **coeffs)
{
  double eps = 1.0/TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff;                                  // *param;
  double x, y;
  // double t = TDatabase::TimeDB->CURRENTTIME;

  if (TDatabase::ParamDB->SOLD_PARAMETER_TYPE== FEM_FCT)
    eps = 0.0;
  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    // param = parameters[i];

    x = X[i];
    y = Y[i];

    // diffusion
    coeff[0] = eps;
    // convection in x direction
    coeff[1] = 0.5 - y;
    // convection in y direction
    coeff[2] = x - 0.5;
    // reaction
    coeff[3] = 0;
    // rhs
    coeff[4] = 0;
    // rhs from previous time step
    coeff[5] = 0;
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

/****************************************************************/
/* computes errors in l1 and l2 norm                            */
/****************************************************************/
void ComputeDiscreteErrors(TFESpace2D *fespace,
TFEFunction2D *u,
double *sol,
double *errors,
double *lump_mass)
{
  int i,j,index, N_U, N_Cells, *GlobalNumbers, *BeginIndex, *DOF;
  double x, y, val[4], *err;
  double l1_error, l2_error;
  FE2D CurrentElement;
  TCollection *coll;
  TBaseCell *cell;

  // get arrays with the numbering of the dof
  GlobalNumbers = fespace->GetGlobalNumbers();
  BeginIndex = fespace->GetBeginIndex();
  N_U = fespace->GetN_DegreesOfFreedom();

  // get collection and number of cells
  coll = fespace->GetCollection();
  N_Cells = coll->GetN_Cells();

  // vector for error
  err = new double[N_U];

  // loop over cells
  for(i=0;i<N_Cells;i++)
  {
    // next cell
    cell = coll->GetCell(i);
    // pointer to global indices of dof connected with this cell
    DOF = GlobalNumbers + BeginIndex[i];
    // finds finite element on cell i
    CurrentElement = fespace->GetFE2D(i, cell);
    switch(CurrentElement)
    {
      // P_1
      case C_P1_2D_T_A:
	case C_UL1_2D_T_A:
        for (j=0;j<3;j++)
        {
          // global index
          index = DOF[j];
          // coordinates of vertex
          cell->GetVertex(j)->GetCoords(x, y);
          Exact_LS(x,y,val);
          err[index]=val[0] - sol[index];
        }
        break;
      case C_Q1_2D_Q_A:
      case C_Q1_2D_Q_M:
	case C_UL1_2D_Q_A:
	case C_UL1_2D_Q_M:
        for (j=0;j<4;j++)
        {
          // global index
          if(j<2)
            index = DOF[j];
          if(j==2)
            index = DOF[3];
          if(j==3)
            index = DOF[2];
          // coordinates of vertex
          cell->GetVertex(j)->GetCoords(x, y);
          Exact_LS(x,y,val);
          err[index] = val[0] - sol[index];
        }
        break;
      default:
        OutPut("ComputeDiscreteErrors not implemented for element "
          << CurrentElement << endl);
        exit(4711);
    }
  }

  l1_error = l2_error = 0.0;
  for (i=0;i<N_U;i++)
  {
    l1_error += lump_mass[i] * fabs(err[i]);
    l2_error += lump_mass[i] * err[i]*err[i];

  }
  delete err;

  errors[0] = l1_error;
  errors[1] = sqrt(l2_error);
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

//**************************************************************
//  UltraLocalErrors_TCD
//  computes errors in the vertices of the mesh cells
//  for ultra local projection schemes for some TCD problems
//**************************************************************

void UltraLocalError_TCD(TFEFunction2D *uh, 
        double *val, double *lumpmass)
{
  int i,j, N_Cells, N_Edges;
  double *Values, min_uh = 1e10, max_uh = 1e-10, x, y;
  double val_uh[4], val_u[4];
  TCollection *coll;
  TFESpace2D *fespace;
  TBaseCell *cell;

  // get coefficients of fe function
  fespace = uh->GetFESpace2D();
  Values = uh->GetValues();
  // get collection
  coll = fespace->GetCollection();
  // get number of cells
  N_Cells = coll->GetN_Cells();

  // loop over the mesh cells
  for(i=0;i<N_Cells;i++)
  {
    // get cell
    cell = coll->GetCell(i);
    // number of edges == number of vertices
    N_Edges = cell->GetN_Edges();
    // loop over the vertices
    for (j=0;j< N_Edges;j++)
    {
	// get coordinates of the edges 
	cell->GetVertex(j)->GetCoords(x, y);
	// compute local function value
	uh->FindGradientLocal(cell,i ,x, y, val_uh);
	// check for maximum and minimum
	if (val_uh[0] > max_uh)
	    max_uh = val_uh[0];
	if (val_uh[0] < min_uh)
	    min_uh = val_uh[0];
	// compute exact solutoin 
	Exact_LS(x,y,val_u);
	// error in this node
	//err[index]= val_u[0] - val_uh[0];
    }	  
  }
  val[2] = min_uh;
  val[3] = max_uh;
} // UltraLocalProjection_TCD
