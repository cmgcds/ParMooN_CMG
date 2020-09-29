//
// u(x,y) = unknown
// p(x,y) = unknown

#define __BULK_ACAD_TEST__

void ExampleFile()
{
  int range;

  OutPut(" Example: Bulk_Fallung_Driven_Cavity_Academic.h ") ;
  TDatabase::ParamDB->BULK_u_infty = 1e-2;
  OutPut(" u_infty " << TDatabase::ParamDB->BULK_u_infty);

  TDatabase::ParamDB->BULK_D_P_0 = 1e-9;
  TDatabase::ParamDB->BULK_D_P_MAX = 1e-4;
  TDatabase::ParamDB->BULK_c_C_infty_sat = 1e-4;
  TDatabase::ParamDB->BULK_C_g = 50;

  TDatabase::ParamDB->BULK_D_P_MIN = TDatabase::ParamDB->BULK_D_P_0/TDatabase::ParamDB->BULK_D_P_MAX;
  TDatabase::ParamDB->BULK_c_C_infty = TDatabase::ParamDB->BULK_c_C_infty_sat
      * exp(TDatabase::ParamDB->BULK_C_2/TDatabase::ParamDB->BULK_D_P_0);
  TDatabase::ParamDB->BULK_f_infty = TDatabase::ParamDB->BULK_u_infty/(TDatabase::ParamDB->BULK_C_g
      *TDatabase::ParamDB->BULK_k_g*pow(TDatabase::ParamDB->BULK_D_P_MAX,3)*TDatabase::ParamDB->BULK_l_infty);
  OutPut("BULK d_p_min " << TDatabase::ParamDB->BULK_D_P_MIN  <<
	 " u_infty " << TDatabase::ParamDB->BULK_u_infty <<
	 " c_C_infty " << TDatabase::ParamDB->BULK_c_C_infty <<
	 " f_infty " << TDatabase::ParamDB->BULK_f_infty << endl);

  TDatabase::ParamDB->SAVE_DATA = TRUE;
}
// ========================================================================
// initial solution
// ========================================================================
void InitialU1(double x, double y, double *values)
{
  values[0] = 0;
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
    cond = DIRICHLET;
}

void U1BoundValue(int BdComp, double Param, double &value)
{
    value = 0;
}

void U2BoundValue(int BdComp, double Param, double &value)
{
  value = 0;
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *x, double *y,
               double **parameters, double **coeffs)
{
  int i;
  double *coeff;
  double l_infty = TDatabase::ParamDB->BULK_l_infty;
  double u_infty = TDatabase::ParamDB->BULK_u_infty;
  double density = TDatabase::ParamDB->BULK_density;
  double dynvisc = TDatabase::ParamDB->BULK_dynamic_viscosity;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    coeff[0] = dynvisc/(density*l_infty*u_infty);
    //coeff[0] = 1;
    coeff[1] = 0; // f1
    coeff[2] = 0; // f2
  }
}


// ========================================================================
// definitions for the substance A
// ========================================================================

// initial conditon
void InitialCondition_c_A(double x, double y, double *values)
{
    values[0] = 0;
}

// kind of boundary condition (for FE space needed)
void BoundCondition_c_A(int BdComp, double Param, BoundCond &cond)
{
    cond = NEUMANN;
}

// value of boundary condition
void BoundValue_c_A(int BdComp, double Param, double &value)
{
    value = 0;
}

void NoCoeffs(int n_points, double *X, double *Y,
		    double **parameters, double **coeffs)
{
    return;
}

void BilinearCoeffs(int n_points, double *X, double *Y,
		    double **parameters, double **coeffs)
{
  int i;
  double *coeff, *param;
  double x, y;
  double t = TDatabase::TimeDB->CURRENTTIME;
  double L_infty = TDatabase::ParamDB->BULK_l_infty;
  double U_infty = TDatabase::ParamDB->BULK_u_infty;
  double C_infty = TDatabase::ParamDB->BULK_c_infty;
  double D_A = TDatabase::ParamDB->BULK_D_A;
  double k_r = TDatabase::ParamDB->BULK_k_r;
  double T_infty, eps, c;
  
  T_infty = L_infty/U_infty;
  eps = D_A/(L_infty*U_infty);
  c = k_r*C_infty * L_infty /U_infty;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    param = parameters[i];

    coeff[0] = eps;
    coeff[1] = 0;
    coeff[2] = 0;
    coeff[3] = 0;
    coeff[4] = 0; 
    //OutPut(param[0] << " " << param[1] << " " << param[2] << endl); 
  }
}

/****************************************************************/
/* finds the nodes which are Neumann and should be Dirichlet    */
/* for FEM_FCT schemes                                          */
/****************************************************************/

void CheckWrongNeumannNodes_c_A(TCollection *Coll, TFESpace2D *fespace,
int &N_neum_to_diri, int* &neum_to_diri,
int* &neum_to_diri_bdry,
double* &neum_to_diri_param)
{
  const int max_entries = 4096;
  int i, j, N_, min_val;
  int N_Cells, N_V, diri_counter = 0, found, diri_counter_1 = 0;
  int *global_numbers, *begin_index, *dof;
  int boundary_vertices[4], tmp_diri[max_entries], tmp_bdry[max_entries];
  double x[4], y[4], eps = 1e-8, tmp_param[max_entries];
  TBaseCell *cell;
  TVertex *vertex;
  FE2D CurrentElement;
  int range = (int)TDatabase::ParamDB->P7;

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
      if ((fabs(x[j])<eps) && ((y[j]>=range/32.0-eps)&&(y[j]<=(range+1)/32.0+eps)))
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
            if (fabs(x[j])<eps)
            {
              tmp_bdry[diri_counter] = 3;
              tmp_param[diri_counter] = 1-y[j];
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

  OutPut("CheckNeumannNodesForVelocity c_A: N_neum_to_diri " << diri_counter_1 << endl);
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

// ========================================================================
// definitions for the substance B
// ========================================================================

// initial conditon
void InitialCondition_c_B(double x, double y, double *values)
{
     values[0] = 0;
}

// kind of boundary condition (for FE space needed)
void BoundCondition_c_B(int BdComp, double Param, BoundCond &cond)
{
    cond = NEUMANN;
}

// value of boundary condition
void BoundValue_c_B(int BdComp, double Param, double &value)
{
    value = 0;
}

/****************************************************************/
/* finds the nodes which are Neumann and should be Dirichlet    */
/* for FEM_FCT schemes                                          */
/****************************************************************/

void CheckWrongNeumannNodes_c_B(TCollection *Coll, TFESpace2D *fespace,
int &N_neum_to_diri, int* &neum_to_diri,
int* &neum_to_diri_bdry,
double* &neum_to_diri_param)
{
  const int max_entries = 4096;
  int i, j, N_, min_val;
  int N_Cells, N_V, diri_counter = 0, found, diri_counter_1 = 0;
  int *global_numbers, *begin_index, *dof;
  int boundary_vertices[4], tmp_diri[max_entries], tmp_bdry[max_entries];
  double x[4], y[4], eps = 1e-8, tmp_param[max_entries];
  TBaseCell *cell;
  TVertex *vertex;
  FE2D CurrentElement;
  int range = (int)TDatabase::ParamDB->P8;

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
      if ((fabs(1.0-x[j])<eps) && ((y[j]>=range/32.0-eps)&&(y[j]<=(range+1)/32.0+eps)))
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
            if ((fabs(1.0-x[j])<eps))
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

  OutPut("CheckNeumannNodesForVelocity c_B: N_neum_to_diri " << diri_counter_1 << endl);
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

// ========================================================================
// definitions for the substance C
// ========================================================================

// exact solution
void Exact_c_C(double x, double y, double *values)
{
  double t = TDatabase::TimeDB->CURRENTTIME;

  values[0] = 0.6*(sin(2*Pi*x)*sin(2*Pi*y)+1)*sin(Pi*t/2)*t;
  values[1] = 0.6*(2*Pi*cos(2*Pi*x)*sin(2*Pi*y))*sin(Pi*t/2)*t;
  values[2] = 0.6*(2*Pi*sin(2*Pi*x)*cos(2*Pi*y))*sin(Pi*t/2)*t;
  values[3] = 0;
}

// initial conditon
void InitialCondition_c_C(double x, double y, double *values)
{
  double t = TDatabase::TimeDB->CURRENTTIME;
  values[0] = 0.6*(sin(2*Pi*x)*sin(2*Pi*y)+1)*sin(Pi*t/2)*t;	
}

// kind of boundary condition (for FE space needed)
void BoundCondition_c_C(int BdComp, double Param, BoundCond &cond)
{
   cond = NEUMANN;
}

// value of boundary condition
void BoundValue_c_C(int BdComp, double Param, double &value)
{
  double t = TDatabase::TimeDB->CURRENTTIME;
  value = 0.6*sin(Pi*t/2)*t;
}

// param[2] : c_A
// param[3] : c_B
// param[4] : c_C (old)
// param[5] : r_g
void BilinearCoeffs_Cc(int n_points, double *X, double *Y,
		    double **parameters, double **coeffs)
{
  int i;
  double *coeff, *param;
  double x, y, rhs, c_C;
  double t = TDatabase::TimeDB->CURRENTTIME;
  double L_infty = TDatabase::ParamDB->BULK_l_infty;
  double U_infty = TDatabase::ParamDB->BULK_u_infty;
  double C_infty = TDatabase::ParamDB->BULK_c_infty;
  double c_C_infty_sat =  TDatabase::ParamDB->BULK_c_C_infty_sat;
  double C_g = TDatabase::ParamDB->BULK_C_g;
  double C_nuc = 15.33;
  double C_2 = TDatabase::ParamDB->BULK_C_2;
  double D_A = TDatabase::ParamDB->BULK_D_A;
  double d_p_0 = TDatabase::ParamDB->BULK_D_P_0;
  double d_p_max = TDatabase::ParamDB->BULK_D_P_MAX;
  double k_g = TDatabase::ParamDB->BULK_k_g;
  double k_r = TDatabase::ParamDB->BULK_k_r;
  double k_nuc =  TDatabase::ParamDB->BULK_k_nuc;
  double eps,  B_C_c, T_infty, lambda_chem, lambda_nuc;
  double d_p_min = TDatabase::ParamDB->BULK_D_P_MIN;
  double c_C_infty = TDatabase::ParamDB->BULK_c_C_infty;
  double r_g;

  // compute derived quantities of the model
  T_infty = L_infty/U_infty;
  // compute coefficients of the equation
  eps = D_A/(L_infty*U_infty);
  lambda_chem = k_r*C_infty*C_infty*L_infty /(U_infty*c_C_infty);
  lambda_nuc = C_nuc*k_nuc*d_p_0*d_p_0*d_p_0*L_infty*pow(c_C_infty,4)/U_infty;


  eps = 1e-7;
  lambda_chem = 0;
  lambda_nuc = 5e-4;
  

  //OutPut("lambda_chem " << lambda_chem << "  lambda_nuc " << lambda_nuc << " " << endl);
  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    param = parameters[i];

    x = X[i];
    y = Y[i];

    if (param[4] < 1)
      B_C_c = 0;
    else
      B_C_c = pow(param[4] - 1,5);

    coeff[0] = eps;
    coeff[1] = 2*(2*y-1)*(1-(2*x-1)*(2*x-1)) * (sin(TDatabase::ParamDB->P1*Pi*t)+1)*TDatabase::ParamDB->P9;
    coeff[2] = -2*(2*x-1)*(1-(2*y-1)*(2*y-1)) *(sin(TDatabase::ParamDB->P1*Pi*t)+1)*TDatabase::ParamDB->P9;
    coeff[3] = 0;
    // following script from 05/08/23 for r_g
    r_g = param[4] - 1e-3;
    // rhs from academic
    // c_t
    rhs = 0.6*(sin(2*Pi*x)*sin(2*Pi*y)+1)*(cos(Pi*t/2) * Pi/2*t+sin(Pi*t/2));
    // u \cdot \nabla c
    rhs += coeff[1] * 0.6*(2*Pi*cos(2*Pi*x)*sin(2*Pi*y))*sin(Pi*t/2)*t;
    rhs += coeff[2] * 0.6*(2*Pi*sin(2*Pi*x)*cos(2*Pi*y))*sin(Pi*t/2)*t;
    // diffusion
    rhs += - coeff[0] * 0.6*(-4*Pi*Pi*sin(2*Pi*x)*sin(2*Pi*y))*sin(Pi*t/2)*t;
    rhs += - coeff[0] *  0.6*(-4*Pi*Pi*sin(2*Pi*x)*sin(2*Pi*y))*sin(Pi*t/2)*t;
    // nucleation
    c_C = 0.6*(sin(2*Pi*x)*sin(2*Pi*y)+1)*sin(Pi*t/2)*t;
    if (c_C > 1)
	rhs += lambda_nuc * pow(c_C - 1,5);
    // integral term
    rhs += (c_C-1e-3) 
	*1e3*(2.0/3.0-1.0/6.0-0.25-(2.0*pow(d_p_min,3)/3.0-pow(d_p_min,6)/6-pow(d_p_min,4)/4))
	* (sin(Pi*x)*sin(Pi*y)+1)
	* sin(Pi*t/2);

  // with d_p - dependent definition of G(c_C,d_p)
     coeff[4] = lambda_chem - lambda_nuc*B_C_c 
	- r_g*param[5] + rhs;//r_chem - r_nuc - r_g    
  }
}

/****************************************************************/
/* finds the nodes which are Neumann and should be Dirichlet    */
/* for FEM_FCT schemes                                          */
/****************************************************************/

void CheckWrongNeumannNodes_c_C(TCollection *Coll, TFESpace2D *fespace,
int &N_neum_to_diri, int* &neum_to_diri,
int* &neum_to_diri_bdry,
double* &neum_to_diri_param)
{
  const int max_entries = 4096;
  int i, j, N_, min_val;
  int N_Cells, N_V, diri_counter = 0, found, diri_counter_1 = 0;
  int *global_numbers, *begin_index, *dof;
  int boundary_vertices[4], tmp_diri[max_entries], tmp_bdry[max_entries];
  double x[4], y[4], eps = 1e-8, tmp_param[max_entries];
  TBaseCell *cell;
  TVertex *vertex;
  FE2D CurrentElement;
  int range_3 = (int)TDatabase::ParamDB->P7;
  int range_1 = (int)TDatabase::ParamDB->P8;

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
      // Dirichlet boundary conditions on the whole boundary
      if ((fabs(x[j])<eps) || (fabs(y[j])<eps) ||(fabs(1-x[j])<eps) ||(fabs(1-y[j])<eps)) 
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
            if (fabs(x[j])<eps)
            {
              tmp_bdry[diri_counter] = 3;
              tmp_param[diri_counter] = 1-y[j];
            }
            if (fabs(1.0-x[j])<eps)
            {
              tmp_bdry[diri_counter] = 1;
              tmp_param[diri_counter] = y[j];
            }
            if (fabs(y[j])<eps)
            {
              tmp_bdry[diri_counter] = 0;
              tmp_param[diri_counter] = x[j];
            }
            if (fabs(1.0-y[j])<eps)
            {
              tmp_bdry[diri_counter] = 2;
              tmp_param[diri_counter] = 1-x[j];
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

  OutPut("CheckNeumannNodesForVelocity c_C: N_neum_to_diri " << diri_counter_1 << endl);
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


/******************************************************************************/
// MOM
/******************************************************************************/

void BoundCondition_mom(int BdComp, double Param, BoundCond &cond)
{
    int range;
    double y;
    
    cond = NEUMANN;
    
    switch(BdComp)
    {
	case 1:
	    range = (int)TDatabase::ParamDB->P8;
	    if ((Param>range/32.0)&&(Param<(range+1)/32.0))
	    {
		cond = DIRICHLET;
	    }
	    break;
	case 3: 
	    range = (int)TDatabase::ParamDB->P7;
	    y = 1-Param;
	    if ((y>range/32.0)&&(y<(range+1)/32.0))
	    {
		cond = DIRICHLET;
	    }
	    break;
    }
}

// value of boundary condition
void BoundValue_mom(int BdComp, double Param, double &value)
{
   value = 0;
}

// initial conditon
void InitialCondition_mom(double x, double y, double *values)
{
  values[0] = 0;	
}
void BilinearCoeffs_mom(int n_points, double *X, double *Y,
		    double **parameters, double **coeffs)
{
  int i, k = TDatabase::ParamDB->INTERNAL_MOMENT;
  double *coeff, *param;
  double x, y;
  double t = TDatabase::TimeDB->CURRENTTIME;
  double l_infty = TDatabase::ParamDB->BULK_l_infty;
  double u_infty = TDatabase::ParamDB->BULK_u_infty;
  double c_C_infty = TDatabase::ParamDB->BULK_c_infty;
  double c_C_infty_sat = TDatabase::ParamDB->BULK_c_C_infty_sat;
  double f_infty = TDatabase::ParamDB->BULK_f_infty;
  double d_p_min = TDatabase::ParamDB->BULK_D_P_MIN;
  double d_p_max = TDatabase::ParamDB->BULK_D_P_MAX;
  double k_g = TDatabase::ParamDB->BULK_k_g;
  double k_nuc = TDatabase::ParamDB->BULK_k_nuc;
  double eps, factor_G, factor_G0, G_c_C, c_d_p, B_c_C, f_d_p_min;

  // this is just for testing
  eps = 1e-10;
  // computed model constants
  factor_G = k_g*c_C_infty*l_infty/(u_infty*d_p_max);
  factor_G0 = k_g*c_C_infty*f_infty;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    param = parameters[i];

    if (TDatabase::ParamDB->BULK_GROWTH_RATE==2)
    {
	OutPut("MOM not for BULK_GROWTH_RATE==2 implemented !!!" << endl);
	exit(4711);
    }
    else
    {
       G_c_C = factor_G0*(param[2] - c_C_infty_sat/c_C_infty);
       c_d_p = factor_G*(param[2] - c_C_infty_sat/c_C_infty);
    }
    
    if (G_c_C > 1e-10)
    {
        // compute rate of nucleation
        B_c_C = k_nuc*pow(c_C_infty*(param[2] - 1),5);
        // truncate negative values
        if (B_c_C < 0)
	    B_c_C = 0;
	f_d_p_min = B_c_C / G_c_C;
    }
    else
	f_d_p_min = 0;

    coeff[0] = eps;
    coeff[1] = param[0]; // u1
    coeff[2] = param[1]; // u2
    coeff[3] = 0;
    if (k>0)
	coeff[4] = c_d_p * f_d_p_min * pow(d_p_min,k) +  c_d_p * k * param[3];
    else
	coeff[4] = c_d_p * f_d_p_min;	
  }
}

void ErrorPSD(int N_x, int N_y, int N_z,
	      double *x_coord, double *y_coord, double *z_coord,
	      double *sol)
{
    int i,N3,N2;
    double error = 0, val;
    double t= TDatabase::TimeDB->CURRENTTIME;
    
    N2 = (N_x+1)*(N_y+1);
    N3 = (N_x+1)*(N_y+1)*(N_z+1);

    for (i=0;i<N3;i++)
    {
	val = 1e3
	    *(2.0 - z_coord[i]*z_coord[i]*z_coord[i]-z_coord[i])
	    *(sin(Pi*x_coord[i])*sin(Pi*y_coord[i])+1)
	    * sin(Pi*t/2);
	//val = sin(t)*(sin(2*Pi*x_coord[i])*sin(2*Pi*y_coord[i])*sin(2*Pi*z_coord[i])+1);
	//if ((i>=N2+100)&&(i<N2+120))
	//OutPut(x_coord[i] << " " << y_coord[i] << " " << z_coord[i] << 
	//      " :ex " << val << " :sol " << sol[i] << endl);
	val -= sol[i];
	error += val*val;
    }
    OutPut(t << " psderr " << sqrt(error/N3) << endl);
    //exit(1);
}


//#endif
