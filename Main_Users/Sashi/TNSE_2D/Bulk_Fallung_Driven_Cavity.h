//
// u(x,y) = unknown
// p(x,y) = unknown

void ExampleFile()
{
  int range;

  OutPut(" Example: Bulk_Fallung_Driven_Cavity.h ") ;
  OutPut(" inflow (u_infty)" << TDatabase::ParamDB->BULK_u_infty);
  OutPut(" upper lid " << TDatabase::ParamDB->P5);
  OutPut(" left " << (int)TDatabase::ParamDB->P7);
  OutPut(" right " << (int)TDatabase::ParamDB->P8);
  OutPut(" lower " << (int)TDatabase::ParamDB->P9 << endl);

  range = (int)TDatabase::ParamDB->P7;
  if ((range<1)||(range>30))
  {
      OutPut("left boundary out of range !!!"<< endl);
      exit(4711);
  }

  range = (int)TDatabase::ParamDB->P8;
  if ((range<1)||(range>30))
  {
      OutPut("right boundary out of range !!!"<< endl);
      exit(4711);
  }

  range = (int)TDatabase::ParamDB->P9;
  if ((range<1)||(range>30))
  {
      OutPut("lower boundary out of range !!!"<< endl);
      exit(4711);
  }
  // set some parameters
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
  int lower = (int)TDatabase::ParamDB->P9;

  if ((i==0)&&((t>lower/32.0)&&(t<(lower+2)/32.0)))
      cond = NEUMANN;
  else
      cond = DIRICHLET;
  // cond = DIRICHLET;
  TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 0; 
}

void U1BoundValue(int BdComp, double Param, double &value)
{
    int range;
    double y, fac = 1024; // 32^2
    

  switch(BdComp)
  {
     case 0: 
        value = 0;
        break;
     case 1:
	range = (int)TDatabase::ParamDB->P8;
	if ((Param>range/32.0)&&(Param<(range+1)/32.0))
        {
           y = Param;
           value =  6*(y-range/32.0)*(y-(range+1)/32.0)*fac;
	}
	else
	    value = 0;
	break;
	// upper boundary
    case 2: if(Param<0.00001 || Param>0.99999) 
              value = 0;
            else
               value = TDatabase::ParamDB->P5/TDatabase::ParamDB->BULK_u_infty;
               //value = 0;
            break;
    case 3: 
	range = (int)TDatabase::ParamDB->P7;
        y = 1-Param;
	if ((y>range/32.0)&&(y<(range+1)/32.0))
        {
           value = -6*(y-range/32.0)*(y-(range+1)/32.0)*fac;
        }
	else
	    value = 0;
	break;
    default: cout << "wrong boundary part number" << endl;
  }
}

void U2BoundValue(int BdComp, double Param, double &value)
{
  value = 0;
/*  int lower = (int)TDatabase::ParamDB->P9;

  if ((BdComp==0)&&((Param>lower/32.0)&&(Param<(lower+1)/32.0)))
      value = -2;
  else
      value = 0;
*/
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
   int range;
   double eps=1e-8;
   double t = TDatabase::TimeDB->CURRENTTIME;
   double t0 = TDatabase::TimeDB->T0;

   if (t<t0)
     values[0] = 0;
   else
   {		
     range = (int)TDatabase::ParamDB->P7;
     if ((fabs(x)<1e-7)&&(y>=range/32.0-eps)&&(y<=(range+1)/32.0+eps))
       values[0] = 1;
     else
       values[0] = 0;
   }
}

// kind of boundary condition (for FE space needed)
void BoundCondition_c_A(int BdComp, double Param, BoundCond &cond)
{
   int range;
   double y,eps=1e-8;

   if ((TDatabase::ParamDB->BULK_REACTION_DISC == GALERKIN) && 
      (TDatabase::ParamDB->BULK_SOLD_PARAMETER_TYPE == FEM_FCT_LIN))
    cond = NEUMANN;
   else
   {      
     if (BdComp==3)
     {
      range = (int)TDatabase::ParamDB->P7;
      y = 1-Param;
      if ((y>=range/32.0-eps)&&(y<=(range+1)/32.0+eps))
      {
         cond = DIRICHLET;
      }
      else
	  cond = NEUMANN;
   }
   else
     cond = NEUMANN;
  }
}

// value of boundary condition
void BoundValue_c_A(int BdComp, double Param, double &value)
{
   int range;
   double y,eps=1e-8;
   double t = TDatabase::TimeDB->CURRENTTIME;
   double t0 = TDatabase::TimeDB->T0;
  
   if (t<t0)
     value = 0;
   else
   {
   if (BdComp==3)
   {
      range = (int)TDatabase::ParamDB->P7;
      y = 1-Param;
      if ((y>=range/32.0-eps)&&(y<=(range+1)/32.0+eps))
      {
	  value = 1;
      }
      else
	  value = 0;
   }
   else
     value = 0;
   }
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
    coeff[1] = param[1]; // u1
    coeff[2] = param[2]; // u2
    coeff[3] = c * param[0];
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
   int range;
   double eps = 1e-8;
   double t = TDatabase::TimeDB->CURRENTTIME;
   double t0 = TDatabase::TimeDB->T0;

   if (t<t0)
     values[0] = 0;
    else
    {
    range = (int)TDatabase::ParamDB->P8;
    if ((fabs(1-x) <1e-7)&&(y>=range/32.0-eps)&&(y<=(range+1)/32.0+eps))
	values[0] = 1;
    else
	values[0] = 0;
}	
}

// kind of boundary condition (for FE space needed)
void BoundCondition_c_B(int BdComp, double Param, BoundCond &cond)
{
   int range;
   double eps = 1e-8;

  if ((TDatabase::ParamDB->BULK_REACTION_DISC == GALERKIN) && 
      (TDatabase::ParamDB->BULK_SOLD_PARAMETER_TYPE == FEM_FCT_LIN))
    cond = NEUMANN;
   else
   {
     if (BdComp==1)
     {
       range = (int)TDatabase::ParamDB->P8;
       if ((Param>=range/32.0-eps)&&(Param<=(range+1)/32.0+eps))
       {
         cond = DIRICHLET;
       }
       else
         cond = NEUMANN;
     }
     else
       cond = NEUMANN;
   }
}

// value of boundary condition
void BoundValue_c_B(int BdComp, double Param, double &value)
{
   int range;
   double t = TDatabase::TimeDB->CURRENTTIME;
   double t0 = TDatabase::TimeDB->T0;
  
   if (t<t0)
     value = 0;
   else
   {   
   if (BdComp==1)
   {
      range = (int)TDatabase::ParamDB->P8;
      if ((Param>=range/32.0)&&(Param<=(range+1)/32.0))
      {
	  value = 1;
      }
      else
	  value = 0;
   }
   else
     value = 0;
   }
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

// initial conditon
void InitialCondition_c_C(double x, double y, double *values)
{
  values[0] = 0;	
}

// kind of boundary condition (for FE space needed)
void BoundCondition_c_C(int BdComp, double Param, BoundCond &cond)
{
    int range;
    double y;

   cond = NEUMANN;

/*
   if ((TDatabase::ParamDB->BULK_REACTION_DISC == GALERKIN) && 
      (TDatabase::ParamDB->BULK_SOLD_PARAMETER_TYPE == FEM_FCT_LIN))
     cond = NEUMANN;
   else
   {
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
*/
}

// value of boundary condition
void BoundValue_c_C(int BdComp, double Param, double &value)
{
   value = 0;
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
  double x, y;
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
  lambda_nuc = C_nuc*k_nuc*pow(d_p_0,3)*L_infty*pow(c_C_infty,4)/U_infty;
  //OutPut("lambda_chem " << lambda_chem << "  lambda_nuc " << lambda_nuc << " " << endl);
  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    param = parameters[i];

    if (param[4] < 1)
      B_C_c = 0;
    else
      B_C_c = pow(param[4] - 1,5);

    coeff[0] = eps;
    coeff[1] = param[0]; // u1
    coeff[2] = param[1]; // u2
    coeff[3] = 0;
    // following script from 05/08/23 for r_g
    r_g = param[4] - c_C_infty_sat/c_C_infty;
    // with d_p -- dependent definition of G(c_C,d_p)
    if (TDatabase::ParamDB->BULK_GROWTH_RATE == 2)
       r_g = 1;
    // with d_p - dependent definition of G(c_C,d_p)
     coeff[4] = lambda_chem*param[2]*param[3] - lambda_nuc*B_C_c 
	- r_g*param[5];//r_chem - r_nuc - r_g    
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
      // comment next two if cases
      /*// vertex on the upper lid
        if ((fabs(x[j])<eps) && ((y[j]>=range_3/32.0-eps)&&(y[j]<=(range_3+1)/32.0+eps)))
      {
        boundary_vertices[j] = 1;
        found++;
      }
      if ((fabs(1.0-x[j])<eps) && ((y[j]>=range_1/32.0-eps)&&(y[j]<=(range_1+1)/32.0+eps)))
      {
        boundary_vertices[j] = 1;
        found++;
	}*/
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

//#endif










