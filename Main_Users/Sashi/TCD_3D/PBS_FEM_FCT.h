#define __UREA__
#define __SIMPATURS__
  
#include <Urea_3d4d.h>
#include <MacroCell.h>

void ExampleFile()
{
  // for velocity switch
  TDatabase::ParamDB->INTERNAL_CONVECTION_EQ_VELOFIELD=0;

  TDatabase::ParamDB->N_CELL_LAYERS = 3;
  TDatabase::ParamDB->DRIFT_Z = 1;
  TDatabase::ParamDB->INTERNAL_PROBLEM_IDENTITY = 1356;

  TDatabase::ParamDB->UREA_SOLD_PARAMETER_TYPE = FEM_FCT_LIN;
  TDatabase::ParamDB->UREA_temp_infty = TDatabase::ParamDB->REACTOR_P23;
  TDatabase::ParamDB->UREA_c_infty = TDatabase::ParamDB->REACTOR_P25;
  TDatabase::ParamDB->UREA_f_infty = TDatabase::ParamDB->REACTOR_P20;
  TDatabase::ParamDB->UREA_rho_sat_1 = 35.36400;
  TDatabase::ParamDB->UREA_rho_sat_2= 1.305000; 
  
#define __FEMFCT__ 
  
#ifdef _MPI
 MPI_Comm Comm;
 int rank;

 Comm = TDatabase::ParamDB->Comm;
 MPI_Comm_rank(Comm, &rank);

 if(rank==TDatabase::ParamDB->Par_P0)
#endif
 {  
  OutPut("Example: PBS_FEM_FCT.h " << endl);
    
  OutPut("UREA_REACTION_DISC: " << TDatabase::ParamDB->UREA_REACTION_DISC << endl);
  OutPut("UREA_PB_DISC: " << TDatabase::ParamDB->UREA_PB_DISC << endl);
  OutPut("UREA_PB_DISC_STAB: " << TDatabase::ParamDB->UREA_PB_DISC_STAB<<endl);
  OutPut("UREA_SOLD_PARAMETER_TYPE: "<< TDatabase::ParamDB->UREA_SOLD_PARAMETER_TYPE <<endl);
  OutPut("UREA_MODEL: " << TDatabase::ParamDB->UREA_MODEL << endl);
  OutPut("UREA_CONC_TOL: " << TDatabase::ParamDB->UREA_CONC_TOL << endl);
  OutPut("UREA_CONC_MAXIT: " << TDatabase::ParamDB->UREA_CONC_MAXIT << endl);

  OutPut("UREA_l_infty: " << TDatabase::ParamDB->UREA_l_infty <<endl);
  OutPut("UREA_u_infty: " << TDatabase::ParamDB->UREA_u_infty <<endl);
  OutPut("UREA_c_infty: " << TDatabase::ParamDB->UREA_c_infty <<endl);
  OutPut("UREA_temp_infty: " << TDatabase::ParamDB->UREA_temp_infty <<endl);
  OutPut("UREA_f_infty: " << TDatabase::ParamDB->UREA_f_infty<<endl);
  OutPut("UREA_nu: " << TDatabase::ParamDB->UREA_nu<<endl); 
  OutPut("UREA_rho: " << TDatabase::ParamDB->UREA_rho<<endl);
  OutPut("UREA_c_p: " << TDatabase::ParamDB->UREA_c_p<<endl);
  OutPut("UREA_lambda: " << TDatabase::ParamDB->UREA_lambda<<endl); 
  OutPut("UREA_D_P_0: " << TDatabase::ParamDB->UREA_D_P_0<<endl); 
  OutPut("UREA_D_P_MAX: " << TDatabase::ParamDB->UREA_D_P_MAX <<endl);
  OutPut("UREA_k_v: " << TDatabase::ParamDB->UREA_k_v<<endl); 
  OutPut("UREA_m_mol: " << TDatabase::ParamDB->UREA_m_mol<<endl);
  OutPut("UREA_D_J: " << TDatabase::ParamDB->UREA_D_J<<endl);
  OutPut("UREA_rho_d: " << TDatabase::ParamDB->UREA_rho_d <<endl);
  OutPut("UREA_k_g: " << TDatabase::ParamDB->UREA_k_g<<endl);
  OutPut("UREA_g: " << TDatabase::ParamDB->UREA_g <<endl);
  OutPut("UREA_rho_sat_1: " << TDatabase::ParamDB->UREA_rho_sat_1 <<endl);
  OutPut("UREA_rho_sat_2: " << TDatabase::ParamDB->UREA_rho_sat_2<<endl); 
  OutPut("UREA_beta_nuc: " << TDatabase::ParamDB->UREA_beta_nuc<<endl); 
  OutPut("UREA_alfa_nuc: " << TDatabase::ParamDB->UREA_alfa_nuc<<endl); 
  OutPut("UREA_INFLOW_SCALE: " << TDatabase::ParamDB->UREA_INFLOW_SCALE <<endl);
  OutPut("UREA_inflow_time: " << TDatabase::ParamDB->UREA_inflow_time <<endl);
 }
  // set some parameters
  //TDatabase::ParamDB->GRID_TYPE = 3;
  //OutPut("GRID_TYPE set to " << TDatabase::ParamDB->GRID_TYPE << endl);
}

void Exact(double x, double y, double z, double *values)
{ 
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}

void Exact_Psd_Intl(double x, double y, double z, double l, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}

// kind of boundary condition (for FE space needed)
void BoundCondition_NSE(double x, double y, double z, BoundCond &cond)
{
  double eps = 1e-8;

  cond = DIRICHLET;

  if (fabs(x-210)<eps)
  {
       // outflow 
       cond = NEUMANN;
       //OutPut("neum");
       TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 0;
  }
}

// ========================================================================
// definitions for the temperature
// ========================================================================

// initial conditon
void InitialCondition_temp(double x, double y, double z, double *values)
{    
  values[0] = TDatabase::ParamDB->REACTOR_P23/TDatabase::ParamDB->UREA_temp_infty;
}

// kind of boundary condition (for FE space needed)
void BoundCondition_temp(double x, double y, double z, BoundCond &cond)
{
    // this will be set correctly in CheckWrongNeumannNodes_temp 
    cond = NEUMANN;
}

// value of boundary condition
//void BoundValue_c_A(int BdComp, double Param, double &value)
void BoundValue_temp(double x, double y, double z, double &value)
{
  double eps = 1e-8;
  double T_D = TDatabase::ParamDB->REACTOR_P23;
  double T_W = TDatabase::ParamDB->REACTOR_P24;
  
  value = 0;

    if ((fabs(x)<eps)&& (y>=0.33333333) && (y<=0.66666667) &&
      (z>=0.33333333) && (z<=0.66666667))
    {
      value = T_D/TDatabase::ParamDB->UREA_temp_infty;      
      //printf(" inflow temp \n" );
    } 
    else
    {
      // outlet
      if ((fabs(x-210)<eps)&&(fabs(y)>eps)&&(fabs(z)>eps)&&(fabs(y-1.0)>eps)&&(fabs(z-1.0)>eps))
      {  value = 0; }
      else // wall
      {  value = T_W/TDatabase::ParamDB->UREA_temp_infty; }
    }

}

// ========================================================================
// BilinearCoeffs for Heat 
// ========================================================================
void BilinearCoeffs_Heat(int n_points, double *x, double *y, double *z,
               double **parameters, double **coeffs)
{
  int i;
  double eps, *coeff, *param;
 
  if(TDatabase::ParamDB->REACTOR_P0)
    eps = 1.0/TDatabase::ParamDB->REACTOR_P0;
  else
    eps = 0.;
  
  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    param = parameters[i];    

    coeff[0] = eps;
    if(TDatabase::ParamDB->INTERNAL_CONVECTION_EQ_VELOFIELD)
     {
      coeff[1] = param[0];  // u1
      coeff[2] = param[1];  // u2
      coeff[3] = param[2];  // u3      
//       cout<< "coeff[0] eps " << eps << " u1 " << param[0]  << "u2 " << param[1] << endl;
     }
    else
     {
      coeff[1] = 0;  // u1
      coeff[2] = 0;  // u2
      coeff[3] = 0;  // u2      
     }     
    coeff[4] = 0;
    coeff[5] = 0.;
    coeff[6] = 0;  
  }
}

/****************************************************************/
/* finds the nodes which are Neumann and should be Dirichlet    */
/* for FEM_FCT schemes                                          */
/****************************************************************/

void CheckWrongNeumannNodes_temp(TCollection *Coll, TFESpace3D *fespace,
				int &N_neum_to_diri, int* &neum_to_diri,
				double* &neum_to_diri_x, 
				double* &neum_to_diri_y,
				double* &neum_to_diri_z) 
{
  const int max_entries = 44000;
  int i, j, N_, min_val;
  int N_Cells, N_V, diri_counter = 0, found, diri_counter_1 = 0;
  int *global_numbers, *begin_index, *dof;
  int boundary_vertices[8], tmp_diri[max_entries]; 
  double x[8], y[8], z[8], eps = 1e-8, tmp_x[max_entries], tmp_y[max_entries], tmp_z[max_entries];
  TBaseCell *cell;
  TVertex *vertex;
  FE3D CurrentElement;

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
      vertex->GetCoords(x[j], y[j], z[j]);
      // vertex on the boundary 
      if ((fabs(x[j])<eps) || (fabs(y[j])<eps) || (fabs(z[j])<eps)  || 
	   (fabs(y[j]-1)<eps) || (fabs(z[j]-1)<eps))
      {
        boundary_vertices[j] = 1;
        found++;
      }
    }
    // no cell with face with vertex on the boundary
    if (found<3)
      continue;
    // finite element on the mesh cell
    CurrentElement = fespace->GetFE3D(i, cell);
    // number of basis functions (= number of d.o.f.)
    N_ = TFEDatabase3D::GetN_BaseFunctFromFE3D(CurrentElement);
    // the array which gives the mapping of the local to the global d.o.f.
    dof = global_numbers+begin_index[i];
    switch(CurrentElement)
    {
      // P_1, Q_1
      case C_P1_3D_T_A:
      case C_Q1_3D_H_A:
      case C_Q1_3D_H_M:
        for (j=0;j<N_V;j++)
        {
          // vertex on the boundary
          if (boundary_vertices[j])
          {
	      // tetrahedron
	      if (CurrentElement==C_P1_3D_T_A)
		  tmp_diri[diri_counter] = dof[j];
	      else
	      {
		  switch(j)
		  {
		      case 0:
		      case 1:
		      case 4:
		      case 5:
			  tmp_diri[diri_counter] = dof[j];
			  break;
		      case 2:
			  tmp_diri[diri_counter] = dof[3];
			  break;
		      case 3:
			  tmp_diri[diri_counter] = dof[2];
			  break;
		      case 6:
			  tmp_diri[diri_counter] = dof[7];
			  break;
		      case 7:
			  tmp_diri[diri_counter] = dof[6];
			  break;
		  }
	      }
	      if (diri_counter > max_entries)
	      {
		  OutPut("tmp_diri too short !!!"<<endl);
		  exit(4711);
	      }
	   
	      if ((fabs(x[j])<eps) || (fabs(y[j])<eps) || (fabs(z[j])<eps)  || 
		  (fabs(y[j]-1)<eps) || (fabs(z[j]-1)<eps))
	      {
		  tmp_x[diri_counter] = x[j];
		  tmp_y[diri_counter] = y[j];
		  tmp_z[diri_counter] = z[j];
	      }
	     // OutPut( tmp_diri[diri_counter] << " " <<
	     // 	      tmp_x[diri_counter] << " " << tmp_y[diri_counter] 
	      //	      << " " << tmp_z[diri_counter]  << endl);
	      diri_counter++;
          }
        }
	//OutPut(endl);
        break;
	default:
	    OutPut("CheckWrongNeumannNodes_temp not implemented for element "
		   << CurrentElement << endl);
	    OutPut("code can be run without CheckWrongNeumannNodes_temp, just delete the exit" << endl);
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
  
  OutPut("CheckWrongNeumannNodes_temp: N_neum_to_diri " << diri_counter_1 << endl);
  N_neum_to_diri = diri_counter_1;
  // allocate array for the indices
  neum_to_diri = new int[diri_counter_1];
  // allocate array for the corresponding x coordinate
  neum_to_diri_x = new double[diri_counter_1];
  // allocate array for the corresponding y coordinate
  neum_to_diri_y = new double[diri_counter_1];
  // allocate array for the corresponding z coordinate
  neum_to_diri_z = new double[diri_counter_1];

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
    neum_to_diri_x[i] = tmp_x[found];
    neum_to_diri_y[i] = tmp_y[found];
    neum_to_diri_z[i] = tmp_z[found];
    tmp_diri[found] = -1;
  }

 // for (i=0;i<diri_counter_1;i++)
 // {
  //  OutPut(i << " " << neum_to_diri[i] << " " << neum_to_diri_x[i]  <<
  //    " " << neum_to_diri_y[i]  <<  " " << neum_to_diri_z[i]  << endl);
 // }
}

// ========================================================================
// definitions for the concentration
// ========================================================================

// initial condition
void InitialCondition_conc(double x, double y, double z, double *values)
{
  double val;
  
  val  = TDatabase::ParamDB->REACTOR_P23;
  val = TDatabase::ParamDB->UREA_rho_sat_1 + TDatabase::ParamDB->UREA_rho_sat_2*(val - 273.15); 

//        printf("InitialCondition_conc %f,\n", val);    
  values[0] = val/TDatabase::ParamDB->UREA_c_infty;
}

// kind of boundary condition (for FE space needed)
void BoundCondition_conc(double x, double y, double z, BoundCond &cond)
{
  // this will be set correctly in CheckWrongNeumannNodes_conc
   cond = NEUMANN;
}

// value of boundary condition
void BoundValue_conc(double x, double y, double z, double &value)
{
  double eps = 1e-8, val, temp;

  value = 0;
  
  // inflow
  if (TDatabase::TimeDB->CURRENTTIME >= TDatabase::TimeDB->T1)
  {
    if ((fabs(x)<eps)&& (y>=0.33333333) && (y<=0.66666667) &&
      (z>=0.33333333) && (z<=0.66666667))
      {
       BoundValue_temp(x,y,z,temp);
       temp *=TDatabase::ParamDB->UREA_temp_infty;
       val = TDatabase::ParamDB->UREA_rho_sat_1 + TDatabase::ParamDB->UREA_rho_sat_2*(temp - 273.15); 

       value = val/TDatabase::ParamDB->UREA_c_infty;       
//        printf("BoundValue_conc %f,\n", value);
      }
  }
}

// ========================================================================
// BilinearCoeffs for concentration 
// ========================================================================
void BilinearCoeffs_Conc(int n_points, double *x, double *y, double *z,
               double **parameters, double **coeffs)
{
  int i;
  double eps, *coeff, *param;

  if(TDatabase::ParamDB->REACTOR_P1)
    eps = 1.0/TDatabase::ParamDB->REACTOR_P1;
  else
    eps = 0.;
  
  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    param = parameters[i];    

    coeff[0] = eps;
    if(TDatabase::ParamDB->INTERNAL_CONVECTION_EQ_VELOFIELD)
     {
      coeff[1] = param[0];  // u1
      coeff[2] = param[1];  // u2
      coeff[3] = param[2];  // u3      
//       cout<< "coeff[0] eps " << eps << " u1 " << param[0]  << "u2 " << param[1] << endl;
     }
    else
     {
      coeff[1] = 0;  // u1
      coeff[2] = 0;  // u2
      coeff[3] = 0;  // u2      
     }    
    coeff[4] = 0.;
    coeff[5] = 0.;
    coeff[6] = 0.;  
  }
}


void CheckWrongNeumannNodes_conc(TCollection *Coll, TFESpace3D *fespace,
				int &N_neum_to_diri, int* &neum_to_diri,
				double* &neum_to_diri_x, 
				double* &neum_to_diri_y,
				double* &neum_to_diri_z) 
{
  const int max_entries = 16000;
  int i, j, N_, min_val;
  int N_Cells, N_V, diri_counter = 0, found, diri_counter_1 = 0;
  int *global_numbers, *begin_index, *dof;
  int boundary_vertices[8], tmp_diri[max_entries]; 
  double x[8], y[8], z[8], eps = 1e-8, tmp_x[max_entries], tmp_y[max_entries], tmp_z[max_entries];
  TBaseCell *cell;
  TVertex *vertex;
  FE3D CurrentElement;

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
      vertex->GetCoords(x[j], y[j], z[j]);
      // vertex on the boundary 
      if ((fabs(x[j])<eps)&& (y[j]>=0.33333333) && (y[j]<=0.66666667) &&
	  (z[j]>=0.33333333) && (z[j]<=0.66666667))
      {
        boundary_vertices[j] = 1;
        found++;
      }
    }
    // no cell with face with vertex on the boundary
    if (found<3)
      continue;
    // finite element on the mesh cell
    CurrentElement = fespace->GetFE3D(i, cell);
    // number of basis functions (= number of d.o.f.)
    N_ = TFEDatabase3D::GetN_BaseFunctFromFE3D(CurrentElement);
    // the array which gives the mapping of the local to the global d.o.f.
    dof = global_numbers+begin_index[i];
    switch(CurrentElement)
    {
      // P_1, Q_1
      case C_P1_3D_T_A:
      case C_Q1_3D_H_A:
      case C_Q1_3D_H_M:
        for (j=0;j<N_V;j++)
        {
          // vertex on the boundary
          if (boundary_vertices[j])
          {
	      // tetrahedron
	      if (CurrentElement==C_P1_3D_T_A)
		  tmp_diri[diri_counter] = dof[j];
	      else
	      {
		  switch(j)
		  {
		      case 0:
		      case 1:
		      case 4:
		      case 5:
			  tmp_diri[diri_counter] = dof[j];
			  break;
		      case 2:
			  tmp_diri[diri_counter] = dof[3];
			  break;
		      case 3:
			  tmp_diri[diri_counter] = dof[2];
			  break;
		      case 6:
			  tmp_diri[diri_counter] = dof[7];
			  break;
		      case 7:
			  tmp_diri[diri_counter] = dof[6];
			  break;
		  }
	      }
	      if (diri_counter > max_entries)
	      {
		  OutPut("tmp_diri too short !!!"<<endl);
		  exit(4711);
	      }
	      if ((fabs(x[j])<eps)&& (y[j]>=0.33333333) && (y[j]<=0.66666667) &&
		  (z[j]>=0.33333333) && (z[j]<=0.66666667))
	      {
		  tmp_x[diri_counter] = x[j];
		  tmp_y[diri_counter] = y[j];
		  tmp_z[diri_counter] = z[j];
	      }
	     // OutPut( tmp_diri[diri_counter] << " " <<
	     //	      tmp_x[diri_counter] << " " << tmp_y[diri_counter] 
	      //	      << " " << tmp_z[diri_counter]  << endl);
	      diri_counter++;
          }
        }
	//OutPut(endl);
        break;
	default:
	    OutPut("CheckWrongNeumannNodes_temp not implemented for element "
		   << CurrentElement << endl);
	    OutPut("code can be run without CheckWrongNeumannNodes_temp, just delete the exit" << endl);
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
  
  OutPut("CheckWrongNeumannNodes_conc: N_neum_to_diri " << diri_counter_1 << endl);
  N_neum_to_diri = diri_counter_1;
  // allocate array for the indices
  neum_to_diri = new int[diri_counter_1];
  // allocate array for the corresponding x coordinate
  neum_to_diri_x = new double[diri_counter_1];
  // allocate array for the corresponding y coordinate
  neum_to_diri_y = new double[diri_counter_1];
  // allocate array for the corresponding z coordinate
  neum_to_diri_z = new double[diri_counter_1];

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
    neum_to_diri_x[i] = tmp_x[found];
    neum_to_diri_y[i] = tmp_y[found];
    neum_to_diri_z[i] = tmp_z[found];
    tmp_diri[found] = -1;
  }

 // for (i=0;i<diri_counter_1;i++)
 // {
 //   OutPut(i << " " << neum_to_diri[i] << " " << neum_to_diri_x[i]  <<
 //     " " << neum_to_diri_y[i]  <<  " " << neum_to_diri_z[i]  << endl);
 // }
}


// ========================================================================
// definitions for the PSD
// ========================================================================

// initial condition
void InitialCondition_psd(double x, double y, double z, double *values)
{
  values[0] = 0;
}

void InitialCondition_psd_Intl(int N_Coord, double *X, double *values)
{
  values[0] = 0;
}

// kind of boundary condition (for FE space needed)
void BoundCondition_psd(double x, double y, double z, BoundCond &cond)
{
  double eps = 1e-8; 
  
  // FEFCT not yet implemented for PBE
   cond = NEUMANN;   
  
    if ((fabs(x)<eps)&& (y>=0.33333333) && (y<=0.66666667) &&
      (z>=0.33333333) && (z<=0.66666667))
     { cond = DIRICHLET; }
   
}

// value of boundary condition
void BoundValue_psd(double x, double y, double z, double &value)
{
 int i, found=0;
 double m, t, eps = 1e-8, L_max = 1.691440e-03, a = TDatabase::ParamDB->REACTOR_P29;
 double scale_f_max = 1./TDatabase::ParamDB->REACTOR_P20;
 
  // inflow
   if ((fabs(x)<eps)&& (y>=0.33333333) && (y<=0.66666667) &&
      (z>=0.33333333) && (z<=0.66666667))  
    {
      
  double inlet_coord[100] ={2.500000e-06,1.356050e-05,3.068150e-05,4.780300e-05,6.492450e-05,8.204550e-05,9.916650e-05,1.162875e-04,1.334090e-04,1.505305e-04,1.676515e-04,1.847725e-04,2.018940e-04,2.190155e-04,2.361365e-04,2.532575e-04,2.703785e-04,2.875000e-04,3.046215e-04,3.217425e-04,3.388635e-04,3.559845e-04,3.731060e-04,3.902275e-04,4.073485e-04,4.244695e-04,4.415910e-04,4.587125e-04,4.758335e-04,4.929545e-04,5.100755e-04,5.271970e-04,5.443185e-04,5.614395e-04,5.785605e-04,5.956815e-04,6.128030e-04,6.299245e-04,6.470455e-04,6.641665e-04,6.812875e-04,6.984090e-04,7.155305e-04,7.326515e-04,7.497725e-04,7.668940e-04,7.840155e-04,8.011365e-04,8.182575e-04,8.353785e-04,8.525000e-04,8.696215e-04,8.867425e-04,9.038635e-04,9.209845e-04,9.381060e-04,9.552275e-04,9.723485e-04,9.894695e-04,1.006591e-03,1.023712e-03,1.040833e-03,1.057955e-03,1.075075e-03,1.092197e-03,1.109318e-03,1.126439e-03,1.143561e-03,1.160682e-03,1.177803e-03,1.194925e-03,1.212045e-03,1.229167e-03,1.246287e-03,1.263409e-03,1.280530e-03,1.297651e-03,1.314773e-03,1.331894e-03,1.349016e-03,1.366137e-03,1.383257e-03,1.400379e-03,1.417500e-03,1.434622e-03,1.451742e-03,1.468863e-03,1.485985e-03,1.503106e-03,1.520228e-03,1.537348e-03,1.554470e-03,1.571591e-03,1.588713e-03,1.605834e-03,1.622954e-03,1.640075e-03,1.657197e-03,1.674318e-03,1.691440e-03};

  double inlet_f_L_seed[100] = {0.000000e+00,1.539144e+09,2.082543e+09,1.056886e+09,1.032820e+09,9.174217e+08,7.339395e+08,6.383108e+08,5.630089e+08,4.776040e+08,3.952118e+08,3.244646e+08,2.668657e+08,2.186269e+08,1.774119e+08,1.430176e+08,1.142823e+08,9.072941e+07,7.094795e+07,5.432124e+07,4.189605e+07,3.280686e+07,2.562972e+07,2.013453e+07,1.620424e+07,1.315460e+07,1.040263e+07,7.846039e+06,5.728553e+06,4.495501e+06,3.800058e+06,3.005911e+06,2.229857e+06,1.675695e+06,1.307405e+06,1.003372e+06,7.486895e+05,6.805768e+05,7.435860e+05,7.812756e+05,7.396855e+05,6.082531e+05,3.731120e+05,1.582579e+05,8.976509e+04,1.072953e+05,9.867716e+04,5.302092e+04,3.384872e+04,7.221184e+04,1.380859e+05,1.713877e+05,1.450725e+05,9.890296e+04,6.954686e+04,3.619885e+04,1.297526e+04,1.931888e+04,5.546763e+04,1.066047e+05,1.070156e+05,5.457797e+04,1.540227e+04,3.078500e+03,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00};       
      
      
      
    Dscal(100, scale_f_max, inlet_f_L_seed);     
      
    for (i=100;i>0;i--)
     { 
      if((a*L_max<=inlet_coord[i])&&(a*L_max>inlet_coord[i-1]))
      {
        m=(inlet_f_L_seed[i]-inlet_f_L_seed[i-1])/(inlet_coord[i] - inlet_coord[i-1]);
        value = m *(a*L_max- inlet_coord[i-1]);
        value += inlet_f_L_seed[i-1];
// 	        OutPut(" a ist !!!"<< a *L_max <<" i ist !!!" <<i << " i ist !!!"<<inlet_coord[i]<<endl);
        found=1;
        break;
       }
      }      
    }// if ((fabs(x)<eps)&& (y>=0.33333333) && (y<=0.666666
    
    if(found==0)
      value= 0.;    
}


void BoundCondition_LminLMax(BoundCond &cond_Lmin, BoundCond &cond_Lmax)
{
  cond_Lmin = DIRICHLET;
//   cond_Lmin = NEUMANN;
  cond_Lmax = NEUMANN;
}

void BoundValue_LMin(double x, double y, double z, double *values)
 {  
   //nucleation
   // no nucleation, so f = 0 at L=0
    values[0] = 0.;
    values[1] = 0.;
    values[2] = 0.;
 }


void BoundValue_LMax(double x, double y,  double z,  double *values)
 {
    values[0] = 0.;
    values[1] = 0.;
    values[2] = 0.; 
 }


// ========================================================================
// BilinearCoeffs for PSD 
// ========================================================================
void BilinearCoeffs_Psd(int n_points, double *x, double *y, double *z,
               double **parameters, double **coeffs)
{
  double eps;
  int i;
  double *coeff, *param;

  if(TDatabase::ParamDB->REACTOR_P2)
   { eps = 1.0/TDatabase::ParamDB->REACTOR_P2; }
  else
   { eps = 0.;}
  
  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    param = parameters[i];
    

    coeff[0] = eps;
    
    if(TDatabase::ParamDB->INTERNAL_CONVECTION_EQ_VELOFIELD)
     {
      coeff[1] = param[0];  // u1
      coeff[2] = param[1];  // u2
      coeff[3] = param[2];  // u3      
//       cout<< "coeff[0] eps " << eps << " u1 " << param[0]  << "u2 " << param[1] << endl;
     }
    else
     {
      coeff[1] = 0;  // u1
      coeff[2] = 0;  // u2
      coeff[3] = 0;  // u2      
     }
     
    coeff[4] = 0;
    coeff[5] = 0.;
    coeff[6] = 0;  
  }
}


void BilinearCoeffs_Psd_Intl(int n_points, int N_Dim, double **Coords,
                             double **parameters, double **coeffs)
{
  int i;
  double eps, *coeff, b;
  double t = TDatabase::TimeDB->CURRENTTIME;

  b = -1e8;// negative, so that C will be taken from the PBS growth term
 
  if(TDatabase::ParamDB->REACTOR_P3)
    eps = 1.0/TDatabase::ParamDB->REACTOR_P3;
  else
    eps = 0.;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];  
   
    coeff[0] = eps;  // diffusion in L -direction    
    coeff[1] = b; // convection in L -direction
    coeff[2] = 0.;    // reaction L -direction
    coeff[3] = 0.;   // rhs in L -direction
  }
}

void GetExampleFileData(BoundCondFunct3D **BoundaryConditions, BoundValueFunct3D **BoundValues, 
                        DoubleFunct3D **InitiaValues, CoeffFct3D **BilinearCoeffs, 
                        int &N_PBEqns, int &N_IndepntScalarEqns, int *Disctypes)
{
  
   N_IndepntScalarEqns = 2;
   N_PBEqns = 1;
   #define __PBS__

   BilinearCoeffs[0] = BilinearCoeffs_Heat;
   BilinearCoeffs[1] = BilinearCoeffs_Conc;
   BilinearCoeffs[2] = BilinearCoeffs_Psd;

   BoundaryConditions[0] = BoundCondition_temp;
   BoundaryConditions[1] = BoundCondition_conc;
   BoundaryConditions[2] = BoundCondition_psd;

   BoundValues[0] = BoundValue_temp;
   BoundValues[1] = BoundValue_conc;
   BoundValues[2] = BoundValue_psd;

   InitiaValues[0] = InitialCondition_temp;
   InitiaValues[1] = InitialCondition_conc;
   InitiaValues[2] = InitialCondition_psd;

   Disctypes[0] = GALERKIN;
   Disctypes[1] = GALERKIN;
//    Disctypes[2] = SDFEM;
   Disctypes[2] = GALERKIN;
}

void Generate1DMesh(TDomain *Domain, double Start, double End, int N_Cells)
{
  int i, j, N_V;
  int *Lines;
  double len, h, x, y, z, *X;
  TVertex **Vetrex;
  TJoint *Joint;
  TBaseCell  **CellTree;

  N_V = N_Cells+1;
  X = new double[N_V];

  for(i=0; i<N_V; i++)
   X[i] = 1. + (1. - Start)*(tanh(2.75*(double(i)/double(N_Cells) - 1.)))/tanh(2.75);

//   h = (End - Start)/N_Cells;
//   for(i=1; i<N_V; i++)
//    X[i] =  h*(double)i;

  X[0] = Start;
  X[N_V-1] = End;

//   for(i=0; i<N_V; i++)
//    cout<< " X[i] " << X[i] <<endl;

  Lines = new int[2*N_Cells];
  Vetrex = new TVertex*[N_V]; 

  y=0.;
  z=0.;
  for(i=0; i<N_Cells; i++)
   {
    Lines[2*i]=i;
    Lines[2*i+1]=i+1;
    Vetrex[i] = new TVertex(X[i], y, z);
   }

  Vetrex[N_Cells] = new TVertex(X[N_V-1], y, z);

  CellTree = new TBaseCell*[N_Cells];

   for (i=0;i<N_Cells;i++)
   {
//     Vetrex[ i ]->GetCoords(x, y);
//     cout<< " x " << x<< " y " << y<<endl;
    CellTree[i] = new TMacroCell(TDatabase::RefDescDB[S_Line], 0);
    CellTree[i]->SetVertex(0, Vetrex[ Lines[ 2*i       ]]);
    CellTree[i]->SetVertex(1, Vetrex[ Lines[ 2*i + 1]]);
    ((TMacroCell *) CellTree[i])->SetSubGridID(0);

//  cout<< " x " <<CellTree[i]->GetN_Edges()<<endl;;
//  cout<< " x " <<TDatabase::RefDescDB[S_Line]->GetN_OrigEdges()<<endl;;

   }


//     Vetrex[ i ]->GetCoords(x, y);
//     cout<< " x " << x<< " y " << y<<endl;
//     exit(0);

   Domain->SetTreeInfo(CellTree, N_Cells);

   TDatabase::IteratorDB[It_EQ]->SetParam(Domain);
   TDatabase::IteratorDB[It_LE]->SetParam(Domain);
   TDatabase::IteratorDB[It_Finest]->SetParam(Domain);
   TDatabase::IteratorDB[It_Between]->SetParam(Domain);
   TDatabase::IteratorDB[It_OCAF]->SetParam(Domain);

   // start joint(vertex)
   Joint = new TJointEqN(CellTree[0]);
   CellTree[0]->SetJoint(0, Joint);


   for(i=1;i<N_Cells;i++)
    {
     Joint = new TJointEqN(CellTree[i-1], CellTree[i]);

     CellTree[i-1]->SetJoint(1, Joint);
     CellTree[i]->SetJoint(0, Joint);
   } // for(i=0;i<N_Cells;i++)

   // end joint(vertex)
   Joint = new TJointEqN(CellTree[N_Cells-1]);
   CellTree[N_Cells-1]->SetJoint(1, Joint);

  delete []  Lines;
}



