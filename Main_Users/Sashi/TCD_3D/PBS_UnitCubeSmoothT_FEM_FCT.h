#define __UREA__
#define __SIMPATURS__
  
#include <Urea_3d4d.h>
#include <MacroCell.h>

void ExampleFile()
{
  
  #define __PBSConstT__  
  #define __FEMFCT__ 

  OutPut("Example: PBS_UnitCubeSmoothT_FEM_FCT.h " << endl);
  
  TDatabase::ParamDB->UREA_SOLD_PARAMETER_TYPE = FEM_FCT_LIN;
  
  
  // for velocity switch
  TDatabase::ParamDB->INTERNAL_CONVECTION_EQ_VELOFIELD=0;

  TDatabase::ParamDB->N_CELL_LAYERS = 3;
  TDatabase::ParamDB->DRIFT_Z = 1;
  TDatabase::ParamDB->INTERNAL_PROBLEM_IDENTITY = 1356;
 
  
#ifdef _MPI
 MPI_Comm Comm;
 int rank;

 Comm = TDatabase::ParamDB->Comm;
 MPI_Comm_rank(Comm, &rank);

 if(rank==TDatabase::ParamDB->Par_P0)
#endif
 {  
  OutPut("Example: PBS.h " << endl);
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

void BoundCondition_FEMFCT(double x, double y, double z, BoundCond &cond)
{
 cond = NEUMANN;
}

//already defined in mai
// void BoundValue_FEM_FCT(double x, double y, double z, double &value)
// {
//  value = 0;   
// }

// ========================================================================
// definitions for the temperature
// ========================================================================

void Exact(double x, double y, double z, double *values)
{
 double t = TDatabase::TimeDB->CURRENTTIME;
 double k = 0.1;
  
  values[0] = (exp(-k*t))*sin(Pi*x)*cos(Pi*y)*cos(Pi*z);
  values[1] = Pi*(exp(-k*t))*cos(Pi*x)*cos(Pi*y)*cos(Pi*z);
  values[2] = -Pi*(exp(-k*t))*cos(Pi*x)*sin(Pi*y)*cos(Pi*z);
  values[3] = -Pi*(exp(-k*t))*cos(Pi*x)*cos(Pi*y)*sin(Pi*z);
  values[4] = -3.*Pi*Pi*(exp(-k*t))*sin(Pi*x)*cos(Pi*y)*cos(Pi*z);
}

// initial conditon
void InitialExact(double x, double y, double z, double *values)
{ 
 double t = 0;
 double k = 0.1;
  
  values[0] = (exp(-k*t))*sin(Pi*x)*cos(Pi*y)*cos(Pi*z);
  values[1] = Pi*(exp(-k*t))*cos(Pi*x)*cos(Pi*y)*cos(Pi*z);
  values[2] = -Pi*(exp(-k*t))*cos(Pi*x)*sin(Pi*y)*cos(Pi*z);
  values[3] = -Pi*(exp(-k*t))*cos(Pi*x)*cos(Pi*y)*sin(Pi*z);
  values[4] = -3.*Pi*Pi*(exp(-k*t))*sin(Pi*x)*cos(Pi*y)*cos(Pi*z);  
}


// initial conditon
void InitialCondition_temp(double x, double y, double z, double *values)
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

// kind of boundary condition (for FE space needed)
void BoundCondition_temp(double x, double y, double z, BoundCond &cond)
{
 cond = DIRICHLET;
}

// value of boundary condition
//void BoundValue_c_A(int BdComp, double Param, double &value)
void BoundValue_temp(double x, double y, double z, double &value)
{
 double t= TDatabase::TimeDB->CURRENTTIME;
 double k = 0.1;
 
 value = (exp(-k*t))*sin(Pi*x)*cos(Pi*y)*cos(Pi*z);   
}


// ========================================================================
// BilinearCoeffs for Heat 
// ========================================================================
void BilinearCoeffs_Heat(int n_points, double *x, double *y, double *z,
               double **parameters, double **coeffs)
{
  double eps;
  int i;
  double *coeff, *param;
  double t = TDatabase::TimeDB->CURRENTTIME;
 
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
      coeff[3] = param[2];  // u2      
//       cout<< "coeff[0] eps " << eps << " u1 " << param[0]  << "u2 " << param[1] << endl;
     }
    else
     {
      coeff[1] = 0;  // u1
      coeff[2] = 0;  // u2
      coeff[3] = 0;  // u2      
     }
     
    // reaction
    coeff[4] = 0;
     // rhs
    coeff[5] = (3.*eps*Pi*Pi - 0.1)*(exp(-0.1*t))*sin(Pi*x[i])*cos(Pi*y[i])*cos(Pi*z[i]); // f
    coeff[6] = 0;  
//     cout<< "coeff[0] eps " << eps << " u1 " << param[0]  << "u2 " << param[1] << endl;
  }
}
// ========================================================================
// definitions for the PSD
// ========================================================================

void Exact_Psd_Intl(double x, double y, double z, double l, double *values)
{
 double k = 0.1;
 double t= TDatabase::TimeDB->CURRENTTIME; 
 
  values[0] = (exp(-k*t))*sin(Pi*x)*cos(Pi*y)*cos(Pi*z)*sin(Pi*l);   
  values[1] = Pi*(exp(-k*t))*cos(Pi*x)*cos(Pi*y)*cos(Pi*z)*sin(Pi*l);    
  values[2] = -Pi*(exp(-k*t))*sin(Pi*x)*sin(Pi*y)*cos(Pi*z)*sin(Pi*l);
  values[3] = -Pi*(exp(-k*t))*sin(Pi*x)*cos(Pi*y)*sin(Pi*z)*sin(Pi*l); 
  values[4] =  Pi*(exp(-k*t))*sin(Pi*x)*cos(Pi*y)*cos(Pi*z)*sin(Pi*l);   
  values[5] = -4*Pi*Pi*(exp(-k*t))*sin(Pi*x)*cos(Pi*y)*cos(Pi*z)*sin(Pi*l);    
}

void Exact_psd(int N_Coord, double *X, double *values)
{
 double t = TDatabase::TimeDB->CURRENTTIME; 
 
 if(N_Coord!=4)
  {
   printf("N_Coord!=4 InitialCondition_psd_Intl !!!\n");  
#ifdef _MPI
     MPI_Finalize();
#endif  
     exit(0);   
  }
 
  values[0] =  (exp(-0.1*t))*sin(Pi*X[0])*cos(Pi*X[1])*cos(Pi*X[2])*sin(Pi*X[3]);
}



// initial condition
void InitialCondition_psd(double x, double y, double z, double *values)
{
 double t = 0;  
 double k = 0.1;  
 double l = TDatabase::ParamDB->REACTOR_P29;
  
 values[0] = (exp(-k*t))*sin(Pi*x)*cos(Pi*y)*cos(Pi*z)*sin(Pi*l);
}

void InitialCondition_psd_Intl(int N_Coord, double *X, double *values)
{
 double t = 0;
  
 if(N_Coord!=4)
  {
   printf("N_Coord!=4 InitialCondition_psd_Intl !!!\n");  
#ifdef _MPI
     MPI_Finalize();
#endif  
     exit(0);   
  }
 
  values[0] = (exp(-0.1*t))*sin(Pi*X[0])*cos(Pi*X[1])*cos(Pi*X[2])*sin(Pi*X[3]);
}



// kind of boundary condition (for FE space needed)
void BoundCondition_psd(double x, double y, double z, BoundCond &cond)
{
   cond = DIRICHLET;  
}

// value of boundary condition
void BoundValue_psd(double x, double y, double z, double &value)
{
  double k = 0.1, temp;
  double l = TDatabase::ParamDB->REACTOR_P29;
  double t= TDatabase::TimeDB->CURRENTTIME;    

  value = (exp(-k*t))*sin(Pi*x)*cos(Pi*y)*cos(Pi*z)*sin(Pi*l);  
}


void BoundCondition_LminLMax(BoundCond &cond_Lmin, BoundCond &cond_Lmax)
{
  cond_Lmin = DIRICHLET;
//   cond_Lmin = NEUMANN;
  cond_Lmax = DIRICHLET;
}

void BoundValue_LMin(double x, double y, double z, double *values)
 {
  double k = 0.1;
  double l = TDatabase::ParamDB->REACTOR_P12;
  double t= TDatabase::TimeDB->CURRENTTIME;
 
  values[0] = (exp(-k*t))*sin(Pi*x)*cos(Pi*y)*cos(Pi*z)*sin(Pi*l); 
 }


void BoundValue_LMax(double x, double y,  double z,  double *values)
 {
  double k = 0.1;
  double l = TDatabase::ParamDB->REACTOR_P13;
  double t= TDatabase::TimeDB->CURRENTTIME;
 
  values[0] = (exp(-k*t))*sin(Pi*x)*cos(Pi*y)*cos(Pi*z)*sin(Pi*l); 
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
  double t = TDatabase::TimeDB->CURRENTTIME;
  double l = TDatabase::ParamDB->REACTOR_P29;
  
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
      coeff[3] = param[2];  // u2      
     }
    else
     {
      coeff[1] = 0;  // u1
      coeff[2] = 0;  // u2
      coeff[3] = 0;  // u2      
     }
     
    // reaction
    coeff[4] = 0;
     // rhs
    coeff[5] = (3.*eps*Pi*Pi - 0.1)*(exp(-0.1*t))*sin(Pi*x[i])*cos(Pi*y[i])*cos(Pi*z[i])*sin(Pi*l); // f
    coeff[6] = 0;  
  }
}


void BilinearCoeffs_Psd_Intl(int n_points, int N_Dim, double **Coords,
                             double **parameters, double **coeffs)
{
  int i;
  double eps, *coeff; // *param;
  double x, y, z, L, c, a[3], b, s[3], h;
  double t = TDatabase::TimeDB->CURRENTTIME;

//   b = -1e8;// negative, so that C will be taken from the PBS growth term
  b = 0;
  c = 0; 

  if(TDatabase::ParamDB->REACTOR_P3)
    eps = 1.0/TDatabase::ParamDB->REACTOR_P3;
  else
    eps = 0.;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    x = Coords[0][i];
    y = Coords[1][i];
    z = Coords[2][i];
    L = Coords[3][i];    
   
    // diffusion
    coeff[0] = eps;   
    
    // convection in z direction
    coeff[1] = b;
    // reaction term
    coeff[2] = c;
    coeff[3] = eps*Pi*Pi*(exp(-0.1*t))*sin(Pi*x)*cos(Pi*y)*cos(Pi*z)*sin(Pi*L); // rhs
  }
}

void GetExampleFileData(BoundCondFunct3D **BoundaryConditions, BoundValueFunct3D **BoundValues, 
                        DoubleFunct3D **InitiaValues, CoeffFct3D **BilinearCoeffs, 
                        int &N_PBEqns, int &N_IndepntScalarEqns, int *Disctypes)
{
  
   N_IndepntScalarEqns = 1;
   N_PBEqns = 0;
   
//      #define __PBS__

   BilinearCoeffs[0] = BilinearCoeffs_Heat;
//    BilinearCoeffs[0] = BilinearCoeffs_Psd;

   BoundaryConditions[0] = BoundCondition_FEMFCT;
//    BoundaryConditions[0] = BoundCondition_psd;

   BoundValues[0] = BoundValue_temp;
//    BoundValues[0] = BoundValue_psd;

   InitiaValues[0] = InitialExact; 
//    InitiaValues[0] = InitialCondition_psd;

   Disctypes[0] = GALERKIN;
 
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


/****************************************************************/
/* finds the nodes which are Neumann and should be Dirichlet    */
/* for FEM_FCT schemes                                          */
/****************************************************************/

// assumed that all boundary faces are of Dirichlet type in Unit cube
void CheckWrongNeumannNodes(TCollection *Coll, TFESpace3D *fespace,
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
      if ( (fabs(x[j])<eps) || (fabs(y[j])<eps) || (fabs(z[j])<eps)  || 
	   (fabs(x[j]-1)<eps) ||(fabs(y[j]-1)<eps) || (fabs(z[j]-1)<eps))
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
               if ( (fabs(x[j])<eps) || (fabs(y[j])<eps) || (fabs(z[j])<eps)  || 
	            (fabs(x[j]-1)<eps) ||(fabs(y[j]-1)<eps) || (fabs(z[j]-1)<eps))
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
	    OutPut("CheckWrongNeumannNodes not implemented for element "
		   << CurrentElement << endl);
	    OutPut("code can be run without CheckWrongNeumannNodes, just delete the exit" << endl);
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
  
  OutPut("CheckWrongNeumannNodes: N_neum_to_diri " << diri_counter_1 << endl);
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

