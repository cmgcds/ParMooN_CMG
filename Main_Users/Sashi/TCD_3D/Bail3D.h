// ==========================================================================
// instationary problem
// ==========================================================================

//===========================================================================
// example file
// =========================================================================
#define __BAIL_3D__

void ExampleFile()
{
  OutPut("Example: Bail3D.h" << endl);
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

// initial conditon
void InitialCondition(double x, double y, double z, double *values)
{
    values[0] = 0;
}

// kind of boundary condition
void BoundCondition(double x, double y, double z, BoundCond &cond)
{
    double eps = 1e-6;

  if (TDatabase::ParamDB->SOLD_PARAMETER_TYPE== FEM_FCT)
    cond = NEUMANN;
  else
  {
      if ((fabs(x-1)<eps)&& (y>=3.0/8.0-eps) && (y<=0.5+eps)
	  && (z>=0.5-eps) && (z<=5.0/8.0+eps))
	  cond = NEUMANN;
      else
	  cond = DIRICHLET;

/*
      if ((fabs(x)<eps) || (fabs(y)<eps) ||  (fabs(z)<eps) ||
	  (fabs(1-y)<eps) || (fabs(1-z)<eps))  
	  cond = DIRICHLET;
      else
	  // at x = 1
	  cond = NEUMANN;
*/
  }
}

// value of boundary condition
void BoundValue(double x, double y, double z, double &value)
{

  double t = TDatabase::TimeDB->CURRENTTIME;
    double eps = 1e-6;
  
    /*if (TDatabase::ParamDB->SOLD_PARAMETER_TYPE== FEM_FCT)
    {
      value = 0.0;
     }
    else*/
 
    {
    if ((x==0) && (y>=5.0/8.0-eps) && (y<=0.75+eps) &&
	(z>= 5.0/8.0-eps) && (z<= 0.75 +eps))
    {
	if ((t>=0) && t<=1)
	    value =  sin(Pi*t/2.0)*1.0;
	else
	{
	    if ((t>1)&&t<=2)
		value = 1.0;
	    else
	    {
		if ((t>2)&&(t<=3))
		    value = sin(Pi*(t-1)/2.0)*1.0;
		else
		{
		    value = 0;
		}
	    }
	}
    }
    else
	value = 0.0;
   }
}

void BilinearCoeffs(int n_points, double *X, double *Y, double *Z,
        double **parameters, double **coeffs)
{
  double eps = 1.0/TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff;                                  // *param;
  double x, y, z, c, a[3], b[3], s[3], h;
  // double t = TDatabase::TimeDB->CURRENTTIME;
  
  // norm of convection vector
  c = 1.0/16 + 1.0/64 + 1;
  c = sqrt(c);

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    // param = parameters[i];

    x = X[i];
    y = Y[i];
    z = Z[i];

    // diffusion
    coeff[0] = eps;
    // convection in x direction
    coeff[1] = 1;
    // convection in y direction
    coeff[2] = -0.25;
    // convection in z direction
    coeff[3] = -0.125;
    // reaction
    // vector from center of outflow (1,7/16,9/16) to (x,y,z)
    a[0] = x-1;
    a[1] = y - 7.0/16;
    a[2] = z - 9.0/ 16;
    // vector from center of inflow (0,11/16,11/16) to (x,y,z)
    b[0] = x;
    b[1] = y - 11.0/16;
    b[2] = z - 11.0/16;
    // cross product
    s[0] = a[1] * b[2] - a[2] * b[1];
    s[1] = a[2] * b[0] - a[0] * b[2];
    s[2] = a[0] * b[1] - a[1] * b[0];
    // norm of cross product = 2 * length of convection * height 
    // area of parallelogram = 2 * area of triangle = 2 * c * h /2 = c * h
    h = sqrt(s[0] * s[0] + s[1] * s[1] + s[2] * s[2]);
    // 2 * height
    h /= c;

    if (h<=0.1)
    {	coeff[4] = 1;
    OutPut(x << " " << y << " " << z << " " << h << " " << coeff[4] << endl);
    }
    else
	coeff[4] = 0;
    // rhs
    coeff[5] = 0;
    // rhs from previous time step
    coeff[6] = 0;
  }
}

/****************************************************************/
//
// for FEM_TVD
//
/****************************************************************/
void CheckWrongNeumannNodes(TCollection *Coll, TFESpace3D *fespace,
int &N_neum_to_diri, int* &neum_to_diri,
			    double* &neum_to_diri_x, double* &neum_to_diri_y, double* &neum_to_diri_z)
{
   const int max_entries = 50000;  
  int i, j, N_, min_val;
  int N_Cells, N_V, diri_counter = 0, found, diri_counter_1 = 0;
  int *global_numbers, *begin_index, *dof;
  int boundary_vertices[8], tmp_diri[max_entries];
  double x[8], y[8], z[8], eps = 1e-6, tmp_x[max_entries], tmp_y[max_entries], tmp_z[max_entries];
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
      // vertex on the upper lid
      if ((fabs(x[j])<eps)||(fabs(y[j])<eps)||(fabs(x[j]-1)<eps)||(fabs(y[j]-1)<eps)||(fabs(z[j])<eps)||(fabs(z[j]-1)<eps))
       {
          if ((fabs(x[j]-1)<eps)&& (y[j]>=3.0/8.0-eps) && (y[j]<=0.5+eps)
	      && (z[j]>=0.5-eps) && (z[j]<=5.0/8.0+eps))
           continue;
{
           // Dirichlet boundary
	   boundary_vertices[j] = 1;
	   found++;}
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
		    if (CurrentElement==C_P1_3D_T_A)
			tmp_diri[diri_counter] = dof[j];
		    else
		    {
			if ((j==0)||(j==1)||(j==4)||(j==5))
                         {
		   			    tmp_diri[diri_counter] = dof[j];
		  			}
			else
			{
			    if (j==2)
				tmp_diri[diri_counter] = dof[3];
			    if (j==3)
				tmp_diri[diri_counter] = dof[2];
			    if (j==6)
				tmp_diri[diri_counter] = dof[7];
			    if (j==7)
				tmp_diri[diri_counter] = dof[6];
			}
		    }
		    if (diri_counter > max_entries)
		    {
			OutPut("tmp_diri too short !!!"<<endl);
			exit(4711);
		    }
                    tmp_x[diri_counter] = x[j];
                    tmp_y[diri_counter] = y[j];
                    tmp_z[diri_counter] = z[j];
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
  // allocate array for the corresponding boundary coordinates
  neum_to_diri_x = new double[diri_counter_1];
  neum_to_diri_y = new double[diri_counter_1];
  neum_to_diri_z = new double[diri_counter_1];
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
      neum_to_diri_x[i] = tmp_x[found];
      neum_to_diri_y[i] = tmp_y[found];
      neum_to_diri_z[i] = tmp_z[found];
      tmp_diri[found] = -1;
  }

  for (i=0;i<diri_counter_1;i++)
  {
      OutPut(i << " " << neum_to_diri[i] << " " << neum_to_diri_x[i] <<
	     " " << neum_to_diri_y[i] << " " << neum_to_diri_z[i] <<  endl);
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
   values[1] = max;
}
