// ======================================================================
// two boundary layer problem, from [JMT97]
// ======================================================================
// #include <ConvDiff2D.h>
#define __TWO_INTERIOR_LAYERS__

void ExampleFile()
{
  OutPut("Example: TwoBoundaryLayers.h" << endl) ;
  TDatabase::ParamDB->INTERNAL_QUAD_RULE = 99;
}

// exact solution
void Exact(double x, double y, double *values)
{
  double Reynolds = TDatabase::ParamDB->RE_NR;
  double t29,t33,t37,t41,t43,t44,t45,t46,t48;

  values[0]= x*y*y -y*y*exp(-2.0*(1.0-x)*Reynolds)
       -x*exp(-3.0*(1.0-y)*Reynolds)
       +exp(-(5.0-2.0*x-3.0*y)*Reynolds);

  t29 = y*y;
  t33 = exp(-2.0*(1.0-x)*Reynolds);
  t37 = exp(-3.0*(1.0-y)*Reynolds);
  t41 = exp(-(5.0-2.0*x-3.0*y)*Reynolds);
  t43 = t29*Reynolds*t33;
  t44 = Reynolds*t41;
  values[1]  = t29-2.0*t43-t37+2.0*t44;
  t45 = x*y;
  t46 = y*t33;
  t48 = x*Reynolds*t37;
  values[2]  = 2.0*t45-2.0*t46-3.0*t48+3.0*t44;  

  values[3] = 0; // laplace
}

// exact solution
void ExactAlternative(double x, double y, double *values)
{
  double eps = 1.0/TDatabase::ParamDB->RE_NR;
 
  values[0]= x*y*y 
      - y*y*exp(2.0*(x - 1.0)/eps)
      - x*exp(3.0*(y - 1.0)/eps)
      + exp((-5.0+2.0*x+3.0*y)/eps);

  values[1] = y*y 
      - y*y*exp(2.0*(x - 1.0)/eps) * 2/eps
      - exp(3.0*(y - 1.0)/eps)
      + exp((-5.0+2.0*x+3.0*y)/eps) * 2/eps;

  values[2] = 2*x*y 
      - 2*y*exp(2.0*(x - 1.0)/eps) 
      - x*exp(3.0*(y - 1.0)/eps) * 3/eps
      + exp((-5.0+2.0*x+3.0*y)/eps) * 3/eps;

  values[3] = 0; // laplace
} 

// kind of boundary condition (for FE space needed)
void BoundCondition(int i, double t, BoundCond &cond)
{
  cond = DIRICHLET;
}

// value of boundary condition
void BoundValue(int BdComp, double Param, double &value)
{
  static double Reynolds = TDatabase::ParamDB->RE_NR;
  double x;

  switch(BdComp)
  {
    case 0: value = -Param*exp(-3*Reynolds) + exp((-5+2*Param)*Reynolds);
            break;
    case 1: value=0;
            break;
    case 2: value=0;
            break;
    case 3: x=1.0-Param;
            value = -exp(-2.0*Reynolds)+2.0*exp(-2.0*Reynolds)*Param
                    -exp(-2.0*Reynolds)*Param*Param
                    +exp(-1.0*(2.0+3.0*Param)*Reynolds);
            break;
    default: cout << "wrong boundary part number" << endl;
            break;
  }
}

void BilinearCoeffs(int n_points, double *X, double *Y,
        double **parameters, double **coeffs)
{
  static double eps=1/TDatabase::ParamDB->RE_NR;
  static double Reynolds=TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff, *param, val[4];
  double x, y;
  double u,t29,t33,t37,t41,t43,t44,t45,t46,t48,t49,t51,t52,t54;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    param = parameters[i];

    coeff[0] = eps; 
    coeff[1] = 2; // b1
    coeff[2] = 3; // b2
    coeff[3] = 1; // c

    x = X[i];
    y = Y[i];

    if (x<0) x=0;
    if (y<0) y=0;
    if (x>1) x=1;
    if (y>1) y=1;

    t29 = y*y;
    t33 = exp(-2.0*(1.0-x)*Reynolds);
    t37 = exp(-3.0*(1.0-y)*Reynolds);
    t41 = exp(-(5.0-2.0*x-3.0*y)*Reynolds);
//    cout << x << t33 << " " << t37 << " " << t41 << endl;
    u = x*t29-t29*t33-x*t37+t41;
    t43 = t29*Reynolds*t33;
    t44 = Reynolds*t41;
    t45 = x*y;
    t46 = y*t33;
    t48 = x*Reynolds*t37;
    t49 = Reynolds*Reynolds;
    t51 = t29*t49*t33;
    t52 = t49*t41;
    t54 = x*t49*t37;

    // right-hand side
    coeff[4]= -(-4.0*t51+13.0*t52+2.0*x-2.0*t33-9.0*t54)/Reynolds
      +2.0*t29-4.0*t43-2.0*t37+13.0*t44+6.0*t45-6.0*t46-9.0*t48+u;

    coeff[5] = sqrt(13); //   \|b\|_infty
      
     // values of solution of continuous problem, just for testing
    Exact(x, y, val);
    coeff[10] = val[0];
    coeff[11] = val[1];
    coeff[12] = val[2];    
  }
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

void CheckWrongNeumannNodes(TCollection *Coll, TFESpace2D *fespace,
			    int &N_neum_to_diri, int* &neum_to_diri,
			    int* &neum_to_diri_bdry, 
			    double* &neum_to_diri_param)
{
    const int max_entries = 4096;  
    int i, j, N_, min_val, type;
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
	// P_1, Q_1
	case C_P2_2D_T_A:
	case C_Q2_2D_Q_A:
	case C_Q2_2D_Q_M:
            // loop over the edges
 	    for (j=0;j<N_V;j++)
	    {
              // check of edge j is on boundary  
              if (boundary_vertices[j] && boundary_vertices[(j+1)%N_V])
              {
		// check if this is a boundary edge
		type = cell->GetJoint(j)->GetType();
		if (!((type == BoundaryEdge)||(type == IsoBoundEdge)))
		  continue;
	        switch(j)
                {
                   case 0:
                     tmp_diri[diri_counter] = dof[0];
                     tmp_diri[diri_counter+1] = dof[1];
                     tmp_diri[diri_counter+2] = dof[2];
                   break;
                  case 1:
                     if (N_V==3)
                     {
                       tmp_diri[diri_counter] = dof[2];
                       tmp_diri[diri_counter+1] = dof[4];
                       tmp_diri[diri_counter+2] = dof[5];
                     }
                     else
                     {
                       tmp_diri[diri_counter] = dof[2];
                       tmp_diri[diri_counter+1] = dof[5];
                       tmp_diri[diri_counter+2] = dof[8];
                     }
                   break;
                  case 2:
                     if (N_V==3)
                     {
                       tmp_diri[diri_counter] = dof[5];
                       tmp_diri[diri_counter+1] = dof[3];
                       tmp_diri[diri_counter+2] = dof[0];
                     }
                     else
                     {
                       tmp_diri[diri_counter] = dof[8];
                       tmp_diri[diri_counter+1] = dof[7];
                       tmp_diri[diri_counter+2] = dof[6];
                     }
                   break;
                   case 3:
                     tmp_diri[diri_counter] = dof[6];
                     tmp_diri[diri_counter+1] = dof[3];
                     tmp_diri[diri_counter+2] = dof[0];
                   break;

                }
              
		if (diri_counter+2 > max_entries)
		{
			OutPut("tmp_diri too short !!!"<<endl);
			exit(4711);
		}


		    if ((fabs(y[j])<eps)&&(fabs(y[(j+1)%N_V])<eps))
		    {
			tmp_bdry[diri_counter] = 0;
			tmp_bdry[diri_counter+1] = 0;
			tmp_bdry[diri_counter+2] = 0;
			tmp_param[diri_counter] = x[j];
			tmp_param[diri_counter+1] = (x[j] + x[(j+1)%N_V])/2.0;
			tmp_param[diri_counter+2] = x[(j+1)%N_V];
		    }
		    if ((fabs(1-y[j])<eps)&&(fabs(1-y[(j+1)%N_V])<eps))
		    {
			tmp_bdry[diri_counter] = 2;
			tmp_bdry[diri_counter+1] = 2;
			tmp_bdry[diri_counter+2] = 2;
			tmp_param[diri_counter] = 1-x[j];
			tmp_param[diri_counter+1] = (1-x[j] + 1-x[(j+1)%N_V])/2.0;
			tmp_param[diri_counter+2] = 1-x[(j+1)%N_V];
		    }
		    if ((fabs(x[j])<eps)&&(fabs(x[(j+1)%N_V])<eps))
		    {
			tmp_bdry[diri_counter] = 3;
			tmp_bdry[diri_counter+1] = 3;
			tmp_bdry[diri_counter+2] = 3;
			tmp_param[diri_counter] = 1-y[j];
			tmp_param[diri_counter+1] = (1-y[j] + 1-y[(j+1)%N_V])/2.0;
			tmp_param[diri_counter+2] = 1-y[(j+1)%N_V];
		    }
		    if ((fabs(1-x[j])<eps)&&(fabs(1-x[(j+1)%N_V])<eps))
		    {
			tmp_bdry[diri_counter] = 1;
			tmp_bdry[diri_counter+1] = 1;
			tmp_bdry[diri_counter+2] = 1;
			tmp_param[diri_counter] = y[j];
			tmp_param[diri_counter+1] = (y[j] + y[(j+1)%N_V])/2.0;
			tmp_param[diri_counter+2] = y[(j+1)%N_V];
		    }
		    diri_counter +=3;
		}
	    }
	    break;
	// P_3, Q_3
	case C_P3_2D_T_A:
	case C_Q3_2D_Q_A:
	case C_Q3_2D_Q_M:
            // loop over the edges
 	    for (j=0;j<N_V;j++)
	    {
              // check of edge j is on boundary  
              if (boundary_vertices[j] && boundary_vertices[(j+1)%N_V])
              {
		// check if this is a boundary edge
		type = cell->GetJoint(j)->GetType();
		if (!((type == BoundaryEdge)||(type == IsoBoundEdge)))
		  continue;

               // P3: local dof 0, 1, 2, 3 are on the boundary
               // Q3: local dof 0, 1, 2, 3 are on the boundary
	        switch(j)
                {
                   case 0:
                     tmp_diri[diri_counter] = dof[0];
                     tmp_diri[diri_counter+1] = dof[1];
                     tmp_diri[diri_counter+2] = dof[2];
		     tmp_diri[diri_counter+3] = dof[3];
                   break;
                  case 1:
                     if (N_V==3)
                     {
                       tmp_diri[diri_counter] = dof[3];
                       tmp_diri[diri_counter+1] = dof[6];
                       tmp_diri[diri_counter+2] = dof[8];
		       tmp_diri[diri_counter+3] = dof[9];
                     }
                     else
                     {
                       tmp_diri[diri_counter] = dof[3];
                       tmp_diri[diri_counter+1] = dof[7];
                       tmp_diri[diri_counter+2] = dof[11];
		       tmp_diri[diri_counter+3] = dof[15];
                     }
                   break;
                  case 2:
                     if (N_V==3)
                     {
                       tmp_diri[diri_counter] = dof[9];
                       tmp_diri[diri_counter+1] = dof[7];
                       tmp_diri[diri_counter+2] = dof[4];
                       tmp_diri[diri_counter+3] = dof[0];
		     }
                     else
                     {
                       tmp_diri[diri_counter] = dof[15];
                       tmp_diri[diri_counter+1] = dof[14];
                       tmp_diri[diri_counter+2] = dof[13];
			tmp_diri[diri_counter+3] = dof[12];
                     }
                   break;
                   case 3:
                     tmp_diri[diri_counter] = dof[12];
                     tmp_diri[diri_counter+1] = dof[8];
                     tmp_diri[diri_counter+2] = dof[4];
		     tmp_diri[diri_counter+3] = dof[0];
                   break;
                }
              
		if (diri_counter+3 > max_entries)
		{
			OutPut("tmp_diri too short !!!"<<endl);
			exit(4711);
		}

		if ((fabs(y[j])<eps)&&(fabs(y[(j+1)%N_V])<eps))
		{
		    tmp_bdry[diri_counter] = 0;
		    tmp_bdry[diri_counter+1] = 0;
		    tmp_bdry[diri_counter+2] = 0;
		    tmp_bdry[diri_counter+3] = 0;
		    tmp_param[diri_counter] = x[j];
		    tmp_param[diri_counter+1] = 2*x[j]/3.0 + x[(j+1)%N_V]/3.0;
		    tmp_param[diri_counter+2] = x[j]/3.0 + 2*x[(j+1)%N_V]/3.0; 
		    tmp_param[diri_counter+3]= x[(j+1)%N_V];
		}
		if ((fabs(1-y[j])<eps)&&(fabs(1-y[(j+1)%N_V])<eps))
		{
		    tmp_bdry[diri_counter] = 2;
		    tmp_bdry[diri_counter+1] = 2;
		    tmp_bdry[diri_counter+2] = 2;
		    tmp_bdry[diri_counter+3] = 2;
		    tmp_param[diri_counter] = 1-x[j];
		    tmp_param[diri_counter+1] = 2*(1-x[j])/3.0 + (1-x[(j+1)%N_V])/3.0;
		    tmp_param[diri_counter+2] = (1-x[j])/3.0 + 2*(1-x[(j+1)%N_V])/3.0;
		    tmp_param[diri_counter+3] = 1-x[(j+1)%N_V];
		}
		    if ((fabs(x[j])<eps)&&(fabs(x[(j+1)%N_V])<eps))
		    {
			tmp_bdry[diri_counter] = 3;
			tmp_bdry[diri_counter+1] = 3;
			tmp_bdry[diri_counter+2] = 3;
			tmp_bdry[diri_counter+3] = 3;
			tmp_param[diri_counter] = 1-y[j];
			tmp_param[diri_counter+1] = 2*(1-y[j])/3.0 + (1-y[(j+1)%N_V])/3.0;
			tmp_param[diri_counter+2] = (1-y[j])/3.0 + (1-y[(j+1)%N_V])/3.0;
			tmp_param[diri_counter+3] = 1-y[(j+1)%N_V];
		    }
		    if ((fabs(1-x[j])<eps)&&(fabs(1-x[(j+1)%N_V])<eps))
		    {
			tmp_bdry[diri_counter] = 1;
			tmp_bdry[diri_counter+1] = 1;
			tmp_bdry[diri_counter+2] = 1;
			tmp_bdry[diri_counter+3] = 1;
			tmp_param[diri_counter] = y[j];
			tmp_param[diri_counter+1] = 2*y[j]/3.0 + y[(j+1)%N_V]/3.0;
			tmp_param[diri_counter+2] = y[j]/3.0 + 2*y[(j+1)%N_V]/3.0;
			tmp_param[diri_counter+3] = y[(j+1)%N_V];
		    }
		    diri_counter +=4;
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

void ComputeExtremalValues(int N, double *sol, double  *values)
{
   int i;
   double max, min,val[4];

   min = 1e10;
   max = -1e10;
   
   for(i=0;i<N;i++)
   {
      if(sol[i] > max)
         max = sol[i];
      if(sol[i] < min)
         min = sol[i];
   }
   // maximal value
   Exact(0.5,0.5,val);

   values[0] = min;
   values[1] = max;
   //OutPut("undershoots " << errors[0] << " overshoots " << errors[1] << endl);
}
