// ======================================================================
// instationary problem
// ======================================================================
#include <TimeConvDiff2D.h>

void ExampleFile()
{
  
#define __LEVELSET__ 

  TDatabase::ParamDB->INTERNAL_CONVECTION_EQ_VELOFIELD=0;
    
  OutPut("Example: Zalesak rotating disk example Levelset.h" << endl);
}

// exact solution
void Exact(double x, double y, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

// exact solution
void InitialCondition(double x, double y, double *values)
{
//  double phi =  sqrt(x*x +(y-25.)*(y-25.)) - 15.;
//   
//   if(phi<=0)
//    { 
//     if(fabs(x)<=2.5 && y<=35)
//      { values[0] = 1.;}
//     else if(phi==0) 
//      { values[0] = 0.;}  
//     else
//     { values[0] = -1.;}
//    }
//   else
//    { values[0] = 1.;}   


  values[0] =  25. - (x*x + (y-40.)*(y-40.));
//  if(phi<0)
//   { values[0] = -1.;}
//   else if(phi==0) 
//    { values[0] = 0.;} 
//   else
//    { values[0] = 1.;}  
//  values[0] = phi;
}

// kind of boundary condition (for FE space needed)
void BoundCondition(int BdComp, double t, BoundCond &cond)
{
   cond = NEUMANN;
//     switch(BdComp)
//      {
//       case 0 : 
//       case 1:  
//       case 2 :        
//       case 3:
//          cond = NEUMANN;
//       break;  
//       
//        default:
//             OutPut("Unknown BoundCondition" << endl);
//             exit(4711);;
//      }
}

// value of boundary condition
void BoundValue(int BdComp, double Param, double &value)
{
  value = 0;
}

void BilinearCoeff(int n_points, double *X, double *Y,
        double **parameters, double **coeffs)
{
  double eps=1./TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff, *param, R;
  double x, y;
  double k0, alpha,  char_L;
  
  alpha = TDatabase::ParamDB->P1;
  char_L = TDatabase::ParamDB->P4;
  k0 = TDatabase::ParamDB->P5;
  R =  Pi/314.0;
     
  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    param = parameters[i];

    x = X[i];
    y = Y[i];

//     R = sqrt(x*x + y*y);
   
    coeff[0] = 0; //diffusion
    
//     if(fabs(R)<1.e-8)
//      {
//       coeff[1] = 0;
//       coeff[2] = 0;
//      }
//     else
//      {
      coeff[1] = R*(-y);
      coeff[2] = R*(x);
//      }      
//       coeff[1] = 0.;
//       coeff[2] = -1.;
      
//     cout << " u1 : " << coeff[1] << " u2 : " << coeff[2] << endl;
    coeff[3] = 0;
    coeff[4] = 0;
  }
}




void GetExampleFileData(BoundCondFunct2D **BoundaryConditions, BoundValueFunct2D **BoundValues, 
                        DoubleFunct2D **InitiaValues, CoeffFct2D **BilinearCoeffs, 
                        int &N_PBEqns, int &N_IndepntScalarEqns, int *Disctypes)
{

  TDatabase::ParamDB->INTERNAL_CONVECTION_EQ_VELOFIELD=0;


   N_IndepntScalarEqns = 1;
   N_PBEqns = 0;

   BilinearCoeffs[0] =  BilinearCoeff;
   BoundaryConditions[0] = BoundCondition;
   BoundValues[0] = BoundValue;
   InitiaValues[0] = InitialCondition;
   Disctypes[0] = GALERKIN;
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
  const int max_entries = 500000;
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
      if ( ( (fabs(x[j]) - 50.)<eps)||( (fabs(y[j]) - 50.)<eps) )
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
            if (fabs(y[j] + 50.)<eps)
            {
              tmp_bdry[diri_counter] = 0;
              tmp_param[diri_counter] = (50.+ x[j])/100.;
            }
            if (fabs(50.-y[j])<eps)
            {
              tmp_bdry[diri_counter] = 2;
              tmp_param[diri_counter] = (50.-x[j])/100.;
            }
            if (fabs(x[j] + 50.)<eps)
            {
              tmp_bdry[diri_counter] = 3;
              tmp_param[diri_counter] = (50.-y[j])/100.;
            }
            if (fabs(50.-x[j])<eps)
            {
              tmp_bdry[diri_counter] = 1;
              tmp_param[diri_counter] = (50. + y[j])/100.;
            }
// 	    OutPut( tmp_diri[diri_counter] << " " <<
// 		    tmp_bdry[diri_counter] << " " << tmp_param[diri_counter] << endl);
            diri_counter++;
          }
        }
// 	OutPut(endl);
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

//   for (i=0;i<diri_counter_1;i++)
//   {
//     OutPut(i << " " << neum_to_diri[i] << " " << neum_to_diri_bdry[i] <<
//       " " << neum_to_diri_param[i] <<  endl);
//   }
}


