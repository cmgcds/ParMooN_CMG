// Navier-Stokes problem, Driven cavity
//
// u(x,y) = unknown
// p(x,y) = unknown
#include <Joint.h>
#include <BoundFace.h>

#define __CALOTTE__

void ExampleFile()
{
  OutPut("Example: Calotte.h, rotation " << TDatabase::ParamDB->P7
	 << ", " << TDatabase::ParamDB->P8 << " % noise (only for U3 !!!)"
	 << ", inflow scale with " << TDatabase::ParamDB->P9 
	 << endl);
}


void InitialU1(double x, double y, double z, double *values)
{
  double w, t1, t2, n,eps=1e-3;

  values[0] = 0;
}


void InitialU2(double x, double y, double z, double *values)
{
  values[0] = 0;
}


void InitialU3(double x, double y, double z, double *values)
{
  double w, t1, t2, n;
  values[0] = 0;
/*  if ((x-5*z/28) * (x-5*z/28) + y*y<= 6.25)
  {
    values[0] = -1;
  }
  if ((x+5*z/28) * (x+5*z/28) + y*y<= 6.25)
  {
    values[0] = -1;
    }*/
}


void InitialP(double x, double y,  double z, double *values)
{
  values[0] = 0;
}


// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y,  double z, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}


void ExactU2(double x, double y,  double z, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}


void ExactU3(double x, double y,  double z, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}


void ExactP(double x, double y,  double z, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}


// kind of boundary condition (for FE space needed)
void BoundCondition(double x, double y, double z, BoundCond &cond)
{
  double eps = 1e-6;

  cond = DIRICHLET;

  // this is for the nodes after the first refinement
  if ((z<10)&&(fabs(x*x+y*y)<=3.465*6.25))
  {
      cond = NEUMANN;
      TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 0;
  }
}


// value of boundary condition
void U1BoundValue(double x, double y, double z, double &value)
{
  value = 0;
}


// value of boundary condition
void U2BoundValue(double x, double y, double z, double &value)
{
  value = 0;
}


// value of boundary condition
void U3BoundValue(double x, double y, double z, double &value)
{
  double eps = 1e-2;
  double noise = TDatabase::ParamDB->P8/100.0;

  value = 0;
  if (fabs(z-140)<eps)
  {
    if ((x+25)*(x+25)+y*y <= 6.25+eps)
      value = (-1.0+ noise * ((double)rand()/RAND_MAX-0.5))*TDatabase::ParamDB->P9;
    if ((x-25)*(x-25)+y*y <= 6.25+eps)
      value = (-1.0+ noise * ((double)rand()/RAND_MAX-0.5))*TDatabase::ParamDB->P9;
  }
  //if ((fabs(z)<10)&& (x*x+y*y <= 2*6.25+eps))
  //    value = -1.0;
}


// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *X, double *Y, double *Z,
double **parameters, double **coeffs)
{
  static double eps = 1/TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff, x, y, z, r, angle, angle1, phi;
  double t = TDatabase::TimeDB->CURRENTTIME;

  // current angle of the rotor
  // note: the rotor has two blades
  angle = -2*Pi*t*TDatabase::ParamDB->P7;
  // normalize to (-Pi,Pi]
  while (angle<=-Pi)
  {
      angle += 2*Pi;
  }
  if (angle<=0)
  {
      angle1 = angle + Pi;
  }
  else
  {
      angle1 = angle - Pi;
  }

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    x = X[i];
    y = Y[i];
    z = Z[i];

    coeff[0] = eps;
    // f1
    coeff[1] = 0;
    // f2
    coeff[2] = 0;
    // f3
    coeff[3] = 0;

    // compute angle of the point 
    r = sqrt(x*x+y*y);
    if ((z<=75)&&(z>=65)&&(r<=25))
    {
	// -pi < angle_1 <= pi 
	phi = atan2(y,x);
	// two blades, the jump of angle_1 at the negative x-axis has
	// to be taken into account
	// change right hand side
	// driving force b=2*Pi*P7*(y,-x,0)
	// the terms on the rhs are (b*\nabla)b
	if ((fabs(angle-phi)<0.175)|| (fabs(angle1-phi)<0.175)
	    || (fabs(angle-phi+2*Pi)<0.175)|| (fabs(angle1-phi+2*Pi)<0.175) 
	    || (fabs(angle-phi-2*Pi)<0.175)|| (fabs(angle1-phi-2*Pi)<0.175) )
	{
	    //OutPut("rot " << x << " " << y << " " << z << endl);
	    coeff[1] = y*2*Pi*TDatabase::ParamDB->P7;
	    coeff[2] = -x*2*Pi*TDatabase::ParamDB->P7;
	}
    }
  }
}

// ========================================================================
//
// check if there are Neumann nodes which should be Diriclet nodes
//
// ========================================================================

void CheckNeumannNodesForVelocity(TCollection *Coll, TFESpace3D *fespace,
				  int &N_neum_to_diri, int* &neum_to_diri, 
				  int &N_neum_bdry, int* &neum_bdry)
{
  int i, j, N_Cells, N_Active, N_V, found, N_, diri[MaxN_BaseFunctions2D];
  int *global_numbers, *begin_index, *dof,  *N_BaseFunct, diri_counter=0;
  int tmp_diri[10000], diri_counter_1 = 0, min_val, face_lid[8], no_face_lid;
  int face_lid_index[4], neum_bdry_counter=0, neum_bdry_counter_1=0, tmp_neum_bdry[10000];
  double x[8],y[8],z[8], eps = 1e-6, sx, sy, sz, r2 = TDatabase::ParamDB->P6*TDatabase::ParamDB->P6;
  TBaseCell *cell;
  TVertex *vertex;
  FE3D CurrentElement;
  TJoint *joint;

  // nodes on the circle should be treated as NEUMANN
  r2+=eps;
  // number of mesh cells
  N_Cells = Coll->GetN_Cells();
  // array with global numbers of d.o.f.
  global_numbers = fespace->GetGlobalNumbers();
  // array which points to the beginning of the global numbers in
  // global_numbers for each mesh cell
  begin_index = fespace->GetBeginIndex();
  // get number of active dof
  N_Active =  fespace->GetActiveBound();
  diri_counter = 0;
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    N_V = cell->GetN_Vertices();
    found = 0;
    for (j=0;j<N_V;j++)
    {
      // read coordinates of the mesh cell
      face_lid[j] = 0;
      vertex = cell->GetVertex(j);
      vertex->GetCoords(x[j], y[j], z[j]);
      // vertex on the upper lid
      //if ((z[j]<1e+0) && (fabs(x[j]*x[j]+y[j]*y[j]+(z[j]-100)*(z[j]-100)-10000)<1e-0))
      if ((z[j]<0.5e+0) && (fabs(x[j]*x[j]+y[j]*y[j]+(z[j]-100)*(z[j]-100)-10000)<10e-0))
      {
        face_lid[j] = 1;
        found++;
	//OutPut(" "<< (fabs(x[j]*x[j]+y[j]*y[j]+(z[j]-100)*(z[j]-100)-10000)));
      }
    }
    // no cell with face with vertex in the outflow region
    if (found < 3)
      continue;
    if (found==4)
      {
	OutPut("4 vertices found"<<endl);
	exit(4711);
      }
    // finite element on the mesh cell
    CurrentElement = fespace->GetFE3D(i, cell);
    // number of basis functions (= number of d.o.f.)
    N_ = TFEDatabase3D::GetN_BaseFunctFromFE3D(CurrentElement);
    // the array which gives the mapping of the local to the global d.o.f.
    dof = global_numbers+begin_index[i];
    switch(CurrentElement)
    {
      // P_2
      case C_P2_3D_T_A:
     // P_2_bubble
      case C_B2_3D_T_A:
        found = 0;
        for (j=0;j<N_V;j++)
        {
          if (face_lid[j])
          {
            face_lid_index[found] = j;
            found++;
          }
        }

        //for (j=0;j<3;j++)
        //	OutPut(x[j] << " " << y[j] << " " << z[j] << endl);
        // loop over all nodal degrees of freedom
        for (j=0;j<10;j++)
        {
          // vertex 3 is not on the top face
          if ((face_lid_index[0]!=3)&& (face_lid_index[1]!=3) && (face_lid_index[2]!=3))
          {
            no_face_lid  = 0; 
	    // check for boundary face
	    joint = cell->GetJoint(0);
	    if (!((joint->GetType() == BoundaryFace ||
		   joint->GetType() == IsoBoundFace)))
	      continue;
           // all local dof with index >5 are not on the top face
            if (j>5)
              continue;
          }
          // vertex 0 is not on the top face
          if ((face_lid_index[0]!=0)&& (face_lid_index[1]!=0) && (face_lid_index[2]!=0))
          {
            no_face_lid  = 2; 
	    // check for boundary face
	    joint = cell->GetJoint(2);
	    if (!((joint->GetType() == BoundaryFace ||
		   joint->GetType() == IsoBoundFace)))
	      continue;
            // the local dofs which are not on the top face
            if ((j==0)||(j==1)||(j==3)||(j==6))
              continue;
          }
          // vertex 1 is not on the top face
          if ((face_lid_index[0]!=1)&& (face_lid_index[1]!=1) && (face_lid_index[2]!=1))
          {
            no_face_lid  = 3; 
	    // check for boundary face
	    joint = cell->GetJoint(3);
	    if (!((joint->GetType() == BoundaryFace ||
		   joint->GetType() == IsoBoundFace)))
	      continue;
            // the local dofs which are not on the top face
            if ((j==1)||(j==2)||(j==4)||(j==7))
              continue;
          }
          // vertex 2 is not on the top face
          if ((face_lid_index[0]!=2)&& (face_lid_index[1]!=2) && (face_lid_index[2]!=2))
          {
            no_face_lid  = 1; 
	    // check for boundary face
	    joint = cell->GetJoint(1);
	    if (!((joint->GetType() == BoundaryFace ||
		   joint->GetType() == IsoBoundFace)))
	      continue;
            // the local dofs which are not on the top face
            if ((j==3)||(j==4)||(j==5)||(j==8))
              continue;
          }
          if (dof[j] < N_Active)
          {
            switch(j)
            {
              case 0:
                sx = x[face_lid_index[0]];
                sy = y[face_lid_index[0]];
		sz = z[face_lid_index[0]];
                break;
              case 1:
                sx = (x[face_lid_index[0]]+x[face_lid_index[1]])/2.0;
                sy = (y[face_lid_index[0]]+y[face_lid_index[1]])/2.0;
                sz = (z[face_lid_index[0]]+z[face_lid_index[1]])/2.0;
                break;
              case 2:
                sx = x[face_lid_index[1]];
                sy = y[face_lid_index[1]];
		sz = z[face_lid_index[1]];
                break;
              case 3:
                sx = (x[face_lid_index[0]]+x[face_lid_index[2]])/2.0;
                sy = (y[face_lid_index[0]]+y[face_lid_index[2]])/2.0;
                sz = (z[face_lid_index[0]]+z[face_lid_index[2]])/2.0;
                break;
              case 4:
                sx = (x[face_lid_index[1]]+x[face_lid_index[2]])/2.0;
                sy = (y[face_lid_index[1]]+y[face_lid_index[2]])/2.0;
                sz = (z[face_lid_index[1]]+z[face_lid_index[2]])/2.0;
                break;
              case 5:
                sx = x[face_lid_index[2]];
                sy = y[face_lid_index[2]];
                sz = z[face_lid_index[2]];
                break;
            }
            //if (sz>6.254e-2+eps)
	    // nodes on the outflow circle will be treated as Dirichlet
	    if (sz>(6.254e-2+eps)*0.876)
	    {
              //OutPut(sz << " " << 6.25195e-2 - sz  << " wrong Neumann node " << dof[j] << endl);
              tmp_diri[diri_counter] = dof[j];
              diri_counter++;
              // set face bubble to Dirichlet
	      /* if (CurrentElement==C_B2_3D_T_A)
             {
                 switch(no_face_lid)
                {
               case 0:
                   tmp_diri[diri_counter] = dof[10];
                   diri_counter++;                          
                   break;
               case 1:
                   tmp_diri[diri_counter] = dof[11];
                   diri_counter++;                          
                   break;
               case 2:
                   tmp_diri[diri_counter] = dof[12];
                   diri_counter++;                          
                   break;
               case 3:
                   tmp_diri[diri_counter] = dof[13];
                   diri_counter++;                          
                   break;
               }
		 //OutPut(sz << " " << 6.254e-2 - sz  << " wrong Neumann node on face " << dof[j] << endl);
		 }*/
            }
	    else
	      {
		OutPut("correct Neumann node " << " " << sz  << " " <<  dof[j] << endl);
		//if (sz<=6.254e-2+eps)
		// {
		    tmp_neum_bdry[neum_bdry_counter] = dof[j];
		    neum_bdry_counter++;
		    // }
	      }
            if (diri_counter > 10000)
            {
              OutPut("tmp_diri too short !!!"<<endl);
              exit(4711);
            }
          }
          else
          {
            switch(j)
            {
              case 0:
                sx = x[face_lid_index[0]];
                sy = y[face_lid_index[0]];
                sz = z[face_lid_index[0]];
                break;
              case 1:
                sx = (x[face_lid_index[0]]+x[face_lid_index[1]])/2.0;
                sy = (y[face_lid_index[0]]+y[face_lid_index[1]])/2.0;
                sz = (z[face_lid_index[0]]+z[face_lid_index[1]])/2.0;
                break;
              case 2:
                sx = x[face_lid_index[1]];
                sy = y[face_lid_index[1]];
                sz = z[face_lid_index[1]];
                break;
              case 3:
                sx = (x[face_lid_index[0]]+x[face_lid_index[2]])/2.0;
                sy = (y[face_lid_index[0]]+y[face_lid_index[2]])/2.0;
                sz = (z[face_lid_index[0]]+z[face_lid_index[2]])/2.0;
                break;
              case 4:
                sx = (x[face_lid_index[1]]+x[face_lid_index[2]])/2.0;
                sy = (y[face_lid_index[1]]+y[face_lid_index[2]])/2.0;
                sz = (z[face_lid_index[1]]+z[face_lid_index[2]])/2.0;
                break;
              case 5:
                sx = x[face_lid_index[2]];
                sy = y[face_lid_index[2]];
                sz = z[face_lid_index[2]];
                break;
            }
            if (sz<6.254e-2+eps)
            {
              OutPut(sz << " wrong Dirichlet node " << dof[j] << endl);
              OutPut("This case cannot be handled !!!" << endl);
              exit(4711);
            }
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
    tmp_diri[found] = -1;
  }


  // condense
  for (i=0;i<neum_bdry_counter;i++)
  {
    if (tmp_neum_bdry[i] == -1)
      continue;
    neum_bdry_counter_1++;
    for (j=i+1;j<neum_bdry_counter;j++)
    {
      if (tmp_neum_bdry[i] == tmp_neum_bdry[j])
      {
        tmp_neum_bdry[j] = -1;
      }
    }
  }

  OutPut("CheckNeumannNodesForVelocity: N_neum on outflow boundary " << neum_bdry_counter_1 << endl);

  N_neum_bdry = neum_bdry_counter_1;
  // allocate array for the indices
  neum_bdry = new int[neum_bdry_counter_1];
  // fill array and sort
  for (i=0;i<neum_bdry_counter_1;i++)
  {
    min_val = tmp_neum_bdry[0];
    found = 0;
    for (j=1;j<neum_bdry_counter;j++)
    {
      if ((tmp_neum_bdry[j]>0) && ((tmp_neum_bdry[j] < min_val) ||
        (min_val == -1)))
      {
        min_val =  tmp_neum_bdry[j];
        found = j;
      }
    }
    neum_bdry[i] = tmp_neum_bdry[found];
    //OutPut(i << " " <<  neum_bdry[i] << endl);
    tmp_neum_bdry[found] = -1;
  }
}


// ========================================================================
//
// manipulate matrices and rhs such that all Dirichlet nodes are fix
//
// ========================================================================
void SetDirichletNodesFromNeumannNodes(TSquareMatrix3D **SQMATRICES,
TMatrix3D **MATRICES,
double *rhs,
int N_U,
int N_neum_to_diri,
int *neum_to_diri)
{
  TSquareMatrix3D *MatrixA11;
  TMatrix3D *Matrix_B;
  double *Entries_A;
  int i, j, l, l0, l1, index, *RowPtr_A, *KCol_A;

  switch(TDatabase::ParamDB->NSTYPE)
  {
    case 2:
      MatrixA11 = SQMATRICES[0];
      RowPtr_A      = MatrixA11->GetRowPtr();
      KCol_A        = MatrixA11->GetKCol();
      Entries_A     = MatrixA11->GetEntries();
      for (i=0;i<N_neum_to_diri;i++)
      {
        index = neum_to_diri[i];
        l0 = RowPtr_A[index];
        l1 = RowPtr_A[index+1];
        for (l=l0;l<l1;l++)
        {
          // diagonal entry
          if (KCol_A[l]==index)
            Entries_A[l] = 1;
          else
            Entries_A[l] = 0;
        }
      }
      break;
    case 4:
      MatrixA11 = SQMATRICES[0];
      RowPtr_A      = MatrixA11->GetRowPtr();
      KCol_A        = MatrixA11->GetKCol();
      Entries_A     = MatrixA11->GetEntries();
      for (i=0;i<N_neum_to_diri;i++)
      {
        index = neum_to_diri[i];
        l0 = RowPtr_A[index];
        l1 = RowPtr_A[index+1];
        for (l=l0;l<l1;l++)
        {
          // diagonal entry
          if (KCol_A[l]==index)
            Entries_A[l] = 1;
          else
            Entries_A[l] = 0;
        }
      }
      // A12
      MatrixA11 = SQMATRICES[1];
      RowPtr_A      = MatrixA11->GetRowPtr();
      KCol_A        = MatrixA11->GetKCol();
      Entries_A     = MatrixA11->GetEntries();
      for (i=0;i<N_neum_to_diri;i++)
      {
        index = neum_to_diri[i];
        l0 = RowPtr_A[index];
        l1 = RowPtr_A[index+1];
        for (l=l0;l<l1;l++)
        {
          Entries_A[l] = 0;
        }
      }
      // A13
      MatrixA11 = SQMATRICES[2];
      RowPtr_A      = MatrixA11->GetRowPtr();
      KCol_A        = MatrixA11->GetKCol();
      Entries_A     = MatrixA11->GetEntries();
      for (i=0;i<N_neum_to_diri;i++)
      {
        index = neum_to_diri[i];
        l0 = RowPtr_A[index];
        l1 = RowPtr_A[index+1];
        for (l=l0;l<l1;l++)
        {
          Entries_A[l] = 0;
        }
      }
      // A21
      MatrixA11 = SQMATRICES[3];
      RowPtr_A      = MatrixA11->GetRowPtr();
      KCol_A        = MatrixA11->GetKCol();
      Entries_A     = MatrixA11->GetEntries();
      for (i=0;i<N_neum_to_diri;i++)
      {
        index = neum_to_diri[i];
        l0 = RowPtr_A[index];
        l1 = RowPtr_A[index+1];
        for (l=l0;l<l1;l++)
        {
          Entries_A[l] = 0;
        }
      }
      // A22
      MatrixA11 = SQMATRICES[4];
      RowPtr_A      = MatrixA11->GetRowPtr();
      KCol_A        = MatrixA11->GetKCol();
      Entries_A     = MatrixA11->GetEntries();
      for (i=0;i<N_neum_to_diri;i++)
      {
        index = neum_to_diri[i];
        l0 = RowPtr_A[index];
        l1 = RowPtr_A[index+1];
        for (l=l0;l<l1;l++)
        {
          // diagonal entry
          if (KCol_A[l]==index)
            Entries_A[l] = 1;
          else
            Entries_A[l] = 0;
        }
      }
      // A23
      MatrixA11 = SQMATRICES[5];
      RowPtr_A      = MatrixA11->GetRowPtr();
      KCol_A        = MatrixA11->GetKCol();
      Entries_A     = MatrixA11->GetEntries();
      for (i=0;i<N_neum_to_diri;i++)
      {
        index = neum_to_diri[i];
        l0 = RowPtr_A[index];
        l1 = RowPtr_A[index+1];
        for (l=l0;l<l1;l++)
        {
          Entries_A[l] = 0;
        }
      }
      // A31
      MatrixA11 = SQMATRICES[6];
      RowPtr_A      = MatrixA11->GetRowPtr();
      KCol_A        = MatrixA11->GetKCol();
      Entries_A     = MatrixA11->GetEntries();
      for (i=0;i<N_neum_to_diri;i++)
      {
        index = neum_to_diri[i];
        l0 = RowPtr_A[index];
        l1 = RowPtr_A[index+1];
        for (l=l0;l<l1;l++)
        {
          Entries_A[l] = 0;
        }
      }
      // A32
      MatrixA11 = SQMATRICES[7];
      RowPtr_A      = MatrixA11->GetRowPtr();
      KCol_A        = MatrixA11->GetKCol();
      Entries_A     = MatrixA11->GetEntries();
      for (i=0;i<N_neum_to_diri;i++)
      {
        index = neum_to_diri[i];
        l0 = RowPtr_A[index];
        l1 = RowPtr_A[index+1];
        for (l=l0;l<l1;l++)
        {
          Entries_A[l] = 0;
        }
      }
      // A33
      MatrixA11 = SQMATRICES[8];
      RowPtr_A      = MatrixA11->GetRowPtr();
      KCol_A        = MatrixA11->GetKCol();
      Entries_A     = MatrixA11->GetEntries();
      for (i=0;i<N_neum_to_diri;i++)
      {
        index = neum_to_diri[i];
        l0 = RowPtr_A[index];
        l1 = RowPtr_A[index+1];
        for (l=l0;l<l1;l++)
        {
          // diagonal entry
          if (KCol_A[l]==index)
            Entries_A[l] = 1;
          else
            Entries_A[l] = 0;
        }
      }
      break;
    case 1:
    case 3:
      OutPut("SetDirichletNodesFromNeumannNodes does not work for NSTYPE: "<<
        TDatabase::ParamDB->NSTYPE << endl);
      exit(4711);
    default:
      OutPut("SetDirichletNodesFromNeumannNodes not implemented"<<endl);
      exit(4711);
  }

  for (j=0;j<3;j++)
  {
    Matrix_B = MATRICES[j];
    RowPtr_A      = Matrix_B->GetRowPtr();
    Entries_A     = Matrix_B->GetEntries();
    for (i=0;i<N_neum_to_diri;i++)
    {
      index = neum_to_diri[i];
      l0 = RowPtr_A[index];
      l1 = RowPtr_A[index+1];
      for (l=l0;l<l1;l++)
        Entries_A[l] = 0;
    }
  }
  for (j=0;j<3;j++)
  {
    l0 = j* N_U;
    for (i=0;i<N_neum_to_diri;i++)
    {
      index = neum_to_diri[i] + l0;
      rhs[index] = 0;
    }
  }
}

// ========================================================================
//
// move grid points to obtain circular inlet and outlet
//
// ========================================================================

void MakeInAndOutflowCircular(TCollection *Coll)
{
  int i,k, j, N_Cells, N_Active, N_V, found, N_;
  int  face_lid[8];
  int face_lid_index[4];
  double x[8],y[8],z[8],eps = 1e-6;
  TBaseCell *cell;
  TVertex *vertex;
  double d,x1,y1,z1,dist_max,d1,xdir,ydir,zdir,dneu;
  int index, index1,index_max;
  double  r=2.5; 
  int found1;
 
  OutPut("MakeInAndOutflowCircular"<<endl);

  /********************************************/
  /**       Round inlet                      **/ 
  /********************************************/
  // number of mesh cells
  N_Cells = Coll->GetN_Cells();
  // loop over the mesh cells
  for(i=0;i<N_Cells;i++)
  {
      cell = Coll->GetCell(i);
      cell->SetClipBoard(0);
      //number of vertices per element
      N_V = cell->GetN_Vertices();
      //counters for 3 vertices, that forms a face, of the top
      found = found1 = 0;
      // loop over the number of vetices per element
      for (j=0;j<N_V;j++)
      {
	  // read coordinates of the mesh cell
	  face_lid[j] = 0;
	  vertex = cell->GetVertex(j);
	  vertex->GetCoords(x[j], y[j], z[j]);
	  //vertices of elements on the top and interior of the left inlet
	  if((fabs(z[j]-140)<eps) && ((x[j]+25)*(x[j]+25)+y[j]*y[j]<6.25+eps) )
	  {
	      //count the vertices with this condition
	      face_lid_index[found] = j;
	      k=j;
	      face_lid[j]=1;
	      found++;
	      //OutPut(" bun   "<< x[j] <<"    "<<y[j]<< "    "<<z[j] <<endl );
	  }
	  //analog for the right outlet
	  if((fabs(z[j]-140)<eps) && ((x[j]-25)*(x[j]-25)+y[j]*y[j]<6.25+eps) )
	  {
	      face_lid_index[found1] = j;
	      k=j;
	      face_lid[j]=1;
	      found1++;
	      OutPut(i << " " <<j <<" nou  "<< x[j] <<"    "<<y[j]<< "    "<<z[j] <<endl );
	  }
	  
      }
      // for the left inlet observe only the three consecutiv founded vertices (a face = a triangle)
      if(found==3)
      {
	  // OutPut("i " << i <<endl );
	  index=0;
	  dist_max = -1;
	  for (k=0;k<3;k++)
	  {
	      //save the 3 coordinates of a founded faces and calculate the distance between every face vertices and the middle point of the inlet 
	      x1=x[face_lid_index[k]];
	      y1=y[face_lid_index[k]];
	      z1=z[face_lid_index[k]];
	      OutPut("    "<< x1 );  
	      OutPut("    "<< y1 );
	      OutPut("    "<< z1 );
	      d=sqrt((x1+25.0)*(x1+25.0)+(y1)*(y1)+(z1-140.0)*(z1-140.0));
	      OutPut("    "<< d <<endl );  
	      //vertices on the inlet circle are counted
	      if (fabs(d-2.5)<eps)
	      {
		  index++;
		  continue;
	      }
	      //for the 2 face vertices at the same distance of the inlet, observe only one
	      //observe a face with exactly one point of the inlet circle  and the other 2 interior of the inlet  but not at the same distance of the inlet center
	      if (fabs(d-dist_max) < eps)
	      {
		  index = -2;
	      }
	      else
	      {
		  if (d > dist_max)
		  {
		      dist_max = d;
		      index_max = face_lid_index[k];
		  }
                  }
	  } 
	  if (index==1)
	  {
	      //translate the vertex with this condition on the inlet circle 
	      cell->SetClipBoard(1);
	  }               
      }
      //analog for the right outlet
      if(found1==3)
      {
	  index=0;
	  dist_max = -1;
	  for (k=0;k<3;k++)
	  {
	      
	      x1=x[face_lid_index[k]];
	      y1=y[face_lid_index[k]];
	      z1=z[face_lid_index[k]];
	      d=sqrt((x1-25.0)*(x1-25.0)+(y1)*(y1)+(z1-140.0)*(z1-140.0));
	      if (fabs(d-2.5)<eps)
	      {
		  index++;
		  continue;
	      }
	      if (fabs(d-dist_max) < eps)
	      {
		      index = -2;
	      }
	      else
	      {
		  if (d > dist_max)
		  {
			  dist_max = d;
			  index_max = face_lid_index[k];
		  }
	      }
	  } 
	  if (index==1)
	  {
	      cell->SetClipBoard(1);
	  }               
      }
  }
  /********************************************/
  /**Round second time the inlet, outlet     **/ 
  /********************************************/
  for(i=0;i<N_Cells;i++)
      //for(i=0;i<2;i++)
  {
      cell = Coll->GetCell(i);
      if (cell->GetClipBoard() != 1)
	  continue;
      N_V = cell->GetN_Vertices();
      // OutPut("N_V="  << N_V << endl);
      //OutPut("i"  << i << endl);
      found = found1 = 0;
      for (j=0;j<N_V;j++)
      {
	  // read coordinates of the mesh cell
	  face_lid[j] = 0;
	  vertex = cell->GetVertex(j);
	  vertex->GetCoords(x[j], y[j], z[j]);
	  if((fabs(z[j]-140)<eps) && ((x[j]+25)*(x[j]+25)+y[j]*y[j]<6.25+eps) )
	  {
	      face_lid_index[found] = j;
	      k=j;
	      face_lid[j]=1;
	      found++;
	      //OutPut(" bun   "<< x[j] <<"    "<<y[j]<< "    "<<z[j] <<endl );
	  }
	if((fabs(z[j]-140)<eps) && ((x[j]-25)*(x[j]-25)+y[j]*y[j]<6.25+eps) )
	{
	    face_lid_index[found1] = j;
	    k=j;
	    face_lid[j]=1;
	    found1++;
	    OutPut(i << " " <<j <<" nou  "<< x[j] <<"    "<<y[j]<< "    "<<z[j] <<endl );
	}
      }
      
      if(found==3)
      {
	  OutPut("i " << i <<endl );
	  index=0;
	  dist_max = -1;
	  for (k=0;k<3;k++)
	  {
	      
	      x1=x[face_lid_index[k]];
	      y1=y[face_lid_index[k]];
	      z1=z[face_lid_index[k]];
	      OutPut("    "<< x1 );  
	      OutPut("    "<< y1 );
	      OutPut("    "<< z1 );
	      d=sqrt((x1+25.0)*(x1+25.0)+(y1)*(y1)+(z1-140.0)*(z1-140.0));
	      OutPut("    "<< d <<endl );  
	      if (fabs(d-2.5)<eps)
	      {
		  index++;
		  continue;
	      }
	      if (d > dist_max)
	      {
		  dist_max = d;
		  index_max = face_lid_index[k];
	      }
	  } 
	  if ((index==1))//&&(dist_max >1.75))
	  {
	      x1=x[index_max];
	      y1=y[index_max];
	      z1=z[index_max];
	      d1=sqrt((x1+25.0)*(x1+25.0)+(y1)*(y1)+(z1-140.0)*(z1-140.0));
	      OutPut("    "<< x1 );  
	      OutPut("    "<< y1 );
	      OutPut("    "<< z1 );
	      OutPut("    "<< d1);
	      xdir=-25.+(x1+25.)*r/d1;
	      ydir=y1*r/d1;
	      zdir=140.+(z1-140.)*r/d1;
	      
	      OutPut("    "<< xdir );  
	      OutPut("    "<< ydir );
	      OutPut("    "<< zdir << endl );
	      vertex = cell->GetVertex(index_max);
	      vertex->SetCoords(xdir,ydir,zdir);
	      //OutPut("to move "  <<index_max << endl); 
	  }               
      }
      
      if(found1==3)
      {
	  index=0;
	  dist_max = -1;
	  for (k=0;k<3;k++)
	  {
	      x1=x[face_lid_index[k]];
	      y1=y[face_lid_index[k]];
	      z1=z[face_lid_index[k]];
	      //OutPut(face_lid_index[k] << "  "<< x1 );  
	      //OutPut("    "<< y1 );
	      //OutPut("    "<< z1 );
	      d=sqrt((x1-25.0)*(x1-25.0)+(y1)*(y1)+(z1-140.0)*(z1-140.0));
	      //OutPut("    "<< d <<endl );  
	      if (fabs(d-2.5)<eps)
	      {
		  index++;
		  continue;
	      }
	      if (d > dist_max)
	      {
		  dist_max = d;
		  index_max = face_lid_index[k];
	      }
	  } 
	  if ((index==1))//&&(dist_max >1.75))
	  {
	      x1=x[index_max];
	      y1=y[index_max];
	      z1=z[index_max];
	      d1=sqrt((x1-25.0)*(x1-25.0)+(y1)*(y1)+(z1-140.0)*(z1-140.0));
	      //OutPut("vor    "<< x1 );  
	      //OutPut("    "<< y1 );
	      // OutPut("    "<< z1 );
	      //OutPut("    "<< d1);
	      xdir=25.+(x1-25.)*r/d1;
	      ydir=y1*r/d1;
	      zdir=140.+(z1-140.)*r/d1;
	      
	      //OutPut("hint    "<< xdir );  
	      //OutPut("    "<< ydir );
	      //OutPut("    "<< zdir << endl );
	      vertex = cell->GetVertex(index_max);
	      vertex->SetCoords(xdir,ydir,zdir);
	      OutPut("to move "  <<index_max << endl); 
	  }               
      }
  }
}


