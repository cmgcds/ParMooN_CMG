// ======================================================================
// Example from P.W. Hemker 1996, extended to 3D
// ======================================================================
#define __HEMKER1996__

#include <ConvDiff3D.h>
#include <PeriodicJoint.h>
void ExampleFile()
{
  OutPut("Example: Hemker1996_3D.h ") ;
  TDatabase::ParamDB->DRIFT_Z = 3.0;
  OutPut("DRIFT_Z set to  " << TDatabase::ParamDB->DRIFT_Z  << endl) ;  
}

// exact solution
void Exact(double x, double y,  double z, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}

void BoundCondition(double x, double y, double z, BoundCond &cond)
{
    double eps = 1e-6;

  if (TDatabase::ParamDB->SOLD_PARAMETER_TYPE != FEM_TVD)
  {
      if ((fabs(y-3)<eps) || (fabs(y+3)<eps) || (fabs(x-9)<eps))
	  cond = NEUMANN;
      else
	  cond = DIRICHLET;
      }
  else
      cond = NEUMANN;
}


// value of boundary condition
void BoundValue(double x, double y, double z, double &value)
{
    double eps = 1e-6;

    if ((fabs(y-3)<eps) || (fabs(y+3)<eps) || (fabs(x-9)<eps) || (fabs(x+3)<eps) )
	value = 0;
    else
	value = 1;
}

void BilinearCoeffs(int n_points, double *x, double *y, double *z,
        double **parameters, double **coeffs)
{
  double eps=1/TDatabase::ParamDB->RE_NR;
  double angle = 0, v1, v2;
  int i;
  double *coeff, *param;

  v1 = cos(angle);
  v2 = sin(angle);

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    coeff[0] = eps;
    coeff[1] = v1;
    coeff[2] = v2;
    coeff[3] = 0;
    coeff[4] = 0;
    coeff[5] = 0;

    coeff[6] = 0;
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
 
// periodic b.c. :
//  z = 0 and z = 3
void SetPeriodicFaceJoints(TCollection *Coll)
{
  int i, j, N_Cells, N_Faces, l1, jj, ii, N_Faces1;
  double X[MaxN_QuadPoints_3D], Y[MaxN_QuadPoints_3D];
  double Z[MaxN_QuadPoints_3D];
  double X1[MaxN_QuadPoints_3D], Y1[MaxN_QuadPoints_3D];
  double Z1[MaxN_QuadPoints_3D];
  double x,y,z,x1,y1,z1;
  TBaseCell *cell, *cell1;
  TJoint *joint, *joint1;
  const int *TmpFV, *TmpLen;
  int MaxLen, MaxLen1, z_vert_0, z_vert_2;
  int found, z1_vert_0, z1_vert_2;
  double RE=TDatabase::ParamDB->RE_NR;
  double z_bottom, z_top;

  z_bottom = 0;
  z_top = TDatabase::ParamDB->DRIFT_Z; 

  N_Cells = Coll->GetN_Cells();

  OutPut("SetPeriodicJoints " << endl);
  // first loop over cells
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    N_Faces = cell->GetN_Faces();
    // check if face on boundary
    for (j=0;j<N_Faces;j++)
    {
      joint=cell->GetJoint(j);
      //not on boundary
      if ((!(joint->GetType() == BoundaryFace)) &&
        (!(joint->GetType() == IsoBoundFace)))
        continue;
      // find vertices on the face
      cell->GetShapeDesc()->GetFaceVertex(TmpFV, TmpLen, MaxLen);
      //OutPut("vertices on face "<< TmpLen[j]<< endl);
      z_vert_0 =  z_vert_2 = 0;
      // compute coordinates of vertices
      for (l1=0;l1<TmpLen[j];l1++)
      {
        cell->GetVertex(TmpFV[j*MaxLen+l1])->GetCoords(X[l1], Y[l1], Z[l1]);
        // check if vertex on z = z_bottom
        if (fabs(Z[l1]-z_bottom)<1e-5)
          z_vert_0++;
        // check if vertex on z = z_top
        if (fabs(Z[l1]-z_top)<1e-5)
          z_vert_2++;
      }
      found = 0;
      if (z_vert_0==TmpLen[j])
        found++;
      if (z_vert_2==TmpLen[j])
        found++;
      // face not at the interesting boundaries
      if (!found)
        continue;
      // compute barycentre
      x = y = z = 0;
      for (l1=0;l1<TmpLen[j];l1++)
      {
        x+=X[l1];
        y+=Y[l1];
        z+=Z[l1];
      }
      x /= TmpLen[j];
      y /= TmpLen[j];
      z /= TmpLen[j];
      //OutPut("bary " << x << " " << y << " " << z << endl);
      //for (l1=0;l1<TmpLen[j];l1++)
      // {
      // OutPut("face " << X[l1] << " " << Y[l1] << " " << Z[l1] << endl);
      //}
      // inner loop over the cells
      for(ii=i+1;ii<N_Cells;ii++)
      {
        cell1 = Coll->GetCell(ii);
        N_Faces1 = cell1->GetN_Faces();
        // check if face on boundary
        for (jj=0;jj<N_Faces1;jj++)
        {
          joint1=cell1->GetJoint(jj);
          //not on boundary
          if (!(joint1->GetType() == BoundaryFace))
            continue;
          // find vertices on the face
          cell1->GetShapeDesc()->GetFaceVertex(TmpFV, TmpLen, MaxLen1);

          x1 = y1 = z1 = 0;
          z1_vert_0 =  z1_vert_2 = 0;
          // compute coordinates of vertices
          // count number of vertices which fulfill the criteria
          for (l1=0;l1<TmpLen[jj];l1++)
          {
            cell1->GetVertex(TmpFV[jj*MaxLen1+l1])->
              GetCoords(X1[l1], Y1[l1], Z1[l1]);
            // check if vertex on z = z_bottom
            if (fabs(Z1[l1]-z_bottom)<1e-5)
              z1_vert_0++;
            // check if vertex on z = z_top
            if (fabs(Z1[l1]-z_top)<1e-5)
              z1_vert_2++;
            x1+=X1[l1];
            y1+=Y1[l1];
            z1+=Z1[l1];
          }

          // bary center
          x1 /= TmpLen[jj];
          y1 /= TmpLen[jj];
          z1 /= TmpLen[jj];
          found = 0;
          //  OutPut("baryopp " << x1 << " " << y1 << " " << z1 << endl);

          // all vertices of original face are on z_bottom
          // and all vertices of second face are on z_top
          if ((z_vert_0==TmpLen[j])&&(z1_vert_2==TmpLen[jj]))
          {
            // check if the x,y - coordinates of the barycenters are the same
            if  ((fabs(x-x1)>1e-5)||(fabs(y-y1)>1e-5))
              continue;
            // the face match
            found++;
          }
          // do the same vice versa
          if ((z_vert_2==TmpLen[j])&&(z1_vert_0==TmpLen[jj]))
          {
            if  ((fabs(x-x1)>1e-5)||(fabs(y-y1)>1e-5))
              continue;
            found++;
          }
          if (!found)
            continue;
          // for (l1=0;l1<TmpLen[jj];l1++)
          // {
          //  OutPut("opp " << X1[l1] << " " << Y1[l1] << " " << Z1[l1] << endl);
          // }
          // delete old joints
          delete joint;
          delete joint1;
          // make new joint
          joint = new TPeriodicJoint(cell, cell1);
          // set joint
          cell->SetJoint(j,joint);
          cell1->SetJoint(jj,joint);
          // find opposite vertex to local vertex zero of face
          if (((z_vert_0==TmpLen[j])&&(z1_vert_2==TmpLen[jj]))
            || ((z_vert_2==TmpLen[j])&&(z1_vert_0==TmpLen[jj])))
          {
            found = -1;
            for (l1=0;l1<TmpLen[jj];l1++)
            {
              if ((fabs(X[0]-X1[l1])<1e-5)&& (fabs(Y[0]-Y1[l1])<1e-5))
              {
                found = l1;
                break;
              }
            }
          }
          //OutPut("opposite to zero vertex " << found << endl);
          joint->SetMapType(found);
        }
      }
    }
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
   const int max_entries = 100000;  
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
      OutPut(i<<endl);
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
      if ((fabs(x[j]+3)<eps)||(sqrt(x[j]*x[j]+y[j]*y[j])<=1+eps))
      {
              // Dirichlet boundary
	   boundary_vertices[j] = 1;
	   found++;
      }
    }
    OutPut("f" << found<<endl);
    // no cell with face with vertex on the boundary
    if (found<3) 
	continue;
    // finite element on the mesh cell
    CurrentElement = fespace->GetFE3D(i, cell);
    // number of basis functions (= number of d.o.f.)
    N_ = TFEDatabase3D::GetN_BaseFunctFromFE3D(CurrentElement);
    // the array which gives the mapping of the local to the global d.o.f.
    dof = global_numbers+begin_index[i];
    OutPut("ff" << dof[0]<<endl);
    switch(CurrentElement)
    {
	// P_1, Q_1
	case C_P1_3D_T_A:
	case C_Q1_3D_H_A:
	case C_Q1_3D_H_M:
	    for (j=0;j<N_V;j++)
	    {
		OutPut("j" << j<<endl);

		// vertex on the boundary
		if (boundary_vertices[j])
		{
		    if (CurrentElement==C_P1_3D_T_A)
			tmp_diri[diri_counter] = dof[j];
		    else
		    {
			if ((j==0)||(j==1)||(j==4)||(j==5))
                         {
		OutPut("j" << x[j] << " " << y[j] << " " << z[j]<< endl);
			     tmp_diri[diri_counter] = dof[j];
		OutPut("j" << j<<endl);
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
		    OutPut("d"<< diri_counter<<endl);
                    tmp_x[diri_counter] = x[j];
                    tmp_y[diri_counter] = y[j];
                    tmp_z[diri_counter] = z[j];
		    OutPut("e"<< diri_counter<<endl);
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
    OutPut("c"<<endl);
 
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
