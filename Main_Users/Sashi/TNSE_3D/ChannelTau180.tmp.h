// Navier-Stokes problem
// channel flow in 3D
//

#include <PeriodicJoint.h>

#define U_INFTY 1
#define __CHANNEL_TOBIAS__
#define __CHANNEL_TAU180__

// ========================================================================
// example file
// ========================================================================

void ExampleFile()
{
  OutPut("Example: ChannelTau180.h"<<endl);
  TDatabase::ParamDB->INTERNAL_PROBLEM_IDENTITY = 180;
}


// ========================================================================
// exact solution
// ========================================================================
double DNS_profile(double zz)
{
  int i;
  double value, val;

  double z[65] =  {0.0000e-00, 3.0118e-04, 1.2045e-03, 2.7095e-03, 4.8153e-03, 7.5205e-03, 1.0823e-02, 1.4722e-02, 1.9215e-02, 2.4298e-02, 2.9969e-02, 3.6224e-02,  
	           4.3060e-02, 5.0472e-02, 5.8456e-02, 6.7007e-02, 7.6120e-02, 8.5790e-02, 9.6011e-02, 1.0678e-01, 1.1808e-01, 1.2991e-01, 1.4227e-01, 1.5515e-01,  
	           1.6853e-01, 1.8242e-01, 1.9679e-01, 2.1165e-01, 2.2699e-01, 2.4279e-01, 2.5905e-01, 2.7575e-01, 2.9289e-01, 3.1046e-01, 3.2844e-01, 3.4683e-01,  
	           3.6561e-01, 3.8477e-01, 4.0430e-01, 4.2419e-01, 4.4443e-01, 4.6500e-01, 4.8590e-01, 5.0710e-01, 5.2860e-01, 5.5039e-01, 5.7244e-01, 5.9476e-01,  
	           6.1732e-01, 6.4010e-01, 6.6311e-01, 6.8632e-01, 7.0972e-01, 7.3329e-01, 7.5702e-01, 7.8090e-01, 8.0491e-01, 8.2904e-01, 8.5327e-01, 8.7759e-01,  
	           9.0198e-01, 9.2644e-01, 9.5093e-01, 9.7546e-01, 1.0000e-00};  

  double  Umean[65] = {0.0000e+00, 5.3639e-02, 2.1443e-01, 4.8197e-01, 8.5555e-01, 1.3339e+00, 1.9148e+00, 2.5939e+00, 3.3632e+00, 4.2095e+00, 5.1133e+00,
         	                     6.0493e+00, 6.9892e+00, 7.9052e+00, 8.7741e+00, 9.5790e+00, 1.0311e+01, 1.0967e+01, 1.1550e+01, 1.2066e+01, 1.2520e+01, 1.2921e+01,  
		   1.3276e+01, 1.3590e+01, 1.3870e+01, 1.4121e+01, 1.4349e+01, 1.4557e+01, 1.4750e+01, 1.4931e+01, 1.5101e+01, 1.5264e+01, 1.5419e+01,  
		   1.5569e+01, 1.5714e+01, 1.5855e+01, 1.5993e+01, 1.6128e+01, 1.6260e+01, 1.6389e+01, 1.6515e+01, 1.6637e+01, 1.6756e+01, 1.6872e+01,  
	                     1.6985e+01, 1.7094e+01, 1.7200e+01, 1.7302e+01, 1.7400e+01, 1.7494e+01, 1.7585e+01, 1.7672e+01, 1.7756e+01, 1.7835e+01, 1.7911e+01,  
      	                     1.7981e+01, 1.8045e+01, 1.8103e+01, 1.8154e+01, 1.8198e+01, 1.8235e+01, 1.8264e+01, 1.8285e+01, 1.8297e+01, 1.8301e+01}; 

  if (zz<=1)
    {
      val = zz;
    }
  else
    {
      val = 2-zz;
    }

  for(i=0;i<64;i++)
    {
      if(val>=z[i] && val<=z[i+1])
	{
	  value = val*(Umean[i+1]-Umean[i])/(z[i+1]-z[i]) + (Umean[i]*z[i+1]-Umean[i+1]*z[i])/(z[i+1]-z[i]); 
	}
    }

  return(value);
}

void InitialU1(double x, double y, double z, double *values)
{
  // in contrast to the literature, the coordinates y and z are interchanged!!
  // the initial setup is as in V. Gravemeier, J. Comput. Phys. (2006)
  // with 10% random noise (positive and negative)
  double  noise = 0.1;
  double scale, RE=TDatabase::ParamDB->RE_NR;
  
  if (RE==180)
  {
      scale = 25*0.9378;
      TDatabase::ParamDB->INTERNAL_BULK_BEFORE = 15.6803;
      TDatabase::ParamDB->INTERNAL_BULK_AFTER = 15.6803;
  }
  else
  {
    if (RE==590)
	scale = 28; 
    else
    {
	// arithmetic mean
	scale = 3.0/410.0 * (RE-180) + 25;  
    }
  }

  // values[0] = scale*z*(2-z)+noise*2*scale*(2*(double)rand()/RAND_MAX-1)/3;
  values[0] = DNS_profile(z)+noise*2*scale*(2*(double)rand()/RAND_MAX-1)/3;
}


void InitialU2(double x, double y, double z, double *values)
{
  double  noise = 0.1, scale;
  double RE=TDatabase::ParamDB->RE_NR;
  if (RE==180)
  {
      scale = 25*0.9378;
  }
  else
  {
    if (RE==590)
	scale = 28; 
    else
    {
	// arithmetic mean
	scale = 3.0/410.0 * (RE-180) + 25;  
    }
  }

  values[0] = noise*2*scale*(2*(double)rand()/RAND_MAX-1)/3;
}


void InitialU3(double x, double y, double z, double *values)
{
  double  noise = 0.1, scale;
  double RE=TDatabase::ParamDB->RE_NR;
  if (RE==180)
  {
      scale = 25*0.9378;
  }
  else
  {
    if (RE==590)
	scale = 28; 
    else
    {
	// arithmetic mean
	scale = 3.0/410.0 * (RE-180) + 25;  
    }
  }
  values[0] = noise*2*scale*(2*(double)rand()/RAND_MAX-1)/3;
}


void InitialP(double x, double y,  double z, double *values)
{
  values[0] = 0;
  TDatabase::ParamDB->INTERNAL_PERIODIC_IDENTITY = 4;
}


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
  cond = DIRICHLET;
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
  value = 0;
}


// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *X, double *Y, double *Z,
double **parameters, double **coeffs)
{
  double eps = 1/TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff, u1, u2, dt;

  u1 = TDatabase::ParamDB->INTERNAL_BULK_BEFORE;
  u2 = TDatabase::ParamDB->INTERNAL_BULK_AFTER;
  dt = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    coeff[0] = eps;
    coeff[1] = 1 + (u1 - u2)/ dt ;
    coeff[2] = 0;
    coeff[3] = 0;
  }
}


// periodic b.c. :
//  x = -1 and x = 1
//  z = 0 and z = 2
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
  int MaxLen, MaxLen1, x_vert_m1, x_vert_1, y_vert_0, y_vert_2;
  int found, x1_vert_m1, x1_vert_1, y1_vert_0, y1_vert_2;
  double y_bottom = -2*Pi/3, y_top = 2*Pi/3, x_bottom = -2*Pi, x_top = 2*Pi;

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
      y_vert_0 =  y_vert_2 = x_vert_m1 = x_vert_1 = 0;
      // compute coordinates of vertices
      for (l1=0;l1<TmpLen[j];l1++)
      {
        cell->GetVertex(TmpFV[j*MaxLen+l1])->GetCoords(X[l1], Y[l1], Z[l1]);
        // check if vertex on y = y_bottom
        if (fabs(Y[l1]-y_bottom)<1e-5)
          y_vert_0++;
        // check if vertex on y = y_top
        if (fabs(Y[l1]-y_top)<1e-5)
          y_vert_2++;
        // check if vertex on x = x_bottom
        if (fabs(X[l1]-x_bottom)<1e-5)
          x_vert_m1++;
        // check if vertex on x = x_top
        if (fabs(X[l1]-x_top)<1e-5)
          x_vert_1++;
      }
      found = 0;
      if (y_vert_0==TmpLen[j])
        found++;
      if (y_vert_2==TmpLen[j])
        found++;
      if (x_vert_m1==TmpLen[j])
        found++;
      if (x_vert_1==TmpLen[j])
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
          y1_vert_0 =  y1_vert_2 = x1_vert_m1 = x1_vert_1 = 0;
          // compute coordinates of vertices
          // count number of vertices which fulfill the criteria
          for (l1=0;l1<TmpLen[jj];l1++)
          {
            cell1->GetVertex(TmpFV[jj*MaxLen1+l1])->
              GetCoords(X1[l1], Y1[l1], Z1[l1]);
            // check if vertex on y = y_bottom
            if (fabs(Y1[l1]-y_bottom)<1e-5)
              y1_vert_0++;
            // check if vertex on y = y_top
            if (fabs(Y1[l1]-y_top)<1e-5)
              y1_vert_2++;
            // check if vertex on x = x_bottom
            if (fabs(X1[l1]-x_bottom)<1e-5)
              x1_vert_m1++;
            // check if vertex on x = x_top
            if (fabs(X1[l1]-x_top)<1e-5)
              x1_vert_1++;
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

          // all vertices of original face are on y_bottom
          // and all vertices of second face are on y_top
          if ((y_vert_0==TmpLen[j])&&(y1_vert_2==TmpLen[jj]))
          {
            // check if the x,z - coordinates of the barycenters are the same
            if  ((fabs(x-x1)>1e-5)||(fabs(z-z1)>1e-5))
              continue;
            // the face match
            found++;
          }
          // do the same vice versa
          if ((y_vert_2==TmpLen[j])&&(y1_vert_0==TmpLen[jj]))
          {
            if  ((fabs(x-x1)>1e-5)||(fabs(z-z1)>1e-5))
              continue;
            found++;
          }
          // all vertices of original face are on x_bottom
          // and all vertices of second face are on x_top
          if ((x_vert_m1==TmpLen[j])&&(x1_vert_1==TmpLen[jj]))
          {
            // check if the y,z - coordinates of the barycenters are the same
            if  ((fabs(z-z1)>1e-5)||(fabs(y-y1)>1e-5))
              continue;
            found++;
          }
          // do the same vice versa
          if ((x_vert_1==TmpLen[j])&&(x1_vert_m1==TmpLen[jj]))
          {
            if  ((fabs(z-z1)>1e-5)||(fabs(y-y1)>1e-5))
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
          if (((y_vert_0==TmpLen[j])&&(y1_vert_2==TmpLen[jj]))
            || ((y_vert_2==TmpLen[j])&&(y1_vert_0==TmpLen[jj])))
          {
            found = -1;
            for (l1=0;l1<TmpLen[jj];l1++)
            {
              if ((fabs(X[0]-X1[l1])<1e-5)&& (fabs(Z[0]-Z1[l1])<1e-5))
              {
                found = l1;
                break;
              }
            }
          }
          if (((x_vert_m1==TmpLen[j])&&(x1_vert_1==TmpLen[jj]))
            || ((x_vert_1==TmpLen[j])&&(x1_vert_m1==TmpLen[jj])))
          {
            found = -1;
            for (l1=0;l1<TmpLen[jj];l1++)
            {
              if ((fabs(Z[0]-Z1[l1])<1e-5)&& (fabs(Y[0]-Y1[l1])<1e-5))
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


// periodic b.c. :
//  x = -1 and x = 1
//  z = 0 and z = 2
void SetZCoordinates(TCollection *Coll, int level)
{
  int i, j, k, N_Cells, N_Layers, N_V, layer_int;
  TBaseCell *cell;
  TVertex *vertex;
  double x,y,z, layer;

  N_Cells = Coll->GetN_Cells();
  // number of layers on initial grid * 2^level
  N_Layers = TDatabase::ParamDB->N_CELL_LAYERS*(int)(pow(2.0,level));
  OutPut("number of layers " << N_Layers << endl);

  // first loop over cells
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    N_V = cell->GetN_Vertices();
    for (j=0;j<N_V;j++)
    {
      // read coordinates
      vertex = cell->GetVertex(j);
      vertex->GetCoords(x, y, z);
      // check if z coordinate fits to the prescribed distribution
      // compute inverse mapping
      //layer = acos(1-z)*N_Layers/Pi;
      layer = (atanh(tanh(2.75)*(z-1))/2.75+1)*N_Layers/2.0;
      layer_int = (int)(layer+1e-7);
      //  OutPut("z " << z << " " << layer << " ");
      // not the correct position
      if (fabs(layer-layer_int)>1e-5)
      {
        //OutPut(fabs(layer-layer_int) << " ");
        // find closest odd integer
        if (fabs(layer_int/2.0 - (int)(layer_int/2.0+1e-7))>0.1)
          k = layer_int;
        else
          k = layer_int+1;
        // compute new z coordinate
        //z = 1-cos(Pi*k/N_Layers);
        z = 1 + tanh(2.75*(2.0*k/N_Layers -1))/tanh(2.75);
       // set new z coordinate
        //OutPut(" k " << k << " z " << z << endl);
        vertex->SetCoords(x,y,z);
      }
    }
  }
  OutPut(endl);
}


void CheckZCoordinates(TCollection *Coll, int level)
{
  int i, j, k, N_Cells, N_Layers, N_V, found;
  TBaseCell *cell;
  TVertex *vertex;
  double x,y,z, layer, *zcoor;

  N_Cells = Coll->GetN_Cells();
  // number of layers on initial grid * 2^level
  N_Layers = TDatabase::ParamDB->N_CELL_LAYERS*(int)(pow(2.0,level));
  zcoor = new double[N_Layers+1];
  for (i=0;i<=N_Layers;i++)
      // zcoor[i] = 1-cos(Pi*i/N_Layers);
    zcoor[i] = 1 + tanh(2.75*(2.0*i/N_Layers -1))/tanh(2.75);

  // first loop over cells
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    N_V = cell->GetN_Vertices();
    for (j=0;j<N_V;j++)
    {
      // read coordinates
      vertex = cell->GetVertex(j);
      vertex->GetCoords(x, y, z);
      // check if z coordinate fits to the prescribed distribution
      // compute inverse mapping
      found = 0;
      for (k=0;k<=N_Layers;k++)
      {
        if (fabs(zcoor[k] - z)<1e-6)
        {
          found++;
          break;
        }
      }
      if (found)
        continue;
      OutPut("coordinate " << z << " not in list !!!");
      exit(4711);
    }
  }
}


// routine for computing the coordinates of the d.o.f. for a finite
// element space
// input: TCollection *Coll    -- collection of mesh cells
//        TFESpace3D *fespace  -- finite element space
// output: double *x_dof, double *y_dof, double *z_dof --
//         arrays for the coordinates (initialized with -4711 in each entry)
//         *N_z_layers : number of dof-layers in z-direction
//         *coord_z_layers: array with coordinates of layers in z-direction
void GetCoordinatesOfDof(TCollection *Coll, TFESpace3D *fespace,
double *x_dof, double *y_dof, double *z_dof,
int *N_z_layers, double* &coord_z_layers)
{
  int i,j, k, l, m, n, N_Cells, N_V, *global_numbers, *begin_index, *dof, N_;
  int *N_BaseFunct, max_z_layers = 1000, N_z_layers_help = 0, change = 0;
  double x[8],y[8],z[8], *coord_z_layers_help;
  TBaseCell *cell;
  TVertex *vertex;
  FE3D CurrentElement;

  // M->23.11.05
  double lx,ly,lz,help;
  int gdof;

  // allocate coord_z_layers_help
  coord_z_layers_help = new double[max_z_layers];

  // initialize coord_z_layers_help
  for(l=0;l<max_z_layers;l++)
  {
    coord_z_layers_help[l] = -4711.0;
  }

  // number of mesh cells
  N_Cells = Coll->GetN_Cells();
  // array with global numbers of d.o.f.
  global_numbers = fespace->GetGlobalNumbers();
  // array which points to the beginning of the global numbers in
  // global_numbers for each mesh cell
  begin_index = fespace->GetBeginIndex();
  // information on the number of basis functions for the available fe
  N_BaseFunct = TFEDatabase3D::GetN_BaseFunctFromFE3D();

  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    N_V = cell->GetN_Vertices();
    for (j=0;j<N_V;j++)
    {
      // read coordinates of the mesh cell
      vertex = cell->GetVertex(j);
      vertex->GetCoords(x[j], y[j], z[j]);
    }
    // finite element on the mesh cell
    CurrentElement = fespace->GetFE3D(i, cell);
    if ((CurrentElement!=C_Q2_3D_H_A)&&(CurrentElement!=C_Q2_3D_H_M)
	&&(CurrentElement!=C_Q3_3D_H_A)&&(CurrentElement!=C_Q3_3D_H_M)
	&&(CurrentElement!=C_Q4_3D_H_A)&&(CurrentElement!=C_Q4_3D_H_M))
    {
      OutPut("GetCoordinatesOfDof for finite element " << CurrentElement
        << " not implemented " << endl);
      exit(4711);
    }
    // number of basis functions (= number of d.o.f.)
    N_ = N_BaseFunct[CurrentElement];

    // the array which gives the mapping of the local to the global d.o.f.
    dof = global_numbers+begin_index[i];

    // loop over the local d.o.f. --> N_ = [ (k+1)^2-1 + k(k+1)^2 ]+1
    for (j=0;j<N_;j++)
    {
      switch(N_)                 // which kind of Q_k element --> k=2,3,4 ?  => N_ = 27,64,125
      {
        case 27:                 // for the case k=2 --> Q_2 elements
        {
          switch(j)              // compute coordinates of the local d.o.f. depending on their position
          {
            case 0:
              lx = x[0]; ly = y[0]; lz = z[0];
              break;
            case 1:
              lx = 0.5*(x[0]+x[1]); ly =0.5*(y[0]+y[1]); lz = 0.5*(z[0]+z[1]);
              break;
            case 2:
              lx = x[1]; ly = y[1]; lz = z[1];
              break;
            case 3:
              lx = 0.5*(x[0]+x[3]); ly = 0.5*(y[0]+y[3]); lz = 0.5*(z[0]+z[3]);
              break;
            case 4:
              lx = 0.25*(x[0]+x[1]+x[2]+x[3]); ly = 0.25*(y[0]+y[1]+y[2]+y[3]); lz = 0.25*(z[0]+z[1]+z[2]+z[3]);
              break;
            case 5:
              lx = 0.5*(x[1]+x[2]); ly = 0.5*(y[1]+y[2]); lz = 0.5*(z[1]+z[2]);
              break;
            case 6:
              lx = x[3]; ly = y[3]; lz = z[3];
              break;
            case 7:
              lx = 0.5*(x[2]+x[3]); ly = 0.5*(y[2]+y[3]); lz = 0.5*(z[2]+z[3]);
              break;
            case 8:
              lx = x[2]; ly = y[2]; lz = z[2];
              break;
            case 9:
              lx = 0.5*(x[0]+x[4]); ly = 0.5*(y[0]+y[4]); lz = 0.5*(z[0]+z[4]);
              break;
            case 10:
              lx = 0.25*(x[0]+x[1]+x[4]+x[5]); ly = 0.25*(y[0]+y[1]+y[4]+y[5]); lz = 0.25*(z[0]+z[1]+z[4]+z[5]);
              break;
            case 11:
              lx = 0.5*(x[1]+x[5]); ly = 0.5*(y[1]+y[5]); lz = 0.5*(z[1]+z[5]);
              break;
            case 12:
              lx = 0.25*(x[0]+x[3]+x[4]+x[7]); ly = 0.25*(y[0]+y[3]+y[4]+y[7]); lz = 0.25*(z[0]+z[3]+z[4]+z[7]);
              break;
            case 13:
              lx = 0.125*(x[0]+x[1]+x[2]+x[3]+x[4]+x[5]+x[6]+x[7]);
              ly = 0.125*(y[0]+y[1]+y[2]+y[3]+y[4]+y[5]+y[6]+y[7]);
              lz = 0.125*(z[0]+z[1]+z[2]+z[3]+z[4]+z[5]+z[6]+z[7]);
              break;
            case 14:
              lx = 0.25*(x[1]+x[2]+x[5]+x[6]); ly = 0.25*(y[1]+y[2]+y[5]+y[6]); lz = 0.25*(z[1]+z[2]+z[5]+z[6]);
              break;
            case 15:
              lx = 0.5*(x[3]+x[7]); ly = 0.5*(y[3]+y[7]); lz = 0.5*(z[3]+z[7]);
              break;
            case 16:
              lx = 0.25*(x[2]+x[3]+x[6]+x[7]); ly = 0.25*(y[2]+y[3]+y[6]+y[7]); lz = 0.25*(z[2]+z[3]+z[6]+z[7]);
              break;
            case 17:
              lx = 0.5*(x[2]+x[6]); ly = 0.5*(y[2]+y[6]); lz = 0.5*(z[2]+z[6]);
              break;
            case 18:
              lx = x[4]; ly = y[4]; lz = z[4];
              break;
            case 19:
              lx = 0.5*(x[4]+x[5]); ly = 0.5*(y[4]+y[5]); lz = 0.5*(z[4]+z[5]);
              break;
            case 20:
              lx = x[5]; ly = y[5]; lz = z[5];
              break;
            case 21:
              lx = 0.5*(x[4]+x[7]); ly = 0.5*(y[4]+y[7]); lz = 0.5*(z[4]+z[7]);
              break;
            case 22:
              lx = 0.25*(x[4]+x[5]+x[6]+x[7]); ly = 0.25*(y[4]+y[5]+y[6]+y[7]); lz = 0.25*(z[4]+z[5]+z[6]+z[7]);
              break;
            case 23:
              lx = 0.5*(x[5]+x[6]); ly = 0.5*(y[5]+y[6]); lz = 0.5*(z[5]+z[6]);
              break;
            case 24:
              lx = x[7]; ly = y[7]; lz = z[7];
              break;
            case 25:
              lx = 0.5*(x[6]+x[7]); ly = 0.5*(y[6]+y[7]); lz = 0.5*(z[6]+z[7]);
              break;
            case 26:
              lx = x[6]; ly = y[6]; lz = z[6];
	      break;
          }
         break;
        }
        case 64:                 // for the case k=3 --> Q_3 elements
        {
          switch(j)              // compute coordinates of the local d.o.f. depending on their position
          {
            case  0:
              lx = x[0];  ly = y[0];  lz = z[0];
              break;
            case  1:
              lx = (x[0]*2+x[1])/3;  ly = (y[0]*2+y[1])/3;  lz = (z[0]*2+z[1])/3;
              break;
            case  2:
              lx = (x[0]+x[1]*2)/3;  ly = (y[0]+y[1]*2)/3;  lz = (z[0]+z[1]*2)/3;
              break;
            case  3:
              lx = x[1];  ly = y[1];  lz = z[1];
              break;
            case  4:
              lx = (x[0]*2+x[3])/3;  ly = (y[0]*2+y[3])/3;  lz = (z[0]*2+z[3])/3;
              break;
            case  5:
              lx = (((x[0]*2+x[3])/3)*2+(x[1]*2+x[2])/3)/3;
              ly = (((y[0]*2+y[3])/3)*2+(y[1]*2+y[2])/3)/3;
              lz = (((z[0]*2+z[3])/3)*2+(z[1]*2+z[2])/3)/3;
              break;
            case  6:
              lx = ((x[0]*2+x[3])/3+2*(x[1]*2+x[2])/3)/3;
              ly = ((y[0]*2+y[3])/3+2*(y[1]*2+y[2])/3)/3;
              lz = ((z[0]*2+z[3])/3+2*(z[1]*2+z[2])/3)/3;
              break;
            case  7:
              lx = (x[1]*2+x[2])/3;  ly = (y[1]*2+y[2])/3;  lz = (z[1]*2+z[2])/3;
              break;
            case  8:
              lx = (x[0]+x[3]*2)/3;  ly = (y[0]+y[3]*2)/3;  lz = (z[0]+z[3]*2)/3;
              break;
            case  9:
              lx = (2*(x[0]+x[3]*2)/3+(x[1]+x[2]*2)/3)/3;
              ly = (2*(y[0]+y[3]*2)/3+(y[1]+y[2]*2)/3)/3;
              lz = (2*(z[0]+z[3]*2)/3+(z[1]+z[2]*2)/3)/3;
              break;
            case 10:
              lx = ((x[0]+x[3]*2)/3+2*(x[1]+x[2]*2)/3)/3;
              ly = ((y[0]+y[3]*2)/3+2*(y[1]+y[2]*2)/3)/3;
              lz = ((z[0]+z[3]*2)/3+2*(z[1]+z[2]*2)/3)/3;
              break;
            case 11:
              lx = (x[1]+x[2]*2)/3;  ly = (y[1]+y[2]*2)/3;  lz = (z[1]+z[2]*2)/3;
              break;
            case 12:
              lx = x[3];  ly = y[3];  lz = z[3];
              break;
            case 13:
              lx = (x[3]*2+x[2])/3;  ly = (y[3]*2+y[2])/3;  lz = (z[3]*2+z[2])/3;
              break;
            case 14:
              lx = (x[3]+x[2]*2)/3;  ly = (y[3]+y[2]*2)/3;  lz = (z[3]+z[2]*2)/3;
              break;
            case 15:
              lx = x[2];  ly = y[2];  lz = z[2];
              break;
            case 16:
              lx = (x[0]*2+x[4])/3;  ly = (y[0]*2+y[4])/3;  lz = (z[0]*2+z[4])/3;
              break;
            case 17:
              lx = (2*(x[0]*2+x[4])/3+(x[1]*2+x[5])/3)/3;
              ly = (2*(y[0]*2+y[4])/3+(y[1]*2+y[5])/3)/3;
              lz = (2*(z[0]*2+z[4])/3+(z[1]*2+z[5])/3)/3;
              break;
            case 18:
              lx = ((x[0]*2+x[4])/3+2*(x[1]*2+x[5])/3)/3;
              ly = ((y[0]*2+y[4])/3+2*(y[1]*2+y[5])/3)/3;
              lz = ((z[0]*2+z[4])/3+2*(z[1]*2+z[5])/3)/3;
              break;
            case 19:
              lx = (x[1]*2+x[5])/3;  ly = (y[1]*2+y[5])/3;  lz = (z[1]*2+z[5])/3;
              break;
            case 20:
              lx = (2*(x[0]*2+x[4])/3+(x[3]*2+x[7])/3)/3;
              ly = (2*(y[0]*2+y[4])/3+(y[3]*2+y[7])/3)/3;
              lz = (2*(z[0]*2+z[4])/3+(z[3]*2+z[7])/3)/3;
              break;
            case 21:
              lx = (2*((2*(x[0]*2+x[4])/3+(x[3]*2+x[7])/3)/3)+((2*(x[1]*2+x[5])/3+(x[2]*2+x[6])/3)/3))/3;
              ly = (2*((2*(y[0]*2+y[4])/3+(y[3]*2+y[7])/3)/3)+((2*(y[1]*2+y[5])/3+(y[2]*2+y[6])/3)/3))/3;
              lz = (2*((2*(z[0]*2+z[4])/3+(z[3]*2+z[7])/3)/3)+((2*(z[1]*2+z[5])/3+(z[2]*2+z[6])/3)/3))/3;
              break;
            case 22:
              lx = (2*((2*(x[1]*2+x[5])/3+(x[2]*2+x[6])/3)/3)+((2*(x[0]*2+x[4])/3+(x[3]*2+x[7])/3)/3))/3;
              ly = (2*((2*(y[1]*2+y[5])/3+(y[2]*2+y[6])/3)/3)+((2*(y[0]*2+y[4])/3+(y[3]*2+y[7])/3)/3))/3;
              lz = (2*((2*(z[1]*2+z[5])/3+(z[2]*2+z[6])/3)/3)+((2*(z[0]*2+z[4])/3+(z[3]*2+z[7])/3)/3))/3;
              break;
            case 23:
              lx = (2*(x[1]*2+x[5])/3+(x[2]*2+x[6])/3)/3;
              ly = (2*(y[1]*2+y[5])/3+(y[2]*2+y[6])/3)/3;
              lz = (2*(z[1]*2+z[5])/3+(z[2]*2+z[6])/3)/3;
              break;
            case 24:
              lx = (2*(x[3]*2+x[7])/3+(x[0]*2+x[4])/3)/3;
              ly = (2*(y[3]*2+y[7])/3+(y[0]*2+y[4])/3)/3;
              lz = (2*(z[3]*2+z[7])/3+(z[0]*2+z[4])/3)/3;
              break;
            case 25:
              lx = (2*((2*(x[3]*2+x[7])/3+(x[0]*2+x[4])/3)/3)+((2*(x[1]*2+x[5])/3+(x[2]*2+x[6])/3)/3))/3;
              ly = (2*((2*(y[3]*2+y[7])/3+(y[0]*2+y[4])/3)/3)+((2*(y[1]*2+y[5])/3+(y[2]*2+y[6])/3)/3))/3;
              lz = (2*((2*(z[3]*2+z[7])/3+(z[0]*2+z[4])/3)/3)+((2*(z[1]*2+z[5])/3+(z[2]*2+z[6])/3)/3))/3;
              break;
            case 26:
              lx = (2*((2*(x[1]*2+x[5])/3+(x[2]*2+x[6])/3)/3)+((2*(x[3]*2+x[7])/3+(x[0]*2+x[4])/3)/3))/3;
              ly = (2*((2*(y[1]*2+y[5])/3+(y[2]*2+y[6])/3)/3)+((2*(y[3]*2+y[7])/3+(y[0]*2+y[4])/3)/3))/3;
              lz = (2*((2*(z[1]*2+z[5])/3+(z[2]*2+z[6])/3)/3)+((2*(z[3]*2+z[7])/3+(z[0]*2+z[4])/3)/3))/3;
              break;
            case 27:
              lx = ((x[1]*2+x[5])/3+2*(x[2]*2+x[6])/3)/3;
              ly = ((y[1]*2+y[5])/3+2*(y[2]*2+y[6])/3)/3;
              lz = ((z[1]*2+z[5])/3+2*(z[2]*2+z[6])/3)/3;
              break;
            case 28:
              lx = (x[3]*2+x[7])/3;  ly = (y[3]*2+y[7])/3;  lz = (z[3]*2+z[7])/3;
              break;
            case 29:
              lx = (2*(x[3]*2+x[7])/3+(x[2]*2+x[6])/3)/3;
              ly = (2*(y[3]*2+y[7])/3+(y[2]*2+y[6])/3)/3;
              lz = (2*(z[3]*2+z[7])/3+(z[2]*2+z[6])/3)/3;
              break;
            case 30:
              lx = (2*(x[2]*2+x[6])/3+(x[3]*2+x[7])/3)/3;
              ly = (2*(y[2]*2+y[6])/3+(y[3]*2+y[7])/3)/3;
              lz = (2*(z[2]*2+z[6])/3+(z[3]*2+z[7])/3)/3;
              break;
            case 31:
              lx = (x[2]*2+x[6])/3;  ly = (y[2]*2+y[6])/3;  lz = (z[2]*2+z[6])/3;
              break;
            case 32:
              lx = (2*x[4]+x[0])/3;  ly = (2*y[4]+y[0])/3;  lz = (2*z[4]+z[0])/3;
              break;
            case 33:
              lx = (2*(2*x[4]+x[0])/3+(2*x[5]+x[1])/3)/3;
              ly = (2*(2*y[4]+y[0])/3+(2*y[5]+y[1])/3)/3;
              lz = (2*(2*z[4]+z[0])/3+(2*z[5]+z[1])/3)/3;
              break;
            case 34:
              lx = (2*(2*x[5]+x[1])/3+(2*x[4]+x[0])/3)/3;
              ly = (2*(2*y[5]+y[1])/3+(2*y[4]+y[0])/3)/3;
              lz = (2*(2*z[5]+z[1])/3+(2*z[4]+z[0])/3)/3;
              break;
            case 35:
              lx = (2*x[5]+x[1])/3;  ly = (2*y[5]+y[1])/3;  lz = (2*z[5]+z[1])/3;
              break;
            case 36:
              lx = (2*(2*x[4]+x[0])/3+(2*x[7]+x[3])/3)/3;
              ly = (2*(2*y[4]+y[0])/3+(2*y[7]+y[3])/3)/3;
              lz = (2*(2*z[4]+z[0])/3+(2*z[7]+z[3])/3)/3;
              break;
            case 37:
              lx = (2*((2*(2*x[4]+x[0])/3+(2*x[7]+x[3])/3)/3)+((2*(2*x[5]+x[1])/3+(2*x[6]+x[2])/3)/3))/3;
              ly = (2*((2*(2*y[4]+y[0])/3+(2*y[7]+y[3])/3)/3)+((2*(2*y[5]+y[1])/3+(2*y[6]+y[2])/3)/3))/3;
              lz = (2*((2*(2*z[4]+z[0])/3+(2*z[7]+z[3])/3)/3)+((2*(2*z[5]+z[1])/3+(2*z[6]+z[2])/3)/3))/3;
              break;
            case 38:
              lx = (2*((2*(2*x[5]+x[1])/3+(2*x[6]+x[2])/3)/3)+((2*(2*x[4]+x[0])/3+(2*x[7]+x[3])/3)/3))/3;
              ly = (2*((2*(2*y[5]+y[1])/3+(2*y[6]+y[2])/3)/3)+((2*(2*y[4]+y[0])/3+(2*y[7]+y[3])/3)/3))/3;
              lz = (2*((2*(2*z[5]+z[1])/3+(2*z[6]+z[2])/3)/3)+((2*(2*z[4]+z[0])/3+(2*z[7]+z[3])/3)/3))/3;
              break;
            case 39:
              lx = (2*(2*x[5]+x[1])/3+(2*x[6]+x[2])/3)/3;
              ly = (2*(2*y[5]+y[1])/3+(2*y[6]+y[2])/3)/3;
              lz = (2*(2*z[5]+z[1])/3+(2*z[6]+z[2])/3)/3;
              break;
            case 40:
              lx = (2*(2*x[7]+x[3])/3+(2*x[4]+x[0])/3)/3;
              ly = (2*(2*y[7]+y[3])/3+(2*y[4]+y[0])/3)/3;
              lz = (2*(2*z[7]+z[3])/3+(2*z[4]+z[0])/3)/3;
              break;
            case 41:
              lx = (2*((2*(2*x[7]+x[3])/3+(2*x[4]+x[0])/3)/3)+((2*(2*x[6]+x[2])/3+(2*x[5]+x[1])/3)/3))/3;
              ly = (2*((2*(2*y[7]+y[3])/3+(2*y[4]+y[0])/3)/3)+((2*(2*y[6]+y[2])/3+(2*y[5]+y[1])/3)/3))/3;
              lz = (2*((2*(2*z[7]+z[3])/3+(2*z[4]+z[0])/3)/3)+((2*(2*z[6]+z[2])/3+(2*z[5]+z[1])/3)/3))/3;
              break;
            case 42:
              lx = (2*((2*(2*x[6]+x[2])/3+(2*x[5]+x[1])/3)/3)+((2*(2*x[7]+x[3])/3+(2*x[4]+x[0])/3)/3))/3;
              ly = (2*((2*(2*y[6]+y[2])/3+(2*y[5]+y[1])/3)/3)+((2*(2*y[7]+y[3])/3+(2*y[4]+y[0])/3)/3))/3;
              lz = (2*((2*(2*z[6]+z[2])/3+(2*z[5]+z[1])/3)/3)+((2*(2*z[7]+z[3])/3+(2*z[4]+z[0])/3)/3))/3;
              break;
            case 43:
              lx = (2*(2*x[6]+x[2])/3+(2*x[5]+x[1])/3)/3;
              ly = (2*(2*y[6]+y[2])/3+(2*y[5]+y[1])/3)/3;
              lz = (2*(2*z[6]+z[2])/3+(2*z[5]+z[1])/3)/3;
              break;
            case 44:
              lx = (2*x[7]+x[3])/3;  ly = (2*y[7]+y[3])/3;  lz = (2*z[7]+z[3])/3;
              break;
            case 45:
              lx = (2*(2*x[7]+x[3])/3+(2*x[6]+x[2])/3)/3;
              ly = (2*(2*y[7]+y[3])/3+(2*y[6]+y[2])/3)/3;
              lz = (2*(2*z[7]+z[3])/3+(2*z[6]+z[2])/3)/3;
              break;
            case 46:
              lx = (2*(2*x[6]+x[2])/3+(2*x[7]+x[3])/3)/3;
              ly = (2*(2*y[6]+y[2])/3+(2*y[7]+y[3])/3)/3;
              lz = (2*(2*z[6]+z[2])/3+(2*z[7]+z[3])/3)/3;
              break;
            case 47:
              lx = (2*x[6]+x[2])/3;  ly = (2*y[6]+y[2])/3;  lz = (2*z[6]+z[2])/3;
              break;
            case 48:
              lx = x[4];  ly = y[4];  lz = z[4];
              break;
            case 49:
              lx = (2*x[4]+x[5])/3;  ly = (2*y[4]+y[5])/3;  lz = (2*z[4]+z[5])/3;
              break;
            case 50:
              lx = (2*x[5]+x[4])/3;  ly = (2*y[5]+y[4])/3;  lz = (2*z[5]+z[4])/3;
              break;
            case 51:
              lx = x[5];  ly = y[5];  lz = z[5];
              break;
            case 52:
              lx = (2*x[4]+x[7])/3;  ly = (2*y[4]+y[7])/3;  lz = (2*z[4]+z[7])/3;
              break;
            case 53:
              lx = (2*(2*x[4]+x[7])/3+(2*x[5]+x[6])/3)/3;
              ly = (2*(2*y[4]+y[7])/3+(2*y[5]+y[6])/3)/3;
              lz = (2*(2*z[4]+z[7])/3+(2*z[5]+z[6])/3)/3;
              break;
            case 54:
              lx = (2*(2*x[5]+x[6])/3+(2*x[4]+x[7])/3)/3;
              ly = (2*(2*y[5]+y[6])/3+(2*y[4]+y[7])/3)/3;
              lz = (2*(2*z[5]+z[6])/3+(2*z[4]+z[7])/3)/3;
              break;
            case 55:
              lx = (2*x[5]+x[6])/3;  ly = (2*y[5]+y[6])/3;  lz = (2*z[5]+z[6])/3;
              break;
            case 56:
              lx = (2*x[7]+x[4])/3;  ly = (2*y[7]+y[4])/3;  lz = (2*z[7]+z[4])/3;
              break;
            case 57:
              lx = (2*(2*x[7]+x[4])/3+(2*x[6]+x[5])/3)/3;
              ly = (2*(2*y[7]+y[4])/3+(2*y[6]+y[5])/3)/3;
              lz = (2*(2*z[7]+z[4])/3+(2*z[6]+z[5])/3)/3;
              break;
            case 58:
              lx = (2*(2*x[6]+x[5])/3+(2*x[7]+x[4])/3)/3;
              ly = (2*(2*y[6]+y[5])/3+(2*y[7]+y[4])/3)/3;
              lz = (2*(2*z[6]+z[5])/3+(2*z[7]+z[4])/3)/3;
              break;
            case 59:
              lx = (2*x[6]+x[5])/3;  ly = (2*y[6]+y[5])/3;  lz = (2*z[6]+z[5])/3;
              break;
            case 60:
              lx = x[7];  ly = y[7];  lz = z[7];
              break;
            case 61:
              lx = (2*x[7]+x[6])/3;  ly = (2*y[7]+y[6])/3;  lz = (2*z[7]+z[6])/3;
              break;
            case 62:
              lx = (2*x[6]+x[7])/3;  ly = (2*y[6]+y[7])/3;  lz = (2*z[6]+z[7])/3;
              break;
            case 63:
              lx = x[6];  ly = y[6];  lz = z[6];
              break;
          }
	  break;
        }
        case 125:                // for the case k=4 --> Q_4 elements
        {
          switch(j)
          {
            case   0:
              lx = x[0];  ly = y[0];  lz = z[0];
              break;
            case   1:
              lx = (3*x[0]+x[1])/4;  ly = (3*y[0]+y[1])/4;  lz = (3*z[0]+z[1])/4;
              break;
            case   2:
              lx = (2*x[0]+2*x[1])/4;  ly = (2*y[0]+2*y[1])/4;  lz = (2*z[0]+2*z[1])/4;
              break;
            case   3:
              lx = (x[0]+3*x[1])/4;  ly = (y[0]+3*y[1])/4;  lz = (z[0]+3*z[1])/4;
              break;
            case   4:
              lx = x[1];  ly = y[1];  lz = z[1];
              break;
            case   5:
              lx = (3*x[0]+x[3])/4;  ly = (3*y[0]+y[3])/4;  lz = (3*z[0]+z[3])/4;
              break;
            case   6:
              lx = (3*((3*x[0]+x[3])/4)+((3*x[1]+x[2])/4))/4;
              ly = (3*((3*y[0]+y[3])/4)+((3*y[1]+y[2])/4))/4;
              lz = (3*((3*z[0]+z[3])/4)+((3*z[1]+z[2])/4))/4;
              break;
            case   7:
              lx = (2*((3*x[0]+x[3])/4)+2*((3*x[1]+x[2])/4))/4;
              ly = (2*((3*y[0]+y[3])/4)+2*((3*y[1]+y[2])/4))/4;
              lz = (2*((3*z[0]+z[3])/4)+2*((3*z[1]+z[2])/4))/4;
              break;
            case   8:
              lx = (((3*x[0]+x[3])/4)+3*((3*x[1]+x[2])/4))/4;
              ly = (((3*y[0]+y[3])/4)+3*((3*y[1]+y[2])/4))/4;
              lz = (((3*z[0]+z[3])/4)+3*((3*z[1]+z[2])/4))/4;
              break;
            case   9:
              lx = (3*x[1]+x[2])/4;  ly = (3*y[1]+y[2])/4;  lz = (3*z[1]+z[2])/4;
              break;
            case  10:
              lx = (2*x[0]+2*x[3])/4;  ly = (2*y[0]+2*y[3])/4;  lz = (2*z[0]+2*z[3])/4;
              break;
            case  11:
              lx = (3*((2*x[0]+2*x[3])/4)+((2*x[1]+2*x[2])/4))/4;
              ly = (3*((2*y[0]+2*y[3])/4)+((2*y[1]+2*y[2])/4))/4;
              lz = (3*((2*z[0]+2*z[3])/4)+((2*z[1]+2*z[2])/4))/4;
              break;
            case  12:
              lx = (2*((2*x[0]+2*x[3])/4)+2*((2*x[1]+2*x[2])/4))/4;
              ly = (2*((2*y[0]+2*y[3])/4)+2*((2*y[1]+2*y[2])/4))/4;
              lz = (2*((2*z[0]+2*z[3])/4)+2*((2*z[1]+2*z[2])/4))/4;
              break;
            case  13:
              lx = (((2*x[0]+2*x[3])/4)+3*((2*x[1]+2*x[2])/4))/4;
              ly = (((2*y[0]+2*y[3])/4)+3*((2*y[1]+2*y[2])/4))/4;
              lz = (((2*z[0]+2*z[3])/4)+3*((2*z[1]+2*z[2])/4))/4;
              break;
            case  14:
              lx = (2*x[1]+2*x[2])/4;  ly = (2*y[1]+2*y[2])/4;  lz = (2*z[1]+2*z[2])/4;
              break;
            case  15:
              lx = (x[0]+3*x[3])/4;  ly = (y[0]+3*y[3])/4;  lz = (z[0]+3*z[3])/4;
              break;
            case  16:
              lx = (3*((x[0]+3*x[3])/4)+((x[1]+3*x[2])/4))/4;
              ly = (3*((y[0]+3*y[3])/4)+((y[1]+3*y[2])/4))/4;
              lz = (3*((z[0]+3*z[3])/4)+((z[1]+3*z[2])/4))/4;
              break;
            case  17:
              lx = (2*((x[0]+3*x[3])/4)+2*((x[1]+3*x[2])/4))/4;
              ly = (2*((y[0]+3*y[3])/4)+2*((y[1]+3*y[2])/4))/4;
              lz = (2*((z[0]+3*z[3])/4)+2*((z[1]+3*z[2])/4))/4;
              break;
            case  18:
              lx = (((x[0]+3*x[3])/4)+3*((x[1]+3*x[2])/4))/4;
              ly = (((y[0]+3*y[3])/4)+3*((y[1]+3*y[2])/4))/4;
              lz = (((z[0]+3*z[3])/4)+3*((z[1]+3*z[2])/4))/4;
              break;
            case  19:
              lx = (x[1]+3*x[2])/4;  ly = (y[1]+3*y[2])/4;  lz = (z[1]+3*z[2])/4;
              break;
            case  20:
              lx = x[3];  ly = y[3];  lz = z[3];
              break;
            case  21:
              lx = (3*x[3]+x[2])/4;  ly = (3*y[3]+y[2])/4;  lz = (3*z[3]+z[2])/4;
              break;
            case  22:
              lx = (2*x[3]+2*x[2])/4;  ly = (2*y[3]+2*y[2])/4;  lz = (2*z[3]+2*z[2])/4;
              break;
            case  23:
              lx = (x[3]+3*x[2])/4;  ly = (y[3]+3*y[2])/4;  lz = (z[3]+3*z[2])/4;
              break;
            case  24:
              lx = x[2];  ly = y[2];  lz = z[2];
              break;
            case  25:
              lx = (3*x[0]+x[4])/4;  ly = (3*y[0]+y[4])/4;  lz = (3*z[0]+z[4])/4;
              break;
            case  26:
              lx = (3*((3*x[0]+x[4])/4)+((3*x[1]+x[5])/4))/4;
              ly = (3*((3*y[0]+y[4])/4)+((3*y[1]+y[5])/4))/4;
              lz = (3*((3*z[0]+z[4])/4)+((3*z[1]+z[5])/4))/4;
              break;
            case  27:
              lx = (2*((3*x[0]+x[4])/4)+2*((3*x[1]+x[5])/4))/4;
              ly = (2*((3*y[0]+y[4])/4)+2*((3*y[1]+y[5])/4))/4;
              lz = (2*((3*z[0]+z[4])/4)+2*((3*z[1]+z[5])/4))/4;
              break;
            case  28:
              lx = (((3*x[0]+x[4])/4)+3*((3*x[1]+x[5])/4))/4;
              ly = (((3*y[0]+y[4])/4)+3*((3*y[1]+y[5])/4))/4;
              lz = (((3*z[0]+z[4])/4)+3*((3*z[1]+z[5])/4))/4;
              break;
            case  29:
              lx = (3*x[1]+x[5])/4;  ly = (3*y[1]+y[5])/4;  lz = (3*z[1]+z[5])/4;
              break;
            case  30:
              lx = (3*((3*x[0]+x[4])/4)+((3*x[3]+x[7])/4))/4;
              ly = (3*((3*y[0]+y[4])/4)+((3*y[3]+y[7])/4))/4;
              lz = (3*((3*z[0]+z[4])/4)+((3*z[3]+z[7])/4))/4;
              break;
            case  31:
              lx = (3*((3*((3*x[0]+x[4])/4)+((3*x[3]+x[7])/4))/4)+((3*((3*x[1]+x[5])/4)+((3*x[2]+x[6])/4))/4))/4;
              ly = (3*((3*((3*y[0]+y[4])/4)+((3*y[3]+y[7])/4))/4)+((3*((3*y[1]+y[5])/4)+((3*y[2]+y[6])/4))/4))/4;
              lz = (3*((3*((3*z[0]+z[4])/4)+((3*z[3]+z[7])/4))/4)+((3*((3*z[1]+z[5])/4)+((3*z[2]+z[6])/4))/4))/4;
              break;
            case  32:
              lx = (2*((3*((3*x[0]+x[4])/4)+((3*x[3]+x[7])/4))/4)+2*((3*((3*x[1]+x[5])/4)+((3*x[2]+x[6])/4))/4))/4;
              ly = (2*((3*((3*y[0]+y[4])/4)+((3*y[3]+y[7])/4))/4)+2*((3*((3*y[1]+y[5])/4)+((3*y[2]+y[6])/4))/4))/4;
              lz = (2*((3*((3*z[0]+z[4])/4)+((3*z[3]+z[7])/4))/4)+2*((3*((3*z[1]+z[5])/4)+((3*z[2]+z[6])/4))/4))/4;
              break;
            case  33:
              lx = (((3*((3*x[0]+x[4])/4)+((3*x[3]+x[7])/4))/4)+3*((3*((3*x[1]+x[5])/4)+((3*x[2]+x[6])/4))/4))/4;
              ly = (((3*((3*y[0]+y[4])/4)+((3*y[3]+y[7])/4))/4)+3*((3*((3*y[1]+y[5])/4)+((3*y[2]+y[6])/4))/4))/4;
              lz = (((3*((3*z[0]+z[4])/4)+((3*z[3]+z[7])/4))/4)+3*((3*((3*z[1]+z[5])/4)+((3*z[2]+z[6])/4))/4))/4;
              break;
            case  34:
              lx = (3*((3*x[1]+x[5])/4)+((3*x[2]+x[6])/4))/4;
              ly = (3*((3*y[1]+y[5])/4)+((3*y[2]+y[6])/4))/4;
              lz = (3*((3*z[1]+z[5])/4)+((3*z[2]+z[6])/4))/4;
              break;
            case  35:
              lx = (2*((3*x[0]+x[4])/4)+2*((3*x[3]+x[7])/4))/4;
              ly = (2*((3*y[0]+y[4])/4)+2*((3*y[3]+y[7])/4))/4;
              lz = (2*((3*z[0]+z[4])/4)+2*((3*z[3]+z[7])/4))/4;
              break;
            case  36:
              lx = (3*((2*((3*x[0]+x[4])/4)+2*((3*x[3]+x[7])/4))/4)+((2*((3*x[1]+x[5])/4)+2*((3*x[2]+x[6])/4))/4))/4;
              ly = (3*((2*((3*y[0]+y[4])/4)+2*((3*y[3]+y[7])/4))/4)+((2*((3*y[1]+y[5])/4)+2*((3*y[2]+y[6])/4))/4))/4;
              lz = (3*((2*((3*z[0]+z[4])/4)+2*((3*z[3]+z[7])/4))/4)+((2*((3*z[1]+z[5])/4)+2*((3*z[2]+z[6])/4))/4))/4;
              break;
            case  37:
              lx = (2*((2*((3*x[0]+x[4])/4)+2*((3*x[3]+x[7])/4))/4)+2*((2*((3*x[1]+x[5])/4)+2*((3*x[2]+x[6])/4))/4))/4;
              ly = (2*((2*((3*y[0]+y[4])/4)+2*((3*y[3]+y[7])/4))/4)+2*((2*((3*y[1]+y[5])/4)+2*((3*y[2]+y[6])/4))/4))/4;
              lz = (2*((2*((3*z[0]+z[4])/4)+2*((3*z[3]+z[7])/4))/4)+2*((2*((3*z[1]+z[5])/4)+2*((3*z[2]+z[6])/4))/4))/4;
              break;
            case  38:
              lx = (((2*((3*x[0]+x[4])/4)+2*((3*x[3]+x[7])/4))/4)+3*((2*((3*x[1]+x[5])/4)+2*((3*x[2]+x[6])/4))/4))/4;
              ly = (((2*((3*y[0]+y[4])/4)+2*((3*y[3]+y[7])/4))/4)+3*((2*((3*y[1]+y[5])/4)+2*((3*y[2]+y[6])/4))/4))/4;
              lz = (((2*((3*z[0]+z[4])/4)+2*((3*z[3]+z[7])/4))/4)+3*((2*((3*z[1]+z[5])/4)+2*((3*z[2]+z[6])/4))/4))/4;
              break;
            case  39:
              lx = (2*((3*x[1]+x[5])/4)+2*((3*x[2]+x[6])/4))/4;
              ly = (2*((3*y[1]+y[5])/4)+2*((3*y[2]+y[6])/4))/4;
              lz = (2*((3*z[1]+z[5])/4)+2*((3*z[2]+z[6])/4))/4;
              break;
            case  40:
              lx = (((3*x[0]+x[4])/4)+3*((3*x[3]+x[7])/4))/4;
              ly = (((3*y[0]+y[4])/4)+3*((3*y[3]+y[7])/4))/4;
              lz = (((3*z[0]+z[4])/4)+3*((3*z[3]+z[7])/4))/4;
              break;
            case  41:
              lx = (3*((((3*x[0]+x[4])/4)+3*((3*x[3]+x[7])/4))/4)+((((3*x[1]+x[5])/4)+3*((3*x[2]+x[6])/4))/4))/4;
              ly = (3*((((3*y[0]+y[4])/4)+3*((3*y[3]+y[7])/4))/4)+((((3*y[1]+y[5])/4)+3*((3*y[2]+y[6])/4))/4))/4;
              lz = (3*((((3*z[0]+z[4])/4)+3*((3*z[3]+z[7])/4))/4)+((((3*z[1]+z[5])/4)+3*((3*z[2]+z[6])/4))/4))/4;
              break;
            case  42:
              lx = (2*((((3*x[0]+x[4])/4)+3*((3*x[3]+x[7])/4))/4)+2*((((3*x[1]+x[5])/4)+3*((3*x[2]+x[6])/4))/4))/4;
              ly = (2*((((3*y[0]+y[4])/4)+3*((3*y[3]+y[7])/4))/4)+2*((((3*y[1]+y[5])/4)+3*((3*y[2]+y[6])/4))/4))/4;
              lz = (2*((((3*z[0]+z[4])/4)+3*((3*z[3]+z[7])/4))/4)+2*((((3*z[1]+z[5])/4)+3*((3*z[2]+z[6])/4))/4))/4;
              break;
            case  43:
              lx = (((((3*x[0]+x[4])/4)+3*((3*x[3]+x[7])/4))/4)+3*((((3*x[1]+x[5])/4)+3*((3*x[2]+x[6])/4))/4))/4;
              ly = (((((3*y[0]+y[4])/4)+3*((3*y[3]+y[7])/4))/4)+3*((((3*y[1]+y[5])/4)+3*((3*y[2]+y[6])/4))/4))/4;
              lz = (((((3*z[0]+z[4])/4)+3*((3*z[3]+z[7])/4))/4)+3*((((3*z[1]+z[5])/4)+3*((3*z[2]+z[6])/4))/4))/4;
              break;
            case  44:
              lx = (((3*x[1]+x[5])/4)+3*((3*x[2]+x[6])/4))/4;
              ly = (((3*y[1]+y[5])/4)+3*((3*y[2]+y[6])/4))/4;
              lz = (((3*z[1]+z[5])/4)+3*((3*z[2]+z[6])/4))/4;
              break;
            case  45:
              lx = (3*x[3]+x[7])/4;  ly = (3*y[3]+y[7])/4;  lz = (3*z[3]+z[7])/4;
              break;
            case  46:
              lx = (3*((3*x[3]+x[7])/4)+((3*x[2]+x[6])/4))/4;
              ly = (3*((3*y[3]+y[7])/4)+((3*y[2]+y[6])/4))/4;
              lz = (3*((3*z[3]+z[7])/4)+((3*z[2]+z[6])/4))/4;
              break;
            case  47:
              lx = (2*((3*x[3]+x[7])/4)+2*((3*x[2]+x[6])/4))/4;
              ly = (2*((3*y[3]+y[7])/4)+2*((3*y[2]+y[6])/4))/4;
              lz = (2*((3*z[3]+z[7])/4)+2*((3*z[2]+z[6])/4))/4;
              break;
            case  48:
              lx = (((3*x[3]+x[7])/4)+3*((3*x[2]+x[6])/4))/4;
              ly = (((3*y[3]+y[7])/4)+3*((3*y[2]+y[6])/4))/4;
              lz = (((3*z[3]+z[7])/4)+3*((3*z[2]+z[6])/4))/4;
              break;
            case  49:
              lx = (3*x[2]+x[6])/4;  ly = (3*y[2]+y[6])/4;  lz = (3*z[2]+z[6])/4;
              break;
            case  50:
              lx = (2*x[0]+2*x[4])/4;  ly = (2*y[0]+2*y[4])/4;  lz = (2*z[0]+2*z[4])/4;
              break;
            case  51:
              lx = (3*((2*x[0]+2*x[4])/4)+((2*x[1]+2*x[5])/4))/4;
              ly = (3*((2*y[0]+2*y[4])/4)+((2*y[1]+2*y[5])/4))/4;
              lz = (3*((2*z[0]+2*z[4])/4)+((2*z[1]+2*z[5])/4))/4;
              break;
            case  52:
              lx = (2*((2*x[0]+2*x[4])/4)+2*((2*x[1]+2*x[5])/4))/4;
              ly = (2*((2*y[0]+2*y[4])/4)+2*((2*y[1]+2*y[5])/4))/4;
              lz = (2*((2*z[0]+2*z[4])/4)+2*((2*z[1]+2*z[5])/4))/4;
              break;
            case  53:
              lx = (((2*x[0]+2*x[4])/4)+3*((2*x[1]+2*x[5])/4))/4;
              ly = (((2*y[0]+2*y[4])/4)+3*((2*y[1]+2*y[5])/4))/4;
              lz = (((2*z[0]+2*z[4])/4)+3*((2*z[1]+2*z[5])/4))/4;
              break;
            case  54:
              lx = (2*x[1]+2*x[5])/4;  ly = (2*y[1]+2*y[5])/4;  lz = (2*z[1]+2*z[5])/4;
              break;
            case  55:
              lx = (3*((2*x[0]+2*x[4])/4)+((2*x[3]+2*x[7])/4))/4;
              ly = (3*((2*y[0]+2*y[4])/4)+((2*y[3]+2*y[7])/4))/4;
              lz = (3*((2*z[0]+2*z[4])/4)+((2*z[3]+2*z[7])/4))/4;
              break;
            case  56:
              lx = (3*((3*((2*x[0]+2*x[4])/4)+((2*x[3]+2*x[7])/4))/4)+((3*((2*x[1]+2*x[5])/4)+((2*x[2]+2*x[6])/4))/4))/4;
              ly = (3*((3*((2*y[0]+2*y[4])/4)+((2*y[3]+2*y[7])/4))/4)+((3*((2*y[1]+2*y[5])/4)+((2*y[2]+2*y[6])/4))/4))/4;
              lz = (3*((3*((2*z[0]+2*z[4])/4)+((2*z[3]+2*z[7])/4))/4)+((3*((2*z[1]+2*z[5])/4)+((2*z[2]+2*z[6])/4))/4))/4;
              break;
            case  57:
              lx = (2*((3*((2*x[0]+2*x[4])/4)+((2*x[3]+2*x[7])/4))/4)+2*((3*((2*x[1]+2*x[5])/4)+((2*x[2]+2*x[6])/4))/4))/4;
              ly = (2*((3*((2*y[0]+2*y[4])/4)+((2*y[3]+2*y[7])/4))/4)+2*((3*((2*y[1]+2*y[5])/4)+((2*y[2]+2*y[6])/4))/4))/4;
              lz = (2*((3*((2*z[0]+2*z[4])/4)+((2*z[3]+2*z[7])/4))/4)+2*((3*((2*z[1]+2*z[5])/4)+((2*z[2]+2*z[6])/4))/4))/4;
              break;
            case  58:
              lx = (((3*((2*x[0]+2*x[4])/4)+((2*x[3]+2*x[7])/4))/4)+3*((3*((2*x[1]+2*x[5])/4)+((2*x[2]+2*x[6])/4))/4))/4;
              ly = (((3*((2*y[0]+2*y[4])/4)+((2*y[3]+2*y[7])/4))/4)+3*((3*((2*y[1]+2*y[5])/4)+((2*y[2]+2*y[6])/4))/4))/4;
              lz = (((3*((2*z[0]+2*z[4])/4)+((2*z[3]+2*z[7])/4))/4)+3*((3*((2*z[1]+2*z[5])/4)+((2*z[2]+2*z[6])/4))/4))/4;
              break;
            case  59:
              lx = (3*((2*x[1]+2*x[5])/4)+((2*x[2]+2*x[6])/4))/4;
              ly = (3*((2*y[1]+2*y[5])/4)+((2*y[2]+2*y[6])/4))/4;
              lz = (3*((2*z[1]+2*z[5])/4)+((2*z[2]+2*z[6])/4))/4;
              break;
            case  60:
              lx = (2*((2*x[0]+2*x[4])/4)+2*((2*x[3]+2*x[7])/4))/4;
              ly = (2*((2*y[0]+2*y[4])/4)+2*((2*y[3]+2*y[7])/4))/4;
              lz = (2*((2*z[0]+2*z[4])/4)+2*((2*z[3]+2*z[7])/4))/4;
              break;
            case  61:
              lx = (3*((2*((2*x[0]+2*x[4])/4)+2*((2*x[3]+2*x[7])/4))/4)+((2*((2*x[1]+2*x[5])/4)+2*((2*x[2]+2*x[6])/4))/4))/4;
              ly = (3*((2*((2*y[0]+2*y[4])/4)+2*((2*y[3]+2*y[7])/4))/4)+((2*((2*y[1]+2*y[5])/4)+2*((2*y[2]+2*y[6])/4))/4))/4;
              lz = (3*((2*((2*z[0]+2*z[4])/4)+2*((2*z[3]+2*z[7])/4))/4)+((2*((2*z[1]+2*z[5])/4)+2*((2*z[2]+2*z[6])/4))/4))/4;
              break;
            case  62:
              lx = (2*((2*((2*x[0]+2*x[4])/4)+2*((2*x[3]+2*x[7])/4))/4)+2*((2*((2*x[1]+2*x[5])/4)+2*((2*x[2]+2*x[6])/4))/4))/4;
              ly = (2*((2*((2*y[0]+2*y[4])/4)+2*((2*y[3]+2*y[7])/4))/4)+2*((2*((2*y[1]+2*y[5])/4)+2*((2*y[2]+2*y[6])/4))/4))/4;
              lz = (2*((2*((2*z[0]+2*z[4])/4)+2*((2*z[3]+2*z[7])/4))/4)+2*((2*((2*z[1]+2*z[5])/4)+2*((2*z[2]+2*z[6])/4))/4))/4;
              break;
            case  63:
              lx = (((2*((2*x[0]+2*x[4])/4)+2*((2*x[3]+2*x[7])/4))/4)+3*((2*((2*x[1]+2*x[5])/4)+2*((2*x[2]+2*x[6])/4))/4))/4;
              ly = (((2*((2*y[0]+2*y[4])/4)+2*((2*y[3]+2*y[7])/4))/4)+3*((2*((2*y[1]+2*y[5])/4)+2*((2*y[2]+2*y[6])/4))/4))/4;
              lz = (((2*((2*z[0]+2*z[4])/4)+2*((2*z[3]+2*z[7])/4))/4)+3*((2*((2*z[1]+2*z[5])/4)+2*((2*z[2]+2*z[6])/4))/4))/4;
              break;
            case  64:
              lx = (2*((2*x[1]+2*x[5])/4)+2*((2*x[2]+2*x[6])/4))/4;
              ly = (2*((2*y[1]+2*y[5])/4)+2*((2*y[2]+2*y[6])/4))/4;
              lz = (2*((2*z[1]+2*z[5])/4)+2*((2*z[2]+2*z[6])/4))/4;
              break;
            case  65:
              lx = (((2*x[0]+2*x[4])/4)+3*((2*x[3]+2*x[7])/4))/4;
              ly = (((2*y[0]+2*y[4])/4)+3*((2*y[3]+2*y[7])/4))/4;
              lz = (((2*z[0]+2*z[4])/4)+3*((2*z[3]+2*z[7])/4))/4;
              break;
            case  66:
              lx = (3*((((2*x[0]+2*x[4])/4)+3*((2*x[3]+2*x[7])/4))/4)+((((2*x[1]+2*x[5])/4)+3*((2*x[2]+2*x[6])/4))/4))/4;
              ly = (3*((((2*y[0]+2*y[4])/4)+3*((2*y[3]+2*y[7])/4))/4)+((((2*y[1]+2*y[5])/4)+3*((2*y[2]+2*y[6])/4))/4))/4;
              lz = (3*((((2*z[0]+2*z[4])/4)+3*((2*z[3]+2*z[7])/4))/4)+((((2*z[1]+2*z[5])/4)+3*((2*z[2]+2*z[6])/4))/4))/4;
              break;
            case  67:
              lx = (2*((((2*x[0]+2*x[4])/4)+3*((2*x[3]+2*x[7])/4))/4)+2*((((2*x[1]+2*x[5])/4)+3*((2*x[2]+2*x[6])/4))/4))/4;
              ly = (2*((((2*y[0]+2*y[4])/4)+3*((2*y[3]+2*y[7])/4))/4)+2*((((2*y[1]+2*y[5])/4)+3*((2*y[2]+2*y[6])/4))/4))/4;
              lz = (2*((((2*z[0]+2*z[4])/4)+3*((2*z[3]+2*z[7])/4))/4)+2*((((2*z[1]+2*z[5])/4)+3*((2*z[2]+2*z[6])/4))/4))/4;
              break;
            case  68:
              lx = (((((2*x[0]+2*x[4])/4)+3*((2*x[3]+2*x[7])/4))/4)+3*((((2*x[1]+2*x[5])/4)+3*((2*x[2]+2*x[6])/4))/4))/4;
              ly = (((((2*y[0]+2*y[4])/4)+3*((2*y[3]+2*y[7])/4))/4)+3*((((2*y[1]+2*y[5])/4)+3*((2*y[2]+2*y[6])/4))/4))/4;
              lz = (((((2*z[0]+2*z[4])/4)+3*((2*z[3]+2*z[7])/4))/4)+3*((((2*z[1]+2*z[5])/4)+3*((2*z[2]+2*z[6])/4))/4))/4;
              break;
            case  69:
              lx = (((2*x[1]+2*x[5])/4)+3*((2*x[2]+2*x[6])/4))/4;
              ly = (((2*y[1]+2*y[5])/4)+3*((2*y[2]+2*y[6])/4))/4;
              lz = (((2*z[1]+2*z[5])/4)+3*((2*z[2]+2*z[6])/4))/4;
              break;
            case  70:
              lx = (2*x[3]+2*x[7])/4;  ly = (2*y[3]+2*y[7])/4;  lz = (2*z[3]+2*z[7])/4;
              break;
            case  71:
              lx = (3*((2*x[3]+2*x[7])/4)+((2*x[2]+2*x[6])/4))/4;
              ly = (3*((2*y[3]+2*y[7])/4)+((2*y[2]+2*y[6])/4))/4;
              lz = (3*((2*z[3]+2*z[7])/4)+((2*z[2]+2*z[6])/4))/4;
              break;
            case  72:
              lx = (2*((2*x[3]+2*x[7])/4)+2*((2*x[2]+2*x[6])/4))/4;
              ly = (2*((2*y[3]+2*y[7])/4)+2*((2*y[2]+2*y[6])/4))/4;
              lz = (2*((2*z[3]+2*z[7])/4)+2*((2*z[2]+2*z[6])/4))/4;
              break;
            case  73:
              lx = (((2*x[3]+2*x[7])/4)+3*((2*x[2]+2*x[6])/4))/4;
              ly = (((2*y[3]+2*y[7])/4)+3*((2*y[2]+2*y[6])/4))/4;
              lz = (((2*z[3]+2*z[7])/4)+3*((2*z[2]+2*z[6])/4))/4;
              break;
            case  74:
              lx = (2*x[2]+2*x[6])/4;  ly = (2*y[2]+2*y[6])/4;  lz = (2*z[2]+2*z[6])/4;
              break;
            case  75:
              lx = (x[0]+3*x[4])/4;  ly = (y[0]+3*y[4])/4;  lz = (z[0]+3*z[4])/4;
              break;
            case  76:
              lx = (3*((x[0]+3*x[4])/4)+((x[1]+3*x[5])/4))/4;
              ly = (3*((y[0]+3*y[4])/4)+((y[1]+3*y[5])/4))/4;
              lz = (3*((z[0]+3*z[4])/4)+((z[1]+3*z[5])/4))/4;
              break;
            case  77:
              lx = (2*((x[0]+3*x[4])/4)+2*((x[1]+3*x[5])/4))/4;
              ly = (2*((y[0]+3*y[4])/4)+2*((y[1]+3*y[5])/4))/4;
              lz = (2*((z[0]+3*z[4])/4)+2*((z[1]+3*z[5])/4))/4;
              break;
            case  78:
              lx = (((x[0]+3*x[4])/4)+3*((x[1]+3*x[5])/4))/4;
              ly = (((y[0]+3*y[4])/4)+3*((y[1]+3*y[5])/4))/4;
              lz = (((z[0]+3*z[4])/4)+3*((z[1]+3*z[5])/4))/4;
              break;
            case  79:
              lx = (x[1]+3*x[5])/4;  ly = (y[1]+3*y[5])/4;  lz = (z[1]+3*z[5])/4;
              break;
            case  80:
              lx = (3*((x[0]+3*x[4])/4)+((x[3]+3*x[7])/4))/4;
              ly = (3*((y[0]+3*y[4])/4)+((y[3]+3*y[7])/4))/4;
              lz = (3*((z[0]+3*z[4])/4)+((z[3]+3*z[7])/4))/4;
              break;
            case  81:
              lx = (3*((3*((x[0]+3*x[4])/4)+((x[3]+3*x[7])/4))/4)+((3*((x[1]+3*x[5])/4)+((x[2]+3*x[6])/4))/4))/4;
              ly = (3*((3*((y[0]+3*y[4])/4)+((y[3]+3*y[7])/4))/4)+((3*((y[1]+3*y[5])/4)+((y[2]+3*y[6])/4))/4))/4;
              lz = (3*((3*((z[0]+3*z[4])/4)+((z[3]+3*z[7])/4))/4)+((3*((z[1]+3*z[5])/4)+((z[2]+3*z[6])/4))/4))/4;
              break;
            case  82:
              lx = (2*((3*((x[0]+3*x[4])/4)+((x[3]+3*x[7])/4))/4)+2*((3*((x[1]+3*x[5])/4)+((x[2]+3*x[6])/4))/4))/4;
              ly = (2*((3*((y[0]+3*y[4])/4)+((y[3]+3*y[7])/4))/4)+2*((3*((y[1]+3*y[5])/4)+((y[2]+3*y[6])/4))/4))/4;
              lz = (2*((3*((z[0]+3*z[4])/4)+((z[3]+3*z[7])/4))/4)+2*((3*((z[1]+3*z[5])/4)+((z[2]+3*z[6])/4))/4))/4;
              break;
            case  83:
              lx = (((3*((x[0]+3*x[4])/4)+((x[3]+3*x[7])/4))/4)+3*((3*((x[1]+3*x[5])/4)+((x[2]+3*x[6])/4))/4))/4;
              ly = (((3*((y[0]+3*y[4])/4)+((y[3]+3*y[7])/4))/4)+3*((3*((y[1]+3*y[5])/4)+((y[2]+3*y[6])/4))/4))/4;
              lz = (((3*((z[0]+3*z[4])/4)+((z[3]+3*z[7])/4))/4)+3*((3*((z[1]+3*z[5])/4)+((z[2]+3*z[6])/4))/4))/4;
              break;
            case  84:
              lx = (3*((x[1]+3*x[5])/4)+((x[2]+3*x[6])/4))/4;
              ly = (3*((y[1]+3*y[5])/4)+((y[2]+3*y[6])/4))/4;
              lz = (3*((z[1]+3*z[5])/4)+((z[2]+3*z[6])/4))/4;
              break;
            case  85:
              lx = (2*((x[0]+3*x[4])/4)+2*((x[3]+3*x[7])/4))/4;
              ly = (2*((y[0]+3*y[4])/4)+2*((y[3]+3*y[7])/4))/4;
              lz = (2*((z[0]+3*z[4])/4)+2*((z[3]+3*z[7])/4))/4;
              break;
            case  86:
              lx = (3*((2*((x[0]+3*x[4])/4)+2*((x[3]+3*x[7])/4))/4)+((2*((x[1]+3*x[5])/4)+2*((x[2]+3*x[6])/4))/4))/4;
              ly = (3*((2*((y[0]+3*y[4])/4)+2*((y[3]+3*y[7])/4))/4)+((2*((y[1]+3*y[5])/4)+2*((y[2]+3*y[6])/4))/4))/4;
              lz = (3*((2*((z[0]+3*z[4])/4)+2*((z[3]+3*z[7])/4))/4)+((2*((z[1]+3*z[5])/4)+2*((z[2]+3*z[6])/4))/4))/4;
              break;
            case  87:
              lx = (2*((2*((x[0]+3*x[4])/4)+2*((x[3]+3*x[7])/4))/4)+2*((2*((x[1]+3*x[5])/4)+2*((x[2]+3*x[6])/4))/4))/4;
              ly = (2*((2*((y[0]+3*y[4])/4)+2*((y[3]+3*y[7])/4))/4)+2*((2*((y[1]+3*y[5])/4)+2*((y[2]+3*y[6])/4))/4))/4;
              lz = (2*((2*((z[0]+3*z[4])/4)+2*((z[3]+3*z[7])/4))/4)+2*((2*((z[1]+3*z[5])/4)+2*((z[2]+3*z[6])/4))/4))/4;
              break;
            case  88:
              lx = (((2*((x[0]+3*x[4])/4)+2*((x[3]+3*x[7])/4))/4)+3*((2*((x[1]+3*x[5])/4)+2*((x[2]+3*x[6])/4))/4))/4;
              ly = (((2*((y[0]+3*y[4])/4)+2*((y[3]+3*y[7])/4))/4)+3*((2*((y[1]+3*y[5])/4)+2*((y[2]+3*y[6])/4))/4))/4;
              lz = (((2*((z[0]+3*z[4])/4)+2*((z[3]+3*z[7])/4))/4)+3*((2*((z[1]+3*z[5])/4)+2*((z[2]+3*z[6])/4))/4))/4;
              break;
            case  89:
              lx = (2*((x[1]+3*x[5])/4)+2*((x[2]+3*x[6])/4))/4;
              ly = (2*((y[1]+3*y[5])/4)+2*((y[2]+3*y[6])/4))/4;
              lz = (2*((z[1]+3*z[5])/4)+2*((z[2]+3*z[6])/4))/4;
              break;
            case  90:
              lx = (((x[0]+3*x[4])/4)+3*((x[3]+3*x[7])/4))/4;
              ly = (((y[0]+3*y[4])/4)+3*((y[3]+3*y[7])/4))/4;
              lz = (((z[0]+3*z[4])/4)+3*((z[3]+3*z[7])/4))/4;
              break;
            case  91:
              lx = (3*((((x[0]+3*x[4])/4)+3*((x[3]+3*x[7])/4))/4)+((((x[1]+3*x[5])/4)+3*((x[2]+3*x[6])/4))/4))/4;
              ly = (3*((((y[0]+3*y[4])/4)+3*((y[3]+3*y[7])/4))/4)+((((y[1]+3*y[5])/4)+3*((y[2]+3*y[6])/4))/4))/4;
              lz = (3*((((z[0]+3*z[4])/4)+3*((z[3]+3*z[7])/4))/4)+((((z[1]+3*z[5])/4)+3*((z[2]+3*z[6])/4))/4))/4;
              break;
            case  92:
              lx = (2*((((x[0]+3*x[4])/4)+3*((x[3]+3*x[7])/4))/4)+2*((((x[1]+3*x[5])/4)+3*((x[2]+3*x[6])/4))/4))/4;
              ly = (2*((((y[0]+3*y[4])/4)+3*((y[3]+3*y[7])/4))/4)+2*((((y[1]+3*y[5])/4)+3*((y[2]+3*y[6])/4))/4))/4;
              lz = (2*((((z[0]+3*z[4])/4)+3*((z[3]+3*z[7])/4))/4)+2*((((z[1]+3*z[5])/4)+3*((z[2]+3*z[6])/4))/4))/4;
              break;
            case  93:
              lx = (((((x[0]+3*x[4])/4)+3*((x[3]+3*x[7])/4))/4)+3*((((x[1]+3*x[5])/4)+3*((x[2]+3*x[6])/4))/4))/4;
              ly = (((((y[0]+3*y[4])/4)+3*((y[3]+3*y[7])/4))/4)+3*((((y[1]+3*y[5])/4)+3*((y[2]+3*y[6])/4))/4))/4;
              lz = (((((z[0]+3*z[4])/4)+3*((z[3]+3*z[7])/4))/4)+3*((((z[1]+3*z[5])/4)+3*((z[2]+3*z[6])/4))/4))/4;
              break;
            case  94:
              lx = (((x[1]+3*x[5])/4)+3*((x[2]+3*x[6])/4))/4;
              ly = (((y[1]+3*y[5])/4)+3*((y[2]+3*y[6])/4))/4;
              lz = (((z[1]+3*z[5])/4)+3*((z[2]+3*z[6])/4))/4;
              break;
            case  95:
              lx = (x[3]+3*x[7])/4;  ly = (y[3]+3*y[7])/4;  lz = (z[3]+3*z[7])/4;
              break;
            case  96:
              lx = (3*((x[3]+3*x[7])/4)+((x[2]+3*x[6])/4))/4;
              ly = (3*((y[3]+3*y[7])/4)+((y[2]+3*y[6])/4))/4;
              lz = (3*((z[3]+3*z[7])/4)+((z[2]+3*z[6])/4))/4;
              break;
            case  97:
              lx = (2*((x[3]+3*x[7])/4)+2*((x[2]+3*x[6])/4))/4;
              ly = (2*((y[3]+3*y[7])/4)+2*((y[2]+3*y[6])/4))/4;
              lz = (2*((z[3]+3*z[7])/4)+2*((z[2]+3*z[6])/4))/4;
              break;
            case  98:
              lx = (((x[3]+3*x[7])/4)+3*((x[2]+3*x[6])/4))/4;
              ly = (((y[3]+3*y[7])/4)+3*((y[2]+3*y[6])/4))/4;
              lz = (((z[3]+3*z[7])/4)+3*((z[2]+3*z[6])/4))/4;
              break;
            case  99:
              lx = (x[2]+3*x[6])/4;  ly = (y[2]+3*y[6])/4;  lz = (z[2]+3*z[6])/4;
              break;
            case 100:
              lx = x[4];  ly = y[4];  lz = z[4];
              break;
            case 101:
              lx = (3*x[4]+x[5])/4;  ly = (3*y[4]+y[5])/4;  lz = (3*z[4]+z[5])/4;
              break;
            case 102:
              lx = (2*x[4]+2*x[5])/4;  ly = (2*y[4]+2*y[5])/4;  lz = (2*z[4]+2*z[5])/4;
              break;
            case 103:
              lx = (x[4]+3*x[5])/4;  ly = (y[4]+3*y[5])/4;  lz = (z[4]+3*z[5])/4;
              break;
            case 104:
              lx = x[5];  ly = y[5];  lz = z[5];
              break;
            case 105:
              lx = (3*x[4]+x[7])/4;  ly = (3*y[4]+y[7])/4;  lz = (3*z[4]+z[7])/4;
              break;
            case 106:
              lx = (3*((3*x[4]+x[7])/4)+((3*x[5]+x[6])/4))/4;
              ly = (3*((3*y[4]+y[7])/4)+((3*y[5]+y[6])/4))/4;
              lz = (3*((3*z[4]+z[7])/4)+((3*z[5]+z[6])/4))/4;
              break;
            case 107:
              lx = (2*((3*x[4]+x[7])/4)+2*((3*x[5]+x[6])/4))/4;
              ly = (2*((3*y[4]+y[7])/4)+2*((3*y[5]+y[6])/4))/4;
              lz = (2*((3*z[4]+z[7])/4)+2*((3*z[5]+z[6])/4))/4;
              break;
            case 108:
              lx = (((3*x[4]+x[7])/4)+3*((3*x[5]+x[6])/4))/4;
              ly = (((3*y[4]+y[7])/4)+3*((3*y[5]+y[6])/4))/4;
              lz = (((3*z[4]+z[7])/4)+3*((3*z[5]+z[6])/4))/4;
              break;
            case 109:
              lx = (3*x[5]+x[6])/4;  ly = (3*y[5]+y[6])/4;  lz = (3*z[5]+z[6])/4;
              break;
            case 110:
              lx = (2*x[4]+2*x[7])/4;  ly = (2*y[4]+2*y[7])/4;  lz = (2*z[4]+2*z[7])/4;
              break;
            case 111:
              lx = (3*((2*x[4]+2*x[7])/4)+((2*x[5]+2*x[6])/4))/4;
              ly = (3*((2*y[4]+2*y[7])/4)+((2*y[5]+2*y[6])/4))/4;
              lz = (3*((2*z[4]+2*z[7])/4)+((2*z[5]+2*z[6])/4))/4;
              break;
            case 112:
              lx = (2*((2*x[4]+2*x[7])/4)+2*((2*x[5]+2*x[6])/4))/4;
              ly = (2*((2*y[4]+2*y[7])/4)+2*((2*y[5]+2*y[6])/4))/4;
              lz = (2*((2*z[4]+2*z[7])/4)+2*((2*z[5]+2*z[6])/4))/4;
              break;
            case 113:
              lx = (((2*x[4]+2*x[7])/4)+3*((2*x[5]+2*x[6])/4))/4;
              ly = (((2*y[4]+2*y[7])/4)+3*((2*y[5]+2*y[6])/4))/4;
              lz = (((2*z[4]+2*z[7])/4)+3*((2*z[5]+2*z[6])/4))/4;
              break;
            case 114:
              lx = (2*x[5]+2*x[6])/4;  ly = (2*y[5]+2*y[6])/4;  lz = (2*z[5]+2*z[6])/4;
              break;
            case 115:
              lx = (x[4]+3*x[7])/4;  ly = (y[4]+3*y[7])/4;  lz = (z[4]+3*z[7])/4;
              break;
            case 116:
              lx = (3*((x[4]+3*x[7])/4)+((x[5]+3*x[6])/4))/4;
              ly = (3*((y[4]+3*y[7])/4)+((y[5]+3*y[6])/4))/4;
              lz = (3*((z[4]+3*z[7])/4)+((z[5]+3*z[6])/4))/4;
              break;
            case 117:
              lx = (2*((x[4]+3*x[7])/4)+2*((x[5]+3*x[6])/4))/4;
              ly = (2*((y[4]+3*y[7])/4)+2*((y[5]+3*y[6])/4))/4;
              lz = (2*((z[4]+3*z[7])/4)+2*((z[5]+3*z[6])/4))/4;
              break;
            case 118:
              lx = (((x[4]+3*x[7])/4)+3*((x[5]+3*x[6])/4))/4;
              ly = (((y[4]+3*y[7])/4)+3*((y[5]+3*y[6])/4))/4;
              lz = (((z[4]+3*z[7])/4)+3*((z[5]+3*z[6])/4))/4;
              break;
            case 119:
              lx = (x[5]+3*x[6])/4;  ly = (y[5]+3*y[6])/4;  lz = (z[5]+3*z[6])/4;
              break;
            case 120:
              lx = x[7];  ly = y[7];  lz = z[7];
              break;
            case 121:
              lx = (3*x[7]+x[6])/4;  ly = (3*y[7]+y[6])/4;  lz = (3*z[7]+z[6])/4;
              break;
            case 122:
              lx = (2*x[7]+2*x[6])/4;  ly = (2*y[7]+2*y[6])/4;  lz = (2*z[7]+2*z[6])/4;
              break;
            case 123:
	      lx = (3*x[6]+x[7])/4;  ly =(3*y[6]+y[7])/4;  lz = (3*z[6]+z[7])/4;
              break;
            case 124:
              lx = x[6];  ly = y[6];  lz = z[6];
              break;
          }
	  break;
        }
      }
      gdof = dof[j];
      // peridoic b.c. are set only from one side
      if ((fabs(lx-2*Pi)<1e-6)||(fabs(ly-2*Pi/3)<1e-6))
        continue;
      //OutPut(i << " loc " << j << " glob " << gdof << " x " << lx << " y " <<  ly
      //<< " z " << lz << endl);
      // check if the entries in the vectors of coordinates are already filled
      // and compare to computed coordinates
      // (assume that the vectors of coordinates are initialized with -4711)
      if (x_dof[gdof] != -4711)
      {
        if((fabs(lx - x_dof[gdof])>1e-5)||(fabs(ly - y_dof[gdof])>1e-5)||(fabs(lz - z_dof[gdof])>1e-5))
        {
          OutPut("Error in GetCoordinatesOfDof "<< j << " " << gdof << endl);
          OutPut("x " << x_dof[gdof] << " " << lx << endl);
          OutPut("y " << y_dof[gdof] << " " << ly << endl);
          OutPut("z " << z_dof[gdof] << " " << lz << endl);
          exit(1);
        }
      }
      // fill the vectors of the coordinates
      x_dof[gdof] = lx;
      y_dof[gdof] = ly;
      z_dof[gdof] = lz;

      // check if lz is already in the list and compute N_z_layers_help
      for (n=0;n<max_z_layers;n++)
      {
        if(fabs(coord_z_layers_help[n]-lz)<1e-6)
        {
          break;
        }
        if(coord_z_layers_help[n] == -4711.)
        {
          coord_z_layers_help[n] = lz;
          N_z_layers_help++;
          break;
        }
      }
    }
  }
  //  OutPut("N_z_layers_help : "<<N_z_layers_help<<endl);

  // set N_z_layers
  *N_z_layers = N_z_layers_help;

  // allocate array with z coordinates
  coord_z_layers = new double[N_z_layers_help];

  // order z-coordinates and fill coord_z_layers
  while(change!=1)
  {
    change = 1;

    for (k=0;k<N_z_layers_help-1;k++)
    {
      if(coord_z_layers_help[k+1]<coord_z_layers_help[k])
      {
        help = coord_z_layers_help[k];
        coord_z_layers_help[k] = coord_z_layers_help[k+1];
        coord_z_layers_help[k+1] = help;

        change = 0;
      }

    }
  }

  // fill coord_z_layers
  for(m=0;m<N_z_layers_help;m++)
  {
    coord_z_layers[m] = coord_z_layers_help[m];
    //OutPut(m << " " <<  coord_z_layers[m] << endl);
  }

  delete coord_z_layers_help;
  //exit(1);
}


void  ComputeMeanVelocity(int N_z_layers,double *coord_z_layers,
double *mean_velocity, double *rms_velocity,
TFEFunction3D *U1, double *x_dof,
double *y_dof, double *z_dof)
{
  double T = TDatabase::TimeDB->CURRENTTIME;
  double dt =  TDatabase::TimeDB->TIMESTEPLENGTH;

  double N_U, *u1, *sum_values_per_layer, bulk_velocity, *rms_values_per_layer;
  int *number_of_sums, i, j;

  N_U = U1->GetLength();
  u1 = U1->GetValues();
  OutPut(N_U << endl);

  // initialisation(with 0) of the vectors for the
  // average determination

  sum_values_per_layer = new double[N_z_layers];
  memset(sum_values_per_layer,0,N_z_layers*SizeOfDouble);

  number_of_sums = new int[N_z_layers];
  memset(number_of_sums,0,N_z_layers*SizeOfInt);

  rms_values_per_layer = new double[N_z_layers];
  memset(rms_values_per_layer,0,N_z_layers*SizeOfDouble);

  // summation of all values per layer
  for (i=0;i<N_U;i++)
  {
    for (j=0;j<N_z_layers;j++)
    {
      if (fabs(z_dof[i]-coord_z_layers[j])<1e-6)
      {
        sum_values_per_layer[j] = sum_values_per_layer[j] + u1[i];
        number_of_sums[j]++;
        break;
      }
    }
  }

  // determination of the average
  for (i=0;i<N_z_layers;i++)
  {
    sum_values_per_layer[i] = sum_values_per_layer[i]/number_of_sums[i];
    //OutPut("av " << number_of_sums[i] << " " << sum_values_per_layer[i] << endl);
  }

  // summation over all values for the calculation of rms_values_per_layer
  for (i=0;i<N_U;i++)
  {
      for (j=0;j<N_z_layers;j++)
      {
	  if (fabs(z_dof[i]-coord_z_layers[j])<1e-6)
	  {
	      rms_values_per_layer[j] += (u1[i]-sum_values_per_layer[j])*(u1[i]-sum_values_per_layer[j]);
	      break;
	  }
      }
  }

  // calculation of rms_velocity
  for (i=0;i<N_z_layers;i++)
  {
      rms_values_per_layer[i] = sqrt(rms_values_per_layer[i]/N_z_layers); 
  }

  if (T>0)
  {
      for (i=0;i<N_z_layers;i++)
      {
	  mean_velocity[i] = mean_velocity[i] + (dt/T)*(sum_values_per_layer[i]-mean_velocity[i]);
	  rms_velocity[i] = rms_velocity[i]+(dt/T)*(rms_values_per_layer[i] - rms_velocity[i] );
	  OutPut("mean velo : time " << T << " " << setw(10) << coord_z_layers[i] << " " << 
		 setw(10) <<  180*(1-fabs(1-coord_z_layers[i])) << " mean  " <<
		 setw(10) <<  mean_velocity[i] << " rms " << setw(10) <<  rms_velocity[i] << endl);
      }
  }
  else
  {
      for (i=0;i<N_z_layers;i++)
      {
	  mean_velocity[i] = sum_values_per_layer[i];
	  rms_velocity[i] = rms_values_per_layer[i];
	  OutPut("mean velo : time " << T << " " << setw(10) << coord_z_layers[i] << " " << 
		 setw(10) <<  180*(1-fabs(1-coord_z_layers[i])) << " " <<
		 setw(10) <<  mean_velocity[i] << endl);
      }
  }   
  // determination of the bulk velocity
  bulk_velocity = 0.0;
  for (i=0;i<N_z_layers-1;i++) 
    {
     bulk_velocity += 0.5*( mean_velocity[i+1]+mean_velocity[i])*(coord_z_layers[i+1]-coord_z_layers[i]);      
    }
  bulk_velocity /= 2.0;
  
  OutPut("bulk velo : time " << T << " " << bulk_velocity << endl);
  TDatabase::ParamDB->INTERNAL_BULK_AFTER = bulk_velocity;
  delete rms_values_per_layer;
  delete sum_values_per_layer;
  delete number_of_sums;

  //    OutPut("mean " << N_U << " " << u1[0] << endl);
}
