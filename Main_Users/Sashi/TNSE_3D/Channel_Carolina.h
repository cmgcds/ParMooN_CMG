// Navier-Stokes problem
// channel flow in 3D
//

#include <PeriodicJoint.h>

#define U_INFTY 1   
#define __CHANNEL_CAROLINA__
#define __CHANNEL_TOBIAS__

// ========================================================================
// example file
// ========================================================================

void ExampleFile()
{
  OutPut("Example: Channel_Carolina.h with velocity of top wall " <<  
	 TDatabase::ParamDB->P9 << 
	 " NOISE " << TDatabase::ParamDB->P7 << endl);
  TDatabase::ParamDB->INTERNAL_PROBLEM_IDENTITY = 6;
}
// ========================================================================
// exact solution
// ========================================================================
void InitialU1(double x, double y, double z, double *values)
{
// interpolation of the boundary condition + noise
    values[0] = TDatabase::ParamDB->P9*sqrt(z)*(1+TDatabase::ParamDB->P7*sin(10*Pi*z)); 
    values[0] = TDatabase::ParamDB->P7*sin(10*Pi*z); 
}

void InitialU2(double x, double y, double z, double *values)
{
    values[0] = TDatabase::ParamDB->P7*sin(8*Pi*z);
}

void InitialU3(double x, double y, double z, double *values)
{
   values[0] = TDatabase::ParamDB->P7*sin(14*Pi*z); 
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
    if (fabs(z-1) < 1e-5)
	value = TDatabase::ParamDB->P9;
    else
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
  double *coeff;

  TDatabase::ParamDB->RE_NR = 1.0;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    coeff[0] = 1;
    coeff[1] = 0;
    coeff[2] = 0;
    coeff[3] = 0;
  }
}

// periodic b.c. : This is the x-y version
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
  double y_bottom = 0, y_top = 1, x_bottom = 0, x_top = 1; 

  N_Cells = Coll->GetN_Cells();

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

void SetZCoordinates(TCollection *Coll, int level)
{
  int i, j, k, N_Cells, N_Layers, N_V, layer_int, grid_type;
  TBaseCell *cell;
  TVertex *vertex;
  double x,y,z, layer, omega = 4.5;

  N_Cells = Coll->GetN_Cells();
  // number of layers on initial grid * 2^level
  N_Layers = TDatabase::ParamDB->N_CELL_LAYERS*(int)(pow(2.0,level));
  grid_type = TDatabase::ParamDB->GRID_TYPE;
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
      layer = (atanh(tanh(omega)*z)/omega)*N_Layers;
      layer_int = (int)(layer+1e-7);
      //OutPut(layer << " ");
      // not the correct position
      if (fabs(layer-layer_int)>1e-5)
      {
        //OutPut(fabs(layer-layer_int) << " ");
        // find closest odd integer
        if (fabs(layer_int/2.0 - (int)(layer_int/2.0+1e-7))>0.1)
          k = layer_int;
        else
          k = layer_int+1;
        //OutPut(" k " << k << " ");
        // compute new z coordinate
	z = tanh(omega*k/(N_Layers))/tanh(omega);
        // set new z coordinate
        vertex->SetCoords(x,y,z);
      }
    }
  }
  OutPut(endl);
}


void CheckZCoordinates(TCollection *Coll, int level)
{
  int i, j, k, N_Cells, N_Layers, N_V, found, grid_type;
  TBaseCell *cell;
  TVertex *vertex;
  double x,y,z, layer, *zcoor, omega = 4.5;

  N_Cells = Coll->GetN_Cells();
  grid_type = TDatabase::ParamDB->GRID_TYPE;
  // number of layers on initial grid * 2^level
  N_Layers = TDatabase::ParamDB->N_CELL_LAYERS*(int)(pow(2.0,level));
  zcoor = new double[N_Layers+1];
  for (i=0;i<=N_Layers;i++)
  {
      zcoor[i] =  tanh(omega*i/(N_Layers))/tanh(omega);
      OutPut( zcoor[i] << endl);
  }
  
  OutPut(1-zcoor[N_Layers-1] << ", first line condition : " << 1.0/(5*TDatabase::ParamDB->P9) << endl ); 

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

