// Navier-Stokes problem
// channel flow in 3D
//

#include <PeriodicJoint.h>

#define U_INFTY 1   
#define __CHANNEL_TOBIAS__

// ========================================================================
// example file
// ========================================================================

void ExampleFile()
{
  OutPut("Example: Channel_Tobias.nonzero_rhs.h, tau_omega (P8) = " << 
         TDatabase::ParamDB->P8 << endl);
}
// ========================================================================
// exact solution
// ========================================================================
void InitialU1(double x, double y, double z, double *values)
{
   values[0] = 0;
}

void InitialU2(double x, double y, double z, double *values)
{
   values[0] = 0;
}

void InitialU3(double x, double y, double z, double *values)
{
   values[0] = 0; 
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
  double *coeff;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    coeff[0] = eps;
    coeff[1] = 0;
    coeff[2] = 0;
    coeff[3] = TDatabase::ParamDB->P8/Pi;
  }
}

// periodic b.c. :
//  x = 0 and x = 2 Pi
//  z = 0 and z = Pi
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
  int MaxLen, MaxLen1, x_vert_m1, x_vert_1, z_vert_0, z_vert_2;
  int found, x1_vert_m1, x1_vert_1, z1_vert_0, z1_vert_2;

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
      z_vert_0 =  z_vert_2 = x_vert_m1 = x_vert_1 = 0;
      // compute coordinates of vertices
      for (l1=0;l1<TmpLen[j];l1++)
      {
        cell->GetVertex(TmpFV[j*MaxLen+l1])->GetCoords(X[l1], Y[l1], Z[l1]);
        // check if vertex on z = 0
        if (fabs(Z[l1])<1e-5)
          z_vert_0++;
        // check if vertex on z = Pi
        if (fabs(Z[l1]-Pi)<1e-5)
          z_vert_2++;
        // check if vertex on x = 0
        if (fabs(X[l1])<1e-5)
          x_vert_m1++;
        // check if vertex on x = 2 Pi        
        if (fabs(X[l1]-2*Pi)<1e-5)
          x_vert_1++;
      }
      found = 0;
      if (z_vert_0==TmpLen[j])
        found++;
      if (z_vert_2==TmpLen[j])
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
      OutPut("bary " << x << " " << y << " " << z << endl);
      for (l1=0;l1<TmpLen[j];l1++)
       {
        OutPut("face " << X[l1] << " " << Y[l1] << " " << Z[l1] << endl);
      }
       // inner loop over the cells
      //for(ii=i+1;ii<N_Cells;ii++)
      for(ii=i;ii<N_Cells;ii++)
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
          z1_vert_0 =  z1_vert_2 = x1_vert_m1 = x1_vert_1 = 0;
          // compute coordinates of vertices
          for (l1=0;l1<TmpLen[jj];l1++)
          {
            cell1->GetVertex(TmpFV[jj*MaxLen1+l1])->
              GetCoords(X1[l1], Y1[l1], Z1[l1]);
            // check if vertex on z = 0
            if (fabs(Z1[l1])<1e-5)
              z1_vert_0++;
            // check if vertex on z = Pi
            if (fabs(Z1[l1]-Pi)<1e-5)
              z1_vert_2++;
            // check if vertex on x = 0
            if (fabs(X1[l1])<1e-5)
              x1_vert_m1++;
            // check if vertex on x = 2 Pi        
            if (fabs(X1[l1]-2*Pi)<1e-5)
              x1_vert_1++;
            x1+=X1[l1];
            y1+=Y1[l1];
            z1+=Z1[l1];
          }
          x1 /= TmpLen[jj];
          y1 /= TmpLen[jj];
          z1 /= TmpLen[jj];
          found = 0;
          //OutPut("baryopp " << x1 << " " << y1 << " " << z1 << endl);
          if ((z_vert_0==TmpLen[j])&&(z1_vert_2==TmpLen[jj]))
          {
            if  ((fabs(x-x1)>1e-5)||(fabs(y-y1)>1e-5))
              continue;
            found++;
          }
          if ((z_vert_2==TmpLen[j])&&(z1_vert_0==TmpLen[jj]))
          {
            if  ((fabs(x-x1)>1e-5)||(fabs(y-y1)>1e-5))
              continue;
            found++;
          }
          if ((x_vert_m1==TmpLen[j])&&(x1_vert_1==TmpLen[jj]))
          {
            if  ((fabs(z-z1)>1e-5)||(fabs(y-y1)>1e-5))
              continue;
             found++;
         }
          if ((x_vert_1==TmpLen[j])&&(x1_vert_m1==TmpLen[jj]))
          {
            if  ((fabs(z-z1)>1e-5)||(fabs(y-y1)>1e-5))
              continue;
            found++;
          }
          if (!found)
            continue;
           for (l1=0;l1<TmpLen[jj];l1++)
           {
            OutPut("opp " << X1[l1] << " " << Y1[l1] << " " << Z1[l1] << endl);
           }
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

