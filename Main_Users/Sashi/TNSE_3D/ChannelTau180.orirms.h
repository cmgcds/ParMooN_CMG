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
double DNS_profile_180(double zz)
{
  int i;
  double value, val;

  double z[65] =
  {
    0.0000e-00, 3.0118e-04, 1.2045e-03, 2.7095e-03, 4.8153e-03, 7.5205e-03, 1.0823e-02, 1.4722e-02, 1.9215e-02, 2.4298e-02, 2.9969e-02, 3.6224e-02,
    4.3060e-02, 5.0472e-02, 5.8456e-02, 6.7007e-02, 7.6120e-02, 8.5790e-02, 9.6011e-02, 1.0678e-01, 1.1808e-01, 1.2991e-01, 1.4227e-01, 1.5515e-01,
    1.6853e-01, 1.8242e-01, 1.9679e-01, 2.1165e-01, 2.2699e-01, 2.4279e-01, 2.5905e-01, 2.7575e-01, 2.9289e-01, 3.1046e-01, 3.2844e-01, 3.4683e-01,
    3.6561e-01, 3.8477e-01, 4.0430e-01, 4.2419e-01, 4.4443e-01, 4.6500e-01, 4.8590e-01, 5.0710e-01, 5.2860e-01, 5.5039e-01, 5.7244e-01, 5.9476e-01,
    6.1732e-01, 6.4010e-01, 6.6311e-01, 6.8632e-01, 7.0972e-01, 7.3329e-01, 7.5702e-01, 7.8090e-01, 8.0491e-01, 8.2904e-01, 8.5327e-01, 8.7759e-01,
    9.0198e-01, 9.2644e-01, 9.5093e-01, 9.7546e-01, 1.0000e-00
  };

  double  Umean[65] =
  {
    0.0000e+00, 5.3639e-02, 2.1443e-01, 4.8197e-01, 8.5555e-01, 1.3339e+00, 1.9148e+00, 2.5939e+00, 3.3632e+00, 4.2095e+00, 5.1133e+00,
    6.0493e+00, 6.9892e+00, 7.9052e+00, 8.7741e+00, 9.5790e+00, 1.0311e+01, 1.0967e+01, 1.1550e+01, 1.2066e+01, 1.2520e+01, 1.2921e+01,
    1.3276e+01, 1.3590e+01, 1.3870e+01, 1.4121e+01, 1.4349e+01, 1.4557e+01, 1.4750e+01, 1.4931e+01, 1.5101e+01, 1.5264e+01, 1.5419e+01,
    1.5569e+01, 1.5714e+01, 1.5855e+01, 1.5993e+01, 1.6128e+01, 1.6260e+01, 1.6389e+01, 1.6515e+01, 1.6637e+01, 1.6756e+01, 1.6872e+01,
    1.6985e+01, 1.7094e+01, 1.7200e+01, 1.7302e+01, 1.7400e+01, 1.7494e+01, 1.7585e+01, 1.7672e+01, 1.7756e+01, 1.7835e+01, 1.7911e+01,
    1.7981e+01, 1.8045e+01, 1.8103e+01, 1.8154e+01, 1.8198e+01, 1.8235e+01, 1.8264e+01, 1.8285e+01, 1.8297e+01, 1.8301e+01
  };

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

// ========================================================================
// exact solution
// ========================================================================
double DNS_profile_395(double zz)
{
  int i;
  double value, val;

  double z[129] =
  {
    0.0000e-00,  7.5298e-05,  3.0118e-04,  6.7762e-04,  1.2045e-03,  1.8819e-03,  2.7095e-03,  3.6874e-03,  4.8153e-03,  6.0930e-03,  7.5205e-03,  9.0974e-03,  1.0823e-02,
    1.2699e-02,  1.4722e-02,  1.6895e-02,  1.9215e-02,  2.1683e-02,  2.4298e-02,  2.7060e-02,  2.9969e-02,  3.3024e-02,  3.6224e-02,  3.9569e-02,  4.3060e-02,  4.6694e-02,
    5.0472e-02,  5.4393e-02,  5.8456e-02,  6.2661e-02,  6.7007e-02,  7.1494e-02,  7.6120e-02,  8.0886e-02,  8.5790e-02,  9.0832e-02,  9.6011e-02,  1.0133e-01,  1.0678e-01,
    1.1236e-01,  1.1808e-01,  1.2393e-01,  1.2991e-01,  1.3603e-01,  1.4227e-01,  1.4864e-01,  1.5515e-01,  1.6178e-01,  1.6853e-01,  1.7541e-01,  1.8242e-01,  1.8954e-01,
    1.9679e-01,  2.0416e-01,  2.1165e-01,  2.1926e-01,  2.2699e-01,  2.3483e-01,  2.4279e-01,  2.5086e-01,  2.5905e-01,  2.6735e-01,  2.7575e-01,  2.8427e-01,  2.9289e-01,
    3.0162e-01,  3.1046e-01,  3.1940e-01,  3.2844e-01,  3.3758e-01,  3.4683e-01,  3.5617e-01,  3.6561e-01,  3.7514e-01,  3.8477e-01,  3.9449e-01,  4.0430e-01,  4.1420e-01,
    4.2419e-01,  4.3427e-01,  4.4443e-01,  4.5467e-01,  4.6500e-01,  4.7541e-01,  4.8590e-01,  4.9646e-01,  5.0710e-01,  5.1782e-01,  5.2860e-01,  5.3946e-01,  5.5039e-01,
    5.6138e-01,  5.7244e-01,  5.8357e-01,  5.9476e-01,  6.0601e-01,  6.1732e-01,  6.2868e-01,  6.4010e-01,  6.5158e-01,  6.6311e-01,  6.7469e-01,  6.8632e-01,  6.9799e-01,
    7.0972e-01,  7.2148e-01,  7.3329e-01,  7.4513e-01,  7.5702e-01,  7.6894e-01,  7.8090e-01,  7.9289e-01,  8.0491e-01,  8.1696e-01,  8.2904e-01,  8.4114e-01,  8.5327e-01,
    8.6542e-01,  8.7759e-01,  8.8978e-01,  9.0198e-01,  9.1420e-01,  9.2644e-01,  9.3868e-01,  9.5093e-01,  9.6319e-01,  9.7546e-01,  9.8773e-01,  1.0000e-00
  };

  double  Umean[129] =
  {
    0.0000e+00,   2.9538e-02,  1.1811e-01,   2.6562e-01,   4.7198e-01,   7.3701e-01,   1.0605e+00,  1.4417e+00,  1.8799e+00,  2.3729e+00,  2.9177e+00,  3.5093e+00,
    4.1409e+00,  4.8032e+00,  5.4854e+00,  6.1754e+00,  6.8611e+00,  7.5309e+00,  8.1754e+00,  8.7870e+00,  9.3607e+00,  9.8937e+00,  1.0385e+01,  1.0836e+01,
    1.1248e+01,  1.1624e+01,  1.1966e+01,  1.2278e+01,  1.2563e+01,  1.2822e+01,  1.3060e+01,  1.3278e+01,  1.3479e+01,  1.3664e+01,  1.3837e+01,  1.3998e+01,
    1.4148e+01,  1.4290e+01,  1.4425e+01,  1.4552e+01,  1.4673e+01,  1.4790e+01,  1.4902e+01,  1.5011e+01,  1.5117e+01,  1.5221e+01,  1.5322e+01,  1.5421e+01,
    1.5518e+01,  1.5614e+01,  1.5707e+01,  1.5799e+01,  1.5890e+01,  1.5979e+01,  1.6067e+01,  1.6153e+01,  1.6239e+01,  1.6324e+01,  1.6409e+01,  1.6493e+01,
    1.6576e+01,  1.6659e+01,  1.6741e+01,  1.6823e+01,  1.6903e+01,  1.6984e+01,  1.7063e+01,  1.7141e+01,  1.7218e+01,  1.7294e+01,  1.7369e+01,  1.7443e+01,
    1.7517e+01,  1.7590e+01,  1.7664e+01,  1.7738e+01,  1.7812e+01,  1.7886e+01,  1.7960e+01,  1.8034e+01,  1.8108e+01,  1.8182e+01,  1.8254e+01,  1.8326e+01,
    1.8396e+01,  1.8466e+01,  1.8535e+01,  1.8603e+01,  1.8669e+01,  1.8734e+01,  1.8797e+01,  1.8859e+01,  1.8919e+01,  1.8978e+01,  1.9035e+01,  1.9092e+01,
    1.9148e+01,  1.9202e+01,  1.9256e+01,  1.9308e+01,  1.9359e+01,  1.9408e+01,  1.9456e+01,  1.9503e+01,  1.9548e+01,  1.9593e+01,  1.9636e+01,  1.9678e+01,
    1.9719e+01,  1.9758e+01,  1.9796e+01,  1.9832e+01,  1.9865e+01,  1.9897e+01,  1.9927e+01,  1.9955e+01,  1.9981e+01,  2.0004e+01,  2.0026e+01,  2.0046e+01,
    2.0064e+01,  2.0080e+01,  2.0094e+01,  2.0106e+01,  2.0116e+01,  2.0123e+01,  2.0129e+01,  2.0132e+01,  2.0133e+01
  };

  if (zz<=1)
  {
    val = zz;
  }
  else
  {
    val = 2-zz;
  }

  for(i=0;i<128;i++)
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
  double RE=TDatabase::ParamDB->RE_NR;

  if (RE==180)
  {
    TDatabase::ParamDB->INTERNAL_BULK_MEAN = 15.6803;
    TDatabase::ParamDB->INTERNAL_BULK_SIMULATION = 15.6803;
    values[0] = DNS_profile_180(z)
      +noise*TDatabase::ParamDB->INTERNAL_BULK_SIMULATION*(2*(double)rand()/RAND_MAX-1);
  }
  else
  {
    if (RE==395)
    {
	TDatabase::ParamDB->INTERNAL_BULK_MEAN = 17.5452;
	TDatabase::ParamDB->INTERNAL_BULK_SIMULATION = 17.5452;
	values[0] = DNS_profile_395(z)
	    +noise*TDatabase::ParamDB->INTERNAL_BULK_SIMULATION*(2*(double)rand()/RAND_MAX-1);
    }
    else
    {
	OutPut("wrong reynolds number " << endl);
	exit(4711);
    }
  }
  if ((fabs(z)<1e-6)||(fabs(2-z)<1e-6))
      values[0] = 0;
  // values[0] = scale*z*(2-z)+noise*2*scale*(2*(double)rand()/RAND_MAX-1)/3;
}

void InitialU2(double x, double y, double z, double *values)
{
  double  noise = 0.1;
 
  values[0] = noise*TDatabase::ParamDB->INTERNAL_BULK_SIMULATION*(2*(double)rand()/RAND_MAX-1);
  if ((fabs(z)<1e-6)||(fabs(2-z)<1e-6))
      values[0] = 0;
}

void InitialU3(double x, double y, double z, double *values)
{
  double  noise = 0.1;

  values[0] = noise*TDatabase::ParamDB->INTERNAL_BULK_SIMULATION*(2*(double)rand()/RAND_MAX-1);
  if ((fabs(z)<1e-6)||(fabs(2-z)<1e-6))
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
  double *coeff, u1, u2, dt;

  u1 = TDatabase::ParamDB->INTERNAL_BULK_MEAN;
  u2 = TDatabase::ParamDB->INTERNAL_BULK_SIMULATION;
  dt = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    coeff[0] = eps;
    coeff[1] = 1 + (u1 - u2)/dt ;
    //coeff[1] = 1;
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
  double RE=TDatabase::ParamDB->RE_NR;
  double y_bottom, y_top, x_bottom, x_top;

  if (RE==180)
  {
       y_bottom = -2*Pi/3;
       y_top = 2*Pi/3; 
       x_bottom = -2*Pi; 
       x_top = 2*Pi;
  }
  else
  {
       y_bottom = -Pi/2;
       y_top = Pi/2; 
       x_bottom = -Pi; 
       x_top = Pi;
  }

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
  int i, j, k, N_Cells, N_Layers, N_V, layer_int, grid_type;
  TBaseCell *cell;
  TVertex *vertex;
  double x,y,z, layer;

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
      if (grid_type == 0)
	  layer = (atanh(tanh(2.75)*(z-1))/2.75+1)*N_Layers/2.0;
      else
	  layer = acos(1-z)*N_Layers/Pi;
      layer_int = (int)(layer+1e-7);
      //   OutPut(layer << " ");
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
	if (grid_type == 0)
	    z = 1 + tanh(2.75*(2.0*k/N_Layers -1))/tanh(2.75);
	else
	    z = 1-cos(Pi*k/N_Layers);
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
  double x,y,z, layer, *zcoor;

  N_Cells = Coll->GetN_Cells();
  grid_type = TDatabase::ParamDB->GRID_TYPE;
  // number of layers on initial grid * 2^level
  N_Layers = TDatabase::ParamDB->N_CELL_LAYERS*(int)(pow(2.0,level));
  zcoor = new double[N_Layers+1];
  for (i=0;i<=N_Layers;i++)
  {
      if (grid_type==0)
	  zcoor[i] = 1 + tanh(2.75*(2.0*i/N_Layers -1))/tanh(2.75);
      else
	  zcoor[i] = 1-cos(Pi*i/N_Layers);
  }
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
  double x[8],y[8],z[8], *coord_z_layers_help, per_x, per_y;
  double RE=TDatabase::ParamDB->RE_NR;
 TBaseCell *cell;
  TVertex *vertex;
  FE3D CurrentElement;
  double lx,ly,lz,help;
  int gdof;

  if (RE==180)
  { 
      per_x = 2*Pi;
      per_y = 2*Pi/3;
  }
  else
  {
      per_x = Pi;
      per_y = Pi/2;
  }

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
      if ((fabs(lx-per_x)<1e-6)||(fabs(ly-per_y)<1e-6))
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

/*******************************************************************************/
/*                                                                             */
/* COMPUTATIONS OF MEAN VALUES FOR COMPARISON WITH [MKM99]                     */
/*                                                                             */
/*******************************************************************************/

void  ComputeMeanVelocity(int N_z_layers,   
			  double *coord_z_layers,   
			  double *rms_velocity,
			  double *mean_velocity, 
			  double *mean_velocity_u2, 
			  double *mean_velocity_u3,
			  double *dmean_velocity, 
			  double *R_xx, 
			  double *R_xy, 
			  double *R_xz, 
			  double *R_yy, 
			  double *R_yz, 
			  double *R_zz,
			  TFEFunction3D *U1,  
			  TFEFunction3D *U2,  
			  TFEFunction3D *U3,
			  double *x_dof,  
			  double *y_dof,  
			  double *z_dof)
{
  double T  = TDatabase::TimeDB->CURRENTTIME;
  double dt = TDatabase::TimeDB->TIMESTEPLENGTH;
  double RE = TDatabase::ParamDB->RE_NR;

  double N_U, *u1, *u2, *u3, *sum_values_per_layer, bulk_velocity, *rms_values_per_layer;
  double a, b, D;
  double *d_mean_velocity, *sum_values_per_layer_u2, *sum_values_per_layer_u3;
  double *sum_values_per_layer_xx,  *sum_values_per_layer_yy,  *sum_values_per_layer_zz;
  double *sum_values_per_layer_xy,  *sum_values_per_layer_xz,  *sum_values_per_layer_yz;
  int *number_of_sums, i, j;

  // length of one component of the solution vector
  N_U = U1->GetLength();

  // get the components of the solution
  u1 = U1->GetValues();
  u2 = U2->GetValues();
  u3 = U3->GetValues();

  // initialisation(with 0) of the vectors for the  average determination
  sum_values_per_layer = new double[11*N_z_layers];
  memset(sum_values_per_layer,0,11*N_z_layers*SizeOfDouble);

  // set pointers for the individual components 
  sum_values_per_layer_u2 = sum_values_per_layer + N_z_layers;
  sum_values_per_layer_u3 = sum_values_per_layer_u2 +  N_z_layers;
  rms_values_per_layer = sum_values_per_layer_u3 + N_z_layers;
  d_mean_velocity = rms_values_per_layer + N_z_layers;
  sum_values_per_layer_xx = d_mean_velocity + N_z_layers;
  sum_values_per_layer_xy = sum_values_per_layer_xx + N_z_layers;
  sum_values_per_layer_xz = sum_values_per_layer_xy + N_z_layers;
  sum_values_per_layer_yy = sum_values_per_layer_xz + N_z_layers;
  sum_values_per_layer_yz = sum_values_per_layer_yy + N_z_layers;
  sum_values_per_layer_zz = sum_values_per_layer_yz + N_z_layers;

  number_of_sums = new int[N_z_layers];
  memset(number_of_sums,0,N_z_layers*SizeOfInt);

  // summation of all values per layer
  // loop over all velo d.o.f.
  for (i=0;i<N_U;i++)
  {
      // find the layer for the velo d.o.f
    for (j=0;j<N_z_layers;j++)
    {
	// check if coordinate of d.o.f. 
      if (fabs(z_dof[i]-coord_z_layers[j])<1e-6)
      {
	  //u1, u2, u3
        sum_values_per_layer[j] = sum_values_per_layer[j] + u1[i];
        sum_values_per_layer_u2[j] = sum_values_per_layer_u2[j] + u2[i];
        sum_values_per_layer_u3[j] = sum_values_per_layer_u3[j] + u3[i];

	// ui*uj , i=1,2,3
        sum_values_per_layer_xx[j] = sum_values_per_layer_xx[j] + (u1[i]*u1[i]);
        sum_values_per_layer_yy[j] = sum_values_per_layer_yy[j] + (u2[i]*u2[i]);
        sum_values_per_layer_zz[j] = sum_values_per_layer_zz[j] + (u3[i]*u3[i]);
        sum_values_per_layer_xy[j] = sum_values_per_layer_xy[j] + (u1[i]*u2[i]);
        sum_values_per_layer_xz[j] = sum_values_per_layer_xz[j] + (u1[i]*u3[i]);
        sum_values_per_layer_yz[j] = sum_values_per_layer_yz[j] + (u3[i]*u2[i]);

        number_of_sums[j]++;

        break;
      }
    }
  }

  // determination of the average
  for (i=0;i<N_z_layers;i++)
  {
    sum_values_per_layer[i] = sum_values_per_layer[i]/number_of_sums[i];
    sum_values_per_layer_u2[i] = sum_values_per_layer_u2[i]/number_of_sums[i];
    sum_values_per_layer_u3[i] = sum_values_per_layer_u3[i]/number_of_sums[i];
    sum_values_per_layer_xx[i] = sum_values_per_layer_xx[i]/number_of_sums[i];
    sum_values_per_layer_yy[i] = sum_values_per_layer_yy[i]/number_of_sums[i];
    sum_values_per_layer_zz[i] = sum_values_per_layer_zz[i]/number_of_sums[i];
    sum_values_per_layer_xy[i] = sum_values_per_layer_xy[i]/number_of_sums[i];
    sum_values_per_layer_xz[i] = sum_values_per_layer_xz[i]/number_of_sums[i];
    sum_values_per_layer_yz[i] = sum_values_per_layer_yz[i]/number_of_sums[i];
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
  
  // compute rms values
  for (i=0;i<N_z_layers;i++)
  {
      rms_values_per_layer[i] = sqrt(rms_values_per_layer[i]/number_of_sums[i]);
  }

  // update the mean values
  if (T>0)
  {
    for (i=0;i<N_z_layers;i++)
    {
      mean_velocity[i] = mean_velocity[i] + (dt/(T+dt))*(sum_values_per_layer[i]-mean_velocity[i]);
      mean_velocity_u2[i] = mean_velocity_u2[i] + (dt/(T+dt))*(sum_values_per_layer_u2[i] - mean_velocity_u2[i]);
      mean_velocity_u3[i] = mean_velocity_u3[i] + (dt/(T+dt))*(sum_values_per_layer_u3[i] - mean_velocity_u3[i]);

      rms_velocity[i] = rms_velocity[i]+(dt/(T+dt))*(rms_values_per_layer[i] - rms_velocity[i] );

      R_xx[i]  =  R_xx[i] + (dt/(T+dt)) * (sum_values_per_layer_xx[i] - R_xx[i]);
      R_yy[i]  =  R_yy[i] + (dt/(T+dt)) * (sum_values_per_layer_yy[i] - R_yy[i]);
      R_zz[i]  =  R_zz[i] + (dt/(T+dt)) * (sum_values_per_layer_zz[i] - R_zz[i]);
      R_xy[i]  =  R_xy[i] + (dt/(T+dt)) * (sum_values_per_layer_xy[i] - R_xy[i]);
      R_xz[i]  =  R_xz[i] + (dt/(T+dt)) * (sum_values_per_layer_xz[i] - R_xz[i]);
      R_yz[i]  =  R_yz[i] + (dt/(T+dt)) * (sum_values_per_layer_yz[i] - R_yz[i]);
    }
  }
  else
  {
    // initial time
    for (i=0;i<N_z_layers;i++)
    {
      mean_velocity[i] = sum_values_per_layer[i];
      mean_velocity_u2[i] = sum_values_per_layer_u2[i];
      mean_velocity_u3[i] = sum_values_per_layer_u3[i];

      rms_velocity[i] = rms_values_per_layer[i];

      R_xx[i] = sum_values_per_layer_xx[i];
      R_yy[i] = sum_values_per_layer_yy[i];
      R_zz[i] = sum_values_per_layer_zz[i];
      R_xy[i] = sum_values_per_layer_xy[i];
      R_xz[i] = sum_values_per_layer_xz[i];
      R_yz[i] = sum_values_per_layer_yz[i];
    }
  }

  // derivative of the current mean velocity profile
  a=0;
  b=0;
  D=0;
  for (i=0;i<N_z_layers;i++)
  {
    // lower boundary, one-sided difference
    if (i==0)
    {
      d_mean_velocity[i] = mean_velocity[i+1]/coord_z_layers[i+1];
    }
    else
    {
      // upper boundary, one-sided difference
      if (i==N_z_layers-1)
      {
        d_mean_velocity[i] = (mean_velocity[i-1]-mean_velocity[i])/(coord_z_layers[i-1]-coord_z_layers[i]);
      }
      else
      {
	// interior points, two-sided difference
        D = (coord_z_layers[i]-coord_z_layers[i-1]) * (coord_z_layers[i+1]*coord_z_layers[i+1]
						       -coord_z_layers[i-1]*coord_z_layers[i-1])
          -(coord_z_layers[i+1]-coord_z_layers[i-1])*(coord_z_layers[i]*coord_z_layers[i]
						      -coord_z_layers[i-1]*coord_z_layers[i-1]);

        b = ( (mean_velocity[i]-mean_velocity[i-1])* (coord_z_layers[i+1]*coord_z_layers[i+1]
						      -coord_z_layers[i-1]*coord_z_layers[i-1])
          - (mean_velocity[i+1]-mean_velocity[i-1])*(coord_z_layers[i]*coord_z_layers[i]
						     -coord_z_layers[i-1]*coord_z_layers[i-1]) )/D;

        a = ( (coord_z_layers[i]-coord_z_layers[i-1])*(mean_velocity[i+1]-mean_velocity[i-1])
          - (coord_z_layers[i+1]-coord_z_layers[i-1])*(mean_velocity[i]-mean_velocity[i-1]) )/D;
        d_mean_velocity[i] = 2*a*coord_z_layers[i]+b;
      }
    }
  }

  // compute time average derivative of mean velocity
  for (i=0;i<N_z_layers;i++)
  {
      dmean_velocity[i] = dmean_velocity[i] + (dt/(T+dt))*(d_mean_velocity[i]-dmean_velocity[i]);
  }


  // output of the results, velocity profiles
  for (i=0;i<N_z_layers;i++)
  {
    OutPut("mean velo : time " << T << 
	   " " << setw(10) << coord_z_layers[i] <<
	   " " << setw(10) <<  RE*(1-fabs(1-coord_z_layers[i])) <<
	   " mean  " << setw(10) <<  mean_velocity[i] <<
	   " rms " << setw(10) <<  rms_velocity[i] <<
	   " d_mean " <<setw(10)<< d_mean_velocity[i]<< 
	   " d_mean_av " <<setw(10)<< dmean_velocity[i]<< endl);
  }
  
  // output of the results, Reynolds stresses
  for (i=0;i<N_z_layers;i++)
  {
      OutPut("time " << T << 
	   " " << setw(10) << coord_z_layers[i] <<
	   " " << setw(10) <<  RE*(1-fabs(1-coord_z_layers[i])) <<
	   " R_uu "<<setw(8)<<R_xx[i]- mean_velocity[i] * mean_velocity[i] <<
	   " R_vv "<<setw(8)<<R_zz[i]- mean_velocity_u3[i]* mean_velocity_u3[i]  <<
	   " R_ww "<<setw(8)<<R_yy[i]- mean_velocity_u2[i]* mean_velocity_u2[i]  <<
	   " R_uv "<<setw(8)<<R_xz[i]- mean_velocity[i]* mean_velocity_u3[i] <<
	   " R_uw "<<setw(8)<<R_xy[i]- mean_velocity[i]* mean_velocity_u2[i] <<
	   " R_vw "<<setw(8)<<R_yz[i]- mean_velocity_u2[i]* mean_velocity_u3[i]  <<
      endl);
  }

  // determination of the bulk velocity, by trapezoidal rule
  bulk_velocity = 0.0;
  for (i=0;i<N_z_layers-1;i++)
  {
    bulk_velocity += 0.5*( mean_velocity[i+1]+mean_velocity[i])*(coord_z_layers[i+1]-coord_z_layers[i]);
  }
  bulk_velocity /= 2.0;

  OutPut("bulk velo : time " << T << " " << setw(10) <<bulk_velocity << " " 
	 << " tau " << (TDatabase::ParamDB->INTERNAL_BULK_MEAN - bulk_velocity) + 1
	 << endl);
  TDatabase::ParamDB->INTERNAL_BULK_SIMULATION = bulk_velocity;

  delete sum_values_per_layer;
  delete number_of_sums;
}
