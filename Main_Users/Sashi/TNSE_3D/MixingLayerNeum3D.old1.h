// Navier-Stokes problem
// channel flow in 3D
//

#include <PeriodicJoint.h>
#define __MIXINGLAYERSLIP3D__

#define U_INFTY 1   


// ========================================================================
// example file
// ========================================================================

void ExampleFile()
{
  OutPut("Example: MixingLayerNeum3D.h, U_INFTY "<< U_INFTY<< " NOISE is P7 !!!" 
         << " ini corrected" << endl);
  OutPut("initial noise waves "  << TDatabase::ParamDB->P5 << " " <<  
         TDatabase::ParamDB->P6 << endl);
}
// ========================================================================
// exact solution
// ========================================================================
void InitialU1(double x, double y, double z, double *values)
{
  double w, sigma;
  double NOISE = TDatabase::ParamDB->P7;
  sigma = 1/TDatabase::ParamDB->P8;
  double wave1 = TDatabase::ParamDB->P5;
  double wave2 = TDatabase::ParamDB->P6;
  w = 2*y/sigma;
 
  if (w>=0)
    values[0] = U_INFTY * (1-exp(-2*w))/(1+exp(-2*w));
  else
    values[0] = U_INFTY * (exp(2*w)-1)/(exp(2*w)+1);    
  values[0] -= NOISE* U_INFTY *exp(-w*w)*8*y*
     (cos(4*Pi*x)+cos(10*Pi*x)+cos(wave1*Pi*(z-1))+cos(wave2*Pi*(z-1)))/(sigma*sigma);
  values[0] -= NOISE* U_INFTY *exp(-w*w)*cos(wave1*Pi*(z-1))*wave1*Pi;
  values[0] -= NOISE* U_INFTY *exp(-w*w)*cos(wave2*Pi*(z-1))*wave2*Pi;
 
}

void InitialU2(double x, double y, double z, double *values)
{
  double w, sigma;
  double  NOISE = TDatabase::ParamDB->P7;
  sigma = 1/TDatabase::ParamDB->P8;
  w = 2*y/sigma;
  values[0] = NOISE*U_INFTY*exp(-w*w)*sin(4*Pi*x)*4*Pi;
  values[0] += NOISE*U_INFTY*exp(-w*w)*sin(10*Pi*x)*10*Pi;
}

void InitialU3(double x, double y, double z, double *values)
{
  double w, sigma;
  double NOISE = TDatabase::ParamDB->P7;
  sigma = 1/TDatabase::ParamDB->P8;
  w = 2*y/sigma;
  values[0] = NOISE*U_INFTY*exp(-w*w)*sin(4*Pi*x)*4*Pi;
  values[0] += NOISE*U_INFTY*exp(-w*w)*sin(10*Pi*x)*10*Pi;
}

void InitialP(double x, double y,  double z, double *values)
{
  values[0] = 0;
  TDatabase::ParamDB->INTERNAL_PERIODIC_IDENTITY = 1;
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
  if ((fabs(y-1)<1e-6)||(fabs(y+1)<1e-6))
  {
    cond = SLIP_FRICTION_PENETRATION_RESISTANCE;
    TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION = 1;
  }
  else
    cond = NEUMANN;
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
    coeff[3] = 0;
  }
}


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
      if (!(joint->GetType() == BoundaryFace))
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
        // check if vertex on z = 2
        if (fabs(Z[l1]-2)<1e-5)
          z_vert_2++;
        // check if vertex on x = -1
        if (fabs(X[l1]+1)<1e-5)
          x_vert_m1++;
        // check if vertex on x = 1        
        if (fabs(X[l1]-1)<1e-5)
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
          z1_vert_0 =  z1_vert_2 = x1_vert_m1 = x1_vert_1 = 0;
          // compute coordinates of vertices
          for (l1=0;l1<TmpLen[jj];l1++)
          {
            cell1->GetVertex(TmpFV[jj*MaxLen1+l1])->
              GetCoords(X1[l1], Y1[l1], Z1[l1]);
            // check if vertex on z = 0
            if (fabs(Z1[l1])<1e-5)
              z1_vert_0++;
            // check if vertex on z = 2
            if (fabs(Z1[l1]-2)<1e-5)
              z1_vert_2++;
            // check if vertex on x = -1
            if (fabs(X1[l1]+1)<1e-5)
              x1_vert_m1++;
            // check if vertex on x = 1        
            if (fabs(X1[l1]-1)<1e-5)
              x1_vert_1++;
            x1+=X1[l1];
            y1+=Y1[l1];
            z1+=Z1[l1];
          }
          x1 /= TmpLen[jj];
          y1 /= TmpLen[jj];
          z1 /= TmpLen[jj];
          found = 0;
          //  OutPut("baryopp " << x1 << " " << y1 << " " << z1 << endl);
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


// thickness of spanwise vorticity
void ComputeMomentumThickness3D(TFEFunction3D *Vorticity, double *thickness)
{
  const int MAX_Y_COORD = 1000;
  int i,j,k,l,ll,mm,found,max_lines,N_Cells,N_Faces,N_Vort;
  int middle[8],middle_plane;
  TCollection *Coll;
  TBaseCell *cell;
  double X[MaxN_QuadPoints_3D], Y[MaxN_QuadPoints_3D];
  double Z[MaxN_QuadPoints_3D];
  double l1,l2,x,z,val0[4],val1[4],val2[4],val3[4],val4[4];
  double integral = 0;
  TFESpace3D *vorticity_space;
  TJoint *joint, *joint1;
  const int *TmpFV, *TmpLen;
  int  MaxLen1, in_xz_plane;
  double global_y_coord[MAX_Y_COORD], aver_vort[MAX_Y_COORD], h;

  vorticity_space=Vorticity->GetFESpace3D();

  // get pointer to set of mesh cells which define the fe space
  Coll = vorticity_space->GetCollection();
  // get number of mesh cells
  N_Cells = Coll->GetN_Cells();

  max_lines = 0;
  // loop over all mesh cells
  for(i=0;i<N_Cells;i++)
  {
    // get current mesh cell
    cell = Coll->GetCell(i);  
    N_Faces = cell->GetN_Faces();
    // loop over all faces
    for (j=0;j<N_Faces;j++)
    {
       joint1=cell->GetJoint(j);
       // find vertices on the face
       cell->GetShapeDesc()->GetFaceVertex(TmpFV, TmpLen, MaxLen1);
       // compute coordinates of vertices
       for (ll=0;ll<TmpLen[j];ll++)
       {
         cell->GetVertex(TmpFV[j*MaxLen1+ll])->
            GetCoords(X[ll], Y[ll], Z[ll]);
       }
       // check if the face is orthogonal to (0,1,0)
       in_xz_plane = 0;
       for (ll = 0; ll < TmpLen[j]; ll++)
          for (mm = ll+1; mm < TmpLen[j]; mm++)
             if (fabs(Y[ll]-Y[mm])<1e-5)
              in_xz_plane++;
       if (in_xz_plane < 6)
          continue;
       // compute integral on face using two times the edge midpoint
       // rule
       // for (ll = 0; ll < TmpLen[j]; ll++)
       //  OutPut(X[ll] << " " << Y[ll] << " " << Z[ll] << endl);
        
       // compute area of the face (l1 * l2)
       l1 = sqrt((X[0]-X[1])*(X[0]-X[1])+(Z[0]-Z[1])*(Z[0]-Z[1]));
       l2 = sqrt((X[2]-X[1])*(X[2]-X[1])+(Z[2]-Z[1])*(Z[2]-Z[1]));
       // two times edge midpoint rule
       x = (X[0]+X[1])/2.0;
       z = (Z[0]+Z[1])/2.0;
       Vorticity->FindGradientLocal(cell, i,x,0,z,val0);
       x = (X[0]+X[3])/2.0;
       z = (Z[0]+Z[3])/2.0;
       Vorticity->FindGradientLocal(cell, i,x,0,z,val1);
       x = (X[2]+X[3])/2.0;
       z = (Z[2]+Z[3])/2.0;
       Vorticity->FindGradientLocal(cell, i,x,0,z,val2);
       x = (X[2]+X[1])/2.0;
       z = (Z[2]+Z[1])/2.0;
       Vorticity->FindGradientLocal(cell, i,x,0,z,val3);
       // opposite vertices
       if (fabs(X[2]-X[0])<1e-5 || fabs(Z[2]-Z[0])<1e-5)
       {
         OutPut("No opposite vertices in ComputeMomentumthickness !!!"<<endl);
         exit(4711);
       }
       x = (X[2]+X[0])/2.0;
       z = (Z[2]+Z[0])/2.0;
       Vorticity->FindGradientLocal(cell, i,x,0,z,val4);
    
       integral = l1*l2*(val0[0]*val0[0]+val1[0]*val1[0]+val2[0]*val2[0]+val3[0]*val3[0]
                         +2*val4[0]*val4[0])/6;
       
       found = 0;
       for (ll=0;ll<max_lines;ll++)
         if (fabs(global_y_coord[ll] - Y[0])<1e-5)
         {
           found++;
           aver_vort[ll] += integral;
         }
       if (!found)
       {
          global_y_coord[max_lines] = Y[0];
          aver_vort[max_lines] = integral;
          max_lines++;
       }
    }
  }

  // compute averages
  // scaling factor 0.25 from averaging over the plane
  // scaling facter 0.5 for inner planes (y \neq -1, 1) since they 
  // are computed twice,
  // inner plane are summed up with factor one in the next step, 
  // outer plane only with factor 1/2 -> apply scaling already here
  for (i=0;i<max_lines;i++)
  {
     //if (fabs(global_y_coord[i] -1 ) < 1e-5 || fabs(global_y_coord[i]  +1 ) < 1e-5)
     // aver_vort[i] /= 8.0;
     //else
      aver_vort[i] /= 8.0;
  }
  //for (i=0;i<max_lines;i++)
  //{
  //   OutPut(global_y_coord[i] << " " << aver_vort[i] << endl);
  //}

  // compute line integral (equidistant case)
  h = 2.0/(max_lines-1);
  integral = 0;
  for (i=0;i<max_lines;i++)
  {
    integral += aver_vort[i];
  }
  integral *= h;
  thickness[0] = 0.5 - integral/(4.0*U_INFTY*U_INFTY);

  return;


}
// thickness of spanwise vorticity
void ComputeVorticityThickness3D(TFEFunction3D *Vorticity, double *thickness)
{
  int i,j,k,l,found,max_lines,N_Cells,N_Faces,N_Vort;
  int middle[8],middle_plane;
  TCollection *Coll;
  TBaseCell *cell;
  double X[MaxN_QuadPoints_3D], Y[MaxN_QuadPoints_3D];
  double Z[MaxN_QuadPoints_3D];
  double x,z,l1,l2,val0[4],val1[4],val2[4],val3[4],val4[4];
  double integral = 0;
  TFESpace3D *vorticity_space;
    

  vorticity_space=Vorticity->GetFESpace3D();

  // get pointer to set of mesh cells which define the fe space
  Coll = vorticity_space->GetCollection();
  // get number of mesh cells
  N_Cells = Coll->GetN_Cells();

  max_lines = 0;
  // loop over all mesh cells
  for(i=0;i<N_Cells;i++)
  {
    // get current mesh cell
    cell = Coll->GetCell(i);  
    middle_plane = 0;
    for (j=0;j< 8;j++)
    {
      middle[j] = 0;
    }
    // for all vertices
    for (j=0;j< 8;j++)
    {
      cell->GetVertex(j)->GetCoords(X[j], Y[j], Z[j]);
      if (fabs(Y[j])<1e-5)
      {
        middle_plane++;
        middle[j] = 1;
      }
    }
    // no face on y = 0
    if (middle_plane != 4)
      continue;

    if ((middle[0])&&(middle[1])&&(middle[2]))
    {
      l1 = sqrt((X[0]-X[1])*(X[0]-X[1])+(Z[0]-Z[1])*(Z[0]-Z[1]));
      l2 = sqrt((X[2]-X[1])*(X[2]-X[1])+(Z[2]-Z[1])*(Z[2]-Z[1]));
      // two times edge midpoint rule
      x = (X[0]+X[1])/2.0;
      z = (Z[0]+Z[1])/2.0;
      Vorticity->FindGradientLocal(cell, i,x,0,z,val0);
      x = (X[0]+X[3])/2.0;
      z = (Z[0]+Z[3])/2.0;
      Vorticity->FindGradientLocal(cell, i,x,0,z,val1);
      x = (X[2]+X[3])/2.0;
      z = (Z[2]+Z[3])/2.0;
      Vorticity->FindGradientLocal(cell, i,x,0,z,val2);
      x = (X[2]+X[1])/2.0;
      z = (Z[2]+Z[1])/2.0;
      Vorticity->FindGradientLocal(cell, i,x,0,z,val3);
      x = (X[2]+X[0])/2.0;
      z = (Z[2]+Z[0])/2.0;
      Vorticity->FindGradientLocal(cell, i,x,0,z,val4);
    }
    if ((middle[0])&&(middle[1])&&(middle[4]))
    {
      l1 = sqrt((X[0]-X[1])*(X[0]-X[1])+(Z[0]-Z[1])*(Z[0]-Z[1]));
      l2 = sqrt((X[0]-X[4])*(X[0]-X[4])+(Z[0]-Z[4])*(Z[0]-Z[4]));
      // two times edge midpoint rule
      x = (X[0]+X[1])/2.0;
      z = (Z[0]+Z[1])/2.0;
      Vorticity->FindGradientLocal(cell, i,x,0,z,val0);
      x = (X[0]+X[4])/2.0;
      z = (Z[0]+Z[4])/2.0;
      Vorticity->FindGradientLocal(cell, i,x,0,z,val1);
      x = (X[1]+X[5])/2.0;
      z = (Z[1]+Z[5])/2.0;
      Vorticity->FindGradientLocal(cell, i,x,0,z,val2);
      x = (X[4]+X[5])/2.0;
      z = (Z[4]+Z[5])/2.0;
      Vorticity->FindGradientLocal(cell, i,x,0,z,val3);
      x = (X[5]+X[0])/2.0;
      z = (Z[5]+Z[0])/2.0;
      Vorticity->FindGradientLocal(cell, i,x,0,z,val4);
    }
    // face 0 4 3 7
    if ((middle[0])&&(middle[3])&&(middle[4]))
    {
      l1 = sqrt((X[0]-X[3])*(X[0]-X[3])+(Z[0]-Z[3])*(Z[0]-Z[3]));
      l2 = sqrt((X[0]-X[4])*(X[0]-X[4])+(Z[0]-Z[4])*(Z[0]-Z[4]));
      // two times edge midpoint rule
      x = (X[0]+X[3])/2.0;
      z = (Z[0]+Z[3])/2.0;
      Vorticity->FindGradientLocal(cell, i,x,0,z,val0);
      x = (X[0]+X[4])/2.0;
      z = (Z[0]+Z[4])/2.0;
      Vorticity->FindGradientLocal(cell, i,x,0,z,val1);
      x = (X[3]+X[7])/2.0;
      z = (Z[3]+Z[7])/2.0;
      Vorticity->FindGradientLocal(cell, i,x,0,z,val2);
      x = (X[4]+X[7])/2.0;
      z = (Z[4]+Z[7])/2.0;
      Vorticity->FindGradientLocal(cell, i,x,0,z,val3);
      x = (X[7]+X[0])/2.0;
      z = (Z[7]+Z[0])/2.0;
      Vorticity->FindGradientLocal(cell, i,x,0,z,val4);
    }
    // face 2 3 6 7
    if ((middle[2])&&(middle[3])&&(middle[6]))
    {
      l1 = sqrt((X[3]-X[2])*(X[3]-X[2])+(Z[3]-Z[2])*(Z[3]-Z[2]));
      l2 = sqrt((X[2]-X[6])*(X[2]-X[6])+(Z[2]-Z[6])*(Z[2]-Z[6]));
      // two times edge midpoint rule
      x = (X[3]+X[2])/2.0;
      z = (Z[3]+Z[2])/2.0;
      Vorticity->FindGradientLocal(cell, i,x,0,z,val0);
      x = (X[2]+X[6])/2.0;
      z = (Z[2]+Z[6])/2.0;
      Vorticity->FindGradientLocal(cell, i,x,0,z,val1);
      x = (X[3]+X[7])/2.0;
      z = (Z[3]+Z[7])/2.0;
      Vorticity->FindGradientLocal(cell, i,x,0,z,val2);
      x = (X[7]+X[6])/2.0;
      z = (Z[7]+Z[6])/2.0;
      Vorticity->FindGradientLocal(cell, i,x,0,z,val3);
      x = (X[2]+X[7])/2.0;
      z = (Z[2]+Z[7])/2.0;
      Vorticity->FindGradientLocal(cell, i,x,0,z,val4);
    }
    // face 1 2 5 6
    if ((middle[1])&&(middle[2])&&(middle[5]))
    {
      l1 = sqrt((X[1]-X[2])*(X[1]-X[2])+(Z[1]-Z[2])*(Z[1]-Z[2]));
      l2 = sqrt((X[1]-X[5])*(X[1]-X[5])+(Z[1]-Z[5])*(Z[1]-Z[5]));
      // two times edge midpoint rule
      x = (X[1]+X[2])/2.0;
      z = (Z[1]+Z[2])/2.0;
      Vorticity->FindGradientLocal(cell, i,x,0,z,val0);
      x = (X[1]+X[5])/2.0;
      z = (Z[1]+Z[5])/2.0;
      Vorticity->FindGradientLocal(cell, i,x,0,z,val1);
      x = (X[2]+X[6])/2.0;
      z = (Z[2]+Z[6])/2.0;
      Vorticity->FindGradientLocal(cell, i,x,0,z,val2);
      x = (X[5]+X[6])/2.0;
      z = (Z[5]+Z[6])/2.0;
      Vorticity->FindGradientLocal(cell, i,x,0,z,val3);
      x = (X[1]+X[6])/2.0;
      z = (Z[1]+Z[6])/2.0;
      Vorticity->FindGradientLocal(cell, i,x,0,z,val4);
    }
    // face 4 5 6 7
    if ((middle[4])&&(middle[4])&&(middle[6]))
    {
      l1 = sqrt((X[4]-X[5])*(X[4]-X[5])+(Z[4]-Z[5])*(Z[4]-Z[5]));
      l2 = sqrt((X[4]-X[7])*(X[4]-X[7])+(Z[4]-Z[7])*(Z[4]-Z[7]));
      // two times edge midpoint rule
      x = (X[4]+X[5])/2.0;
      z = (Z[4]+Z[5])/2.0;
      Vorticity->FindGradientLocal(cell, i,x,0,z,val0);
      x = (X[4]+X[7])/2.0;
      z = (Z[4]+Z[7])/2.0;
      Vorticity->FindGradientLocal(cell, i,x,0,z,val1);
      x = (X[5]+X[6])/2.0;
      z = (Z[5]+Z[6])/2.0;
      Vorticity->FindGradientLocal(cell, i,x,0,z,val2);
      x = (X[7]+X[6])/2.0;
      z = (Z[7]+Z[6])/2.0;
      Vorticity->FindGradientLocal(cell, i,x,0,z,val3);
      x = (X[4]+X[6])/2.0;
      z = (Z[4]+Z[6])/2.0;
      Vorticity->FindGradientLocal(cell, i,x,0,z,val4);
    }
    //OutPut(l1 << " " << l2 << " " << val0[0] << " " <<  val1[0] << " " << val2[0] << " " 
    //      <<  val3[0] << " " << val4[0] << endl);
    // update integral 
    integral += l1*l2*(val0[0]+val1[0]+val2[0]+val3[0]+2*val4[0])/6;
  }
 
  integral /= 2; // computed twice 
  integral /= 4; //mean value

  thickness[0] =  2*U_INFTY/fabs(integral); 
  return;

}
