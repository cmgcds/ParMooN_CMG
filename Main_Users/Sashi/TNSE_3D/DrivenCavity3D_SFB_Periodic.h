// Navier-Stokes problem
// channel flow in 3D
//

#include <PeriodicJoint.h>
#define __DC3D_SFB_PERIODIC__


// ========================================================================
// example file
// ========================================================================

void ExampleFile()
{
  int range;

  OutPut("Example: DrivenCavity3D_SFB_Periodic.h ") ;
  OutPut("inflow " << TDatabase::ParamDB->P6);
  OutPut(" upper lid " << TDatabase::ParamDB->P5);
  OutPut(" left " << (int)TDatabase::ParamDB->P7);
  OutPut(" right " << (int)TDatabase::ParamDB->P8);
  OutPut(" lower " << (int)TDatabase::ParamDB->P9 << endl);

  range = (int)TDatabase::ParamDB->P7;
  if ((range<1)||(range>30))
  {
      OutPut("left boundary out of range !!!"<< endl);
      exit(4711);
  }

  range = (int)TDatabase::ParamDB->P8;
  if ((range<1)||(range>30))
  {
      OutPut("right boundary out of range !!!"<< endl);
      exit(4711);
  }

  range = (int)TDatabase::ParamDB->P9;
  if ((range<1)||(range>30))
  {
      OutPut("lower boundary out of range !!!"<< endl);
      exit(4711);
  }
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
   int lower = (int)TDatabase::ParamDB->P9;

  if ((fabs(y)<1e-6)&&((x>lower/32.0)&&(x<(lower+2)/32.0)))
      cond = NEUMANN;
  else
      cond = DIRICHLET;
  TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 0; 
}

// value of boundary condition
void U1BoundValue(double x, double y, double z, double &value)
{
    int range;
/*    double y, fac = 1024; // 32^2
    

  switch(BdComp)
  {
     case 0: 
        value = 0;
        break;
     case 1:
	range = (int)TDatabase::ParamDB->P8;
	if ((Param>range/32.0)&&(Param<(range+1)/32.0))
        {
           y = Param;
           value =  6*(y-range/32.0)*(y-(range+1)/32.0)*fac;
           //  OutPut(value << " ");
        }
	else
	    value = 0;
	break;
    case 2: if(Param<0.00001 || Param>0.99999) 
              value = 0;
            else
               value = TDatabase::ParamDB->P5/TDatabase::ParamDB->P6;
               //value = 0;
            break;
    case 3: 
	range = (int)TDatabase::ParamDB->P7;
        y = 1-Param;
	if ((y>range/32.0)&&(y<(range+1)/32.0))
        {
           value = -6*(y-range/32.0)*(y-(range+1)/32.0)*fac;
        }
	else
	    value = 0;
	break;
    default: cout << "wrong boundary part number" << endl;
  }
*/
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
  int i;
  double *coeff;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    coeff[0] = 1e-6/TDatabase::ParamDB->P6;
    coeff[1] = 0;
    coeff[2] = 0;
    coeff[3] = 0;
  }
}


// ========================================================================
// periodic b.c. for z=0 and z=0.5
// ========================================================================
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
      z_vert_0 =  z_vert_2 = 0;
      // compute coordinates of vertices
      for (l1=0;l1<TmpLen[j];l1++)
      {
        cell->GetVertex(TmpFV[j*MaxLen+l1])->GetCoords(X[l1], Y[l1], Z[l1]);
        // check if vertex on z = 0
        if (fabs(Z[l1])<1e-5)
          z_vert_0++;
        // check if vertex on z = 0.5
        if (fabs(Z[l1]-0.5)<1e-5)
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
          for (l1=0;l1<TmpLen[jj];l1++)
          {
            cell1->GetVertex(TmpFV[jj*MaxLen1+l1])->
              GetCoords(X1[l1], Y1[l1], Z1[l1]);
            // check if vertex on z = 0
            if (fabs(Z1[l1])<1e-5)
              z1_vert_0++;
            // check if vertex on z = 0.5
            if (fabs(Z1[l1]-0.5)<1e-5)
              z1_vert_2++;
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
