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

void  FindPeriodicXDOFs(TFESpace3D *USpace, int dof_length,
                       int *left_dof, int *right_dof, 
                        int *direction, int *new_length)
{
  int i, j, N_Cells, left[8], right[8], rightno, leftno;
  int left_curr=0, right_curr=0, found;
  int *UGlobalNumbers, *UBeginIndex, *DOF;
  TCollection *Coll;
  TBaseCell *cell;
  TJoint *joint;
  TBoundFace *boundface;
  TBoundComp3D *BoundComp;
  double *left_y,*left_z, *right_y,*right_z, tmp;
  double X[MaxN_QuadPoints_3D], Y[MaxN_QuadPoints_3D];
  double Z[MaxN_QuadPoints_3D];

  memset(left_dof, 0, dof_length*SizeOfInt);
  memset(right_dof, 0, dof_length*SizeOfInt);
  left_y = new double[dof_length];
  left_z = new double[dof_length];
  right_y = new double[dof_length];
  right_z = new double[dof_length];

  UGlobalNumbers = USpace->GetGlobalNumbers();
  UBeginIndex = USpace->GetBeginIndex();
  
  Coll = USpace->GetCollection();
  N_Cells = Coll->GetN_Cells();

  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    for (j=0;j<8;j++)
      left[j] = right[j] = 0;

    // get coordinates for all vertices
    for (j=0;j<8;j++)
    {
      cell->GetVertex(j)->GetCoords(X[j], Y[j], Z[j]);
      // check if a face is on x=-1 or x=1
      if (fabs(X[j]+1)<1e-5)
        left[j] = 1;
      if (fabs(X[j]-1)<1e-5)
        right[j] = 1;
    }
    
    leftno = rightno = 0;
    for (j=0;j<8;j++)
    {
      leftno+=left[j];
      rightno += right[j];
    }
    if ( (leftno<4)&&(rightno<4) )
      continue;
    
    DOF = UGlobalNumbers + UBeginIndex[i];
    // face on x=-1
    if (leftno==4)
    {
      // face 0 1 2 3
      if ((left[0]==1)&&(left[1]==1)&&(left[2]==1))
      {
        OutPut("periodic 1 not implemented !!!"<< endl);
        exit(4711);
       }
      // face 0 1 4 5
      if ((left[0]==1)&&(left[1]==1)&&(left[4]==1))
      {
        left_dof[left_curr] = DOF[0];
        left_y[left_curr] = Y[0];
        left_z[left_curr] = Z[0];
        left_curr++;
        left_dof[left_curr] = DOF[1];
        left_y[left_curr] = (Y[0]+Y[1])/2;
        left_z[left_curr] = (Z[0]+Z[1])/2;
        left_curr++;
        left_dof[left_curr] = DOF[2];
        left_y[left_curr] = Y[1];
        left_z[left_curr] = Z[1];
        left_curr++;
        left_dof[left_curr] = DOF[9];
        left_y[left_curr] = (Y[0]+Y[4])/2;
        left_z[left_curr] = (Z[0]+Z[4])/2;
        left_curr++;
        left_dof[left_curr] = DOF[10];
        left_y[left_curr] = (Y[0]+Y[5])/2;
        left_z[left_curr] = (Z[0]+Z[5])/2;
        left_curr++;
        left_dof[left_curr] = DOF[11];
        left_y[left_curr] = (Y[1]+Y[5])/2;
        left_z[left_curr] = (Z[1]+Z[5])/2;
        left_curr++;
        left_dof[left_curr] = DOF[18];
        left_y[left_curr] = Y[4];
        left_z[left_curr] = Z[4];
        left_curr++;
        left_dof[left_curr] = DOF[19];
        left_y[left_curr] = (Y[4]+Y[5])/2;
        left_z[left_curr] = (Z[4]+Z[5])/2;
        left_curr++;
        left_dof[left_curr] = DOF[20];
        left_y[left_curr] = Y[5];
        left_z[left_curr] = Z[5];
        left_curr++;       
      }
      // face 0 3 4 7
      if ((left[0]==1)&&(left[3]==1)&&(left[4]==1))
      {
        left_dof[left_curr] = DOF[0];
        left_y[left_curr] = Y[0];
        left_z[left_curr] = Z[0];
        left_curr++;
        left_dof[left_curr] = DOF[3];
        left_y[left_curr] = (Y[0]+Y[3])/2;
        left_z[left_curr] = (Z[0]+Z[3])/2;
        left_curr++;
        left_dof[left_curr] = DOF[6];
        left_y[left_curr] = Y[3];
        left_z[left_curr] = Z[3];
        left_curr++;
        left_dof[left_curr] = DOF[9];
        left_y[left_curr] = (Y[0]+Y[4])/2;
        left_z[left_curr] = (Z[0]+Z[4])/2;
        left_curr++;
        left_dof[left_curr] = DOF[12];
        left_y[left_curr] = (Y[0]+Y[7])/2;
        left_z[left_curr] = (Z[0]+Z[7])/2;
        left_curr++;
        left_dof[left_curr] = DOF[15];
        left_y[left_curr] = (Y[3]+Y[7])/2;
        left_z[left_curr] = (Z[3]+Z[7])/2;
        left_curr++;
        left_dof[left_curr] = DOF[18];
        left_y[left_curr] = Y[4];
        left_z[left_curr] = Z[4];
        left_curr++;
        left_dof[left_curr] = DOF[21];
        left_y[left_curr] = (Y[4]+Y[7])/2;
        left_z[left_curr] = (Z[4]+Z[7])/2;
        left_curr++;
        left_dof[left_curr] = DOF[24];
        left_y[left_curr] = Y[7];
        left_z[left_curr] = Z[7];
        left_curr++;       
      }
      // face 4 5 6 7
      if ((left[4]==1)&&(left[5]==1)&&(left[6]==1))
      {
        OutPut("periodic 2 not implemented !!!"<< endl);
        exit(4711);
      }
      // face 1 2 5 6
      if ((left[1]==1)&&(left[2]==1)&&(left[5]==1))
      {
        OutPut("periodic 3 not implemented !!!"<< endl);
        exit(4711);
       }
      // face 2 3 6 7
      if ((left[2]==1)&&(left[3]==1)&&(left[6]==1))
      {
        OutPut("periodic 4 not implemented !!!"<< endl);
        exit(4711);
       }
    }
    if (rightno==4)
    {
      // face 0 1 2 3
      if ((right[0]==1)&&(right[1]==1)&&(right[2]==1))
      {
        OutPut("periodic 5 not implemented !!!"<< endl);
        exit(4711);
       }
      // face 0 1 4 5
      if ((right[0]==1)&&(right[1]==1)&&(right[4]==1))
      {
        right_dof[right_curr] = DOF[0];
        right_y[right_curr] = Y[0];
        right_z[right_curr] = Z[0];
        right_curr++;
        right_dof[right_curr] = DOF[1];
        right_y[right_curr] = (Y[0]+Y[1])/2;
        right_z[right_curr] = (Z[0]+Z[1])/2;
        right_curr++;
        right_dof[right_curr] = DOF[2];
        right_y[right_curr] = Y[1];
        right_z[right_curr] = Z[1];
        right_curr++;
        right_dof[right_curr] = DOF[9];
        right_y[right_curr] = (Y[0]+Y[4])/2;
        right_z[right_curr] = (Z[0]+Z[4])/2;
        right_curr++;
        right_dof[right_curr] = DOF[10];
        right_y[right_curr] = (Y[0]+Y[5])/2;
        right_z[right_curr] = (Z[0]+Z[5])/2;
        right_curr++;
        right_dof[right_curr] = DOF[11];
        right_y[right_curr] = (Y[1]+Y[5])/2;
        right_z[right_curr] = (Z[1]+Z[5])/2;
        right_curr++;
        right_dof[right_curr] = DOF[18];
        right_y[right_curr] = Y[4];
        right_z[right_curr] = Z[4];
        right_curr++;
        right_dof[right_curr] = DOF[19];
        right_y[right_curr] = (Y[4]+Y[5])/2;
        right_z[right_curr] = (Z[4]+Z[5])/2;
        right_curr++;
        right_dof[right_curr] = DOF[20];
        right_y[right_curr] = Y[5];
        right_z[right_curr] = Z[5];
        right_curr++;       
      }
      // face 0 3 4 7
      if ((right[0]==1)&&(right[3]==1)&&(right[4]==1))
      {
        right_dof[right_curr] = DOF[0];
        right_y[right_curr] = Y[0];
        right_z[right_curr] = Z[0];
        right_curr++;
        right_dof[right_curr] = DOF[3];
        right_y[right_curr] = (Y[0]+Y[3])/2;
        right_z[right_curr] = (Z[0]+Z[3])/2;
        right_curr++;
        right_dof[right_curr] = DOF[6];
        right_y[right_curr] = Y[3];
        right_z[right_curr] = Z[3];
        right_curr++;
        right_dof[right_curr] = DOF[9];
        right_y[right_curr] = (Y[0]+Y[4])/2;
        right_z[right_curr] = (Z[0]+Z[4])/2;
        right_curr++;
        right_dof[right_curr] = DOF[12];
        right_y[right_curr] = (Y[0]+Y[7])/2;
        right_z[right_curr] = (Z[0]+Z[7])/2;
        right_curr++;
        right_dof[right_curr] = DOF[15];
        right_y[right_curr] = (Y[3]+Y[7])/2;
        right_z[right_curr] = (Z[3]+Z[7])/2;
        right_curr++;
        right_dof[right_curr] = DOF[18];
        right_y[right_curr] = Y[4];
        right_z[right_curr] = Z[4];
        right_curr++;
        right_dof[right_curr] = DOF[21];
        right_y[right_curr] = (Y[4]+Y[7])/2;
        right_z[right_curr] = (Z[4]+Z[7])/2;
        right_curr++;
        right_dof[right_curr] = DOF[24];
        right_y[right_curr] = Y[7];
        right_z[right_curr] = Z[7];
        right_curr++;       
      }
      // face 4 5 6 7
      if ((right[4]==1)&&(right[5]==1)&&(right[6]==1))
      {
        OutPut("periodic 6 not implemented !!!"<< endl);
        exit(4711);
      }
      // face 1 2 5 6
      if ((right[1]==1)&&(right[2]==1)&&(right[5]==1))
      {
        OutPut("periodic 7 not implemented !!!"<< endl);
        exit(4711);
       }
      // face 2 3 6 7
      if ((right[2]==1)&&(right[3]==1)&&(right[6]==1))
      {
        OutPut("periodic 8 not implemented !!!"<< endl);
        exit(4711);
       }
    }
      
  } // endfor i (cells)

  OutPut("left " << left_curr << " right " << right_curr << endl);

  leftno =0;
  for (i=0;i<left_curr;i++)
  {
    found = 0;
    // check if already in list
    for (j=0;j<leftno-1;j++)
    {     
      if (left_dof[i]==left_dof[j])
      {
        found = 1;
        break;
      }
    }
    if (found)
      continue;
    found = 0;
    // find right counterpart
    for (j=leftno;j< right_curr; j++)
    {
      if ((left_y[i]==right_y[j])&&(left_z[i]==right_z[j]))
      {
         left_dof[leftno] = left_dof[i];
         left_y[leftno] = left_y[i];
         left_z[leftno] = left_z[i];
         rightno =  right_dof[leftno];
         right_dof[leftno] = right_dof[j];
         right_dof[j] = rightno;
         tmp = right_y[leftno];
         right_y[leftno] = right_y[j];
         right_y[j] = tmp;
         tmp = right_z[leftno];
         right_z[leftno] = right_z[j];
         right_z[j] = tmp;
         found=1;
         leftno++;
      }
      if (found)
        break;
    }
    if (!found)
    {
      OutPut("no counterpart found !!!" << endl);
      exit(4711);
    }
  }

  for (i=0;i<leftno;i++)
  {
    OutPut(i << " " << left_dof[i] << " " <<  left_y[i] << " " << left_z[i] << endl);
    OutPut("     " << " " << right_dof[i] << " " <<  right_y[i] << " " << right_z[i] << endl);
    if (left_y[i]<0)
      direction[i] = -1;
    if (left_y[i]>0)
      direction[i] = 1;
    
  }

  new_length[0] = leftno;

  delete left_y;
  delete left_z;
  delete right_y;
  delete right_z;
} 
void  FindPeriodicZDOFs(TFESpace3D *USpace, int dof_length,
                        int *left_dof, int *right_dof, 
                        int *direction, int *new_length)
{
  int i, j, N_Cells, left[8], right[8], rightno, leftno;
  int left_curr=0, right_curr=0, found;
  int *UGlobalNumbers, *UBeginIndex, *DOF;
  TCollection *Coll;
  TBaseCell *cell;
  TJoint *joint;
  TBoundFace *boundface;
  TBoundComp3D *BoundComp;
  double *left_y,*left_x, *right_y,*right_x, tmp;
  double X[MaxN_QuadPoints_3D], Y[MaxN_QuadPoints_3D];
  double Z[MaxN_QuadPoints_3D];

  memset(left_dof, 0, dof_length*SizeOfInt);
  memset(right_dof, 0, dof_length*SizeOfInt);
  left_y = new double[dof_length];
  left_x = new double[dof_length];
  right_y = new double[dof_length];
  right_x = new double[dof_length];

  UGlobalNumbers = USpace->GetGlobalNumbers();
  UBeginIndex = USpace->GetBeginIndex();
  
  Coll = USpace->GetCollection();
  N_Cells = Coll->GetN_Cells();

  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    for (j=0;j<8;j++)
      left[j] = right[j] = 0;

    // get coordinates for all vertices
    for (j=0;j<8;j++)
    {
      cell->GetVertex(j)->GetCoords(X[j], Y[j], Z[j]);
      // check if a face is on z=0 or z=2
      if (fabs(Z[j])<1e-5)
      {
        left[j] = 1;
      }
      if (fabs(Z[j]-2)<1e-5)
      {
        right[j] = 1;
      }
    }
    
    leftno = rightno = 0;
    for (j=0;j<8;j++)
    {
      leftno+=left[j];
      rightno += right[j];
    }
    if ( (leftno<4)&&(rightno<4) )
      continue;
    
    DOF = UGlobalNumbers + UBeginIndex[i];
    // face on z=0
    if (leftno==4)
    {
      // face 0 1 2 3
      if ((left[0]==1)&&(left[1]==1)&&(left[2]==1))
      {
        left_dof[left_curr] = DOF[0];
        left_y[left_curr] = Y[0];
        left_x[left_curr] = X[0];
        left_curr++;
        left_dof[left_curr] = DOF[1];
        left_y[left_curr] = (Y[0]+Y[1])/2;
        left_x[left_curr] = (X[0]+X[1])/2;
        left_curr++;
        left_dof[left_curr] = DOF[2];
        left_y[left_curr] = Y[1];
        left_x[left_curr] = X[1];
        left_curr++;
        left_dof[left_curr] = DOF[3];
        left_y[left_curr] = (Y[0]+Y[3])/2;
        left_x[left_curr] = (X[0]+X[3])/2;
        left_curr++;
        left_dof[left_curr] = DOF[4];
        left_y[left_curr] = (Y[0]+Y[2])/2;
        left_x[left_curr] = (X[0]+X[2])/2;
        left_curr++;
        left_dof[left_curr] = DOF[5];
        left_y[left_curr] = (Y[1]+Y[2])/2;
        left_x[left_curr] = (X[1]+X[2])/2;
        left_curr++;
        left_dof[left_curr] = DOF[6];
        left_y[left_curr] = Y[3];
        left_x[left_curr] = X[3];
        left_curr++;
        left_dof[left_curr] = DOF[7];
        left_y[left_curr] = (Y[3]+Y[2])/2;
        left_x[left_curr] = (X[3]+X[2])/2;
        left_curr++;
        left_dof[left_curr] = DOF[8];
        left_y[left_curr] = Y[2];
        left_x[left_curr] = X[2];
        left_curr++;       
       }
      // face 0 1 4 5
      if ((left[0]==1)&&(left[1]==1)&&(left[4]==1))
      {
        left_dof[left_curr] = DOF[0];
        left_y[left_curr] = Y[0];
        left_x[left_curr] = X[0];
        left_curr++;
        left_dof[left_curr] = DOF[1];
        left_y[left_curr] = (Y[0]+Y[1])/2;
        left_x[left_curr] = (X[0]+X[1])/2;
        left_curr++;
        left_dof[left_curr] = DOF[2];
        left_y[left_curr] = Y[1];
        left_x[left_curr] = X[1];
        left_curr++;
        left_dof[left_curr] = DOF[9];
        left_y[left_curr] = (Y[0]+Y[4])/2;
        left_x[left_curr] = (X[0]+X[4])/2;
        left_curr++;
        left_dof[left_curr] = DOF[10];
        left_y[left_curr] = (Y[0]+Y[5])/2;
        left_x[left_curr] = (X[0]+X[5])/2;
        left_curr++;
        left_dof[left_curr] = DOF[11];
        left_y[left_curr] = (Y[1]+Y[5])/2;
        left_x[left_curr] = (X[1]+X[5])/2;
        left_curr++;
        left_dof[left_curr] = DOF[18];
        left_y[left_curr] = Y[4];
        left_x[left_curr] = X[4];
        left_curr++;
        left_dof[left_curr] = DOF[19];
        left_y[left_curr] = (Y[4]+Y[5])/2;
        left_x[left_curr] = (X[4]+X[5])/2;
        left_curr++;
        left_dof[left_curr] = DOF[20];
        left_y[left_curr] = Y[5];
        left_x[left_curr] = X[5];
        left_curr++;       
      }
      // face 0 3 4 7
      if ((left[0]==1)&&(left[3]==1)&&(left[4]==1))
      {
        left_dof[left_curr] = DOF[0];
        left_y[left_curr] = Y[0];
        left_x[left_curr] = X[0];
        left_curr++;
        left_dof[left_curr] = DOF[3];
        left_y[left_curr] = (Y[0]+Y[3])/2;
        left_x[left_curr] = (X[0]+X[3])/2;
        left_curr++;
        left_dof[left_curr] = DOF[6];
        left_y[left_curr] = Y[3];
        left_x[left_curr] = X[3];
        left_curr++;
        left_dof[left_curr] = DOF[9];
        left_y[left_curr] = (Y[0]+Y[4])/2;
        left_x[left_curr] = (X[0]+X[4])/2;
        left_curr++;
        left_dof[left_curr] = DOF[12];
        left_y[left_curr] = (Y[0]+Y[7])/2;
        left_x[left_curr] = (X[0]+X[7])/2;
        left_curr++;
        left_dof[left_curr] = DOF[15];
        left_y[left_curr] = (Y[3]+Y[7])/2;
        left_x[left_curr] = (X[3]+X[7])/2;
        left_curr++;
        left_dof[left_curr] = DOF[18];
        left_y[left_curr] = Y[4];
        left_x[left_curr] = X[4];
        left_curr++;
        left_dof[left_curr] = DOF[21];
        left_y[left_curr] = (Y[4]+Y[7])/2;
        left_x[left_curr] = (X[4]+X[7])/2;
        left_curr++;
        left_dof[left_curr] = DOF[24];
        left_y[left_curr] = Y[7];
        left_x[left_curr] = X[7];
        left_curr++;       
      }
      // face 4 5 6 7
      if ((left[4]==1)&&(left[5]==1)&&(left[6]==1))
      {
        OutPut("periodic 2 not implemented !!!"<< endl);
        exit(4711);
      }
      // face 1 2 5 6
      if ((left[1]==1)&&(left[2]==1)&&(left[5]==1))
      {
        OutPut("periodic 3 not implemented !!!"<< endl);
        exit(4711);
       }
      // face 2 3 6 7
      if ((left[2]==1)&&(left[3]==1)&&(left[6]==1))
      {
        OutPut("periodic 4 not implemented !!!"<< endl);
        exit(4711);
       }
    }
    if (rightno==4)
    {
      // face 0 1 2 3
      {
        right_dof[right_curr] = DOF[0];
        right_y[right_curr] = Y[0];
        right_x[right_curr] = X[0];
        right_curr++;
        right_dof[right_curr] = DOF[1];
        right_y[right_curr] = (Y[0]+Y[1])/2;
        right_x[right_curr] = (X[0]+X[1])/2;
        right_curr++;
        right_dof[right_curr] = DOF[2];
        right_y[right_curr] = Y[1];
        right_x[right_curr] = X[1];
        right_curr++;
        right_dof[right_curr] = DOF[3];
        right_y[right_curr] = (Y[0]+Y[3])/2;
        right_x[right_curr] = (X[0]+X[3])/2;
        right_curr++;
        right_dof[right_curr] = DOF[4];
        right_y[right_curr] = (Y[0]+Y[2])/2;
        right_x[right_curr] = (X[0]+X[2])/2;
        right_curr++;
        right_dof[right_curr] = DOF[5];
        right_y[right_curr] = (Y[1]+Y[2])/2;
        right_x[right_curr] = (X[1]+X[2])/2;
        right_curr++;
        right_dof[right_curr] = DOF[6];
        right_y[right_curr] = Y[3];
        right_x[right_curr] = X[3];
        right_curr++;
        right_dof[right_curr] = DOF[7];
        right_y[right_curr] = (Y[2]+Y[3])/2;
        right_x[right_curr] = (X[2]+X[3])/2;
        right_curr++;
        right_dof[right_curr] = DOF[8];
        right_y[right_curr] = Y[2];
        right_x[right_curr] = X[2];
        right_curr++;       
      }
      // face 0 1 4 5
      if ((right[0]==1)&&(right[1]==1)&&(right[4]==1))
      {
        right_dof[right_curr] = DOF[0];
        right_y[right_curr] = Y[0];
        right_x[right_curr] = X[0];
        right_curr++;
        right_dof[right_curr] = DOF[1];
        right_y[right_curr] = (Y[0]+Y[1])/2;
        right_x[right_curr] = (X[0]+X[1])/2;
        right_curr++;
        right_dof[right_curr] = DOF[2];
        right_y[right_curr] = Y[1];
        right_x[right_curr] = X[1];
        right_curr++;
        right_dof[right_curr] = DOF[9];
        right_y[right_curr] = (Y[0]+Y[4])/2;
        right_x[right_curr] = (X[0]+X[4])/2;
        right_curr++;
        right_dof[right_curr] = DOF[10];
        right_y[right_curr] = (Y[0]+Y[5])/2;
        right_x[right_curr] = (X[0]+X[5])/2;
        right_curr++;
        right_dof[right_curr] = DOF[11];
        right_y[right_curr] = (Y[1]+Y[5])/2;
        right_x[right_curr] = (X[1]+X[5])/2;
        right_curr++;
        right_dof[right_curr] = DOF[18];
        right_y[right_curr] = Y[4];
        right_x[right_curr] = X[4];
        right_curr++;
        right_dof[right_curr] = DOF[19];
        right_y[right_curr] = (Y[4]+Y[5])/2;
        right_x[right_curr] = (X[4]+X[5])/2;
        right_curr++;
        right_dof[right_curr] = DOF[20];
        right_y[right_curr] = Y[5];
        right_x[right_curr] = X[5];
        right_curr++;       
      }
      // face 0 3 4 7
      if ((right[0]==1)&&(right[3]==1)&&(right[4]==1))
      {
        right_dof[right_curr] = DOF[0];
        right_y[right_curr] = Y[0];
        right_x[right_curr] = X[0];
        right_curr++;
        right_dof[right_curr] = DOF[3];
        right_y[right_curr] = (Y[0]+Y[3])/2;
        right_x[right_curr] = (X[0]+X[3])/2;
        right_curr++;
        right_dof[right_curr] = DOF[6];
        right_y[right_curr] = Y[3];
        right_x[right_curr] = X[3];
        right_curr++;
        right_dof[right_curr] = DOF[9];
        right_y[right_curr] = (Y[0]+Y[4])/2;
        right_x[right_curr] = (X[0]+X[4])/2;
        right_curr++;
        right_dof[right_curr] = DOF[12];
        right_y[right_curr] = (Y[0]+Y[7])/2;
        right_x[right_curr] = (X[0]+X[7])/2;
        right_curr++;
        right_dof[right_curr] = DOF[15];
        right_y[right_curr] = (Y[3]+Y[7])/2;
        right_x[right_curr] = (X[3]+X[7])/2;
        right_curr++;
        right_dof[right_curr] = DOF[18];
        right_y[right_curr] = Y[4];
        right_x[right_curr] = X[4];
        right_curr++;
        right_dof[right_curr] = DOF[21];
        right_y[right_curr] = (Y[4]+Y[7])/2;
        right_x[right_curr] = (X[4]+X[7])/2;
        right_curr++;
        right_dof[right_curr] = DOF[24];
        right_y[right_curr] = Y[7];
        right_x[right_curr] = X[7];
        right_curr++;       
      }
      // face 4 5 6 7
      if ((right[4]==1)&&(right[5]==1)&&(right[6]==1))
      {
        OutPut("periodic 6 not implemented !!!"<< endl);
        exit(4711);
      }
      // face 1 2 5 6
      if ((right[1]==1)&&(right[2]==1)&&(right[5]==1))
      {
        OutPut("periodic 7 not implemented !!!"<< endl);
        exit(4711);
       }
      // face 2 3 6 7
      if ((right[2]==1)&&(right[3]==1)&&(right[6]==1))
      {
        OutPut("periodic 8 not implemented !!!"<< endl);
        exit(4711);
       }
    }
      
  } // endfor i (cells)

  OutPut("left " << left_curr << " right " << right_curr << endl);

  leftno =0;
  for (i=0;i<left_curr;i++)
  {
    found = 0;
    // check if already in list
    for (j=0;j<leftno-1;j++)
    {     
      if (left_dof[i]==left_dof[j])
      {
        found = 1;
        break;
      }
    }
    if (found)
      continue;
    found = 0;
    // find right counterpart
    for (j=leftno;j< right_curr; j++)
    {
      if ((left_y[i]==right_y[j])&&(left_x[i]==right_x[j]))
      {
         left_dof[leftno] = left_dof[i];
         left_y[leftno] = left_y[i];
         left_x[leftno] = left_x[i];
         rightno =  right_dof[leftno];
         right_dof[leftno] = right_dof[j];
         right_dof[j] = rightno;
         tmp = right_y[leftno];
         right_y[leftno] = right_y[j];
         right_y[j] = tmp;
         tmp = right_x[leftno];
         right_x[leftno] = right_x[j];
         right_x[j] = tmp;
         found=1;
         leftno++;
      }
      if (found)
        break;
    }
    if (!found)
    {
      OutPut("no counterpart found !!!" << endl);
      exit(4711);
    }
  }

  for (i=0;i<leftno;i++)
  {
    OutPut("Z " << i << " " << left_dof[i] << " " <<  left_x[i] << " " << left_y[i] << endl);
    OutPut("     " << " " << right_dof[i] << " " <<  right_x[i] << " " << right_y[i] << endl);
    if (left_y[i]<0)
      direction[i] = -1;
    if (left_y[i]>0)
      direction[i] = 1;
  }

  new_length[0] = leftno;

  delete left_y;
  delete left_x;
  delete right_y;
  delete right_x;
} 

void EnsurePeriodicity(double *sol,int N_U,int dof_length_new, 
                       int *left_dof,int *right_dof, int *direction)
{
  int i, left, right, dir;
  double tmp,max, leftdof, rightdof;
  
  max = 0;
  leftdof = 0;
  for (i=0;i<dof_length_new;i++)
  {
    left = left_dof[i];
    right = right_dof[i];
    if (fabs(sol[left]-sol[right])>max)
    {
      max = fabs(sol[left]-sol[right]);
      leftdof = left;
      rightdof = right;
      dir = direction[i];
    }
    if (direction[i]==0)
    {
      tmp = sol[left]+sol[right];
      tmp/=2;
    }
    if (direction[i]==1)
    {
      tmp = sol[right];
    }
    if (direction[i]== -1)
    {
      tmp = sol[left];
    }    
    sol[left] = tmp;
    sol[right] = tmp;
    left = left_dof[i]+N_U;
    right = right_dof[i]+N_U;
    if (fabs(sol[left]-sol[right])>max)
    {
      max = fabs(sol[left]-sol[right]);
      leftdof = left;
      rightdof = right;
      dir = direction[i];
    }
    if (direction[i]==0)
    {
      tmp = sol[left]+sol[right];
      tmp/=2;
    }
    if (direction[i]==1)
    {
      tmp = sol[right];
    }
    if (direction[i]== -1)
    {
      tmp = sol[left];
    }    
    sol[left] = tmp;
    sol[right] = tmp;
    left = left_dof[i]+2*N_U;
    right = right_dof[i]+2*N_U;
    if (fabs(sol[left]-sol[right])>max)
    {
      max = fabs(sol[left]-sol[right]);
      leftdof = left;
      rightdof = right;
      dir = direction[i];
    }
    if (direction[i]==0)
    {
      tmp = sol[left]+sol[right];
      tmp/=2;
    }
    if (direction[i]==1)
    {
      tmp = sol[right];
    }
    if (direction[i]== -1)
    {
      tmp = sol[left];
    }    
    sol[left] = tmp;
    sol[right] = tmp;
  }
  OutPut("periodicity: maximal difference " << max << " " << leftdof << " " << rightdof << 
         " " << dir << endl);
}

// thickness of spanwise vorticity
void ComputeVorticiyThickness3D(TFEFunction3D *Vorticity, double *thickness)
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
