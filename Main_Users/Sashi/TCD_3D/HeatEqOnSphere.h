
#include <TimeUtilities.h>


void ExampleFile()
{
  OutPut("Example: EvolvingEllipsoid.h" << endl) ;
}

// ========================================================================
// initial condition
// ========================================================================
void InitialU1(double x, double y,  double z, double *values)
{

  values[0] = 0;
}

void InitialU2(double x, double y,  double z, double *values)
{

  values[0] = 0;
}

void InitialU3(double x, double y,  double z, double *values)
{
  values[0] = 0;
}

void InitialP(double x, double y,  double z, double *values)
{
  values[0] = 0;
}

// ========================================================================
// exact solution
// ========================================================================

// kind of boundary condition (for FE space needed)
void BoundCondition(double x, double y, double z, BoundCond &cond)
{
     TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 0;
     cond  = FREESURF;
}

// ========================================================================
// declaration for grid handling
// ========================================================================
void GridBoundCondition(double x, double y, double z, BoundCond &cond)
{
  cond = DIRICHLET;
}

// value of boundary condition
void GridBoundValue(double x, double y, double z, double &value)
{
  value = 0;
}

void GridCoeffs(int n_points, double *x, double *y, double *Z,
        double **parameters, double **coeffs)
{
  int i;
  double *coeff;
  double r2;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    coeff[0] = 1.;
    coeff[1] = 0;
    coeff[2] = 0;
    coeff[3] = 0;
    coeff[4] = 0;
  }
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *X, double *Y, double *Z,
               double **parameters, double **coeffs)
{
  static double eps = 1.;
  int i;
  double *coeff, x, y, z;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    coeff[0] = eps;
    coeff[1] = 0;
    coeff[2] = 0;
    coeff[3] = 0;

  }
}



// kind of boundary condition (for FE space needed)
void SurfBoundCondition(int BdComp, double t, BoundCond &cond)
{
//   cond = DIRICHLET;
//  each edge of all triangles on a 3D surface is an interior edge
  cout << "Error!! each edge of all triangles on a 3D surface is an interior edge"<<endl;
  cout << "Error!! Check the Example file "<<endl;
  exit(0);

}

// value of boundary condition
void SurfBoundValue(int BdComp, double Param, double &value)
{
//   value = 0;
//  each edge of all triangles on a 3D surface is an interior edge
  cout << "Error!! each edge of all surface triangles on a 3D surface is an interior edge"<<endl;
  cout << "Error!! Check the Example file "<<endl;
  exit(0);
}


void InitialS(double x, double y, double *values)
{
   double h, r;

//  if(fabs(x)<1e-5  && fabs(y)<1e-5 )
//  cout << "InitialS : x : " << x << " y : " << y <<endl;

// if(x*y>0)
  values[0] = x*y;
// else
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;

}
void ExactS(double x, double y, double *values)
{

 double t = TDatabase::TimeDB->CURRENTTIME;
  values[0] = exp(-6.*t)*x*y;
  values[1] = exp(-6.*t)*y;
  values[2] = exp(-6.*t)*x;
  values[3] = 0;
  values[4] = 0;

}

void InitialSall(double x, double y,  double z, double *values)
{
  values[0] = x*y;
}

void ExactSall(double x, double y,  double z, double *values)
{
 double t = TDatabase::TimeDB->CURRENTTIME;
 double v, R;

  R = 1.;
  v =  sqrt(x*x + y*y + z*z) / R;
  x /= v;
  y /= v;
  z /= v;

  values[0] = x*y*exp(-6.*t);
  values[1] = exp(-6.*t)*y;
  values[2] = exp(-6.*t)*x;
  values[3] = 0;
  values[4] = 0;
}

void SurfAllCoeffs(int n_points, double *X, double *Y, double *Z,
               double **parameters, double **coeffs)
{
  static double eps = 1.;
  int i;
  double *coeff, x, y, z;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    coeff[0] = eps;
    coeff[1] = 0;
    coeff[2] = 0;
    coeff[3] = 0;

  }
}


void SurfCoeffs(int n_points, double *x, double *y,
        double **parameters, double **coeffs)
{
  int i;
  double *coeff;
  double r2;
  static double eps = 1.;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    coeff[0] = eps;   // eps

    if(TDatabase::ParamDB->FR_NR == 0)
       coeff[1] = 0;
    else
       coeff[1] = TDatabase::ParamDB->FR_NR;

  }
}

void GetRhs(double t, double x1, double x2, double x3, double &rhs)
{
 rhs =  0.;
}

void  ProjectBDPoints(TDomain *Domain)
{

  int i, j, k, l, n, N_Cells;
  double x, y, z, v, R;
  TJoint *joint;
  TShapeDesc *ShapeDesc;
  const int *TmpFV,  *TmpLen;
  int MaxLen;
  TBaseCell *Me;
  TCollection *Coll;


 R = 1.0 ;// radius of the sphere

  Coll = Domain->GetCollection(It_Finest, 0);
  N_Cells = Coll->GetN_Cells();
  for(i=0;i<N_Cells;i++)
  {
    Me = Coll->GetCell(i);
    k = Me->GetN_Faces();
    for(l=0;l<k;l++)
     {
       joint = Me->GetJoint(l);
       if(joint->GetType() == BoundaryFace || joint->GetType() == IsoBoundFace )
        {
           Me->GetShapeDesc()->GetFaceVertex(TmpFV, TmpLen, MaxLen);
           for(n=0;n<TmpLen[l];n++)
            {
              Me->GetVertex(TmpFV[l*MaxLen+n])->GetCoords(x,y,z);
              v= sqrt(x*x + y*y + z*z)/R;
//              cout<< "v "  << v <<endl;
              x /=v;
              y /=v;
              z /=v;
              Me->GetVertex(TmpFV[l*MaxLen+n])->SetCoords(x,y,z);

//              cout<< "R  "  << sqrt(x*x + y*y + z*z) <<endl;
           }
        }
     } // endfor l
  } // endfor i

}

