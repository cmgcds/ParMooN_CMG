
#include <TimeUtilities.h>


void ExampleFile()
{
  OutPut("Example: ExpandingSphere.h" << endl) ;
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
  static double eps = 1./TDatabase::ParamDB->RE_NR;
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
  cout << "Error!! Check the Example file "<<endl;
  exit(0);
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

 double t = TDatabase::TimeDB->CURRENTTIME;
//  if(fabs(x)<1e-5  && fabs(y)<1e-5 )
//  cout << "InitialS : x : " << x << " y : " << y <<endl;

  values[0] = exp(-6.*t)*x*y;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
  cout << "Error!! Check the Example file "<<endl;
  exit(0);

}
void ExactS(double x, double y, double *values)
{


  cout << "Error!! Check the Example file "<<endl;
  exit(0);

 double t = TDatabase::TimeDB->CURRENTTIME;
  values[0] = exp(-6.*t)*x*y;
  values[1] = exp(-6.*t)*y;
  values[2] = exp(-6.*t)*x;
  values[3] = 0;
  values[4] = 0;

}

void InitialSall(double x, double y,  double z, double *values)
{
 double t = TDatabase::TimeDB->CURRENTTIME;
  values[0] = x*y*exp(-6.*t);
}

void ExactSall_oldStep(double x, double y,  double z, double *values)
{
 double oldt = TDatabase::ParamDB->P14; // old time step should be storen in main prgram
 double t = TDatabase::TimeDB->CURRENTTIME;
 double a, amp=TDatabase::ParamDB->P15;

  a=1. + amp*sin(t);

  values[0] = x*y*exp(-6.*oldt)/a;
  values[1] = exp(-6.*oldt)*y/a;
  values[2] = exp(-6.*oldt)*x/a;
  values[3] = 0;
  values[4] = 0;
}


void ExactSall(double x, double y,  double z, double *values)
{
 double t = TDatabase::TimeDB->CURRENTTIME;
 double a, amp=TDatabase::ParamDB->P15;

  a=1. + amp*sin(t);

  values[0] = x*y*exp(-6.*t)/a;
  values[1] = exp(-6.*t)*y/a;
  values[2] = exp(-6.*t)*x/a;
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

  cout << "Error!! Check the Example file "<<endl;
  exit(0);

}


void SurfCoeffs(int n_points, double *x, double *y,
        double **parameters, double **coeffs)
{
  int i;
  double *coeff;
  double  r2;
  static double eps = 1./TDatabase::ParamDB->RE_NR;

  for(i=0;i<n_points;i++)
   {
    coeff = coeffs[i];
    coeff[0] = eps;   // eps

    if(TDatabase::ParamDB->FR_NR == 0)
       coeff[1] = 0;
    else
       coeff[1] = TDatabase::ParamDB->FR_NR;

   }
  cout << "Error!! Check the Example file "<<endl;
  exit(0);
}

void GetRhs(double t, double x1, double x2, double x3, double &rhs)
{
 double D, amp=TDatabase::ParamDB->P15;
 double alph=exp(-6.*t);
 double a=1.+amp*sin(t);
 double da=amp*cos(t);

 D = x1*x1 + x2*x2 + x3*x3;

// Eulerian
//  rhs = x1*x2*alph*(2.*da/a - 6. + 6./D );

// Lagrangian manner
//  rhs = x1*x2*alph*(da/a - 6. + 6./D );

// no tang-divergence term
//  rhs = x1*x2*alph*(- 6. + 6./D );

// no  tang-divergence and laplace terms
//  rhs = x1*x2*alph*(- 6.);
  rhs = (1./a - 1.)*6*x1*x2*alph/a;

}

// // only on reference sphere
// void GetUTangU(double y1, double y2, double y3, double &U1, double &U2, double &U3, double &TangDivU)
// {
//  double v, amp=TDatabase::ParamDB->P15;
//  double t=TDatabase::TimeDB->CURRENTTIME;
//  double a=1.+amp*sin(t);
//  double da=amp*cos(t);
// 
//   TangDivU =da /sqrt(a);
//   v =  da /(2.*sqrt(a));
//   U1 = v*y1;
//   U2 = v*y2;
//   U3 = v*y3;
// }

void  ProjectBDPoints(TDomain *Domain)
{

  int i, j, k, l, n, N_Cells;
  int MaxLen;
  const int *TmpFV,  *TmpLen;

  double x, y, z, v, R, t, a;
  double  amp=TDatabase::ParamDB->P15;

  TJoint *joint;
  TShapeDesc *ShapeDesc;
  TBaseCell *Me;
  TCollection *Coll;

  t = TDatabase::TimeDB->CURRENTTIME;
  a = sqrt(1. + amp*sin(t));
  R = a;// radius of the sphere

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
//               v= sqrt( (x/a)*(x/a) + (y/b)*(y/b) + (z/c)*(z/c))/R;
              v= sqrt( x*x + y*y + z*z )/R;
              x /=v;
              y /=v;
              z /=v;
//              cout<< "v "  << v <<endl;
              Me->GetVertex(TmpFV[l*MaxLen+n])->SetCoords(x,y,z);
//              cout<< "R  "  << sqrt((x/a)*(x/a) + (y/b)*(y/b) + (z/c)*(z/c)) <<endl;
           }
        }
     } // endfor l
  } // endfor i

}


// calculates the W^{n+1} of the given domain \Omega_t^{n+1}
void GetCurrentGridVelo(TFEVectFunct3D *RefGridPos, TFEVectFunct3D *GridPos, TFEVectFunct3D *AuxGridPos, 
                        TFEVectFunct3D *GridVelo, double  t, double dt,
                        TSquareMatrix3D **Mat, int *KCol, int *RowPtr,
                        double **Entries, TParDirectSolver *Grid_SMPSolver)
{
 TCollection *Coll;
 TFESpace3D *FESpace;
 TBaseCell *cell;
 FE3D FEId;
 TFE3D *Element;
 TJoint *joint;
 TIsoBoundFace *isojoint;
 TVertex **Vertices;
 TFEDesc3D *FEDesc;

 int i, j, k, l, m, n, N_U, N_Cells, N_Faces, N_LocalDOFs, *JointDOF;
 int N_Inner, N_BoundaryNodes;
 int *BeginIndex, *GlobalNumbers, *DOF;

 double x, y, z, t1, t2, Mult, t_new;
 double *ValuesV, IsoX, IsoY, IsoZ;
 double *RefX, *RefY, *RefZ;
 double *ValuesNewX, *ValuesNewY, *ValuesNewZ;
 double *RefPosX, *RefPosY, *RefPosZ;
 double *PosX, *PosY, *PosZ;
 double a, a_new, theta, phi, *sol, *rhs;
 double *Isotheta, *Isophi;
 double  amp=TDatabase::ParamDB->P15, da;

 ValuesV = GridVelo->GetValues();
 FESpace = GridVelo->GetFESpace3D(); // assume fespacces in all vecfunct is same
 Coll = FESpace->GetCollection();
 N_Cells = Coll->GetN_Cells();
 BeginIndex = FESpace->GetBeginIndex();
 GlobalNumbers = FESpace->GetGlobalNumbers();

 N_U = GridVelo->GetLength();
 N_Inner = FESpace->GetN_Inner();
 N_BoundaryNodes = N_U - N_Inner;

 PosX = GridPos->GetValues();
 PosY = PosX + N_U;
 PosZ = PosX + 2*N_U;

 RefPosX = RefGridPos->GetValues();
 RefPosY = RefPosX + N_U;
 RefPosZ = RefPosX + 2*N_U;

 ValuesNewX = AuxGridPos->GetValues();
 ValuesNewY = ValuesNewX + N_U;
 ValuesNewZ = ValuesNewX + 2*N_U;

 memcpy(ValuesNewX, RefPosX, 3*N_U*SizeOfDouble);

 t_new = t + dt;
 a_new = 1. + amp*sin(t_new);
 Mult = sqrt(a_new);

 rhs = new double[3*N_U];
 sol = new double[3*N_U];

   // determine new position of boundary vertices
  for(i=0;i<N_Cells;i++)
   {
    cell  = Coll->GetCell(i);
    N_Faces = cell->GetN_Faces();
    for(j=0;j<N_Faces;j++)
     {
      if( (cell->GetJoint(j))->GetType() == IsoBoundFace )
       {
        DOF = GlobalNumbers + BeginIndex[i];

        FEId = FESpace->GetFE3D(i, cell);
        Element = TFEDatabase3D::GetFE3D(FEId);
        N_LocalDOFs = Element->GetN_DOF();

        for(k=0;k<N_LocalDOFs;k++)
         {
          m = DOF[k];
          if(m>=N_Inner)
           {
            ValuesNewX[m] = Mult*RefPosX[m];
            ValuesNewY[m] = Mult*RefPosY[m];
            ValuesNewZ[m] = Mult*RefPosZ[m];
// only at t=0
//  cout<< "u1 "<< (ValuesNewX[m] -RefPosX[m] )/dt << 
//         " Exact " << (sqrt(a_new) - sqrt(a))*RefPosX[m]/dt <<endl;
//  cout<< "u2 "<< (ValuesNewY[m] -RefPosY[m] )/dt << 
//         " Exact " << (sqrt(a_new)*RefPosY[m] -  sqrt(a)*RefPosY[m])/dt<<endl;
//  cout<< "u3 "<< (ValuesNewZ[m] -RefPosZ[m] )/dt << 
//         " Exact " << (sqrt(a_new)*RefPosZ[m] - sqrt(a)*RefPosZ[m])/dt <<endl;
          } //  if(m>=N_Inner)
         } //    for(k=0;k<N_LocalDOFs;k++)
       }   //   if( !(cell->GetJoint(j)->InnerJoint()) )
      } // endfor j
     } //   for(i=0;i<N_Cells;i++)


  memset(sol, 0, SizeOfDouble*3*N_U);
  memset(rhs, 0, SizeOfDouble*3*N_U);
  memset(ValuesV, 0, SizeOfDouble*3*N_U);

  memcpy(sol+N_Inner, ValuesNewX+N_Inner, N_BoundaryNodes*SizeOfDouble);
  memcpy(sol+(N_U+N_Inner), ValuesNewY+N_Inner, N_BoundaryNodes*SizeOfDouble);
  memcpy(sol+(2*N_U+N_Inner), ValuesNewZ+N_Inner, N_BoundaryNodes*SizeOfDouble);

  Daxpy(N_BoundaryNodes, -1., PosX+N_Inner, sol+N_Inner);
  Daxpy(N_BoundaryNodes, -1., PosY+N_Inner, sol+(N_U+N_Inner));
  Daxpy(N_BoundaryNodes, -1., PosZ+N_Inner, sol+(2*N_U+N_Inner));

  memcpy(rhs+N_Inner, sol+N_Inner, N_BoundaryNodes*SizeOfDouble);
  memcpy(rhs+(N_U+N_Inner), sol+(N_U+N_Inner), N_BoundaryNodes*SizeOfDouble);
  memcpy(rhs+(2*N_U+N_Inner), sol+(2*N_U+N_Inner), N_BoundaryNodes*SizeOfDouble);

// solve the system
  SolveGridEquation3D(Entries, sol, rhs, KCol, RowPtr, N_U);


//    OutPut("Grid_SMPSolver start : "<< endl);
//    Grid_SMPSolver->Solve(Mat[0], Mat[1], Mat[2],
//                          Mat[3], Mat[4], Mat[5],
//                          Mat[6], Mat[7], Mat[8],
//                          sol, rhs);

  memcpy(ValuesV, sol, 3*N_U*SizeOfDouble);
  Dscal(3*N_U, 1./dt, ValuesV);

 delete [] sol;
 delete [] rhs;
}

void MoveGrid(TFEVectFunct3D *GridPos, TFEVectFunct3D *GridVelo,
              double t, double tau)
{
 TCollection *Coll;
 TFESpace3D *FESpace;
 TBaseCell *cell;
 TJoint *joint;
 TIsoBoundFace *isojoint;
 TVertex **Vertices, **RefVertices;

 int i, j, k, N_W, N_Cells, N_Faces, l;

 double Mult, IsoX, IsoY, IsoZ, theta, phi;
 double a, da, *W, *Grid_Pos;
 double *Isotheta, *Isophi;
 double  amp=TDatabase::ParamDB->P15;

 FESpace = GridVelo->GetFESpace3D(); // assume fespacces in all vecfunct is same
 Coll = FESpace->GetCollection();
 N_Cells = Coll->GetN_Cells();

  N_W = GridVelo->GetLength();
  W = GridVelo->GetValues();
  Grid_Pos = GridPos->GetValues();

  Daxpy(3*N_W, tau, W, Grid_Pos);
  GridPos->DataToGrid();


// update the isoparametric points also
  a = 1. + amp*sin(t);
  Mult = sqrt(a);

//   Mult = sqrt(1. + amp*sin(t));

  if(TDatabase::ParamDB->USE_ISOPARAMETRIC)
   {
    for(i=0;i<N_Cells;i++)
     {
      cell  = Coll->GetCell(i);
      N_Faces = cell->GetN_Faces();
      for(j=0;j<N_Faces;j++)
       {
        joint = cell->GetJoint(j);
        if( (cell->GetJoint(j))->GetType() == IsoBoundFace )
         {
          isojoint = (TIsoBoundFace *)joint;
          k = isojoint->GetN_Vertices();
          Vertices = isojoint->GetVertices();
          RefVertices = isojoint->GetRefVertices();

         for(l=0;l<k;l++)
          {
           RefVertices[l]->GetCoords(IsoX, IsoY,  IsoZ);
           IsoX = Mult*IsoX;
           IsoY = Mult*IsoY;
           IsoZ = Mult*IsoZ;

           Vertices[l]->SetCoords(IsoX, IsoY,  IsoZ);
          } // endfor l
       } // endif
    } // endfor j
  } // endfor i
 }// if(TDatabase::ParamDB->USE_ISOPARAMETRIC)
}


// void MoveCurrentDomainToRefDomain(double t, double old_t,
//                                   TFEVectFunct3D *GridPos, TFEVectFunct3D *AuxGridPos)
// {
//  TCollection *Coll;
//  TFESpace3D *FESpace;
//  TBaseCell *cell;
//  FE3D FEId;
//  TFE3D *Element;
//  TJoint *joint;
//  TIsoBoundFace *isojoint;
//  TVertex **Vertices;
//  TFEDesc3D *FEDesc;
// 
//  TCollection *Coll_low;
//  TBaseCell *Me_low;
//  FE3D Velo_FEId;
//  FE2D FEId_low;
//  TFE2D *Element_low;
// 
//  BaseFunct3D LocBF[N_BaseFuncts3D];
//  BaseFunct3D *BaseFuncts;
//  BF3DRefElements RefElement;
//  RefTrans3D RefTrans;
//  TRefTrans3D *F_K;
//  QuadFormula2D LineQuadFormula;
//  QuadFormula3D QF3;
//  TQuadFormula2D *qf;
//  QuadFormula2D QuadFormula;
//  TFEDesc3D *Velo_FeDesc, *FeDesc;
//  TFEDesc2D *FeDesc_low;
// 
// 
//  int i, j, k, l, m, n, N_U, N_Cells, N_Faces, N_LocalDOFs, *JointDOF;
//  int N_Inner, N_BoundaryNodes;
//  int *BeginIndex, *GlobalNumbers, *DOF;
//  int N_Cells_low, N, N_LocalUsedElements, local_i, local_j, ORDER;
//  int N_BaseFunct_low,  N_Points, N_JointDOF, begin, end, *N_BaseFuncts, N_BaseFunct;
//  int *BeginIndex_low, *GlobalNumbers_low, *DOF_LOW, TestDOF, AnsatzDOF, IJoint;
// 
//  double x, y, z, t1, t2, Mult, t_new;
//  double IsoX, IsoY, IsoZ;
//  double *RefX, *RefY, *RefZ;
//  double *ValuesNewX, *ValuesNewY, *ValuesNewZ;
//  double *PreviousPosX, *PreviousPosY, *PreviousPosZ;
//  double a, a_old, theta, phi;
//  double *Isotheta, *Isophi;
//  double  amp=TDatabase::ParamDB->P15, da;
//  double *ValuesM, *Weights, *p1, *p2, **uref;
//  double LocMatrixM[MaxN_BaseFunctions3D*MaxN_BaseFunctions3D];
//  double X[MaxN_QuadPoints_3D], Y[MaxN_QuadPoints_3D], Z[MaxN_QuadPoints_3D];
//  double a1, a2, a3, b1, b2, b3, n1, n2, n3, len;
//  double test000, ansatz000, val, surfacearea=0.;
// 
//  Isotheta = new double[10];
//  Isophi = new double[10];
// 
//  FESpace = GridPos->GetFESpace3D(); // assume fespacces in all vecfunct is same
//  Coll = FESpace->GetCollection();
//  N_Cells = Coll->GetN_Cells();
//  BeginIndex = FESpace->GetBeginIndex();
//  GlobalNumbers = FESpace->GetGlobalNumbers();
// 
//  N_U = GridPos->GetLength();
//  N_Inner = FESpace->GetN_Inner();
//  N_BoundaryNodes = N_U - N_Inner;
// 
//  PreviousPosX = GridPos->GetValues();
//  PreviousPosY = PreviousPosX + N_U;
//  PreviousPosZ = PreviousPosX + 2*N_U;
// 
//  ValuesNewX = AuxGridPos->GetValues();
//  ValuesNewY = ValuesNewX + N_U;
//  ValuesNewZ = ValuesNewX + 2*N_U;
// 
//  memcpy(ValuesNewX, PreviousPosX, 3*N_U*SizeOfDouble);
// 
//  a = 1. + amp*sin(t);
//  a_old =  1. + amp*sin(old_t);
//  Mult = sqrt(a_old/a);
// //  Mult = 1.;
// //  cout<< " t " << t << " old_t " << old_t << " Mult " << Mult << endl;
// 
//    // determine new position of boundary vertices
//   for(i=0;i<N_Cells;i++)
//    {
//     cell  = Coll->GetCell(i);
//     N_Faces = cell->GetN_Faces();
//     for(j=0;j<N_Faces;j++)
//      {
//       if( (cell->GetJoint(j))->GetType() == IsoBoundFace )
//        {
//         DOF = GlobalNumbers + BeginIndex[i];
// 
//         FEId = FESpace->GetFE3D(i, cell);
//         Element = TFEDatabase3D::GetFE3D(FEId);
//         N_LocalDOFs = Element->GetN_DOF();
// 
//         for(k=0;k<N_LocalDOFs;k++)
//          {
//           m = DOF[k];
//           if(m>=N_Inner)
//            {
//             ValuesNewX[m] =  Mult*PreviousPosX[m];
//             ValuesNewY[m] =  Mult*PreviousPosY[m];
//             ValuesNewZ[m] =  Mult*PreviousPosZ[m];
//            } //  if(m>=N_Inner)
//          } //    for(k=0;k<N_LocalDOFs;k++)
//        }   //   if( !(cell->GetJoint(j)->InnerJoint()) )
//       } // endfor j
//      } //   for(i=0;i<N_Cells;i++)
// 
//   AuxGridPos->DataToGrid();
// 
//   if(TDatabase::ParamDB->USE_ISOPARAMETRIC)
//    {
//     for(i=0;i<N_Cells;i++)
//      {
//       cell  = Coll->GetCell(i);
//       N_Faces = cell->GetN_Faces();
//       for(j=0;j<N_Faces;j++)
//        {
//         joint = cell->GetJoint(j);
//         if( (cell->GetJoint(j))->GetType() == IsoBoundFace )
//          {
//           isojoint = (TIsoBoundFace *)joint;
//           k = isojoint->GetN_Vertices();
//           Vertices = isojoint->GetVertices();
//           isojoint->GetParameters(Isotheta, Isophi);
// 
//          for(l=0;l<k;l++)
//           {
//            Vertices[l]->GetCoords(IsoX, IsoY,  IsoZ);
//            IsoX = Mult*IsoX;
//            IsoY = Mult*IsoY;
//            IsoZ = Mult*IsoZ;
// 
//            Vertices[l]->SetCoords(IsoX, IsoY,  IsoZ);
//           } // endfor l
//        } // endif
//     } // endfor j
//   } // endfor i
//  }// if(TDatabase::ParamDB->USE_ISOPARAMETRIC)
// 
// }
// 
// 
// void RestoreToCurrentDomain(double t, double old_t,
//                             TFEVectFunct3D *GridPos, TFEVectFunct3D *AuxGridPos)
// {
//  int i, j, k, l, m, n, N_Cells, N_Faces;
// 
//  double a, a_old, Mult, amp=TDatabase::ParamDB->P15, da;
//  double IsoX, IsoY, IsoZ;
// 
//  TCollection *Coll;
//  TFESpace3D *FESpace;
//  TBaseCell *cell;
//  TJoint *joint;
//  TIsoBoundFace *isojoint;
//  TVertex **Vertices;
// 
//  FESpace = GridPos->GetFESpace3D(); // assume fespacces in all vecfunct is same
//  Coll = FESpace->GetCollection();
//  N_Cells = Coll->GetN_Cells();
// 
// // put the domain to the original position
//  GridPos->DataToGrid();
// 
//  a = 1. + amp*sin(t);
//  a_old =  1. + amp*sin(old_t);
//  Mult =sqrt(a/a_old);
// 
// 
//   if(TDatabase::ParamDB->USE_ISOPARAMETRIC)
//    {
//     for(i=0;i<N_Cells;i++)
//      {
//       cell  = Coll->GetCell(i);
//       N_Faces = cell->GetN_Faces();
//       for(j=0;j<N_Faces;j++)
//        {
//         joint = cell->GetJoint(j);
//         if( (cell->GetJoint(j))->GetType() == IsoBoundFace )
//          {
//           isojoint = (TIsoBoundFace *)joint;
//           k = isojoint->GetN_Vertices();
//           Vertices = isojoint->GetVertices();
// 
//          for(l=0;l<k;l++)
//           {
//            Vertices[l]->GetCoords(IsoX, IsoY,  IsoZ);
//            IsoX = Mult*IsoX;
//            IsoY = Mult*IsoY;
//            IsoZ = Mult*IsoZ;
// 
//            Vertices[l]->SetCoords(IsoX, IsoY,  IsoZ);
//           } // endfor l
//        } // endif
//     } // endfor j
//   } // endfor i
//  }// if(TDatabase::ParamDB->USE_ISOPARAMETRIC)
// 
// 
// // cout<< " surfacearea " << surfacearea <<endl;
// }
// 









