// =======================================================================
// 
// Purpose:     Main program for rigid body motion
//
// Author   :     Bhanu Teja, Sashikumaar Ganesan
//  date    :     29.05.2014
// modified :   
// ======================================================================= 
#include <Domain.h>
#include <FEDatabase3D.h>
#include <SquareStructure3D.h>
#include <Structure3D.h>
#include <AuxParam3D.h>
#include <DirectSolver.h>
#include <Assemble3D.h>
#include <Output3D.h>
#include <DiscreteForm3D.h>
#include <QuadAffin.h>
#include <QuadBilinear.h>
#include <QuadIsoparametric.h>
#include <TetraAffin.h>
#include <TetraIsoparametric.h>
#include <HexaAffin.h>
#include <HexaIsoparametric.h>
#include <HexaTrilinear.h>
#include <TriaAffin.h>
#include <TriaIsoparametric.h>
#include <IsoBoundEdge.h>
#include <IsoInterfaceJoint.h>
#include <FESpace3D.h>
#include <Database.h>
#include <string.h>
#include <stdlib.h>
// #include <malloc.h>
#include <math.h>
#include <stdio.h>
#include <LinAlg.h>
#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <sys/types.h>
#include <BoundFace.h>
#include <InterfaceJoint3D.h>
#include <tetgen.h>
#include <GridCell.h>
#include <MacroCell.h>
#include <BdPlane.h>
#include <BdSphere.h>
#include <IsoBoundFace.h>
#include <InterfaceJoint3D.h>
#include <IsoInterfaceJoint3D.h>

#include <LinAlg.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

// #define PI 3.14159265
// ============================================================
//  include example file
// #include "../../Examples/TNSE_3D/RotatingBody.h"
#include "../../Examples/RB_3D/aerofoil_expt.h"
// ============================================================

extern "C" {

void dgetrf_(int *m, int *n, double *A, int *lda, int *ipiv, int *info);

void dgetrs_(char *trans, int *n, int *nrhs, double *A, int *lda, int *ipiv,
             double *B, int *ldb, int *info);

void dsptrd_(char *uplo, int *n,  double *AP, double *D, double *E, 
             double *TAU, int *info);

void dopgtr_(char *uplo, int *n,  double *AP, double *TAU, double *Q, int *LDQ, double *work,int *info);

void  dsteqr_(char *compz, int *N, double *D,  double *E, double *Z, int *LDZ, double *work, int *info);
 
}
// GetMOI(Coll, gridpos, CGx_Body, CGy_Body, CGz_Body, MOItensor, mass);
void GetMOI(TCollection *Coll, double *gridpos, int N_GridDOFs, double &CGx_Body, double &CGy_Body, double &CGz_Body, double *MOItensor, double &mass)
{
  int i,j,k;
  int N_Cells, N_Joints;
  TBaseCell *cell;
  TJoint *joint;
  JointType jointtype;
  boolean IsIsoparametric;
  RefTrans3D RefTrans;
  TRefTrans3D *rt;
  QuadFormula3D QuadFormula;
  TQuadFormula3D *qf2;
  int N_Points;
  double *weights, *xi, *eta, *zeta, Izz=0., Iyy=0., Ixx=0., Ixy=0., Iyz=0., Izx=0.;
  double Izzcell, Iyycell, Ixxcell, Ixycell, Iyzcell, Izxcell;
  double absdetjk[MaxN_QuadPoints_3D];
  double X[MaxN_QuadPoints_3D], Y[MaxN_QuadPoints_3D], Z[MaxN_QuadPoints_3D];
  double vol, locvol;
  int MaxApproxOrder = 2;
  
  CGx_Body = CGy_Body = CGz_Body = 0;
  
    for(i=0;i<N_GridDOFs;i++)
  {
    CGx_Body +=gridpos[i];
    CGy_Body +=gridpos[i+N_GridDOFs];
    CGz_Body +=gridpos[i+2*N_GridDOFs];
  }
  CGx_Body /=(double)N_GridDOFs;
  CGy_Body /=(double)N_GridDOFs; 
  CGz_Body /=(double)N_GridDOFs;  // end of CG 
  cout <<"CGx is " <<CGx_Body <<endl;
  cout <<"CGy is " <<CGy_Body <<endl;
  cout <<"CGz is " <<CGz_Body <<endl;
    
  N_Cells = Coll->GetN_Cells();
    for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    N_Joints = cell->GetN_Joints();
    IsIsoparametric = FALSE;
    for(j=0;j<N_Joints;j++)
    {
      joint = cell->GetJoint(j);
      switch(cell->GetType())
      {
	case Tetrahedron:
	  RefTrans = TetraAffin;
	  QuadFormula = TFEDatabase3D::GetQFTetraFromDegree(3);
	break;
	
	case Brick:
	  RefTrans = HexaAffin;
	  QuadFormula = TFEDatabase3D::GetQFHexaFromDegree(2);
	break;

	case Hexahedron:
	  RefTrans = HexaTrilinear;
	  QuadFormula = TFEDatabase3D::GetQFHexaFromDegree(2);
	break;
      } // endswitch GetType

//       jointtype = joint->GetType();
//       if(jointtype == BoundaryFace)
//       {
//         if( ((TBoundFace*)joint)->GetBoundComp()->GetType() != Plane)
// // 	if( (joint)->GetBoundComp()->GetType() != Plane)
//           IsIsoparametric = TRUE;
//       }
// 
//       if(jointtype == InterfaceJoint3D)
//       {
//         if( ((TInterfaceJoint3D*)joint)->GetBoundComp()->GetType() != Plane)
// // 	if( (joint)->GetBoundComp()->GetType() != Plane)
//           IsIsoparametric = TRUE;
//       }
// 
//       if(jointtype == IsoInterfaceJoint3D)
//         IsIsoparametric = TRUE;
// 
//       if(jointtype == IsoBoundFace)
//         IsIsoparametric = TRUE;
// 
//       if(jointtype == IsoJointEqN)
//         IsIsoparametric = TRUE;
    } // endfor j

    if(IsIsoparametric)
    {
      switch(RefTrans)
      {
	case HexaAffin:

	case HexaTrilinear:
	RefTrans = HexaIsoparametric;
	break;

	case TetraAffin:
	RefTrans = TetraIsoparametric;
	break;
      }                                                     // endswitch
    }                                                       // endif IsIsoparametric

    qf2 = TFEDatabase3D::GetQuadFormula3D(QuadFormula);
    qf2->GetFormulaData(N_Points, weights, xi, eta, zeta);

    rt = TFEDatabase3D::GetRefTrans3D(RefTrans);
    switch(RefTrans)
    {
      case TetraAffin:
        // cout << "TetraAffin" << endl;
        ((TTetraAffin *)rt)->SetCell(cell);
        ((TTetraAffin *)rt)->GetOrigFromRef(N_Points, xi, eta, zeta, X, Y, Z, absdetjk);
      break;
//       case TetraIsoparametric:
//         cout << "TetraIsoparametric" << endl;
//         ((TTetraIsoparametric *)rt)->SetApproximationOrder(MaxApproxOrder);
//         ((TTetraIsoparametric *)rt)->SetQuadFormula(QuadFormula);
//         ((TTetraIsoparametric *)rt)->SetCell(cell);
//         ((TTetraIsoparametric *)rt)->GetOrigFromRef(N_Points, xi, eta, zeta, X, Y, Z, absdetjk);
//       break;
      case HexaAffin:
        // cout << "HexaAffin" << endl;
        ((THexaAffin *)rt)->SetCell(cell);
        ((THexaAffin *)rt)->GetOrigFromRef(N_Points, xi, eta, zeta, X, Y, Z, absdetjk);
      break;
      case HexaTrilinear:
        // cout << "HexaTrilinear" << endl;
        ((THexaTrilinear *)rt)->SetCell(cell);
        ((THexaTrilinear *)rt)->GetOrigFromRef(N_Points, xi, eta, zeta, X, Y, Z, absdetjk);
      break;
//       case HexaIsoparametric:
//         // cout << "HexaIsoparametric" << endl;
//         ((THexaIsoparametric *)rt)->SetApproximationOrder(MaxApproxOrder);
//         ((THexaIsoparametric *)rt)->SetQuadFormula(QuadFormula);
//         ((THexaIsoparametric *)rt)->SetCell(cell);
//         ((THexaIsoparametric *)rt)->GetOrigFromRef(N_Points, xi, eta, zeta, X, Y, Z, absdetjk);
//       break;
      
      default:
        cout<<"Wrong transformation, Isoparameteric not yet implemented/tested!! "<<endl;
	exit(0);
      break;      
    } // endswitch
    
      locvol = 0., Izzcell = 0., Iyycell = 0., Ixxcell = 0.; Ixycell = 0.; Iyzcell = 0.; Izxcell = 0.;
    for(j=0;j<N_Points;j++)
    {
      locvol += weights[j]*absdetjk[j];
      Izzcell+= weights[j]*absdetjk[j]*(pow((X[j] - CGx_Body),2) + pow((Y[j] - CGy_Body),2));
      Iyycell+= weights[j]*absdetjk[j]*(pow((X[j] - CGx_Body),2) + pow((Z[j] - CGz_Body),2));
      Ixxcell+= weights[j]*absdetjk[j]*(pow((Y[j] - CGy_Body),2) + pow((Z[j] - CGz_Body),2));
      Ixycell+= weights[j]*absdetjk[j]*(X[j] - CGx_Body)*(Y[j] - CGy_Body);
      Iyzcell+= weights[j]*absdetjk[j]*(Y[j] - CGy_Body)*(Z[j] - CGz_Body);
      Izxcell+= weights[j]*absdetjk[j]*(Z[j] - CGz_Body)*(X[j] - CGx_Body);      
    }
   mass +=locvol;
   Izz+= Izzcell;
   Iyy+= Iyycell;
   Ixx+= Ixxcell;
   Ixy+= Ixycell;
   Iyz+= Iyzcell;
   Izx+= Izxcell;   
     } // end of i < N_Cells

  //multiply with the density
  Izz *=TDatabase::ParamDB->P0;
  Iyy *=TDatabase::ParamDB->P0;
  Ixx *=TDatabase::ParamDB->P0;
  Ixy *=TDatabase::ParamDB->P0;
  Iyz *=TDatabase::ParamDB->P0;
  Izx *=TDatabase::ParamDB->P0;
  mass*=TDatabase::ParamDB->P0;
     cout<< "Izz is "<< Izz<<endl;
     cout<< "Iyy is "<< Iyy<<endl;
     cout<< "Ixx is "<< Ixx<<endl;
     cout<< "Ixy is "<< Ixy<<endl;
     cout<< "Iyz is "<< Iyz<<endl;
     cout<< "Izx is "<< Izx<<endl;
     cout<< "mass is "<< mass<<endl;
//      MOItensor[0] = Izz; MOItensor[1] = -Iyz; MOItensor[2] = Iyy; //packed MOItensor columnwise for LinAlg.C
//      MOItensor[3] = -Izx; MOItensor[4] = -Ixy; MOItensor[5] = Ixx; 
     MOItensor[0] = Ixx; MOItensor[1] = -Ixy; MOItensor[2] = Iyy; //packed MOItensor columnwise for LinAlg.C
     MOItensor[3] = -Izx; MOItensor[4] = -Iyz; MOItensor[5] = Izz; 
} // end of GetMOI


// UpdateMOIs(N_Cells, CGxyz, CGx_Body, CGy_Body, CGz_Body, locvol, Izzt, Iyyt, Ixxt);
void UpdateMOIs(int N_Cells, double *CGxyz, double CGx_Body, double CGy_Body, double CGz_Body, double locvol, double &Izzt, double &Iyyt, double &Ixxt)
{
  int i;
  double x2y2, z2x2, y2z2, Izz, Iyy, Ixx;
  for(i=0;i<N_Cells;i++)
  {  
   x2y2 = (CGxyz[i] - CGx_Body)*(CGxyz[i] - CGx_Body) + (CGxyz[i+N_Cells] - CGy_Body)*(CGxyz[i+N_Cells] - CGy_Body);
   z2x2 = (CGxyz[i+2*N_Cells] - CGz_Body)*(CGxyz[i+2*N_Cells] - CGz_Body) + (CGxyz[i] - CGx_Body)*(CGxyz[i] - CGx_Body);
   y2z2 = (CGxyz[i+N_Cells] - CGy_Body)*(CGxyz[i+N_Cells] - CGy_Body) + (CGxyz[i+2*N_Cells] - CGz_Body)*(CGxyz[i+2*N_Cells] - CGz_Body);   
   Izz += locvol*x2y2;
   Iyy += locvol*z2x2;
   Ixx += locvol*y2z2;
  }
    //multiply with the density
  Izz *=TDatabase::ParamDB->P0;
  Iyy *=TDatabase::ParamDB->P0;
  Ixx *=TDatabase::ParamDB->P0;  
  Izzt = Izz; Iyyt = Iyy; Ixxt = Ixx;
}

void vecscal(int n, double alpha, double *x, double *y)
{
  register int i;
  register double scal;
  register double *a;

  scal = alpha;
  a = x;
  for(i=0; i<n; i++)
  {
    *y = *a*scal;
    a++;
    y++;
  }
}
/** x := alpha*x */
void vecscal(int n, double alpha, double *x)
{
  register int i;
  register double scal;
  register double *a;

  scal = alpha;
  a = x;
  for(i=0; i<n; i++)
  {
    *a *=scal;
    a++;
  }
}
/** x := alpha+x */
void vecadd(int n, double alpha, double *x)
{
  register int i;
  register double scal;
  register double *a;

  scal = alpha;
  a = x;
  for(i=0; i<n; i++)
  {
    *a += scal;
    a++;
  }
}
/** y := alpha+x */
void vecadd(int n, double alpha, double *x, double *y)
{
  register int i;
  register double scal;
  register double *a, *b;

  scal = alpha;
  a = x;
  b = y;
  for(i=0; i<n; i++)
  {
    *b = *a + scal;
    a++;
    b++;
  }
}
/** y := alpha*x + y */
void vecsum(int n, double alpha, double *x, double *y)
{
  register int i;
  register double *a, *b;
  register double scal;

  a = x;
  b = y;
  scal = alpha;
  for(i=0;i<n;i++)
  {
    *b += scal * *a;
    a++;
    b++;
  }
}
/** y := x + y */
void vecsum(int n, double *x, double *y)
{
  register int i;
  register double *a, *b;

  a = x;
  b = y;
  for(i=0;i<n;i++)
  {
    *b += *a;
    a++;
    b++;
  }
}
// z = x + y
void vecsum(int n, double *x, double *y, double *z)
{
  register int i;
  register double *a, *b, *c;

  a = x;
  b = y;
  c = z;
  for(i=0;i<n;i++)
  {
    *c = *a + *b;
    a++;
    b++;
    c++;
  }
}
//  generate F vector given the input U vector
//  generateFvec(Coll, finputkn, Fvecn, Ixx, Iyy, Izz, mass, Un0);
void generateFvec(TCollection *Coll, double *ffinputkn, double *fFvec, double fIxx, double fIyy, double fIzz,  double fmass, double *fUn0)
{
  double CGx, CGy, CGz, totFx=0, totFy=0, totFz=0, totmomentx=10, totmomenty=10, totmomentz=10, *inpforCGx, *inpforCGy, *inpforCGz;
  
  fFvec[0] = ffinputkn[6];
  fFvec[1] = ffinputkn[7];
  fFvec[2] = ffinputkn[8];
  fFvec[3] = ffinputkn[9];
  fFvec[4] = ffinputkn[10];
  fFvec[5] = ffinputkn[11];

  fFvec[6] = totFx/fmass;
  fFvec[7] = totFy/fmass;
  fFvec[8] = totFz/fmass;
  fFvec[9] = (totmomentx - (fIzz - fIyy)*ffinputkn[10]*ffinputkn[11])/fIxx;
  fFvec[10] = (totmomenty - (fIxx - fIzz)*ffinputkn[11]*ffinputkn[9])/fIyy;
  fFvec[11] = (totmomentz - (fIyy - fIxx)*ffinputkn[9]*ffinputkn[10])/fIzz; 
  
//     int i;  
//   for(i=0;i<12;i++)
//     cout <<"rhs["<<i<<"] "<<fFvec[i]<<endl; 
//   
//   cout<<endl;
//   exit(0);
}
// FixtoPaxes(gridposPaxes, gridpos_CM, Paxes);
void FixtoPaxes(double *gridposPaxes, double *gridpos_CM, double *Paxes, int N_GridDOFs)
{
  int i;
  double x,y,z, xproj=0, yproj=0, zproj=0;
  for(i=0; i<N_GridDOFs; i++)
  {
    x = gridpos_CM[i];
    y = gridpos_CM[i+N_GridDOFs];
    z = gridpos_CM[i+2*N_GridDOFs];
    
    xproj = x*Paxes[0] + y*Paxes[1] + z*Paxes[2];
    yproj = x*Paxes[3] + y*Paxes[4] + z*Paxes[5];
    zproj = x*Paxes[6] + y*Paxes[7] + z*Paxes[8];
    
    gridposPaxes[i] = xproj;
    gridposPaxes[i+N_GridDOFs] = yproj;
    gridposPaxes[i+2*N_GridDOFs] = zproj;
  }  
}
// stdAxeswrtPaxes(Paxes, stdaxeswrtpaxes);
void stdAxeswrtPaxes(double *Paxes, double *stdaxeswrtpaxes)
{
  stdaxeswrtpaxes[0] = Paxes[0];
  stdaxeswrtpaxes[1] = Paxes[3];
  stdaxeswrtpaxes[2] = Paxes[6];
  stdaxeswrtpaxes[3] = Paxes[1];
  stdaxeswrtpaxes[4] = Paxes[4];
  stdaxeswrtpaxes[5] = Paxes[7];
  stdaxeswrtpaxes[6] = Paxes[2];
  stdaxeswrtpaxes[7] = Paxes[5];
  stdaxeswrtpaxes[8] = Paxes[8];
//   for(int i=0;i<9;i++)
//   cout << "stdaxeswrtpaxes["<< i<<"] "<< stdaxeswrtpaxes[i]<<endl;
    
  double norm1 = pow(Paxes[0]*Paxes[0] + Paxes[3]*Paxes[3] + Paxes[6]*Paxes[6],0.5);
  double norm2 = pow(Paxes[1]*Paxes[1] + Paxes[4]*Paxes[4] + Paxes[7]*Paxes[7],0.5);
  double norm3 = pow(Paxes[2]*Paxes[2] + Paxes[5]*Paxes[5] + Paxes[8]*Paxes[8],0.5);
  
  vecscal(3, 1./norm1, stdaxeswrtpaxes);
  vecscal(3, 1./norm2, stdaxeswrtpaxes+3);
  vecscal(3, 1./norm3, stdaxeswrtpaxes+6);
  cout<<endl;
}

// shifttostdaxes(rotvecPaxes, stdaxeswrtpaxes, rotvecstd);
void shifttostdaxes(double *rotvecPaxes, double *stdaxeswrtpaxes, double *rotvecstd)
{
  rotvecstd[0] = rotvecPaxes[0]*stdaxeswrtpaxes[0] + rotvecPaxes[1]*stdaxeswrtpaxes[1] + rotvecPaxes[2]*stdaxeswrtpaxes[2];
  rotvecstd[1] = rotvecPaxes[0]*stdaxeswrtpaxes[3] + rotvecPaxes[1]*stdaxeswrtpaxes[4] + rotvecPaxes[2]*stdaxeswrtpaxes[5];
  rotvecstd[2] = rotvecPaxes[0]*stdaxeswrtpaxes[6] + rotvecPaxes[1]*stdaxeswrtpaxes[7] + rotvecPaxes[2]*stdaxeswrtpaxes[8];
}

//RotatePaxes(Paxes, rotvecstd);
void RotatePaxes(double* Paxes, double* rotvecstd)
{
    int j;
    double newPaxes[9];
    double alpha1=rotvecstd[0], beta1=rotvecstd[1], gamma1=rotvecstd[2];
    for(j=0;j<3;j++)  // rotation of x, y and z coordinates of grid
    {           
      newPaxes[3*j] = Paxes[3*j]*cos(gamma1)*cos(beta1) - Paxes[3*j+1]*sin(gamma1)*cos(beta1) + Paxes[3*j+2]*sin(beta1);
      newPaxes[3*j+1] = Paxes[3*j]*(sin(gamma1)*cos(alpha1) + sin(alpha1)*sin(beta1)*cos(gamma1)) + Paxes[3*j+1]*(cos(alpha1)*cos(gamma1) - sin(gamma1)*sin(beta1)*sin(alpha1)) - Paxes[3*j+2]*sin(alpha1)*cos(beta1);
      newPaxes[3*j+2] = Paxes[3*j]*(sin(alpha1)*sin(gamma1) - cos(alpha1)*sin(beta1)*cos(gamma1)) + Paxes[3*j+1]*(cos(gamma1)*sin(alpha1) + cos(alpha1)*sin(beta1)*sin(gamma1)) + Paxes[3*j+2]*cos(beta1)*cos(alpha1);         
    }
    memcpy(Paxes, newPaxes, 9*sizeof(double));
}

/** calculate the eigenvalue of the system using Lapack routines*/
void FindEigenValuesss(double *ap, int N_Eqn, char &COMPZ, double *d, double *z)
{
// Arguments:
//  ap         double precision array which contains the packed upper triangular matrix column wise
//              a[i,j] = a[i +(j-1)*j/2]
//  UPLO    (input) CHARACTER*1
//          = 'U':  Upper triangle of A is packed;
//          = 'L':  Lower triangle of A is packed./**/
//  N_Eqn    : order of the matrix ap 
//  COMPZ   (input) CHARACTER*1
//           = 'N':  Compute eigenvalues only.
//           = 'V':  Compute eigenvalues and eigenvectors of the original
//                   symmetric matrix.  On entry, Z must contain the
//                   orthogonal matrix used to reduce the original matrix
//                   to tridiagonal form.
//           = 'I':  Compute eigenvalues and eigenvectors of the
//                   tridiagonal matrix.  Z is initialized to the identity
//                   matrix.
//  d       (output) DOUBLE PRECISION array, dimension (LDZ, N_Eqn) eigen values
//  z       (output) DOUBLE PRECISION array, dimension (LDZ, N)
//           COMPZ = 'V', Z contains the
//           orthonormal eigenvectors of the original symmetric matrix,
//           and if COMPZ = 'I', Z contains the orthonormal eigenvectors
//           of the symmetric tridiagonal matrix.
//           If COMPZ = 'N', then Z is not referenced.


 int info, i;
 char uplo='U';
 char compz='N';
//  double *d = new double[N_Eqn];
 double *e = new double[N_Eqn-1];
 double *tau = new double[N_Eqn-1];
 double *q = new double[N_Eqn*N_Eqn];
 double *work =  new double[2*N_Eqn-1];
 
 dsptrd_(&uplo, &N_Eqn,  ap, d, e, tau, &info);
  
 dopgtr_(&uplo, &N_Eqn,  ap, tau, q, &N_Eqn, work, &info);

 dsteqr_(&COMPZ, &N_Eqn, d, e, q, &N_Eqn, work, &info);

 memcpy(z, q, 9*sizeof(double));
 
  for(i=0; i<N_Eqn; i++)
   cout<< " Eigen values " << d[i] << endl;
  
  for (i=0;i<N_Eqn*N_Eqn;i++)
    cout<< "q["<<i<<"]="<<q[i]<<endl;
//  
}

//===========================================================
int main(int argc, char* argv[])
{ 
  TDomain *Domain = new TDomain();
  TDatabase *Database = new TDatabase();
  TFEDatabase3D *FEDatabase = new TFEDatabase3D();  
  TCollection *Coll;
  TBaseCell *Cell;
  TFESpace3D *Grid_Space, *Cell_Space;
  TFEVectFunct3D *GridPos;
  TOutput3D *Output;
  
  int i, j, k, l, N_Cells, ret, img=1;
  int N_GridDOFs, m=0;  
  double *gridpos, *gridpos_CM, *gridposPaxes;
  double CGx_Body, CGy_Body, CGz_Body, Ixx=0, Iyy=0, Izz=0, mass;
//   double locvol = 5.0;
  char *PRM, *GEO;
  char ReadinDat[] = "readin.dat";
  char GridString[] = "X";
  char Description[] = "description";
  char NameString[] = "VMS";
  char WString[] = "w";
  const char vtkdir[] = "VTK";
  char *VtkBaseName;
  std::ostringstream os;
  os << " ";
//======================================================================
// read parameter file
//======================================================================
   if(argc>=2)
    { ret=Domain->ReadParam(argv[1]); }
    else
      { ret=Domain->ReadParam(ReadinDat); }  
  if(ret==-1)
   {
    exit(-1);
   }
  OpenFiles();
  OutFile.setf(std::ios::scientific);
  Database->WriteParamDB(argv[0]);
  Database->WriteTimeDB();
  ExampleFile();
  
  VtkBaseName = TDatabase::ParamDB->VTKBASENAME;
  mkdir(vtkdir, 0777);
//======================================================================
// read boundary parameterization and initialize coarse grid
//======================================================================
   PRM = TDatabase::ParamDB->BNDFILE;
   GEO = TDatabase::ParamDB->GEOFILE;
   Domain->Init(PRM, GEO);
//     TetraGen(Domain);
   for(i=0;i<TDatabase::ParamDB->UNIFORM_STEPS;i++)
     Domain->RegRefineAll();
   Coll = Domain->GetCollection(It_Finest, 0);
//    N_Cells = Coll->GetN_Cells();
   //=========================================================================
   // construct all finite element spaces
   //=========================================================================
//                       TFESpace3D(coll, NameString, WString, GridBoundCondition, 1);  
   Grid_Space  =  new TFESpace3D(Coll, NameString, WString, GridBoundCondition, 1);   //1 -> p1 finite element space
   N_GridDOFs = Grid_Space->GetN_DegreesOfFreedom();  
   //=========================================================================
   // memory allocate all vectors and construction of all fefunction
   //=========================================================================      
   gridpos = new double[3*N_GridDOFs];  //original coordinates of the grid(mesh)
   gridpos_CM = new double[3*N_GridDOFs]; //coordinates of the mesh WRT centre of mass
//    gridposPaxes = new double[3*N_GridDOFs];
   
   memset(gridpos, 0, 3*N_GridDOFs*SizeOfDouble);
   GridPos = new TFEVectFunct3D(Grid_Space, GridString, GridString, gridpos, N_GridDOFs, 3); 
   GridPos->GridToData();
    //=========================================================================
  Output = new TOutput3D(1, 1, 1, 1, Domain, Coll, vtkdir);
  Output->AddFEVectFunct(GridPos);
  os.seekp(std::ios::beg);
  Output->AddParameter(TDatabase::TimeDB->CURRENTTIME,os.str().c_str());
//=========================================================================
  //rigid body location at T = 0
     if(TDatabase::ParamDB->WRITE_VTK)
      {
       os.seekp(std::ios::beg);
        if(img<10) os <<  "VTK/"<<VtkBaseName<<".0000"<<img<<".vtk" << ends;
         else if(img<100) os <<  "VTK/"<<VtkBaseName<<".000"<<img<<".vtk" << ends;
          else if(img<1000) os <<  "VTK/"<<VtkBaseName<<".00"<<img<<".vtk" << ends;
           else if(img<10000) os <<  "VTK/"<<VtkBaseName<<".0"<<img<<".vtk" << ends;
            else  os <<  "VTK/"<<VtkBaseName<<"."<<img<<".vtk" << ends;
       Output->WriteVtk(os.str().c_str());
       img++;
      }    

   memcpy(gridpos_CM, gridpos, 3*N_GridDOFs*SizeOfDouble);  //gridpos_CM for gridpos WRT CM. changed in later step (line 418, 419,420)

    double MOItensor[6]; // six independent components of MOI tensor
    double PMOI[3]; //for principle moments of inertia
    double Paxes[9]; //for principle axes
    double stdaxeswrtpaxes[9]; // standard axes WRT priciple axes
   
    memset(MOItensor, 0, 6*SizeOfDouble);
    memset(PMOI, 0, 3*SizeOfDouble);
    memset(Paxes, 0, 9*SizeOfDouble);
    memset(stdaxeswrtpaxes, 0, 9*SizeOfDouble);  //standard axes unit vectors WRT principle axes
    char COMPZ = 'V'; // needed for FindEigenValues
// void GetMOI(TCollection *Coll, double *CGxyz, double &CGx_Body, double &CGy_Body, double &CGz_Body, double *MOItensor, double &mass)
    GetMOI(Coll, gridpos, N_GridDOFs, CGx_Body, CGy_Body, CGz_Body, MOItensor, mass); // returns CGx_Body, CGy_Body, CGz_Body, MOItensor and mass(density of solid phase is taken from DAT file)
// void FindEigenValues(double *ap, int N_Eqn, char &COMPZ, double *d, double *z)
    FindEigenValuesss(MOItensor, 3, COMPZ, PMOI, Paxes);
    for(i=0;i<9;i++)
      cout <<"Paxes["<<i<<"] "<<Paxes[i]<<endl;
    for(i=0;i<3;i++)
      cout <<"PMOI["<<i<<"] "<<PMOI[i]<<endl;

  Ixx = PMOI[0]; Iyy = PMOI[1]; Izz = PMOI[2];
  // x,y,z of grid WRT CM
  vecadd(N_GridDOFs, -CGx_Body, gridpos_CM);              // Remove CMX from grid x
  vecadd(N_GridDOFs, -CGy_Body, gridpos_CM+N_GridDOFs);   // Remove CMY from grid y
  vecadd(N_GridDOFs, -CGz_Body, gridpos_CM+2*N_GridDOFs); // Remove CMZ from grid z
  //now gridpos_CM is ready
  memcpy(gridpos, gridpos_CM, 3*N_GridDOFs*SizeOfDouble);
  double x=0, y=0, z=0, xdot=0, ydot=0, zdot=0, alpha = 0, beta = 0, gamma = 0, alphadot=0, betadot=0, gammadot=0;
//   around x = alpha; around y = beta ;around z = gamma;
  double k1[12], k2[12], k3[12], k4[12], Un0[12], Un1[12];
  Un0[0] = x; Un0[1] = y; Un0[2] = z; Un0[3] = alpha; Un0[4] = beta; Un0[5] = gamma;
  Un0[6] = xdot; Un0[7] = ydot; Un0[8] = zdot; Un0[9] = alphadot; Un0[10] = betadot; Un0[11] = gammadot;
  double ksum[12], finputk1[12], Fvec1[12], finputk2[12], Fvec2[12], finputk3[12], Fvec3[12], finputk4[12], Fvec4[12];
  memset(ksum, 0, 12*SizeOfDouble);
  memset(finputk1, 0, 12*SizeOfDouble);
  memset(finputk2, 0, 12*SizeOfDouble);
  memset(finputk3, 0, 12*SizeOfDouble);
  memset(finputk4, 0, 12*SizeOfDouble);
  double *rotgridx, *rotgridy, *rotgridz, rotvecPaxes[3], rotvecstd[3];
  memset(rotvecPaxes, 0, 3*SizeOfDouble);
  memset(rotvecstd, 0, 3*SizeOfDouble);
  rotgridx = new double[N_GridDOFs];
  rotgridy = new double[N_GridDOFs];
  rotgridz = new double[N_GridDOFs];
  double t = 0, deltat, timeend;  
  double alpha1=0, beta1=0, gamma1=0;
  double newPaxes[9];
  double dispx=0, dispy=0, dispz=0;
  memset(newPaxes, 0, 9*SizeOfDouble);
  deltat = TDatabase::TimeDB->TIMESTEPLENGTH;
  timeend = TDatabase::TimeDB->ENDTIME;

//=========================================================================   
  //time stepping and RK4
  while(t<timeend)
  {
   m++;
   TDatabase::TimeDB->CURRENTTIME += deltat;
   
   OutPut(endl << "CURRENT TIME: ");
   OutPut(TDatabase::TimeDB->CURRENTTIME << endl);  

    memcpy(finputk1, Un0, 12*sizeof(double));
    generateFvec(Coll, finputk1, Fvec1, Ixx, Iyy, Izz, mass, Un0);
    
    vecscal(12, deltat, Fvec1);  // Fvec1 is now k1
    memcpy(ksum, Fvec1, 12*sizeof(double)); // ksum has k1 now
    vecscal(12, 0.5, Fvec1); //Fvec1 is now k1/2
    memcpy(finputk2, Fvec1, 12*sizeof(double)); //finputk2 has k1/2 (Fvec1) now
    vecsum(12, Un0, finputk2); //finputk2 has Un0 + k1/2
    generateFvec(Coll, finputk2, Fvec2, Ixx, Iyy, Izz, mass, Un0);

    vecscal(12, deltat, Fvec2); // Fvec2 is now k2
    memcpy(finputk3, Fvec2, 12*sizeof(double)); //finputk3 is k2 (Fvec2) now
    vecscal(12, 2, Fvec2); //Fvec2 is now 2*k2
    vecsum(12, Fvec2, ksum); //ksum has k1 + 2*k2    
    vecscal(12, 0.5, finputk3); //finputk3 has k2/2 now
    vecsum(12, Un0, finputk3); //finputk3 is now Un0 + k2/2
    generateFvec(Coll, finputk3, Fvec3, Ixx, Iyy, Izz, mass, Un0);

    vecscal(12, deltat, Fvec3); // Fvec3 is now k3
    memcpy(finputk4, Fvec3, 12*sizeof(double)); //finputk4 has k3 (Fvec3) now
    vecscal(12, 2, Fvec3); //Fvec3 is now 2*k3
    vecsum(12, Fvec3, ksum); //ksum has k1 + 2*k2 + 2*k3 
    vecsum(12, Un0, finputk4); //finputk4 is now Un0 + k3
    generateFvec(Coll, finputk4, Fvec4, Ixx, Iyy, Izz, mass, Un0);
    
    vecscal(12, deltat, Fvec4); // Fvec4 is now k4
    vecsum(12, Fvec4, ksum); //ksum has k1 + 2*k2 + 2*k3 + k4
    vecscal(12, 1.0/6.0, ksum); // ksum is now (k1 + 2*k2 + 2*k3 + k4)/6
    vecsum(12, Un0, ksum); // ksum is now Un0 + (k1 + 2*k2 + 2*k3 + k4)/6

    memcpy(Un1, ksum, 12*sizeof(double));
    
//     exit(0);
// ===================================================================    
//     Un1[0] is x; Un1[1] is y; Un1[2] is z; Un1[3] is alpha; Un1[4] is beta; Un1[5] is gamma;
//     alpha is around x axis(Ixx); beta is around y axis(Iyy); gamma is around z axis(Izz);
    
      
    memcpy(rotvecPaxes, Un1+3, 3*sizeof(double)); // copied rotation solution i..e \delta alpha, beta, gamma to rotvecPaxes. rotation vector WRT Paxes
//     cout <<"sol[3] "<< Un1[3] << endl;
//     cout <<"sol[4] "<< Un1[4] <<endl;
//     cout <<"sol[5] "<< Un1[5] <<endl;
//     
    
    rotvecPaxes[0] = rotvecPaxes[0] - Un0[3];
    rotvecPaxes[1] = rotvecPaxes[1] - Un0[4];
    rotvecPaxes[2] = rotvecPaxes[2] - Un0[5];
    
//     cout <<"rotvecPaxes[0] "<< rotvecPaxes[0] << endl;
//     cout <<"rotvecPaxes[1] "<< rotvecPaxes[1] <<endl;
//     cout <<"rotvecPaxes[2] "<< rotvecPaxes[2] <<endl;
    
//     exit(0);
    
    stdAxeswrtPaxes(Paxes, stdaxeswrtpaxes); // obtained standard axes WRT Paxes      
    shifttostdaxes(rotvecPaxes, stdaxeswrtpaxes, rotvecstd); // projected rotvecPaxes to standard axes.    
    alpha1 = rotvecstd[0];
    beta1  = rotvecstd[1];
    gamma1 = rotvecstd[2];
      
    cout <<"sol[0] "<< rotvecstd[0] << endl;
    cout <<"sol[1] "<< rotvecstd[1] <<endl;
    cout <<"sol[2] "<< rotvecstd[2] <<endl;    
    
    for(j=0;j<N_GridDOFs;j++)  // rotation of x, y and z coordinates of grid
    {           
      rotgridx[j] = gridpos[j]*cos(gamma1)*cos(beta1) - gridpos[j+N_GridDOFs]*sin(gamma1)*cos(beta1) + gridpos[j+2*N_GridDOFs]*sin(beta1);
      rotgridy[j] = gridpos[j]*(sin(gamma1)*cos(alpha1) + sin(alpha1)*sin(beta1)*cos(gamma1)) + gridpos[j+N_GridDOFs]*(cos(alpha1)*cos(gamma1) - sin(gamma1)*sin(beta1)*sin(alpha1)) - gridpos[j+2*N_GridDOFs]*sin(alpha1)*cos(beta1);
      rotgridz[j] = gridpos[j]*(sin(alpha1)*sin(gamma1) - cos(alpha1)*sin(beta1)*cos(gamma1)) + gridpos[j+N_GridDOFs]*(cos(gamma1)*sin(alpha1) + cos(alpha1)*sin(beta1)*sin(gamma1)) + gridpos[j+2*N_GridDOFs]*cos(beta1)*cos(alpha1);         
    }

//     cout << "Pos " << Ddot(3*N_GridDOFs, gridpos, gridpos)<<endl;
    
    RotatePaxes(Paxes, rotvecstd);
    
    memcpy(gridpos, rotgridx, N_GridDOFs*sizeof(double));
    memcpy(gridpos+N_GridDOFs, rotgridy, N_GridDOFs*sizeof(double));
    memcpy(gridpos+2*N_GridDOFs, rotgridz, N_GridDOFs*sizeof(double));

    dispx = Un1[0];
    dispy = Un1[1];
    dispz = Un1[2];
    
//     cout << "dispx "<< Un1[0]<<endl;
//     cout << "dispy "<< Un1[1]<<endl;
//     cout << "dispz "<< Un1[2]<<endl;
    
    vecadd(N_GridDOFs, dispx+CGx_Body, gridpos);    // translate gridpos
    vecadd(N_GridDOFs, dispy+CGy_Body, gridpos+N_GridDOFs);    
    vecadd(N_GridDOFs, dispz+CGz_Body, gridpos+2*N_GridDOFs);  

    memcpy(Un0, Un1, 12*sizeof(double));
    GridPos->DataToGrid();
    
// cout << "Pos " << Ddot(3*N_GridDOFs, gridpos, gridpos)<<endl;
//    exit(0);    
    
  if(m==1 || m % TDatabase::TimeDB->STEPS_PER_IMAGE == 0) 
     if(TDatabase::ParamDB->WRITE_VTK)
      {
       os.seekp(std::ios::beg);
        if(img<10) os <<  "VTK/"<<VtkBaseName<<".0000"<<img<<".vtk" << ends;
         else if(img<100) os <<  "VTK/"<<VtkBaseName<<".000"<<img<<".vtk" << ends;
          else if(img<1000) os <<  "VTK/"<<VtkBaseName<<".00"<<img<<".vtk" << ends;
           else if(img<10000) os <<  "VTK/"<<VtkBaseName<<".0"<<img<<".vtk" << ends;
            else  os <<  "VTK/"<<VtkBaseName<<"."<<img<<".vtk" << ends;
       Output->WriteVtk(os.str().c_str());
       img++;
      }    
     t+=deltat;
     
    vecadd(N_GridDOFs, -dispx-CGx_Body, gridpos);    // again shift to CM for next rotation
    vecadd(N_GridDOFs, -dispy-CGy_Body, gridpos+N_GridDOFs);    
    vecadd(N_GridDOFs, -dispz-CGz_Body, gridpos+2*N_GridDOFs);  
  }  
  delete[] gridpos_CM, gridpos, rotgridx, rotgridy, rotgridz;
}
