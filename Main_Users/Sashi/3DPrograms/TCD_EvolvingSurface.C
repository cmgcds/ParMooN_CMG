// =======================================================================
//
// Purpose:     Main program for computing scalars on evolving surface
//              No Multigrid !
// Author:      Sashikumaar Ganesan    26/06/09
//
// =======================================================================

#include <omp.h>

#include <Domain.h>
#include <Database.h>
#include <FEDatabase3D.h>
#include <FEDatabase2D.h>
#include <FESpace3D.h>
#include <FESpace2D.h>
#include <SquareStructure3D.h>
#include <Structure3D.h>
#include <SquareStructure2D.h>
#include <Structure2D.h>
#include <AuxParam3D.h>
#include <AuxParam2D.h>
#include <Aux2D3D.h>
#include <QuadAffin.h>
#include <Solver.h>
#include <Assemble3D.h>
#include <Assemble2D.h>
#include <Output3D.h>
#include <Output2D.h>
#include <DiscreteForm3D.h>
#include <DiscreteForm2D.h>
#include <LinAlg.h>
#include <TNSE3D_ParamRout.h>
#include <BoundFace.h>
#include <Collection.h>

#include <string.h>
#include <sstream>
#include <MooNMD_Io.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>

#include <DirectSolver.h>
#include <ParDirectSolver.h>

#include <MainUtilities.h>
#include <TimeUtilities.h>

#define AMG 0
#define GMG 1

#include <FreeSurface3D.h>
#include <IsoBoundFace.h>
#include <HexaAffin.h>
#include <HexaIsoparametric.h>
#include <TetraAffin.h>
#include <TetraIsoparametric.h>

#include <tetgen.h>
#include <GridCell.h>
#include <MacroCell.h>
#include <BdPlane.h>
#include <BdSphere.h>
#include <IsoBoundFace.h>
#include <InterfaceJoint3D.h>
#include <IsoInterfaceJoint3D.h>

// =======================================================================
// include current example
// =======================================================================
// #include "Examples/TNSE_3D/Drop3D.h"
 #include "Examples/TCD_3D/EvolvingEllipsoid.h"
//  #include "Examples/TCD_3D/HeatEqOnSphere.h"
//  #include "Examples/TCD_3D/ExpandingSphere.h"
//  #include "Examples/TCD_3D/ExpandingSphere_ConstMass.h"
//  #include "Examples/TCD_3D/ExpandingSphere_ConstMassNoDiff.h"

void  AssembleSurf2D(int n_fespaces, TFESpace3D **fespaces, TFEFunction3D **fefunctions, 
                     int N_FESpaces_low, TFESpace2D **fespaces_low, int N_SquareMatrices,
                     TSquareMatrix2D **sqmatrices_low, int N_Rhs, double **rhs, 
                     TFESpace2D **ferhs_low, int *Cell_array, int *Joint_array, double time, double dt, double scale)
{
  TCollection *Coll, *Coll_low;
  TBaseCell *Me, *Me_low;
  FE3D Velo_FEId, FEId;
  FE2D FEId_low;
  TFE2D *Element;

  BaseFunct3D LocBF[N_BaseFuncts3D];
  BaseFunct3D *BaseFuncts;
  BF3DRefElements RefElement;
  RefTrans3D RefTrans;
  TRefTrans3D *F_K;
  QuadFormula2D LineQuadFormula;
  QuadFormula3D QF3;
  TQuadFormula2D *qf;
  QuadFormula2D QuadFormula;
  TFEDesc3D *Velo_FeDesc, *FeDesc;
  TFEDesc2D *FeDesc_low;

  int i, j, k, l, m, n, N_Cells_low, N, N_LocalUsedElements, local_i, local_j, ORDER;
  int Velo_N_BaseFunct, N_BaseFunct_low,  N_Points, N_JointDOF, begin, end, *N_BaseFuncts;
  int LocN_BF[N_BaseFuncts3D],  *KCol, *RowPtr, *JointDOF;
  int *BeginIndex_low, *GlobalNumbers_low, *Velo_DOF, *DOF_LOW, TestDOF, AnsatzDOF, IJoint ;
  int *Velo_BeginIndex, *Velo_GlobalNumbers, Velo_N_JointDOF;
  int maxlen, index1, index2, ORDER_LOW, N_BaseFunct;
  const int *TmpFV, *TmpFV_aux, *TmpLen, *TmpLen_aux;


  double x0, y0, x1, y1, t0, t1, n0, n1, normn;
  double AbsDetjk[MaxN_QuadPoints_3D], Mult;
  double *weights, *xi, *eta;
  double **uref, **uxiref, **uetaref, **uzetaref;
  double **Velo_uref, **Velo_uxiref, **Velo_uetaref, **Velo_uzetaref;
  double *Weights, *p1, *p2, *zeta, *ValuesA, *ValuesM, *ValuesS;
  double X[MaxN_QuadPoints_3D], Y[MaxN_QuadPoints_3D], Z[MaxN_QuadPoints_3D];
  double uorig[MaxN_BaseFunctions3D], uxorig[MaxN_BaseFunctions3D];
  double uyorig[MaxN_BaseFunctions3D], uzorig[MaxN_BaseFunctions3D];
  double val, theta, ngrad_ansatz, ngrad_test, TangDivU, UTgrad;
  double  X_B[100], Y_B[100],  d1, d2, e1, e2;
  double LocMatrixA[MaxN_BaseFunctions3D*MaxN_BaseFunctions3D];
  double LocMatrixM[MaxN_BaseFunctions3D*MaxN_BaseFunctions3D];
  double LocMatrixS[MaxN_BaseFunctions3D*MaxN_BaseFunctions3D];
  double LocRhs[MaxN_BaseFunctions3D];
  double ansatz000, ansatz100, ansatz010, ansatz001;
  double test000, test100, test010,  test001, *u1, *u2, *u3, u1x, u2x, u3x, u1y, u2y, u3y;
  double u1z, u2z, u3z, U1, U2, U3;
  double a1, a2, a3, b1, b2, b3;
  double n2, n3, len,  X1, Y1,  Z1;
  double  D, d3, ngrad, surfacearea=0.;
  double  e3, ngrad2, error=0.;
  double c0, Re = TDatabase::ParamDB->RE_NR, r2, r;
  double Pr = TDatabase::ParamDB->Pr_NR, rhsval;
  double amp=TDatabase::ParamDB->P15;
  double t = TDatabase::TimeDB->CURRENTTIME, v;
  double a=1.+amp*sin(t);
  double da=amp*cos(t);

  c0 =  Re;
//   c0 =  1.0;
  BaseFuncts = TFEDatabase3D::GetBaseFunct3D_IDFromFE3D();
  N_BaseFuncts = TFEDatabase3D::GetN_BaseFunctFromFE3D();

  Coll = fespaces[0]->GetCollection(); // all spaces use same Coll

  Velo_BeginIndex = fespaces[1]->GetBeginIndex();
  Velo_GlobalNumbers = fespaces[1]->GetGlobalNumbers();
  u1 = fefunctions[0]->GetValues();
  u2 = fefunctions[1]->GetValues();
  u3 = fefunctions[2]->GetValues();

  Coll_low = fespaces_low[0]->GetCollection(); // all low spaces use same Coll
  N_Cells_low = Coll_low->GetN_Cells();
  BeginIndex_low =  fespaces_low[0]->GetBeginIndex();
  GlobalNumbers_low =  fespaces_low[0]->GetGlobalNumbers();

  RowPtr = sqmatrices_low[0]->GetRowPtr();
  KCol = sqmatrices_low[0]->GetKCol();

  ValuesA = sqmatrices_low[0]->GetEntries();
  ValuesM = sqmatrices_low[1]->GetEntries();
  ValuesS = sqmatrices_low[2]->GetEntries();

  N_LocalUsedElements = n_fespaces;

// ########################################################################
// loop over all low space cells
// ########################################################################
  for(i=0;i<N_Cells_low;i++)
  {
    N = Cell_array[i];
    IJoint = Joint_array[i];

    Me = Coll->GetCell(N);
    FEId = fespaces[0]->GetFE3D(N, Me);  // scalar space in the entire domain
    FeDesc = TFEDatabase3D::GetFEDesc3DFromFE3D(FEId);
    N_BaseFunct = FeDesc->GetN_DOF();
    N_JointDOF = FeDesc->GetN_JointDOF();
    JointDOF = FeDesc->GetJointDOF(IJoint);
    RefElement = TFEDatabase3D::GetRefElementFromFE3D(FEId);

    DOF_LOW = GlobalNumbers_low + BeginIndex_low[i];
    Me_low = Coll_low->GetCell(i);
    FEId_low = fespaces_low[0]->GetFE2D(i, Me_low);
    Element = TFEDatabase2D::GetFE2D(FEId_low);
    N_BaseFunct_low = Element->GetN_DOF();

    if(N_JointDOF != N_BaseFunct_low )
     {
      cout<< N_JointDOF << " N_JointDOF != N_BaseFunct_low  " << N_BaseFunct_low <<endl;
      exit(0);
    }

    Velo_FEId = fespaces[1]->GetFE3D(N, Me);  // velocity space in the entire domain
    Velo_FeDesc = TFEDatabase3D::GetFEDesc3DFromFE3D(Velo_FEId);
    Velo_N_BaseFunct = Velo_FeDesc->GetN_DOF();
    Velo_N_JointDOF = FeDesc->GetN_JointDOF();
    Velo_DOF = Velo_GlobalNumbers + Velo_BeginIndex[N];

    memset(LocMatrixA, 0, N_BaseFunct_low*N_BaseFunct_low*SizeOfDouble);
    memset(LocMatrixM, 0, N_BaseFunct_low*N_BaseFunct_low*SizeOfDouble);
    memset(LocMatrixS, 0, N_BaseFunct_low*N_BaseFunct_low*SizeOfDouble);
    memset(LocRhs, 0, N_BaseFunct_low*SizeOfDouble);

    l = TFEDatabase3D::GetPolynomialDegreeFromFE3D(FEId);
    ORDER = TFEDatabase3D::GetAccuracyFromFE3D(FEId);


    switch(RefElement)
    {
      case BFUnitHexahedron:

        QuadFormula = TFEDatabase3D::GetQFQuadFromDegree(3*l);
        qf = TFEDatabase3D::GetQuadFormula2D(QuadFormula);
        qf->GetFormulaData(N_Points, Weights, p1, p2);

        RefTrans = HexaIsoparametric;
        F_K = TFEDatabase3D::GetRefTrans3D(RefTrans);

        QF3 = TFEDatabase3D::GetQFHexaFromDegree(3*l);
        ((THexaIsoparametric *)F_K)->SetQuadFormula(QF3);
        ((THexaIsoparametric *)F_K)->SetApproximationOrder(ORDER);
        ((THexaIsoparametric *)F_K)->SetCell(Me);

      break;

      case BFUnitTetrahedron:

       QuadFormula = TFEDatabase3D::GetQFTriaFromDegree(3*l);
       qf = TFEDatabase3D::GetQuadFormula2D(QuadFormula);
       qf->GetFormulaData(N_Points, Weights, p1, p2);

       if(TDatabase::ParamDB->USE_ISOPARAMETRIC)
       {
        RefTrans = TetraIsoparametric;
        F_K = TFEDatabase3D::GetRefTrans3D(RefTrans);
        QF3 = TFEDatabase3D::GetQFTetraFromDegree(3*l);
        ((TTetraIsoparametric *)F_K)->SetQuadFormula(QF3);
        ((TTetraIsoparametric *)F_K)->SetApproximationOrder(ORDER);
        ((TTetraIsoparametric *)F_K)->SetCell(Me);
        ((TTetraIsoparametric *)F_K)->GetOrigBoundFromRef(IJoint, N_Points, p1, p2,
                                                          X,  Y,  Z);
        }
        else
        {
        RefTrans = TetraAffin;
        F_K = TFEDatabase3D::GetRefTrans3D(RefTrans);
        QF3 = TFEDatabase3D::GetQFTetraFromDegree(3*l);
        ((TTetraAffin *)F_K)->SetCell(Me);
        ((TTetraAffin *)F_K)->GetOrigBoundFromRef(IJoint, N_Points, p1, p2, X,  Y,  Z);
        }

      break;
    } // endswitch

// cout << " QuadFormula " << QuadFormula <<endl;

    TFEDatabase3D::GetBaseFunct3DFromFE3D(Velo_FEId)->MakeRefElementData(QuadFormula);

    Velo_uref = TFEDatabase3D::GetJointValues3D(BaseFuncts[Velo_FEId], QuadFormula, IJoint);
    Velo_uxiref = TFEDatabase3D::GetJointDerivatives3D(BaseFuncts[Velo_FEId],
              QuadFormula, IJoint, D100);
    Velo_uetaref = TFEDatabase3D::GetJointDerivatives3D(BaseFuncts[Velo_FEId],
              QuadFormula, IJoint, D010);
    Velo_uzetaref = TFEDatabase3D::GetJointDerivatives3D(BaseFuncts[Velo_FEId],
              QuadFormula, IJoint, D001);

//     for(k=0;k<1;k++)
//      {
//      for(l=0;l<Velo_N_BaseFunct;l++)
//       cout<< "  " << Velo_uref[k][l] << " " << Velo_uxiref[k][l] << " " << Velo_uetaref[k][l] << " " <<  Velo_uzetaref[k][l]<<endl;
//       cout<<endl;
//      }
// // exit(0);

    TFEDatabase3D::GetBaseFunct3DFromFE3D(FEId)->MakeRefElementData(QuadFormula);

    uref = TFEDatabase3D::GetJointValues3D(BaseFuncts[FEId],
              QuadFormula, IJoint);
    uxiref = TFEDatabase3D::GetJointDerivatives3D(BaseFuncts[FEId],
              QuadFormula, IJoint, D100);
    uetaref = TFEDatabase3D::GetJointDerivatives3D(BaseFuncts[FEId],
              QuadFormula, IJoint, D010);
    uzetaref = TFEDatabase3D::GetJointDerivatives3D(BaseFuncts[FEId],
              QuadFormula, IJoint, D001);

//     for(k=0;k<1;k++)
//      {
//      for(l=0;l<N_BaseFunct;l++)
//       cout<< "  " << uref[k][l] << " " << uxiref[k][l] << " " << uetaref[k][l] << " " <<  uzetaref[k][l]<<endl;
//       cout<<endl;
//      }
// exit(0);

    for(k=0;k<N_Points;k++)
    {
      switch(RefElement)
      {
        case BFUnitHexahedron:
          ((THexaIsoparametric *)F_K)->GetTangentVectors(
                IJoint, p1[k], p2[k], a1, a2, a3, b1, b2, b3);
        break;

       case BFUnitTetrahedron:
       if(TDatabase::ParamDB->USE_ISOPARAMETRIC)
          ((TTetraIsoparametric *)F_K)->GetTangentVectors(
                IJoint, p1[k], p2[k], a1, a2, a3, b1, b2, b3);
       else
          ((TTetraAffin *)F_K)->GetTangentVectors(
                IJoint, p1[k], p2[k], a1, a2, a3, b1, b2, b3);
        break;
      } // endswitch


      n1 = a2*b3 - a3*b2;
      n2 = a3*b1 - a1*b3;
      n3 = a1*b2 - a2*b1;

//       OutPut("xi: " << p1[k] << " eta: " << p2[k] << endl);
//       OutPut("point: " << k << " t1: " << a1 << " " << a2 << " " << a3 << endl);
//       OutPut("point: " << k << " t2: " << b1 << " " << b2 << " " << b3 << endl);

      len = sqrt(n1*n1 + n2*n2 + n3*n3);
      n1 /= len;
      n2 /= len;
      n3 /= len;

//     OutPut("point: " << k << " normal: " << n1 << " " << n2 << " " << n3 << endl);

//       get velocity and its derivative at this quad point
      switch(RefElement)
       {
        case BFUnitHexahedron:
          ((THexaIsoparametric *)F_K)->GetOrigValues(
                IJoint, p1[k], p2[k], Velo_N_BaseFunct,
                Velo_uref[k], Velo_uxiref[k], Velo_uetaref[k], Velo_uzetaref[k],
                uorig, uxorig, uyorig, uzorig);
        break;


        case BFUnitTetrahedron:
        if(TDatabase::ParamDB->USE_ISOPARAMETRIC)
          ((TTetraIsoparametric *)F_K)->GetOrigValues(
                IJoint, p1[k], p2[k], Velo_N_BaseFunct,
                Velo_uref[k], Velo_uxiref[k], Velo_uetaref[k], Velo_uzetaref[k],
                uorig, uxorig, uyorig, uzorig);
        else
          ((TTetraAffin *)F_K)->GetOrigValues(
                IJoint, p1[k], p2[k], Velo_N_BaseFunct,
                Velo_uref[k], Velo_uxiref[k], Velo_uetaref[k], Velo_uzetaref[k],
                uorig, uxorig, uyorig, uzorig);
        break;
       } // endswitch

          u1x=0.;  u1y=0.; u1z=0.;
          u2x=0.;  u2y=0.; u2z=0.;
          u3x=0.;  u3y=0.; u3z=0.;
          U1=0.; U2=0.; U3=0.;

          for(l=0;l<Velo_N_BaseFunct;l++)
            {
             m = Velo_DOF[l];
             U1 +=u1[m]*uorig[l];
             u1x += u1[m]*uxorig[l];
             u1y += u1[m]*uyorig[l];
             u1z += u1[m]*uzorig[l];

             U2 +=u2[m]*uorig[l];
             u2x += u2[m]*uxorig[l];
             u2y += u2[m]*uyorig[l];
             u2z += u2[m]*uzorig[l];

             U3 +=u3[m]*uorig[l];
             u3x += u3[m]*uxorig[l];
             u3y += u3[m]*uyorig[l];
             u3z += u3[m]*uzorig[l];

            }
//          cout <<" U.n:  " << n1*U1 + n2*U2 + n3*U3 <<" exact:  " <<  v*sqrt(D)  <<endl;

//          v=n1*U1 + n2*U2 + n3*U3;
//          U1 -=v*n1;
//          U2 -=v*n2;
//          U3 -=v*n3;

//          D = X[k]*X[k]/(a*a) + Y[k]*Y[k] + Z[k]*Z[k];
//          v=da*X[k]*X[k]/(2.*D*a*a);
//       cout <<" U1:  " << U1 <<" n1 exact:  " << da/(2.*a)*X[k] - v*X[k]/a  <<endl;
//       cout <<" U2:  " << U2 <<" n2 exact:  " <<  v*Y[k] <<endl;
//       cout <<" U3:  " << U3 <<" n3 exact:  " <<  v*Z[k] <<endl;

// check on sphere
//       v = sqrt(X[k]*X[k] + Y[k]*Y[k] + Z[k]*Z[k]);

//       cout <<" n1x:  " << n1 <<" n1 exact:  " <<  X[k]/(a*sqrt(D))  <<endl;
//       cout <<" n2x:  " << n2 <<" n2 exact:  " <<  Y[k]/sqrt(D)   <<endl;
//       cout <<" n3x:  " << n3 <<" n3 exact:  " <<  Z[k]/sqrt(D)  <<endl;

         TangDivU =   u1x - (u1x*n1 + u1y*n2 + u1z*n3)*n1 +
                      u2y - (u2x*n1 + u2y*n2 + u2z*n3)*n2 +
                      u3z - (u3x*n1 + u3y*n2 + u3z*n3)*n3;

// //          TangDivU =  (amp*cos(t)) /( 1. + amp*sin(t) ); // expanding sphere

//           v = X[k]*X[k]/(a*a);
//          TangDivU =  ( da/(2*a) )*( (Y[k]*Y[k] + Z[k]*Z[k])/(v+Y[k]*Y[k] +Z[k]*Z[k]) ) ;  // evolving ellipsoid
//           cout <<" U1:  " << U1 << " TangDivU :  " <<TangDivU<<  endl;
// Expanding Sphere

//           cout <<" U1:  " << U1 << " TangDivU :  " <<TangDivU<< endl;
//           cout <<" U1:  " << U1 << " da/2a x:  " <<( (amp*cos(t)) /(2.*( 1. + amp*sin(t) ) ))* X[k] << endl;
//           cout <<" u1x:  " << u1x <<" u1y:  " << u1y <<" u1z:  " << u1z <<endl;
//           cout <<" u2x:  " << u2x <<" u2y:  " << u2y<<" u2z:  " << u2z  <<endl;
//           cout <<" u3x:  " << u3x <<" u3y:  " << u3y<<" u3z:  " << u3z <<endl;
//           cout <<" TangDivU:  " << TangDivU << " da/a:  " << (amp*cos(t)) /( 1. + amp*sin(t) )  <<endl;

// evolving ellipsoid

//           cout <<" TangDivU:  " << TangDivU << " da/a:  " << ( da/(2*a) )*( (Y[k]*Y[k] + Z[k]*Z[k])/(v+Y[k]*Y[k] +Z[k]*Z[k]) )  <<endl;

// exit(0);


      switch(RefElement)
       {
        case BFUnitHexahedron:
          ((THexaIsoparametric *)F_K)->GetOrigValues(
                IJoint, p1[k], p2[k], N_BaseFunct,
                uref[k], uxiref[k], uetaref[k], uzetaref[k],
                uorig, uxorig, uyorig, uzorig);
        break;

        case BFUnitTetrahedron:
        if(TDatabase::ParamDB->USE_ISOPARAMETRIC)
          ((TTetraIsoparametric *)F_K)->GetOrigValues(
                IJoint, p1[k], p2[k], N_BaseFunct,
                uref[k], uxiref[k], uetaref[k], uzetaref[k],
                uorig, uxorig, uyorig, uzorig);
        else
          ((TTetraAffin *)F_K)->GetOrigValues(
                IJoint, p1[k], p2[k], N_BaseFunct,
                uref[k], uxiref[k], uetaref[k], uzetaref[k],
                uorig, uxorig, uyorig, uzorig);
        break;
       } // endswitch


          Mult = Weights[k]* len;
//           surfacearea +=Mult;
          for(l=0;l<N_JointDOF;l++)
           {
              local_j   = JointDOF[l];

              test000 = uorig[local_j];
              test100 = uxorig[local_j];
              test010 = uyorig[local_j];
              test001 = uzorig[local_j];

              ngrad = test100*n1 + test010*n2 +test001*n3;

              d1 = test100 - ngrad*n1;
              d2 = test010 - ngrad*n2;
              d3 = test001 - ngrad*n3;

// // testing


// //        for torus example
//              r =  sqrt(X[k]* X[k] + (Y[k]-1.)*(Y[k]-1.) + Z[k]*Z[k]);
//             if(r<0.25)
//                c1 = 100.;
//             else
//                c1 = 0.;
// 

             GetRhs(time, X[k], Y[k],  Z[k], rhsval);
             LocRhs[l] += rhsval*Mult*test000;

             for(m=0;m<N_JointDOF;m++)
              {
               local_i   = JointDOF[m];

               ansatz000 = uorig[local_i];
               ansatz100 = uxorig[local_i];
               ansatz010 = uyorig[local_i];
               ansatz001 = uzorig[local_i];

               ngrad2 = ansatz100*n1 +ansatz010*n2 + ansatz001*n3;

               e1 = ansatz100 - ngrad2*n1;
               e2 = ansatz010 - ngrad2*n2;
               e3 = ansatz001 - ngrad2*n3;

               val = c0*(d1*e1 + d2*e2 + d3*e3);     // diffusive term
//                val += TangDivU*test000*ansatz000;   // tangential divergence term
               val *= Mult;
               LocMatrixA[l*N_JointDOF+m] += val;

               val = TangDivU*test000*ansatz000;   // tangential divergence term
               val *= Mult;
               LocMatrixS[l*N_JointDOF+m] += val;

//                cout << " val " << val << endl;
               val  = test000*ansatz000;
               val *= Mult;
               LocMatrixM[l*N_JointDOF+m] += val;
             } // for(m=0;m<N_BaseFun
          } //  for(l=0;l<N_Joint
   } //  for(k=0;k<N_Points

//   add to global matrices
    for(l=0;l<N_JointDOF;l++)
     {
      TestDOF = DOF_LOW[l];
      begin = RowPtr[TestDOF];
      end = RowPtr[TestDOF+1];

      rhs[0][TestDOF]  +=LocRhs[l];
//           cout << " rhs[0][TestDOF] " <<rhs[0][TestDOF] << " LocRhs[l] " << LocRhs[l] <<endl;
      for(n=begin;n<end;n++)
       {
       for(m=0;m<N_JointDOF;m++)
        {
         if(KCol[n] == DOF_LOW[m])
          {
           ValuesA[n] +=LocMatrixA[l*N_JointDOF+m];
           ValuesM[n] +=LocMatrixM[l*N_JointDOF+m];
           ValuesS[n] +=LocMatrixS[l*N_JointDOF+m];
//           cout << " ValuesA[n] " <<ValuesA[n] <<endl;
           break;
          }
        } // for(m=0;m<N_JointDOF
      } // for(n=begin;n<end;n++)
     } // for(l=0;l<N_JointDOF

  } //  for(i=0;i<N_Cells_low
// cout<< " surfacearea Assemble2D " << surfacearea <<endl;
}

// interpolation for the  fefunction only on the surface
// quadrature points are projected to original surface.
void SurfInterpolate(TFEFunction3D *fefunction, TFEFunction2D *fefunct_low, 
                     int *Cell_array, int *Joint_No, DoubleFunct3D *Exact)
{
  int i, j,k,l, N_Cells_low, N;
  TBaseCell *cell;
  TCollection *Coll, *Coll_low;;
  FE3D FEId;
  TFE3D *Element;
  BaseFunct3D BF;
  TNodalFunctional3D *nf;
  int N_Cells, IJoint, *JointDOF;
  int N_DOFs, N_JointDOFs;
  int *BeginIndex, *GlobalNumbers;
  int N_, N_Points;
  double s, *xi, *eta, *zeta;
  double Val[MaxN_BaseFunctions3D];
  double OutVal[MaxN_BaseFunctions3D];
  int *DOF, Index;
  RefTrans3D F_K;
  TRefTrans3D *rt;
  double X[MaxN_PointsForNodal3D], Y[MaxN_PointsForNodal3D];
  double Z[MaxN_PointsForNodal3D];
  double AbsDetjk[MaxN_PointsForNodal3D];
  double PointValues[MaxN_PointsForNodal3D];
  double FunctionalValues[MaxN_PointsForNodal3D];
  double FctVal[5], *Values,  v;
  int PolynomialDegree, ApproxOrder;
  QuadFormula3D QuadFormula;
  boolean IsIsoparametric;
  TJoint *joint;
  JointType jointtype;
  BoundTypes bdtype;
  int N_Faces;
  BF3DRefElements RefElement;
  RefTrans3D RefTrans, *RefTransArray;
  TFESpace3D *FESpace3D;
  TFESpace2D *fespace_low;
  TFEDesc3D *FEDesc;
  // begin code

  fespace_low=fefunct_low->GetFESpace2D();
  Coll_low = fespace_low->GetCollection();
  N_Cells_low = Coll_low->GetN_Cells();

  FESpace3D=fefunction->GetFESpace3D();
  Values = fefunction->GetValues();
  Coll = FESpace3D->GetCollection();
  N_Cells = Coll->GetN_Cells();
  BeginIndex = FESpace3D->GetBeginIndex();
  GlobalNumbers = FESpace3D->GetGlobalNumbers();
  N_DOFs = FESpace3D->GetN_DegreesOfFreedom();

  memset(Values, 0, SizeOfDouble*N_DOFs);
  RefTransArray = TFEDatabase3D::GetRefTrans3D_IDFromFE3D();

  for(i=0;i<N_Cells_low;i++)
  {
    N = Cell_array[i];
    IJoint = Joint_No[i];

    cell = Coll->GetCell(N);
    FEId = FESpace3D->GetFE3D(N, cell);
    FEDesc = TFEDatabase3D::GetFEDesc3DFromFE3D(FEId);
    JointDOF = FEDesc->GetJointDOF(IJoint);
    N_JointDOFs =FEDesc->GetN_JointDOF();

    Element = TFEDatabase3D::GetFE3D(FEId);
    nf = Element->GetNodalFunctional3D();
    nf->GetPointsForFace(IJoint, N_Points, xi, eta, zeta);

    PolynomialDegree = TFEDatabase3D::GetPolynomialDegreeFromFE3D(FEId);
    ApproxOrder = TFEDatabase3D::GetAccuracyFromFE3D(FEId);

    RefElement = Element->GetBaseFunct3D()->GetRefElement();
    switch(RefElement)
    {
      case BFUnitHexahedron:
        QuadFormula = TFEDatabase3D::GetQFHexaFromDegree
                         (3*PolynomialDegree);
        N_Faces = 6;
      break;

      case BFUnitTetrahedron:
        QuadFormula = TFEDatabase3D::GetQFTetraFromDegree
                         (3*PolynomialDegree);
        N_Faces = 4;
      break;
    }

    RefTrans = RefTransArray[FEId];

    IsIsoparametric = FALSE;
    if (TDatabase::ParamDB->USE_ISOPARAMETRIC)
    {
      for(j=0;j<N_Faces;j++)
      {
        joint = cell->GetJoint(j);
        jointtype = joint->GetType();
        if(jointtype == BoundaryFace)
        {
          bdtype = ((TBoundFace *)(joint))->GetBoundComp()->GetType();
          if(bdtype != Plane)
            IsIsoparametric = TRUE;
        }
        if(jointtype == IsoBoundFace)
          IsIsoparametric = TRUE;
      }
    } // endif 
  
    if(IsIsoparametric)
    {
      switch(RefElement)
      {
        case BFUnitHexahedron:
          RefTrans = HexaIsoparametric;
        break;
  
        case BFUnitTetrahedron:
          RefTrans = TetraIsoparametric;
        break;
      }
    } // endif IsIsoparametric
//     cout << "IsIsoparametric: " << IsIsoparametric << endl;

    switch(RefTrans)
    {
      case HexaAffin:
        rt = TFEDatabase3D::GetRefTrans3D(HexaAffin);
        ((THexaAffin *)rt)->SetCell(cell);
        F_K = HexaAffin;
        break;
/*      case HexaTrilinear:
        rt = TFEDatabase3D::GetRefTrans3D(HexaTrilinear);
        ((THexaTrilinear *)rt)->SetCell(cell);
        F_K = HexaTrilinear;
        break;*/
      case HexaIsoparametric:
        rt = TFEDatabase3D::GetRefTrans3D(HexaIsoparametric);
        ((THexaIsoparametric *)rt)->SetApproximationOrder(ApproxOrder);
        ((THexaIsoparametric *)rt)->SetQuadFormula(QuadFormula);
        ((THexaIsoparametric *)rt)->SetCell(cell);
        F_K = HexaIsoparametric;
        break;
      case TetraAffin:
        rt = TFEDatabase3D::GetRefTrans3D(TetraAffin);
        ((TTetraAffin *)rt)->SetCell(cell);
        F_K = TetraAffin;
        break;
      case TetraIsoparametric:
        rt = TFEDatabase3D::GetRefTrans3D(TetraIsoparametric);
        ((TTetraIsoparametric *)rt)->SetApproximationOrder(ApproxOrder);
        ((TTetraIsoparametric *)rt)->SetQuadFormula(QuadFormula);
        ((TTetraIsoparametric *)rt)->SetCell(cell);
        F_K = TetraIsoparametric;
        break;
    }
    TFEDatabase3D::GetOrigFromRef(F_K, N_Points, xi, eta, zeta,
                                X, Y, Z, AbsDetjk);

    // cout << "----------------" << endl;
    for(j=0;j<N_Points;j++)
     {
      Exact(X[j], Y[j], Z[j], FctVal);
      PointValues[j] = FctVal[0];
     }

    nf->GetFaceFunctionals(PointValues,FunctionalValues);

    DOF = GlobalNumbers+BeginIndex[N];

    for(j=0;j<N_JointDOFs;j++)
      Values[ DOF[ JointDOF[j] ] ] = FunctionalValues[j];
  }

}


// error in the fefunction only on the surface
// quadrature points are projected to original surface.
void GetSurfErrors(TFEFunction3D *fefunction, TFEFunction2D *fefunct_low, double *errors,
                   int *Cell_array, int *Joint_array, DoubleFunct3D *Exact)
{
  int i, j, k, l, m, n, N, N_Cells_low, N_LocalUsedElements, local_i, local_dof, ORDER;
  int N_BaseFunct,  N_Points, N_JointDOF, begin, end, *N_BaseFuncts, IJoint;
  TCollection *Coll, *Coll_low;
  TBaseCell *Me, *Me_low;
  FE3D FEId;
  TFE3D *ele;
  FE2D FEId_low;
  TFE2D *Element;
  TFESpace3D *fespace;
  TFESpace2D *fespace_low;
  BaseFunct3D LocBF[N_BaseFuncts3D];
  BaseFunct3D *BaseFuncts;
  int *DOF, *JointDOF, *BeginIndex, *GlobalNumbers, N_BaseFunct_low;
  TFEDesc3D *FeDesc;
  BF3DRefElements RefElement;
  RefTrans3D RefTrans;
  TRefTrans3D *F_K;
  QuadFormula3D QF3;
  TQuadFormula2D *qf;
  QuadFormula2D QuadFormula;
  double *u, **uref, **uxiref, **uetaref, **uzetaref;
  double *Weights, *p1, *p2, *zeta, LocL2U, LocH1U;
  double *weights, *xi, *eta;
  double  X_B[100], Y_B[100], Values[6], v, U, ux, uy, uz, X_P, Y_P, Z_P;
  double X[MaxN_QuadPoints_3D], Y[MaxN_QuadPoints_3D], Z[MaxN_QuadPoints_3D];
  double uorig[MaxN_BaseFunctions3D], uxorig[MaxN_BaseFunctions3D];
  double uyorig[MaxN_BaseFunctions3D], uzorig[MaxN_BaseFunctions3D];
  double a1, a2, a3, b1, b2, b3, ulx, uly, ulz;
  double n1, n2, n3, AbsDetjk, Mult, SurfArea=0.;
  double test000, test100, test010,  test001, d1, d2, d3, ngrad, nugrad;

  BaseFuncts = TFEDatabase3D::GetBaseFunct3D_IDFromFE3D();
  N_BaseFuncts = TFEDatabase3D::GetN_BaseFunctFromFE3D();

  fespace=fefunction->GetFESpace3D();
  Coll = fespace->GetCollection();
  BeginIndex = fespace->GetBeginIndex();
  GlobalNumbers = fespace->GetGlobalNumbers();
  u = fefunction->GetValues();

  fespace_low=fefunct_low->GetFESpace2D();
  Coll_low = fespace_low->GetCollection();
  N_Cells_low = Coll_low->GetN_Cells();

  errors[0] =0.; // L2-norm
  errors[1] =0.; // H1-seminorm
// ########################################################################
// loop over all surf cells
// ########################################################################
  for(i=0;i<N_Cells_low;i++)
  {
    LocL2U = 0.;
    LocH1U = 0.;

    N = Cell_array[i];
    Me = Coll->GetCell(N);
    IJoint = Joint_array[i];

    FEId = fespace->GetFE3D(N, Me);
    ele = TFEDatabase3D::GetFE3D(FEId);

    RefElement = TFEDatabase3D::GetRefElementFromFE3D(FEId);
    FeDesc = TFEDatabase3D::GetFEDesc3DFromFE3D(FEId);
    N_JointDOF =  FeDesc->GetN_JointDOF();
    JointDOF = FeDesc->GetJointDOF(IJoint);
    N_BaseFunct = FeDesc->GetN_DOF();
    DOF = GlobalNumbers + BeginIndex[N];

    Me_low = Coll_low->GetCell(i);
    FEId_low = fespace_low->GetFE2D(i, Me_low);
    Element = TFEDatabase2D::GetFE2D(FEId_low);
    N_BaseFunct_low = Element->GetN_DOF();

    l = TFEDatabase3D::GetPolynomialDegreeFromFE3D(FEId);
    ORDER = TFEDatabase3D::GetAccuracyFromFE3D(FEId);
//     ORDER = 2;
    switch(RefElement)
    {
      case BFUnitHexahedron:

        QuadFormula = TFEDatabase3D::GetQFQuadFromDegree(3*l);
        qf = TFEDatabase3D::GetQuadFormula2D(QuadFormula);
        qf->GetFormulaData(N_Points, Weights, p1, p2);

        RefTrans = HexaIsoparametric;
        F_K = TFEDatabase3D::GetRefTrans3D(RefTrans);

        QF3 = TFEDatabase3D::GetQFHexaFromDegree(3*l);
        ((THexaIsoparametric *)F_K)->SetQuadFormula(QF3);
        ((THexaIsoparametric *)F_K)->SetApproximationOrder(ORDER);
        ((THexaIsoparametric *)F_K)->SetCell(Me);

      break;

      case BFUnitTetrahedron:

       QuadFormula = TFEDatabase3D::GetQFTriaFromDegree(3*l);
       qf = TFEDatabase3D::GetQuadFormula2D(QuadFormula);
       qf->GetFormulaData(N_Points, Weights, p1, p2);

       if(TDatabase::ParamDB->USE_ISOPARAMETRIC)
       {
        RefTrans = TetraIsoparametric;
        F_K = TFEDatabase3D::GetRefTrans3D(RefTrans);
        QF3 = TFEDatabase3D::GetQFTetraFromDegree(3*l);
        ((TTetraIsoparametric *)F_K)->SetQuadFormula(QF3);
        ((TTetraIsoparametric *)F_K)->SetApproximationOrder(ORDER);
        ((TTetraIsoparametric *)F_K)->SetCell(Me);
        ((TTetraIsoparametric *)F_K)->GetOrigBoundFromRef(IJoint, N_Points, p1, p2,
                                                          X,  Y,  Z);
        }
        else
        {
        RefTrans = TetraAffin;
        F_K = TFEDatabase3D::GetRefTrans3D(RefTrans);
        QF3 = TFEDatabase3D::GetQFTetraFromDegree(3*l);
        ((TTetraAffin *)F_K)->SetCell(Me);
        ((TTetraAffin *)F_K)->GetOrigBoundFromRef(IJoint, N_Points, p1, p2, X,  Y,  Z);
        }

      break;
    } // endswitch

// cout << " QuadFormula " << QuadFormula <<endl;

    TFEDatabase3D::GetBaseFunct3DFromFE3D(FEId)->MakeRefElementData(QuadFormula);

    for(k=0;k<N_Points;k++)
    {
      switch(RefElement)
      {
        case BFUnitHexahedron:
          ((THexaIsoparametric *)F_K)->GetTangentVectors(
                IJoint, p1[k], p2[k], a1, a2, a3, b1, b2, b3);
        break;

       case BFUnitTetrahedron:
       if(TDatabase::ParamDB->USE_ISOPARAMETRIC)
          ((TTetraIsoparametric *)F_K)->GetTangentVectors(
                IJoint, p1[k], p2[k], a1, a2, a3, b1, b2, b3);
       else
          ((TTetraAffin *)F_K)->GetTangentVectors(
                IJoint, p1[k], p2[k], a1, a2, a3, b1, b2, b3);
        break;
      } // endswitch

      n1 = a2*b3 - a3*b2;
      n2 = a3*b1 - a1*b3;
      n3 = a1*b2 - a2*b1;
      AbsDetjk = sqrt(n1*n1 + n2*n2 + n3*n3);
      n1 /= AbsDetjk;
      n2 /= AbsDetjk;
      n3 /= AbsDetjk;

      Mult = Weights[k]* AbsDetjk;

      uref = TFEDatabase3D::GetJointValues3D(BaseFuncts[FEId],
                QuadFormula, IJoint);
      uxiref = TFEDatabase3D::GetJointDerivatives3D(BaseFuncts[FEId],
                QuadFormula, IJoint, D100);
      uetaref = TFEDatabase3D::GetJointDerivatives3D(BaseFuncts[FEId],
                QuadFormula, IJoint, D010);
      uzetaref = TFEDatabase3D::GetJointDerivatives3D(BaseFuncts[FEId],
                QuadFormula, IJoint, D001);


      switch(RefElement)
      {
        case BFUnitHexahedron:
          ((THexaIsoparametric *)F_K)->GetOrigValues(
                IJoint, p1[k], p2[k], N_BaseFunct,
                uref[k], uxiref[k], uetaref[k], uzetaref[k],
                uorig, uxorig, uyorig, uzorig);
        break;

        case BFUnitTetrahedron:
        if(TDatabase::ParamDB->USE_ISOPARAMETRIC)
          ((TTetraIsoparametric *)F_K)->GetOrigValues(
                IJoint, p1[k], p2[k], N_BaseFunct,
                uref[k], uxiref[k], uetaref[k], uzetaref[k],
                uorig, uxorig, uyorig, uzorig);
        else
          ((TTetraAffin *)F_K)->GetOrigValues(
                IJoint, p1[k], p2[k], N_BaseFunct,
                uref[k], uxiref[k], uetaref[k], uzetaref[k],
                uorig, uxorig, uyorig, uzorig);
        break;
      } // endswitch

//       get solution and its gradients
     U=0.; ux=0.;  uy=0.;  uz=0.;
     for(l=0;l<N_BaseFunct_low;l++)
      {
       local_dof   = JointDOF[l];

       test000 = uorig[local_dof];
       test100 = uxorig[local_dof];
       test010 = uyorig[local_dof];
       test001 = uzorig[local_dof];
       ngrad = test100*n1 + test010*n2 +test001*n3;

       d1 = test100 - ngrad*n1;
       d2 = test010 - ngrad*n2;
       d3 = test001 - ngrad*n3;

       m = DOF[local_dof];

       U  += u[m]*test000;
       ux += u[m]*d1;
       uy += u[m]*d2;
       uz += u[m]*d3;
      }

//  Projection of all points to the free boundary in normal direction 
//       v =  sqrt(X[k]* X[k] + Y[k]*Y[k] + Z[k]*Z[k] );// /r;
//       X[k] /= v;
//       Y[k] /=  v;
//       Z[k] /=  v;
// projected in the Example file, not here

      Exact(X[k], Y[k], Z[k], Values);

      nugrad = Values[1]*n1 + Values[2]*n2 + Values[3]*n3;
      ulx = Values[1] - nugrad*n1;
      uly = Values[2] - nugrad*n2;
      ulz = Values[3] - nugrad*n3;

      errors[0] +=Mult*(Values[0]-U)*(Values[0]-U);
      errors[1] +=Mult*(ux-ulx)*(ux-ulx);
      errors[1] +=Mult*(uy-uly)*(uy-uly);
      errors[1] +=Mult*(uz-ulz)*(uz-ulz);

      SurfArea +=Mult;
// cout << " " << uz - ulz<< " "<<endl;
    } // for(k=0;k<N_Points
 } //   for(i=0;i<N_Cells_low;i++)
   errors[0] = sqrt(errors[0]);
   errors[1] = sqrt(errors[1]);

//   double r, amp=TDatabase::ParamDB->P15;
//   double t=TDatabase::TimeDB->CURRENTTIME;
//    r = sqrt(1.+amp*sin(t));
//    cout << " surface area : " << SurfArea <<" Exact : " <<  4.*Pi*r*r<<  endl;
}


void GetSurfMass(TFESpace3D *fespace, TFEFunction2D *fefunct_low, int *Cell_array, int *Joint_array, double *Mass)
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

 TCollection *Coll_low;
 TBaseCell *Me_low;
 FE3D Velo_FEId;
 FE2D FEId_low;
 TFE2D *Element_low;
 TFESpace2D *fespace_low;

 BaseFunct3D LocBF[N_BaseFuncts3D];
 BaseFunct3D *BaseFuncts;
 BF3DRefElements RefElement;
 RefTrans3D RefTrans;
 TRefTrans3D *F_K;
 QuadFormula2D LineQuadFormula;
 QuadFormula3D QF3;
 TQuadFormula2D *qf;
 QuadFormula2D QuadFormula;
 TFEDesc3D *Velo_FeDesc, *FeDesc;
 TFEDesc2D *FeDesc_low;

 int i, j, k, l, m, n, N_U, N_Cells, N_Faces, N_LocalDOFs, *JointDOF, dof;
 int N_Inner, N_BoundaryNodes;
//  int *BeginIndex, *GlobalNumbers, *DOF;
 int N_Cells_low, N, N_LocalUsedElements, local_i, local_j, ORDER;
 int N_BaseFunct_low,  N_Points, N_JointDOF, begin, end, *N_BaseFuncts, N_BaseFunct;
 int *BeginIndex_low, *GlobalNumbers_low, *DOF_LOW, TestDOF, AnsatzDOF, IJoint;

 double x, y, z, t1, t2, Mult, t_new;
 double IsoX, IsoY, IsoZ;
 double *RefX, *RefY, *RefZ;
 double *ValuesNewX, *ValuesNewY, *ValuesNewZ;
 double *PreviousPosX, *PreviousPosY, *PreviousPosZ;
 double a, a_old, theta, phi;
 double *Isotheta, *Isophi;
 double  amp=TDatabase::ParamDB->P15, da;
 double *values, *Weights, *p1, *p2, **uref;
 double LocMatrixM[MaxN_BaseFunctions3D*MaxN_BaseFunctions3D];
 double X[MaxN_QuadPoints_3D], Y[MaxN_QuadPoints_3D], Z[MaxN_QuadPoints_3D];
 double a1, a2, a3, b1, b2, b3, n1, n2, n3, len;
 double test000, ansatz000, val;


  BaseFuncts = TFEDatabase3D::GetBaseFunct3D_IDFromFE3D();
  N_BaseFuncts = TFEDatabase3D::GetN_BaseFunctFromFE3D();


  Coll = fespace->GetCollection();
//   BeginIndex = fespace->GetBeginIndex();
//   GlobalNumbers = fespace->GetGlobalNumbers();

  fespace_low=fefunct_low->GetFESpace2D();
  Coll_low = fespace_low->GetCollection();
  N_Cells_low = Coll_low->GetN_Cells();
  BeginIndex_low =  fespace_low->GetBeginIndex();
  GlobalNumbers_low =  fespace_low->GetGlobalNumbers();
  values = fefunct_low->GetValues();

  Mass[0] = 0.;
  Mass[1] = 0.;

// ########################################################################
// loop over all low space cells
// ########################################################################
  for(i=0;i<N_Cells_low;i++)
   {
    N = Cell_array[i];
    IJoint = Joint_array[i];

    cell = Coll->GetCell(N);
    FEId = fespace->GetFE3D(N, cell);  // scalar space in the entire domain
    FeDesc = TFEDatabase3D::GetFEDesc3DFromFE3D(FEId);
    N_BaseFunct = FeDesc->GetN_DOF();
    N_JointDOF = FeDesc->GetN_JointDOF();
    JointDOF = FeDesc->GetJointDOF(IJoint);
    RefElement = TFEDatabase3D::GetRefElementFromFE3D(FEId);

    DOF_LOW = GlobalNumbers_low + BeginIndex_low[i];
    Me_low = Coll_low->GetCell(i);
    FEId_low = fespace_low->GetFE2D(i, Me_low);
    Element_low = TFEDatabase2D::GetFE2D(FEId_low);
    N_BaseFunct_low = Element_low->GetN_DOF();

    if(N_JointDOF != N_BaseFunct_low )
     {
      cout<< N_JointDOF << " N_JointDOF != N_BaseFunct_low  " << N_BaseFunct_low <<endl;
      exit(0);
    }

    memset(LocMatrixM, 0, N_BaseFunct_low*N_BaseFunct_low*SizeOfDouble);

    l = TFEDatabase3D::GetPolynomialDegreeFromFE3D(FEId);
    ORDER = TFEDatabase3D::GetAccuracyFromFE3D(FEId);

    switch(RefElement)
    {
      case BFUnitHexahedron:

        QuadFormula = TFEDatabase3D::GetQFQuadFromDegree(3*l);
        qf = TFEDatabase3D::GetQuadFormula2D(QuadFormula);
        qf->GetFormulaData(N_Points, Weights, p1, p2);

        RefTrans = HexaIsoparametric;
        F_K = TFEDatabase3D::GetRefTrans3D(RefTrans);

        QF3 = TFEDatabase3D::GetQFHexaFromDegree(3*l);
        ((THexaIsoparametric *)F_K)->SetQuadFormula(QF3);
        ((THexaIsoparametric *)F_K)->SetApproximationOrder(ORDER);
        ((THexaIsoparametric *)F_K)->SetCell(cell);

      break;

      case BFUnitTetrahedron:

       QuadFormula = TFEDatabase3D::GetQFTriaFromDegree(3*l);
       qf = TFEDatabase3D::GetQuadFormula2D(QuadFormula);
       qf->GetFormulaData(N_Points, Weights, p1, p2);

       if(TDatabase::ParamDB->USE_ISOPARAMETRIC)
       {
        RefTrans = TetraIsoparametric;
        F_K = TFEDatabase3D::GetRefTrans3D(RefTrans);
        QF3 = TFEDatabase3D::GetQFTetraFromDegree(3*l);
        ((TTetraIsoparametric *)F_K)->SetQuadFormula(QF3);
        ((TTetraIsoparametric *)F_K)->SetApproximationOrder(ORDER);
        ((TTetraIsoparametric *)F_K)->SetCell(cell);
        ((TTetraIsoparametric *)F_K)->GetOrigBoundFromRef(IJoint, N_Points, p1, p2,
                                                          X,  Y,  Z);
        }
       else
        {
        RefTrans = TetraAffin;
        F_K = TFEDatabase3D::GetRefTrans3D(RefTrans);
        QF3 = TFEDatabase3D::GetQFTetraFromDegree(3*l);
        ((TTetraAffin *)F_K)->SetCell(cell);
        ((TTetraAffin *)F_K)->GetOrigBoundFromRef(IJoint, N_Points, p1, p2, X,  Y,  Z);
        }

      break;
    } // endswitch

// cout << " QuadFormula " << QuadFormula <<endl;

    TFEDatabase3D::GetBaseFunct3DFromFE3D(FEId)->MakeRefElementData(QuadFormula);

    uref = TFEDatabase3D::GetJointValues3D(BaseFuncts[FEId], QuadFormula, IJoint);

    for(k=0;k<N_Points;k++)
     {
      switch(RefElement)
      {
        case BFUnitHexahedron:
          ((THexaIsoparametric *)F_K)->GetTangentVectors(
                IJoint, p1[k], p2[k], a1, a2, a3, b1, b2, b3);
        break;

       case BFUnitTetrahedron:
       if(TDatabase::ParamDB->USE_ISOPARAMETRIC)
          ((TTetraIsoparametric *)F_K)->GetTangentVectors(
                IJoint, p1[k], p2[k], a1, a2, a3, b1, b2, b3);
       else
          ((TTetraAffin *)F_K)->GetTangentVectors(
                IJoint, p1[k], p2[k], a1, a2, a3, b1, b2, b3);
        break;
      } // endswitch


      n1 = a2*b3 - a3*b2;
      n2 = a3*b1 - a1*b3;
      n3 = a1*b2 - a2*b1;

//       OutPut("xi: " << p1[k] << " eta: " << p2[k] << endl);
//       OutPut("point: " << k << " t1: " << a1 << " " << a2 << " " << a3 << endl);
//       OutPut("point: " << k << " t2: " << b1 << " " << b2 << " " << b3 << endl);

      len = sqrt(n1*n1 + n2*n2 + n3*n3);
      n1 /= len;
      n2 /= len;
      n3 /= len;

     Mult = Weights[k]* len;
     Mass[1] +=Mult;

     for(l=0;l<N_JointDOF;l++)
      {
       local_j   = JointDOF[l];
       dof  =  DOF_LOW[l];
        Mass[0] += Mult*values[dof]*uref[k][local_j];

      } // for(l=0;l<N_JointDOF;l++)

    } // for(k=0;k<N_Points;k++)
  } //   for(i=0;i<N_Cells_low;i++)

 OutPut("Mass " << Mass[0] <<endl);
 OutPut("SurfaceArea " << Mass[1] <<endl);
}


void TetrameshGen(TDomain *Domain)
{
 //======================================================================
 // Tetgen for grid generation begin
 //======================================================================
  int i, j, k, l, N_Coord, N_FVert, *N_FVerts, *Facets;
  double *Vertices;
  int N, N_RootCells, N_Cells, CurrVertex, N_Vertices, ID, N_Faces, N_G, RefLevel=0;
  int CurrNeib, len1, len2, len3, maxEpV = 0, a, b, c, Neib[2], Neighb_tmp, CurrComp;
  int *Tetrahedrals, *PartMarker, *PointNeighb;
  double *Coordinates, N_x, N_y, N_z; 
  tetgenio In, Out;
  double Xmin = 1e10, Xmax = -1e10, Ymin = 1e10, Ymax = -1e10;  tetgenio::facet *F;
  double Zmin = 1e10, Zmax = -1e10, T[4]={0,0,0,0}, S[4]={0,0,0,0};  tetgenio::polygon *P;  
  double StartX, StartY, StartZ, BoundX, BoundY, BoundZ;  char *SMESH, line[100];
  TBaseCell **CellTree,  **SurfCellTree;  std::ostringstream opts;
  TGridCell **DelCell;
  TVertex **VertexDel, **NewVertices, **NewSurfVertices;  opts << " ";
  TBoundPart *BoundPart;
  TBdPlane **UpdateFaceParams;  SMESH = TDatabase::ParamDB->SMESHFILE;
  TJoint *Joint;
  TBoundComp3D *bdcomp; 
  TCollection *coll, *SurfColl;
  TBaseCell *cell;
  TBdSphere *UpdateParam;


 std::ifstream dat(SMESH);

  if (!dat)
  {
    cerr << "cannot open '" << SMESH << "' for input" << endl;
    exit(-1);
  }

  dat.getline (line, 99);

  // determine the number of vettices and alocate memory
  dat >> N_Vertices >> N_Coord;
  cerr << "N_Vertices: " << N_Vertices << endl;
  if(N_Coord!=3)  
   {
    cerr << "N_Vertices: " << N_Vertices << endl;
    cerr << "N_Coord must be 3 but it has: " << N_Coord << endl;
    exit(-1);
  }

  Vertices = new double[3*N_Vertices];

  for(i=0;i<N_Vertices; i++)
   {
    dat.getline (line, 99);
    dat >> j >>Vertices[3*i] >> Vertices[3*i+1] >> Vertices[3*i+2];
//     cout<< i << " vert X: " <<Vertices[3*i] <<endl;
   }

  dat.getline (line, 99);
  // determine the number of vettices and alocate memory
  dat >> N_Faces >> N_FVert;
  cerr << "N_Faces: " << N_Faces << endl;
//   N_FVert =4;
  Facets = new int[N_Faces*N_FVert];
  N_FVerts = new int[N_Faces];

 if(N_FVert==3)
  for(i=0;i<N_Faces; i++)
   {
    dat.getline (line, 99);
//     dat >> N_FVerts[i] >>Facets[3*i] >> Facets[3*i+1] >> Facets[3*i+2];
    dat >> N_FVerts[i] >>Facets[3*i] >> Facets[3*i+1] >> Facets[3*i+2];
//     cout<< i << " Facets X: " <<Facets[3*i] <<endl;
   }
  else if(N_FVert==4)
   for(i=0;i<N_Faces; i++)
    {
     dat.getline (line, 99);
//     dat >> N_FVerts[i] >>Facets[3*i] >> Facets[3*i+1] >> Facets[3*i+2];
     dat >> N_FVerts[i] >>Facets[4*i] >> Facets[4*i+1] >> Facets[4*i+2]>> Facets[4*i+3];
//     cout<< i << " Facets X: " <<Facets[3*i] <<endl;
    }
  else
   {
    cout<< "Number of face vertices should be 3 or 4, see TetrameshGen(TDomain *Domain) !!!!!!!!! "<< endl;
    exit(0); 
   }

  dat.close();

    BoundPart = Domain->GetBdPart(0);

    N = BoundPart->GetN_BdComps();
//     UpdateFaceParams = new TBdPlane *[N];
//     for(i=0; i<N; i++)
//      {
//       UpdateFaceParams[i] = (TBdPlane*)BoundPart->GetBdComp(i);
//      }
//    cout<<"N"<<N <<endl;

    UpdateParam = (TBdSphere*)BoundPart->GetBdComp(0);  // sphere domain
    UpdateParam->SetParams(0.0, 0.0, 0.0, 1.0);

    opts.seekp(std::ios::beg);
    opts<<'p'; // Tetrahedralize the PLC. Switches are chosen to read a PLC (p)
    opts<<"q"<<1.25; // quality mesh generation(q) with a specified quality bound
//     opts<<"a"<<0.1; // maximum volume constraint
//     opts<<'i'; // Inserts a list of additional points into mesh.
    opts<<'z'; // numbers all output items starting from zero
//     opts<<'d'; // Detect intersections of PLC facets.
    opts<<'f'; // Outputs all  faces (including non-boundary) 
    opts<<'e'; // Outputs a list of edges of the triangulation
//     opts<<'I'; // Suppresses mesh iteration numbers.
    opts<<'C'; // Checks the consistency of the final mesh.
//     opts<<'Q'; // Quiet: No terminal output except errors.
//     opts<<'g'; // Outputs mesh to .mesh file for viewing by Medit
    opts<<'Y'; // Suppresses boundary facets/segments splitting
    opts<<'V';  //verbose mode
    opts<<ends;


    In.numberofpoints = N_Vertices;
    In.pointlist = new double[3*In.numberofpoints];
    for(i=0;i<3*N_Vertices; i++)
     In.pointlist[i]= Vertices[i];

    In.numberoffacets = N_Faces;
    In.facetlist = new tetgenio::facet[In.numberoffacets];
    In.facetmarkerlist = new int[In.numberoffacets];
// cout<< " test main 1 " <<endl;
    for(i=0;i<N_Faces; i++)
     {
      F = &In.facetlist[i];
      F->numberofpolygons = 1;
      F->polygonlist = new tetgenio::polygon[F->numberofpolygons];
      F->numberofholes = 0;
      F->holelist = NULL;
      P = &F->polygonlist[0];
      P->numberofvertices = N_FVert;
      P->vertexlist = new int[P->numberofvertices];
      for(j=0;j<N_FVert;j++)
        P->vertexlist[j] = Facets[N_FVert*i + j];
      In.facetmarkerlist[i] = 1;
     }
// cout<< " test main " <<endl;

//     for(i=0;i<In.numberofpoints;i++)
//       OutPut(i<<" (x, y, z) =  "<<
//        In.pointlist[3*i]<<' '<<In.pointlist[3*i+1]<<' '<<In.pointlist[3*i+2]<<endl);

 // Calling  tetrahedralize function of 3dtetgen mesh generator
    tetrahedralize((char*)opts.str().c_str(), &In, &Out);

 //   output: coordinates of all vertices
//  for(i=0;i<Out.numberofpoints;i++)
//   OutPut(" (x, y, z) =  "<< Out.pointlist[3*i]<<' '<<Out.pointlist[3*i+1]<<' '<<Out.pointlist[3*i+2]<<endl);

    Domain->GetTreeInfo(CellTree,N_RootCells);
    coll = Domain->GetCollection(It_Finest, 0);
    N_Cells = coll->GetN_Cells();

//     cout<<"N_RootCells: "<<N_RootCells<<endl;
   // remove all existing vertices and joints
    VertexDel = new TVertex*[8*N_RootCells];
    CurrVertex = 0;

   for(i=0;i<N_Cells;i++)
     {
       cell = coll->GetCell(i);
       N_Faces = cell->GetN_Faces();
       N_Vertices = cell->GetN_Vertices();
       for(j=0;j<N_Faces;j++)
         {
          if(CurrVertex==0)
           {
               VertexDel[CurrVertex++] = cell->GetVertex(j);
            }
           else
            {
             ID = 0;
             for(k=0;k<CurrVertex;k++)
             if(VertexDel[k]==cell->GetVertex(j))
              {
               ID = 1; break;
              }
             if(ID!=1)
              {
               VertexDel[CurrVertex++] = cell->GetVertex(j);
              }
          } // else if(CurrVertex==0)

           ID = 0;
           for(k=0;k<CurrVertex;k++)
           if(VertexDel[k]==cell->GetVertex((j+1)%N_Vertices))
            {
             ID = 1; break;
            }
             if(ID!=1)
            {
             VertexDel[CurrVertex++] = cell->GetVertex((j+1)%N_Vertices);
            }
 
           ID = 0;
           for(k=0;k<CurrVertex;k++)
           if(VertexDel[k]==cell->GetVertex((j+2)%N_Vertices))
            {
             ID = 1; break;
            }
             if(ID!=1)
            {
             VertexDel[CurrVertex++] = cell->GetVertex((j+2)%N_Vertices);
            }
         if(N_Faces==6) // If old cell is hexahedrol
          {  
           ID = 0;
           for(k=0;k<CurrVertex;k++)
           if(VertexDel[k]==cell->GetVertex((j+4)%N_Vertices))
            {
             ID = 1; break;
            }
             if(ID!=1)
            {
             VertexDel[CurrVertex++] = cell->GetVertex((j+4)%N_Vertices);
            }        
        }
      } // for j
    } // for i

    for(i=0;i<CurrVertex;i++)
     delete VertexDel[i];
     delete [] VertexDel;
     OutPut(CurrVertex<<" vertices were deleted"<<endl);

  // remove all existing cells and joints

   for(i=0;i<N_RootCells;i++)
    delete (TGridCell*)CellTree[i];
    delete [] CellTree;
    OutPut(N_RootCells<<" cells were deleted"<<endl);


   N_RootCells = Out.numberoftetrahedra;


  // allocate auxillary fields
   Coordinates = Out.pointlist;
   Tetrahedrals = Out.tetrahedronlist;

  // generate new vertices
   N_G = Out.numberofpoints;
   NewVertices = new TVertex*[N_G];

   for (i=0;i<N_G;i++)
    {
      NewVertices[i] = new TVertex(Coordinates[3*i], Coordinates[3*i+1], Coordinates[3*i+2]);

      // set bounding box
      if (Coordinates[3*i] > Xmax) Xmax = Coordinates[3*i];
      if (Coordinates[3*i] < Xmin) Xmin = Coordinates[3*i];
      if (Coordinates[3*i+1] > Ymax) Ymax = Coordinates[3*i+1];
      if (Coordinates[3*i+1] < Ymin) Ymin = Coordinates[3*i+1];
      if (Coordinates[3*i+2] > Zmax) Zmax = Coordinates[3*i+2];
      if (Coordinates[3*i+2] < Zmin) Zmin = Coordinates[3*i+2];
   }

   // set bounding box
    StartX = Xmin;
    StartY = Ymin;
    StartZ = Zmin;
    BoundX = Xmax - Xmin;
    BoundY = Ymax - Ymin;
    BoundZ = Zmax - Zmin;


   Domain->SetBoundBox(StartX, StartY, StartZ, BoundX, BoundY, BoundZ);
//        cout<<Xmin <<"  "<<Ymin <<"  "<<Zmin<<endl;
//        cout<<Xmax <<"  "<<Ymax <<"  "<<Zmax<<endl;

   CellTree = new TBaseCell*[N_RootCells];
 //  output of each tetraheron vertex indices (four vertices for each)
//   for (i=0;i<N_RootCells;i++)
//        cout<< Tetrahedrals[4*i]<<"  "<<Tetrahedrals[4*i + 1]<<"  "
//          <<Tetrahedrals[4*i + 2]<<"  "<<Tetrahedrals[4*i + 3]<<endl;

   for (i=0;i<N_RootCells;i++)
   {
     CellTree[i] = new TMacroCell(TDatabase::RefDescDB[Tetrahedron],
                                    RefLevel);

     CellTree[i]->SetVertex(0, NewVertices[Tetrahedrals[4*i    ]]);
     CellTree[i]->SetVertex(1, NewVertices[Tetrahedrals[4*i + 1]]);
     CellTree[i]->SetVertex(2, NewVertices[Tetrahedrals[4*i + 2]]);
     CellTree[i]->SetVertex(3, NewVertices[Tetrahedrals[4*i + 3]]);

    CellTree[i]->SetClipBoard(i);
     ((TMacroCell *) CellTree[i])->SetSubGridID(0);
   }

   Domain->SetTreeInfo(CellTree, N_RootCells);


   // initialize iterators
   TDatabase::IteratorDB[It_EQ]->SetParam(Domain);
   TDatabase::IteratorDB[It_LE]->SetParam(Domain);
   TDatabase::IteratorDB[It_Finest]->SetParam(Domain);
   TDatabase::IteratorDB[It_Between]->SetParam(Domain);
   TDatabase::IteratorDB[It_OCAF]->SetParam(Domain);


//     search neighbours
   N_G = Out.numberofpoints;
   PointNeighb = new int[N_G];
   cout<<"numberofpoints "<<N_G<<endl;
   memset(PointNeighb, 0, N_G *SizeOfInt);

     for (i=0;i<4*N_RootCells;i++)
     PointNeighb[Tetrahedrals[i]]++;

   for (i=0;i<N_G;i++)
     if (PointNeighb[i] > maxEpV) maxEpV = PointNeighb[i];
   delete [] PointNeighb;

   cout<<"maxEpV "<< maxEpV<<endl;

   PointNeighb = new int[++maxEpV * N_G];

   memset(PointNeighb, 0, maxEpV*N_G*SizeOfInt);

    // every vertex contains "maxEpV" columns
    // for every vertex at first colomn contains the number of cells containing this vertex
    // at further columns we set the index of corresponding cells containing this vertex
    // cout<<"maxEpV*N_G "<<maxEpV*N_G<<endl;

   for(i=0;i<4*N_RootCells;i++)
    {
     j = Tetrahedrals[i]*maxEpV;
     PointNeighb[j]++;
     //cout<<"j + PointNeighb[j] " << j <<endl;
     PointNeighb[j + PointNeighb[j]] = i / 4;
    }
 //  output of PointNeighb columns for each point
//   for (i=0;i<N_G;i++)
//    {
//     for (j=0;j<maxEpV;j++)
//     cout<<"  "<< PointNeighb[i*maxEpV+j];
//     cout<<endl;
//    }

   // generate new faces 
   N_G =  Out.numberoftrifaces;
   cout<<"numberoftrifaces "<<N_G<<endl;
   for (i=0;i<N_G;i++)
   {
     a = Out.trifacelist[3*i];
     b = Out.trifacelist[3*i+1];
     c = Out.trifacelist[3*i+2];

//      cout<<"  "<< a<<"  "<< b<<"  "<< c<<endl;

     Neib[0] = -1;
     Neib[1] = -1;
     CurrNeib = 0;

     len1 = PointNeighb[a*maxEpV];
     len2 = PointNeighb[b*maxEpV];
     len3 = PointNeighb[c*maxEpV];

   // find the index of the cells containing current face with point indices a,b,c
    for (j=1;j<=len1;j++)
     {
       Neighb_tmp = PointNeighb[a*maxEpV + j];
        for (k=1;k<=len2;k++)
         {
          if (Neighb_tmp == PointNeighb[b*maxEpV + k])
           {
            for (l=1;l<=len3;l++)
             if (Neighb_tmp == PointNeighb[c*maxEpV + l])
             {
              Neib[CurrNeib++] = Neighb_tmp;
              break;
             }
           }
          }
       if (CurrNeib == 2) break;
     }
 //   cout<<"CurrNeib " << CurrNeib <<endl;
// cout<<"Out.trifacemarkerlist[i] : "<<Out.trifacemarkerlist[i]<<endl;
     if (Out.trifacemarkerlist[i]) // 0 for inner edges and Boundcomp+1 for Boundedge respect
      {

       CurrComp = Out.trifacemarkerlist[i] - 1;
//        cout<<"Boundary face CurrComp: "<<CurrComp<<endl;

       bdcomp = Domain->GetBdPart(0)->GetBdComp(CurrComp);

       if(bdcomp->GetTSofXYZ(NewVertices[a]->GetX(), NewVertices[a]->GetY(),
                             NewVertices[a]->GetY(), T[1], S[1]) ||
          bdcomp->GetTSofXYZ(NewVertices[b]->GetX(), NewVertices[b]->GetY(),
                             NewVertices[b]->GetY(), T[2], S[2]) ||
          bdcomp->GetTSofXYZ(NewVertices[c]->GetX(), NewVertices[c]->GetY(),
                             NewVertices[c]->GetY(), T[3], S[3])    )
         {
          cerr<<"Error: could not set parameter values"<<endl;
          OutPut(NewVertices[a]<<endl);
          OutPut(NewVertices[b]<<endl);
          OutPut(NewVertices[c]<<endl);
          exit(0);
         }

       if (CurrNeib == 2)
        {
         if(bdcomp->IsFreeBoundary())
           Joint = new TIsoInterfaceJoint3D(bdcomp, T, S,
                       CellTree[Neib[0]], CellTree[Neib[1]] );
          else
           Joint = new TInterfaceJoint3D(bdcomp, T, S, 
                       CellTree[Neib[0]], CellTree[Neib[1]] );
        }
       else
        {
         if(bdcomp->IsFreeBoundary())
           Joint = new TIsoBoundFace(bdcomp, T, S);
         else
           Joint = new TBoundFace(bdcomp, T, S);
        }
      }
     else
      {
// //       cout<<"Inner face"<<endl;
       if (CurrNeib != 2)
       cerr << "Error !!!!!!!! not enough neighbours!" << endl;

       Joint = new TJointEqN(CellTree[Neib[0]], CellTree[Neib[1]]);

      }

    // First element containing the current face
    // find the local index for the point 'a' on the cell
    for (j=0;j<4;j++)
      if (Tetrahedrals[4*Neib[0]+j] == a) break;

    // find the local index for the point 'b' on the cell
    for (k=0;k<4;k++)
      if (Tetrahedrals[4*Neib[0]+k] == b) break;

       // find the local index for the point 'c' on the cell
    for (l=0;l<4;l++)
      if (Tetrahedrals[4*Neib[0]+l] == c) break;   

     l = l*100 + k*10 + j;  

//      cout<<""<< l <<endl;

     switch (l) // j will contain the local index for the current face
      {
        case 210: case 21: case 102:
        case 120: case 12: case 201:
          j = 0;
          break;  
        case 310: case 31: case 103:
        case 130: case 13: case 301:
          j = 1;
          break;  
        case 321: case 132: case 213:
        case 231: case 123: case 312:
          j = 2;
          break;  
        case 230: case 23: case 302:
        case 320: case 32: case 203:
          j = 3;
          break; 

      default:
       Error("Unable to set the face !!!!!!!!!!!!" << endl);
       exit(0);
     }
      CellTree[Neib[0]]->SetJoint(j, Joint);

   if (Neib[1] != -1) // second element containing the current face
    {
          // find the local index for the point 'a' on the cell
    for (j=0;j<4;j++)
      if (Tetrahedrals[4*Neib[1]+j] == a) break;

    // find the local index for the point 'b' on the cell
    for (k=0;k<4;k++)
      if (Tetrahedrals[4*Neib[1]+k] == b) break;

       // find the local index for the point 'c' on the cell
    for (l=0;l<4;l++)
      if (Tetrahedrals[4*Neib[1]+l] == c) break;   

     l = l*100 + k*10 + j;  

//      cout<<""<< l <<endl;

     switch (l) // j will contain the local index for the current face
      {
        case 210: case 21: case 102:
        case 120: case 12: case 201:
          j = 0;
          break;  
        case 310: case 31: case 103:
        case 130: case 13: case 301:
          j = 1;
          break;  
        case 321: case 132: case 213:
        case 231: case 123: case 312:
          j = 2;
          break;  
        case 230: case 23: case 302:
        case 320: case 32: case 203:
          j = 3;
          break; 

      default:
       Error("Unable to set the face !!!!!!!!!!!!" << endl);
       exit(0);
      }
      CellTree[Neib[1]]->SetJoint(j, Joint);
     }

  if (Joint->GetType() == InterfaceJoint3D ||
      Joint->GetType() == IsoInterfaceJoint3D)
      {
        ((TInterfaceJoint3D*)Joint)->SetMapType();
        ((TInterfaceJoint3D*)(Joint))->CheckOrientation();
      }
      else 
       if (Joint->GetType() == JointEqN)
           Joint->SetMapType();

  }

  delete [] PointNeighb;


cout<<"Out.numberofpoints "<< Out.numberofpoints <<endl;
// cout<<"cout main tetra "<< N_Faces <<endl;
// exit(0);
 //======================================================================
 // Tetgen for grid generation end
 //======================================================================


delete [] Vertices;

}


void  Domain3DSurface(TDomain *Domain, TDomain *SurfDomain, int *Cell_No, int *Joint_No)
{
  int i,j,k,l,m,n,m1, N,N_G,N_SurfCells, maxlen, N_BoundVert, *EdgeList, *PointNeighb;
  int N_Cells, N_Joints, *EdgeListALL, *Triangles,ID,N_Edges, begin,end, n1, *Bd_Part;
  int maxEpV, a, b, Neib[2], CurrNeib, len1, len2, Neighb_tmp;
  TBaseCell *Me;
  TJoint *joint;
//   TIsoJointEqN *isojoint;
  TCollection *Coll;
  TVertex **BoundVetrex;
  TBaseCell  **SurfCellTree;
  const int *TmpFV,  *TmpLen;
  boolean check=FALSE; 
//   TBoundComp *BoundComp;
  std::ostringstream os;
      os << " ";
  double x0, y0, z0;


  Coll = Domain->GetCollection(It_Finest, 0);
  N_SurfCells = 0;
  N_Cells = Coll->GetN_Cells();
  for(i=0;i<N_Cells;i++)
  {
    Me = Coll->GetCell(i);
    k = Me->GetN_Faces();
    for(l=0;l<k;l++)
     {
       joint = Me->GetJoint(l);
       if(joint->GetType() == BoundaryFace || joint->GetType() == IsoBoundFace ||
           joint->GetType() == IsoInterfaceJoint3D )
         N_SurfCells++;
     } // endfor l
  } // endfor i

  OutPut("N_SurfCells: " << N_SurfCells << endl);

  if(k!=4)
   {
     cout << "Domain3DSurface  implemented only for tetrahedral !!!" << endl;
     exit(0);
    }

//    cout << " N_Cells " << N_SurfCells << endl;

//   Cell_No = new int[N_SurfCells];
//   Joint_No = new int[N_SurfCells];

  Bd_Part = new int[N_SurfCells];
  Triangles = new int[3*N_SurfCells];
  EdgeListALL = new int[6*N_SurfCells];
  BoundVetrex = new TVertex*[3*N_SurfCells];

  m1 = 0;
  N    = 0; 
  N_BoundVert = 0;
  for(i=0;i<N_Cells;i++)
   {
     Me = Coll->GetCell(i);
     k = Me->GetN_Faces();
     for(l=0;l<k;l++)
      {
        joint = Me->GetJoint(l);
        if(joint->GetType() == BoundaryFace || joint->GetType() == IsoBoundFace ||
           joint->GetType() == IsoInterfaceJoint3D )
         {
           Cell_No[N] = i;
           Joint_No[N] = l;
           Bd_Part[N++] =(((TBoundFace *)joint)->GetBoundComp())->GetID();
           Me->GetShapeDesc()->GetFaceVertex(TmpFV, TmpLen, maxlen);
 // getting freesurf vertices ---begin
           for(n=0;n<TmpLen[l];n++)
            {
              ID = 0;
              for(m=0;m<N_BoundVert;m++)
                if(BoundVetrex[m]==Me->GetVertex(TmpFV[maxlen*l+n]))
                 {
                   ID = 1;
                   Triangles[m1] = m;
                   m1++;
                   break; // vertex is already in the collection
                 }

              if(ID==0)
               {
                 Triangles[m1] = N_BoundVert;
                  m1++;
                 BoundVetrex[N_BoundVert++] = Me->GetVertex(TmpFV[maxlen*l+n]);
               }
          }
// getting freesurf vertices ---end
       }  //  if(joint->GetType() ==
     } // endfor l
   } // endfor i

//    cout << " N_vertices:  " << N_BoundVert << endl;
//    cout << " N_Cells " << N << endl;
//   for(i=0;i<N_SurfCells;i++)
//   cout << " Triangles:  " <<Triangles[3*i] << "    " << Triangles[3*i+1] << "    " << Triangles[3*i+2] << endl;
//   for(i=0;i<N_SurfCells;i++)
//  cout << " Cell_No:  " <<Cell_No[level][i] << " Joint_No  " << Joint_No[level][i] << endl;
// // write .smesh for tetgen
//       os.seekp(std::ios::beg);
// 
//       os << "Boundary.smesh" << ends;
//       std::ofstream dat(os.str().c_str());
//       if (!dat)
//        {
//         cerr << "cannot open file for output" << endl;
// //         exit(0);
//        }
//        dat << "# Surface mesh created for droplet by MooNMD" << endl;
// 
//       dat << N_BoundVert << " " <<  3<< endl;
//       for (i=0;i<N_BoundVert;i++)
//        {
//         BoundVetrex[i]->GetCoords(x0, y0, z0);
//         dat<< i << "  " <<  "  " << x0 << " " << "  " <<  y0 << " " << "  " <<  z0 <<endl;
//         }
// 
//       dat << N_SurfCells << " " <<  0<< endl;
//       for (i=0;i<N_SurfCells;i++)
//        {
//          dat<< 3 << "  " <<  "  " << Triangles[3*i    ] << " " << "  " <<  Triangles[3*i + 1] << " " << "  " <<  Triangles[3*i + 2] <<endl;
//         }
// 
//       dat <<  0<< endl;
//       dat <<  0<< endl;
//       dat.close();
// //       cout << endl;
//     cout << "Boundary wrote output into file " << endl;

   for (i=0;i<N_SurfCells;i++)
   {
    EdgeListALL[6*i    ] = Triangles[3*i    ];
    EdgeListALL[6*i + 1] = Triangles[3*i + 1];

    EdgeListALL[6*i + 2] = Triangles[3*i + 1];
    EdgeListALL[6*i + 3] = Triangles[3*i + 2];

    EdgeListALL[6*i + 4] = Triangles[3*i + 2];
    EdgeListALL[6*i + 5] = Triangles[3*i    ];
   }

// first Triangle
  N_Edges = 3;
  for(i=1;i<N_SurfCells;i++)
   {
    for(j=0;j<3;j++)
     {
       begin = EdgeListALL[6*i + 2*j ];
      if(j==2)
       end   = EdgeListALL[6*i    ];
      else
       end   = EdgeListALL[6*i +  2*j +1];

      for(k=0;k<i;k++)
       { 
        for(l=0;l<3;l++)
           {
             if(l==2)
              n1= 0;
             else  
              n1 = 2*l+1;

             if( (begin==EdgeListALL[6*k + 2*l]  &&   end == EdgeListALL[6*k + n1]) || 
                 (end  ==EdgeListALL[6*k + 2*l]  && begin == EdgeListALL[6*k + n1])   )
              {
               check = TRUE;
               break;
              }
           }// for l
         if(check)
          break;
        } // for k

         if(check)
          check= FALSE;
         else
           N_Edges++;

       } // for j
    }  // for i

//   cout << " Number of surface edges " << N_Edges << endl;
//  starting and ending vertices in each edge
  EdgeList = new int[2*N_Edges];
  N = 0;
  EdgeList[N++] = EdgeListALL[0];
  EdgeList[N++] = EdgeListALL[1];

  EdgeList[N++] = EdgeListALL[2];
  EdgeList[N++] = EdgeListALL[3];

  EdgeList[N++] = EdgeListALL[4];
  EdgeList[N++] = EdgeListALL[5];


  for(i=1;i<N_SurfCells;i++)
   {
    for(j=0;j<3;j++)
     {
       begin = EdgeListALL[6*i + 2*j];
      if(j==2)
       end   = EdgeListALL[6*i    ];
      else
       end   = EdgeListALL[6*i + 2*j +1];

      for(k=0;k<i;k++)
       { 
        for(l=0;l<3;l++)
           {
             if(l==2)
              n1= 0;
             else  
              n1 = 2*l+1;
             if( (begin==EdgeListALL[6*k + 2*l]  &&   end == EdgeListALL[6*k + n1]) || 
                 (end  ==EdgeListALL[6*k + 2*l]  && begin == EdgeListALL[6*k + n1])   )
              {
               check = TRUE;
               break;
              }
           }// for l
         if(check)
          break;
        } // for k

         if(check)
          check= FALSE;
         else
          {
           EdgeList[N++] = begin;
           EdgeList[N++] = end;
          }

       } // for j
    }  // for i

//   cout << " Number of surface edges " << N/2 << endl;

   SurfCellTree = new TBaseCell*[N_SurfCells];
   for (i=0;i<N_SurfCells;i++)
   {
    SurfCellTree[i] = new TMacroCell(TDatabase::RefDescDB[Triangle], 0);

    SurfCellTree[i]->SetVertex(0, BoundVetrex[Triangles[3*i    ]]);
    SurfCellTree[i]->SetVertex(1, BoundVetrex[Triangles[3*i + 1]]);
    SurfCellTree[i]->SetVertex(2, BoundVetrex[Triangles[3*i + 2]]);

    ((TMacroCell *) SurfCellTree[i])->SetSubGridID(0);
    ((TBaseCell *) SurfCellTree[i])->SetBd_Part(Bd_Part[i] );

   }

   SurfDomain->SetTreeInfo(SurfCellTree, N_SurfCells);

   TDatabase::IteratorDB[It_EQ]->SetParam(SurfDomain);
   TDatabase::IteratorDB[It_LE]->SetParam(SurfDomain);
   TDatabase::IteratorDB[It_Finest]->SetParam(SurfDomain);
   TDatabase::IteratorDB[It_Between]->SetParam(SurfDomain);
   TDatabase::IteratorDB[It_OCAF]->SetParam(SurfDomain);

  // search neighbours
  N_G = N_BoundVert;
  PointNeighb = new int[N_G];
  memset(PointNeighb, 0, N_G *SizeOfInt);

// degree of each vertex
  for (i=0;i<3*N_SurfCells;i++)
    PointNeighb[Triangles[i]]++;

  maxEpV = 0;
  for (i=0;i<N_G;i++)
    if (PointNeighb[i] > maxEpV) maxEpV = PointNeighb[i];

  delete [] PointNeighb;

//   cout << " Maximum degree of a vertex " << maxEpV << endl;
  PointNeighb = new int[++maxEpV * N_G];

  memset(PointNeighb, 0, maxEpV * N_G *SizeOfInt);

    // every vertex contains "maxEpV" columns
    // for every vertex at first colomn contains the number of cells containing this vertex
    // at further columns we set the index of the corresponding cells containing this vertex
 
  for(i=0;i<3*N_SurfCells;i++)
   {
     j = Triangles[i]*maxEpV;
     PointNeighb[j]++;
     PointNeighb[j + PointNeighb[j]] = i / 3;
   }

  // generate new edges
//   N_G = Out.numberofedges;
  N_G = N_Edges;

  for (i=0;i<N_G;i++)
   {
    a = EdgeList[2*i];
    b = EdgeList[2*i+1];
    Neib[0] = -1;
    Neib[1] = -1;
    CurrNeib = 0;

    len1 = PointNeighb[a*maxEpV];
    len2 = PointNeighb[b*maxEpV];

  //  cout<< a<< " test  " << b<<endl;
  //  cout<< len1<< " len1  " << len2<<endl;

  // find indexes of cells containing the current edge
   for (j=1;j<=len1;j++)
    {
      Neighb_tmp = PointNeighb[a*maxEpV + j];
       for (k=1;k<=len2;k++)
        if (Neighb_tmp == PointNeighb[b*maxEpV + k])
        {
          Neib[CurrNeib++] = Neighb_tmp;
          break;
        }
      if (CurrNeib == 2) break;
    }

//  cout << a << " Neib "<< b<<endl;

 // each edge of all triangles on a 3D surface is an interior edge
   if (CurrNeib != 2)
      cerr << "Error!!!!!!!! not enough neighbours!" << endl;

   joint = new TJointEqN(SurfCellTree[Neib[0]], SurfCellTree[Neib[1]]);

   // find the local index for the point 'a' on the cell
   for (j=0;j<3;j++)
     if (Triangles[3*Neib[0]+j] == a) break;

    // find the local index for the point 'b' on the cell
   for (k=0;k<3;k++)
     if (Triangles[3*Neib[0]+k] == b) break;

//   cout<< a <<  " test main " <<b <<endl;

   k = k*10 + j;

    switch (k)
    {
      case  1:
      case 10:
        j = 0;
        break;
      case 12:
      case 21:
        j = 1;
        break;
      case  2:
      case 20:
        j = 2;
        break;
    }

   SurfCellTree[Neib[0]]->SetJoint(j, joint);

      // find the local index for the point 'a' on the cell
      for (j=0;j<3;j++)
        if (Triangles[3*Neib[1]+j] == a) break;

      // find the local index for the point 'b' on the cell
      for (k=0;k<3;k++)
        if (Triangles[3*Neib[1]+k] == b) break;

      k = k*10 + j;

      switch (k) // j will contain the local index for the current
      {
        case  1:
        case 10:
          j = 0;
          break;
        case 12:
        case 21:
          j = 1;
          break;
        case  2:
        case 20:
          j = 2;
          break;
      }

    SurfCellTree[Neib[1]]->SetJoint(j, joint);
   }
  delete [] PointNeighb;
  delete []  Bd_Part;
  delete []  Triangles;
  delete []  EdgeListALL;
  delete []  EdgeList;
//   delete []  BoundVetrex;

//    cout <<   " Joint_No " <<Joint_No[0] <<endl;
//    cout << " Cell_No " <<  Cell_No[0] <<endl;

}


int main(int argc, char* argv[])
{
  TDomain *Domain = new TDomain();
  TDomain *SurfDomain = new TDomain();
  TDatabase *Database = new TDatabase();
  TFEDatabase3D *FEDatabase = new TFEDatabase3D();
  TFEDatabase2D *FEDatabase_surface = new TFEDatabase2D();
  TCollection *coll, *mortarcoll = NULL, *Surf_Coll;
  TDiscreteForm3D *DiscreteFormGalerkin;
  TDiscreteForm3D *DiscreteFormUpwind;
  TDiscreteForm3D *DiscreteFormGrid;
  TDiscreteForm2D *DiscreteFormSurfactant;
  TDiscreteForm3D *DiscreteFormNLGalerkin;
  TDiscreteForm3D *DiscreteFormNLUpwind;
  TDiscreteForm3D *DiscreteFormRHS;
  TFESpace2D *Surf_ScalarSpace;
  TFESpace3D *ScalarSpace, *VeloSpace, *GridSpace;
  TFESpace3D *fesp[4], *ferhs[6];
  TFESpace2D *fesp_low[4];
  TBaseCell *Me;
  TJoint *joint;
  TVertex **FreeVertices, **FreeIsoVertices, **IsoVertices;
  TIsoBoundFace **FreeBoundFace;
  TSquareStructure3D *sqGridStructure;
  TSquareStructure2D *sqSurfStructure;
  TSquareMatrix2D *MatricesSurf_Scalar_A, *MatricesSurf_Scalar_A_Old;
  TSquareMatrix2D *MatricesSurf_Scalar_M, *MatricesSurf_Scalar_M_Old;
  TSquareMatrix2D *MatricesSurf_Scalar_TDivU, *MatricesSurf_Scalar_TDivU_Old;
  TSquareMatrix2D *SQMATRICES_SURF[3];
  TFEVectFunct3D *GridVelocity, *RefGridPos, *AuxGridPos, *GridPos;
  TFEFunction2D *Surf_U1, *Surf_U2, *Surf_U3, *Surf_Scalar, *Surf_Scalar_Exact, *fefct_low[12];
  TFEFunction3D *Scalar, *Scalar_Exact, *U1, *U2, *U3, *fefct[6], *Grid_U1, *Grid_U2, *Grid_U3;
  TOutput3D *Output;
  TSquareMatrix3D *MatricesG11, *MatricesG12, *MatricesG13;
  TSquareMatrix3D *MatricesG21, *MatricesG22, *MatricesG23;
  TSquareMatrix3D *MatricesG31, *MatricesG32, *MatricesG33, *SQMATRICES_GRID[9];
  TAuxParam3D *aux;
  TFESpace2D *surfferhs[1];
  BoundCondFunct3D *GridBoundConditions[1];
  BoundValueFunct3D *GridBoundValues[1];
  TParDirectSolver *Grid_SMPSolver, *Scalar_SMPSolver;
  BoundCondFunct2D *SurfBoundConditions[1];
  BoundValueFunct2D *SurfBoundValues[1];
  MultiIndex3D AllDerivatives[4] = { D000, D100, D010, D001};

  bool MOVING_GRID=TRUE, Initialize_ScalarSolver=TRUE, Initialize_GridSolver=TRUE;

  int i, j, k, l, m, n, N_Cells, ret, N_S, N_Surf_S, maxlen, ID, ORDER;
  int **cellAndJointNo, *Cell_No, *Joint_No, VSP, N_SubSteps, N_G;
  int N_IsoVert, N_FreeVert, m1, m2, N_Vertices, methods, N_Rhs, len;
  int *BCN, *FJN, N_FreeJoints, N_Surf_U, img=1, time_discs, Eq_pos, Eq_pos1, Eq_pos2, Eq_pos3;
  int very_first_time=0, N_U, N_FESpaces, N_FESpaces_low, N_SquareMatrices;
  int *GridKCol, *GridRowPtr, N_GInner, N_GBoundaryNodes;
  int N_Eqn, N_Entries, N_MatEntries, *KCol, *RowPtr, begin, end;
  int *ScalarKCol, *ScalarRowPtr;
  const int *TmpFV, *TmpFV_aux, *TmpLen, *TmpLen_aux;

  double errors[9], RE_NR, **LinComb, *X, *Y, *Z, oldtau, end_time;
  double total_time, *Surf_rhs, *Surf_rhs_Old, *Surf_B, *U;
  double *Entries[9],t1, t2, *refpos, *auxpos, *pos;
  double *Surf_sol, *Surf_sol_Exact, *sol, *sol_Exact, tau, old_time, *Grid_U, *Grid_U_Old, *Surf_defect;
  double *SRHSs[1], gamma, *MatEntries, *Entries_Scalar;
  double olderror = 0, l_infty_l_2 = 0, l_infty_l_2_time=-4711.0, l_2_l_2u=0, l_2_H_1u=0;
  double olderror_l_2_l_2u=0, olderror_l_2_H_1u=0., hmin, hmax, Mass[2], scale=1., amp;

  char *PRM, *GEO;
  char *PsBaseName, *VtkBaseName;
  // Strings
  char ReadinDat[] = "readin.dat";
  char Name[] = "Surface";
  char PsiString[] = "Surfactant";
  char ExString[] = "Surfactant_Exact";
  char Surf_UString[] = "S_U";
  char UString[] = "u";
  char WString[] = "w";
  char refposString[] = "refpos";
  char auxposString[] = "auxpos";
  char posString[] = "pos";

  std::ostringstream os;

  N_SubSteps = GetN_SubSteps();
  oldtau = 1;
  end_time = TDatabase::TimeDB->ENDTIME;
  if (TDatabase::TimeDB->TIMESTEPLENGTH_CONTROL)
    time_discs = 2;
  else
    time_discs = 1;

//======================================================================
// read parameter file
//======================================================================
  total_time = GetTime();
  if(argc>=2)
    ret=Domain->ReadParam(argv[1]);
  else
    ret=Domain->ReadParam(ReadinDat);

  RE_NR=TDatabase::ParamDB->RE_NR;

  if(ret==-1)
  {
    exit(-1);
  }

  OpenFiles();
  OutFile.setf(std::ios::scientific);

  Database->CheckParameterConsistencyNSE();
  Database->WriteParamDB();
  Database->WriteTimeDB();
  ExampleFile();


  if( (TDatabase::ParamDB->DISCTYPE==2) )
  {
    OutPut("SDFEM does not work!" << endl);
    Error("SDFEM does not work!" << endl);
    exit(4711);
  }
  if( (TDatabase::ParamDB->DISCTYPE==5) )
  {
    OutPut("DISCTYPE 5 NOT IMPLEMENTED!" << endl);
    Error("DISCTYPE 5 NOT IMPLEMENTED!" << endl);
    exit(4711);
  }

  PRM = TDatabase::ParamDB->BNDFILE;
  GEO = TDatabase::ParamDB->GEOFILE;

  PsBaseName = TDatabase::ParamDB->PSBASENAME;
  VtkBaseName = TDatabase::ParamDB->VTKBASENAME;

//======================================================================
// creating discrete forms
//======================================================================

 InitializeDiscreteFormsForFreeSurface(
  DiscreteFormGalerkin,
  DiscreteFormUpwind,
  DiscreteFormNLGalerkin,
  DiscreteFormNLUpwind,
  DiscreteFormRHS,
  DiscreteFormGrid,
  DiscreteFormSurfactant,
  LinCoeffs, GridCoeffs, SurfCoeffs,
 TDatabase::ParamDB->NSTYPE);

//======================================================================
// read boundary parameterization and initialize coarse grid
//======================================================================
//   Domain->SurfGrid(PRM, GEO);
  Domain->Init(PRM, GEO);
//   smeshcreat();
//   SphereSmeshCreat();
  TetrameshGen(Domain);
//   ProjectBDPoints(Domain);

  for(i=0;i<TDatabase::ParamDB->UNIFORM_STEPS;i++)
   {
    Domain->RegRefineAll();
    ProjectBDPoints(Domain);
   }

  GridBoundConditions[0] = GridBoundCondition;
  GridBoundValues[0] = GridBoundValue;

  SurfBoundConditions[0] = SurfBoundCondition;
  SurfBoundValues[0] = SurfBoundValue;

  TDatabase::TimeDB->CURRENTTIME = TDatabase::TimeDB->STARTTIME;


  ORDER = 0;

  VSP = TDatabase::ParamDB->ANSATZ_ORDER;

  if (abs(VSP) > 20)
    ORDER = abs(VSP) - 20;
  else if ( abs(VSP) > 10)
    ORDER = abs(VSP) - 10;
  else
   ORDER = abs(VSP);

//======================================================================
// Generating FE spaces and allocating memory for their matrices
//======================================================================

   coll=Domain->GetCollection(It_Finest, 0);
   N_Cells = coll->GetN_Cells();

//     SmeshCreatQuad(Domain, SurfDomain, Cell_No_array, Joint_No_array, i);
   Cell_No = new int[N_Cells];
   Joint_No = new int[6*N_Cells];

   Domain3DSurface(Domain, SurfDomain, Cell_No, Joint_No);

   Surf_Coll = SurfDomain->GetCollection(It_Finest, 0);

   TDatabase::IteratorDB[It_EQ]->SetParam(Domain);
   TDatabase::IteratorDB[It_LE]->SetParam(Domain);
   TDatabase::IteratorDB[It_Finest]->SetParam(Domain);
   TDatabase::IteratorDB[It_Between]->SetParam(Domain);
   TDatabase::IteratorDB[It_OCAF]->SetParam(Domain);

   OutPut("cells " << N_Cells<< endl);
   cout << endl << endl;

   VeloSpace = new TFESpace3D(coll , Name, PsiString, GridBoundCondition,
                              ContP_USpace, TDatabase::ParamDB->ANSATZ_ORDER);
   N_U = VeloSpace->GetN_DegreesOfFreedom();

   GridSpace = new TFESpace3D(coll, Name, WString, GridBoundCondition, ContP_USpace, 1);
   N_G = GridSpace->GetN_DegreesOfFreedom();
   N_GInner = GridSpace->GetN_Inner();
   N_GBoundaryNodes = N_G - N_GInner;

   ScalarSpace = new TFESpace3D(coll , Name, PsiString, GridBoundCondition,
                                ContP_USpace, TDatabase::ParamDB->ANSATZ_ORDER);
   N_S =  ScalarSpace->GetN_DegreesOfFreedom();


   Surf_ScalarSpace = new TFESpace2D(Surf_Coll, Name, PsiString, SurfBoundCondition,
                                     TDatabase::ParamDB->ANSATZ_ORDER, NULL);
   N_Surf_S = Surf_ScalarSpace->GetN_DegreesOfFreedom();

   // For Freesurfint
 if(ORDER==2)
  {
   BCN = new int[6*N_Cells]; // BOUNDARY CELL NUMBER
   FJN = new int[6*N_Cells]; // FACE NUMBER IN BOUNDARY CELL
   N_FreeJoints = 0;

   for(j=0;j<N_Cells;j++)
    {
      Me = coll->GetCell(j);
      k = Me->GetN_Faces();
      for(l=0;l<k;l++)
       {
        joint = Me->GetJoint(l);
        if(joint->GetType() == IsoBoundFace)
         {
          BCN[N_FreeJoints] = j;
          FJN[N_FreeJoints] = l;
          N_FreeJoints++;
         }
      } // endfor l
    } // endfor

   N_IsoVert = 5*N_FreeJoints;

   FreeVertices = new TVertex*[4*N_FreeJoints];
   FreeIsoVertices = new TVertex*[N_IsoVert];
   FreeBoundFace = new TIsoBoundFace*[N_FreeJoints];
   m1=0;
   m2=0;
   N_FreeVert=0;


// Second order: one face bubble vertex (not needed for isoparametric approx.) and 3 edge mid vertices
// Generating 3 edge mid vertices for Second order
   if(ORDER==2)
    {
     LinComb =  new double*[5];
     for(n=0;n<5;n++)
      LinComb[n] =  new double[4];

     X = new double[4];
     Y = new double[4];
     Z = new double[4];

     if(k == 4) //tetrehedrol
      {
       N_Vertices = 4; // without face vertices, see IsoBoundFace.C
       LinComb[0][0]=0.5; LinComb[1][0]=0.5; LinComb[2][0]=0.0; LinComb[3][0]=1.0/3.0;
       LinComb[0][1]=0.5; LinComb[1][1]=0.0; LinComb[2][1]=0.5; LinComb[3][1]=1.0/3.0;
       LinComb[0][2]=0.0; LinComb[1][2]=0.5; LinComb[2][2]=0.5; LinComb[3][2]=1.0/3.0;
      }
     else //hexahedrol
      {
       N_Vertices = 5;
       LinComb[0][0]=0.5; LinComb[1][0]=0.5; LinComb[2][0]=0.25; LinComb[3][0]=0.0; LinComb[4][0]=0.0;
       LinComb[0][1]=0.0; LinComb[1][1]=0.5; LinComb[2][1]=0.25; LinComb[3][1]=0.0; LinComb[4][1]=0.5;
       LinComb[0][2]=0.0; LinComb[1][2]=0.0; LinComb[2][2]=0.25; LinComb[3][2]=0.5; LinComb[4][2]=0.5;
       LinComb[0][3]=0.5; LinComb[1][3]=0.0; LinComb[2][3]=0.25; LinComb[3][3]=0.5; LinComb[4][3]=0.0;
      }
     }

    for(j=0;j<N_Cells;j++)
     {

      Me = coll->GetCell(j);
      k = Me->GetN_Faces();
      for(l=0;l<k;l++)
       {
        joint = Me->GetJoint(l);
        if(joint->GetType() == IsoBoundFace)
         {
          (Me->GetShapeDesc())->GetFaceVertex(TmpFV, TmpLen, maxlen);

          for(n=0;n<TmpLen[l];n++)
           {
            if(N_FreeVert==0)
             FreeVertices[N_FreeVert++] = Me->GetVertex(TmpFV[maxlen*l+n]);

            ID = 0;
            for(m=0;m<N_FreeVert;m++)
             {
              if(FreeVertices[m]==Me->GetVertex(TmpFV[maxlen*l+n]))
               {
                ID = 1;
                break;
               }
             }

           if(ID!=1)
            FreeVertices[N_FreeVert++] = Me->GetVertex(TmpFV[maxlen*l+n]);
          }

         FreeBoundFace[m1] = (TIsoBoundFace *)joint;
         if(ORDER==2)
          {
           for(n=0;n<TmpLen[l];n++)
            Me->GetVertex(TmpFV[maxlen*l+n])->GetCoords(X[n], Y[n], Z[n]);

           if(k == 4) //tetrehedrol
            FreeBoundFace[m1]->GenVert(N_Vertices, TmpLen[l], LinComb, X, Y, Z);
           else //hexahedrol
            FreeBoundFace[m1]->GenHexaVert(N_Vertices, TmpLen[l], LinComb, X, Y, Z);
          }

         IsoVertices = FreeBoundFace[m1]->GetVertices();
         N_IsoVert = FreeBoundFace[m1]->GetN_Vertices(); // iso vertices in each face
	 m1++;
         for(n=0;n<N_IsoVert;n++)
          FreeIsoVertices[m2++] = IsoVertices[n];
        }
      }
    } //  for(j=0;j<N_Cells;j++)

//    for(n=0;n<5;n++)
//     delete [] LinComb[n];

   delete [] LinComb;

   delete [] X;
   delete [] Y;
   delete [] Z;

   N_IsoVert = m2;  // total Iso vertices
   cout << "N_FreeJoints: " << N_FreeJoints << endl;
   cout << "N_FreeVertices: " << N_FreeVert << endl;
   cout << "N_IsoVertices: " << N_IsoVert << endl; //each freesurf face contains 4 vertices
 } // if ORDER==2


   sqSurfStructure = new TSquareStructure2D(Surf_ScalarSpace);
   sqSurfStructure->Sort();
   sqSurfStructure->Sort_ForDirectSolver(); // necessary for ParDirectSolver

   MatricesSurf_Scalar_A = new TSquareMatrix2D(sqSurfStructure);
   MatricesSurf_Scalar_A_Old = new TSquareMatrix2D(sqSurfStructure);

   MatricesSurf_Scalar_M = new TSquareMatrix2D(sqSurfStructure);
   MatricesSurf_Scalar_M_Old = new TSquareMatrix2D(sqSurfStructure);

   MatricesSurf_Scalar_TDivU = new TSquareMatrix2D(sqSurfStructure);
   MatricesSurf_Scalar_TDivU_Old = new TSquareMatrix2D(sqSurfStructure);

  #ifdef __MOVINGGRID__
   sqGridStructure = new TSquareStructure3D(GridSpace);
   sqGridStructure->Sort_ForDirectSolver(); // necessary for ParDirectSolver

   MatricesG11 = new TSquareMatrix3D(sqGridStructure);
   MatricesG12 = new TSquareMatrix3D(sqGridStructure);
   MatricesG13 = new TSquareMatrix3D(sqGridStructure);
   MatricesG21 = new TSquareMatrix3D(sqGridStructure);
   MatricesG22 = new TSquareMatrix3D(sqGridStructure);
   MatricesG23 = new TSquareMatrix3D(sqGridStructure);
   MatricesG31 = new TSquareMatrix3D(sqGridStructure);
   MatricesG32 = new TSquareMatrix3D(sqGridStructure);
   MatricesG33 = new TSquareMatrix3D(sqGridStructure);
  #endif

   OutPut("dof Surf_Scalar : "<< setw(20) << N_Surf_S << endl);
   OutPut("dof Scalar : "<< setw(20) << N_S << endl);
   OutPut("Grid Velo  dof : "<< setw(20) <<3* N_G << endl);
//    OutPut("Velo dof : "<< setw(20) << 3*N_U << endl);

   Surf_sol = new double[N_Surf_S];
   Surf_sol_Exact = new double[N_Surf_S];
   Surf_rhs = new double[N_Surf_S];
   Surf_rhs_Old = new double[N_Surf_S];
   Surf_B = new double [N_Surf_S];  // working rhs array
   Surf_defect = new double[N_Surf_S];
   sol = new double[N_S];
   sol_Exact = new double[N_S];
   U = new double[3*N_U]; // 3D velocity in domain
   Grid_U = new double[3*N_G];
   Grid_U_Old = new double[3*N_G];

   refpos = new double[3*N_G];
   auxpos = new double[3*N_G];
   pos = new double[3*N_G];

   memset(Surf_sol, 0, N_Surf_S*SizeOfDouble);
   memset(Surf_sol_Exact, 0, N_Surf_S*SizeOfDouble);
   memset(Surf_rhs, 0, N_Surf_S*SizeOfDouble);
   memset(Surf_rhs_Old, 0, N_Surf_S*SizeOfDouble);
   memset(sol, 0, N_S*SizeOfDouble);
   memset(sol_Exact, 0, N_S*SizeOfDouble);
   memset(U, 0, N_U*SizeOfDouble);
   memset(Grid_U, 0, 3*N_G*SizeOfDouble);
   memset(Grid_U_Old, 0, 3*N_G*SizeOfDouble);
   memset(refpos, 0, 3*N_G*SizeOfDouble);
   memset(auxpos, 0, 3*N_G*SizeOfDouble);
   memset(pos, 0, 3*N_G*SizeOfDouble);

//  scalar
   Surf_Scalar = new TFEFunction2D(Surf_ScalarSpace, PsiString, PsiString, Surf_sol, N_Surf_S);
   Surf_Scalar_Exact = new TFEFunction2D(Surf_ScalarSpace, ExString, ExString, Surf_sol_Exact, N_Surf_S);

//  scalar in the whole domain for output in VTK
   Scalar = new TFEFunction3D(ScalarSpace, PsiString, PsiString, sol, N_S);
   Scalar_Exact = new TFEFunction3D(ScalarSpace, ExString, ExString, sol_Exact, N_S);
   SurfInterpolate(Scalar, Surf_Scalar, Cell_No, Joint_No, InitialSall);
   MapDomainToSurf(Scalar, Surf_Scalar, Cell_No, Joint_No);
   memset(sol, 0, N_S*SizeOfDouble);

   SurfInterpolate(Scalar_Exact, Surf_Scalar_Exact, Cell_No, Joint_No, InitialSall);
   MapDomainToSurf(Scalar_Exact, Surf_Scalar_Exact, Cell_No, Joint_No);
   memset(sol_Exact, 0, N_S*SizeOfDouble);

// velocity in the entire domain
   U1 = new TFEFunction3D(VeloSpace, UString, UString, U,       N_U);
   U2 = new TFEFunction3D(VeloSpace, UString, UString, U+N_U,   N_U);
   U3 = new TFEFunction3D(VeloSpace, UString, UString, U+2*N_U, N_U);

// Grid velocity in the entire domain
   GridVelocity = new TFEVectFunct3D(GridSpace, WString, WString, Grid_U, N_G, 3);
   Grid_U1 = GridVelocity->GetComponent(0);
   Grid_U2 = GridVelocity->GetComponent(1);
   Grid_U3 = GridVelocity->GetComponent(2);

   GridPos = new TFEVectFunct3D(GridSpace, posString, WString, pos, N_G, 3);
   AuxGridPos = new TFEVectFunct3D(GridSpace, auxposString, WString, auxpos, N_G, 3);
   RefGridPos = new TFEVectFunct3D(GridSpace, refposString, WString, refpos, N_G, 3);


   GridPos->GridToData();
//    RefGridPos->GridToData();
//    AuxGridPos->GridToData();
   memcpy(auxpos, pos, 3*N_G*SizeOfDouble);
   memcpy(refpos, pos, 3*N_G*SizeOfDouble);

   Entries[0] = MatricesG11->GetEntries();
   Entries[1] = MatricesG12->GetEntries();
   Entries[2] = MatricesG13->GetEntries();
   Entries[3] = MatricesG21->GetEntries();
   Entries[4] = MatricesG22->GetEntries();
   Entries[5] = MatricesG23->GetEntries();
   Entries[6] = MatricesG31->GetEntries();
   Entries[7] = MatricesG32->GetEntries();
   Entries[8] = MatricesG33->GetEntries();

   GridKCol = sqGridStructure->GetKCol();
   GridRowPtr = sqGridStructure->GetRowPtr();


// *****************************************************************************************************
// Set Fespaces for output
//  ****************************************************************************************************

   if(TDatabase::ParamDB->WRITE_VTK )
     {
      // output only velocity and pressure
        Output = new TOutput3D(3, 3 , 1, 1, Domain);

  #ifdef __MOVINGGRID__
//           Output->AddFEVectFunct(GridVelocity);
  #endif
        Output->AddFEFunction(Scalar);
//         Output->AddFEFunction(Scalar_Exact);

        os.seekp(std::ios::beg);
        Output->AddParameter(TDatabase::TimeDB->CURRENTTIME,os.str().c_str());
      }

// assemble matrices for grid movement
// reference domain is the initial domain
  #ifdef __MOVINGGRID__
        t1 = GetTime();
    // assemble matrix for grid moving
        fesp[0] = GridSpace;
        SQMATRICES_GRID[0] = MatricesG11;
        SQMATRICES_GRID[0]->Reset();
        SQMATRICES_GRID[1] = MatricesG12;
        SQMATRICES_GRID[1]->Reset();
        SQMATRICES_GRID[2] = MatricesG13;
        SQMATRICES_GRID[2]->Reset();
        SQMATRICES_GRID[3] = MatricesG21;
        SQMATRICES_GRID[3]->Reset();
        SQMATRICES_GRID[4] = MatricesG22;
        SQMATRICES_GRID[4]->Reset();
        SQMATRICES_GRID[5] = MatricesG23;
        SQMATRICES_GRID[5]->Reset();
        SQMATRICES_GRID[6] = MatricesG31;
        SQMATRICES_GRID[6]->Reset();
        SQMATRICES_GRID[7] = MatricesG32;
        SQMATRICES_GRID[7]->Reset();
        SQMATRICES_GRID[8] = MatricesG33;
        SQMATRICES_GRID[8]->Reset();

        aux = new TAuxParam3D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL,
                              0, NULL);

        Assemble3D(1, fesp,
                   9, SQMATRICES_GRID,
                   0, NULL,
                   0, NULL, NULL,
                   DiscreteFormGrid,
                   GridBoundConditions,
                   GridBoundValues,
                   aux);

         delete aux;

         memset(Entries[1] + GridRowPtr[N_G-N_GBoundaryNodes], 0,
                N_GBoundaryNodes*SizeOfDouble);
         memset(Entries[2] + GridRowPtr[N_G-N_GBoundaryNodes], 0,
                N_GBoundaryNodes*SizeOfDouble);

         memset(Entries[3] + GridRowPtr[N_G-N_GBoundaryNodes], 0,
                N_GBoundaryNodes*SizeOfDouble);
         memset(Entries[5] + GridRowPtr[N_G-N_GBoundaryNodes], 0,
                N_GBoundaryNodes*SizeOfDouble);

         memset(Entries[6] + GridRowPtr[N_G-N_GBoundaryNodes], 0,
                N_GBoundaryNodes*SizeOfDouble);
         memset(Entries[7] + GridRowPtr[N_G-N_GBoundaryNodes], 0,
                N_GBoundaryNodes*SizeOfDouble);

         t2 = GetTime();
//          cout << "Grid assembling done"<< endl;
//          OutPut("Time for Grid assembling: " << t2-t1 << endl);


         /** Construct SMP Direct Solver for grid velocity **/
          // Sort_ForDirectSolver() should be called before this construction direct solver!
          // collect data of the global system
//             OutPut("memory after : " << setw(10) << GetMemory() << endl);

         if(Initialize_GridSolver)
          {
            N_Entries = GridRowPtr[N_G];

            N_Eqn = 3*N_G;
            N_MatEntries = 9*GridRowPtr[N_G];
            MatEntries = new double[N_MatEntries];

            KCol = new int[N_MatEntries];
            RowPtr = new int[N_Eqn+1];

            RowPtr[0] = 0 + 1; // fortran style for ParDiso direct solver

            for(i=0;i<N_G;i++)
             {
              begin = GridRowPtr[i];
              end = GridRowPtr[i+1];
              len = end - begin;

              Eq_pos1 = 3*begin;

              for(j=begin;j<end;j++)
               {
                MatEntries[Eq_pos1] = Entries[0][j];
                KCol[Eq_pos1] = GridKCol[j] + 1; // fortran style

                MatEntries[len + Eq_pos1] = Entries[1][j];
                KCol[len + Eq_pos1] = N_G + GridKCol[j] + 1; // fortran style

                MatEntries[2*len + Eq_pos1] = Entries[2][j];
                KCol[2*len + Eq_pos1] = 2*N_G + GridKCol[j] + 1; // fortran style

                Eq_pos1++;
               }
              RowPtr[i+1] = 3*begin + 3*len + 1; // fortran style for ParDiso direct solver

              Eq_pos2 = 3*N_Entries + 3*begin;
              for(j=begin;j<end;j++)
               {
                MatEntries[Eq_pos2] = Entries[3][j];
                KCol[Eq_pos2] = GridKCol[j] + 1; // fortran style

                MatEntries[len + Eq_pos2] = Entries[4][j];
                KCol[len + Eq_pos2] = N_G + GridKCol[j] + 1; // fortran style

                MatEntries[2*len + Eq_pos2] = Entries[5][j];
                KCol[2*len + Eq_pos2] = 2*N_G + GridKCol[j] + 1; // fortran style

                Eq_pos2++;
               }
              RowPtr[N_G + i+1] = 3*N_Entries + 3*begin + 3*len + 1;

              Eq_pos3 = 6*N_Entries + 3*begin;
              for(j=begin;j<end;j++)
               {
                MatEntries[Eq_pos3] = Entries[6][j];
                KCol[Eq_pos3] = GridKCol[j] + 1; // fortran style

                MatEntries[len + Eq_pos3] = Entries[7][j];
                KCol[len + Eq_pos3] = N_G + GridKCol[j] + 1; // fortran style

                MatEntries[2*len + Eq_pos3] = Entries[8][j];
                KCol[2*len + Eq_pos3] = 2*N_G + GridKCol[j] + 1; // fortran style

                Eq_pos3++;
               }
              RowPtr[2*N_G + i+1] = 6*N_Entries + 3*begin + 3*len + 1;

             } // for(i=0;i<N_G;i++)

//             cout<< " RowPtr[2*N_G + i+1] " << RowPtr[2*N_G + i] <<endl;
//             cout<< " N_MatEntries " << N_MatEntries <<endl;
//             cout<< " N_Entries " << N_Entries <<endl;

            Grid_SMPSolver = new TParDirectSolver(N_Eqn, RowPtr, KCol, N_MatEntries);
            Grid_SMPSolver->SymbolicFactorize(MatEntries);

            delete [] RowPtr;
            delete [] KCol;
            delete [] MatEntries;

           Initialize_GridSolver = FALSE;

           t1 = GetTime();
           OutPut("Time for Grid_SMPSolver initialize: " << t1-t2 << endl);

          }
/** Construct SMP Direct Solver for grid velocity - end**/


// to get the initial matrix and rhs
//        W^{n+1} of the domain \Omega_t^{n+1}
       GetCurrentGridVelo(RefGridPos, GridPos, AuxGridPos, GridVelocity,
                          TDatabase::TimeDB->STARTTIME, TDatabase::TimeDB->TIMESTEPLENGTH,
                          SQMATRICES_GRID, GridKCol, GridRowPtr,
                          Entries, Grid_SMPSolver);


  #endif


       oldtau = TDatabase::TimeDB->TIMESTEPLENGTH;
//         cout << "Surfactant assembling begin"<< endl;
        N_FESpaces = 2;
        fesp[0] = ScalarSpace;  // scalar space in entire domain
        fesp[1] = GridSpace;  // velocity space

        fefct[0] = Grid_U1;  // u1 = Grid_U1
        fefct[1] = Grid_U2;  //u2 = Grid_U2
        fefct[2] = Grid_U3;  //u3 = Grid_U2

        N_FESpaces_low = 1;
        fesp_low[0] = Surf_ScalarSpace;
        fefct_low[0] = Surf_Scalar;

        N_Rhs =1;
        surfferhs[0] = Surf_ScalarSpace;
        SRHSs[0] = Surf_rhs;
        memset(SRHSs[0], 0, N_Surf_S*SizeOfDouble);

        N_SquareMatrices = 3;
        SQMATRICES_SURF[0] = MatricesSurf_Scalar_A;
        SQMATRICES_SURF[0]->Reset();
        SQMATRICES_SURF[1] = MatricesSurf_Scalar_M;
        SQMATRICES_SURF[1]->Reset();
        SQMATRICES_SURF[2] = MatricesSurf_Scalar_TDivU;
        SQMATRICES_SURF[2]->Reset();

        AssembleSurf2D(N_FESpaces, fesp, fefct, N_FESpaces_low, fesp_low,
                       N_SquareMatrices, SQMATRICES_SURF, N_Rhs, SRHSs, 
                       surfferhs, Cell_No, Joint_No,  TDatabase::TimeDB->STARTTIME,
                       TDatabase::TimeDB->TIMESTEPLENGTH, scale);

        memcpy(Surf_rhs_Old, Surf_rhs, N_Surf_S*SizeOfDouble);

        MatricesSurf_Scalar_A_Old->Reset();
        MatAdd(MatricesSurf_Scalar_A_Old, MatricesSurf_Scalar_A, 1.);
        MatricesSurf_Scalar_M_Old->Reset();
        MatAdd(MatricesSurf_Scalar_M_Old, MatricesSurf_Scalar_M, 1.);
        MatricesSurf_Scalar_TDivU_Old->Reset();
        MatAdd(MatricesSurf_Scalar_TDivU_Old, MatricesSurf_Scalar_TDivU, 1.);

   if(TDatabase::ParamDB->WRITE_VTK)
    {
//        MapDomainToSurf(U1Array[mg_level-1],SU1Array[mg_level-1], Cell_No_array[mg_level-1], Joint_No_array[mg_level-1]);
//        MapSurfToDomain(SU1Array[mg_level-1], Surfactant, Cell_No_array[mg_level-1], Joint_No_array[mg_level-1]);
       MapSurfToDomain(Surf_Scalar, Scalar, Cell_No, Joint_No);

       os.seekp(std::ios::beg);
       if(img<10) os << VtkBaseName<<".0000"<<img<<".vtk" << ends;
       else if(img<100) os << VtkBaseName<<".000"<<img<<".vtk" << ends;
       else if(img<1000) os << VtkBaseName<<".00"<<img<<".vtk" << ends;
       else if(img<10000) os << VtkBaseName<<".0"<<img<<".vtk" << ends;
       else  os << VtkBaseName<<"."<<img<<".vtk" << ends;
       Output->WriteVtk(os.str().c_str());

       img++;
      }

  N_SubSteps = GetN_SubSteps();
//   oldtau = 1;
  end_time = TDatabase::TimeDB->ENDTIME;
  if (TDatabase::TimeDB->TIMESTEPLENGTH_CONTROL)
    time_discs = 2;
  else
    time_discs = 1;

  m=0;

  Surf_Coll->GetHminHmax(&hmin,&hmax);
  OutPut("Surface h_min : " << hmin << "Surface h_max : " << hmax << endl);


//   TDatabase::TimeDB->TIMESTEPLENGTH= hmin*hmin;
//   OutPut("TIMESTEPLENGTH : " <<hmax*hmax << endl);

//   coll->GetHminHmax(&hmin,&hmax);
//   OutPut("h_min : " << hmin << " h_max : " << hmax << endl);


      OutPut(endl << "CURRENT TIME: ");
      OutPut(TDatabase::TimeDB->CURRENTTIME << endl);
      GetSurfMass(ScalarSpace, Surf_Scalar, Cell_No, Joint_No, Mass);

      MapSurfToDomain(Surf_Scalar, Scalar, Cell_No, Joint_No);
      GetSurfErrors(Scalar, Surf_Scalar, errors, Cell_No, Joint_No, ExactSall);
      OutPut( "L2: " << errors[0] << endl);
      OutPut( "H1-semi: " << errors[1] << endl);

 amp=TDatabase::ParamDB->P15;
//======================================================================
// start of time cycle
//======================================================================
  while(TDatabase::TimeDB->CURRENTTIME< end_time)
  {                              // time cycle
    m++;
    TDatabase::TimeDB->INTERNAL_STARTTIME = TDatabase::TimeDB->CURRENTTIME;

    for (methods=0;methods<time_discs;methods++)
    {
      if (time_discs==2)
      {
        if (methods==0) // fractional-step-theta-scheme
          TDatabase::TimeDB->TIME_DISC = 3;

        else           // crank nicolson scheme
        {              // take solution of first scheme as initial iterate
          TDatabase::TimeDB->TIME_DISC = 2;
          TDatabase::TimeDB->CURRENTTIME = TDatabase::TimeDB->INTERNAL_STARTTIME;
        }
        N_SubSteps = GetN_SubSteps();
      }

      for(l=0;l<N_SubSteps;l++)      // sub steps of fractional step theta
       {
        SetTimeDiscParameters();

        if (m==1)
        {
          OutPut("Theta1: " << TDatabase::TimeDB->THETA1<< endl);
          OutPut("Theta2: " << TDatabase::TimeDB->THETA2<< endl);
          OutPut("Theta3: " << TDatabase::TimeDB->THETA3<< endl);
          OutPut("Theta4: " << TDatabase::TimeDB->THETA4<< endl);
        }

        TDatabase::ParamDB->P14 = TDatabase::TimeDB->CURRENTTIME;
        old_time = TDatabase::TimeDB->CURRENTTIME;

        tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
        TDatabase::TimeDB->CURRENTTIME += tau;

        OutPut(endl << "CURRENT TIME: ");
        OutPut(TDatabase::TimeDB->CURRENTTIME << endl);

/* moving boundary is know, so the computation is performed on the new domain \Omega_t^{n+1}*/
/* Move the domain \Omega_t^{n+1} =  \Omega_t^n + dtW^n */

        #ifdef __MOVINGGRID__

       MoveGrid(GridPos, GridVelocity, TDatabase::TimeDB->CURRENTTIME, tau);
//     get the current time step domain velocity w^{n+1}

       GetCurrentGridVelo(RefGridPos, GridPos, AuxGridPos, GridVelocity,
                          TDatabase::TimeDB->CURRENTTIME, tau,
                          SQMATRICES_GRID, GridKCol, GridRowPtr,
                          Entries, Grid_SMPSolver);


        // scaling due the old sol and old rhs are in the previous timestep domain
        GetScale(TDatabase::TimeDB->CURRENTTIME, old_time, scale);
//         cout<< " scale " << scale << endl;

//       get the old solution u^n on the new domain \Omega_t^{n+1}
//       memset(sol_Exact, 0, N_S*SizeOfDouble);
//       memset(Surf_sol_Exact, 0, N_Surf_S*SizeOfDouble);
//       SurfInterpolate(Scalar_Exact, Surf_Scalar_Exact, Cell_No, Joint_No, ExactSall_oldStep);
//       MapDomainToSurf(Scalar_Exact, Surf_Scalar_Exact, Cell_No, Joint_No);
//       memcpy(Surf_sol,  Surf_sol_Exact, N_Surf_S*SizeOfDouble);


// get the previous domain solution on new domain
        MatAdd(MatricesSurf_Scalar_M_Old, MatricesSurf_Scalar_TDivU_Old, tau);
        memset(Surf_defect, 0, N_Surf_S*SizeOfDouble);
        MatVectActive(MatricesSurf_Scalar_M_Old, Surf_sol, Surf_defect);
        MatAdd(MatricesSurf_Scalar_M_Old, MatricesSurf_Scalar_TDivU_Old, -tau);

        DirectSolver(MatricesSurf_Scalar_M_Old, Surf_defect, Surf_sol);

//       for(i=0; i<N_Surf_S; i++)
//         cout<< " Orig scale " <<  Surf_sol_Exact[i]/Surf_sol[i]  << endl;
        #endif

// *********************************************************************************************
// assembling rhs and matrices
// *********************************************************************************************
        // working array for rhs is Surf_B, initialize Surf_B
        memset(Surf_B, 0, N_Surf_S*SizeOfDouble);
        Daxpy(N_Surf_S, tau*TDatabase::TimeDB->THETA3, Surf_rhs_Old, Surf_B);

        #ifdef __MOVINGGRID__
        // scaling due the old sol and old rhs are in the previous timestep domain
        // scale is the determinant of the jacobian on 3D surface, i.e, 2D
        Dscal(N_Surf_S, scale, Surf_B);
        #endif

        N_FESpaces = 2;
        fesp[0] = ScalarSpace;  // scalar space in entire domain
        fesp[1] = GridSpace;  // velocity space

        fefct[0] = Grid_U1;  // u1 = Grid_U1
        fefct[1] = Grid_U2;  //u2 = Grid_U2
        fefct[2] = Grid_U3;  //u3 = Grid_U2

        N_FESpaces_low = 1;
        fesp_low[0] = Surf_ScalarSpace;
        fefct_low[0] = Surf_Scalar;

        N_Rhs =1;
        surfferhs[0] = Surf_ScalarSpace;
        SRHSs[0] = Surf_rhs;
        memset(SRHSs[0], 0, N_Surf_S*SizeOfDouble);

        N_SquareMatrices = 3;
        SQMATRICES_SURF[0] = MatricesSurf_Scalar_A;
        SQMATRICES_SURF[0]->Reset();
        SQMATRICES_SURF[1] = MatricesSurf_Scalar_M;
        SQMATRICES_SURF[1]->Reset();
        SQMATRICES_SURF[2] = MatricesSurf_Scalar_TDivU;
        SQMATRICES_SURF[2]->Reset();

        AssembleSurf2D(N_FESpaces, fesp, fefct, N_FESpaces_low, fesp_low,
                       N_SquareMatrices, SQMATRICES_SURF, N_Rhs, SRHSs, 
                       surfferhs, Cell_No, Joint_No, TDatabase::TimeDB->CURRENTTIME,
                       tau, scale);

//         add rhs from current sub time step to rhs array B
//         Daxpy(N_Surf_S, tau, Surf_rhs, Surf_B);
        Daxpy(N_Surf_S, tau*TDatabase::TimeDB->THETA4, Surf_rhs, Surf_B);

        // add the tangential divergence matrix to the rhs mat
        MatAdd(MatricesSurf_Scalar_A_Old, MatricesSurf_Scalar_TDivU_Old, 1.);

        MatAdd(MatricesSurf_Scalar_M_Old, MatricesSurf_Scalar_A_Old, -tau*TDatabase::TimeDB->THETA2);
        gamma = 0.;
        memset(Surf_defect, 0, N_Surf_S*SizeOfDouble);
        MatVectActive(MatricesSurf_Scalar_M_Old, Surf_sol, Surf_defect);

        #ifdef __MOVINGGRID__
        // scaling due the old sol and old rhs are in the previous timestep domain
        Daxpy(N_Surf_S, scale, Surf_defect, Surf_B);
        #else
        Daxpy(N_Surf_S, 1., Surf_defect, Surf_B);
        #endif

        memcpy(Surf_rhs_Old, Surf_rhs, N_Surf_S*SizeOfDouble);
        MatricesSurf_Scalar_A_Old->Reset();
        MatAdd(MatricesSurf_Scalar_A_Old, MatricesSurf_Scalar_A, 1.);
        MatricesSurf_Scalar_M_Old->Reset();
        MatAdd(MatricesSurf_Scalar_M_Old, MatricesSurf_Scalar_M, 1.);
        MatricesSurf_Scalar_TDivU_Old->Reset();
        MatAdd(MatricesSurf_Scalar_TDivU_Old, MatricesSurf_Scalar_TDivU, 1.);

       //=====================================================================
       // assembling of system matrix
       //=====================================================================
        // add the tangential divergence matrix to the lhs mat
        MatAdd(MatricesSurf_Scalar_A, MatricesSurf_Scalar_TDivU, 1.);
        MatAdd(MatricesSurf_Scalar_M, MatricesSurf_Scalar_A,  tau*TDatabase::TimeDB->THETA1);

        if(Initialize_ScalarSolver)
         {
          /** Construct SMP Direct Solver for the surface scalar **/
           // Sort_ForDirectSolver() should be called before this construction direct solver!
           // collect data of the global system

           Entries_Scalar = MatricesSurf_Scalar_M->GetEntries();
           N_Eqn = N_Surf_S;
           ScalarKCol = sqSurfStructure->GetKCol();
           ScalarRowPtr = sqSurfStructure->GetRowPtr();
           N_Entries = ScalarRowPtr[N_Surf_S];
           KCol = new int[N_Entries];
           RowPtr = new int[N_Eqn+1];

           RowPtr[0] = 0 + 1; // fortran style for ParDiso direct solver
           Eq_pos = 0;
           for(i=0;i<N_Surf_S;i++)
            {
             begin = ScalarRowPtr[i];
             end = ScalarRowPtr[i+1];
             for(j=begin;j<end;j++)
              {
               KCol[Eq_pos] = ScalarKCol[j] + 1; // fortran style
               Eq_pos++;
              }

            RowPtr[i+1] = Eq_pos + 1; // fortran style for ParDiso direct solver
           }

          Scalar_SMPSolver = new TParDirectSolver(N_Eqn, RowPtr, KCol, N_Entries);
          Scalar_SMPSolver->SymbolicFactorize(Entries_Scalar);

          Initialize_ScalarSolver = FALSE;
          delete [] KCol; delete [] RowPtr;
         }

        Scalar_SMPSolver->FactorizeAndSolve(MatricesSurf_Scalar_M, Surf_B, Surf_sol);
//         Scalar_SMPSolver->Solve(MatricesSurf_Scalar_M, Surf_B, Surf_sol);

//         DirectSolver(MatricesSurf_Scalar_M, Surf_B, Surf_sol);

      } //  for(l=0;l<N_SubSteps;l++) 
   } // if (time_discs==2)

   if(m==1 || (m % TDatabase::TimeDB->STEPS_PER_IMAGE) == 0)
   if(TDatabase::ParamDB->WRITE_VTK)
    {
       MapSurfToDomain(Surf_Scalar, Scalar, Cell_No, Joint_No);

       os.seekp(std::ios::beg);
       if(img<10) os << VtkBaseName<<".0000"<<img<<".vtk" << ends;
       else if(img<100) os << VtkBaseName<<".000"<<img<<".vtk" << ends;
       else if(img<1000) os << VtkBaseName<<".00"<<img<<".vtk" << ends;
       else if(img<10000) os << VtkBaseName<<".0"<<img<<".vtk" << ends;
       else  os << VtkBaseName<<"."<<img<<".vtk" << ends;
       Output->WriteVtk(os.str().c_str());

       img++;
    }

    if(TDatabase::ParamDB->MEASURE_ERRORS)
     {
      MapSurfToDomain(Surf_Scalar, Scalar, Cell_No, Joint_No);

      GetSurfErrors(Scalar, Surf_Scalar, errors, Cell_No, Joint_No, ExactSall);

      OutPut( "L2: " << errors[0] << endl);
      OutPut( "H1-semi: " << errors[1] << endl);

      // error in L^infty(0,t,L^2)
      if (sqrt(errors[0]*errors[0]) > l_infty_l_2)
       {
        l_infty_l_2 = sqrt(errors[0]*errors[0]);
        l_infty_l_2_time =  TDatabase::TimeDB->CURRENTTIME;
       }
      OutPut( l_infty_l_2_time <<  " l_infty(L2(u)) " << l_infty_l_2 << endl);

      // error in L^2(0,t,L^2)
      l_2_l_2u += (errors[0]*errors[0] + olderror_l_2_l_2u)*TDatabase::TimeDB->TIMESTEPLENGTH/2.0;
      l_2_H_1u += (errors[1]*errors[1] + olderror_l_2_H_1u)*TDatabase::TimeDB->TIMESTEPLENGTH/2.0;

      OutPut( TDatabase::TimeDB->CURRENTTIME << " " );
      OutPut( "L2(0,t,L2)(u) : " <<  sqrt(l_2_l_2u) << endl);
      OutPut( TDatabase::TimeDB->CURRENTTIME << " " );
      OutPut( "L2(0,t,H1-semi)(u) : " <<  sqrt(l_2_H_1u) << endl);

      olderror_l_2_l_2u = errors[0]*errors[0];
      olderror_l_2_H_1u = errors[1]*errors[1];
     }

  } // while

// ======================================================================
// end of time cycle
// ======================================================================


//       OutPut( "T_inftyL2: " << l_infty_l_2 << endl);
//       OutPut( "T_L2L2: " << l_2_l_2u << endl);

//       OutPut( "T_inftyH1: " << T_inftyH1 << endl);
//       OutPut( "T_L2H1: " << sqrt(Terrors[1]) << endl);
//   OutPut("TIMESTEPLENGTH : " <<hmax*hmax << endl);
  CloseFiles();
  OutPut("used time: " << GetTime() << endl);
  return 0;
}


void smeshcreat()
 {

int i, j, k, l, m, n, *p;
double dt, ds, t, s, R, r;
std::ostringstream os;
os << " ";

m=35;
n =50;
R = 1.;
r = 0.25;


 p = new int[m*n];

ds=2.*Pi/double(m);
dt=2.*Pi/double(n);

cout<< " number of vertices " <<m*n <<endl;
os.seekp(std::ios::beg);
os << "torus.smesh" << ends;
std::ofstream dat(os.str().c_str());
if (!dat)
  {
   cerr << "cannot open file for output" << endl;
   exit(0);
  }
 dat << "# Surface mesh created for torus by MooNMD" << endl;
 dat << m*n << " " <<  3<< endl;

k=0;
for(i=0;i<m;i++)
 {
  s = double(i)*ds;
  for(j=0;j<n;j++)
   {
     t = double(j)*dt;
     p[i*n + j] = k;
     k++;
    dat<< 3 << "  " <<  "  " << (R + r*cos(s))* cos(t)<< " " << "  " <<  (R + r*cos(s))* sin(t) << " " << "  " <<   r*sin(s) <<endl;
   }
 }

dat << m*n << " " <<  0<< endl;
for(i=0;i<m;i++)
 {
  for(j=0;j<n;j++)
   {
    if(i!=m-1 && j!=n-1)
        dat<< 4 << "  " <<  "  " << p[i*n + j] << " " << "  " <<  p[i*n+(j+1)] << " " << "  " 
                                 <<  p[(i+1)*n + (j+1)]<< " " << "  " <<  p[(i+1)*n + j    ] <<endl;
    else if(i==m-1 && j==n-1)
        dat<< 4 << "  " <<  "  " << p[i*n + j] << " " << "  " <<  p[i*n+0] << " " << "  " 
                                  <<  p[(0)*n + (0)]<< " " << "  " <<  p[(0)*n + j    ] <<endl;
    else if(i==m-1)
        dat<< 4 << "  " <<  "  " << p[i*n + j] << " " << "  " <<  p[i*n+(j+1)] << " " << "  " 
                                 <<  p[(0)*n + (j+1)]<< " " << "  " <<  p[(0)*n + j    ] <<endl;
    else if(j==n-1)
        dat<< 4 << "  " <<  "  " << p[i*n + j] << " " << "  " <<  p[i*n+(0)] << " " << "  " 
                                 <<  p[(i+1)*n + (0)]<< " " << "  " <<  p[(i+1)*n + j    ] <<endl;
   }
 }
      dat <<  0<< endl;
      dat <<  0<< endl;
      dat.close();
//       cout << endl;
    cout << "Boundary wrote output into file " << endl;
delete [] p;

exit(0);

}


void  SmeshCreatQuad(TDomain *Domain, TDomain *SurfDomain, int **Cell_No, int **Joint_No, int level)
{
  int i,j,k,l,m,n,m1, N,N_G,N_SurfCells, maxlen, N_BoundVert, *EdgeList, *PointNeighb;
  int N_Cells, N_Joints, *EdgeListALL, *SCells,ID,N_Edges, begin,end, n1, *Bd_Part;
  int maxEpV, a, b, Neib[2], CurrNeib, len1, len2, Neighb_tmp;
  TBaseCell *Me;
  TJoint *joint;
//   TIsoJointEqN *isojoint;
  TCollection *Coll;
  TVertex **BoundVetrex;
  TBaseCell  **SurfCellTree;
  const int *TmpFV,  *TmpLen;
  boolean check=FALSE; 
//   TBoundComp *BoundComp;
  std::ostringstream os;
      os << " ";
  double x0, y0, z0;


  Coll = Domain->GetCollection(It_Finest, 0);
  N_SurfCells = 0;
  N_Cells = Coll->GetN_Cells();
  for(i=0;i<N_Cells;i++)
  {
    Me = Coll->GetCell(i);
    k = Me->GetN_Faces();
    for(l=0;l<k;l++)
     {
       joint = Me->GetJoint(l);
       if(joint->GetType() == BoundaryFace || joint->GetType() == IsoBoundFace ||
           joint->GetType() == IsoInterfaceJoint3D )
         N_SurfCells++;
     } // endfor l
  } // endfor i

  OutPut("N_SurfCells: " << N_SurfCells << endl);


  Cell_No[level] = new int[N_SurfCells];
  Joint_No[level] = new int[N_SurfCells];
  Bd_Part = new int[N_SurfCells];
  SCells = new int[4*N_SurfCells];
//   EdgeListALL = new int[6*N_SurfCells];
  BoundVetrex = new TVertex*[4*N_SurfCells];

  m1 = 0;
  N    = 0; 
  N_BoundVert = 0;
  for(i=0;i<N_Cells;i++)
   {
     Me = Coll->GetCell(i);
     k = Me->GetN_Faces();
     for(l=0;l<k;l++)
      {
        joint = Me->GetJoint(l);
        if(joint->GetType() == BoundaryFace || joint->GetType() == IsoBoundFace ||
           joint->GetType() == IsoInterfaceJoint3D )
         {
//            Cell_No[level][N] = i;
//            Joint_No[level][N] = l;
           Bd_Part[N++] =(((TBoundFace *)joint)->GetBoundComp())->GetID();
           Me->GetShapeDesc()->GetFaceVertex(TmpFV, TmpLen, maxlen);
 // getting freesurf vertices ---begin
           for(n=0;n<TmpLen[l];n++)
            {
              ID = 0;
              for(m=0;m<N_BoundVert;m++)
                if(BoundVetrex[m]==Me->GetVertex(TmpFV[maxlen*l+n]))
                 {
                   ID = 1;
                   SCells[m1] = m;
                   m1++;
                   break; // vertex is already in the collection
                 }

              if(ID==0)
               {
                 SCells[m1] = N_BoundVert;
                  m1++;
                 BoundVetrex[N_BoundVert++] = Me->GetVertex(TmpFV[maxlen*l+n]);
               }
          }
// getting freesurf vertices ---end
       }  //  if(joint->GetType() ==
     } // endfor l
   } // endfor i

   cout << " N_vertices:  " << N_BoundVert << endl;
//   for(i=0;i<N_SurfCells;i++)
//   cout << " SCells:  " <<SCells[3*i] << "    " << SCells[3*i+1] << "    " << SCells[3*i+2] << endl;
//   for(i=0;i<N_SurfCells;i++)
//  cout << " Cell_No:  " <<Cell_No[level][i] << " Joint_No  " << Joint_No[level][i] << endl;
// write .smesh for tetgen
      os.seekp(std::ios::beg);

      os << "Boundary.smesh" << ends;
      std::ofstream dat(os.str().c_str());
      if (!dat)
       {
        cerr << "cannot open file for output" << endl;
//         exit(0);
       }
       dat << "# Surface mesh created for droplet by MooNMD" << endl;

      dat << N_BoundVert << " " <<  3<< endl;
      for (i=0;i<N_BoundVert;i++)
       {
        BoundVetrex[i]->GetCoords(x0, y0, z0);
        dat<< i << "  " <<  "  " << x0 << " " << "  " <<  y0 << " " << "  " <<  z0 <<endl;
        }

      dat << N_SurfCells << " " <<  0<< endl;
      for (i=0;i<N_SurfCells;i++)
       {
         dat<< 4 << "  " <<  "  " << SCells[4*i    ] << " " << "  " <<  SCells[4*i + 1] << " " << "  " <<  SCells[4*i + 2]<< " " << "  " <<  SCells[4*i + 3] <<endl;
        }

      dat <<  0<< endl;
      dat <<  0<< endl;
      dat.close();
//       cout << endl;
    cout << "Boundary wrote output into file " << endl;


  delete []  SCells;
//   delete []  EdgeListALL;

exit(0);
}



void GetSphericalCoords(TFEVectFunct3D *CurrentGridPos, TFEVectFunct3D *AuxGridPos)
{
 TCollection *Coll;
 TFESpace3D *FESpace;
 TBaseCell *cell;
 FE3D FEId;
 TFE3D *Element;

 int i, j, k, m, n, N_U, N_Cells, N_Faces, N_LocalDOFs;
 int N_Inner, N_BoundaryNodes;
 int *BeginIndex, *GlobalNumbers, *DOF;

 double x, y, z, r;
 double *ValuesX, *ValuesY, *ValuesZ;
 double *AuxX, *AuxY;

 FESpace = CurrentGridPos->GetFESpace3D(); // assume fespacces in all vecfunct is same
 Coll = FESpace->GetCollection();
 N_Cells = Coll->GetN_Cells();

 N_U = CurrentGridPos->GetLength();
 BeginIndex = FESpace->GetBeginIndex();
 GlobalNumbers = FESpace->GetGlobalNumbers();

 ValuesX = CurrentGridPos->GetValues();
 ValuesY = ValuesX + N_U;
 ValuesZ = ValuesX + 2*N_U;

 AuxX = AuxGridPos->GetValues();
 AuxY = AuxX + N_U;

   // determine new position of boundary vertices
  for(i=0;i<N_Cells;i++)
   {
    cell  = Coll->GetCell(i);
    N_Faces = cell->GetN_Faces();
    for(j=0;j<N_Faces;j++)
     {
      if( (cell->GetJoint(j))->GetType() == IsoBoundFace ||
          (cell->GetJoint(j))->GetType() == BoundaryFace   )
       {
        DOF = GlobalNumbers + BeginIndex[i];

        FEId = FESpace->GetFE3D(i, cell);
        Element = TFEDatabase3D::GetFE3D(FEId);
        N_LocalDOFs = Element->GetN_DOF();

        for(k=0;k<N_LocalDOFs;k++)
         {
          m = DOF[k];;
          if(m>=N_Inner)
           {
            x = ValuesX[m]; // should be created before calling this routine
            y = ValuesY[m];
            z = ValuesZ[m];

            r = sqrt(x*x + y*y + z*z);
            AuxX[m] = atan2(y, x);
            AuxY[m] = acos(z/r);
          } //  if(m>=N_Inner)
         } //    for(k=0;k<N_LocalDOFs;k++)
       }   //   if( !(cell->GetJoint(j)->InnerJoint()) )
      } // endfor j
     } //   for(i=0;i<N_Cells;i++)



}







