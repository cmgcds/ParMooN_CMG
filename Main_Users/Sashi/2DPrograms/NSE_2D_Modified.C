// =======================================================================
//
// Purpose:     main program
//
// Author:      Volker Behns  22.07.97
//              Volker John   Jan. 2000
//
// =======================================================================

#include <Domain.h>
#include <Database.h>
#include <DiscreteForm2D.h>
#include <FEDatabase2D.h>
#include <FESpace2D.h>
#include <SquareStructure2D.h>
#include <Structure2D.h>
#include <AuxParam2D.h>
#include <QuadAffin.h>
#include <Solver.h>
#include <DirectSolver.h>
#include <Assemble2D.h>
#include <Output2D.h>
#include <LinAlg.h>
#include <NSE2D_ParamRout.h>
#include <TNSE2D_ParamRout.h>
#include <string.h>
#include <sstream>
#include <MooNMD_Io.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>

double bound = 0;

// ======================================================================
// utilities for main program
// ======================================================================
#include <MainUtilities.h>
#include <Upwind.h>

#include <NSE_MultiGrid.h>
#include <NSE_MGLevel1.h>
#include <NSE_MGLevel2.h>
#include <NSE_MGLevel3.h>
#include <NSE_MGLevel4.h>
#include <ItMethod.h>
#include <FixedPointIte.h>
#include <MultiGridIte.h>
#include <FgmresIte.h>
#include <MultiGrid2D.h>
#include <MGLevel2D.h>
#include <MultiGridScaIte.h>
#include <InterfaceJoint.h>
#include <TriaIsoparametric.h>
#include <QuadIsoparametric.h>
#include <QuadBilinear.h>
#include <Remesh2D.h>

#define AMG 0
#define GMG 1

// =======================================================================
// include current example
// =======================================================================
// #include "Examples/NSE_2D/Linear.h"
// #include "Examples/NSE_2D/Const.h"
// #include "Examples/NSE_2D/FSHabil.h"
// #include "Examples/NSE_2D/FSHabil_slip.h"
// #include "Examples/NSE_2D/DC2.h"
// #include "Examples/NSE_2D/DrivenCavity.h"
// #include "../Examples/NSE_2D/Benchmark.h"
// #include "Examples/NSE_2D/Benchmark_Neum.h"
// #include "Examples/NSE_2D/Frequence.h"
// #include "Examples/NSE_2D/Poiseuille.h"
// #include "Examples/NSE_2D/Poiseuille2.h"
// #include "Examples/NSE_2D/Einfach.h"
// #include "Examples/NSE_2D/SinCos.h"
// #include "Examples/NSE_2D/ChannelStep.h"
// #include "Examples/NSE_2D/BraessSarazin.h"
// #include "Examples/NSE_2D/BraessSarazinGeneralized.h"
// #include "Examples/NSE_2D/ZeroSolution.h"
// #include "Examples/NSE_2D/FerroMagnet.h"
// #include "../Examples/TNSE_2D/ChannelStep.h"
// #include "../Examples/NSE_2D/PressureGradient.h"
// #include "../Examples/NSE_2D/constant_velo_free_slip.h"
// #include "Examples/TNSE_2D/TrapezCavity.h"
// #include "Examples/TNSE_2D/StepCavity.h"
// #include "../Examples/TNSE_2D/ChannelStep.h"
// #include "Examples/NSE_2D/BenchmarkQuadNeum.h"
// #include "Examples/NSE_2D/atife_diri.h"
// #include "Examples/NSE_2D/atife_slip.h"
// #include "Examples/NSE_2D/Channel30_diri.h"
// #include "Examples/NSE_2D/Channel30_slip.h"
// #include "Examples/NSE_2D/atife_diri_sincos.h"
// #include "Examples/NSE_2D/atife_slip_sincos.h"
// #include "Examples/NSE_2D/Brennstoffzelle.h"
// #include "BOSCH/data/bosch_0433175329_00.h"
// #include "Examples/NSE_2D/NonSymmUnitSq.h"
// #include "Examples/NSE_2D/ZeroVeloLinPress.h"
// #include "Examples/NSE_2D/GleitRB_HalfCircle.h" 
// #include "Examples/NSE_2D/Example_sashikumar.h"
// #include "Examples/NSE_2D/Example_sashikumar.01.h"
// #include "Examples/NSE_2D/Example_sashikumar.02.h"
// #include "Examples/NSE_2D/Example_sashikumar.03.h"
// #include "Examples/NSE_2D/NoFlow.h"
// #include "Examples/TNSE_2D/UnitSquareHoleDrag.h"
// #include "Examples/TNSE_2D/SFBCavity.h"
// #include "Examples/TNSE_2D/SFBCavity_2.h"
// #include "Examples/TNSE_2D/SFBCavity_3.h"
#include "Examples/NSE_2D/Jump.h"
// using exact curvature on non-adapted mesh
void CalculateIntegrals(TFESpace2D *velocity_space, double *f1, double *f2)
{
  int i,j,k,l,m;
  TCollection *coll;
  TBaseCell *cell;
  int N_Cells, N_Joints;
  double x1, y1, x2, y2, dx, dy, x, y, xi, eta;
  double xmin, xmax, ymin, ymax;
  const double r=1.0; // circle radius
  const double r2=r*r; // square of circle radius
  double s1, s2;
  int N_Points;
  double X[4], Y[4], T[4], t;
  FE2D FEid;
  TBaseFunct2D *bf;
  TQuadFormula1D *qf;
  int N_QFPoints;
  double *zeta, *weights;
  int *BeginIndex, *GlobalNumbers, *DOF;
  RefTrans2D RFid;
  double values[MaxN_BaseFunctions2D];
  double v1[MaxN_BaseFunctions2D];
  double v2[MaxN_BaseFunctions2D];
  int N_BF;

  coll = velocity_space->GetCollection();
  N_Cells = coll->GetN_Cells();

  double We = TDatabase::ParamDB->WB_NR;
  BeginIndex = velocity_space->GetBeginIndex();
  GlobalNumbers = velocity_space->GetGlobalNumbers();

  qf = TFEDatabase2D::GetQuadFormula1D(Gauss11Line);
  qf->GetFormulaData(N_QFPoints, weights, zeta);

  for(i=0;i<N_Cells;i++)
  {
    N_Points = 0;
    cell = coll->GetCell(i);
    N_Joints = cell->GetN_Joints();
    for(j=0;j<N_Joints;j++)
    {
      cell->GetVertex(j)->GetCoords(x1, y1);
      cell->GetVertex((j+1)%N_Joints)->GetCoords(x2, y2);
      dx = x2-x1;
      dy = y2-y1;
      xmin = (x1<x2)?x1:x2;
      xmax = (x1>x2)?x1:x2;
      ymin = (y1<y2)?y1:y2;
      ymax = (y1>y2)?y1:y2;

      if(fabs(dy)<1e-8)
      {
        // horizontal joint
        if(y1*y1<=r2)
        {
          s1 = sqrt(r2-y1*y1);
          s2 = -s1;
          if(xmin<=s1 && s1<=xmax)
          {
            X[N_Points] = s1;
            Y[N_Points] = y1;
            t = atan2(Y[N_Points], X[N_Points]);
            if(t<0) t+=2*Pi;
            T[N_Points] = t;
            // cout << "cell: " << i << " joint: " << j << endl;
            // cout << X[N_Points] << " " << Y[N_Points] << " " << T[N_Points]
            // << endl;
            N_Points++;
          }
          if(xmin<=s2 && s2<=xmax)
          {
            X[N_Points] = s2;
            Y[N_Points] = y1;
            t = atan2(Y[N_Points], X[N_Points]);
            if(t<0) t+=2*Pi;
            T[N_Points] = t;
            // cout << "cell: " << i << " joint: " << j << endl;
            // cout << X[N_Points] << " " << Y[N_Points] << " " << T[N_Points]
            // << endl;
            N_Points++;
          }
        }
      }
      else
      {
        if(fabs(dx)<1e-8)
        {
          // vertical joint
          if(x1*x1<=r2)
          {
            s1 = sqrt(r2-x1*x1);
            s2 = -s1;
            if(ymin<=s1 && s1<=ymax)
            {
              X[N_Points] = x1;
              Y[N_Points] = s1;
              t = atan2(Y[N_Points], X[N_Points]);
              if(t<0) t+=2*Pi;
              T[N_Points] = t;
              // cout << "cell: " << i << " joint: " << j << endl;
              // cout << X[N_Points] << " " << Y[N_Points] << " " << T[N_Points]
              // << endl;
              N_Points++;
            }
            if(ymin<=s2 && s2<=ymax)
            {
              X[N_Points] = x1;
              Y[N_Points] = s2;
              t = atan2(Y[N_Points], X[N_Points]);
              if(t<0) t+=2*Pi;
              T[N_Points] = t;
              // cout << "cell: " << i << " joint: " << j << endl;
              // cout << X[N_Points] << " " << Y[N_Points] << " " << T[N_Points]
              // << endl;
              N_Points++;
            }
          }
        }
        else
        {
          // diagonal joint
          cout << "Diagonal joints have not been implemented yet! " << endl;
          exit(-1);
        }
      }
    } // endfor j

    if(N_Points == 2)
    {
      // if(N_Points) cout << "cell: " << i << endl;

      // correction around 0 angle
      if(fabs(T[0]-T[1])>Pi)
      {
        if(T[0]>T[1])
          T[0]-=2*Pi;
        else
          T[1]-=2*Pi;
      }

      // sort such that T is increasing
      if(T[1]<T[0])
      {
        t = T[0];    x = X[0];    y = Y[0];
        T[0] = T[1]; X[0] = X[1]; Y[0] = Y[1];
        T[1] = t;    X[1] = x;    Y[1] = y;
      }
      // for(j=0;j<N_Points;j++)
      // {
      //   cout << "X: " << X[j] << " Y: " << Y[j];
      //   cout << " T: " << T[j]*180/Pi << endl;
      // }

      // find basis functions on cell i
      FEid = velocity_space->GetFE2D(i, cell);
      bf = TFEDatabase2D::GetBaseFunct2DFromFE2D(FEid);
      N_BF = bf->GetDimension();
      for(j=0;j<N_BF;j++)
      {
        v1[j] = 0;
        v2[j] = 0;
      }

      DOF = GlobalNumbers + BeginIndex[i];

      switch(cell->GetType())
      {
        case Triangle:
          TFEDatabase2D::SetCellForRefTrans(cell, TriaAffin);
          RFid = TriaAffin;
        break;

        case Quadrangle:
          TFEDatabase2D::SetCellForRefTrans(cell, QuadBilinear);
          RFid = QuadBilinear;
        break;

        case Parallelogram:
        case Rectangle:
          TFEDatabase2D::SetCellForRefTrans(cell, QuadAffin);
          RFid = QuadAffin;
        break;

        default:
          cout << "Illegal cell shape" << endl;
          exit(-1);
      }

      for(j=0;j<N_QFPoints;j++)
      {
        t = T[0]+(T[1]-T[0])*(zeta[j]+1)/2;
        x = r*cos(t);
        y = r*sin(t);
        TFEDatabase2D::GetRefFromOrig(RFid, x, y, xi, eta);

        bf->GetDerivatives(D00, xi, eta, values);

        for(k=0;k<N_BF;k++)
        {
          v1[k] += (weights[j]*(T[1]-T[0])*0.5*values[k]*x)/We;
          v2[k] += (weights[j]*(T[1]-T[0])*0.5*values[k]*y)/We;
        }
      } // endfor j

      for(k=0;k<N_BF;k++)
      {
        f1[DOF[k]] -= v1[k];
        f2[DOF[k]] -= v2[k];
      }
    } 
    else
      if(N_Points != 0)
      {
        cout << "Only two intersection points can be handled" << endl;
        exit(-1);
      }
  } // endfor i
} // CalculateIntegrals


// using CSF for curvature on non-adapted mesh
void CSF_Int(TFESpace2D *velocity_space, double *f1, double *f2, double h_max)
{
  int i,j,k,l,m;
  TCollection *coll;
  TBaseCell *cell;
  int N_Cells, N_Joints;
  double x1, y1, x2, y2, dx, dy, x, y, *xi, *eta;
  double xmin, xmax, ymin, ymax;
  const double r=1.0; // circle radius
  const double r2=r*r; // square of circle radius
  double s1, s2;
  int N_Points;
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D], T[4], t;
  FE2D FEId;
  TFE2D *ele;
  TBaseFunct2D *bf;
  QuadFormula2D QuadFormula;
  TQuadFormula2D *qf2;
  int N_QFPoints;
  double *zeta, *weights, Mult, limit;
  int *BeginIndex, *GlobalNumbers, *DOF;
  RefTrans2D RefTrans;
  BF2DRefElements RefElement;
  double values[MaxN_QuadPoints_2D][MaxN_BaseFunctions2D];
  double v1[MaxN_BaseFunctions2D];
  double v2[MaxN_BaseFunctions2D];
  int N_Vertices;
  boolean  CSF_CELL;
  double AbsDetjk[MaxN_QuadPoints_2D], delta, radius, local_r, dist, n1, n2, theta, dirac_delta;
  TRefTrans2D *F_K;
  int N_BaseFunct, *N_BaseFuncts;

// transition region length
  limit = 1.5*h_max;
//   delta = 2.*limit;


// cout << " tehta " << (180./Pi)*atan2(-1.,0.) <<endl;
// exit(0);

  coll = velocity_space->GetCollection();
  N_Cells = coll->GetN_Cells();

  double We = TDatabase::ParamDB->WB_NR;
  BeginIndex = velocity_space->GetBeginIndex();
  GlobalNumbers = velocity_space->GetGlobalNumbers();
  N_BaseFuncts = TFEDatabase2D::GetN_BaseFunctFromFE2D();

  for(i=0;i<N_Cells;i++)
  {
    cell = coll->GetCell(i);
    N_Vertices = cell->GetN_Vertices();
    CSF_CELL = FALSE;
    for(k=0;k<N_Vertices;k++)
     {
      cell->GetVertex(k)->GetCoords(x, y);
      radius = sqrt(x*x+y*y);
      if( (radius >(r - limit)) && (radius< (r + limit))  )
       {
        CSF_CELL = TRUE;
        break;
       }
     } // endfor k

    if(CSF_CELL)
     {
      
      cell->SetPhase_No(5);
//  part of the cell i is in CSF region
      //cout << "Cell " << i << " has free surface." << endl;
      FEId = velocity_space->GetFE2D(i, cell);
      bf = TFEDatabase2D::GetBaseFunct2DFromFE2D(FEId);
      RefElement = TFEDatabase2D::GetRefElementFromFE2D(FEId);
      l = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);

      switch(RefElement)
       {
        case BFUnitSquare:

          QuadFormula = TFEDatabase2D::GetQFQuadFromDegree(5*l);
          qf2 = TFEDatabase2D::GetQuadFormula2D(QuadFormula);
          qf2->GetFormulaData(N_Points, weights, xi, eta);
          RefTrans = QuadBilinear;
          F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
          ((TQuadBilinear *)F_K)->SetCell(cell);
          ((TQuadBilinear *)F_K)->GetOrigFromRef(N_Points, xi, eta,
                                         X, Y, AbsDetjk);

         break;
        default:
          cout << " CSF not implemented for Trangular element" << endl;
       } // endswitch
 

     DOF = GlobalNumbers + BeginIndex[i];
     N_BaseFunct = N_BaseFuncts[FEId];

     for(k=0;k<N_Points;k++)
      {
        bf->GetDerivatives(D00, xi[k], eta[k], values[k]);
      }
     for(k=0;k<N_Points;k++)
      {
      x = X[k];
      y = Y[k];
//       cout <<' ' << " x " << x << ' ' << " y " << y <<endl;
      local_r = sqrt(x*x* + y*y);
     if((local_r > (r - limit)) && (local_r < (r + limit)))
      {

//    cout <<' ' << " x " << x << ' ' << " y " << y <<' ' << local_r << endl;
       dist = -(r - local_r);
       theta = atan2(y,x);
       if(theta<0.) theta +=2.*Pi;

// cout<< "theta : " << (180./Pi)*theta << endl;
       n1 = r*cos(theta);
       n2 = r*sin(theta);
//    dist is a discrete zero level-set function
//    dist is negative inside and positive outside varies from -limit to +limit
//    make  dist as -1 to +1 in the transition region   
       dirac_delta = (1./(2.*limit))*(1.+ cos(Pi*dist/limit))/0.5;
       Mult = weights[k]*AbsDetjk[k];
 OutPut(" limit " <<  limit << " x local_r " <<  local_r << " dirac_delta " << dirac_delta<< " cos dist " <<cos(Pi*dist) <<endl);
//        for(l=0;l<N_BaseFunct;l++)
//         {
//          v1[l] = 0.;
//          v2[l] = 0.;
//         }

       for(l=0;l<N_BaseFunct;l++)
        {
         f1[DOF[l]] += Mult*values[k][l]*n1*dirac_delta;
         f2[DOF[l]] += Mult*values[k][l]*n2*dirac_delta; 
        } // for(l
//        for(l=0;l<N_BaseFunct;l++)
//         {
//          f1[DOF[k]] += v1[l];
//          f2[DOF[k]] += v2[l];
//         } // for(l

       } // if((local_r >
      }  //  for(k=0
     } //     if(CSF_CELL)
  } // endfor i
} // CalculateIntegrals


// using exact curvature on adapted mesh
void IsoIntegrals(TFESpace2D *velocity_space, double *f1, double *f2)
{
  int i,j,k,l, l2;
  int N_Cells, N_Joints, N_Vertices;
  TBaseCell *cell;
  TCollection *coll;
  int *BeginIndex, *GlobalNumbers, *DOF, *Neib_DOF;
  TJoint *joint;
  double xm, ym, x, y;
  const double r=1;
  const double r2=r*r;
  RefTrans2D RFid;
  TRefTrans2D *F_K;
  double ta1, ta2;
  FE2D FEid;
  TBaseFunct2D *bf;
  TQuadFormula1D *qf;
  int N_QFPoints;
  double *zeta, *weights;
  double *values[1];
  double val[MaxN_BaseFunctions2D];
  double v1[MaxN_BaseFunctions2D];
  double v2[MaxN_BaseFunctions2D];
  int N_BF;
  double t, s, len, n1, n2, Kappa;
  double We = TDatabase::ParamDB->WB_NR;
  TInterfaceJoint *intface_joint;
  TBaseCell *Neib_cell, *Me;
  int Neib_Cell_ID, j1, i1;
  TJoint *Neib_joint;
  boolean update;

  Kappa = 1.0; // curvature
  values[0] = val;

  coll = velocity_space->GetCollection();
  N_Cells = coll->GetN_Cells();

  BeginIndex = velocity_space->GetBeginIndex();
  GlobalNumbers = velocity_space->GetGlobalNumbers();

  qf = TFEDatabase2D::GetQuadFormula1D(Gauss11Line);
  qf->GetFormulaData(N_QFPoints, weights, zeta);



  for(i=0;i<N_Cells;i++)
  {
    cell = coll->GetCell(i);

    N_Joints = cell->GetN_Joints();
    for(j=0;j<N_Joints;j++)
    {
      joint = cell->GetJoint(j);
      if(joint->GetType() == InterfaceJoint)
      {
        xm = 0;
        ym = 0;
        N_Vertices = cell->GetN_Vertices();
        if(N_Vertices != 4)
        {
          cout << "ONLY quads are allowed" << endl;
          exit(-1);
        }
        for(k=0;k<N_Vertices;k++)
        {
          cell->GetVertex(k)->GetCoords(x, y);
          xm += x;
          ym += y;
        } // endfor k
        xm /= N_Vertices;
        ym /= N_Vertices;
        if(xm*xm+ym*ym<r2)
        {
          // cell is inside the circle
          // cout << "cell inside" << endl;
          if(TDatabase::ParamDB->USE_ISOPARAMETRIC)
            RFid = QuadIsoparametric;
          else
            RFid = QuadBilinear;
          F_K = TFEDatabase2D::GetRefTrans2D(RFid);
          TFEDatabase2D::SetCellForRefTrans(cell, RFid);

          // find basis functions on cell i
          FEid = velocity_space->GetFE2D(i, cell);
          bf = TFEDatabase2D::GetBaseFunct2DFromFE2D(FEid);
          N_BF = bf->GetDimension();
          for(k=0;k<N_BF;k++)
          {
            v1[k] = 0;
            v2[k] = 0;
          }

          DOF = GlobalNumbers + BeginIndex[i];

          for(k=0;k<N_QFPoints;k++)
          {
            t = zeta[k];
            F_K->GetTangent(j, t, ta1, ta2);
            len = sqrt(ta1*ta1+ta2*ta2);
            n1 =  ta2/len;
            n2 = -ta1/len;
            bf->GetValues(1, zeta+k, j, values);
            // cout << "tangent: " << ta1 << " " << ta2 << endl;

            s = weights[k]*len/We;
            for(l=0;l<N_BF;l++)
            {
              v1[l] += val[l]*s*n1*Kappa;
              v2[l] += val[l]*s*n2*Kappa;
            }
          }
          for(k=0;k<N_BF;k++)
          {
            f1[DOF[k]] += v1[k];
            f2[DOF[k]] += v2[k];
          }

// // cout<< " test " <<endl;
// // for nonconforming element 
// //  normally entry values should be zero even for nonconforming
// 
//     intface_joint = (TInterfaceJoint *)joint;
//     Neib_cell = intface_joint->Get_Neighb(cell);
//     Neib_Cell_ID = Neib_cell->GetPhase_No();
// //          cout<< i<< " Neib_Cell_ID "<<Neib_Cell_ID<<endl;
// 
//      if(Neib_Cell_ID != 1)
//       {
//        cout<< " two cells are in same phase"<<endl;
//          exit(-1);
//       }
// // 
//     for(i1=0;i1<N_Cells;i1++)
//      {
//       Me = coll->GetCell(i1);
//       if(Me==Neib_cell)
//        {
//       for(j1=0;j1<N_Joints;j1++)
//        {
//         Neib_joint = Neib_cell->GetJoint(j1);
//         if(Neib_joint == joint)
//          {
// 
//            if(TDatabase::ParamDB->USE_ISOPARAMETRIC)
//             RFid = QuadIsoparametric;
//           else
//             RFid = QuadBilinear;
// 
//           F_K = TFEDatabase2D::GetRefTrans2D(RFid);
//           TFEDatabase2D::SetCellForRefTrans(Me, RFid);
// //                cout << " test cell ID " << Neib_Cell_ID << endl;
//           // find basis functions on cell i1
//           FEid = velocity_space->GetFE2D(i1, Me);
//           bf = TFEDatabase2D::GetBaseFunct2DFromFE2D(FEid);
//           N_BF = bf->GetDimension();
//           for(k=0;k<N_BF;k++)
//           {
//             v1[k] = 0.;
//             v2[k] = 0.;
//           }
// 
//           Neib_DOF = GlobalNumbers + BeginIndex[i1];
// 
//           for(k=0;k<N_QFPoints;k++)
//           {
//             t = zeta[k];
//             F_K->GetTangent(j1, t, ta1, ta2);
//             len = sqrt(ta1*ta1+ta2*ta2);
//             n1 =  ta2/len;
//             n2 = -ta1/len;
// // to make outward normal w.r.t inner phase
// //             n1 = -n1 ;
// //             n2 = -n2;
// 
//             bf->GetValues(1, zeta+k, j1, values);
//             // cout << "tangent: " << ta1 << " " << ta2 << endl;
// 
//             s = 0.5*weights[k]*len/We;
//             for(l=0;l<N_BF;l++)
//             {
//             update = TRUE;
//             for(l2=0;l2<N_BF;l2++)
//              {
//               if(Neib_DOF[l]==DOF[l2] )
//                {
// //              already updated w.r.t. inner phase
//                 update = FALSE;
//                 break;
//                }
//              }
// 
//             if(update)
//              {
//               v1[l] += val[l]*s*n1*Kappa;
//               v2[l] += val[l]*s*n2*Kappa;
// 
//              }
//             }
//           }
//           for(k=0;k<N_BF;k++)
//            {
//             f1[Neib_DOF[k]] += v1[k];
//             f2[Neib_DOF[k]] += v2[k];
// //             cout<< l<< " v1[k]: " << v1[k]<< " v2[k]: " << v2[k]<<endl;
//            }
// 
// 
//           } // if(Neib_joint == joint)
//          } //  for(j1
//         } //   if(Me==Neib_cell)
//       } //  for(i1=0;i1<N_Cells;i1++)
// 

        } // innercell
      } // endif
    } // endfor j
  } // endfor i
}


// using Laplace Beltrami operator
void IsoInt(TFESpace2D *velocity_space, double *rhs1, double *rhs2)
{
  int i,j,k,l,m, i1, j1, k1;
  TBaseCell *cell, *Me, *Neib_cell;
  TCollection *Coll;
  int N_Cells, N_Vertices, N_Edges;
  TJoint *joint, *joint_Neib;
  TInterfaceJoint *intface_joint;
  TIsoBoundEdge *isoboundedge;
  TBoundComp *BoundComp;
  int comp, Neib_Cell_ID, N_Joints_Neib;
  double t0, t1, n0, n1, normn;
  BoundCond Cond0, Cond1;
  int JointNumbers[MAXN_JOINTS], N_IsoJoints;
  FE2D FEId;
  TFE2D *ele;
  RefTrans2D RefTrans;
  TRefTrans2D *F_K;
  BF2DRefElements RefElement;
  QuadFormula1D LineQuadFormula;
  TQuadFormula1D *qf1;
  int N_LinePoints, update;
  double *LineWeights, *zeta;
  double x0, y0, x1, y1, xi, eta;
  int N_BaseFunct, N_BaseFunct_Neib, *N_BaseFuncts, N_BaseFuncts_Neib;
  double **uref, **uxiref, **uetaref;
  double uorig[MaxN_BaseFunctions2D], uxorig[MaxN_BaseFunctions2D];
  double uyorig[MaxN_BaseFunctions2D];
  BaseFunct2D *BaseFuncts;
  double r2, r;
  int *KCol, *RowPtr;
  double *ValuesA11, *ValuesA12, *ValuesA21, *ValuesA22;
  int *BeginIndex, *GlobalNumbers, *DOF, *DOF_Neib, TestDOF, AnsatzDOF;
  int index1, index2, Cell_ID;
  double val, theta;

   double We = TDatabase::ParamDB->WB_NR;
//   double We = 0.1;

  BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
  N_BaseFuncts = TFEDatabase2D::GetN_BaseFunctFromFE2D();

  Coll = velocity_space->GetCollection();
  N_Cells = Coll->GetN_Cells();

  BeginIndex = velocity_space->GetBeginIndex();
  GlobalNumbers = velocity_space->GetGlobalNumbers();

  for(i=0;i<N_Cells;i++)
  {
    // cout << endl << "CELL number: " << i << endl;
   cell = Coll->GetCell(i);
   Cell_ID = cell->GetPhase_No();
   if(Cell_ID==0)
    {
     N_Edges = cell->GetN_Edges();
     for(j=0;j<N_Edges;j++)
     {
      joint = cell->GetJoint(j);
      if(joint->GetType() == InterfaceJoint)
      {
     // cout << " test cell ID " << Cell_ID << endl;
     // cout << "Cell " << i << " has interface" << endl;
     // cout << "joint number: " << j << endl;
      FEId = velocity_space->GetFE2D(i, cell);
      ele = TFEDatabase2D::GetFE2D(FEId);
      RefElement = TFEDatabase2D::GetRefElementFromFE2D(FEId);

      switch(RefElement)
       {
        case BFUnitSquare:
	  if(TDatabase::ParamDB->USE_ISOPARAMETRIC)
           {
	    RefTrans = QuadIsoparametric;
            F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
            ((TQuadIsoparametric *)F_K)->SetCell(cell);
           }
	  else
	   {
	    RefTrans = QuadBilinear;
            F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
            ((TQuadBilinear *)F_K)->SetCell(cell);
//             TFEDatabase2D::SetCellForRefTrans(cell, RefTrans);
	   }
        break;

        case BFUnitTriangle:
          RefTrans = TriaIsoparametric;
          F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
          ((TTriaIsoparametric *)F_K)->SetCell(cell);
          //cout << "Cell RefElement" << RefElement  << endl;
        break;
       } // endswitch

      l = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
      LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(20);
      qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
      qf1->GetFormulaData(N_LinePoints, LineWeights, zeta);
      TFEDatabase2D::GetBaseFunct2DFromFE2D(FEId)
                  ->MakeRefElementData(LineQuadFormula);

      DOF = GlobalNumbers + BeginIndex[i];
      N_BaseFunct = N_BaseFuncts[FEId];

//       cell->GetVertex(j)->GetCoords(x0, y0);
//       cell->GetVertex((j+1) % N_Edges)->GetCoords(x1, y1);
//     //    cout<< " y0= " <<y0<<" y1= "<<y1<<"  x0= "<<x0<<"  x1= "<<x1<<endl;

      uref = TFEDatabase2D::GetJointValues2D(BaseFuncts[FEId],
                       LineQuadFormula, j);
      uxiref = TFEDatabase2D::GetJointDerivatives2D(BaseFuncts[FEId],
                       LineQuadFormula, j, D10);
      uetaref = TFEDatabase2D::GetJointDerivatives2D(BaseFuncts[FEId],
                       LineQuadFormula, j, D01);

      for(k=0;k<N_LinePoints;k++)
       {
        F_K->GetTangent(j, zeta[k], t0, t1);
        normn = sqrt(t0*t0+t1*t1);
        n0 =  t1/normn;
        n1 = -t0/normn;
        // cout << "zeta: " << zeta[k] << endl;
        // cout << "k= " << k << "  tangent: " << t0 << " " << t1 << endl;
        // cout << "length: " << sqrt(t0*t0+t1*t1) << endl;

        switch(RefElement)
         {
          case BFUnitSquare:
             if(TDatabase::ParamDB->USE_ISOPARAMETRIC)
               ((TQuadIsoparametric *)F_K)->GetOrigValues(j,  zeta[k],
                        N_BaseFunct, uref[k], uxiref[k], uetaref[k],
                        uorig, uxorig, uyorig);
             else
	      {
                switch(j)
                {
                 case 0:
                      xi  = zeta[k]; eta = -1;
                      break;
                 case 1:
                      xi = 1; eta = zeta[k];
                      break;
                 case 2:
                     xi  = -zeta[k]; eta = 1;
                     break;
                 case 3:
                    xi = -1 ; eta = -zeta[k];
                    break;
                }
                ((TQuadBilinear *)F_K)->GetOrigValues(xi, eta,
                        N_BaseFunct, uref[k], uxiref[k], uetaref[k],
                        uorig, uxorig, uyorig);

               }
          break;

          case BFUnitTriangle:
               ((TTriaIsoparametric *)F_K)->GetOrigValues(j,  zeta[k],
                        N_BaseFunct, uref[k], uxiref[k], uetaref[k],
                        uorig, uxorig, uyorig);
          break;
         } // endswitch

          // modify matrices
        r2 = 1./(t0*t0+t1*t1);
        r = sqrt(t0*t0+t1*t1)/We;
        for(l=0;l<N_BaseFunct;l++)
         {
          // updating rhs
          val = -r2* t0 * (uxorig[l]*t0+uyorig[l]*t1);
          val *= LineWeights[k]*r;
          rhs1[DOF[l]] += val;

          val = -r2* t1 * (uxorig[l]*t0+uyorig[l]*t1);
          val *= LineWeights[k]*r;
          rhs2[DOF[l]] += val;
         } // endfor l
        } // endfor k
 
// //  begin neighbour cell *****************************
//        intface_joint = (TInterfaceJoint *)joint;
//        Neib_cell = intface_joint->Get_Neighb(cell);
//        Neib_Cell_ID = Neib_cell->GetPhase_No();
//        //   cout<< i<< " Neib_Cell_ID "<<Neib_Cell_ID<<endl;
// 
//        if(Neib_Cell_ID != 1)
//         {
//       cout<< " two cells are in same phase"<<endl;
//          exit(-1);
//         }
// 
//        for(i1=0;i1<N_Cells;i1++)
//         {
//          Me = Coll->GetCell(i1);
//          if(Me==Neib_cell)
//           {
//            N_Joints_Neib = Me->GetN_Joints();
//            for(j1=0;j1<N_Joints_Neib;j1++)
//             {
//              joint_Neib = Me->GetJoint(j1);
//              if(joint_Neib == joint)
//               {
// //                cout << " test cell ID " << Neib_Cell_ID << endl;
//                FEId = velocity_space->GetFE2D(i1, Me);
//                ele = TFEDatabase2D::GetFE2D(FEId);
//                RefElement = TFEDatabase2D::GetRefElementFromFE2D(FEId);
// 
//             switch(RefElement)
//              {
//               case BFUnitSquare:
// 	         if(TDatabase::ParamDB->USE_ISOPARAMETRIC)
//                   {
// 		   RefTrans = QuadIsoparametric;
//                    F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
//                    ((TQuadIsoparametric *)F_K)->SetCell(Me);
//                   }
//                 else
//                   {
//                    RefTrans = QuadBilinear;
//                    F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
//                    ((TQuadBilinear *)F_K)->SetCell(Me);
//                   }
//               break;
// 
//               case BFUnitTriangle:
//                    RefTrans = TriaIsoparametric;
//                    F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
//                   ((TTriaIsoparametric *)F_K)->SetCell(Me);
//               break;
//              } // endswitch
// 
//       l = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
//       LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(2*l);
//       qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
//       qf1->GetFormulaData(N_LinePoints, LineWeights, zeta);
//       TFEDatabase2D::GetBaseFunct2DFromFE2D(FEId)
//                   ->MakeRefElementData(LineQuadFormula);
// 
//       DOF_Neib = GlobalNumbers + BeginIndex[i1];
// //       N_BaseFunct_Neib = N_BaseFuncts[FEId];
// 
//       uref = TFEDatabase2D::GetJointValues2D(BaseFuncts[FEId],
//                        LineQuadFormula, j1);
//       uxiref = TFEDatabase2D::GetJointDerivatives2D(BaseFuncts[FEId],
//                        LineQuadFormula, j1, D10);
//       uetaref = TFEDatabase2D::GetJointDerivatives2D(BaseFuncts[FEId],
//                        LineQuadFormula, j1, D01);
// 
//       for(k=0;k<N_LinePoints;k++)
//        {
//         F_K->GetTangent(j1, zeta[k], t0, t1);
//         normn = sqrt(t0*t0+t1*t1);
// //      tangent w.r.t inner domain
//         t0 = -t0;
//         t1 = -t1;
// //         n0 =  t1/normn;
// //         n1 = -t0/normn;
//         // cout << "zeta: " << zeta[k] << endl;
//         // cout << "k= " << k << "  tangent: " << t0 << " " << t1 << endl;
//         // cout << "length: " << sqrt(t0*t0+t1*t1) << endl;
// 
//         switch(RefElement)
//          {
//           case BFUnitSquare:
//              if(TDatabase::ParamDB->USE_ISOPARAMETRIC)
//               {
//                ((TQuadIsoparametric *)F_K)->GetOrigValues(j1,  zeta[k],
//                         N_BaseFunct, uref[k], uxiref[k], uetaref[k],
//                         uorig, uxorig, uyorig);
//               }
//              else
// 	      {
//                 switch(j1)
//                 {
//                  case 0:
//                       xi  = zeta[k]; eta = -1;
//                       break;
//                  case 1:
//                       xi = 1; eta = zeta[k];
//                       break;
//                  case 2:
//                      xi  = -zeta[k]; eta = 1;
//                      break;
//                  case 3:
//                     xi = -1 ; eta = -zeta[k];
//                     break;
//                 }
//                 ((TQuadBilinear *)F_K)->GetOrigValues(xi, eta,
//                         N_BaseFunct, uref[k], uxiref[k], uetaref[k],
//                         uorig, uxorig, uyorig);
// 
//                }
//           break;
// 
//           case BFUnitTriangle:
//                ((TTriaIsoparametric *)F_K)->GetOrigValues(j1,  zeta[k],
//                         N_BaseFunct, uref[k], uxiref[k], uetaref[k],
//                         uorig, uxorig, uyorig);
//           break;
//          } // endswitch
// 
//           // modify matrices
//         r2 = 1./(t0*t0+t1*t1);
//         r = sqrt(t0*t0+t1*t1)/We;
//         for(l=0;l<N_BaseFunct;l++)
//          {
//            update = 1;
//            for(k1=0;k1<N_BaseFunct;k1++)
//             if(DOF_Neib[l] == DOF[k1])
//              {
// //              already updated from the inner doamin cell
//               update = 0;
//               break;
//              }
//          if(update)
// 	 {
//           // updating rhs
//           val = -r2* t0 * (uxorig[l]*t0+uyorig[l]*t1);
//           val *= LineWeights[k]*r;
//           rhs1[DOF_Neib[l]] += val;
// //           cout << " rhs1: " << val;
// 
//           val = -r2* t1 * (uxorig[l]*t0+uyorig[l]*t1);
//           val *= LineWeights[k]*r;
//           rhs2[DOF_Neib[l]] += val;
// //           cout <<" rhs2: " << val<<endl;
// 	  }
//          } // endfor l
//         } // endfor k
// 
//        } // if joint
//       }  // for j1
//      } // if Me
//     } // for i1
// //  end neighbour cell *****************************

     } // if interfacejoint
    } // for j
   } // if cellID==0
  } // endfor i
 }

void IsoInt_spline(int N_E, TBaseCell **cell, int *EdgeNo,
                   TFESpace2D *velocity_space, double *f1, double *f2)
{
  int i,j,k,l, m, s1;
  int N_Cells, N_Joints, N_Vertices;
  TBaseCell *Cell;
  TCollection *coll;
  int *BeginIndex, *GlobalNumbers, *DOF;
  TJoint *joint;
  TInterfaceJoint *interface;
  TBoundComp2D *comp;
  double xm, ym, X, Y, T, *Local_T;
  const double r=1.;
  const double r2=r*r;
  RefTrans2D RFid;
  TRefTrans2D *F_K;
  double ta1, ta2, dx, dy, dxx, dyy, dt;
  FE2D FEid;
  TBaseFunct2D *bf;
  TQuadFormula1D *qf;
  int N_QFPoints;
  double *zeta, *weights, **cvr;
  double *values[1];
  double val[MaxN_BaseFunctions2D];
  double v1[MaxN_BaseFunctions2D];
  double v2[MaxN_BaseFunctions2D];
  int N_BF;
  double s, len, n1, n2, Kappa, t0, t1;

  int  ISpline, N_Splines, N_V, VSP, i3;
  double *h, *t;
  double *a, *b, *c, *x, *y, teta;
  double *rhs_Spline, *Mx, *My, *Params, *Param9;
  double phi1, phi2, phi3, phi4;
  double dx0, dy0, dx1, dy1, tx, ty;
  TIsoInterfaceJoint *isojoint;
  TVertex **IsoVertices;
  TJoint *Joint;
  TBaseCell *Me;
  double We = TDatabase::ParamDB->WB_NR;
  Kappa = -1.0; // curvature
  values[0] = val;

  if(TDatabase::ParamDB->USE_ISOPARAMETRIC)
   N_V = 2*N_E + 1;
  else
   N_V = N_E + 1; // star and end vertex same 2times

  N_Splines = N_V-1;
  h = new double[N_V+1];
  t = new double[N_V+1];
  a = new double[N_V+1];
  b = new double[N_V+1];
  c = new double[N_V+1];
  rhs_Spline = new double[N_V+1];
  Mx = new double[N_V+1];
  My = new double[N_V+1];
//   Params = new double [10*N_V];
//   Param9 = new double [N_V+1];

  x = new double[N_V+1];
  y = new double[N_V+1];

//   dt =2*Pi/(N_E);
//   T = 0.0;

  m = 0;
  for(i=0;i<N_E;i++)
   {
    Me = cell[i];

    Me->GetVertex(EdgeNo[i])->GetCoords(x[m], y[m]);
    m++;

    if(TDatabase::ParamDB->USE_ISOPARAMETRIC)
    {
//   for iso mid point
     joint = cell[i]->GetJoint(EdgeNo[i]);
     interface = (TInterfaceJoint *)(joint);
     comp = interface->GetBoundComp();
     if(Me == interface->GetNeighbour(0))
        interface->GetParameters(t0, t1); // correct order
     else
        interface->GetParameters(t1, t0); // backward order

//  only second order implemented
      T = (t1-t0)/2.;
      T += t0;
      comp->GetXYofT(T, x[m], y[m]);
      m++;
    }
   }
// cout<< " test spline " <<endl;
//   end vertex is the starting vertex
  x[m] = x[0];
  y[m] = y[0];
  x[m+1] = x[1];
  y[m+1] = y[1];

  t[0] = 0.; h[0] = 0.;
  for(i=1;i<=N_Splines;i++)
   {
    h[i] = sqrt((x[i]-x[i-1])*(x[i]-x[i-1])+(y[i]-y[i-1])*(y[i]-y[i-1]));
    t[i] = t[i-1] + h[i];
   }

  h[N_Splines+1] = h[1];

  for(i=1;i<=N_Splines;i++)
   {
    a[i-1] = 2.0;
    c[i-1] = h[i]/(h[i]+h[i+1]);         // b = lambda
    b[i-1] = h[i+1]/(h[i]+h[i+1]);       // c = mu
    rhs_Spline[i-1] = 3*( b[i-1]*((x[i]-x[i-1])/h[i]) + c[i-1]*((x[i+1]-x[i])/h[i+1]) );
   }

  Solver_3dia_1(N_Splines, a, b, c, rhs_Spline, Mx+1);
  Mx[0] = Mx[N_Splines];

  for(i=1;i<=N_Splines;i++)
   {
    rhs_Spline[i-1] = 3*(  b[i-1]*((y[i]-y[i-1])/h[i]) + c[i-1]*((y[i+1]-y[i])/h[i+1]) );
   }
  Solver_3dia_1(N_Splines, a, b, c, rhs_Spline, My+1);
  My[0] = My[N_Splines];

  coll = velocity_space->GetCollection();
  N_Cells = coll->GetN_Cells();

  BeginIndex = velocity_space->GetBeginIndex();
  GlobalNumbers = velocity_space->GetGlobalNumbers();

  qf = TFEDatabase2D::GetQuadFormula1D(Gauss11Line);
  qf->GetFormulaData(N_QFPoints, weights, zeta);

  if(TDatabase::ParamDB->USE_ISOPARAMETRIC)
   T =  zeta[5] - zeta[0];
  else
   T =  zeta[N_QFPoints-1] - zeta[0];
  Local_T = new double[N_QFPoints];
  cvr= new double *[N_E];

  for(k=0;k<N_QFPoints;k++)
   {
    if(TDatabase::ParamDB->USE_ISOPARAMETRIC)
     {
//   works only for Gauss11Line 1D QuadFormula
      if(k<5)
        Local_T[k] = (zeta[k]- zeta[0])/T;
      else
        Local_T[k] = (zeta[k]- zeta[5])/T;
     }
    else
     Local_T[k] = (zeta[k]- zeta[0])/T;
//  cout<< "k " <<k << " Local_T " <<Local_T[k] <<endl;
// cout<< "k " <<k << " zeta[k] " <<zeta[k] <<endl;
   }

//    cout << "N_QFPoints " << N_QFPoints << endl;
  m = 0;
  for(i=0;i<N_E;i++)
   {
    cvr[i]= new double [N_QFPoints];
    ISpline = m*10;
    for(j=0;j<N_QFPoints;j++)
     {
      T = Local_T[j];

      if((TDatabase::ParamDB->USE_ISOPARAMETRIC) && (j ==5))
       m++;

      phi1 =  6.*T*(T-1)/h[m+1];
      phi2 = -6.*T*(T-1)/h[m+1];
      phi3 =  (3.*T*T-4.*T+1);
      phi4 =  (T*(3.*T-2));

      dx = x[m]*phi1 + x[m+1]*phi2 + Mx[m]*phi3 +Mx[m+1]*phi4;
      dy = y[m]*phi1 + y[m+1]*phi2 + My[m]*phi3 + My[m+1]*phi4;
      // cout<< m << ' '<< dx << ' ' << dy<<endl;

      phi1 =  (6.*(2.*T-1))/(h[m+1]*h[m+1]);
      phi2 = -(6.*(2.*T-1))/(h[m+1]*h[m+1]);
      phi3 = (6.*T - 4.)/(h[m+1]);
      phi4 = (6.*T - 2.)/(h[m+1]);

      dxx = x[m]*phi1 + x[m+1]*phi2 + Mx[m]*phi3 +Mx[m+1]*phi4;
      dyy = y[m]*phi1 + y[m+1]*phi2 + My[m]*phi3 + My[m+1]*phi4;
//   cout<< m << ' '<<dxx << ' ' << dyy<<endl;

      cvr[i][j] = -(dx*dyy - dxx*dy)/(pow((dx*dx + dy*dy), 1.5));
//       cvr[i][j] = - 1.0
//       cout << ' '<<cvr[i][j] << endl;
// exit(0);
     }
    m++;
// cout<<endl;
   }

// exit(0);

  for(i=0;i<N_Cells;i++)
  {
   Cell = coll->GetCell(i);
   if(Cell->GetPhase_No()==0)
     {
         // inner cells
     N_Joints = Cell->GetN_Joints();
     for(j=0;j<N_Joints;j++)
      {
       joint = Cell->GetJoint(j);
       if(joint->GetType() == InterfaceJoint)
       {
        xm = 0;
        ym = 0;
        N_Vertices = Cell->GetN_Vertices();
        if(N_Vertices != 4)
         {
          cout << "ONLY quads are allowed" << endl;
          exit(-1);
        }
        for(k=0;k<N_Vertices;k++)
        {
          Cell->GetVertex(k)->GetCoords(X, Y);
          xm += X;
          ym += Y;
        } // endfor k
        xm /= N_Vertices;
        ym /= N_Vertices;
        if(sqrt(xm*xm+ym*ym)<r2)
        {
         for(s1=0;s1<N_E;s1++)
          {
          if(Cell == cell[s1])
           break;
          if(s1 == N_E-1)
           {
            cout<< " joint not matching "<<endl;
            exit(-1);
            }
          }
//           cout<< N_Splines <<" joint " << s1 <<endl;
          // cell is inside the circle
          // cout << "cell inside" << endl;
          if(TDatabase::ParamDB->USE_ISOPARAMETRIC)
            RFid = QuadIsoparametric;
          else
            RFid = QuadBilinear;
          F_K = TFEDatabase2D::GetRefTrans2D(RFid);
          TFEDatabase2D::SetCellForRefTrans(Cell, RFid);

          // find basis functions on cell i
          FEid = velocity_space->GetFE2D(i, Cell);
          bf = TFEDatabase2D::GetBaseFunct2DFromFE2D(FEid);
          N_BF = bf->GetDimension();
          for(k=0;k<N_BF;k++)
           {
            v1[k] = 0;
            v2[k] = 0;
           }

          DOF = GlobalNumbers + BeginIndex[i];

          for(k=0;k<N_QFPoints;k++)
           {
//             t = zeta[k];
            F_K->GetTangent(j, zeta[k], ta1, ta2);
            len = sqrt(ta1*ta1+ta2*ta2);
            n1 =  ta2/len;
            n2 = -ta1/len;
            bf->GetValues(1, zeta+k, j, values);
            // cout <1< "tangent: " << ta1 << " " << ta2 << endl;

            s = weights[k]*len/We;
            for(l=0;l<N_BF;l++)
            {
              v1[l] += val[l]*s*n1*cvr[s1][k];
              v2[l] += val[l]*s*n2*cvr[s1][k];
            }
          }
          for(k=0;k<N_BF;k++)
          {
            f1[DOF[k]] += v1[k];
            f2[DOF[k]] += v2[k];
          }
        }
      } // endif
    } // endfor j
   }
  } // endfor i

   delete [] h; delete [] t; delete[] a;  delete [] b;  
   delete [] c;
   delete [] rhs_Spline;
   delete [] Mx; delete [] My;
   delete [] x; delete [] y;

   delete [] Local_T;
   for(i=0;i<N_E;i++)
     delete [] cvr[i];
   delete [] cvr;
}


void Intface_U(TFEVectFunct2D *Velocity, double *Int_intfaceU)
 {
  int i,j,k,l;
  int N_Cells, N_Joints, N_Vertices;
  TBaseCell *cell;
  TCollection *coll;
  int *BeginIndex, *GlobalNumbers, *DOF;
  TJoint *joint;
  double xm, ym, x, y;
  const double r=1;
  const double r2=r*r;
  RefTrans2D RFid;
  TRefTrans2D *F_K;
  double ta1, ta2;
  FE2D FEid;
  TBaseFunct2D *bf;
  TQuadFormula1D *qf;
  int N_QFPoints, Cell_ID;
  double *zeta, *weights;
  double *values[1], *ValuesVX, *ValuesVY;
  double val[MaxN_BaseFunctions2D];
  double v1[MaxN_BaseFunctions2D];
  double v2[MaxN_BaseFunctions2D];
  int N_BF;
  double t, s, len, n1, n2, Kappa;
  TFESpace2D *VelocitySpace;

  VelocitySpace = Velocity->GetFESpace2D();
  BeginIndex = VelocitySpace->GetBeginIndex();
  GlobalNumbers = VelocitySpace->GetGlobalNumbers();
  ValuesVX = Velocity->GetValues();
  ValuesVY = ValuesVX + Velocity->GetLength();
  values[0] = val;

  coll = VelocitySpace->GetCollection();
  N_Cells = coll->GetN_Cells();

  qf = TFEDatabase2D::GetQuadFormula1D(Gauss11Line);
  qf->GetFormulaData(N_QFPoints, weights, zeta);

  for(i=0;i<N_Cells;i++)
  {
    cell = coll->GetCell(i);
    N_Joints = cell->GetN_Joints();
    for(j=0;j<N_Joints;j++)
    {
      joint = cell->GetJoint(j);
      if(joint->GetType() == InterfaceJoint)
      {
        xm = 0;
        ym = 0;
        N_Vertices = cell->GetN_Vertices();
        if(N_Vertices != 4)
        {
          cout << "ONLY quads are allowed" << endl;
          exit(-1);
        }
        for(k=0;k<N_Vertices;k++)
        {
          cell->GetVertex(k)->GetCoords(x, y);
          xm += x;
          ym += y;
        } // endfor k
        xm /= N_Vertices;
        ym /= N_Vertices;
        if(xm*xm+ym*ym<r2)
        {
          // cell is inside the circle
          // cout << "cell inside" << endl;
          if(TDatabase::ParamDB->USE_ISOPARAMETRIC)
            RFid = QuadIsoparametric;
          else
            RFid = QuadBilinear;
          F_K = TFEDatabase2D::GetRefTrans2D(RFid);
          TFEDatabase2D::SetCellForRefTrans(cell, RFid);

          // find basis functions on cell i
          FEid = VelocitySpace->GetFE2D(i, cell);
          bf = TFEDatabase2D::GetBaseFunct2DFromFE2D(FEid);
          N_BF = bf->GetDimension();

// 	  cout<<N_BF<<" N_BF " <<endl;
          for(k=0;k<N_BF;k++)
          {
            v1[k] = 0;
            v2[k] = 0;
// 	   cout<< k<<N_BF<<" N_BF " <<endl;
          }

          DOF = GlobalNumbers + BeginIndex[i];

          for(k=0;k<N_QFPoints;k++)
          {
            t = zeta[k];
            F_K->GetTangent(j, t, ta1, ta2);
            len = sqrt(ta1*ta1+ta2*ta2);
            n1 =  ta2/len;
            n2 = -ta1/len;
            bf->GetValues(1, zeta+k, j, values);
//             cout << "tangent: " << ta1 << " " << ta2 << endl;

            s = weights[k]*len;
            for(l=0;l<N_BF;l++)
            {
	     // u.n
             Int_intfaceU[0] += (val[l]*s*n1*ValuesVX[DOF[l]]
	                         + val[l]*s*n2*ValuesVY[DOF[l]]);
	    // u
             Int_intfaceU[1] += (val[l]*s*ValuesVX[DOF[l]]
	                         + val[l]*s*ValuesVY[DOF[l]]);
            }
          }
        }
      } // endif
    } // endfor j
  } // endfor i
//   cout<< " test "<< endl;
 }



void Cal_Intface_U(TFEVectFunct2D *Velocity, double *Int_intfaceU)
{
  int i,j,k,l,m;
  TCollection *coll;
  TBaseCell *cell;
  int N_Cells, N_Joints;
  double x1, y1, x2, y2, dx, dy, x, y, xi, eta;
  double xmin, xmax, ymin, ymax;
  const double r=1.0; // circle radius
  const double r2=r*r; // square of circle radius
  double s1, s2;
  int N_Points;
  double X[4], Y[4], T[4], t;
  FE2D FEid;
  TBaseFunct2D *bf;
  TQuadFormula1D *qf;
  int N_QFPoints;
  double *zeta, *weights, *ValuesVX, *ValuesVY;
  int *BeginIndex, *GlobalNumbers, *DOF;
  RefTrans2D RFid;
  double values[MaxN_BaseFunctions2D];
  double v1[MaxN_BaseFunctions2D];
  double v2[MaxN_BaseFunctions2D];
  int N_BF;
  TFESpace2D *velocity_space;

  velocity_space = Velocity->GetFESpace2D();
  BeginIndex = velocity_space->GetBeginIndex();
  GlobalNumbers = velocity_space->GetGlobalNumbers();
  ValuesVX = Velocity->GetValues();
  ValuesVY = ValuesVX + Velocity->GetLength();

  coll = velocity_space->GetCollection();
  N_Cells = coll->GetN_Cells();

  BeginIndex = velocity_space->GetBeginIndex();
  GlobalNumbers = velocity_space->GetGlobalNumbers();

  qf = TFEDatabase2D::GetQuadFormula1D(Gauss11Line);
  qf->GetFormulaData(N_QFPoints, weights, zeta);

  for(i=0;i<N_Cells;i++)
  {
    N_Points = 0;
    cell = coll->GetCell(i);
    N_Joints = cell->GetN_Joints();
    for(j=0;j<N_Joints;j++)
    {
      cell->GetVertex(j)->GetCoords(x1, y1);
      cell->GetVertex((j+1)%N_Joints)->GetCoords(x2, y2);
      dx = x2-x1;
      dy = y2-y1;
      xmin = (x1<x2)?x1:x2;
      xmax = (x1>x2)?x1:x2;
      ymin = (y1<y2)?y1:y2;
      ymax = (y1>y2)?y1:y2;

      if(fabs(dy)<1e-8)
      {
        // horizontal joint
        if(y1*y1<=r2)
        {
          s1 = sqrt(r2-y1*y1);
          s2 = -s1;
          if(xmin<=s1 && s1<=xmax)
          {
            X[N_Points] = s1;
            Y[N_Points] = y1;
            t = atan2(Y[N_Points], X[N_Points]);
            if(t<0) t+=2*Pi;
            T[N_Points] = t;
            // cout << "cell: " << i << " joint: " << j << endl;
            // cout << X[N_Points] << " " << Y[N_Points] << " " << T[N_Points]
            // << endl;
            N_Points++;
          }
          if(xmin<=s2 && s2<=xmax)
          {
            X[N_Points] = s2;
            Y[N_Points] = y1;
            t = atan2(Y[N_Points], X[N_Points]);
            if(t<0) t+=2*Pi;
            T[N_Points] = t;
            // cout << "cell: " << i << " joint: " << j << endl;
            // cout << X[N_Points] << " " << Y[N_Points] << " " << T[N_Points]
            // << endl;
            N_Points++;
          }
        }
      }
      else
      {
        if(fabs(dx)<1e-8)
        {
          // vertical joint
          if(x1*x1<=r2)
          {
            s1 = sqrt(r2-x1*x1);
            s2 = -s1;
            if(ymin<=s1 && s1<=ymax)
            {
              X[N_Points] = x1;
              Y[N_Points] = s1;
              t = atan2(Y[N_Points], X[N_Points]);
              if(t<0) t+=2*Pi;
              T[N_Points] = t;
              // cout << "cell: " << i << " joint: " << j << endl;
              // cout << X[N_Points] << " " << Y[N_Points] << " " << T[N_Points]
              // << endl;
              N_Points++;
            }
            if(ymin<=s2 && s2<=ymax)
            {
              X[N_Points] = x1;
              Y[N_Points] = s2;
              t = atan2(Y[N_Points], X[N_Points]);
              if(t<0) t+=2*Pi;
              T[N_Points] = t;
              // cout << "cell: " << i << " joint: " << j << endl;
              // cout << X[N_Points] << " " << Y[N_Points] << " " << T[N_Points]
              // << endl;
              N_Points++;
            }
          }
        }
        else
        {
          // diagonal joint
          cout << "Diagonal joints have not been implemented yet! " << endl;
          exit(-1);
        }
      }
    } // endfor j

    if(N_Points == 2)
    {
      // if(N_Points) cout << "cell: " << i << endl;

      // correction around 0 angle
      if(fabs(T[0]-T[1])>Pi)
      {
        if(T[0]>T[1])
          T[0]-=2*Pi;
        else
          T[1]-=2*Pi;
      }

      // sort such that T is increasing
      if(T[1]<T[0])
      {
        t = T[0];    x = X[0];    y = Y[0];
        T[0] = T[1]; X[0] = X[1]; Y[0] = Y[1];
        T[1] = t;    X[1] = x;    Y[1] = y;
      }
      // for(j=0;j<N_Points;j++)
      // {
      //   cout << "X: " << X[j] << " Y: " << Y[j];
      //   cout << " T: " << T[j]*180/Pi << endl;
      // }

      // find basis functions on cell i
      FEid = velocity_space->GetFE2D(i, cell);
      bf = TFEDatabase2D::GetBaseFunct2DFromFE2D(FEid);
      N_BF = bf->GetDimension();
      for(j=0;j<N_BF;j++)
      {
        v1[j] = 0;
        v2[j] = 0;
      }

      DOF = GlobalNumbers + BeginIndex[i];

      switch(cell->GetType())
      {
        case Triangle:
          TFEDatabase2D::SetCellForRefTrans(cell, TriaAffin);
          RFid = TriaAffin;
        break;

        case Quadrangle:
          TFEDatabase2D::SetCellForRefTrans(cell, QuadBilinear);
          RFid = QuadBilinear;
        break;

        case Parallelogram:
        case Rectangle:
          TFEDatabase2D::SetCellForRefTrans(cell, QuadAffin);
          RFid = QuadAffin;
        break;

        default:
          cout << "Illegal cell shape" << endl;
          exit(-1);
      }

      for(j=0;j<N_QFPoints;j++)
      {
        t = T[0]+(T[1]-T[0])*(zeta[j]+1)/2;
        x = r*cos(t);
        y = r*sin(t);
        TFEDatabase2D::GetRefFromOrig(RFid, x, y, xi, eta);

        bf->GetDerivatives(D00, xi, eta, values);

        for(k=0;k<N_BF;k++)
        {
         // u.n
         Int_intfaceU[0] += (weights[j]*(T[1]-T[0])*0.5*values[k]*x*ValuesVX[DOF[k]]
                                 + weights[j]*(T[1]-T[0])*0.5*values[k]*y*ValuesVY[DOF[k]]);
         // u
          Int_intfaceU[1] += (weights[j]*(T[1]-T[0])*0.5*values[k]*ValuesVX[DOF[k]]
                                 + weights[j]*(T[1]-T[0])*0.5*values[k]*ValuesVY[DOF[k]]);
        }
      } // endfor j

//       for(k=0;k<N_BF;k++)
//       {
//         f1[DOF[k]] -= v1[k];
//         f2[DOF[k]] -= v2[k];
//       }
    }
    else
      if(N_Points != 0)
      {
        cout << "Only two intersection points can be handled" << endl;
        exit(-1);
      }
  } // endfor i
} // CalculateIntegrals

int main(int argc, char* argv[])
{
  TDomain *Domain = new TDomain();
  TDatabase *Database = new TDatabase();
  TFEDatabase2D *FEDatabase = new TFEDatabase2D();
  TCollection *coll, *mortarcoll = NULL;
  TBaseCell *cell, ***IntFace_Cell;
  TVertex ***IntFace_Vert;
  TFESpace2D *velocity_space, *pressure_space, *streamfunction_space;
  TFESpace2D *velocity_space_low, *pressure_space_low;
  TFESpace2D *old_u_space, *old_p_space;
  TFESpace2D *pressure_separation_space;
  TFESpace2D **USpaces, **PSpaces, **PsiSpaces;
  TOutput2D *Output;
  TJoint *joint;

  double *rhs, *sol, *oldsol, tol, tolerance, *psi, *defect, *fesol, *soldiff;
  double *rhs_low, *sol_low, *old_sol, *itmethod_sol, *itmethod_rhs;
  double *rhs_high, *nosep_p;
  int i,j,k,l, N_, Len, low, ii, N_Cells, N_Joints, N_Vertices, m1, *N_IntFace;
  int N_Rows, N_Columns, N_U, N_P, N_Unknowns, N_V, Cell_ID;
  int N_Active, N_NonActive, N_U_low,  N_P_low, N_Unknowns_low;
  double *l2u1, *l2u2, *h1u1, *h1u2;
  double *l2p, *h1p, *sd, *l_inf, vect[3], exactvect[3];
  int which;
  double DiffL2, DiffH1, xm, ym;
  char *PRM, *GEO, *MAP;
  int LEVELS, BASELEVEL;
  int ret, pde, img=1;
  double negPower, **Int_intfaceU;
  double x,y,max,min,sum;
  double RE_NR;
  double tau1, tau2;
  double errors[4], p1, p2, errors_mg[4],velo_l2;
  double t1, t2, res, res2, oldres, solver_time,residual;
  double impuls_residual,limit;
  int N_LinIter, **Edge_No;

  std::ostringstream os;
  char *PsBaseName, *GrapeBaseName, *GnuBaseName, *ReadGrapeBaseName;
  char *VtkBaseName, *MatlabBaseName;

  double *val;
  TFEFunction2D *u1, *u2, *p, *fefct[2], *StreamFct;
  TFEFunction2D *u1_low, *u2_low, *p_low, *old_p, *soldiff_fe1,*soldiff_fe2;
  TFEFunction2D **U1Array, **U2Array, **AuxFEFunctArray;
  TFEFunction2D **PArray, *AuxPArray, *separated_pressure_fe_funct;
  TFEFunction2D *separated_pressure_rhs_fe_funct;
  TFEVectFunct2D *u, **UArray, *u_low, *old_u, **AuxFEVectFunctArray;
  TFESpace2D *fesp[2], *ferhs[2];

  TAuxParam2D *aux;

  TSquareStructure2D *sqstructureA, *sqstructurePressSep;
  TStructure2D *structureB, *structureBT;
  TSquareMatrix2D *sqmatrixA, *SQMATRICES[4];
  TSquareMatrix2D *sqmatrixA11, *sqmatrixA12;
  TSquareMatrix2D *sqmatrixA21, *sqmatrixA22;
  TSquareMatrix2D **MatricesA;
  TSquareMatrix2D **MatricesA11, **MatricesA12;
  TSquareMatrix2D **MatricesA21, **MatricesA22;
  TMatrix2D *matrixB1, *matrixB2, *MATRICES[4];
  TMatrix2D *matrixB1T, *matrixB2T;
  TMatrix2D **MatricesB1, **MatricesB2, **MatricesB1T, **MatricesB2T;
  TMatrix **matrices = (TMatrix **)MATRICES;
  TSquareMatrix **sqmatrices = (TSquareMatrix **)SQMATRICES;
  MatVecProc *MatVect;
  DefectProc *Defect;
  TSquareMatrix2D *sqmatrixPressSep;

  TSquareStructure2D *sqstructureA_low;
  TStructure2D *structureB_low, *structureBT_low;
  TSquareMatrix2D *sqmatrixA_low;
  TSquareMatrix2D *sqmatrixA11_low, *sqmatrixA12_low;
  TSquareMatrix2D *sqmatrixA21_low, *sqmatrixA22_low;
  TSquareMatrix2D **MatricesA_low;
  TSquareMatrix2D **MatricesA11_low, **MatricesA12_low;
  TSquareMatrix2D **MatricesA21_low, **MatricesA22_low;
  TMatrix2D *matrixB1_low, *matrixB2_low;
  TMatrix2D *matrixB1T_low, *matrixB2T_low;
  TMatrix2D **MatricesB1_low, **MatricesB2_low, **MatricesB1T_low, **MatricesB2T_low;

  int N_SquareMatrices, N_RectMatrices;
  int N_Rhs, N_FESpaces;

  double **RhsArray;

  TNSE_MGLevel *MGLevel, *MGLevel_low;
  TNSE_MultiGrid *MG;
  
  double *RHSs[3];
  int *N_Uarray, *N_Parray;

  TDiscreteForm2D *DiscreteFormGalerkin;
  TDiscreteForm2D *DiscreteFormSDFEM;
  TDiscreteForm2D *DiscreteFormUpwind;
  TDiscreteForm2D *DiscreteFormSmagorinsky;
  
  TDiscreteForm2D *DiscreteFormNLGalerkin;
  TDiscreteForm2D *DiscreteFormNLSDFEM;
  TDiscreteForm2D *DiscreteFormNLUpwind;
  TDiscreteForm2D *DiscreteFormNLSmagorinsky;

  TDiscreteForm2D *DiscreteFormPressSep;
  TDiscreteForm2D *DiscreteFormAuxProbPressSep;

  TDiscreteForm2D *DiscreteForm;
  
  BoundCondFunct2D *BoundaryConditions[2], *BoundaryConditionsPressureSeparation[1];
  BoundValueFunct2D *BoundValues[2], *BoundaryValuesPressureSeparation[1];
  double average, hmin, hmax, *h_max;

  TItMethod *itmethod, *prec, *Auxprec, *Auxitmethod;
  int Max_It, FirstSolve;
  double omega, alpha[2], alpha_fine[2], cd,cl,dp;
  int N_Parameters=2,n_aux, **downwind;
  double Parameters[4],delta0,delta1, reatt_pt,  reatt_point[3];
  double convergence_speed, residuals[10], firsterror,lasterror;
  double firsterrorl2,lasterrorl2,p3,p4;
  int last_digit_ite,slow_conv,last_sq,mg_level,mg_type, pre_calculation=1;
  int calculations=1,ll, zerostart;
  int velocity_space_code, pressure_space_code;
  int pressure_separation = 0, N_P_sep;
  double *separated_pressure_array, *separated_pressure_aux_array;
  double *pressure_aux_array;
  double *rhsPressSep;

  TMultiGrid2D *AuxMG;
  TMGLevel2D *AuxMGLevel;
  double *Auxitmethod_sol, *Auxitmethod_rhs;
  int number_pressep;

  char Readin[] = "readin.dat";
  char Name[] = "name";
  char UString[] = "u";
  char PString[] = "p";
  char PsiString[] = "psi";
  char DString[] = "d";
  char PsepString[] = "psep";

#ifdef __BENCH__
  double Cd, Cl, dP1[3], dP2[3];
#endif

  os << " ";

//======================================================================
// read parameter file
//======================================================================
  if (argc >= 2)
    ret = Domain->ReadParam(argv[1]);
  else
    ret = Domain->ReadParam(Readin);

  if (ret == -1)
  {
    exit(-1);
  }

  RE_NR = TDatabase::ParamDB->RE_NR;
  delta0 = TDatabase::ParamDB->DELTA0;
  delta1 = TDatabase::ParamDB->DELTA1;
  convergence_speed = TDatabase::ParamDB->SC_DIV_FACTOR;

  OpenFiles();
  OutFile.setf(std::ios::scientific);

  Database->CheckParameterConsistencyNSE();
  Database->WriteParamDB();
  ExampleFile();

//======================================================================
// copy read parameters into local variables
//======================================================================

  PRM = TDatabase::ParamDB->BNDFILE;
  GEO = TDatabase::ParamDB->GEOFILE;
  MAP = TDatabase::ParamDB->MAPFILE;

  PsBaseName = TDatabase::ParamDB->PSBASENAME;
  GrapeBaseName = TDatabase::ParamDB->GRAPEBASENAME;
  GnuBaseName = TDatabase::ParamDB->GNUBASENAME;
  ReadGrapeBaseName = TDatabase::ParamDB->READGRAPEBASENAME;
  VtkBaseName = TDatabase::ParamDB->VTKBASENAME;
  MatlabBaseName = TDatabase::ParamDB->MATLABBASENAME;

  mg_type = TDatabase::ParamDB->SC_MG_TYPE_SADDLE;
  if (mg_type)
    mg_level = 0;
  else
    mg_level = -1;
  LEVELS = TDatabase::ParamDB->LEVELS;
  BASELEVEL = TDatabase::ParamDB->UNIFORM_STEPS;
  l2u1 = new double[LEVELS+1];
  l2u2 = new double[LEVELS+1];
  l2p = new double[LEVELS+1];
  h1u1 = new double[LEVELS+1];
  h1u2 = new double[LEVELS+1];
  h1p = new double[LEVELS+1];
  sd = new double[LEVELS+1];
  l_inf = new double[LEVELS+1];
  h_max = new double[LEVELS+1];

  U1Array = new TFEFunction2D*[LEVELS+1];
  U2Array = new TFEFunction2D*[LEVELS+1];
  PArray = new TFEFunction2D*[LEVELS+1];
  UArray = new TFEVectFunct2D*[LEVELS+1];

  RhsArray = new double* [LEVELS+1];
  N_Uarray = new int[LEVELS+1];
  N_Parray = new int[LEVELS+1];

  USpaces = new TFESpace2D*[LEVELS+1];
  PSpaces = new TFESpace2D*[LEVELS+1];
  PsiSpaces = new TFESpace2D*[LEVELS+1];
  OutPut("before matrices " << endl);

  N_IntFace = new int[LEVELS+1];
  IntFace_Cell = new TBaseCell**[LEVELS+1];
  IntFace_Vert = new TVertex**[LEVELS+1];
  Edge_No      = new int *[LEVELS+1];
  Int_intfaceU = new double* [LEVELS+1];

  switch(TDatabase::ParamDB->NSTYPE)
  {
    case 1:
      MatricesA = new TSquareMatrix2D*[LEVELS+1];

      MatricesB1 = new TMatrix2D*[LEVELS+1];
      MatricesB2 = new TMatrix2D*[LEVELS+1];
      MatVect = MatVect_NSE1;
      Defect = Defect_NSE1;
    break;

    case 2:
      MatricesA = new TSquareMatrix2D*[LEVELS+1];

      MatricesB1 = new TMatrix2D*[LEVELS+1];
      MatricesB2 = new TMatrix2D*[LEVELS+1];
      MatricesB1T = new TMatrix2D*[LEVELS+1];
      MatricesB2T = new TMatrix2D*[LEVELS+1];
      MatVect = MatVect_NSE2;
      Defect = Defect_NSE2;
    break;

    case 3:
      MatricesA11 = new TSquareMatrix2D*[LEVELS+1];
      MatricesA12 = new TSquareMatrix2D*[LEVELS+1];
      MatricesA21 = new TSquareMatrix2D*[LEVELS+1];
      MatricesA22 = new TSquareMatrix2D*[LEVELS+1];

      MatricesB1 = new TMatrix2D*[LEVELS+1];
      MatricesB2 = new TMatrix2D*[LEVELS+1];
      MatVect = MatVect_NSE3;
      Defect = Defect_NSE3;
    break;

    case 4:
      MatricesA11 = new TSquareMatrix2D*[LEVELS+1];
      MatricesA12 = new TSquareMatrix2D*[LEVELS+1];
      MatricesA21 = new TSquareMatrix2D*[LEVELS+1];
      MatricesA22 = new TSquareMatrix2D*[LEVELS+1];

      MatricesB1 = new TMatrix2D*[LEVELS+1];
      MatricesB2 = new TMatrix2D*[LEVELS+1];
      MatricesB1T = new TMatrix2D*[LEVELS+1];
      MatricesB2T = new TMatrix2D*[LEVELS+1];
      MatVect = MatVect_NSE4;
      Defect = Defect_NSE4;
    break;
  } // endswitch
  downwind = new int*[LEVELS+1];

//======================================================================
// read boundary parameterization and initialize coarse grid
//======================================================================
  Domain->Init(PRM, GEO);
//     if (TDatabase::ParamDB->CONVERT_QUAD_TO_TRI)
//     {
//         Domain->ConvertQuadToTri(TDatabase::ParamDB->CONVERT_QUAD_TO_TRI);

//       OutPut("ConvertQuadToTri" << endl);
// exit(0);
//     }
//======================================================================
// do some special mortar stuff
//======================================================================
#ifdef __MORTAR__
  int N_Mortar;
  TVelocity_SpaceD *fespace_mortar;
  TStructure2D *struct_mortar;
  TMatrix2D *matrix_mortar;

  if (!strcmp(MAP, "NO_MAP_FILE"))
  {
    OutPut("switch off 'MORTAR = -D__MORTAR__' in the makefile !!!" << endl);
    OutPut("set 'MORTAR = ' in the makefile !!!" << endl);      
    exit(1);
  }
  else
  {
    Domain->ReadMapFile(MAP, Database);
    //Domain->TestMortar();

    if (TDatabase::ParamDB->CONVERT_QUAD_TO_TRI)
    {
      Domain->ConvertQuadToTri();
      OutPut("ConvertQuadToTri" << endl);
// exit(0);
    }

    Domain->RegRefineSub(0);
    Domain->PS("Grid1.ps",It_Finest,0);
    //Domain->PS("Grid2.ps",It_Finest,0);
    //Domain->RegRefineAll();
  }
#endif // __MORTAR__

//======================================================================
// initialize all discrete forms
//======================================================================
 OutPut("before discrete forms " << endl);

  if( (TDatabase::ParamDB->DISCTYPE == 2) &&
     ((TDatabase::ParamDB->NSTYPE == 1) ||
      (TDatabase::ParamDB->NSTYPE == 3) ))
  {
    OutPut("DISCTYPE=2 works only with NSTYPE=2 or NSTYPE=4!" << endl);
    Error("DISCTYPE=2 works only with NSTYPE=2 or NSTYPE=4!" << endl);
    return -1;
  }

  InitializeDiscreteForms(
    DiscreteFormGalerkin, DiscreteFormSDFEM,
    DiscreteFormUpwind, DiscreteFormSmagorinsky,
    DiscreteFormNLGalerkin, DiscreteFormNLSDFEM,
    DiscreteFormNLUpwind, DiscreteFormNLSmagorinsky,
    DiscreteFormPressSep,
    DiscreteFormAuxProbPressSep,
    LinCoeffs, TDatabase::ParamDB->NSTYPE);

 OutPut("after discrete forms " << endl);
  BoundaryConditions[0] = BoundCondition;
  BoundaryConditions[1] = BoundCondition;

  BoundValues[0] = U1BoundValue;
  BoundValues[1] = U2BoundValue;

  BoundaryConditionsPressureSeparation[0] = BoundaryConditionPressSep;
  BoundaryValuesPressureSeparation[0] = BoundaryValuePressSep;
  

// refine up to user defined coarsest level

  for(i=0;i<TDatabase::ParamDB->SC_COARSEST_LEVEL_SADDLE;i++)
    Domain->RegRefineAll();

// initialize solver parameters

  limit = TDatabase::ParamDB->SC_NONLIN_RES_NORM_MIN_SADDLE;
  Max_It = TDatabase::ParamDB->SC_NONLIN_MAXIT_SADDLE;
  omega = TDatabase::ParamDB->SC_NONLIN_DAMP_FACTOR_SADDLE;
  velocity_space_code =   TDatabase::ParamDB->VELOCITY_SPACE;

  Parameters[0] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_SADDLE;
  Parameters[1] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_FINE_SADDLE;
  // Parameters[2] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_COARSE_SADDLE;

  if(TDatabase::ParamDB->SOLVER_TYPE == GMG)
  {
    i=1;
    MG = new TNSE_MultiGrid(i, N_Parameters, Parameters);
  }
  Parameters[2] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_SCALAR;
  Parameters[3] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_FINE_SCALAR;
  if ((TDatabase::ParamDB->PRESSURE_SEPARATION==3)||
      (TDatabase::ParamDB->PRESSURE_SEPARATION==4)||
      (TDatabase::ParamDB->PRESSURE_SEPARATION==5))
  {
    AuxMG = new TMultiGrid2D(i, N_Parameters, Parameters+2);
  }

  FirstSolve = TDatabase::ParamDB->SC_FIRST_SOLUTION_LEVEL_SADDLE;
 OutPut("before loop " << endl);


double r2 = 1.0; // inner domain radius
//======================================================================
// loop over all levels
//======================================================================

  for(i=0;i<LEVELS;i++)
  {
    mg_level++;
    OutPut("*******************************************************" << endl);
    OutPut("******           GEOMETRY  LEVEL ");
    OutPut(i << "              *******" << endl);
    OutPut("******           MULTIGRID LEVEL ");
    OutPut(mg_level << "              *******" << endl);
    OutPut("*******************************************************" << endl);
    solver_time = 0.0;
    N_LinIter = 0;
    OutPut("memory before: " << setw(10) << GetMemory() << endl);
    for (j=0;j<10;j++)
      residuals[j]=1e10;
    slow_conv = 0;

    // refine grid if level is greater than 0
    if (i)
      Domain->RegRefineAll();
//  cout<< "test cell ID 1 "<<endl;
    coll=Domain->GetCollection(It_Finest, 0);
    Output = new TOutput2D(2, 2, 1, 1,Domain);
    cout << endl << endl;

    if(TDatabase::ParamDB->WRITE_PS)
    {
      // write grid into an Postscript file
      os.seekp(std::ios::beg);
      os << PsBaseName << i << ".ps" << ends;
      Domain->PS(os.str().c_str(),It_Finest,0);
    }


// ***************************************************************
//  for calculating isointegrals on interfaces between two doamins
// ***************************************************************
//  Split the domain into two phases
  N_Cells = coll->GetN_Cells();
if(!(TDatabase::ParamDB->REACTOR_P28==12345 || TDatabase::ParamDB->REACTOR_P28==123))
 {
  m1 = 0;
  for(l=0;l<N_Cells;l++)
  {
   cell = coll->GetCell(l);
   xm = 0.;  ym = 0.;
   N_Vertices = cell->GetN_Vertices();
   for(k=0;k<N_Vertices;k++)
    {
     cell->GetVertex(k)->GetCoords(x, y);
     xm += x;
     ym += y;
    } // endfor k
   xm /= N_Vertices;
   ym /= N_Vertices;
   if(sqrt(xm*xm+ym*ym)<r2)
    {
     cell->SetPhase_No(0);  // inner phase
     N_Joints = cell->GetN_Joints();
     for(j=0;j<N_Joints;j++)
      {
       joint = cell->GetJoint(j);
       if(joint->GetType() == InterfaceJoint)
         m1++;
      }
     }
    else
     {
      cell->SetPhase_No(1);  // outer phase
//         cout<< "test cell ID 1 "<<endl;
    }
   }

// data for spline approximation
 N_IntFace[i] = m1;
 cout <<"Number of interface vertices " << N_IntFace[i] << endl;
 IntFace_Cell[i] = new TBaseCell*[m1];
 IntFace_Vert[i] = new TVertex*[m1];
 Edge_No[i]      = new int [m1];

  m1 = 0;
  for(l=0;l<N_Cells;l++)
  {
   cell = coll->GetCell(l);
   Cell_ID = cell->GetPhase_No();
    if(Cell_ID==0) // inner cells
     {
      N_Joints = cell->GetN_Joints();
      for(j=0;j<N_Joints;j++)
       {
        joint = cell->GetJoint(j);
        if(joint->GetType() == InterfaceJoint)
         {
          IntFace_Cell[i][m1] = cell;
          IntFace_Vert[i][m1] = cell->GetVertex(j);
          IntFace_Vert[i][m1]->GetCoords(x, y);
//           cout << x << ' ' << y<< ' '<<endl;
          Edge_No[i][m1++] = j;
         }
       }
     }
   }

   // sort interface vertices in anti-clockwise
   SortIntfaceVert(IntFace_Cell[i], IntFace_Vert[i], Edge_No[i], N_IntFace[i]);
 }
// ***************************************************************
//  end - data for spline
// ***************************************************************




#ifdef __MORTAR__
    mortarcoll = Domain->GetMortarColl(It_Mortar1, MAX_ItLevel);
    Domain->InitMortarJoints(It_Mortar1, MAX_ItLevel, mortarcoll);
#endif

    // get spaces for low order disc on finest geo grid
    if (mg_type==1)
    {
      velocity_space_low = new TFESpace2D(coll,Name,UString,BoundCondition, 
                                    Non_USpace,1, mortarcoll);
      pressure_space_low = new TFESpace2D(coll,Name,PString,BoundCondition, 
                                    DiscP_PSpace,0, mortarcoll);
    }
    // get spaces of high order disc on finest geo grid
    if ((i>=FirstSolve)||(mg_type==0))
      GetVelocityAndPressureSpace(coll,BoundCondition,
                                  mortarcoll, velocity_space,
                                  pressure_space, &pressure_space_code,
                                  TDatabase::ParamDB->VELOCITY_SPACE,
                                  TDatabase::ParamDB->PRESSURE_SPACE);   
        

    TDatabase::ParamDB->INTERNAL_PRESSURE_SPACE = pressure_space_code;	
    
    // fe space for stream function
    streamfunction_space = new TFESpace2D(coll,Name,PsiString,BoundCondition, 
                                          1, mortarcoll);



#ifdef __MORTAR__
    fespace_mortar = new TVelocity_SpaceD(mortarcoll, "mortar space",
                     "mortar space", velocity_space);

    N_Mortar =  fespace_mortar->GetN_DegreesOfFreedom();

    struct_mortar = new TStructure2D(fespace_mortar, velocity_space);
    matrix_mortar = new TMatrix2D(struct_mortar);
#endif

    // build fespace hierarchy
    // set values and pointers for low order fe space
    if (mg_type==1)
    {
      USpaces[i] = velocity_space_low;
      PSpaces[i] = pressure_space_low;
      N_U_low = velocity_space_low->GetN_DegreesOfFreedom();
      N_U = N_U_low;
      N_P_low = pressure_space_low->GetN_DegreesOfFreedom();
      N_P = N_P_low;      
      N_Uarray[i] = velocity_space_low->GetN_DegreesOfFreedom();
      N_Parray[i] = pressure_space_low->GetN_DegreesOfFreedom();
    }
    // set values and pointers for high order fe space
    if ((i>=FirstSolve)||(mg_type==0))
    {
      USpaces[mg_level] = velocity_space;
      PSpaces[mg_level] = pressure_space;    
      N_U = velocity_space->GetN_DegreesOfFreedom();
      N_P = pressure_space->GetN_DegreesOfFreedom();
      N_Uarray[mg_level] = N_U;
      N_Parray[mg_level] = N_P;
      N_Active = velocity_space->GetActiveBound();
      N_NonActive = N_U - N_Active;
      PsiSpaces[i] = streamfunction_space;
      N_V = streamfunction_space->GetN_DegreesOfFreedom();    
    }

    if (TDatabase::ParamDB->PRESSURE_SEPARATION>=1)
       //  ((i>=FirstSolve)||(mg_type==0)))
    {
       if (i<FirstSolve) 
           pressure_space_code=0;
      OutPut("pressure_space_code " << pressure_space_code << endl); 
      // allocate finite element space for separated pressure
      switch (pressure_space_code)
      {
         case 0:
            pressure_separation_space = new TFESpace2D(coll,Name,PsiString,BoundaryConditionPressSep, 
                                                       1, mortarcoll);
            break;
         case 1:
         case -11:
            pressure_separation_space = new TFESpace2D(coll,Name,PsiString,BoundConditionVMM, 
                                                       2, mortarcoll);
            break;
         default:
            OutPut("case for pressure_space_code not implemented" << endl);
            exit(4711);
      }
      N_P_sep = pressure_separation_space->GetN_DegreesOfFreedom();
      separated_pressure_array = new double[N_P_sep];
      memset(separated_pressure_array,0, N_P_sep*SizeOfDouble);
      separated_pressure_aux_array = new double[N_P_sep];
      memset(separated_pressure_aux_array,0, N_P_sep*SizeOfDouble);
      separated_pressure_fe_funct = new TFEFunction2D(pressure_separation_space, 
                                                      PsepString, PsepString,separated_pressure_array,
                                                      N_P_sep);
       if (i<FirstSolve) 
          N_P = N_P_sep;
       nosep_p = new double[N_P];
      if (TDatabase::ParamDB->PRESSURE_SEPARATION==1)
         pressure_separation = 0;
      else
         pressure_separation = 1;
      if ((TDatabase::ParamDB->PRESSURE_SEPARATION==3)||
          (TDatabase::ParamDB->PRESSURE_SEPARATION==4)||
          (TDatabase::ParamDB->PRESSURE_SEPARATION==5))
      {
         // allocate matrices 
         sqstructurePressSep = new TSquareStructure2D(pressure_separation_space);
         sqmatrixPressSep = new TSquareMatrix2D(sqstructurePressSep);
         // rhs for auxiliary problem
         rhsPressSep = new double[N_P_sep];
         memset(rhsPressSep, 0, N_P_sep*SizeOfDouble);
         separated_pressure_rhs_fe_funct = new TFEFunction2D(pressure_separation_space, 
                                                             PsepString, PsepString,rhsPressSep,
                                                             N_P_sep);
         // auxiliary array for prolongation
         pressure_aux_array = new double[N_P];
         memset(pressure_aux_array,0, N_P*SizeOfDouble);
         //OutPut("NP " << N_P << " NPSEP " << N_P_sep << endl);
         // allocate multigrid level
         n_aux = 4;
         AuxMGLevel = new TMGLevel2D(i, sqmatrixPressSep, 
                                     rhsPressSep, separated_pressure_array, 
                                     n_aux, NULL);
         AuxMG->AddLevel(AuxMGLevel);
         if (i <= FirstSolve)
         {
            // the matrix for the auxiliary problem must be build here
            // it is needed on the lower levels of the multigrid method
            OutPut("assemble pressure separation problem"<<endl);
            N_SquareMatrices = 1;
            N_RectMatrices = 0;          
            N_Rhs = 1;
            SQMATRICES[0] = sqmatrixPressSep;
            SQMATRICES[0]->Reset();
            fesp[0] = pressure_separation_space;
            ferhs[0] = pressure_separation_space;
            // form of rhs on lower levels not important
            // use for for PRESSURE_SEPARATION 3
            number_pressep = TDatabase::ParamDB->PRESSURE_SEPARATION;
            TDatabase::ParamDB->PRESSURE_SEPARATION=3;
            N_FESpaces = 1;
            aux = new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL); 
            RHSs[0] = rhsPressSep;
            
            // assemble
            Assemble2D(N_FESpaces, fesp, 
                       N_SquareMatrices, SQMATRICES, 
                       0, NULL, 
                       N_Rhs, RHSs, ferhs,
                       DiscreteFormAuxProbPressSep, 
                       BoundaryConditionsPressureSeparation, 
                       BoundaryValuesPressureSeparation, 
                       aux);
            delete aux;
            TDatabase::ParamDB->PRESSURE_SEPARATION =  number_pressep; 
         }
      }
    }

    // build matrices for high order disc
    if ((i>=FirstSolve)||(mg_type==0))
    {
      // matrix structures
      structureB = new TStructure2D(pressure_space, velocity_space);
      structureBT = new TStructure2D(velocity_space, pressure_space);
      sqstructureA = new TSquareStructure2D(velocity_space);
      sqstructureA->Sort();
      
      // allocate matrices
      switch(TDatabase::ParamDB->NSTYPE)
      {
        case 1:
          matrixB1 = new TMatrix2D(structureB);
          matrixB2 = new TMatrix2D(structureB);
          
          MatricesB1[mg_level] = matrixB1;
          MatricesB2[mg_level] = matrixB2;
          
          sqmatrixA = new TSquareMatrix2D(sqstructureA);
          
          MatricesA[mg_level] = sqmatrixA;
          break;
          
        case 2:
          matrixB1 = new TMatrix2D(structureB);
          matrixB2 = new TMatrix2D(structureB);
          matrixB1T = new TMatrix2D(structureBT);
          matrixB2T = new TMatrix2D(structureBT);
          
          MatricesB1[mg_level] = matrixB1;
          MatricesB2[mg_level] = matrixB2;
          MatricesB1T[mg_level] = matrixB1T;
          MatricesB2T[mg_level] = matrixB2T;
          
          sqmatrixA = new TSquareMatrix2D(sqstructureA);
          
          MatricesA[mg_level] = sqmatrixA;
          break;
          
        case 3:
          matrixB1 = new TMatrix2D(structureB);
          matrixB2 = new TMatrix2D(structureB);
          
          MatricesB1[mg_level] = matrixB1;
          MatricesB2[mg_level] = matrixB2;

          sqmatrixA11 = new TSquareMatrix2D(sqstructureA);
          sqmatrixA12 = new TSquareMatrix2D(sqstructureA);
          sqmatrixA21 = new TSquareMatrix2D(sqstructureA);
          sqmatrixA22 = new TSquareMatrix2D(sqstructureA);
          
          MatricesA11[mg_level] = sqmatrixA11;
          MatricesA12[mg_level] = sqmatrixA12;
          MatricesA21[mg_level] = sqmatrixA21;
          MatricesA22[mg_level] = sqmatrixA22;
          break;
          
        case 4:
          matrixB1 = new TMatrix2D(structureB);
          matrixB2 = new TMatrix2D(structureB);
          matrixB1T = new TMatrix2D(structureBT);
          matrixB2T = new TMatrix2D(structureBT);
          
          MatricesB1[mg_level] = matrixB1;
          MatricesB2[mg_level] = matrixB2;
          MatricesB1T[mg_level] = matrixB1T;
          MatricesB2T[mg_level] = matrixB2T;
          
          sqmatrixA11 = new TSquareMatrix2D(sqstructureA);
          sqmatrixA12 = new TSquareMatrix2D(sqstructureA);
          sqmatrixA21 = new TSquareMatrix2D(sqstructureA);
          sqmatrixA22 = new TSquareMatrix2D(sqstructureA);
          
          MatricesA11[mg_level] = sqmatrixA11;
          MatricesA12[mg_level] = sqmatrixA12;
          MatricesA21[mg_level] = sqmatrixA21;
          MatricesA22[mg_level] = sqmatrixA22;
          break;
      }
    }

    // build matrices for low order disc
    if (mg_type==1)
    {
      // matrix structures
      structureB_low = new TStructure2D(pressure_space_low, velocity_space_low);
      structureBT_low = new TStructure2D(velocity_space_low, pressure_space_low);
      sqstructureA_low = new TSquareStructure2D(velocity_space_low);
      sqstructureA_low->Sort();
      
      switch(TDatabase::ParamDB->NSTYPE)
      {
        case 1:
          matrixB1_low = new TMatrix2D(structureB_low);
          matrixB2_low = new TMatrix2D(structureB_low);
          
          MatricesB1[i] = matrixB1_low;
          MatricesB2[i] = matrixB2_low;
          
          sqmatrixA_low = new TSquareMatrix2D(sqstructureA_low);
          
          MatricesA[i] = sqmatrixA_low;
          break;
          
        case 2:
          matrixB1_low = new TMatrix2D(structureB_low);
          matrixB2_low = new TMatrix2D(structureB_low);
          matrixB1T_low = new TMatrix2D(structureBT_low);
          matrixB2T_low = new TMatrix2D(structureBT_low);
          
          MatricesB1[i] = matrixB1_low;
          MatricesB2[i] = matrixB2_low;
          MatricesB1T[i] = matrixB1T_low;
          MatricesB2T[i] = matrixB2T_low;
          
          sqmatrixA_low = new TSquareMatrix2D(sqstructureA_low);
          
          MatricesA[i] = sqmatrixA_low;
          break;
          
        case 3:
          matrixB1_low = new TMatrix2D(structureB_low);
          matrixB2_low = new TMatrix2D(structureB_low);

          MatricesB1[i] = matrixB1_low;
          MatricesB2[i] = matrixB2_low;
          
          sqmatrixA11_low = new TSquareMatrix2D(sqstructureA_low);
          sqmatrixA12_low = new TSquareMatrix2D(sqstructureA_low);
          sqmatrixA21_low = new TSquareMatrix2D(sqstructureA_low);
          sqmatrixA22_low = new TSquareMatrix2D(sqstructureA_low);
          
          MatricesA11[i] = sqmatrixA11_low;
          MatricesA12[i] = sqmatrixA12_low;
          MatricesA21[i] = sqmatrixA21_low;
          MatricesA22[i] = sqmatrixA22_low;
          break;
          
        case 4:
          matrixB1_low = new TMatrix2D(structureB_low);
          matrixB2_low = new TMatrix2D(structureB_low);
          matrixB1T_low = new TMatrix2D(structureBT_low);
          matrixB2T_low = new TMatrix2D(structureBT_low);
          
          MatricesB1[i] = matrixB1_low;
          MatricesB2[i] = matrixB2_low;
          MatricesB1T[i] = matrixB1T_low;
          MatricesB2T[i] = matrixB2T_low;
          
          sqmatrixA11_low = new TSquareMatrix2D(sqstructureA_low);
          sqmatrixA12_low = new TSquareMatrix2D(sqstructureA_low);
          sqmatrixA21_low = new TSquareMatrix2D(sqstructureA_low);
          sqmatrixA22_low = new TSquareMatrix2D(sqstructureA_low);
          
          MatricesA11[i] = sqmatrixA11_low;
          MatricesA12[i] = sqmatrixA12_low;
          MatricesA21[i] = sqmatrixA21_low;
          MatricesA22[i] = sqmatrixA22_low;
          break;
      }
    } // end if (mg_type==1)

#ifdef __MORTAR__
    N_Unknowns = 2*N_U + N_P + 2* N_Mortar;
    OutPut("dof mortar   : "<< setw(10) << 2*N_Mortar << endl);
#else
    N_Unknowns = 2*N_U + N_P;
    if (mg_type==1)
      N_Unknowns_low = 2*N_U_low + N_P_low;
#endif

    coll->GetHminHmax(&hmin,&hmax);
    h_max[i] = hmax;
    OutPut("h_max : " << hmax << " h_max : " << hmax << " sqrt(h_min) : " << sqrt(hmin)
    << " sqrt(h_min) : " << sqrt(hmin)<<endl);
    OutPut("dof velocity : "<< setw(10) << 2* N_U << endl);
    OutPut("dof pressure : "<< setw(10) << N_P << endl);
    OutPut("dof all      : "<<  setw(10) << N_Unknowns  << endl);
    if (mg_type==1)
    {
      OutPut("dof low order disc     : "<<  setw(10) << N_Unknowns_low  << endl);

      // initialize solver
      // low order disc
      rhs_low = new double[N_Unknowns_low];
      memset(rhs_low, 0, N_Unknowns_low*SizeOfDouble);
      RhsArray[i] = rhs_low;
      sol_low = new double[N_Unknowns_low];
      memset(sol_low, 0, N_Unknowns_low*SizeOfDouble);
    }

    // high order disc
    if ((i>=FirstSolve)||(mg_type==0))
    {
      rhs_high = new double[N_Unknowns];
      memset(rhs_high, 0, N_Unknowns*SizeOfDouble);  
      RhsArray[mg_level] = rhs_high;
      sol = new double[N_Unknowns];
      oldsol = new double[N_Unknowns];
      memset(sol, 0, N_Unknowns*SizeOfDouble);
      memset(oldsol, 0, N_Unknowns*SizeOfDouble);
//        for (k=0;k<N_U;k++)
//         sol[k] = 5;
//       for (k=0;k<N_U;k++)
//       sol[N_U+k] = -4;
    }

    // build multigrid level(s)
    // ( A B' )
    // ( B 0  )
    switch(TDatabase::ParamDB->SOLVER_TYPE)
    {
      case AMG:
        low = mg_level;
      break;

      case GMG:
        // coarsest grid number
        low = 0;
        // determine number of auxiliary arrays
        if ((TDatabase::ParamDB->SC_STEP_LENGTH_CONTROL_ALL_SADDLE)
            || (TDatabase::ParamDB->SC_STEP_LENGTH_CONTROL_FINE_SADDLE))
          n_aux=4;
        else
          n_aux=2;

        if (i==0)
        {
           alpha[0] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_SADDLE;
           alpha[1] = TDatabase::ParamDB->SC_GMG_DAMP_FACTOR_SADDLE;
        }
        else
        {
           alpha[0] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_SADDLE;
           alpha[1] = TDatabase::ParamDB->SC_GMG_DAMP_FACTOR_SADDLE;
        }
        if (mg_type==1)
        {
           alpha_fine[0] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_SADDLE;
           alpha_fine[1] = TDatabase::ParamDB->SC_GMG_DAMP_FACTOR_FINE_SADDLE;
        }
        else
        {
           alpha_fine[0] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_SADDLE;
           alpha_fine[1] = TDatabase::ParamDB->SC_GMG_DAMP_FACTOR_SADDLE;
        }
           
     downwind[i] = new int[coll->GetN_Cells()];
     for (j=0;j<coll->GetN_Cells();j++)
        downwind[i][j] = j;
#ifdef __DOWNWIND__
     DownwindNumberingCells(coll, downwind[i]);
#endif

        // build fe multigrid levels
        switch(TDatabase::ParamDB->NSTYPE)
        {
          case 1:
            // low order disc
            if (mg_type==1)
            {
              MGLevel_low = new TNSE_MGLevel1(i, sqmatrixA_low,
                            matrixB1_low, matrixB2_low,
                            structureBT_low,
                            rhs_low,  sol_low,
                            n_aux, alpha[1], -1, 0, NULL );
              if (i==0)
                MG->AddLevel(MGLevel_low);
              else
                MG->ReplaceLevel(i,MGLevel_low);
            }
            // high order disc
            if ((i>=FirstSolve)||(mg_type==0))
            {
              MGLevel = new TNSE_MGLevel1(mg_level, sqmatrixA,
                            matrixB1, matrixB2,
                            structureBT,
                            rhs_high,  sol,  n_aux, alpha_fine[1],
                            velocity_space_code, pressure_space_code, NULL);
              MG->AddLevel(MGLevel);
            }
          break;

          case 2:
            if (mg_type==1)
            {
              MGLevel_low = new TNSE_MGLevel2(i, sqmatrixA_low,
                            matrixB1_low, matrixB2_low,
                            matrixB1T_low, matrixB2T_low,
                            rhs_low, sol_low,
                            n_aux, alpha[1], -1, 0, NULL);
            if (i==0)
              MG->AddLevel(MGLevel_low);
            else
              MG->ReplaceLevel(i,MGLevel_low);
            }
            // high order disc
            if ((i>=FirstSolve)||(mg_type==0))
            {
              MGLevel = new TNSE_MGLevel2(mg_level, sqmatrixA,
                        matrixB1, matrixB2,
                        matrixB1T, matrixB2T,
                        rhs_high, sol, n_aux, alpha_fine[1],
                        velocity_space_code, pressure_space_code, NULL);
              MG->AddLevel(MGLevel);
            }
            break;

          case 3:
            if (mg_type==1)
            {
              MGLevel_low = new TNSE_MGLevel3(i,
                            sqmatrixA11_low, sqmatrixA12_low,
                            sqmatrixA21_low, sqmatrixA22_low,
                            matrixB1_low, matrixB2_low,
                            structureBT_low,
                            rhs_low, sol_low,
                            n_aux, alpha[1], -1, 0, NULL);
              if (i==0)
                MG->AddLevel(MGLevel_low);
              else
                MG->ReplaceLevel(i,MGLevel_low);
            }
            // high order disc
            if ((i>=FirstSolve)||(mg_type==0))
            {
              MGLevel = new TNSE_MGLevel3(mg_level,
                        sqmatrixA11, sqmatrixA12,
                        sqmatrixA21, sqmatrixA22, matrixB1, matrixB2,
                        structureBT,
                        rhs_high, sol, n_aux, alpha_fine[1],
                        velocity_space_code, pressure_space_code, NULL);
              MG->AddLevel(MGLevel);
            }
          break;

          case 4:
            // low order disc
            if (mg_type==1)
            {
              MGLevel_low = new TNSE_MGLevel4(i,
                            sqmatrixA11_low,  sqmatrixA12_low,
                            sqmatrixA21_low, sqmatrixA22_low,
                            matrixB1_low, matrixB2_low,
                            matrixB1T_low, matrixB2T_low,
                            rhs_low, sol_low,
                            n_aux, alpha[1], -1, 0, NULL);
              if (i==0)
                MG->AddLevel(MGLevel_low);
              else
                MG->ReplaceLevel(i,MGLevel_low);
            }
            // high order disc
            if ((i>=FirstSolve)||(mg_type==0))
            {
              MGLevel = new TNSE_MGLevel4(mg_level, sqmatrixA11,  sqmatrixA12,
                        sqmatrixA21, sqmatrixA22, matrixB1, matrixB2,
                        matrixB1T, matrixB2T,
                        rhs_high, sol, n_aux, alpha_fine[1],
                        velocity_space_code, pressure_space_code, NULL);
              MG->AddLevel(MGLevel);
            }
          break;
        } // end switch(NSTYPE)
      break;
    }

#ifdef __MORTAR__
    Assemble(matrix_mortar);
#endif


    // build new fe functions
    // high order fe space
    if ((i>=FirstSolve)||(mg_type==0))
    {
      u = new TFEVectFunct2D(velocity_space, UString, UString, sol, N_U, 2);
      u1 = u->GetComponent(0);
      u2 = u->GetComponent(1);
      p = new TFEFunction2D(pressure_space, PString, PString, sol+2*N_U, N_P);
      
      U1Array[mg_level] = u1;
      U2Array[mg_level] = u2;
      PArray[mg_level] = p;
      UArray[mg_level] = u;
    }

    // low order fe space
    if (mg_type==1)
    {
      u_low = new TFEVectFunct2D(velocity_space_low, UString, UString, sol_low, N_U_low, 2);
      u1_low = u_low->GetComponent(0);
      u2_low = u_low->GetComponent(1);
      p_low = new TFEFunction2D(pressure_space_low, PString, PString, sol_low+2*N_U_low, N_P_low);
      U1Array[i] = u1_low;
      U2Array[i] = u2_low;
      PArray[i] = p_low;
      UArray[i] = u_low;
    }

    if ((i>=FirstSolve)||(mg_type==0)) // CHECK THIS !!!
    {
      Output->AddFEVectFunct(u);
      Output->AddFEFunction(p);
    }
    // prolongation, to get a good starting iterate 
    if(i && i>FirstSolve)
    {
      Prolongate(old_u_space, USpaces[mg_level],
                 old_u->GetComponent(0)->GetValues(), U1Array[mg_level]->GetValues(),
                 oldsol);
      Prolongate(old_u_space, USpaces[mg_level],
                 old_u->GetComponent(1)->GetValues(), U2Array[mg_level]->GetValues(),
                 oldsol);
      Prolongate(old_p_space, PSpaces[mg_level],
                 old_p->GetValues(), PArray[mg_level]->GetValues(), 
                 oldsol);
      
      if (TDatabase::ParamDB->SOLVER_TYPE==AMG)
      {
        delete USpaces[i-1];
        delete PSpaces[i-1];
        delete U1Array[i-1]->GetValues();
      }
      // copy current solution for assembling the nonlinear terms
      memcpy(oldsol,sol, N_Unknowns*SizeOfDouble);
 
      if (mg_type==1)
      {
        delete old_sol;
        delete old_u;
        delete old_p;
        delete old_u_space;
        delete old_p_space;
      }            
    } // end of prolongate

    if (TDatabase::ParamDB->P9==123456789)
    {
      for (k=0;k<2*N_U;k++)
        sol[k] = 1.0;
      for (k=0;k<N_P;k++)      
        sol[2*N_U+k] = 0.0;
      fesol = new double[N_Unknowns];
      memset(fesol, 0, N_Unknowns*SizeOfDouble);
      soldiff = new double[N_Unknowns];
      memset(soldiff, 0, N_Unknowns*SizeOfDouble);
      soldiff_fe1 = new TFEFunction2D(velocity_space, DString, DString, soldiff,N_U);
      soldiff_fe2 = new TFEFunction2D(velocity_space, DString, DString, soldiff+N_U,N_U);
      pre_calculation = 1;
      calculations = 2;
      TDatabase::ParamDB->SC_LIN_RES_NORM_MIN_SADDLE = 1e-13;
    }

    // restrict solution to all grids 
    if (TDatabase::ParamDB->SOLVER_TYPE==GMG)
      MG->RestrictToAllGrids();      

    // if no solution on this grid, continue
    if(FirstSolve>i) 
      continue;
 
    if(TDatabase::ParamDB->READ_GRAPE_FILE)
    {
      AuxFEFunctArray = new TFEFunction2D*[2];
      AuxPArray =  new TFEFunction2D(velocity_space, PString, PString, oldsol, N_V);
      AuxFEFunctArray[0] = PArray[mg_level];
      AuxFEFunctArray[1] = AuxPArray;
      AuxFEVectFunctArray = new TFEVectFunct2D*[1];
      AuxFEVectFunctArray[0] = UArray[mg_level];
      ReadGrapeFile(ReadGrapeBaseName, 2, 1, AuxFEFunctArray,AuxFEVectFunctArray);
      delete AuxPArray;
      TDatabase::ParamDB->READ_GRAPE_FILE = 0;
      OutPut("u " << Ddot(2*N_U,sol,sol)<< endl);
      // for (ii=0;ii<N_Unknowns;ii++)
      //  OutPut(ii << " " << sol[ii] << endl);
      OutPut("p " << Ddot(N_P,sol+2*N_U,sol+2*N_U)<< endl);
      memcpy(oldsol,sol, N_Unknowns*SizeOfDouble);
      if (TDatabase::ParamDB->SOLVER_TYPE==GMG)
        MG->RestrictToAllGrids();      
     }

    // build the discretizations 
    for(k=low;k<=mg_level;k++)
    {
      rhs = RhsArray[k];
      N_U = N_Uarray[k];
      N_P = N_Parray[k];
      N_Active = USpaces[k]->GetActiveBound();
      N_NonActive = N_U - N_Active;

      RHSs[0] = rhs;
      RHSs[1] = rhs + N_U;
      RHSs[2] = rhs + 2*N_U;
      memset(rhs, 0, (2*N_U+N_P)*SizeOfDouble);

      fesp[0] = USpaces[k];
      fesp[1] = PSpaces[k];

      fefct[0] = U1Array[k];
      fefct[1] = U2Array[k];
      ferhs[0] = USpaces[k];
      ferhs[1] = USpaces[k];

      // find discrete form
// cout<< GALERKIN << " cout " <<TDatabase::ParamDB->DISCTYPE <<endl;
// exit(0);
      if ((mg_type==1) && (k<i+1))
          DiscreteForm = DiscreteFormUpwind;
      else
        switch(TDatabase::ParamDB->DISCTYPE)
        {
          case GALERKIN:
            DiscreteForm = DiscreteFormGalerkin;
// cout<< GALERKIN << " cout " <<TDatabase::ParamDB->DISCTYPE <<endl;
            break;
            
          case SDFEM:
            DiscreteForm = DiscreteFormSDFEM;
            break;
            
          case UPWIND:
            DiscreteForm = DiscreteFormUpwind;
            break;
            
          case SMAGORINSKY:
            DiscreteForm = DiscreteFormSmagorinsky;
            break;
            
          default:
            Error("Unknown DISCTYPE" << endl);
            return -1;
        }
      if (TDatabase::ParamDB->STOKES_PROBLEM)
        DiscreteForm = DiscreteFormUpwind;

      // initialize matrices 
      switch(TDatabase::ParamDB->NSTYPE)
      {
        case 1:
          SQMATRICES[0] = MatricesA[k];
          MATRICES[0] = MatricesB1[k];
          MATRICES[1] = MatricesB2[k];

          SQMATRICES[0]->Reset();
          MATRICES[0]->Reset();
          MATRICES[1]->Reset();

          N_SquareMatrices = 1;
          N_RectMatrices = 2;

          N_Rhs = 2;
          N_FESpaces = 2;
         break;

        case 2:
          SQMATRICES[0] = MatricesA[k];
          MATRICES[0] = MatricesB1[k];
          MATRICES[1] = MatricesB2[k];
          MATRICES[2] = MatricesB1T[k];
          MATRICES[3] = MatricesB2T[k];

          SQMATRICES[0]->Reset();
          MATRICES[0]->Reset();
          MATRICES[1]->Reset();
          MATRICES[2]->Reset();
          MATRICES[3]->Reset();

          N_SquareMatrices = 1;
          N_RectMatrices = 4;

          N_Rhs = 2;
          N_FESpaces = 2;
        break;

        case 3:
          SQMATRICES[0] = MatricesA11[k];
          SQMATRICES[1] = MatricesA12[k];
          SQMATRICES[2] = MatricesA21[k];
          SQMATRICES[3] = MatricesA22[k];
          MATRICES[0] = MatricesB1[k];
          MATRICES[1] = MatricesB2[k];

          SQMATRICES[0]->Reset();
          SQMATRICES[1]->Reset();
          SQMATRICES[2]->Reset();
          SQMATRICES[3]->Reset();
          MATRICES[0]->Reset();
          MATRICES[1]->Reset();

          N_SquareMatrices = 4;
          N_RectMatrices = 2;

          N_Rhs = 2;
          N_FESpaces = 2;
        break;

        case 4:
          SQMATRICES[0] = MatricesA11[k];
          SQMATRICES[1] = MatricesA12[k];
          SQMATRICES[2] = MatricesA21[k];
          SQMATRICES[3] = MatricesA22[k];
          MATRICES[0] = MatricesB1[k];
          MATRICES[1] = MatricesB2[k];
          MATRICES[2] = MatricesB1T[k];
          MATRICES[3] = MatricesB2T[k];

          SQMATRICES[0]->Reset();
          SQMATRICES[1]->Reset();
          SQMATRICES[2]->Reset();
          SQMATRICES[3]->Reset();
          MATRICES[0]->Reset();
          MATRICES[1]->Reset();
          MATRICES[2]->Reset();
          MATRICES[3]->Reset();

          N_SquareMatrices = 4;
          N_RectMatrices = 4;

          N_Rhs = 2;
          N_FESpaces = 2;
        break;
      }
      // get auxiliary values
      // fixed point iteration
      if (TDatabase::ParamDB->SC_NONLIN_ITE_TYPE_SADDLE==0)
        if(TDatabase::ParamDB->DISCTYPE == SMAGORINSKY)
        {
           aux =  new TAuxParam2D(NSN_FESpacesVelo_GradVelo, NSN_FctVelo_GradVelo, 
                               NSN_ParamFctVelo_GradVelo, 
                               NSN_FEValuesVelo_GradVelo, 
                               fesp, fefct, 
                               NSFctVelo_GradVelo, 
                               NSFEFctIndexVelo_GradVelo, NSFEMultiIndexVelo_GradVelo, 
                               NSN_ParamsVelo_GradVelo, NSBeginParamVelo_GradVelo);
        }
        else
        {
          aux =  new TAuxParam2D(NSN_FESpacesVelo, NSN_FctVelo, NSN_ParamFctVelo, 
                               NSN_FEValuesVelo, 
                               fesp, fefct, 
                               NSFctVelo, 
                               NSFEFctIndexVelo, NSFEMultiIndexVelo, 
                               NSN_ParamsVelo, NSBeginParamVelo);
        }
      else // Newton method
         {
           aux =  new TAuxParam2D(NSN_FESpacesVelo_GradVelo, NSN_FctVelo_GradVelo, 
                                NSN_ParamFctVelo_GradVelo, 
                                NSN_FEValuesVelo_GradVelo, 
                                fesp, fefct, 
                                NSFctVelo_GradVelo, 
                                NSFEFctIndexVelo_GradVelo, NSFEMultiIndexVelo_GradVelo, 
                                NSN_ParamsVelo_GradVelo, NSBeginParamVelo_GradVelo);
        }

    // assemble
      Assemble2D(N_FESpaces, fesp, 
               N_SquareMatrices, SQMATRICES, 
               N_RectMatrices, MATRICES, 
               N_Rhs, RHSs, ferhs,
               DiscreteForm, 
               BoundaryConditions, 
               BoundValues, 
               aux);

      if ((DiscreteForm == DiscreteFormUpwind)
          &&(!TDatabase::ParamDB->STOKES_PROBLEM))
      {
        switch(TDatabase::ParamDB->NSTYPE)
        {
          case 1:
          case 2:
            // do upwinding with one matrix
            UpwindForNavierStokes(SQMATRICES[0], U1Array[k], U2Array[k]);
            cout << "UPWINDING DONE : level " << k << endl;
            break;

          case 3:
          case 4:
            // do upwinding with two matrices
            cout << "UPWINDING DONE : level " << k << endl;
            UpwindForNavierStokes(SQMATRICES[0], U1Array[k], U2Array[k]);
            UpwindForNavierStokes(SQMATRICES[3], U1Array[k], U2Array[k]);
            break;
        } // endswitch
      } // endif

    switch((int)(TDatabase::ParamDB->REACTOR_P28+0.5))
    {

      case 123:
        CSF_Int(USpaces[k], RHSs[0], RHSs[1], h_max[k]);
      break;

      case 12345:
        CalculateIntegrals(USpaces[k], RHSs[0], RHSs[1]);
      break;


      case 54321:
        IsoIntegrals(USpaces[k], RHSs[0], RHSs[1]);
      break;
      case 6789:
        IsoInt(USpaces[k], RHSs[0], RHSs[1]);
      break;
      case 9876:
        if (((mg_type==1) && (k<i+1)) || (k == 0))
         IsoInt_spline(N_IntFace[k], IntFace_Cell[k], Edge_No[k],
                       USpaces[k], RHSs[0], RHSs[1]);
        else
         IsoInt_spline(N_IntFace[k-1], IntFace_Cell[k-1], Edge_No[k-1],
                       USpaces[k], RHSs[0], RHSs[1]);
      break; 
     default:
       cout << " Unknown interface integral type " <<endl;
       exit(-1);
      break;
    }


//    U2Array[k]->Interpolate(ExactU2);



    if (TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION >= 1)
    {
      if (TDatabase::ParamDB->NSTYPE <4)
      {
        OutPut("For slip with friction bc NSTYPE 4 is ");
        OutPut("necessary !!!!! " << endl);
        exit(4711);
      }

      // prepare everything for the assembling of slip with friction bc
      // on all levels
      N_FESpaces = 1;
      N_SquareMatrices = 4;
      N_RectMatrices = 2;
      N_Rhs = 2;
      DiscreteForm = NULL;

      SQMATRICES[0] = MatricesA11[k];
      SQMATRICES[1] = MatricesA22[k];
      SQMATRICES[2] = MatricesA12[k];
      SQMATRICES[3] = MatricesA21[k];

      MATRICES[0] = MatricesB1T[k];
      MATRICES[1] = MatricesB2T[k];

      fesp[0] = USpaces[k];
      ferhs[0] = USpaces[k];
      ferhs[1] = USpaces[k];

      RHSs[0] = RhsArray[k];
      RHSs[1] = RhsArray[k]+N_Uarray[k];
      Assemble2DSlipBC(N_FESpaces, fesp,
                       N_SquareMatrices, SQMATRICES,
                       N_RectMatrices, MATRICES,
                       N_Rhs, RHSs, ferhs,
                       DiscreteForm,
                       BoundaryConditions,
                       BoundValues,
                       aux,
                       U1Array[k],U2Array[k]);
      // reset MATRICES for solver
      SQMATRICES[0] = MatricesA11[k];
      SQMATRICES[1] = MatricesA12[k];
      SQMATRICES[2] = MatricesA21[k];
      SQMATRICES[3] = MatricesA22[k];
      MATRICES[0] = MatricesB1[k];
      MATRICES[1] = MatricesB2[k];
      MATRICES[2] = MatricesB1T[k];
      MATRICES[3] = MatricesB2T[k];
    }
    delete aux;

    // pressure separation with Neumann problem
    if (((TDatabase::ParamDB->PRESSURE_SEPARATION==3)||
         (TDatabase::ParamDB->PRESSURE_SEPARATION==4))
        &&(k==mg_level))
    {
       OutPut("assemble pressure separation problem"<<endl);
       N_SquareMatrices = 1;
       N_RectMatrices = 0;          
       N_Rhs = 1;
       SQMATRICES[0] = sqmatrixPressSep;
       SQMATRICES[0]->Reset();
       fesp[0] = pressure_separation_space;
       ferhs[0] = pressure_separation_space;
       if (TDatabase::ParamDB->PRESSURE_SEPARATION==3)
       {
          N_FESpaces = 1;
          aux = new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL); 
       }
       else
       {
          N_FESpaces = 2;
          fesp[1] = USpaces[mg_level];          
          fefct[0] = U1Array[mg_level];
          fefct[1] = U2Array[mg_level];
          aux =  new TAuxParam2D(NSN_FESpacesVelo_GradVelo, NSN_FctVelo_GradVelo, 
                                 NSN_ParamFctVelo_GradVelo, 
                                 NSN_FEValuesVelo_GradVelo, 
                                 fesp+1, fefct, 
                                 NSFctVelo_GradVelo, 
                                 NSFEFctIndexVelo_GradVelo, NSFEMultiIndexVelo_GradVelo, 
                                 NSN_ParamsVelo_GradVelo, NSBeginParamVelo_GradVelo);
       }
       RHSs[0] = rhsPressSep;
       memset(rhsPressSep, 0, N_P_sep*SizeOfDouble);      
       // assemble
       Assemble2D(N_FESpaces, fesp, 
                  N_SquareMatrices, SQMATRICES, 
                  0, NULL, 
                  N_Rhs, RHSs, ferhs,
                  DiscreteFormAuxProbPressSep, 
                  BoundaryConditionsPressureSeparation, 
                  BoundaryValuesPressureSeparation, 
                  aux);
      
       // solve linear system
//        Solver(sqmatrixPressSep, RHSs[0],separated_pressure_array,1);
       Auxprec = new TMultiGridScaIte(MatVect_Scalar, Defect_Scalar, NULL,
                                      0, N_P_sep, AuxMG, 0);
       Auxitmethod = new TFgmresIte(MatVect_Scalar, Defect_Scalar, Auxprec,
                                    0, N_P_sep, 1);
       Auxitmethod_sol = new double[N_P_sep];
       Auxitmethod_rhs = new double[N_P_sep];              
       memcpy(Auxitmethod_sol, separated_pressure_array, N_P_sep*SizeOfDouble);
       memcpy(Auxitmethod_rhs, rhsPressSep, N_P_sep*SizeOfDouble);
       SQMATRICES[0] = sqmatrixPressSep;
       Auxitmethod->Iterate(sqmatrices,NULL,Auxitmethod_sol,Auxitmethod_rhs);
       memcpy(separated_pressure_array, Auxitmethod_sol, N_P_sep*SizeOfDouble);
       memcpy(rhsPressSep, Auxitmethod_rhs, N_P_sep*SizeOfDouble);
       delete Auxitmethod_sol;
       delete Auxitmethod_rhs;
       delete Auxprec;
       delete Auxitmethod;

       // store separated pressure in the original pressure space
       Prolongate(pressure_separation_space, pressure_space,
                  separated_pressure_array, nosep_p,
		  pressure_aux_array);

	// assemble rhs for NSE
	// the gradient of the separated pressure is needed for assembling
	// this has to be said to the assembling routine by an aux object
	fesp[0] = USpaces[mg_level];
	fesp[1] = pressure_separation_space;
	
	fefct[0] = separated_pressure_fe_funct;

	aux =  new TAuxParam2D(NSN_FESpacesPressSep, NSN_FctPressSep, NSN_ParamFctPressSep, 
                               NSN_FEValuesPressSep, 
                               fesp+1, fefct, 
                               NSFctPressSep, 
                               NSFEFctIndexPressSep, NSFEMultiIndexPressSep, 
                               NSN_ParamsPressSep, NSBeginParamPressSep);
	
	// assemble the right hand side
	N_FESpaces = 2;
	N_SquareMatrices = 0;
	N_RectMatrices = 0;
	N_Rhs = 2;
	RHSs[0] = rhs;
	RHSs[1] = rhs + N_U;
	memset(rhs, 0, (2*N_U+N_P)*SizeOfDouble);
	ferhs[0] = USpaces[mg_level];
	ferhs[1] = USpaces[mg_level];
	DiscreteForm = DiscreteFormPressSep;

	Assemble2D(N_FESpaces, fesp, 
		   N_SquareMatrices, SQMATRICES, 
		   N_RectMatrices, MATRICES, 
		   N_Rhs, RHSs, ferhs,
		   DiscreteForm, 
		   BoundaryConditions, 
		   BoundValues, 
		   aux);
	
        // initialize solution array for separated pressure
        memset(sol+2*N_U,0, N_P*SizeOfDouble);
    }    
  } // endfor, assembling done
    
    // set Dirichlet nodes
    memcpy(sol+N_Active, rhs_high+N_Active, N_NonActive*SizeOfDouble);
    memcpy(sol+N_U+N_Active, rhs_high+N_U+N_Active, N_NonActive*SizeOfDouble);

    {
    int N_Hanging, N_Nodes;
    THangingNode *hn, **HangingNodes;
    HNDesc HNDescr;
    THNDesc *HNDescr_Obj;
    double *Coupling;
    int *DOF;
    int FirstHangingNodeNumber;
    double value, value1, value2;
    N_Hanging = velocity_space->GetN_Hanging();
    int i,j,k;
    if(N_Hanging)
    {
      HangingNodes = velocity_space->GetHangingNodes();
      FirstHangingNodeNumber = velocity_space->GetActiveBound();
      for(i=0;i<N_Hanging;i++)
      {
        hn = HangingNodes[i];
        HNDescr = hn->GetType();
        HNDescr_Obj = TFEDatabase2D::GetHNDesc2D(HNDescr);
        N_Nodes = HNDescr_Obj->GetN_Nodes();
        Coupling = HNDescr_Obj->GetCoeff();
        DOF = hn->GetDOF();
  
        // cout << "HN: " << FirstHangingNodeNumber+i << endl;
        
        value1 = 0;
        value2 = 0;
        for(j=0;j<N_Nodes;j++)
        {
          value = Coupling[j];
          k = DOF[j];
          // cout << j << " " << k << " " << Coupling[j] << endl;
          value1 += value * sol[k];
          value2 += value * sol[k+N_U];
        }
        sol[    FirstHangingNodeNumber+i] = value1;
        sol[N_U+FirstHangingNodeNumber+i] = value2;
      } // endfor i
    } // endif N_Hanging
    }


//     cout << "Interpolate" << endl;
//     u1->Interpolate(ExactU1);
//     u2->Interpolate(ExactU2);
//     p->Interpolate(ExactP);

    // compute defect    
    if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE) 
      IntoL20FEFunction(sol+2*N_U, N_P,PSpaces[mg_level],
                        velocity_space_code, pressure_space_code);
    defect = new double[N_Unknowns];
    memset(defect,0,N_Unknowns*SizeOfDouble);
    Defect(sqmatrices,matrices,sol,rhs_high,defect);

    if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE) 
      IntoL20Vector2D(defect+2*N_U, N_P,pressure_space_code);

    /*
    for(j=0;j<N_U;j++)
    {
      cout << setw(3) << j << setw(23) << defect[j];
      cout << setw(23) << defect[j+N_U] << endl;
    }
    for(j=0;j<N_P;j++)
    {
      cout << setw(5) << j << setw(23) << defect[j+2*N_U] << endl;
    }
    exit(-1);
    */
    
    residual =  Ddot(N_Unknowns,defect,defect);
    impuls_residual = Ddot(2*N_U,defect,defect);
    OutPut("nonlinear iteration step   0");
    OutPut(setw(14) << impuls_residual);
    OutPut(setw(14) << residual-impuls_residual);
    OutPut(setw(14) << sqrt(residual) << endl);
 

    // solve system
    switch(TDatabase::ParamDB->SOLVER_TYPE)
    {
      case AMG:
        TDatabase::ParamDB->SC_VERBOSE=1;
        TDatabase::ParamDB->CC_VERBOSE=1;     
        t1 = GetTime();
        switch(TDatabase::ParamDB->NSTYPE)
        {
          case 1:
            Solver(sqmatrixA, matrixB1, matrixB2, rhs_high, sol);
          break;

          case 2:
#ifdef __MORTAR__
            Solver(sqmatrixA, matrixB1T, matrixB2T,
                   matrixB1, matrixB2, matrix_mortar, rhs_high, sol);
#else
            Solver(sqmatrixA, matrixB1T, matrixB2T, 
                   matrixB1, matrixB2, rhs_high, sol);
#endif
          break;

          case 3:
            Solver(sqmatrixA11, sqmatrixA12, sqmatrixA21, 
                   sqmatrixA22, matrixB1, matrixB2, rhs_high, sol);
          break;

          case 4:
            DirectSolver(sqmatrixA11, sqmatrixA12, sqmatrixA21, 
                   sqmatrixA22, matrixB1T, matrixB2T, 
                   matrixB1, matrixB2, rhs_high, sol);
          break;
        }
        t2 = GetTime();
        solver_time += (t2-t1);
      break;

      case GMG:
        t1 = GetTime();
        switch (TDatabase::ParamDB->SC_SOLVER_SADDLE)
        {
          case 11:
            zerostart = 1;
            break;
          case 16:
            zerostart = 0;
            break;
        }
        // build preconditioner
        switch (TDatabase::ParamDB->SC_PRECONDITIONER_SADDLE)
        {
          case 5:
            prec = new TMultiGridIte(MatVect, Defect, NULL,
                                     0, N_Unknowns, MG, zerostart);
            break;
          default:
            OutPut("Unknown preconditioner !!!" << endl);
            exit(4711);
        }
        switch (TDatabase::ParamDB->SC_SOLVER_SADDLE)
        {
          case 11:
            itmethod = new TFixedPointIte(MatVect, Defect, prec,
                                          0, N_Unknowns, 0);
            if (TDatabase::ParamDB->SC_PRECONDITIONER_SADDLE == 5)
            {
              itmethod_sol = new double[N_Unknowns];
              itmethod_rhs = new double[N_Unknowns];              
              memcpy(itmethod_sol, sol, N_Unknowns*SizeOfDouble);
              memcpy(itmethod_rhs, rhs_high, N_Unknowns*SizeOfDouble);
            }
            else
            {
              itmethod_sol = sol;
              itmethod_rhs = rhs_high;            
            }
            break;
          case 16:
            itmethod = new TFgmresIte(MatVect, Defect, prec,
                                      0, N_Unknowns, 0);
            if (TDatabase::ParamDB->SC_PRECONDITIONER_SADDLE == 5)
            {
              itmethod_sol = new double[N_Unknowns];
              itmethod_rhs = new double[N_Unknowns];              
              memcpy(itmethod_sol, sol, N_Unknowns*SizeOfDouble);
              memcpy(itmethod_rhs, rhs_high, N_Unknowns*SizeOfDouble);
            }
            else
            {
              itmethod_sol = sol;
              itmethod_rhs = rhs_high;
            }
            break;
          default:
            OutPut("Unknown solver !!!" << endl);
            exit(4711);
        }
        t2 = GetTime();
        solver_time += (t2-t1);
        for (ll=0;ll<calculations;ll++)
        {
          if(TDatabase::ParamDB->P9 == 123456789)
          {
            fesp[0] = USpaces[mg_level];
            fefct[0] = U1Array[mg_level];
            fefct[1] = U2Array[mg_level];
            aux =  new TAuxParam2D(NSN_FESpacesVelo, NSN_FctVelo, 
                                   NSN_ParamFctVelo, 
                                   NSN_FEValuesVelo, 
                                   fesp, fefct, 
                                   NSFctVelo, 
                                   NSFEFctIndexVelo, NSFEMultiIndexVelo, 
                                   NSN_ParamsVelo, NSBeginParamVelo);
            
            // errors in first velocity component
            U1Array[mg_level]->GetErrors(ExactU1, 3, NSAllDerivatives, 2, 
                                         L2H1Errors, 
                                         NULL, aux, 1, USpaces+mg_level, errors_mg);
            p1 = errors_mg[0] * errors_mg[0];
            
            // errors in second velocity component
            U2Array[mg_level]->GetErrors(ExactU2, 3, NSAllDerivatives, 2, 
                                         L2H1Errors, 
                                         NULL, aux, 1, USpaces+mg_level, errors_mg);
            p1 += errors_mg[0] * errors_mg[0];
            // errors in first velocity component
//          soldiff_fe1->GetErrors(ExactNull, 3, NSAllDerivatives, 2, 
            //                               L2H1Errors, 
            //                     NULL, aux, 1, USpaces+mg_level, errors_mg);
            p1 = errors_mg[0] * errors_mg[0];

            // errors in second velocity component
            //soldiff_fe2->GetErrors(ExactNull, 3, NSAllDerivatives, 2,
            //                     L2H1Errors, 
            //                     NULL, aux, 1, USpaces+mg_level, errors_mg);
            p1 += errors_mg[0] * errors_mg[0];
            firsterror =  p1 = sqrt(p1);
            delete aux;
            
            for(l=0;l<N_Unknowns;l++)          
              soldiff[l] = sol[l]-fesol[l];
            
            p3 = sqrt(Ddot(2*N_U,soldiff,soldiff));
            firsterrorl2 = p3;          
            
            OutPut("iteration -1  L2(u): " <<  p1 << " l2(u) " << p3 << endl);
            p2 = p1;
            p4 = p3;
          } // endif MEASURE_ERRORS

          t1 = GetTime();

//        for (k=0;k<N_U;k++)
//         cout << "u1 " << itmethod_rhs[k]<< " u2 " << itmethod_rhs[k+N_U]<<endl;
// cout<<endl;
          // solve linear system
          N_LinIter+=itmethod->Iterate(sqmatrices,matrices,itmethod_sol,itmethod_rhs);

//        for (k=0;k<N_U;k++)
//         cout << "u1 " << itmethod_sol[k]<< " u2 " << itmethod_sol[k+N_U]<<endl;
// cout<<endl;
 /*
          residual = 1;
          while(residual>1e-10)
          {
            N_LinIter+=prec->Iterate(sqmatrices,matrices,itmethod_sol,itmethod_rhs);
            Defect(sqmatrices,matrices,sol,rhs_high,defect);

            if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE) 
              IntoL20Vector2D(defect+2*N_U, N_P,pressure_space_code);
            residual =  Ddot(N_Unknowns,defect,defect);
            cout << N_LinIter << " res: " << residual << " " << coll->GetN_Cells() << endl;
          }
          */
          switch (TDatabase::ParamDB->SC_SOLVER_SADDLE)
          {
            case 11:
              if (TDatabase::ParamDB->SC_PRECONDITIONER_SADDLE == 5)
              {
                memcpy(sol, itmethod_sol, N_Unknowns*SizeOfDouble);
                memcpy(rhs_high, itmethod_rhs, N_Unknowns*SizeOfDouble);
              }
              break;
            case 16:
              if (TDatabase::ParamDB->SC_PRECONDITIONER_SADDLE == 5)
              {
                memcpy(sol, itmethod_sol, N_Unknowns*SizeOfDouble);
                memcpy(rhs_high, itmethod_rhs, N_Unknowns*SizeOfDouble);
              }
              break;
          }
          t2 = GetTime();
          solver_time += (t2-t1);

          if(TDatabase::ParamDB->P9 == 123456789)
          {
            fesp[0] = USpaces[mg_level];
            fefct[0] = U1Array[mg_level];
            fefct[1] = U2Array[mg_level];
            aux =  new TAuxParam2D(NSN_FESpacesVelo, NSN_FctVelo, 
                                   NSN_ParamFctVelo, 
                                   NSN_FEValuesVelo, 
                                   fesp, fefct, 
                                   NSFctVelo, 
                                   NSFEFctIndexVelo, NSFEMultiIndexVelo, 
                                   NSN_ParamsVelo, NSBeginParamVelo);
            
            // errors in first velocity component
            U1Array[mg_level]->GetErrors(ExactU1, 3, NSAllDerivatives, 2,
                                         L2H1Errors, 
                                         NULL, aux, 1, USpaces+mg_level, errors_mg);
            p1 = errors_mg[0] * errors_mg[0];
            
            // errors in second velocity component
            U2Array[mg_level]->GetErrors(ExactU2, 3, NSAllDerivatives, 2, 
                                         L2H1Errors, 
                                         NULL, aux, 1, USpaces+mg_level, errors_mg);
            // errors in first velocity component
            /* soldiff_fe1->GetErrors(ExactNull, 3, NSAllDerivatives, 2, 
               L2H1Errors, 
               NULL, aux, 1, USpaces+mg_level, errors_mg);
               p1 = errors_mg[0] * errors_mg[0];
               
               // errors in second velocity component
               soldiff_fe2->GetErrors(ExactNull, 3, NSAllDerivatives, 2, 
               L2H1Errors, 
               NULL, aux, 1, USpaces+mg_level, errors_mg);
            */ 
            p1 += errors_mg[0] * errors_mg[0];
            p1 = sqrt(p1);
            delete aux;
            
            for(l=0;l<N_Unknowns;l++)          
              soldiff[l] = sol[l]-fesol[l];          
            p3 = sqrt(Ddot(2*N_U,soldiff,soldiff));

            OutPut("iteration " << j << " L2(u): " <<  p1 << " redu ");
            OutPut(" rateL2 " << pow(p1/firsterror,1.0/N_LinIter));
            OutPut(p1/p2 << " l2(u) " << p3 << " redu " << p3/p4);
            OutPut(" ratel2 " << pow(p3/firsterrorl2,1.0/N_LinIter)<< endl);
            
            p2 = p1;          
            p4 = p3;          
            lasterror = p1;
            lasterrorl2 = p3;
            // if (res/res2 > convergence_speed)
            // {
            //  OutPut("SLOW !!!!!!!!! " << endl);
            //  break;
            // }
          } // endif MEASURE_ERRORS

          for(l=0;l<N_Unknowns;l++)
          {
            p2 = sol[l]-oldsol[l];
            sol[l] = oldsol[l] + omega * p2;
          }
          
          if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE) 
            IntoL20FEFunction(sol+2*N_U, N_P,PSpaces[mg_level], 
                              velocity_space_code, pressure_space_code);
          
          if(TDatabase::ParamDB->P9 == 123456789)
          {
            if (!pre_calculation)
            {
              OutPut("average error reduction rate (L2/l2) " << pow(lasterror/firsterror,1.0/N_LinIter));
              OutPut(" " << pow(lasterrorl2/firsterrorl2,1.0/N_LinIter) << endl);
            }
            if (pre_calculation)
            {
              for (k=0;k<2*N_U+N_P;k++)
              {
                fesol[k] = sol[k];
                sol[k] = 0.0;
              }     
              pre_calculation = 0;  
              N_LinIter=0;
              TDatabase::ParamDB->SC_LIN_RES_NORM_MIN_SADDLE = 1e-10;
            }
          }
        }
        break;
    }
    
//        for (k=0;k<N_U;k++)
//         cout << "u1 " << sol[k]<< " u2 " << sol[k+N_U]<<endl;
// exit(0);
// ************************************************************* //
// end of first nonlinear step                                   
// ************************************************************* //
    // don't know why it does not work
    // if (TDatabase::ParamDB->PRESSURE_SEPARATION<3)
    {
       OutPut("MEMORY: " << setw(10) << GetMemory()/(1048576.0));
       OutPut(" MB" << endl);
    }
// ************************************************************* //
// the nonlinear iteration 
// ************************************************************* //

    for(j=1;j<=Max_It;j++)
    {
      if (TDatabase::ParamDB->SOLVER_TYPE==GMG)
        MG->RestrictToAllGrids();

      memcpy(oldsol, sol, SizeOfDouble*N_Unknowns);
  
      for(k=low;k<=mg_level;k++)
      {
        fesp[0] = USpaces[k];
        fesp[1] = PSpaces[k];
        
        fefct[0] = U1Array[k];
        fefct[1] = U2Array[k];

        if ((k<i+1)&&(mg_type==1))
          DiscreteForm = DiscreteFormNLUpwind;
        else
          switch(TDatabase::ParamDB->DISCTYPE)
          {
            case GALERKIN:
              DiscreteForm = DiscreteFormNLGalerkin;
              break;
              
            case SDFEM:
              DiscreteForm = DiscreteFormNLSDFEM;
              break;
              
            case UPWIND:
              DiscreteForm = DiscreteFormNLUpwind;
              break;
              
            case SMAGORINSKY:
              DiscreteForm = DiscreteFormNLSmagorinsky;
              break;
          } // endswitch
        // this can only happen if pressure separation is applied
        if (TDatabase::ParamDB->STOKES_PROBLEM)
           DiscreteForm = DiscreteFormNLUpwind;
        
        switch(TDatabase::ParamDB->NSTYPE)
        {
          case 1:
            SQMATRICES[0] = MatricesA[k];
            SQMATRICES[0]->Reset();

            N_SquareMatrices = 1;
            N_RectMatrices = 0;

            N_Rhs = 0;
            N_FESpaces = 1;
          break;

          case 2:
            SQMATRICES[0] = MatricesA[k];
            SQMATRICES[0]->Reset();

            N_SquareMatrices = 1;
            if (DiscreteForm == DiscreteFormNLSDFEM)
            {
              N_RectMatrices = 2;
              MATRICES[0] = MatricesB1T[k];
              MATRICES[1] = MatricesB2T[k];

              MATRICES[0]->Reset();
              MATRICES[1]->Reset();

              N_Rhs = 2;
              rhs = RhsArray[k];
              RHSs[0] = rhs;
              RHSs[1] = rhs + N_Uarray[k];
              memset(rhs, 0, (2*N_Uarray[k]+N_Parray[k])*SizeOfDouble);
              N_FESpaces = 2;
              ferhs[0] = USpaces[k];
              ferhs[1] = USpaces[k];
            }
            else
            {
              N_RectMatrices = 0;

              N_Rhs = 0;
              N_FESpaces = 1;
            }
          break;

          case 3:
            if (TDatabase::ParamDB->SC_NONLIN_ITE_TYPE_SADDLE==0)
            {
              if(DiscreteForm == DiscreteFormNLSmagorinsky)
              {
                SQMATRICES[0] = MatricesA11[k];
                SQMATRICES[1] = MatricesA12[k];
                SQMATRICES[2] = MatricesA21[k];
                SQMATRICES[3] = MatricesA22[k];
                SQMATRICES[0]->Reset();
                SQMATRICES[1]->Reset();
                SQMATRICES[2]->Reset();
                SQMATRICES[3]->Reset();
                
                N_SquareMatrices = 4;
                N_RectMatrices = 0;
                last_sq = 3;
              }
              else
              {
                SQMATRICES[0] = MatricesA11[k];
                SQMATRICES[1] = MatricesA22[k];
                SQMATRICES[0]->Reset();
                SQMATRICES[1]->Reset();
                
                N_SquareMatrices = 2;
                N_RectMatrices = 0;
                last_sq = 1;
              }
              N_Rhs = 0;
              N_FESpaces = 1;
            }
            else // Newton method
            {
              SQMATRICES[0] = MatricesA11[k];
              SQMATRICES[1] = MatricesA12[k];
              SQMATRICES[2] = MatricesA21[k];
              SQMATRICES[3] = MatricesA22[k];
              SQMATRICES[0]->Reset();
              SQMATRICES[1]->Reset();
              SQMATRICES[2]->Reset();
              SQMATRICES[3]->Reset();
              
              N_SquareMatrices = 4;
              N_RectMatrices = 0;

              N_Rhs = 2;
              N_FESpaces = 1;             
              rhs = RhsArray[k];
              RHSs[0] = rhs;
              RHSs[1] = rhs + N_Uarray[k];
              memset(rhs, 0, (2*N_Uarray[k]+N_Parray[k])*SizeOfDouble);
              ferhs[0] = USpaces[k];
              ferhs[1] = USpaces[k];
              last_sq = 3;
            }
          break;

          case 4:
            if (TDatabase::ParamDB->SC_NONLIN_ITE_TYPE_SADDLE==0)
            {
              if (DiscreteForm == DiscreteFormNLSDFEM)
              {
                N_SquareMatrices = 2;
                SQMATRICES[0] = MatricesA11[k];
                SQMATRICES[1] = MatricesA22[k];
                SQMATRICES[0]->Reset();
                SQMATRICES[1]->Reset();
                
                N_RectMatrices = 2;
                MATRICES[0] = MatricesB1T[k];
                MATRICES[1] = MatricesB2T[k];
                MATRICES[0]->Reset();
                MATRICES[1]->Reset();
                
                N_Rhs = 2;
                rhs = RhsArray[k];
                RHSs[0] = rhs;
                RHSs[1] = rhs + N_Uarray[k];
                memset(rhs, 0, (2*N_Uarray[k]+N_Parray[k])*SizeOfDouble);
                N_FESpaces = 2;
                
                ferhs[0] = USpaces[k];
                ferhs[1] = USpaces[k];
                last_sq = 1;
              }
              else
              {
                if (DiscreteForm == DiscreteFormNLSmagorinsky)
                {
                  N_RectMatrices = 0;
                  
                  N_SquareMatrices = 4;
                  SQMATRICES[0] = MatricesA11[k];
                  SQMATRICES[1] = MatricesA12[k];
                  SQMATRICES[2] = MatricesA21[k];
                  SQMATRICES[3] = MatricesA22[k];
                  SQMATRICES[0]->Reset();
                  SQMATRICES[1]->Reset();
                  SQMATRICES[2]->Reset();
                  SQMATRICES[3]->Reset();
                  
                  N_Rhs = 0;
                  N_FESpaces = 1;
                  last_sq = 3;
                }
                else // default
                {
                  N_SquareMatrices = 2;
                  SQMATRICES[0] = MatricesA11[k];
                  SQMATRICES[1] = MatricesA22[k];
                  SQMATRICES[0]->Reset();
                  SQMATRICES[1]->Reset();
                  
                  N_RectMatrices = 0;
                  
                  N_Rhs = 0;
                  N_FESpaces = 1;
                  last_sq = 1;
                }
              }
            }
            else // Newton method
            {
              SQMATRICES[0] = MatricesA11[k];
              SQMATRICES[1] = MatricesA12[k];
              SQMATRICES[2] = MatricesA21[k];
              SQMATRICES[3] = MatricesA22[k];
              SQMATRICES[0]->Reset();
              SQMATRICES[1]->Reset();
              SQMATRICES[2]->Reset();
              SQMATRICES[3]->Reset();
              
              N_SquareMatrices = 4;
              N_RectMatrices = 0;

              N_Rhs = 2;
              N_FESpaces = 1;             
              rhs = RhsArray[k];
              RHSs[0] = rhs;
              RHSs[1] = rhs + N_Uarray[k];
              memset(rhs, 0, (2*N_Uarray[k]+N_Parray[k])*SizeOfDouble);
              ferhs[0] = USpaces[k];
              ferhs[1] = USpaces[k];
              last_sq = 3;
              
              if (DiscreteForm == DiscreteFormNLSDFEM)
              {                
                N_RectMatrices = 2;
                MATRICES[0] = MatricesB1T[k];
                MATRICES[1] = MatricesB2T[k];
                MATRICES[0]->Reset();
                MATRICES[1]->Reset();
                N_FESpaces = 2;
               }              
            }
          break;
        } // endswitch

        if (TDatabase::ParamDB->SC_NONLIN_ITE_TYPE_SADDLE==0)
          if (DiscreteForm == DiscreteFormNLSmagorinsky)
          {
            aux =  new TAuxParam2D(NSN_FESpacesVelo_GradVelo, NSN_FctVelo_GradVelo, 
                                 NSN_ParamFctVelo_GradVelo, 
                                 NSN_FEValuesVelo_GradVelo, 
                                 fesp, fefct, 
                                 NSFctVelo_GradVelo, 
                                 NSFEFctIndexVelo_GradVelo, NSFEMultiIndexVelo_GradVelo, 
                                 NSN_ParamsVelo_GradVelo, NSBeginParamVelo_GradVelo);
          }
          else
          {
            aux =  new TAuxParam2D(NSN_FESpacesVelo, NSN_FctVelo, NSN_ParamFctVelo, 
                                 NSN_FEValuesVelo, 
                                 fesp, fefct, 
                                 NSFctVelo, 
                                 NSFEFctIndexVelo, NSFEMultiIndexVelo, 
                                 NSN_ParamsVelo, NSBeginParamVelo);
          }
        else // Newton method
          {
            aux =  new TAuxParam2D(NSN_FESpacesVelo_GradVelo, NSN_FctVelo_GradVelo, 
                                 NSN_ParamFctVelo_GradVelo, 
                                 NSN_FEValuesVelo_GradVelo, 
                                 fesp, fefct, 
                                 NSFctVelo_GradVelo, 
                                 NSFEFctIndexVelo_GradVelo, NSFEMultiIndexVelo_GradVelo, 
                                 NSN_ParamsVelo_GradVelo, NSBeginParamVelo_GradVelo);
          }
        // assembling
        Assemble2D(N_FESpaces, fesp,
                 N_SquareMatrices, SQMATRICES,
                 N_RectMatrices, MATRICES,
                 N_Rhs, RHSs, ferhs,
                 DiscreteForm,
                 BoundaryConditions, 
                 BoundValues, 
                 aux);

        if ((DiscreteForm == DiscreteFormNLUpwind)
            &&(!TDatabase::ParamDB->STOKES_PROBLEM))
        {
          switch(TDatabase::ParamDB->NSTYPE)
          {
            case 1:
            case 2:
              // do upwinding with one matrix
              UpwindForNavierStokes(SQMATRICES[0], U1Array[k], U2Array[k]);
              cout << "UPWINDING DONE : level " << k << endl;
            break;

            case 3:
            case 4:
              // do upwinding with two matrices
              UpwindForNavierStokes(SQMATRICES[0], U1Array[k], U2Array[k]);
              UpwindForNavierStokes(SQMATRICES[last_sq], U1Array[k], U2Array[k]);
              cout << "UPWINDING DONE(2) : level " << k << endl;
            break;
          } // endswitch
        } // endif

        if (TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION >= 1)
        {
       
          // prepare everything for the assembling of slip with friction bc
          // on all levels
          N_FESpaces = 1;
          N_SquareMatrices = 4;
          N_RectMatrices = 0;
          N_Rhs = 2;
          DiscreteForm = NULL;
          
          SQMATRICES[0] = MatricesA11[k];
          SQMATRICES[1] = MatricesA22[k];
          SQMATRICES[2] = MatricesA12[k];
          SQMATRICES[3] = MatricesA21[k];
      
          fesp[0] = USpaces[k];
          ferhs[0] = USpaces[k];
          ferhs[1] = USpaces[k];
          
          RHSs[0] = RhsArray[k];
          RHSs[1] = RhsArray[k]+N_Uarray[k];
          
          Assemble2DSlipBC(N_FESpaces, fesp, 
                           N_SquareMatrices, SQMATRICES, 
                           N_RectMatrices, MATRICES, 
                           N_Rhs, RHSs, ferhs,
                           DiscreteForm, 
                           BoundaryConditions, 
                           BoundValues, 
                           aux,
                           U1Array[k],U2Array[k]);
        }
        delete aux;

      } // endfor k, assembling done

      if ((TDatabase::ParamDB->PRESSURE_SEPARATION==2)&&(j==1))
      {
        OutPut("apply pressure separation"<<endl);
        // save original pressure
        if ((j==1)||(1))
           memcpy(nosep_p,sol+2*N_U, N_P*SizeOfDouble);
        else
        {
           Daxpy(N_P, 1.0, sol+2*N_U, nosep_p);
           memcpy(sol+2*N_U, nosep_p, N_P*SizeOfDouble);
        }
	// assemble separated pressure entry on right hand side
	// first : interpolate discrete normal pressure to the
	//         pressure separation space
	Prolongate(pressure_space, pressure_separation_space,
		   sol+2*N_U, separated_pressure_array,
		   separated_pressure_aux_array);
        
	// second : assemble 
	// the gradient of the separated pressure is needed for assembling
	// this has to be said to the assembling routine by an aux object
	fesp[0] = USpaces[mg_level];
	fesp[1] = pressure_separation_space;
	
	fefct[0] = separated_pressure_fe_funct;

	aux =  new TAuxParam2D(NSN_FESpacesPressSep, NSN_FctPressSep, NSN_ParamFctPressSep, 
                               NSN_FEValuesPressSep, 
                               fesp+1, fefct, 
                               NSFctPressSep, 
                               NSFEFctIndexPressSep, NSFEMultiIndexPressSep, 
                               NSN_ParamsPressSep, NSBeginParamPressSep);
	
	// assemble the right hand side
	N_FESpaces = 2;
	N_SquareMatrices = 0;
	N_RectMatrices = 0;
	N_Rhs = 2;
	RHSs[0] = rhs;
	RHSs[1] = rhs + N_U;
	memset(rhs, 0, (2*N_U+N_P)*SizeOfDouble);
	ferhs[0] = USpaces[mg_level];
	ferhs[1] = USpaces[mg_level];
	DiscreteForm = DiscreteFormPressSep;

	Assemble2D(N_FESpaces, fesp, 
		   N_SquareMatrices, SQMATRICES, 
		   N_RectMatrices, MATRICES, 
		   N_Rhs, RHSs, ferhs,
		   DiscreteForm, 
		   BoundaryConditions, 
		   BoundValues, 
		   aux);
	
        // initialize solution array for separated pressure
        memset(sol+2*N_U,0, N_P*SizeOfDouble);
      }

      if ((TDatabase::ParamDB->PRESSURE_SEPARATION==5)&&(j==1))
      {
         OutPut("assemble pressure separation problem"<<endl);
         N_SquareMatrices = 1;
         N_RectMatrices = 0;          
         N_Rhs = 1;
         SQMATRICES[0] = sqmatrixPressSep;
         SQMATRICES[0]->Reset();
         fesp[0] = pressure_separation_space;
         ferhs[0] = pressure_separation_space;
         N_FESpaces = 2;
         fesp[1] = USpaces[mg_level];          
         fefct[0] = U1Array[mg_level];
         fefct[1] = U2Array[mg_level];
         aux =  new TAuxParam2D(NSN_FESpacesVelo_GradVelo, NSN_FctVelo_GradVelo, 
                                NSN_ParamFctVelo_GradVelo, 
                                NSN_FEValuesVelo_GradVelo, 
                                fesp+1, fefct, 
                                NSFctVelo_GradVelo, 
                                NSFEFctIndexVelo_GradVelo, NSFEMultiIndexVelo_GradVelo, 
                                NSN_ParamsVelo_GradVelo, NSBeginParamVelo_GradVelo);
       
         RHSs[0] = rhsPressSep;
         memset(rhsPressSep, 0, N_P_sep*SizeOfDouble);      
         // assemble
         Assemble2D(N_FESpaces, fesp, 
                    N_SquareMatrices, SQMATRICES, 
                    0, NULL, 
                    N_Rhs, RHSs, ferhs,
                    DiscreteFormAuxProbPressSep, 
                    BoundaryConditionsPressureSeparation, 
                    BoundaryValuesPressureSeparation, 
                    aux);
         
         // solve linear system

         Auxprec = new TMultiGridScaIte(MatVect_Scalar, Defect_Scalar, NULL,
                                        0, N_P_sep, AuxMG, 0);
         Auxitmethod = new TFgmresIte(MatVect_Scalar, Defect_Scalar, Auxprec,
                                      0, N_P_sep, 1);
         Auxitmethod_sol = new double[N_P_sep];
         Auxitmethod_rhs = new double[N_P_sep];              
         memcpy(Auxitmethod_sol, separated_pressure_array, N_P_sep*SizeOfDouble);
         memcpy(Auxitmethod_rhs, rhsPressSep, N_P_sep*SizeOfDouble);
         SQMATRICES[0] = sqmatrixPressSep;
         Auxitmethod->Iterate(sqmatrices,NULL,Auxitmethod_sol,Auxitmethod_rhs);
         memcpy(separated_pressure_array, Auxitmethod_sol, N_P_sep*SizeOfDouble);
         memcpy(rhsPressSep, Auxitmethod_rhs, N_P_sep*SizeOfDouble);
         delete Auxitmethod_sol;
         delete Auxitmethod_rhs;
         delete Auxprec;
         delete Auxitmethod;

         // store separated pressure in the original pressure space
         Prolongate(pressure_separation_space, pressure_space,
                    separated_pressure_array, nosep_p,
                    pressure_aux_array);
         
         // assemble rhs for NSE
         // the gradient of the separated pressure is needed for assembling
         // this has to be said to the assembling routine by an aux object
         fesp[0] = USpaces[mg_level];
         fesp[1] = pressure_separation_space;
         
         fefct[0] = separated_pressure_fe_funct;
         
         aux =  new TAuxParam2D(NSN_FESpacesPressSep, NSN_FctPressSep, NSN_ParamFctPressSep, 
                                NSN_FEValuesPressSep, 
                                fesp+1, fefct, 
                                NSFctPressSep, 
                               NSFEFctIndexPressSep, NSFEMultiIndexPressSep, 
                                NSN_ParamsPressSep, NSBeginParamPressSep);
         
         // assemble the right hand side
         N_FESpaces = 2;
         N_SquareMatrices = 0;
         N_RectMatrices = 0;
         N_Rhs = 2;
         RHSs[0] = rhs;
         RHSs[1] = rhs + N_U;
         memset(rhs, 0, (2*N_U+N_P)*SizeOfDouble);
         ferhs[0] = USpaces[mg_level];
         ferhs[1] = USpaces[mg_level];
         DiscreteForm = DiscreteFormPressSep;

         Assemble2D(N_FESpaces, fesp, 
                    N_SquareMatrices, SQMATRICES, 
                    N_RectMatrices, MATRICES, 
                    N_Rhs, RHSs, ferhs,
                    DiscreteForm, 
                    BoundaryConditions, 
                    BoundValues, 
                    aux);
         
         // initialize solution array for separated pressure
         memset(sol+2*N_U,0, N_P*SizeOfDouble);         
      }

      // end of assembling

       
      // reset MATRICES for solver
      switch(TDatabase::ParamDB->NSTYPE)
      {
        case 1:
          SQMATRICES[0] = MatricesA[mg_level];
          MATRICES[0] = MatricesB1[mg_level];
          MATRICES[1] = MatricesB2[mg_level];
          break;
        case 2:
          SQMATRICES[0] = MatricesA[mg_level];
          MATRICES[0] = MatricesB1[mg_level];
          MATRICES[1] = MatricesB2[mg_level];
          MATRICES[2] = MatricesB1T[mg_level];
          MATRICES[3] = MatricesB2T[mg_level];
          break;
        case 3:
          SQMATRICES[0] = MatricesA11[mg_level];
          SQMATRICES[1] = MatricesA12[mg_level];
          SQMATRICES[2] = MatricesA21[mg_level];
          SQMATRICES[3] = MatricesA22[mg_level];
          MATRICES[0] = MatricesB1[mg_level];
          MATRICES[1] = MatricesB2[mg_level];
          break;
        case 4:
          SQMATRICES[0] = MatricesA11[mg_level];
          SQMATRICES[1] = MatricesA12[mg_level];
          SQMATRICES[2] = MatricesA21[mg_level];
          SQMATRICES[3] = MatricesA22[mg_level];
          MATRICES[0] = MatricesB1[mg_level];
          MATRICES[1] = MatricesB2[mg_level];
          MATRICES[2] = MatricesB1T[mg_level];
          MATRICES[3] = MatricesB2T[mg_level];
          break;
      }

      // set rhs for Dirichlet nodes
      memcpy(sol+N_Active, rhs_high+N_Active, N_NonActive*SizeOfDouble);
      memcpy(sol+N_U+N_Active, rhs_high+N_U+N_Active, N_NonActive*SizeOfDouble);

      {
      int N_Hanging, N_Nodes;
      THangingNode *hn, **HangingNodes;
      HNDesc HNDescr;
      THNDesc *HNDescr_Obj;
      double *Coupling;
      int *DOF;
      int FirstHangingNodeNumber;
      double value, value1, value2;
      N_Hanging = velocity_space->GetN_Hanging();
      int i,j,k;
      if(N_Hanging)
      {
        HangingNodes = velocity_space->GetHangingNodes();
        FirstHangingNodeNumber = velocity_space->GetActiveBound();
        for(i=0;i<N_Hanging;i++)
        {
          hn = HangingNodes[i];
          HNDescr = hn->GetType();
          HNDescr_Obj = TFEDatabase2D::GetHNDesc2D(HNDescr);
          N_Nodes = HNDescr_Obj->GetN_Nodes();
          Coupling = HNDescr_Obj->GetCoeff();
          DOF = hn->GetDOF();
    
          // cout << "HN: " << FirstHangingNodeNumber+i << endl;
          
          value1 = 0;
          value2 = 0;
          for(j=0;j<N_Nodes;j++)
          {
            value = Coupling[j];
            k = DOF[j];
            // cout << j << " " << k << " " << Coupling[j] << endl;
            value1 += value * sol[k];
            value2 += value * sol[k+N_U];
          }
          sol[    FirstHangingNodeNumber+i] = value1;
          sol[N_U+FirstHangingNodeNumber+i] = value2;
        } // endfor i
      } // endif N_Hanging
      }

      // compute defect
      memset(defect,0,N_Unknowns*SizeOfDouble);
      Defect(sqmatrices,matrices,sol,rhs_high,defect);
      if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE) 
        IntoL20Vector2D(defect+2*N_U, N_P,pressure_space_code);
      residual =  Ddot(N_Unknowns,defect,defect);
      impuls_residual = Ddot(2*N_U,defect,defect);
      if (pressure_separation == 1)
  	OutPut("p sep ");
      OutPut("nonlinear iteration step " << setw(3) << j);
      OutPut(setw(14) << impuls_residual);
      OutPut(setw(14) << residual-impuls_residual);
      OutPut(setw(14) << sqrt(residual) << endl);
     
      // check if convergence is too slow
      last_digit_ite = j%10;
      if (sqrt(residual)>convergence_speed*residuals[last_digit_ite])
        slow_conv=1;
      residuals[last_digit_ite]=sqrt(residual);
      
      // check if convergence criterion is fulfilled 
      if ((sqrt(residual)<=limit)||(j==Max_It)||(slow_conv)) 
      {
         if (pressure_separation == 1)
          OutPut("P SEP ");
        OutPut("ITE : " << setw(3) << j);
        OutPut(" (" << setw(3) << N_LinIter << " LINITE)");
        OutPut("  TIME FOR SOLVER : " << solver_time << "s");
        OutPut("  RES : " <<  sqrt(residual));
        if (slow_conv)
          OutPut(" SLOW !!!");
        OutPut(endl);
  	if (TDatabase::ParamDB->PRESSURE_SEPARATION==0) 	    
	    break;
	if (((TDatabase::ParamDB->PRESSURE_SEPARATION==1)||
            (TDatabase::ParamDB->PRESSURE_SEPARATION==2))&&
            (pressure_separation == 1))
        {
          // compute pressure
           memcpy(sol+2*N_U, nosep_p, N_P*SizeOfDouble);  
           //Daxpy(N_P,1.0,nosep_p,sol+2*N_U);
          break;
        }
 	if ((TDatabase::ParamDB->PRESSURE_SEPARATION==3)||
            (TDatabase::ParamDB->PRESSURE_SEPARATION==4)||
            (TDatabase::ParamDB->PRESSURE_SEPARATION==5))
        {
           Daxpy(N_P,1.0,nosep_p,sol+2*N_U);
           break;
        }
       // save original pressure
        memcpy(nosep_p,sol+2*N_U, N_P*SizeOfDouble);
        
	// assemble separated pressure entry on right hand side
	// first : interpolate discrete normal pressure to the
	//         pressure separation space
	Prolongate(pressure_space, pressure_separation_space,
		   sol+2*N_U, separated_pressure_array,
		   separated_pressure_aux_array);
        
        // TESTING: USE INTERPOLATION
        // separated_pressure_fe_funct->Interpolate(ExactP);
	// second : assemble 
	// the gradient of the separated pressure is needed for assembling
	// this has to be said to the assembling routine by an aux object
	fesp[0] = USpaces[mg_level];
	fesp[1] = pressure_separation_space;
	
	fefct[0] = separated_pressure_fe_funct;

	aux =  new TAuxParam2D(NSN_FESpacesPressSep, NSN_FctPressSep, NSN_ParamFctPressSep, 
                               NSN_FEValuesPressSep, 
                               fesp+1, fefct, 
                               NSFctPressSep, 
                               NSFEFctIndexPressSep, NSFEMultiIndexPressSep, 
                               NSN_ParamsPressSep, NSBeginParamPressSep);
	
	// assemble the right hand side
	N_FESpaces = 2;
	N_SquareMatrices = 0;
	N_RectMatrices = 0;
	N_Rhs = 2;
	RHSs[0] = rhs;
	RHSs[1] = rhs + N_U;
	memset(rhs, 0, (2*N_U+N_P)*SizeOfDouble);
	ferhs[0] = USpaces[mg_level];
	ferhs[1] = USpaces[mg_level];
	DiscreteForm = DiscreteFormPressSep;

	Assemble2D(N_FESpaces, fesp, 
		   N_SquareMatrices, SQMATRICES, 
		   N_RectMatrices, MATRICES, 
		   N_Rhs, RHSs, ferhs,
		   DiscreteForm, 
		   BoundaryConditions, 
		   BoundValues, 
		   aux);
	
        // initialize solution array for separated pressure
        memset(sol+2*N_U,0, N_P*SizeOfDouble);
	// compute defect
	memset(defect,0,N_Unknowns*SizeOfDouble);
	Defect(sqmatrices,matrices,sol,rhs_high,defect);
	if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE) 
	    IntoL20Vector2D(defect+2*N_U, N_P,pressure_space_code);
	residual =  Ddot(N_Unknowns,defect,defect);
	impuls_residual = Ddot(2*N_U,defect,defect);
	// reset counter
	j = 0;
	pressure_separation = 1;
 	OutPut("p sep ");
	OutPut("nonlinear iteration step " << setw(3) << j);
	OutPut(setw(14) << impuls_residual);
	OutPut(setw(14) << residual-impuls_residual);
	OutPut(setw(14) << sqrt(residual) << endl);
     }

      switch(TDatabase::ParamDB->SOLVER_TYPE)
      {
        case AMG:
          TDatabase::ParamDB->SC_VERBOSE=0;     
          TDatabase::ParamDB->CC_VERBOSE=0;     
          t1 = GetTime();
          switch(TDatabase::ParamDB->NSTYPE)
          {
            case 1:
              Solver(sqmatrixA, matrixB1, matrixB2, rhs_high, sol);
            break;

            case 2:
#ifdef __MORTAR__
              Solver(sqmatrixA, matrixB1T, matrixB2T,
                     matrixB1, matrixB2, matrix_mortar, rhs_high, sol);
#else
              Solver(sqmatrixA, matrixB1T, matrixB2T, 
                     matrixB1, matrixB2, rhs_high, sol);
#endif
            break;

            case 3:
              Solver(sqmatrixA11, sqmatrixA12, sqmatrixA21, 
                     sqmatrixA22, matrixB1, matrixB2, rhs_high, sol);
            break;

            case 4:
              DirectSolver(sqmatrixA11, sqmatrixA12, sqmatrixA21, 
                     sqmatrixA22, matrixB1T, matrixB2T, 
                     matrixB1, matrixB2, rhs_high, sol);
            break;
          }
          t2 = GetTime();
          solver_time += (t2-t1);
          break;

        case GMG:
          t1 = GetTime();
          switch (TDatabase::ParamDB->SC_SOLVER_SADDLE)
          {
            case 11:
              if (TDatabase::ParamDB->SC_PRECONDITIONER_SADDLE == 5)
              {
                memcpy(itmethod_sol, sol, N_Unknowns*SizeOfDouble);
                memcpy(itmethod_rhs, rhs_high, N_Unknowns*SizeOfDouble);
              }
              break;
            case 16:
              if (TDatabase::ParamDB->SC_PRECONDITIONER_SADDLE == 5)
              {
                memcpy(itmethod_sol, sol, N_Unknowns*SizeOfDouble);
                memcpy(itmethod_rhs, rhs_high, N_Unknowns*SizeOfDouble);
              }
              break;
          }
          N_LinIter+=itmethod->Iterate(sqmatrices,matrices,itmethod_sol,itmethod_rhs);
          switch (TDatabase::ParamDB->SC_SOLVER_SADDLE)
          {
            case 11:
              if (TDatabase::ParamDB->SC_PRECONDITIONER_SADDLE == 5)
              {
                memcpy(sol, itmethod_sol, N_Unknowns*SizeOfDouble);
                memcpy(rhs_high, itmethod_rhs, N_Unknowns*SizeOfDouble);
              }
              break;
            case 16:
              if (TDatabase::ParamDB->SC_PRECONDITIONER_SADDLE == 5)
              {
                memcpy(sol, itmethod_sol, N_Unknowns*SizeOfDouble);
                memcpy(rhs_high, itmethod_rhs, N_Unknowns*SizeOfDouble);
              }
              break;
          }
          t2 = GetTime();    
          solver_time += (t2-t1);
          for(l=0;l<N_Unknowns;l++)
          {
            p2 = sol[l]-oldsol[l];
            sol[l] = oldsol[l] + omega * p2;
          }
          break; 
      }
    } // endfor
    
    if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE) 
      IntoL20FEFunction(sol+2*N_U, N_P,PSpaces[mg_level],
                        velocity_space_code, pressure_space_code);

    psi = new double[N_V];
//     StreamFunction(velocity_space, sol, sol+N_U, streamfunction_space, psi);
//     StreamFct = new TFEFunction2D(streamfunction_space, PsiString, PsiString, psi, N_V);
    // DivU(u1, u2, StreamFct);
//     Output->AddFEFunction(StreamFct);
/*
    for (l=0;l<N_U;l++)
       OutPut("u1(" << l+1 << ") = " << sol[l] << ";" << endl);
    for (l=0;l<N_U;l++)
       OutPut("u2(" << l+1 << ") = " << sol[l+N_U] << ";" << endl);
    for (l=0;l<N_P;l++)
       OutPut("p(" << l+1 << ") = " << sol[2*N_U+l] << ";" << endl);
*/    
    if(TDatabase::ParamDB->WRITE_GRAPE)
    {
      // OutPut("u " << Ddot(2*N_U,sol,sol)<< endl);
      //for (ii=0;ii<N_Unknowns;ii++)
      //  OutPut(ii << " " << sol[ii] << endl);
      // OutPut("p " << Ddot(N_P,sol+2*N_U,sol+2*N_U)<< endl);
      os.seekp(std::ios::beg);
      os << GrapeBaseName << i << ".dat" << ends;
      Output->WriteGrape(os.str().c_str());
    }

    U2Array[mg_level]->Interpolate(InitialU2);

    if(TDatabase::ParamDB->WRITE_VTK)
    {
	os.seekp(std::ios::beg);
	os << VtkBaseName <<".0000"<< img<< ".vtk" << ends;
	Output->WriteVtk(os.str().c_str());
        img++;
    }
    U2Array[mg_level]->Interpolate(ExactU2);

    if(TDatabase::ParamDB->WRITE_GNU)
    {
	os.seekp(std::ios::beg);
	os << GnuBaseName << i << ".gnu" << ends;
	Output->WriteGnuplot(os.str().c_str());
    }
//      if(TDatabase::ParamDB->WRITE_VTK)
//        {
//         os.seekp(std::ios::beg);
//         if(img<10) os << VtkBaseName<<".0000"<<img<<".vtk" << ends;
//          else if(img<100) os << VtkBaseName<<".000"<<img<<".vtk" << ends;
//          else if(img<1000) os << VtkBaseName<<".00"<<img<<".vtk" << ends;
//          else if(img<10000) os << VtkBaseName<<".0"<<img<<".vtk" << ends;
//          else  os << VtkBaseName<<"."<<img<<".vtk" << ends;
//         Output->WriteVtk(os.str().c_str());
//         img++;
//        }

    if(TDatabase::ParamDB->WRITE_MATLAB)
    {
	os.seekp(std::ios::beg);
	os << MatlabBaseName << i << ".m" << ends;
	Output->WriteMatlab(os.str().c_str());
    }

    // measure errors to known solution
    if(TDatabase::ParamDB->MEASURE_ERRORS)
    {
      fesp[0] = USpaces[mg_level];
      fefct[0] = U1Array[mg_level];
      fefct[1] = U2Array[mg_level];
      aux =  new TAuxParam2D(NSN_FESpacesVelo, NSN_FctVelo, 
                           NSN_ParamFctVelo, 
                           NSN_FEValuesVelo,
                           fesp, fefct, 
                           NSFctVelo, 
                           NSFEFctIndexVelo, NSFEMultiIndexVelo, 
                           NSN_ParamsVelo, NSBeginParamVelo);

      // errors in first velocity component
      U1Array[mg_level]->GetErrors(ExactU1, 3, NSAllDerivatives, 2, 
                            L2H1Errors, 
                            NULL, aux, 1, USpaces+mg_level, errors);
      l2u1[i] = errors[0];
      h1u1[i] = errors[1];

      // errors in second velocity component
      U2Array[mg_level]->GetErrors(ExactU2, 3, NSAllDerivatives, 2,
                            L2H1Errors, 
                            NULL, aux, 1, USpaces+mg_level, errors);
      l2u2[i] = errors[0];
      h1u2[i] = errors[1];
    
      // errors in pressure


      PArray[mg_level]->GetErrors(ExactP, 3, NSAllDerivatives, 2, 
                           L2H1Errors, 
                           NULL, aux, 1, PSpaces+mg_level, errors);    
      l2p[i] = errors[0];
      h1p[i] = errors[1];

//       for(j=0;j<100;j++)
//       {
//         x = -2.0+4.0*j/99;
//         y =  0.0;
//         PArray[mg_level]->FindGradient(x, y, errors);
//         cout << "h " << x << " " << y << " " << errors[0] << " ";
//         cout << errors[1] << " " << errors[2] << endl;
//       }
      
//       for(j=0;j<100;j++)
//       {
//         x = -2.0+4.0*j/99;
//         y = -2.0+4.0*j/99;
//         PArray[mg_level]->FindGradient(x, y, errors);
//         cout << "d " << x << " " << y << " " << errors[0] << " ";
//         cout << errors[1] << " " << errors[2] << endl;
//       }
      
      delete aux;

      // output of errors
      // coarsest level



      if(i<=FirstSolve)
      {
        p1 = sqrt(l2u1[i]*l2u1[i]+l2u2[i]*l2u2[i]);
        OutPut("L2(u): " <<  p1 << endl);
        p2 = sqrt(h1u1[i]*h1u1[i]+h1u2[i]*h1u2[i]);
        OutPut("H1-semi(u): " <<  p2 << endl);
        OutPut("L2(p): " <<  l2p[i] << endl);
        OutPut("H1-semi(p): " <<  h1p[i] << endl);
      }
      // not coarsest level
      else
      {
        errors[0] = sqrt(l2u1[i]*l2u1[i]+l2u2[i]*l2u2[i]) ;
        errors[1] = sqrt(l2u1[i-1]*l2u1[i-1]+l2u2[i-1]*l2u2[i-1]);
        OutPut("L2(u): " <<  errors[0] << " order ");
        OutPut(log(errors[1]/errors[0])/ln2 << endl);

        errors[0] = sqrt(h1u1[i]*h1u1[i]+h1u2[i]*h1u2[i]) ;
        errors[1] = sqrt(h1u1[i-1]*h1u1[i-1]+h1u2[i-1]*h1u2[i-1]);
        OutPut("H1-semi(u): " << errors[0] << " order ");
        OutPut(log(errors[1]/errors[0])/ln2 << endl);

        errors[0] = l2p[i];
        errors[1] = l2p[i-1];
        OutPut("L2(p): " <<  errors[0] << " order ");
        OutPut(log(errors[1]/errors[0])/ln2 << endl);

        errors[0] = h1p[i];
        errors[1] = h1p[i-1];
        OutPut("H1-semi(p): " << errors[0] << " order ");
        OutPut(log(errors[1]/errors[0])/ln2 << endl);

      }

      Int_intfaceU[i] = new double[2];

      Int_intfaceU[i][0] = 0.0;
      Int_intfaceU[i][1] = 0.0;

    switch((int)(TDatabase::ParamDB->REACTOR_P28+0.5))
    {
      case 12345:
           Cal_Intface_U(UArray[mg_level], Int_intfaceU[i]);
      break;

      case 54321:
      case  6789:
      case  9876:
           Intface_U(UArray[mg_level], Int_intfaceU[i]);
      break;
    }


//       OutPut( "u.n : " <<Int_intfaceU[i][0] <<  " u : " << Int_intfaceU[i][1]<<endl);

      OutPut("L2(u.n): " <<  Int_intfaceU[i][0] << " order ");
      if(i) OutPut(log(Int_intfaceU[i-1][0]/Int_intfaceU[i][0])/ln2 << endl);
      OutPut( " u : " << Int_intfaceU[i][1]<<endl);

//       OutPut("Int_intface(U): " << Int_intfaceU[i][0] + Int_intfaceU[i][1]<<endl;);

    } // endif MEASURE_ERRORS

#ifdef __BENCH__
    // compute characteristic values (deltaP, Cd, Cl)
    GetCdCl(U1Array[mg_level], U2Array[mg_level], PArray[mg_level],
	    U1Array[mg_level], U2Array[mg_level], Cd, Cl);

    PArray[mg_level]->FindGradient(0.15, 0.2, dP1);
    PArray[mg_level]->FindGradient(0.25, 0.2, dP2);

    OutPut( "C_drag = " << setprecision(16) <<Cd );
    OutPut( " C_lift = " << setprecision(16) << Cl);
    OutPut( " deltaP = " << setprecision(16) << dP1[0] - dP2[0] << endl);
    OutPut( setprecision(7) << endl);
#endif

#ifdef __CHANNELSTEP__
    GetReattachmentPoint(U1Array[mg_level], reatt_pt);
    OutPut( "reattachment: " << reatt_pt<< endl);    
#endif

#ifdef __CHANNEL30__
    GetReattachmentPoint(U1Array[mg_level], reatt_point);
    OutPut( "reattachment: lower " << reatt_point[0] << " left upper " <<
            reatt_point[1] << " right upper " <<  reatt_point[2] << endl);    
#endif

   // remove data which will not be used later
    delete oldsol;
    delete defect;
    delete psi;
   
    if ((mg_type==1)||(TDatabase::ParamDB->SOLVER_TYPE == AMG))
    {
      switch(TDatabase::ParamDB->NSTYPE)
      {
        case 1:
          delete matrixB1;
          delete matrixB2;
          delete sqmatrixA;
          break;
        case 2:
          delete matrixB1;
          delete matrixB2;
          delete matrixB1T;
          delete matrixB2T;
          delete sqmatrixA;
          break;
        case 3:
          delete matrixB1;
          delete matrixB2;
          delete sqmatrixA11;
          delete sqmatrixA12;
          delete sqmatrixA21;
          delete sqmatrixA22;
          break;
        case 4:
          delete matrixB1;
          delete matrixB2;
          delete matrixB1T;
          delete matrixB2T;
          delete sqmatrixA11;
          delete sqmatrixA12;
          delete sqmatrixA21;
          delete sqmatrixA22;
          break;
      }
      delete structureB;
      delete structureBT;
      delete sqstructureA;
      delete rhs_high;
//       delete MGLevel;
    } // end if (mg_type==1)

    old_sol = sol;
    old_u = u;
    old_p = p;
    old_u_space = velocity_space;
    old_p_space = pressure_space;
       
    OutPut("memory after: " << setw(10) << GetMemory() << endl);
  } // endfor i

  OutPut("used time: " << GetTime() << endl);
  OutPut("used bytes: " << GetMemory() << endl);
  CloseFiles();

  return 0;
}
