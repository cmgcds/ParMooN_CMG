// =======================================================================
// 
// Purpose:     Main program for impinging droplet
//
// Author:     Sashikumaar Ganesan
// modified    10.06.2010 
//             20.01.2012
//             27.07.2014 (surfact coupling updated, solution in new timestep domain)
// ======================================================================= 

#include <Domain.h>
#include <Database.h>
#include <DiscreteForm2D.h>
#include <FEDatabase2D.h>
#include <FESpace2D.h>
#include <SquareStructure1D.h>
#include <SquareStructure2D.h>
#include <SquareMatrix1D.h>
#include <FEFunction1D.h>
#include <Structure2D.h>
#include <AuxParam2D.h>
#include <DirectSolver.h>
#include <Assemble2D.h>
#include <Output2D.h>
#include <DiscreteForm2D.h>
#include <LinAlg.h>
#include <TNSE2D_ParamRout.h>
#include <AuxParam2D.h>

#include <Collection.h>
#include <NodalFunctional2D.h>
#include <string.h>
#include <sstream>
#include <MooNMD_Io.h>
#include <stdlib.h>
// #include <malloc.h>

#include <Upwind.h>
#include <NSE_MultiGrid.h>
#include <NSE_MGLevel1.h>
#include <NSE_MGLevel2.h>
#include <NSE_MGLevel3.h>
#include <NSE_MGLevel4.h>
#include <Convolution.h>
#include <ItMethod.h>
#include <FixedPointIte.h>
#include <MultiGridIte.h>
#include <FgmresIte.h>

#include <MultiGrid2D.h>
#include <MGLevel2D.h>
#include <FreeSurface2D.h>

#include <MainUtilities.h>
// #include <TimeUtilities.h>
#include <TimeDiscRout.h>

#include <MacroCell.h>
#include <BoundEdge.h>
#include <IsoBoundEdge.h>
#include <IsoInterfaceJoint.h>
#include <BdLine.h>
#include <BdCircle.h>
#include <GridCell.h>

#include <QuadAffin.h>
#include <QuadBilinear.h>
#include <QuadIsoparametric.h>
#include <TriaAffin.h>
#include <TriaIsoparametric.h>

#  include <sys/stat.h>
#  include <sys/types.h>

// #include "../Examples/TNSE_2D/Droponsolid.h"
// #include "../Examples/TNSE_2D/Drop_Imping_Axial3D.h"
// #include "../Examples/TNSE_2D/ChannelDroplet.h"
#include "../TNSE_2D/HangingDroplet.h"

// extern "C"
// {
//   void triangulate(char*, struct triangulateio*,
//                    struct triangulateio*, struct triangulateio*);
// }



void PrintSurfSurfactant(int N, TVertex **Vertex, TFEFunction2D *Surfact, int &N_BData)
{
 int i, j, k, Cell_No, IJoint, N_DOF_Local;
 double  T_val[3], x1, y1, x2, y2, ArcLength=0.;
 char *VtkBaseName;
 VtkBaseName = TDatabase::ParamDB->VTKBASENAME;

 std::ostringstream os;
 os << " ";

 os.seekp(std::ios::beg);
//       os << "surfact"<< i << ".dat" << ends;
  if(N_BData<10) os << "BDData/"<<VtkBaseName<<"Gamma_0000"<<N_BData<<".data" << ends;
  else if(N_BData<100) os <<"BDData/"<<VtkBaseName<<"Gamma_000"<<N_BData<<".data" << ends;
  else if(N_BData<1000) os <<"BDData/"<<VtkBaseName<<"Gamma_00"<<N_BData<<".data" << ends;
  else if(N_BData<10000) os <<"BDData/"<<VtkBaseName<<"Gamma_0"<<N_BData<<".data" << ends;
  else  os <<"BDData/"<<VtkBaseName<<"Gamma_"<<N_BData<<".data" << ends;

  std::ofstream dat(os.str().c_str());

//   cout << " PrintSurfSurfactant " << N<<endl;
  
  if (!dat)
   {
    cerr << "cannot open file for output" << endl;
    exit(0);
   }
  dat << "%% Surfactant data created by MooNMD" << endl;
  dat << "%% Current Reference Time :" << TDatabase::TimeDB->CURRENTTIME << endl;
  dat << "%% x, y, ArcLength,  Surfact" << endl;

  Vertex[0]->GetCoords(x1, y1);
  for(i=0;i<N;i++)
   {
    Vertex[i]->GetCoords(x2, y2);


    ArcLength += sqrt( (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) );

    T_val[0] = 0.;
    if(Surfact)
     Surfact->FindGradient( x1,  y1, T_val);
    dat << x1<< " " << y1 << " "<< ArcLength << "  " << T_val[0] <<endl;
//     cout<<i<< " Angle of free Vertices "<<(180./Pi)*atan2(y2,x2)<<"  " << T_val[0]<<endl;
    x1=x2; y1=y2;
   }

      dat.close();
      cout << endl;
      cout << "Surfactant data wrote into file " << endl;
 N_BData++;

} // PrintSurfSurfactant


// ====================================================================
// Get the inner angles of the cells in whole domain
// ====================================================================

void MapDomainToSurf(TFEFunction2D *Fe2D, TFEFunction1D *Fe1D,
                     int *Cell_array, int *Joint_array)
{
  int i,j,k,l,n1,n2, N_Cells1D, N_Cells2D, N_DOF, N_DOF1D, N;
  int *GlobalNumbers1D, *GlobalNumbers2D,  *BeginIndex1D, *BeginIndex2D, *JointDOF;
  double  *Values1D, *Values2D;
  TCollection *Coll1D, *Coll2D;
  TFESpace1D *SurfSpace;
  TFESpace2D *FeSpace;
  TBaseCell *Me;
  TJoint *joint;
  FE2D FeId;
  TFEDesc2D *FeDesc;
  FE1D FeId1D;
  TFE1D *Element;

  SurfSpace = Fe1D->GetFESpace1D();
  FeSpace = Fe2D->GetFESpace2D();
  Values1D= Fe1D->GetValues();
  Values2D= Fe2D->GetValues();

  GlobalNumbers1D = SurfSpace->GetGlobalNumbers();
  BeginIndex1D = SurfSpace->GetBeginIndex();

  GlobalNumbers2D = FeSpace->GetGlobalNumbers();
  BeginIndex2D = FeSpace->GetBeginIndex();


  Coll1D = SurfSpace->GetCollection();
  Coll2D = FeSpace->GetCollection();
  N_Cells1D = Coll1D->GetN_Cells();
  N_Cells2D = Coll2D->GetN_Cells();

//   cout <<N_Cells1D<< "  test MapSurfToDomain " << N_Cells2D << endl;
  for(i=0;i<N_Cells1D;i++)
   {
     N = Cell_array[i];
     Me = Coll2D->GetCell(N);
     l = Joint_array[i];
     joint = Me->GetJoint(l);
     FeId = FeSpace->GetFE2D(N, Me);
     FeDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(FeId);
     JointDOF = FeDesc->GetJointDOF(l);
     N_DOF = FeDesc->GetN_JointDOF();

     FeId1D = SurfSpace->GetFE1D(i, Coll1D->GetCell(i));
     Element = TFEDatabase2D::GetFE1D(FeId1D);
     N_DOF1D = Element->GetN_DOF();

     if(N_DOF1D !=N_DOF )
      {
        cout << "Number of degrees of freedon on surface line should be matched with triangular edge" << endl;
        cout << " Check MapDomainToSurf !!!!" << endl;
        cout << "N_DOF1D " << N_DOF1D << " N_DOF2D " << N_DOF << endl;
        exit(0);
      }

     for(j=0;j<N_DOF1D;j++)
      {
       n1 = GlobalNumbers1D[BeginIndex1D[i] + j];
       n2 = GlobalNumbers2D[BeginIndex2D[N]+JointDOF[j]];
       Values1D[n1] = Values2D[n2];
      }// endfor j
   } //  for(i=0;i<N_Cel

//   cout <<N_Cells1D<< "  2 test MapSurfToDomain " << N_Cells2D << endl;

}


 
 void MapSurfToDomain(TFEFunction1D *Fe1D, TFEFunction2D *Fe2D,
                     int *Cell_array, int *Joint_array)
{
  int i,j,k,l,n1,n2, N_Cells1D, N_Cells2D, N_DOF, N_DOF1D, N;
  int *GlobalNumbers1D, *GlobalNumbers2D,  *BeginIndex1D, *BeginIndex2D, *JointDOF;
  double  *Values1D, *Values2D;
  TCollection *Coll1D, *Coll2D;
  TFESpace1D *SurfSpace;
  TFESpace2D *FeSpace;
  TBaseCell *Me;
  TJoint *joint;
  FE2D FeId;
  TFEDesc2D *FeDesc;
  FE1D FeId1D;
  TFE1D *Element;

  SurfSpace = Fe1D->GetFESpace1D();
  FeSpace = Fe2D->GetFESpace2D();
  Values1D= Fe1D->GetValues();
  Values2D= Fe2D->GetValues();

  GlobalNumbers1D = SurfSpace->GetGlobalNumbers();
  BeginIndex1D = SurfSpace->GetBeginIndex();

  GlobalNumbers2D = FeSpace->GetGlobalNumbers();
  BeginIndex2D = FeSpace->GetBeginIndex();

  Coll1D = SurfSpace->GetCollection();
  Coll2D = FeSpace->GetCollection();
  N_Cells1D = Coll1D->GetN_Cells();
  N_Cells2D = Coll2D->GetN_Cells();

//   cout <<N_Cells1D<< "  test MapSurfToDomain " << N_Cells2D << endl;
//   exit(0);
  for(i=0;i<N_Cells1D;i++)
   {
     N = Cell_array[i];
     Me = Coll2D->GetCell(N);
     l = Joint_array[i];
     joint = Me->GetJoint(l);
     FeId = FeSpace->GetFE2D(N, Me);
     FeDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(FeId);
     JointDOF = FeDesc->GetJointDOF(l);
     N_DOF = FeDesc->GetN_JointDOF();

     FeId1D = SurfSpace->GetFE1D(i, Coll1D->GetCell(i));
     Element = TFEDatabase2D::GetFE1D(FeId1D);
     N_DOF1D = Element->GetN_DOF();

     if(N_DOF1D !=N_DOF )
      {
        cout << "Number of degrees of freedon on surface line should be matched with triangular edge" << endl;
        cout << " Chech MapDomainToSurf !!!!" << endl;
        cout << "N_DOF1D " << N_DOF1D << " N_DOF2D " << N_DOF << endl;
        exit(0);
      }

     for(j=0;j<N_DOF1D;j++)
      {
       n1 = GlobalNumbers1D[BeginIndex1D[i] + j];
       n2 = GlobalNumbers2D[BeginIndex2D[N]+JointDOF[j]];
       Values2D[n2] = Values1D[n1];
      }// endfor j


   } //  for(i=0;i<N_Cel

//   cout <<N_Cells1D<< "  test MapSurfToDomain " << N_Cells2D << endl;
//   exit(0);   
}



void Surfact2D_InterfaceInt(TSquareMatrix2D *A, double *rhs, BoundCondFunct2D SurfactBoundCondition, 
                           TFEFunction2D *FeFunct, TFEFunction1D *SurfaceFeFunct, 
                           TFESpace2D *GlobalFE_Space, int *Cell_array, int *Joint_array, double *XmaxVal)
{
  int i, j, k, l, N_Cells, N_Vertices, N_Edges, ORDER;
  int *BeginIndex, *GlobalNumbers, *DOF, *DOF_Surface, TestDOF, AnsatzDOF;
  int *BeginIndex_Surface, *GlobalNumbers_Surface, dof_surface;
  int JointNumbers[MAXN_JOINTS], N_IsoJoints, N_JointDOF;
  int *KCol, *RowPtr, *JointDOF, N_DOF, N, IJoint, local_dof, m;
  int N_LinePoints, ActiveBound, N_BaseFunct_Surface, Global_N;
  int N_BaseFunct, *N_BaseFuncts, index1, Surf_cell_No, N_Cells_Surface;

  double **uref, **uxiref, **uetaref, val, rhsval, *V, *Surfactant;
  double uorig[MaxN_BaseFunctions2D], uxorig[MaxN_BaseFunctions2D];
  double uyorig[MaxN_BaseFunctions2D], U, test00, C;
  double *LineWeights, *zeta, X_B[100], Y_B[100];
  double t0, t1, r_axial, normn, surf, surfx, surfy, surf_flux, CMass, GammaMass;
//   double C_Infty, Gamma_Infty, Bi_a, Bi_d;
  double Bi, Da, beta;
   
  BoundCond Cond0, Cond1;
  FE2D FEId;
  FE1D FEId_Surface;
  TFE2D *ele, *ele_surf;
  TFE1D *Element;
  RefTrans2D RefTrans;
  TRefTrans2D *F_K;
  BF2DRefElements RefElement;
  QuadFormula1D LineQuadFormula;
  TQuadFormula1D *qf1;
  TBaseCell *Me, *Me_Surface;
  TCollection *GlobalColl, *Coll, *Coll_Surface;
  TFESpace2D *fespace;
  TFESpace1D *fespace_surface;
  TJoint *joint;
  TFEDesc2D *FeDesc;
  TFEDesc1D *FeDesc_low;

  GlobalColl = GlobalFE_Space->GetCollection();

  fespace = A->GetFESpace();    // outer phase
  Coll = fespace->GetCollection();
  N_Cells = Coll->GetN_Cells();
  BaseFunct2D *BaseFuncts;


  ActiveBound = fespace->GetActiveBound();
  GlobalNumbers = fespace->GetGlobalNumbers();
  BeginIndex = fespace->GetBeginIndex();
  Surfactant = FeFunct->GetValues();

  BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
  N_BaseFuncts = TFEDatabase2D::GetN_BaseFunctFromFE2D();

  fespace_surface=SurfaceFeFunct->GetFESpace1D();
  Coll_Surface = fespace_surface->GetCollection();
  N_Cells_Surface = Coll_Surface->GetN_Cells();

  GlobalNumbers_Surface = fespace_surface->GetGlobalNumbers();
  BeginIndex_Surface = fespace_surface->GetBeginIndex();
  V = SurfaceFeFunct->GetValues();

//   Gamma_Infty = TDatabase::ParamDB->REACTOR_P13;
//   C_Infty = TDatabase::ParamDB->REACTOR_P14;
//   Bi_a = TDatabase::ParamDB->REACTOR_P15;
//   Bi_d = TDatabase::ParamDB->REACTOR_P16;
  
  Bi = TDatabase::ParamDB->REACTOR_P13; 
  Da = TDatabase::ParamDB->REACTOR_P14;
  beta = TDatabase::ParamDB->REACTOR_P15;
    
  XmaxVal[0] = -1.e10;
  XmaxVal[1] =  0.;
  surf_flux =  0.;

 for(i=0;i<N_Cells_Surface;i++)
  {
    Me_Surface = Coll_Surface->GetCell(i);
    FEId_Surface = fespace_surface->GetFE1D(i, Me_Surface);
    Element = TFEDatabase2D::GetFE1D(FEId_Surface);
    N_BaseFunct_Surface = Element->GetN_DOF();
    DOF_Surface = GlobalNumbers_Surface + BeginIndex_Surface[i];

    Global_N = Cell_array[i];     // outerphase
    Me = GlobalColl->GetCell(Global_N);
//     N = Me->GetLocalCellNo();
    N=Global_N; // free surf flows N=Global_N;
    IJoint = Joint_array[i];

    FEId = fespace->GetFE2D(N, Me);
    ele = TFEDatabase2D::GetFE2D(FEId);

    RefElement = TFEDatabase2D::GetRefElementFromFE2D(FEId);
    FeDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(FEId);
    N_JointDOF =  FeDesc->GetN_JointDOF();
    JointDOF = FeDesc->GetJointDOF(IJoint);
    N_BaseFunct = FeDesc->GetN_DOF();
    DOF = GlobalNumbers + BeginIndex[N];

    l = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
    ORDER = TFEDatabase2D::GetAccuracyFromFE2D(FEId);
    LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(2*l);
    qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
    qf1->GetFormulaData(N_LinePoints, LineWeights, zeta);

    switch(RefElement)
     {
      case BFUnitSquare:
          RefTrans = QuadIsoparametric;
          F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
          ((TQuadIsoparametric *)F_K)->SetCell(Me);
          ((TQuadIsoparametric *)F_K)->GetOrigBoundFromRef(IJoint, N_LinePoints, zeta, X_B, Y_B);
      break;

      case BFUnitTriangle:
          RefTrans = TriaIsoparametric;
          F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
          ((TTriaIsoparametric *)F_K)->SetApproximationOrder(ORDER);  
          ((TTriaIsoparametric *)F_K)->SetCell(Me);
          ((TTriaIsoparametric *)F_K)->GetOrigBoundFromRef(IJoint, N_LinePoints, zeta, X_B, Y_B);
      break;
    } // endswitch

    TFEDatabase2D::GetBaseFunct2DFromFE2D(FEId)
                  ->MakeRefElementData(LineQuadFormula);


    uref = TFEDatabase2D::GetJointValues2D(BaseFuncts[FEId],
                        LineQuadFormula, IJoint);
    uxiref = TFEDatabase2D::GetJointDerivatives2D(BaseFuncts[FEId],
                        LineQuadFormula, IJoint, D10);
    uetaref = TFEDatabase2D::GetJointDerivatives2D(BaseFuncts[FEId],
                        LineQuadFormula, IJoint, D01);

    for(k=0;k<N_LinePoints;k++)
     {

         switch(RefElement)
          {
            case BFUnitSquare:
              ((TQuadIsoparametric *)F_K)->SetCell(Me);
              ((TQuadIsoparametric *)F_K)->GetOrigValues(IJoint,  zeta[k],
                        N_BaseFunct, uref[k], uxiref[k], uetaref[k],
                        uorig, uxorig, uyorig);
              ((TQuadIsoparametric *) F_K)->GetTangent(IJoint, zeta[k], t0, t1);
            break;

            case BFUnitTriangle:
              ((TTriaIsoparametric *)F_K)->SetApproximationOrder(ORDER);
              ((TTriaIsoparametric *)F_K)->SetCell(Me); 
              ((TTriaIsoparametric *)F_K)->GetOrigValues(IJoint,  zeta[k],
                        N_BaseFunct, uref[k], uxiref[k], uetaref[k],
                        uorig, uxorig, uyorig);
              ((TTriaIsoparametric *) F_K)->GetTangent(IJoint, zeta[k], t0, t1);
            break;
          } // endswitch

     normn = sqrt(t0*t0+t1*t1);
     r_axial = fabs(X_B[k]);
     
         if(X_B[k]<=0)
          {
           cout <<"X_B[k] negative in SurfInt  change Quad rule " <<  X_B[k] <<endl;
//         exit(0);
          }

     // get the outer phase surfactant value at this point
     C=0.;
     for(l=0;l<N_BaseFunct;l++)
       C  += Surfactant[DOF[l]]*uorig[l];


     // get surfactant solution  on the interface
     U=0.;
#ifndef __MASSTRANSTEST__
     for(l=0;l<N_BaseFunct_Surface;l++)
      {
       local_dof   = JointDOF[l];
       test00 = uorig[local_dof];
       m = DOF_Surface[l];
       U  += V[m]*test00;
      }
#endif
     
     //U=0.; // for bulk concentration test problem
//     if( XmaxVal[0]< r_axial )
//      {
//       XmaxVal[0] = r_axial;
//       XmaxVal[1] = U;
//      }
     if( XmaxVal[1]< U )
      {
       XmaxVal[0] = r_axial;
       XmaxVal[1] = U;
      }
     rhsval =  -beta*C*(1. - U) + Bi*Da*U;
     rhsval *= normn * LineWeights[k]*r_axial;
     
//   OutPut("Flux rhsval " <<  rhsval <<endl);
  
//         // update rhs for all test functions
     for(l=0;l<N_BaseFunct;l++)
       {
         if((index1 = DOF[l])<ActiveBound)
          {
           rhs[index1] += rhsval*uorig[l];
           surf_flux += rhsval*uorig[l];
//            cout <<"surf_flux " << surf_flux << " rhs "<<  rhsval*uorig[l] <<endl;
          }
       }//  for(l=0;l<N_BaseFunct;l++)
//          cout  <<endl;
    }  // for(k=0;k<N_LinePoints;k
  } // for(i=0;i<N_Cells;

//   OutPut("Flux at the interface " << 2*Pi*surf_flux <<endl);
//   cout << " GammaMass on Interface; " << 2*Pi*GammaMass <<endl;
//   cout << " CMass on Interface; " << 2*Pi*CMass <<endl;
//    cout << " xpos; " << xpos << " XmaxVal; " << XmaxVal << endl;

// exit(0);
}


void Surfact2D_InterfaceInt_Implicit(TSquareMatrix2D *A, double *rhs,
                           BoundCondFunct2D SurfactBoundCondition,
                           TFEFunction2D *FeFunct, TFEFunction1D *SurfaceFeFunct, 
                           TFESpace2D *GlobalFE_Space, int *Cell_array, int *Joint_array, double *XmaxVal)
{
  int i, j, k, l, N_Cells, N_Vertices, N_Edges, ORDER;
  int *BeginIndex, *GlobalNumbers, *DOF, *DOF_Surface, TestDOF, AnsatzDOF;
  int *BeginIndex_Surface, *GlobalNumbers_Surface, dof_surface;
  int JointNumbers[MAXN_JOINTS], N_IsoJoints, N_JointDOF;
  int *KCol, *RowPtr, *JointDOF, N_DOF, N, IJoint, local_dof, m;
  int N_LinePoints, ActiveBound, N_BaseFunct_Surface, Global_N, index2;
  int N_BaseFunct, *N_BaseFuncts, index1, Surf_cell_No, N_Cells_Surface;

  double **uref, **uxiref, **uetaref, val, rhsval, *V, *Surfactant;
  double uorig[MaxN_BaseFunctions2D], uxorig[MaxN_BaseFunctions2D];
  double uyorig[MaxN_BaseFunctions2D], U, test00, C;
  double *LineWeights, *zeta, X_B[100], Y_B[100];
  double t0, t1, r_axial, normn, surf, surfx, surfy, surf_flux, CMass, GammaMass;
//   double C_Infty, Gamma_Infty, Bi_a, Bi_d, *ValuesA;
  double Bi, Da, beta, *ValuesA;
  
  BoundCond Cond0, Cond1;
  FE2D FEId;
  FE1D FEId_Surface;
  TFE2D *ele, *ele_surf;
  TFE1D *Element;
  RefTrans2D RefTrans;
  TRefTrans2D *F_K;
  BF2DRefElements RefElement;
  QuadFormula1D LineQuadFormula;
  TQuadFormula1D *qf1;
  TBaseCell *Me, *Me_Surface;
  TCollection *GlobalColl, *Coll, *Coll_Surface;
  TFESpace2D *fespace;
  TFESpace1D *fespace_surface;
  TJoint *joint;
  TFEDesc2D *FeDesc;
  TFEDesc1D *FeDesc_low;

  GlobalColl = GlobalFE_Space->GetCollection();

  fespace = A->GetFESpace();    // outer phase
  Coll = fespace->GetCollection();
  N_Cells = Coll->GetN_Cells();
  BaseFunct2D *BaseFuncts;

  ActiveBound = fespace->GetActiveBound();
  GlobalNumbers = fespace->GetGlobalNumbers();
  BeginIndex = fespace->GetBeginIndex();
  Surfactant = FeFunct->GetValues();

  ValuesA  = A->GetEntries();
  RowPtr = A->GetRowPtr();
  KCol = A->GetKCol();

 
  BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
  N_BaseFuncts = TFEDatabase2D::GetN_BaseFunctFromFE2D();
  fespace_surface=SurfaceFeFunct->GetFESpace1D();
  Coll_Surface = (SurfaceFeFunct->GetFESpace1D())->GetCollection();
  N_Cells_Surface = Coll_Surface->GetN_Cells();
 
  GlobalNumbers_Surface = fespace_surface->GetGlobalNumbers();
  BeginIndex_Surface = fespace_surface->GetBeginIndex();
  V = SurfaceFeFunct->GetValues();

  Bi = TDatabase::ParamDB->REACTOR_P13; 
  Da = TDatabase::ParamDB->REACTOR_P14;
  beta = TDatabase::ParamDB->REACTOR_P15;
  
  XmaxVal[0] = -1.e10;
  XmaxVal[1] =  0.;

  
 for(i=0;i<N_Cells_Surface;i++)
  {
    Me_Surface = Coll_Surface->GetCell(i);
    FEId_Surface = fespace_surface->GetFE1D(i, Me_Surface);
    Element = TFEDatabase2D::GetFE1D(FEId_Surface);
    N_BaseFunct_Surface = Element->GetN_DOF();
    DOF_Surface = GlobalNumbers_Surface + BeginIndex_Surface[i];
    Global_N = Cell_array[i];     // outerphase
    Me = GlobalColl->GetCell(Global_N);
//     N = Me->GetLocalCellNo();
    N=Global_N; // free surf flows N=Global_N;
    IJoint = Joint_array[i];

    FEId = fespace->GetFE2D(N, Me);
    ele = TFEDatabase2D::GetFE2D(FEId);
    RefElement = TFEDatabase2D::GetRefElementFromFE2D(FEId);
    FeDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(FEId);
    N_JointDOF =  FeDesc->GetN_JointDOF();
    JointDOF = FeDesc->GetJointDOF(IJoint);
    
    N_BaseFunct = FeDesc->GetN_DOF();
    DOF = GlobalNumbers + BeginIndex[N];
    l = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
    ORDER = TFEDatabase2D::GetAccuracyFromFE2D(FEId);
    LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(2*l);
    qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
    qf1->GetFormulaData(N_LinePoints, LineWeights, zeta);

    switch(RefElement)
     {
      case BFUnitSquare:
          RefTrans = QuadIsoparametric;
          F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
          ((TQuadIsoparametric *)F_K)->SetCell(Me);
          ((TQuadIsoparametric *)F_K)->GetOrigBoundFromRef(IJoint,
	                            N_LinePoints, zeta, X_B, Y_B);
      break;

      case BFUnitTriangle:
          RefTrans = TriaIsoparametric;
          F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
          ((TTriaIsoparametric *)F_K)->SetApproximationOrder(ORDER);  
          ((TTriaIsoparametric *)F_K)->SetCell(Me);
          ((TTriaIsoparametric *)F_K)->GetOrigBoundFromRef(IJoint, N_LinePoints, zeta, X_B, Y_B);
      break;
    } // endswitch

    TFEDatabase2D::GetBaseFunct2DFromFE2D(FEId)
                  ->MakeRefElementData(LineQuadFormula);

    uref = TFEDatabase2D::GetJointValues2D(BaseFuncts[FEId],
                        LineQuadFormula, IJoint);
    uxiref = TFEDatabase2D::GetJointDerivatives2D(BaseFuncts[FEId],
                        LineQuadFormula, IJoint, D10);
    uetaref = TFEDatabase2D::GetJointDerivatives2D(BaseFuncts[FEId],
                        LineQuadFormula, IJoint, D01);

    for(k=0;k<N_LinePoints;k++)
     {

         switch(RefElement)
          {
            case BFUnitSquare:
              ((TQuadIsoparametric *)F_K)->SetCell(Me);
              ((TQuadIsoparametric *)F_K)->GetOrigValues(IJoint,  zeta[k],
                        N_BaseFunct, uref[k], uxiref[k], uetaref[k],
                        uorig, uxorig, uyorig);
              ((TQuadIsoparametric *) F_K)->GetTangent(IJoint, zeta[k], t0, t1);
            break;

            case BFUnitTriangle:
              ((TTriaIsoparametric *)F_K)->SetApproximationOrder(ORDER);
              ((TTriaIsoparametric *)F_K)->SetCell(Me); 
              ((TTriaIsoparametric *)F_K)->GetOrigValues(IJoint,  zeta[k],
                        N_BaseFunct, uref[k], uxiref[k], uetaref[k],
                        uorig, uxorig, uyorig);
              ((TTriaIsoparametric *) F_K)->GetTangent(IJoint, zeta[k], t0, t1);
            break;
          } // endswitch

      normn = sqrt(t0*t0+t1*t1);
      r_axial = fabs(X_B[k]);
     
         if(X_B[k]<=0)
          {
           cout <<"X_B[k] negative in SurfInt  change Quad rule " <<  X_B[k] <<endl;
//         exit(0);
          }


     // get surfactant solution  on the interface
     U=0.;
#ifndef __MASSTRANSTEST__
     for(l=0;l<N_BaseFunct_Surface;l++)
      {
       local_dof   = JointDOF[l];
       test00 = uorig[local_dof];
       m = DOF_Surface[l];
       U  += V[m]*test00;
      }
#endif

//      if( XmaxVal[0]< r_axial )
//       {
//        XmaxVal[0] = r_axial;
//        XmaxVal[1] = U;
//       }

     if( XmaxVal[1]< U )
      {
       XmaxVal[0] = r_axial;
       XmaxVal[1] = U;
      }
      
      rhsval = Bi*Da*U;
      rhsval *= normn*LineWeights[k]*r_axial;

//    update the matrix for all test functions
      for(l=0;l<N_BaseFunct;l++)
       {
        TestDOF = DOF[l];

        if(TestDOF<ActiveBound)
          rhs[TestDOF] += rhsval*uorig[l];

        index2 = RowPtr[TestDOF+1];

        for(m=0;m<N_BaseFunct;m++)
         {
          AnsatzDOF = DOF[m];
          // cout << AnsatzDOF << " -- " << TestDOF << endl;
          index1 = RowPtr[TestDOF];

          if(index1+1 == index2) continue;
           while(KCol[index1] != AnsatzDOF) index1++;

          val =  beta*(1. - U)*uorig[m]*uorig[l];
          val *= normn*LineWeights[k]*r_axial;

          ValuesA[index1] += val;
         } // for(m=0;m<N_BaseFunct;m++)
       }//  for(l=0;l<N_BaseFunct;l++)
    }  // for(k=0;k<N_LinePoints;k

  } // for(i=0;i<N_Cells;

//   cout<< " Flux at the interface " << 2*Pi*surf_flux <<endl;
//   cout << " GammaMass on Interface; " << 2*Pi*GammaMass <<endl;
//   cout << " CMass on Interface; " << 2*Pi*CMass <<endl;
//    cout << " XmaxVal; " << XmaxVal[1] << endl;

// exit(0);
}


void  AssembleSurf1D_SolubleSurfact(int n_fespaces, TFESpace2D **fespaces, TFEFunction2D **fefunctions,
                     int N_FESpaces_low, TFESpace1D **fespaces_low, TFEFunction1D *SurfaceFeFunct, int N_SquareMatrices,
                     TSquareMatrix1D **sqmatrices_low, int N_Rhs, double **RHSs, 
                     TFESpace1D **ferhs_low, int *Cell_array, int *Joint_array, double *C_Outer)
{
  int i, j, k, l, m, n, N_Cells_low, N, N_LocalUsedElements, local_i, local_j, ORDER;
  int N_BaseFunct, N_BaseFunct_low,  N_Points, N_JointDOF, begin, end, *N_BaseFuncts;
  int *BeginIndex_low, *GlobalNumbers_low, *DOF, *DOF_LOW, TestDOF, AnsatzDOF, IJoint ;
  int *BeginIndex, *GlobalNumbers, *GlobalNumbers_Outer, *BeginIndex_Outer, N_BaseFunct_Outer;
  int LocN_BF[N_BaseFuncts2D], N_LinePoints, *KCol, *RowPtr, *JointDOF, N_Outer;
  int *DOF_Outer;

  double x0, y0, x1, y1, t0, t1, n0, n1, normn;
  double AbsDetjk[MaxN_QuadPoints_2D], Mult;
  double *weights, *xi, *eta;
  double **uref, **uxiref, **uetaref;
  double **uref_Outer, **uxiref_Outer, **uetaref_Outer,  C, CX, CY;
  double *LineWeights, *zeta, *ValuesA, *ValuesM;
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];
  double uorig[MaxN_BaseFunctions2D], uxorig[MaxN_BaseFunctions2D], uyorig[MaxN_BaseFunctions2D];
  double uorig_Outer[MaxN_BaseFunctions2D], uxorig_Outer[MaxN_BaseFunctions2D];
  double uyorig_Outer[MaxN_BaseFunctions2D];
  double c0, r2;

  double val, rhsval, theta, ngrad_ansatz, ngrad_test, TangDivU;
  double  X_B[100], Y_B[100], r_axial, d1, d2, e1, e2;
  double LocMatrixA[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double LocMatrixM[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double LocRhs[MaxN_BaseFunctions2D];
  double ansatz00, ansatz10, ansatz01 ;
  double test00, test10, test01, *u1, *u2, u1x, u2x, u1y, u2y, U1;
  double *RHS, *V, U, AddedMass, CMass, GammaMass;
  double Bi, Da, beta, Pe_s, NGrad_C;
  
  BaseFunct2D LocBF[N_BaseFuncts2D];
  BaseFunct2D *BaseFuncts;
  boolean *SecondDer;
  BF2DRefElements RefElement;
  RefTrans2D RefTrans;
  TRefTrans2D *F_K;
  QuadFormula1D LineQuadFormula;
  TQuadFormula1D *qf1;
  TFEDesc2D *FeDesc, *FeDesc_Outer;
  TFEDesc1D *FeDesc_low;
  TCollection *Coll, *Coll_low;
  TBaseCell *Me, *Me_low;
  FE2D FEId, FEId_Outer;
  FE1D FEId_low;
  TFE1D *Element;
  TFE2D *ele;

  SecondDer = new boolean[n_fespaces];
// ########################################################################
// store information in local arrays
// ########################################################################
  BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
  N_BaseFuncts = TFEDatabase2D::GetN_BaseFunctFromFE2D();

  Coll = fespaces[0]->GetCollection(); // all spaces use same Coll
  GlobalNumbers = fespaces[0]->GetGlobalNumbers();
  BeginIndex = fespaces[0]->GetBeginIndex();
  u1 = fefunctions[0]->GetValues();
  u2 = fefunctions[1]->GetValues();

  GlobalNumbers_Outer = fespaces[1]->GetGlobalNumbers();
  BeginIndex_Outer = fespaces[1]->GetBeginIndex();
  V = SurfaceFeFunct->GetValues();

  Coll_low = fespaces_low[0]->GetCollection(); // all low spaces use same Coll
  N_Cells_low = Coll_low->GetN_Cells();
  BeginIndex_low =  fespaces_low[0]->GetBeginIndex();
  GlobalNumbers_low =  fespaces_low[0]->GetGlobalNumbers();

  RowPtr = sqmatrices_low[0]->GetRowPtr();
  KCol = sqmatrices_low[0]->GetKCol();

  ValuesA = sqmatrices_low[0]->GetEntries();
  ValuesM = sqmatrices_low[1]->GetEntries();
  RHS = RHSs[0];

  N_LocalUsedElements = n_fespaces;
  for(j=0;j<n_fespaces;j++)
    SecondDer[j]=FALSE;
  
  Bi = TDatabase::ParamDB->REACTOR_P13; 
  Da = TDatabase::ParamDB->REACTOR_P14;
  beta = TDatabase::ParamDB->REACTOR_P15;
  
  Pe_s = TDatabase::ParamDB->REACTOR_P17; // Peclet number
  if(Pe_s==0.)
   c0=0.;
  else
   c0 =  1./Pe_s;
//   AddedMass = 0.;

// ########################################################################
// loop over all low space cells
// ########################################################################
  for(i=0;i<N_Cells_low;i++)
  {
    N = Cell_array[i];
    Me = Coll->GetCell(N);
    IJoint = Joint_array[i];

    FEId = fespaces[0]->GetFE2D(N, Me);  // FEID of velocity space in the outer domain
    ele = TFEDatabase2D::GetFE2D(FEId);
    RefElement = TFEDatabase2D::GetRefElementFromFE2D(FEId);
    FeDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(FEId);
    N_JointDOF = FeDesc->GetN_JointDOF();
    JointDOF = FeDesc->GetJointDOF(IJoint);
    N_BaseFunct = FeDesc->GetN_DOF();
    DOF = GlobalNumbers + BeginIndex[N];

#ifdef __SOLUBLE__ 
//     N_Outer = Me->GetLocalCellNo();
    N_Outer = N; // free surf flows N=N_Outer;
    FEId_Outer = fespaces[1]->GetFE2D(N_Outer, Me);  // FEID of surfactant space in the outer domain
    FeDesc_Outer = TFEDatabase2D::GetFEDesc2DFromFE2D(FEId_Outer);
    N_BaseFunct_Outer = FeDesc_Outer->GetN_DOF();
    DOF_Outer = GlobalNumbers_Outer + BeginIndex_Outer[N_Outer];
#endif     
    DOF_LOW = GlobalNumbers_low + BeginIndex_low[i];
    Me_low = Coll_low->GetCell(i);
    FEId_low = fespaces_low[0]->GetFE1D(i, Me_low);
    Element = TFEDatabase2D::GetFE1D(FEId_low);
    N_BaseFunct_low = Element->GetN_DOF();
    if(N_JointDOF != N_BaseFunct_low )
     {
      cout<< " " << N_JointDOF <<" N_JointDOF != N_BaseFunct_low " << N_BaseFunct_low<<endl;
      exit(0);
    }

    memset(LocMatrixA, 0, N_BaseFunct_low*N_BaseFunct_low*SizeOfDouble);
    memset(LocMatrixM, 0, N_BaseFunct_low*N_BaseFunct_low*SizeOfDouble);
    memset(LocRhs, 0, N_BaseFunct_low*SizeOfDouble);

    l = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
    LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(2*l);
    qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
    qf1->GetFormulaData(N_LinePoints, LineWeights, zeta);

    TFEDatabase2D::GetBaseFunct2DFromFE2D(FEId)->MakeRefElementData(LineQuadFormula);
    TFEDatabase2D::GetBaseFunct2DFromFE2D(FEId_Outer)->MakeRefElementData(LineQuadFormula);

    ORDER = TFEDatabase2D::GetAccuracyFromFE2D(FEId);
    switch(RefElement)
      {
        case BFUnitSquare:
          RefTrans = QuadIsoparametric;
          F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
          ((TQuadIsoparametric *)F_K)->SetCell(Me);
          ((TQuadIsoparametric *)F_K)->GetOrigBoundFromRef(IJoint, N_LinePoints, zeta, X_B, Y_B);
        break;

        case BFUnitTriangle:
          RefTrans = TriaIsoparametric;
          F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
          ((TTriaIsoparametric *)F_K)->SetApproximationOrder(ORDER);  
          ((TTriaIsoparametric *)F_K)->SetCell(Me);
          ((TTriaIsoparametric *)F_K)->GetOrigBoundFromRef(IJoint, N_LinePoints, zeta, X_B, Y_B);
        break;
      } // endswitch

#ifdef __SOLUBLE__ 
      uref_Outer = TFEDatabase2D::GetJointValues2D(BaseFuncts[FEId_Outer],
                     LineQuadFormula, IJoint);
      uxiref_Outer = TFEDatabase2D::GetJointDerivatives2D(BaseFuncts[FEId_Outer],
                        LineQuadFormula, IJoint, D10);
      uetaref_Outer = TFEDatabase2D::GetJointDerivatives2D(BaseFuncts[FEId_Outer],
                        LineQuadFormula, IJoint, D01);
#endif  
      uref = TFEDatabase2D::GetJointValues2D(BaseFuncts[FEId],
                        LineQuadFormula, IJoint);
      uxiref = TFEDatabase2D::GetJointDerivatives2D(BaseFuncts[FEId],
                        LineQuadFormula, IJoint, D10);
      uetaref = TFEDatabase2D::GetJointDerivatives2D(BaseFuncts[FEId],
                        LineQuadFormula, IJoint, D01);

      for(k=0;k<N_LinePoints;k++)
        {
          F_K->GetTangent(IJoint, zeta[k], t0, t1);
          normn = sqrt(t0*t0+t1*t1);
          n0 =  t1/normn;
          n1 = -t0/normn;

          switch(RefElement)
          {
            case BFUnitSquare:
              ((TQuadIsoparametric *)F_K)->GetOrigValues(IJoint,  zeta[k],
                        N_BaseFunct, uref[k], uxiref[k], uetaref[k],
                        uorig, uxorig, uyorig);
#ifdef __SOLUBLE__       
              ((TQuadIsoparametric *)F_K)->GetOrigValues(IJoint,  zeta[k],
                        N_BaseFunct_Outer, uref_Outer[k], uxiref_Outer[k], uetaref_Outer[k],
                        uorig_Outer, uxorig_Outer, uyorig_Outer);
#endif        
            break;

            case BFUnitTriangle:
              ((TTriaIsoparametric *)F_K)->GetOrigValues(IJoint,  zeta[k],
                        N_BaseFunct, uref[k], uxiref[k], uetaref[k],
                        uorig, uxorig, uyorig);
#ifdef __SOLUBLE__ 	      
              ((TTriaIsoparametric *)F_K)->GetOrigValues(IJoint,  zeta[k],
                        N_BaseFunct_Outer, uref_Outer[k], uxiref_Outer[k], uetaref_Outer[k],
                        uorig_Outer, uxorig_Outer, uyorig_Outer);
#endif  

            break;
          } // endswitch

          r_axial = fabs(X_B[k]); // r value in the axial symmetric integral
         if(X_B[k]<=0)
          {
           cout <<"X_B[k] negative in Assemble1D change Quad rule " <<  X_B[k] <<endl;
//         exit(0);
          }
//           cout << " x " << r_axial<< " y " << Y_B[k]<< endl;

//       value of C at this integral point
         C = 0.;
#ifdef __SOLUBLE__  
         for(l=0;l<N_BaseFunct_Outer;l++)
           C += C_Outer[DOF_Outer[l]]*uorig_Outer[l];
#endif  
//           cout <<" C:  " << C <<endl;

//       get surfactant solution  on the interface
        U=0.;
#ifndef __MASSTRANSTEST__
        for(l=0;l<N_BaseFunct_low;l++)
         {
          local_j   = JointDOF[l];
          test00 = uorig[local_j];
          m = DOF_LOW[l];
          U  += V[m]*test00;
         }
#endif
         //         cout <<" U:  out " << U <<endl;

          //       get velocity gradients
          U1 = 0.;  u1x=0.; u2x=0.; u1y=0.; u2y=0.;
          for(l=0;l<N_BaseFunct;l++)
            {
             m = DOF[l];
             U1  += u1[m]*uorig[l];
             u1x += u1[m]*uxorig[l];
             u1y += u1[m]*uyorig[l];
             u2x += u2[m]*uxorig[l];
             u2y += u2[m]*uyorig[l];
            }

          TangDivU =  u1x - (u1x*n0 + u1y*n1)*n0  + U1/r_axial
                    + u2y - (u2x*n0 + u2y*n1)*n1;

//           cout <<" u1x:  " << u1x <<" u1y:  " << u1y <<endl;
//           cout <<" u2x:  " << u2x <<" u2y:  " << u2y <<endl;
//           cout <<" TangDivU:  " << TangDivU <<endl;

          Mult = sqrt(t0*t0+t1*t1)*(LineWeights[k]);

          rhsval = (beta/Da)*C*(1. - U)  -  Bi*U;

          rhsval *= Mult*r_axial;
          U *= Mult*r_axial;
          C *= Mult*r_axial;

          for(l=0;l<N_BaseFunct_low;l++)
           {
            local_j   = JointDOF[l];

            test00  = uorig[local_j];
            test10  = uxorig[local_j];
            test01  = uyorig[local_j];

            ngrad_test= n0*test10 + n1*test01;
            d1 = test10 - ngrad_test*n0;
            d2 = test01 - ngrad_test*n1;

//          rhs
            LocRhs[l] += rhsval*test00; // soluble surfactant relation explicit (in C) form
            
            for(m=0;m<N_BaseFunct_low;m++)
             {
              local_i   = JointDOF[m];

              ansatz00 = uorig[local_i];
              ansatz10 = uxorig[local_i];
              ansatz01 = uyorig[local_i];

//              cout << local_i << " -- " << local_j << endl;
              ngrad_ansatz= n0*ansatz10 + n1*ansatz01;
              e1 = ansatz10 - ngrad_ansatz*n0;
              e2 = ansatz01 - ngrad_ansatz*n1;

//              cout << " Tgrad . n  " << e1*n0 + e2*n1 << endl;
#ifndef __MASSTRANSTEST__
              val  =c0*(d1*e1 + d2*e2);
              val +=TangDivU*test00*ansatz00;
              val *= (Mult*r_axial);
              LocMatrixA[l*N_BaseFunct_low+m] += val;
#endif
              
              val  = test00*ansatz00;
              val *= (Mult*r_axial);
              LocMatrixM[l*N_BaseFunct_low+m] += val;
            }
          } //  for(l=0;l<N_Joint
        } //  for(k=0;k<N_

//       for(l=0;l<N_BaseFunct_low;l++)
//        for(m=0;m<N_BaseFunct_low;m++)
//         cout << " LocMatrixA " << LocMatrixA[l*N_BaseFunct_low+m]<<endl;

//   add to global matrices
    for(l=0;l<N_BaseFunct_low;l++)
     {
      TestDOF = DOF_LOW[l];
      RHS[TestDOF] += LocRhs[l]; // soluble surfactant relation

      begin = RowPtr[TestDOF];
      end = RowPtr[TestDOF+1];
      for(n=begin;n<end;n++)
       {
       for(m=0;m<N_BaseFunct_low;m++)
        {
         if(KCol[n] == DOF_LOW[m])
          {
           ValuesA[n] +=LocMatrixA[l*N_BaseFunct_low+m];
           ValuesM[n] +=LocMatrixM[l*N_BaseFunct_low+m];
           break;
          }
        } // for(m=0;m<N_BaseFunct_low
      } // for(n=begin;n<end;n++)
     } // for(l=0;l<N_BaseFunct_low
    } //  for(i=0;i<N_Cells_low

//   cout << " AddedMass on Interface; " << 2*Pi*AddedMass <<endl;
//   cout << " GammaMass at the Interface; " << 2*Pi*GammaMass <<endl;
//   cout << " CMass at the Interface; " << 2*Pi*CMass <<endl;
 delete [] SecondDer;
// exit(0);
}

void  AssembleSurf1D_SolubleSurfact_Implicit(int n_fespaces, TFESpace2D **fespaces, 
                     TFEFunction2D **fefunctions, int N_FESpaces_low, TFESpace1D **fespaces_low,
                     TFEFunction1D *SurfaceFeFunct, int N_SquareMatrices,
                     TSquareMatrix1D **sqmatrices_low, int N_Rhs, double **RHSs, 
                     TFESpace1D **ferhs_low, int *Cell_array, int *Joint_array, double *C_Outer)
{
  int i, j, k, l, m, n, N_Cells_low, N, N_LocalUsedElements, local_i, local_j, ORDER;
  int N_BaseFunct, N_BaseFunct_low,  N_Points, N_JointDOF, begin, end, *N_BaseFuncts;
  int *BeginIndex_low, *GlobalNumbers_low, *DOF, *DOF_LOW, TestDOF, AnsatzDOF, IJoint ;
  int *BeginIndex, *GlobalNumbers, *GlobalNumbers_Outer, *BeginIndex_Outer, N_BaseFunct_Outer;
  int LocN_BF[N_BaseFuncts2D], N_LinePoints, *KCol, *RowPtr, *JointDOF, N_Outer;
  int *DOF_Outer;

  double x0, y0, x1, y1, t0, t1, n0, n1, normn;
  double AbsDetjk[MaxN_QuadPoints_2D], Mult;
  double *weights, *xi, *eta;
  double **uref, **uxiref, **uetaref;
  double **uref_Outer, **uxiref_Outer, **uetaref_Outer,  C, CX, CY;
  double *LineWeights, *zeta, *ValuesA, *ValuesM;
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];
  double uorig[MaxN_BaseFunctions2D], uxorig[MaxN_BaseFunctions2D], uyorig[MaxN_BaseFunctions2D];
  double uorig_Outer[MaxN_BaseFunctions2D], uxorig_Outer[MaxN_BaseFunctions2D];
  double uyorig_Outer[MaxN_BaseFunctions2D];
  double c0, r2;
  double Pr = TDatabase::ParamDB->PR_NR;
  double val, rhsval, theta, ngrad_ansatz, ngrad_test, TangDivU;
  double  X_B[100], Y_B[100], r_axial, d1, d2, e1, e2;
  double LocMatrixA[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double LocMatrixM[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double LocRhs[MaxN_BaseFunctions2D];
  double ansatz00, ansatz10, ansatz01 ;
  double test00, test10, test01, *u1, *u2, u1x, u2x, u1y, u2y, U1;
  double *RHS, *V, U, AddedMass, CMass, GammaMass;
  double Bi, Da, beta, Pe_s, NGrad_C;
  
  BaseFunct2D LocBF[N_BaseFuncts2D];
  BaseFunct2D *BaseFuncts;
  boolean *SecondDer;
  BF2DRefElements RefElement;
  RefTrans2D RefTrans;
  TRefTrans2D *F_K;
  QuadFormula1D LineQuadFormula;
  TQuadFormula1D *qf1;
  TFEDesc2D *FeDesc, *FeDesc_Outer;
  TFEDesc1D *FeDesc_low;
  TCollection *Coll, *Coll_low;
  TBaseCell *Me, *Me_low;
  FE2D FEId, FEId_Outer;
  FE1D FEId_low;
  TFE1D *Element;
  TFE2D *ele;

  SecondDer = new boolean[n_fespaces];
// ########################################################################
// store information in local arrays
// ########################################################################
  BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
  N_BaseFuncts = TFEDatabase2D::GetN_BaseFunctFromFE2D();

  Coll = fespaces[0]->GetCollection(); // all spaces use same Coll
  GlobalNumbers = fespaces[0]->GetGlobalNumbers();
  BeginIndex = fespaces[0]->GetBeginIndex();
  u1 = fefunctions[0]->GetValues();
  u2 = fefunctions[1]->GetValues();

  GlobalNumbers_Outer = fespaces[1]->GetGlobalNumbers();
  BeginIndex_Outer = fespaces[1]->GetBeginIndex();

  V = SurfaceFeFunct->GetValues();

  Coll_low = fespaces_low[0]->GetCollection(); // all low spaces use same Coll
  N_Cells_low = Coll_low->GetN_Cells();
  BeginIndex_low =  fespaces_low[0]->GetBeginIndex();
  GlobalNumbers_low =  fespaces_low[0]->GetGlobalNumbers();

  RowPtr = sqmatrices_low[0]->GetRowPtr();
  KCol = sqmatrices_low[0]->GetKCol();

  ValuesA = sqmatrices_low[0]->GetEntries();
  ValuesM = sqmatrices_low[1]->GetEntries();
  RHS = RHSs[0];

  N_LocalUsedElements = n_fespaces;
  for(j=0;j<n_fespaces;j++)
    SecondDer[j]=FALSE;
  
  Bi = TDatabase::ParamDB->REACTOR_P13; 
  Da = TDatabase::ParamDB->REACTOR_P14;
  beta = TDatabase::ParamDB->REACTOR_P15;
  
  Pe_s = TDatabase::ParamDB->REACTOR_P17; // Peclet number
  if(Pe_s==0.)
   c0=0.;
  else
   c0 =  1./Pe_s;

//   AddedMass = 0.;
//   CMass =0.;
// ==================================================================================
// loop over all low space cells
// ==================================================================================
  for(i=0;i<N_Cells_low;i++)
  {
    N = Cell_array[i];
    Me = Coll->GetCell(N);
    IJoint = Joint_array[i];

    FEId = fespaces[0]->GetFE2D(N, Me);  // FEID of velocity space in the outer domain
    ele = TFEDatabase2D::GetFE2D(FEId);
    RefElement = TFEDatabase2D::GetRefElementFromFE2D(FEId);
    FeDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(FEId);
    N_JointDOF = FeDesc->GetN_JointDOF();
    JointDOF = FeDesc->GetJointDOF(IJoint);
    N_BaseFunct = FeDesc->GetN_DOF();
    DOF = GlobalNumbers + BeginIndex[N];

    
//     if(N==1)
//     {
//      for(l=0;l<N_BaseFunct;l++)
//      {
//       cout << DOF[l] << " N_BaseFunct " << N_BaseFunct << " u2 " << u2[DOF[l]] ;
//       fespaces[0]->GetDOFPosition(DOF[l], x0, y0);
//        cout << " x " << x0 << " y " << y0 << endl;
//      }
//     }
// exit(0);      
    
//     N_Outer = Me->GetLocalCellNo();
    N_Outer = N; // free surf flows N=N_Outer;    
    FEId_Outer = fespaces[1]->GetFE2D(N_Outer, Me);  // FEID of surfactant space in the outer domain
    FeDesc_Outer = TFEDatabase2D::GetFEDesc2DFromFE2D(FEId_Outer);
    N_BaseFunct_Outer = FeDesc_Outer->GetN_DOF();
    DOF_Outer = GlobalNumbers_Outer + BeginIndex_Outer[N_Outer];

    DOF_LOW = GlobalNumbers_low + BeginIndex_low[i];
    Me_low = Coll_low->GetCell(i);
    FEId_low = fespaces_low[0]->GetFE1D(i, Me_low);
    Element = TFEDatabase2D::GetFE1D(FEId_low);
    N_BaseFunct_low = Element->GetN_DOF();
    if(N_JointDOF != N_BaseFunct_low )
     {
      cout<< " " << N_JointDOF <<" N_JointDOF != N_BaseFunct_low " << N_BaseFunct_low<<endl;
      exit(0);
    }

    memset(LocMatrixA, 0, N_BaseFunct_low*N_BaseFunct_low*SizeOfDouble);
    memset(LocMatrixM, 0, N_BaseFunct_low*N_BaseFunct_low*SizeOfDouble);
    memset(LocRhs, 0, N_BaseFunct_low*SizeOfDouble);

    l = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
    LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(2*l);
    qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
    qf1->GetFormulaData(N_LinePoints, LineWeights, zeta);

    TFEDatabase2D::GetBaseFunct2DFromFE2D(FEId)->MakeRefElementData(LineQuadFormula);
    TFEDatabase2D::GetBaseFunct2DFromFE2D(FEId_Outer)->MakeRefElementData(LineQuadFormula);

    ORDER = TFEDatabase2D::GetAccuracyFromFE2D(FEId);
    switch(RefElement)
      {
        case BFUnitSquare:
          RefTrans = QuadIsoparametric;
          F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
          ((TQuadIsoparametric *)F_K)->SetCell(Me);
          ((TQuadIsoparametric *)F_K)->GetOrigBoundFromRef(IJoint, N_LinePoints, zeta, X_B, Y_B);
        break;

        case BFUnitTriangle:
          RefTrans = TriaIsoparametric;
          F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
          ((TTriaIsoparametric *)F_K)->SetApproximationOrder(ORDER);  
          ((TTriaIsoparametric *)F_K)->SetCell(Me);
          ((TTriaIsoparametric *)F_K)->GetOrigBoundFromRef(IJoint, N_LinePoints, zeta, X_B, Y_B);
        break;
      } // endswitch

#ifdef __SOLUBLE__ 
      uref_Outer = TFEDatabase2D::GetJointValues2D(BaseFuncts[FEId_Outer],
                     LineQuadFormula, IJoint);
      uxiref_Outer = TFEDatabase2D::GetJointDerivatives2D(BaseFuncts[FEId_Outer],
                        LineQuadFormula, IJoint, D10);
      uetaref_Outer = TFEDatabase2D::GetJointDerivatives2D(BaseFuncts[FEId_Outer],
                        LineQuadFormula, IJoint, D01);

 #endif          
      

      uref = TFEDatabase2D::GetJointValues2D(BaseFuncts[FEId],
                        LineQuadFormula, IJoint);
      uxiref = TFEDatabase2D::GetJointDerivatives2D(BaseFuncts[FEId],
                        LineQuadFormula, IJoint, D10);
      uetaref = TFEDatabase2D::GetJointDerivatives2D(BaseFuncts[FEId],
                        LineQuadFormula, IJoint, D01);

      for(k=0;k<N_LinePoints;k++)
        {
          F_K->GetTangent(IJoint, zeta[k], t0, t1);
          normn = sqrt(t0*t0+t1*t1);
          n0 =  t1/normn;
          n1 = -t0/normn;

          switch(RefElement)
          {
            case BFUnitSquare:
              ((TQuadIsoparametric *)F_K)->GetOrigValues(IJoint,  zeta[k],
                        N_BaseFunct, uref[k], uxiref[k], uetaref[k],
                        uorig, uxorig, uyorig);
#ifdef __SOLUBLE__       
              ((TQuadIsoparametric *)F_K)->GetOrigValues(IJoint,  zeta[k],
                        N_BaseFunct_Outer, uref_Outer[k], uxiref_Outer[k], uetaref_Outer[k],
                        uorig_Outer, uxorig_Outer, uyorig_Outer);
#endif      
            break;

            case BFUnitTriangle:
              ((TTriaIsoparametric *)F_K)->GetOrigValues(IJoint,  zeta[k],
                        N_BaseFunct, uref[k], uxiref[k], uetaref[k],
                        uorig, uxorig, uyorig);
#ifdef __SOLUBLE__ 	      
              ((TTriaIsoparametric *)F_K)->GetOrigValues(IJoint,  zeta[k],
                        N_BaseFunct_Outer, uref_Outer[k], uxiref_Outer[k], uetaref_Outer[k],
                        uorig_Outer, uxorig_Outer, uyorig_Outer);
#endif   
            break;
          } // endswitch
         
          r_axial = fabs(X_B[k]); // r value in the axial symmetric integral
         if(X_B[k]<=0)
          {
           cout <<"X_B[k] negative in Assemble1D change Quad rule " <<  X_B[k] <<endl;
//         exit(0);
          }
//           cout << " x " << r_axial<< " y " << Y_B[k]<< endl;

         // value of C at this integral point
         C = 0.;
#ifdef __SOLUBLE__ 
         for(l=0;l<N_BaseFunct_Outer;l++)
           C += C_Outer[DOF_Outer[l]]*uorig_Outer[l];
#endif
          // cout <<" C:  " << C <<endl;

          // get velocity gradients
          U1 = 0.;  u1x=0.; u2x=0.; u1y=0.; u2y=0.;
          for(l=0;l<N_BaseFunct;l++)
            {
             m = DOF[l];
             U1  += u1[m]*uorig[l];
             u1x += u1[m]*uxorig[l];
             u1y += u1[m]*uyorig[l];
             u2x += u2[m]*uxorig[l];
             u2y += u2[m]*uyorig[l];  
            }

          TangDivU =  u1x - (u1x*n0 + u1y*n1)*n0  + U1/r_axial
                    + u2y - (u2x*n0 + u2y*n1)*n1;

//           cout <<" u1x:  " << u1x <<" u1y:  " << u1y <<endl;
//           cout <<" u2x:  " << u2x <<" u2y:  " << u2y <<endl;
//           cout << u2x <<" TangDivU:  " <<    U2 <<endl;


//   	 if(fabs(TangDivU)>1e3)
// 	  {
//            OutPut(" C:  " << C <<" TangDivU:  " << TangDivU <<endl);
//             for(l=0;l<N_BaseFunct;l++)
// 	     {
//               OutPut("Cell " << N << " l " << l << " Dof " <<  DOF[l] << " u1:  " <<u1[DOF[l]] << " u2:  " <<u2[DOF[l]]<< 
//                      " uorig:  " <<uorig[l] << " uxorig:  " << uxorig[l] <<" uyorig:  " << uyorig[l]<<endl);
// 	   
// 	      
// 	     }
	     
// 	     OutPut( " u1x:  " << u1x <<" u1y:  " << u1y 
// 	          <<" u2x:  " << u2x <<" u2y:  " << u2y <<endl);
	     
//            Me->GetVertex(IJoint)->GetCoords(x0, y0);       
//            Me->GetVertex((IJoint+1) % 3)->GetCoords(x1, y1);    
//            cout << IJoint <<  " X  " << x0<<" Y  " << y0<<" X  " << x1<<" Y  " << y1 <<endl;
// 	   exit(0);
	   
// 	  }
    
          Mult = sqrt(t0*t0+t1*t1)*LineWeights[k];

          rhsval =   (beta/Da)*C;
          rhsval *= Mult*r_axial;
//           AddedMass +=rhsval;
//           CMass +=Mult*C;

          for(l=0;l<N_BaseFunct_low;l++)
           {
            local_j   = JointDOF[l];

            test00  = uorig[local_j];
            test10  = uxorig[local_j];
            test01  = uyorig[local_j];

            ngrad_test= n0*test10 + n1*test01;
            d1 = test10 - ngrad_test*n0;
            d2 = test01 - ngrad_test*n1;

//          rhs
            LocRhs[l] += rhsval*test00; 

            for(m=0;m<N_BaseFunct_low;m++)
             {
              local_i   = JointDOF[m];

              ansatz00 = uorig[local_i];
              ansatz10 = uxorig[local_i];
              ansatz01 = uyorig[local_i];

//              cout << local_i << " -- " << local_j << endl;
              ngrad_ansatz= n0*ansatz10 + n1*ansatz01;
              e1 = ansatz10 - ngrad_ansatz*n0;
              e2 = ansatz01 - ngrad_ansatz*n1;

#ifndef __MASSTRANSTEST__
//              cout << " Tgrad . n  " << e1*n0 + e2*n1 << endl;
              val = c0*(d1*e1 + d2*e2);
              val += ((beta/Da)*C  +  Bi)*ansatz00*test00;
              val +=TangDivU*test00*ansatz00;
              val *= (Mult*r_axial);
              LocMatrixA[l*N_BaseFunct_low+m] += val;
#endif
              
              val  = test00*ansatz00;
              val *= (Mult*r_axial);
              LocMatrixM[l*N_BaseFunct_low+m] += val;
            }
          } //  for(l=0;l<N_Joint
        } //  for(k=0;k<N_

//       if(i==18 &&  TDatabase::ParamDB->P15==1)
//       {
//       for(l=0;l<N_BaseFunct_low;l++)
//        for(m=0;m<N_BaseFunct_low;m++)
//         cout << i << " LocMatrixA " << LocMatrixA[l*N_BaseFunct_low+m]<<endl;
//         exit(0);
//       }
       
//   add to global matrices
    for(l=0;l<N_BaseFunct_low;l++)
     {
      TestDOF = DOF_LOW[l];
      RHS[TestDOF] += LocRhs[l]; // soluble surfactant relation

      begin = RowPtr[TestDOF];
      end = RowPtr[TestDOF+1];
      for(n=begin;n<end;n++)
       {
       for(m=0;m<N_BaseFunct_low;m++)
        {
         if(KCol[n] == DOF_LOW[m])
          {
           ValuesA[n] +=LocMatrixA[l*N_BaseFunct_low+m];
           ValuesM[n] +=LocMatrixM[l*N_BaseFunct_low+m];
 
//            cout << TestDOF  << ", " << DOF_LOW[m] <<" M: " 
//            << LocMatrixM[l*N_BaseFunct_low+m] <<endl;   
           break;
          }
        } // for(m=0;m<N_BaseFunct_low
      } // for(n=begin;n<end;n++)
     } // for(l=0;l<N_BaseFunct_low
    } //  for(i=0;i<N_Cells_low

//   cout << " AddedMass on Interface; " << 2*Pi*AddedMass <<endl;
//   cout << " GammaMass at the Interface; " << 2*Pi*GammaMass <<endl;
//   cout << " CMass at the Interface; " << 2*Pi*CMass <<endl;
 delete [] SecondDer;
// exit(0);
}



void Getcellangle(TFESpace2D *Space, double *MinMaxAngle)
{
 int i,j,k,l, N_Cells, N_Edges;
 int found,  N_LinePoints;

 double TX[4], TY[4], hE[4], Theta, tx, ty, Test, MQI=0.;
 TBaseCell *cell;
 FE2D FEId;
 BF2DRefElements RefElement;
 TRefTrans2D *F_K;
 RefTrans2D RefTrans;
 TCollection *Cells;

  MinMaxAngle[0] = 180;  // Min_Angel = 180
  MinMaxAngle[1] = 0;  // Max_Angel = 0
  Cells = Space->GetCollection();
  N_Cells = Cells->GetN_Cells();
     
//      TX      = new double[4];  // Max no edges in 2d
//      TY      = new double[4];  // Max no edges in 2d
//      hE      = new double[4];  // Max no edges in 2d

  for(i=0;i<N_Cells;i++)
   {
     cell    = Cells->GetCell(i);
     N_Edges = cell->GetN_Edges();

     FEId = Space->GetFE2D(i, cell);
     RefElement = TFEDatabase2D::GetRefElementFromFE2D(FEId);

     switch(RefElement)
        {
         case BFUnitTriangle:

            RefTrans = TriaAffin;
            F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
            ((TTriaAffin *)F_K)->SetCell(cell);

          break;

          case BFUnitSquare:

            RefTrans = QuadAffin;
            F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
            ((TQuadAffin *)F_K)->SetCell(cell);

          break;

          default:
            Error("only triangles and quadrilaterals are allowed" << endl);
            Error("file: " << __FILE__ << " line " << __LINE__ << endl);
            exit(-1);
          } // endswitch

     for(j=0;j<N_Edges;j++)
      {
        F_K->GetTangent(j, 0, tx, ty);
        TX[j] = tx;
        TY[j] = ty;
        hE[j] = sqrt(tx*tx+ty*ty);

   // cout <<"cell : " <<i << "  j= " << j << ": " <<TX[j]<< "------ " << TY[j] << endl;
       } // endfor j

//      Test = 0;
      k = N_Edges -1;
      for(j=0;j<N_Edges;j++)
      {
       if(hE[j]==0.0 || hE[k]== 0.0 )
        Theta = 0.0;
       else
        Theta = acos(-(TX[j]*TX[k]+TY[j]*TY[k])/(hE[j]*hE[k]))*(180/3.141592654);

       k = j;
//        Test +=Theta;
       if(MinMaxAngle[0]>Theta) MinMaxAngle[0] = Theta;
       if(MinMaxAngle[1]<Theta) MinMaxAngle[1] = Theta;
//        cout <<"cell : " <<i << "  j= " << j << ": " << " Theta : " << Theta << endl;
//  *****************************************************
//  Grid test

      MQI += (60. - Theta)*(60. - Theta);
//  *****************************************************

     }
//       cout <<"cell : " <<i <<  " sum of 3 angels : " << Test << endl;
     //  cout<<endl;

   } // endfor i

   MQI /=double(3*N_Cells);
   MQI = sqrt(MQI);

// OutPut("Mesh Quality Indicator: "<< MQI<< endl);
//    delete [] TX;
//    delete [] TY;
//    delete [] hE;
 //cout<< " Min_Angel: "<< MinMaxAngle[0]<< "  Max_Angel : "<<MinMaxAngle[1]<< endl;
// exit(0);
}

void Get_KE(TFEVectFunct2D *Velocity, double *parameters)
 {
  int i,j,k,l, polydegree, Phase_No;
  int N_Cells, N_Joints, N_Vertices;
  TBaseCell *cell;
  TCollection *coll;
  int *BeginIndex, *GlobalNumbers, *DOF;
  TJoint *joint;
  double KE=0., volume=0., volume_1=0., volume_2=0.;
  double KE_QP, KE_1=0., x_mass=0., y_mass=0.;
  double u1_rise, u2_rise, U1_Rise=0., U2_Rise=0.;
  TRefTrans2D *F_K;
  double U1, U2;
  FE2D FEid;
  TBaseFunct2D *bf;
  int N_QFPoints;
  double *weights, *xi, *eta;
  double values[MaxN_QuadPoints_2D][MaxN_BaseFunctions2D];
  double AbsDetjk[MaxN_QuadPoints_2D], X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];
  double *ValuesVX, *ValuesVY;
  int N_BF;
  double Mult, r_axial;
  TFESpace2D *VelocitySpace;
  JointType jointtype;
  BoundTypes bdtype;
  RefTrans2D RefTrans;
  boolean IsIsoparametric;
  QuadFormula2D QuadFormula;
  TQuadFormula2D *qf2;

  VelocitySpace = Velocity->GetFESpace2D();
  BeginIndex = VelocitySpace->GetBeginIndex();
  GlobalNumbers = VelocitySpace->GetGlobalNumbers();
  ValuesVX = Velocity->GetValues();
  ValuesVY = ValuesVX + Velocity->GetLength();

  coll = VelocitySpace->GetCollection();
  N_Cells = coll->GetN_Cells();

  for(i=0;i<N_Cells;i++)
   {
    cell = coll->GetCell(i);
    FEid = VelocitySpace->GetFE2D(i, cell);

    RefTrans = TFEDatabase2D::GetRefTrans2D_IDFromFE2D(FEid);
    N_Joints = cell->GetN_Joints();

    IsIsoparametric = FALSE;

    if(TDatabase::ParamDB->USE_ISOPARAMETRIC)
     {
     for(j=0;j<N_Joints;j++)
      {
       joint = cell->GetJoint(j);
       jointtype = joint->GetType();
       if(jointtype == BoundaryEdge)
        {
         bdtype = ((TBoundEdge *)(joint))->GetBoundComp()->GetType();
         if(bdtype != Line)  IsIsoparametric = TRUE;
        }
       if(jointtype == InterfaceJoint)
        {
          bdtype = ((TInterfaceJoint *)(joint))->GetBoundComp()->GetType();
          if(bdtype != Line)
            IsIsoparametric = TRUE;
        }
        if(jointtype == IsoInterfaceJoint || jointtype == IsoBoundEdge)
         IsIsoparametric = TRUE;

      } // for(j=0;j<
     } // if(TDatabase::ParamDB->USE_ISOPARAMETRIC)

   if(IsIsoparametric)
    {
      switch(N_Joints)
      {
        case 4:
          RefTrans = QuadIsoparametric;
        break;

        case 3:
          RefTrans = TriaIsoparametric;
        break;
      }
    } // endif IsIsoparametric

    F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);

    switch(RefTrans)
    {
      case TriaAffin:
        polydegree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEid);
        QuadFormula = TFEDatabase2D::GetQFTriaFromDegree(9);
        qf2 = TFEDatabase2D::GetQuadFormula2D(QuadFormula);
        qf2->GetFormulaData(N_QFPoints, weights, xi, eta);
        ((TTriaAffin *)F_K)->SetCell(cell);
        ((TTriaAffin *)F_K)->GetOrigFromRef(N_QFPoints, xi, eta,
                                         X, Y, AbsDetjk);
      break;

      case TriaIsoparametric:
        polydegree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEid);
        QuadFormula = TFEDatabase2D::GetQFTriaFromDegree(9);
        qf2 = TFEDatabase2D::GetQuadFormula2D(QuadFormula);
        qf2->GetFormulaData(N_QFPoints, weights, xi, eta);
        ((TTriaIsoparametric *)F_K)->SetApproximationOrder(polydegree);
        ((TTriaIsoparametric *)F_K)->SetQuadFormula(QuadFormula);
        ((TTriaIsoparametric *)F_K)->SetCell(cell);
        ((TTriaIsoparametric *)F_K)->GetOrigFromRef(N_QFPoints, xi, eta,
                                         X, Y, AbsDetjk);
      break;

      case QuadAffin:
        polydegree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEid);
        QuadFormula = TFEDatabase2D::GetQFQuadFromDegree(3*polydegree);
        qf2 = TFEDatabase2D::GetQuadFormula2D(QuadFormula);
        qf2->GetFormulaData(N_QFPoints, weights, xi, eta);
        ((TQuadAffin *)F_K)->SetCell(cell);
        ((TQuadAffin *)F_K)->GetOrigFromRef(N_QFPoints, xi, eta,
                                         X, Y, AbsDetjk);
      break;

      case QuadBilinear:
        polydegree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEid);
        QuadFormula = TFEDatabase2D::GetQFQuadFromDegree(3*polydegree);
        qf2 = TFEDatabase2D::GetQuadFormula2D(QuadFormula);
        qf2->GetFormulaData(N_QFPoints, weights, xi, eta);
        ((TQuadBilinear *)F_K)->SetCell(cell);
        ((TQuadBilinear *)F_K)->GetOrigFromRef(N_QFPoints, xi, eta,
                                         X, Y, AbsDetjk);
      break;

      case QuadIsoparametric:
        polydegree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEid);
        QuadFormula = TFEDatabase2D::GetQFQuadFromDegree(3*polydegree);
        qf2 = TFEDatabase2D::GetQuadFormula2D(QuadFormula);
        qf2->GetFormulaData(N_QFPoints, weights, xi, eta);
        ((TQuadIsoparametric *)F_K)->SetApproximationOrder(polydegree);
        ((TQuadIsoparametric *)F_K)->SetQuadFormula(QuadFormula);
        ((TQuadIsoparametric *)F_K)->SetCell(cell);
        ((TQuadIsoparametric *)F_K)->GetOrigFromRef(N_QFPoints, xi, eta,
                                         X, Y, AbsDetjk);
      break;
    }

    // find basis functions on cell i
    bf = TFEDatabase2D::GetBaseFunct2DFromFE2D(FEid);
    N_BF = bf->GetDimension();
    DOF = GlobalNumbers + BeginIndex[i];

    for(k=0;k<N_QFPoints;k++)
     {
      bf->GetDerivatives(D00, xi[k], eta[k], values[k]);
      r_axial = fabs(X[k]);
      if(X[k]<=0)
       {
        cout <<"X[k] negative in Get_KE change Quad rule " <<  X[k] <<endl;
//         exit(0);
       }

      Mult = r_axial*weights[k]*AbsDetjk[k];

      KE_QP = 0.;
      u1_rise = 0.;
      u2_rise = 0.;

      for(l=0;l<N_BF;l++)
       {
        j = DOF[l];
        U1 = ValuesVX[j];
        U2 = ValuesVY[j];


 //       u.u
        KE_QP += (values[k][l]*U1 * values[k][l]*U1
	        + values[k][l]*U2 * values[k][l]*U2 );
        u1_rise += values[k][l]*U1;
        u2_rise += values[k][l]*U2;
       }

     KE += (KE_QP*Mult);
     volume += Mult;
     x_mass   += (X[k]*Mult);
     y_mass   += (Y[k]*Mult);
     U1_Rise  += (u1_rise * Mult);
     U2_Rise  += (u2_rise * Mult);

    }  // for(k=0;k<N_QF
   } // endfor i

    KE      /=volume;
    x_mass  /=volume;
    y_mass  /=volume;
    U1_Rise /=volume;
    U2_Rise /=volume;

 parameters[0] = volume;
 parameters[1] = KE;
 parameters[2] = x_mass;
 parameters[3] = y_mass;
 parameters[4] = U1_Rise;
 parameters[5] = U2_Rise;

//   cout<< " Volume_1: "<< volume_1<< " Volume_2: "<< volume_2
//       << " Volume: "<< volume<< endl;
//   cout<< " x_mass: "<< x_mass << " y_mass: "<< y_mass<< endl;
//   cout<< " KE_1: "<< KE_1<< " KE: "<< KE<< endl;
//   cout<< " U1_Rise: "<< U1_Rise << " U2_Rise: "<< U2_Rise<< endl;
//   exit(0);
 }

void GetSurfactMass(TFEFunction2D *fefunction, TFEFunction1D *fefunct_low,
                   int *Cell_array, int *Joint_array, double *errors)
{
  int i, j, k, l, m, n, N, N_Cells_low, N_LocalUsedElements, local_i, local_dof, ORDER;
  int N_BaseFunct,  N_Points, N_JointDOF, begin, end, *N_BaseFuncts, IJoint;
  TCollection *Coll, *Coll_low;
  TBaseCell *Me, *Me_low;
  FE2D FEId;
  TFE2D *ele;
  FE1D FEId_low;
  TFE1D *Element;
  TFESpace2D *fespace;
  TFESpace1D *fespace_low;
  BaseFunct2D LocBF[N_BaseFuncts2D];
  BaseFunct2D *BaseFuncts;
  int *DOF, *JointDOF, *BeginIndex, *GlobalNumbers, N_BaseFunct_low;
  TFEDesc2D *FeDesc;
  TFEDesc1D *FeDesc_low;
  BF2DRefElements RefElement;
  RefTrans2D RefTrans;
  TRefTrans2D *F_K;
  QuadFormula1D LineQuadFormula;
  TQuadFormula1D *qf1;
  int N_LinePoints;
  double *LineWeights, t0, t1, normn, n0, n1, X_P, Y_P;
  double *u, **uref, **uxiref, **uetaref, **uzetaref;
  double *Weights, *p1, *p2, *zeta, LocL2U, LocH1U,  Surf[100], Surf_DOF[10], Surf_val[10];
  double r_axial;
  double  X_B[MaxN_QuadPoints_2D], Y_B[MaxN_QuadPoints_2D], Exact_Surf[10], v, Mult;
  double uorig[MaxN_BaseFunctions2D], uxorig[MaxN_BaseFunctions2D];
  double uyorig[MaxN_BaseFunctions2D];
  double test00, test10, test01;
  double U, ux, uy, d1, d2, ngrad, h_K =0., h_K_min=1e8, h_K_max=-1e8;

  BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
  N_BaseFuncts = TFEDatabase2D::GetN_BaseFunctFromFE2D();

  fespace=fefunction->GetFESpace2D();
  Coll = fespace->GetCollection();
  BeginIndex = fespace->GetBeginIndex();
  GlobalNumbers = fespace->GetGlobalNumbers();
  u = fefunction->GetValues();

  fespace_low=fefunct_low->GetFESpace1D();
  Coll_low = fespace_low->GetCollection();
  N_Cells_low = Coll_low->GetN_Cells();

  errors[0] =0.; // surfactant mass
  errors[1] =0.; // surface area

// ########################################################################
// loop over all surf cells
// ########################################################################
  for(i=0;i<N_Cells_low;i++)
  {
    h_K =0.;
    Me_low = Coll_low->GetCell(i);
    FEId_low = fespace_low->GetFE1D(i, Me_low);
    Element = TFEDatabase2D::GetFE1D(FEId_low);
    N_BaseFunct_low = Element->GetN_DOF();

    N = Cell_array[i];
    Me = Coll->GetCell(N);
    IJoint = Joint_array[i];

    FEId = fespace->GetFE2D(N, Me);
    ele = TFEDatabase2D::GetFE2D(FEId);

    RefElement = TFEDatabase2D::GetRefElementFromFE2D(FEId);
    FeDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(FEId);
    N_JointDOF =  FeDesc->GetN_JointDOF();
    JointDOF = FeDesc->GetJointDOF(IJoint);
    N_BaseFunct = FeDesc->GetN_DOF();
    DOF = GlobalNumbers + BeginIndex[N];

    l = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
    ORDER = TFEDatabase2D::GetAccuracyFromFE2D(FEId);
    LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(2*l);
    qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
    qf1->GetFormulaData(N_LinePoints, LineWeights, zeta);

    switch(RefElement)
     {
      case BFUnitSquare:
          RefTrans = QuadIsoparametric;
          F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
          ((TQuadIsoparametric *)F_K)->SetCell(Me);
          ((TQuadIsoparametric *)F_K)->GetOrigBoundFromRef(IJoint,
	                            N_LinePoints, zeta, X_B, Y_B);
      break;

      case BFUnitTriangle:
          RefTrans = TriaIsoparametric;
          F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
          ((TTriaIsoparametric *)F_K)->SetApproximationOrder(ORDER);	  
          ((TTriaIsoparametric *)F_K)->SetCell(Me);
          ((TTriaIsoparametric *)F_K)->GetOrigBoundFromRef(IJoint, N_LinePoints, zeta, X_B, Y_B);

      break;
    } // endswitch

    TFEDatabase2D::GetBaseFunct2DFromFE2D(FEId)
                  ->MakeRefElementData(LineQuadFormula);


      uref = TFEDatabase2D::GetJointValues2D(BaseFuncts[FEId],
                        LineQuadFormula, IJoint);
      uxiref = TFEDatabase2D::GetJointDerivatives2D(BaseFuncts[FEId],
                        LineQuadFormula, IJoint, D10);
      uetaref = TFEDatabase2D::GetJointDerivatives2D(BaseFuncts[FEId],
                        LineQuadFormula, IJoint, D01);


    for(k=0;k<N_LinePoints;k++)
     {

         switch(RefElement)
          {
            case BFUnitSquare:
              ((TQuadIsoparametric *)F_K)->SetCell(Me);
              ((TQuadIsoparametric *)F_K)->GetOrigValues(IJoint,  zeta[k],
                        N_BaseFunct, uref[k], uxiref[k], uetaref[k],
                        uorig, uxorig, uyorig);
              ((TQuadIsoparametric *) F_K)->GetTangent(IJoint, zeta[k], t0, t1);
            break;

            case BFUnitTriangle:
              ((TTriaIsoparametric *)F_K)->SetApproximationOrder(ORDER);
              ((TTriaIsoparametric *)F_K)->SetCell(Me); 
              ((TTriaIsoparametric *)F_K)->GetOrigValues(IJoint,  zeta[k],
                        N_BaseFunct, uref[k], uxiref[k], uetaref[k],
                        uorig, uxorig, uyorig);
              ((TTriaIsoparametric *) F_K)->GetTangent(IJoint, zeta[k], t0, t1);
            break;
          } // endswitch


      normn = sqrt(t0*t0+t1*t1);
      n0 =  t1/normn;
      n1 = -t0/normn;

//  //  Projection of all points to the free boundary in normal direction  !! only for circle!!!
//       v =  sqrt(X_B[k]* X_B[k] + Y_B[k]*Y_B[k] );// /r;
      X_P =  X_B[k]  ;
      Y_P =  Y_B[k]  ;
      r_axial = fabs(X_P);

//       get solution 
     U=0.;
     for(l=0;l<N_BaseFunct_low;l++)
      {
       local_dof   = JointDOF[l];

       test00 = uorig[local_dof];

       m = DOF[local_dof];
       U  += u[m]*test00;
      }
// cout << "ux " << ux <<endl;
//       ExactS(X_P, Y_P, Exact_Surf);

      Mult = LineWeights[k]*normn*r_axial;
      h_K +=normn;
//       errors[0] +=Mult*(Exact_Surf[0]-U)*(Exact_Surf[0]-U);
      errors[0] +=U*Mult;
      errors[1] +=Mult;
     } //  for(k=0;k<N_LinePoints;k++)
  } // for(i=0;i<N
//     OutPut("h_K_min and h_K_max of free surface: "<< h_K_min << " " << h_K_max<<endl; );

   errors[0] *=2.*Pi;
   errors[1] *=2.*Pi;

//    OutPut( "Time, Surfactant Mass " <<TDatabase::TimeDB->CURRENTTIME<< " " <<errors[0]<< " "<<endl);
//    OutPut( "Time, Surface area " <<TDatabase::TimeDB->CURRENTTIME<< " " <<errors[1]<< " "<<endl);

// exit(0);
}


void Get_FeFunction2DMass(TFEFunction2D *fefunction, double *parameters)
 {
  int i,j,k,l, polydegree, N_QFPoints, ORDER;
  int N_Cells, N_Joints, N_Vertices;
  int *BeginIndex, *GlobalNumbers, *DOF, N_BF;

  double *U, Mult, r_axial, val, mass, volume, Concentration;
  double *weights, *xi, *eta;
  double values[MaxN_QuadPoints_2D][MaxN_BaseFunctions2D];
  double AbsDetjk[MaxN_QuadPoints_2D], X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];

  TJoint *joint;
  TFESpace2D *FeSpace;
  TBaseCell *cell;
  TCollection *coll;
  JointType jointtype;
  BoundTypes bdtype;
  RefTrans2D RefTrans;
  boolean IsIsoparametric;
  QuadFormula2D QuadFormula;
  TQuadFormula2D *qf2;
  FE2D FEid;
  TBaseFunct2D *bf;
  TRefTrans2D *F_K;

  FeSpace = fefunction->GetFESpace2D();
  BeginIndex = FeSpace->GetBeginIndex();
  GlobalNumbers = FeSpace->GetGlobalNumbers();
  U = fefunction->GetValues();

  coll = FeSpace->GetCollection();
  N_Cells = coll->GetN_Cells();

  mass = 0.;
  volume = 0.;
  Concentration = 0.;
  for(i=0;i<N_Cells;i++)
   {
    cell = coll->GetCell(i);
    FEid = FeSpace->GetFE2D(i, cell);

    RefTrans = TFEDatabase2D::GetRefTrans2D_IDFromFE2D(FEid);
    N_Joints = cell->GetN_Joints();
    IsIsoparametric = FALSE;
    if(TDatabase::ParamDB->USE_ISOPARAMETRIC)
     {
     for(j=0;j<N_Joints;j++)
      {
       joint = cell->GetJoint(j);
       jointtype = joint->GetType();
       if(jointtype == BoundaryEdge)
        {
         bdtype = ((TBoundEdge *)(joint))->GetBoundComp()->GetType();
         if(bdtype != Line)  IsIsoparametric = TRUE;
        }
       if(jointtype == InterfaceJoint)
        {
          bdtype = ((TInterfaceJoint *)(joint))->GetBoundComp()->GetType();
          if(bdtype != Line)
            IsIsoparametric = TRUE;
        }
        if(jointtype == IsoInterfaceJoint || jointtype == IsoBoundEdge)
         IsIsoparametric = TRUE;

      } // for(j=0;j<
     } // if(TDatabase::ParamDB->USE_ISOPARAMETRIC)

   if(IsIsoparametric)
    {
      switch(N_Joints)
      {
        case 4:
          RefTrans = QuadIsoparametric;
        break;

        case 3:
          RefTrans = TriaIsoparametric;
        break;
      }
    } // endif IsIsoparametric

    F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
    ORDER = TFEDatabase2D::GetAccuracyFromFE2D(FEid);
    switch(RefTrans)
    {
      case TriaAffin:
        polydegree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEid);
        QuadFormula = TFEDatabase2D::GetQFTriaFromDegree(9);
	qf2 = TFEDatabase2D::GetQuadFormula2D(QuadFormula);
        qf2->GetFormulaData(N_QFPoints, weights, xi, eta);
        ((TTriaAffin *)F_K)->SetCell(cell);
//         locvol = ((TTriaAffin *)rt)->GetVolume();
        ((TTriaAffin *)F_K)->GetOrigFromRef(N_QFPoints, xi, eta,
                                         X, Y, AbsDetjk);
      break;

      case TriaIsoparametric:
        polydegree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEid);
        QuadFormula = TFEDatabase2D::GetQFTriaFromDegree(9);
	qf2 = TFEDatabase2D::GetQuadFormula2D(QuadFormula);
        qf2->GetFormulaData(N_QFPoints, weights, xi, eta);
        ((TTriaIsoparametric *)F_K)->SetApproximationOrder(ORDER);
        ((TTriaIsoparametric *)F_K)->SetQuadFormula(QuadFormula);
        ((TTriaIsoparametric *)F_K)->SetCell(cell);
//         locvol = ((TTriaIsoparametric *)F_K)->GetVolume();
        ((TTriaIsoparametric *)F_K)->GetOrigFromRef(N_QFPoints, xi, eta,
                                         X, Y, AbsDetjk);
      break;

      case QuadAffin:
        polydegree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEid);
        QuadFormula = TFEDatabase2D::GetQFQuadFromDegree(3*polydegree);
	qf2 = TFEDatabase2D::GetQuadFormula2D(QuadFormula);
        qf2->GetFormulaData(N_QFPoints, weights, xi, eta);
        ((TQuadAffin *)F_K)->SetCell(cell);
//         locvol = ((TQuadAffin *)rt)->GetVolume();
        ((TQuadAffin *)F_K)->GetOrigFromRef(N_QFPoints, xi, eta,
                                         X, Y, AbsDetjk);
      break;

      case QuadBilinear:
        polydegree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEid);
        QuadFormula = TFEDatabase2D::GetQFQuadFromDegree(3*polydegree);
	qf2 = TFEDatabase2D::GetQuadFormula2D(QuadFormula);
        qf2->GetFormulaData(N_QFPoints, weights, xi, eta);
        ((TTriaIsoparametric *)F_K)->SetApproximationOrder(polydegree);
        ((TQuadBilinear *)F_K)->SetCell(cell);
        ((TQuadBilinear *)F_K)->GetOrigFromRef(N_QFPoints, xi, eta,
                                         X, Y, AbsDetjk);
      break;

      case QuadIsoparametric:
        polydegree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEid);
        QuadFormula = TFEDatabase2D::GetQFQuadFromDegree(3*polydegree);
	qf2 = TFEDatabase2D::GetQuadFormula2D(QuadFormula);
        qf2->GetFormulaData(N_QFPoints, weights, xi, eta);
        ((TQuadIsoparametric *)F_K)->SetApproximationOrder(polydegree);
        ((TQuadIsoparametric *)F_K)->SetQuadFormula(QuadFormula);
        ((TQuadIsoparametric *)F_K)->SetCell(cell);
//         locvol = ((TQuadIsoparametric *)rt)->GetVolume();
        ((TQuadIsoparametric *)F_K)->GetOrigFromRef(N_QFPoints, xi, eta,
                                         X, Y, AbsDetjk);
      break;
    }

    // find basis functions on cell i
    bf = TFEDatabase2D::GetBaseFunct2DFromFE2D(FEid);
    N_BF = bf->GetDimension();
    DOF = GlobalNumbers + BeginIndex[i];

    for(k=0;k<N_QFPoints;k++)
     {
      bf->GetDerivatives(D00, xi[k], eta[k], values[k]);
      r_axial = fabs(X[k]);
      if(X[k]<=0)
       {
        cout <<"X[k] negative in Get_KE change Quad rule " <<  X[k] <<endl;
//         exit(0);
       }

      Mult = r_axial*weights[k]*AbsDetjk[k];
      val = 0.;
      for(l=0;l<N_BF;l++)
       {
        j = DOF[l];
        val += U[j]*values[k][l];
       }

     mass += val*Mult;
     volume += Mult;
    } //  for(k=0;k<N_QFPoints;
   } //  for(i=0;i<N_Cells;i++)

   parameters[0] = 2.*Pi*mass;
   parameters[1] = 2.*Pi*volume;
   Concentration = mass/volume;
   parameters[2] = Concentration;


//    OutPut( "Time, C Mass " <<TDatabase::TimeDB->CURRENTTIME<< " " <<parameters[0]<< " "<<endl);
//    OutPut( "Time, C Surface area " <<TDatabase::TimeDB->CURRENTTIME<< " " <<parameters[1]<< " "<<endl);

 }



// ======================================================================
// utilities for main program
// ======================================================================

int main(int argc, char* argv[])
{
  TDomain *Domain = new TDomain();
  TDatabase *Database = new TDatabase();
  TFEDatabase2D *FEDatabase = new TFEDatabase2D();
  TCollection *coll, *NSE_coll, *mortarcoll = NULL;
  TBaseCell *Me, *cell, **Free_Cells, **NSE_Cells;
  TGridCell **DelCell;
  TFESpace2D *velocity_space, *pressure_space, *streamfunction_space, *convolution_space, *fesps;
  TFESpace2D  *velocity_space_output, *pressure_space_output;
  TFESpace2D *Grid_space, *vorticity_space, *surfact_space,*grid_space;

  TOutput2D *Output;
  TFEVectFunct2D *RefGridPos, *AuxGridPos, *GridPos, *ReparamPos;
  TFEFunction2D *fefct[4];
  TFESpace2D *fesp[3], *ferhs_T[3], *ferhs[2];
  TAuxParam2D *aux;
  TDiscreteForm2D *DiscreteFormMatrixT_MRhs;
  TDiscreteForm2D *DiscreteForm;
  TMatrix2D *MATRICES[4];
  TSquareMatrix2D *SQMATRICES[8], *SQMATRICES_GRID[4], *SQMATRICES_SURFACT[4];
  TDiscreteForm2D *DiscreteFormGalerkin;
  TDiscreteForm2D *DiscreteFormNLGalerkin;
  TDiscreteForm2D *DiscreteFormGrid, *DiscreteFormSurfact, *DiscreteFormSurfact_SUPG;
  TSquareMatrix **sqmatrices = (TSquareMatrix **)SQMATRICES;
  TMatrix **matrices = (TMatrix **)MATRICES;
  MatVecProc *MatVect;
  DefectProc *Defect;

  TFESpace2D **FESpaces_All = new TFESpace2D *[7];      
  TFEFunction2D **FEFunctions_All = new TFEFunction2D *[7];    
  TFEVectFunct2D **FEVectFuncts_All = new TFEVectFunct2D*[2];
  TStructure2D **Structure_All = new TStructure2D *[2];
  TSquareStructure2D **SquareStructure_All = new TSquareStructure2D *[3];
  TSquareMatrix2D **SqMat_All = new TSquareMatrix2D *[16];
  TMatrix2D **Mat_All = new TMatrix2D *[4];
  TSquareMatrix1D **SqMat_IFace = new TSquareMatrix1D*[2];
  
  TSquareStructure1D **IFaceStruct = new TSquareStructure1D*[1];
  TFEFunction1D **IFaceFeFunct = new TFEFunction1D*[1];
  
  
#ifdef __SURFACT__ 
  // variables for surfactant
  TDomain *IFaceDomain = new TDomain();  
  
  int **N_List = new int*[4], N_IFaceCells, N_FESpaces_low, N_IActive;   
   
  FE1D *FE1D_List;
  TCollection *IFace_Coll;
  TFESpace1D *IFaceSurfact_space;
  TFESpace1D **IFaceFeSpaces = new TFESpace1D*[2];
  TSquareMatrix1D *SQMATRICES_IFace[2];
  TFESpace1D *IFacefesp[1], *IFaceferhs[1];  
#endif
  
  double total_time,*Coordinates, IntitialTimestep;
  double  t, teta, dt,x,y,gamma, tx,ty,sx,sy, R_Theta[3];
  double left, right, top, bottom,T_a, T_b;
  double x0, y0,x1,y1,hi, residual, impuls_residual, oldresidual, solver_time;
  double end_time, t1, t2, t4, t3;
  double *B, *defect;
  double *RHSs[3], *refpos, *auxpos, *pos, *ReparamDisp, *ReparamMeshVelo;
  double  TX[2], TY[2], solver_time_curr;
  double SLPX, SLPY, *Entries[4], tau, oldtau, limit, *sol_output, InitVolume, CurrVolume;  
  double Lx, Ly, Rx, Ry,  x2, y2, x3, y3, x4, y4, fh, fhlimit, fhtot, fhmin, fhmax;
  double *Angle = new double[2], **FreePts = new double *[2];  
  double **Sol_All = new double *[5], *tmp_Gsol, *tmp_Gd;
  double **Rhs_All = new double *[4], MaxWetD=0., T_MaxWetD=0.;
  double *S_BX, *S_BY, *Intpol_Coord, *Intpol_VeloValues, h_interface;
  double *C_defect, *Csol_old, *oldsol, *C_B, *CRhs_old;
  double *I_defect, *Isol_old, *I_B, *IRhs_old, *Csol_nonlinearstep;
  double  *CRHSs[1], GammaXmaxVal[2];
  double residual_scalar, oldresidual_scalar, limit_scalar, *SRHSs[1];
  double Surf_Mass[2], Params[10],Initial_SurfactMass, Initial_IFaceSurfactMass;
  double tn_2=0., tn_1=0., tn=0.,  ztn_2=0., ztn_1=0., ztn=0.;
  
  int ret,N_RootCells,N_Cells,N_Joints, N_Vertices,N_G, N_NSE_Cells, N_NonNSE_Cells, N_E;
  int N_SlipBound_Vert,  N_FreeBound_Vert,  N_AxialBound_Vert,N_Interf_Vertices;
  int In_Index,CurrComp,CurrVertex, img=1, RemeshImg=1, N_BData=0;
  int ID,CurrNeib,Neib[2], N_SquareMatrices, N_RectMatrices;
  int a,b,i,X,Y,j,k,l,len1, len2,Neighb_tmp,Indextemp;
  int *PartMarker, *Triangles, *NSE_GlobalCllNo;
  int *PointNeighb,maxEpV = 0, Max_It, N_U_output, N_P_output;
  int  N_U, N_P,N_Unknowns,N_pressureDOF,N_Rhs,N_FESpaces;
  int  N_Hori1, N_Hori2,N_Verti,N_Boundary_Vert,N_surfactDOF, N_surfactActive, N_surfactNonActive,N_IsurfactDOF;
  int velocity_space_code, pressure_space_code;
  int m,m0,m1,m2,m3,m4,m5,m6, N_Active, N_LinIterCurr, N_LinIter;
  int ORDER,VSP, N_GActive, N_GBoundaryNodes, N_SubSteps, very_first_time=0;
  int **IsoCellEdgeNos, *GridKCol, *GridRowPtr, *RowPtr, last_sq;
  int  *JointDOF, N_DOF, dof, *DOF, *GlobalNumbers, *BeginIndex, N_ReParam=0, N_Remesh=0, *N_MovVert;
  int surf_couple_var, Max_It_scalar, N_S, IncrTimeStepCount=0, DecrTimeStepCount=0;
   
  char *PRM, *GEO, *PsBaseName, *VtkBaseName;
  char ReadinDat[] = "readin.dat";
  char CString[] = "C";
  char NameString[]  = "name";
  char UString[] = "u";
  char PString[] = "p";
  char WString[] = "w";
  const char vtkdir[] = "VTK";
  const char BDdir[] = "BDData";
  char IFaceSString[] = "C_I";
  
  bool remeshed=FALSE, reparam = FALSE, SurfSurfactReparam=FALSE, ReduceTiemStep=FALSE;

  TBoundPart *BoundPart;
  TBoundComp *BoundComp;
  TBdLine *UpdateAxialBound, *UpdateHoriWallBound;
  TBdCircle *UpdateFreeBound;
  TBaseCell **CellTree;
  TVertex *vert, **VertexDel,**NewVertices;
  TJoint *Joint;
  TBoundEdge *Solid_Joint;
  TBoundEdge ***Bound_Joint = new TBoundEdge**[2];
  TIsoBoundEdge **Free_Joint, *IsoJoint;
  TVertex ***MovBoundVert = new TVertex**[3];
  TVertex *temp_Mov;
  TBoundEdge *tempSlip_Joint;
  FE2D FeId;
  TFEDesc2D *FeDesc;

  BoundCondFunct2D *BoundaryConditions[2];
  BoundValueFunct2D *BoundValues[2];

  BoundCondFunct2D *GridBoundaryConditions[1];
  BoundValueFunct2D *GridBoundValues[1];

#ifdef __SURFACT__
  BoundCondFunct2D *SurfactBoundaryConditions[1];
  BoundValueFunct2D *SurfactBoundValues[1];
#endif

  struct triangulateio In, Out;

  std::ostringstream opts;
  std::ostringstream os;

  os << " ";
  opts << " ";
//======================================================================
// read parameter file
//======================================================================
  total_time = GetTime();
  if(argc>=2)
    ret=Domain->ReadParam(argv[1]);
  else
    ret=Domain->ReadParam(ReadinDat);

  if(ret==-1)
  {
    exit(-1);
  }

  OpenFiles();
  OutFile.setf(std::ios::scientific);

  Database->CheckParameterConsistencyNSE();
  Database->WriteParamDB(argv[0]);
  Database->WriteTimeDB();
  ExampleFile();
//======================================================================
// copy read parameters into local variables
//======================================================================
  if(TDatabase::ParamDB->DISCTYPE==2 )
  {
    OutPut("SDFEM does not work!" << endl);
    Error("SDFEM does not work!" << endl);
    exit(4711);
  }
  if( TDatabase::ParamDB->DISCTYPE==5 )
  {
    OutPut("DISCTYPE 5 NOT IMPLEMENTED!" << endl);
    Error("DISCTYPE 5 NOT IMPLEMENTED!" << endl);
    exit(4711);
  }

  PRM = TDatabase::ParamDB->BNDFILE;
  GEO = TDatabase::ParamDB->GEOFILE;
  
  PsBaseName = TDatabase::ParamDB->PSBASENAME;
  VtkBaseName = TDatabase::ParamDB->VTKBASENAME;

  mkdir(vtkdir, 0777);
  mkdir(BDdir, 0777);
#define __AXIAL3D__ 

//======================================================================
// read boundary parameterization and initialize coarse grid
//======================================================================

   Domain->Init(PRM, GEO);

// write grid into an Postscript file
   os.seekp(std::ios::beg);
   os << "Domain_old" << ".ps" << ends;
   Domain->PS(os.str().c_str(),It_Finest,0);

   N_FreeBound_Vert = int (TDatabase::ParamDB->P6);    //Freesurf except end point
   S_BX = new double[N_FreeBound_Vert];
   S_BY = new double[N_FreeBound_Vert];
   
   dt = TDatabase::ParamDB->P7/TDatabase::ParamDB->P6;
   h_interface = dt;
   
   
   for(i=0;i<N_FreeBound_Vert; i++)
    {
     S_BX[i] = ((double)i)*dt;
     S_BY[i] = 0.;       
     //cout << "x " << S_BX[i] << endl;
    }

   /** Generate new mesh */
   N_MovVert = new int[3];
//    TriaReMeshGen(Domain, N_FreeBound_Vert, S_BX, S_BY);
   TriaReMeshGen_hori(Domain, N_FreeBound_Vert, S_BX, S_BY);
   BoundPart = Domain->GetBdPart(0);   
#ifdef  __WITHCYLINDERTICKNESS__     
//    UpdateHoriWallBound = (TBdLine*)BoundPart->GetBdComp(1);      
//    UpdateAxialBoundTop = (TBdLine*)BoundPart->GetBdComp(4);   
   UpdateAxialBound = (TBdLine*)BoundPart->GetBdComp(4);
#else
   UpdateAxialBound = (TBdLine*)BoundPart->GetBdComp(3);
#endif   
   
   // write grid into an Postscript file
   os.seekp(std::ios::beg);
   os << "Domain" << ".ps" << ends;
   Domain->PS(os.str().c_str(),It_Finest,0);
   
//    exit(0);
//======================================================================
// Triangular for grid generation end
//======================================================================   

   //Initialize DiscreteForms         
  InitializeDiscreteForms_Moving(DiscreteFormGalerkin, DiscreteFormNLGalerkin,
                                 DiscreteFormGrid, LinCoeffs, GridCoeffs);
 
  InitializeDiscreteForms_Moving(DiscreteFormSurfact, DiscreteFormSurfact_SUPG, SurfactCoeffs);   
   
  BoundaryConditions[0] = BoundCondition;
  BoundaryConditions[1] = BoundCondition;

  BoundValues[0] = U1BoundValue;
  BoundValues[1] = U2BoundValue;

  GridBoundaryConditions[0] = GridBoundCondition;
  GridBoundValues[0] = GridBoundValue;

#ifdef __SURFACT__
  SurfactBoundaryConditions[0] = SurfactBoundCondition;
  SurfactBoundValues[0] = SurfactBoundValue;
#endif


//======================================================================
// construct all finite element spaces
//======================================================================
  coll=Domain->GetCollection(It_Finest, 0);
  N_Cells = coll->GetN_Cells();
  
#ifdef __SURFACT__
//  list of outer phase cells containing interafce
  N_List[0] = new int[N_Cells];   // Cell_No
  N_List[1] = new int[N_Cells];   //Joint_No
  Domain2DSurf_2Phase(coll, IFaceDomain, N_List);  
  
  IFace_Coll = IFaceDomain->GetCollection(It_Finest, 0);
  N_IFaceCells= IFace_Coll->GetN_Cells();
  
  OutPut("N_Cells      : "<< setw(10) << N_Cells  << endl);
  OutPut("N_IFaceCells : "<<  setw(10) <<  N_IFaceCells  << endl); 
    
  FE1D_List = new FE1D[N_Cells];
  for(j=0;j<N_Cells;j++)
  FE1D_List[j] = FE1D(TDatabase::ParamDB->ANSATZ_ORDER);

  TDatabase::IteratorDB[It_EQ]->SetParam(Domain);
  TDatabase::IteratorDB[It_LE]->SetParam(Domain);
  TDatabase::IteratorDB[It_Finest]->SetParam(Domain);
  TDatabase::IteratorDB[It_Between]->SetParam(Domain);
  TDatabase::IteratorDB[It_OCAF]->SetParam(Domain);    
# endif
  
  ORDER = 0;
  VSP = TDatabase::ParamDB->VELOCITY_SPACE;

  if(abs(VSP) > 20)
   {ORDER = abs(VSP) - 20;}
  else if (abs(VSP) > 10)
    {ORDER = abs(VSP) - 10;}
  else ORDER = abs(VSP);

//   NSE_coll = new TCollection(N_NSE_Cells, NSE_Cells);
//======================================================================
// construct all finite element spaces
//======================================================================
  // get velocity and pressure spacess
  GetVelocityAndPressureSpace(coll,BoundCondition,
                              mortarcoll, velocity_space,
                              pressure_space, &pressure_space_code,
                              TDatabase::ParamDB->VELOCITY_SPACE,
                              TDatabase::ParamDB->PRESSURE_SPACE);
  velocity_space_code =   TDatabase::ParamDB->VELOCITY_SPACE;
  FESpaces_All[0] = velocity_space;
  FESpaces_All[1] =  pressure_space;  
  N_Active =  FESpaces_All[0]->GetActiveBound();
  N_U = FESpaces_All[0]->GetN_DegreesOfFreedom();
  N_P = FESpaces_All[1]->GetN_DegreesOfFreedom();

  GlobalNumbers = FESpaces_All[0]->GetGlobalNumbers();
  BeginIndex = FESpaces_All[0]->GetBeginIndex();
  
// mesh velocity space 
   grid_space = new TFESpace2D(coll, NameString, CString, GridBoundCondition, 1, NULL);
   FESpaces_All[2] =  grid_space;     
   N_G = FESpaces_All[2]->GetN_DegreesOfFreedom();
   N_GActive = FESpaces_All[2]->GetActiveBound();
   N_GBoundaryNodes = N_G - N_GActive;
   
   OutPut("N_G          : "<< setw(10) << N_G  << endl);
   OutPut("N_GActive    : "<< setw(10) << N_GActive  << endl);
   
// surfact space
#ifdef __SURFACT__
   surfact_space = new TFESpace2D(coll, NameString, CString, SurfactBoundCondition,
                                  TDatabase::ParamDB->ANSATZ_ORDER, NULL);
   FESpaces_All[3] =  surfact_space;  
   N_surfactDOF = FESpaces_All[3]->GetN_DegreesOfFreedom();
   N_surfactActive = FESpaces_All[3]->GetActiveBound();
   N_surfactNonActive = N_surfactDOF - N_surfactActive;
   OutPut("N_surfactDOF    : "<< setw(10) << N_surfactDOF  << endl);
 
   IFaceSurfact_space = new TFESpace1D(IFace_Coll , IFaceSString, IFaceSString, FE1D_List);

   IFaceFeSpaces[0] = IFaceSurfact_space;
   N_IsurfactDOF = IFaceFeSpaces[0]->GetN_DegreesOfFreedom();
   N_IActive = IFaceFeSpaces[0]->GetActiveBound();   
   
   FESpaces_All[4] =  new TFESpace2D(coll, NameString, IFaceSString, SurfactBoundCondition,
                                   ContP_USpace, TDatabase::ParamDB->ANSATZ_ORDER, NULL);

   N_S =  FESpaces_All[4]->GetN_DegreesOfFreedom();
#endif   
   
//======================================================================
// construct all finite element functions
//======================================================================
  N_Unknowns = 2*N_U + N_P;
  Sol_All[0] = new double[N_Unknowns];  
  Rhs_All[0] = new double[N_Unknowns];
  
  OutPut("dof velocity : "<< setw(10) << 2* N_U << endl);
  OutPut("dof pressure : "<< setw(10) << N_P << endl);
  OutPut("dof all      : "<<  setw(10) << N_Unknowns  << endl);  
  
  B = new double[N_Unknowns];
  defect = new double[N_Unknowns];
  
  memset(Sol_All[0], 0, N_Unknowns*SizeOfDouble);
  memset(Rhs_All[0], 0, N_Unknowns*SizeOfDouble);

  //velo vect
  FEVectFuncts_All[0] =  new TFEVectFunct2D(FESpaces_All[0], UString, UString, Sol_All[0], N_U, 2);
  FEFunctions_All[0] = FEVectFuncts_All[0]->GetComponent(0);
  FEFunctions_All[1] = FEVectFuncts_All[0]->GetComponent(1); 
  FEFunctions_All[0]->Interpolate(InitialU1);
  FEFunctions_All[1]->Interpolate(InitialU2);
  
  //pressure FeFunction
  FEFunctions_All[2] = new TFEFunction2D(FESpaces_All[1], PString,  PString,  Sol_All[0]+2*N_U, N_P);

  tmp_Gd = new double[2*N_G];
  tmp_Gsol = new double[2*N_G];
  Sol_All[1] = new double[2*N_G];
  Rhs_All[1] = new double[2*N_G];   
  
  memset(Sol_All[1], 0, 2*N_G*SizeOfDouble);   
  FEVectFuncts_All[1]  = new TFEVectFunct2D(FESpaces_All[2], WString, WString, Sol_All[1], N_G, 2);
  
  FEFunctions_All[3] = FEVectFuncts_All[1]->GetComponent(0);
  FEFunctions_All[4] = FEVectFuncts_All[1]->GetComponent(1);
//======================================================================
// surfact space finite element functions
//======================================================================
#ifdef __SURFACT__
  Sol_All[2] = new double[N_surfactDOF];
  Rhs_All[2] = new double[N_surfactDOF];
  Csol_old = new double[N_surfactDOF];
  CRhs_old = new double[N_surfactDOF];
  C_B = new double[N_surfactDOF];  
  C_defect = new double[N_surfactDOF];  
  Csol_nonlinearstep = new double[N_surfactDOF];
   
  memset(Sol_All[2], 0, N_surfactDOF*SizeOfDouble);
  memset(Csol_old, 0, N_surfactDOF*SizeOfDouble);
  memset(Rhs_All[2], 0, N_surfactDOF*SizeOfDouble);

  // surfact fefunction
  FEFunctions_All[5] = new TFEFunction2D(FESpaces_All[3], CString, CString, Sol_All[2], N_surfactDOF);
  FEFunctions_All[5]->Interpolate(InitialSuract);
 
 
   Sol_All[3] = new double[N_IsurfactDOF];
   Rhs_All[3] = new double[N_IsurfactDOF];

   I_defect = new double[N_IsurfactDOF];
   Isol_old = new double[N_IsurfactDOF];
   I_B = new double[N_IsurfactDOF];
   IRhs_old = new double[N_IsurfactDOF];

   
   memset(Sol_All[3], 0, N_IsurfactDOF*SizeOfDouble);
   memset(Rhs_All[3], 0, N_IsurfactDOF*SizeOfDouble);

   memset(I_defect, 0, N_IsurfactDOF*SizeOfDouble);
   memset(Isol_old, 0, N_IsurfactDOF*SizeOfDouble);
   memset(I_B, 0, N_IsurfactDOF*SizeOfDouble);
 
   
   IFaceFeFunct[0] = new TFEFunction1D(IFaceFeSpaces[0], IFaceSString, IFaceSString, Sol_All[3], N_IsurfactDOF);
   IFaceFeFunct[0]->Interpolate(InitialS);   
    
   //for output 
   Sol_All[4] =  new double[N_S];
   memset(Sol_All[4], 0, N_S*SizeOfDouble);
   FEFunctions_All[6]  = new TFEFunction2D(FESpaces_All[4], IFaceSString,  IFaceSString, Sol_All[4],  N_S);   
 #endif 
  
//======================================================================
// allocate memory for all matrices
//======================================================================
  Structure_All[0] = new TStructure2D(FESpaces_All[1], FESpaces_All[0]);  // B
  Structure_All[1] = new TStructure2D(FESpaces_All[0], FESpaces_All[1]); // BT
  
  //velo
  SquareStructure_All[0] = new TSquareStructure2D(FESpaces_All[0]);  
  SquareStructure_All[0]->Sort();

  // grid 
  SquareStructure_All[1] = new TSquareStructure2D(FESpaces_All[2]); 
  SquareStructure_All[1]->Sort();
  
  //surfact
#ifdef __SURFACT__
  SquareStructure_All[2] = new TSquareStructure2D(FESpaces_All[3]);
  SquareStructure_All[2]->Sort();
#endif
  
#ifdef __SURFACT__ 
  /* interface surfactant matrices */
  IFaceStruct[0] = new TSquareStructure1D(IFaceFeSpaces[0]);
  IFaceStruct[0]->Sort();
#endif   
    
    
  // for NSE
  MatVect = MatVect_NSE4;
  Defect = Defect_NSE4;

  SqMat_All[0] = new TSquareMatrix2D(SquareStructure_All[0]); // M11
  SqMat_All[1] = new TSquareMatrix2D(SquareStructure_All[0]); // M12
  SqMat_All[2] = new TSquareMatrix2D(SquareStructure_All[0]); // M21
  SqMat_All[3] = new TSquareMatrix2D(SquareStructure_All[0]); // M22
  
  SqMat_All[4] = new TSquareMatrix2D(SquareStructure_All[0]); // A11
  SqMat_All[5] = new TSquareMatrix2D(SquareStructure_All[0]); // A12
  SqMat_All[6] = new TSquareMatrix2D(SquareStructure_All[0]); // A21
  SqMat_All[7] = new TSquareMatrix2D(SquareStructure_All[0]); // A22
  SqMat_All[8] = new TSquareMatrix2D(SquareStructure_All[0]); // F11
  SqMat_All[9] = new TSquareMatrix2D(SquareStructure_All[0]); // F22

  Mat_All[0] = new TMatrix2D(Structure_All[0]); // B1
  Mat_All[1] = new TMatrix2D(Structure_All[0]); // B2
  Mat_All[2] = new TMatrix2D(Structure_All[1]); // B1T
  Mat_All[3] = new TMatrix2D(Structure_All[1]); // B2T

  // for mesh
  SqMat_All[10] = new TSquareMatrix2D(SquareStructure_All[1]); // G11
  SqMat_All[11] = new TSquareMatrix2D(SquareStructure_All[1]); // G12
  SqMat_All[12] = new TSquareMatrix2D(SquareStructure_All[1]); // G21
  SqMat_All[13] = new TSquareMatrix2D(SquareStructure_All[1]); // G22

  // for heat
#ifdef __SURFACT__  
  SqMat_All[14]  = new TSquareMatrix2D(SquareStructure_All[2]); // T_M
  SqMat_All[15] = new TSquareMatrix2D(SquareStructure_All[2]); // T_A
#endif     
  IsoCellEdgeNos = new int *[2];
  
#ifdef __SURFACT__   
 /* interface surfactant matrices */
  SqMat_IFace[0] = new TSquareMatrix1D(IFaceStruct[0]); // Gamma_A
  SqMat_IFace[1] = new TSquareMatrix1D(IFaceStruct[0]); // Gamma_M  
#endif   
  
 //  ====================================================================================   

  GetMovingBoundData(coll, N_MovVert, Bound_Joint, MovBoundVert, Free_Joint,
                     Free_Cells, IsoCellEdgeNos, S_BX[0], S_BY[0]);

  delete [] S_BX;
  delete [] S_BY;  
                                        
//  ====================================================================================  
// assemble matrix for grid moving - begin
//  ====================================================================================  
    fesp[0] = FESpaces_All[2];
    SQMATRICES_GRID[0] = SqMat_All[10];
    SQMATRICES_GRID[0]->Reset();
    SQMATRICES_GRID[1] = SqMat_All[11];
    SQMATRICES_GRID[1]->Reset();
    SQMATRICES_GRID[2] = SqMat_All[12];
    SQMATRICES_GRID[2]->Reset();
    SQMATRICES_GRID[3] = SqMat_All[13];
    SQMATRICES_GRID[3]->Reset();
    aux = new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);  
       
    Assemble2D(1, fesp,
             4, SQMATRICES_GRID,
             0, NULL,
             0, NULL, NULL,
             DiscreteFormGrid,
             GridBoundaryConditions,
             GridBoundValues,
             aux);
    delete aux;   
   
     Entries[0] = SqMat_All[10]->GetEntries();
     Entries[1] = SqMat_All[11]->GetEntries();
     Entries[2] = SqMat_All[12]->GetEntries();
     Entries[3] = SqMat_All[13]->GetEntries();

     GridKCol = SquareStructure_All[1]->GetKCol();
     GridRowPtr = SquareStructure_All[1]->GetRowPtr();

  // for Dirichlet rows in off-diagonal matrices
  memset(Entries[1] + GridRowPtr[N_GActive], 0, (GridRowPtr[N_G] - GridRowPtr[N_GActive])*SizeOfDouble);
  memset(Entries[2] + GridRowPtr[N_GActive], 0, (GridRowPtr[N_G] - GridRowPtr[N_GActive])*SizeOfDouble);

  refpos = new double[2*N_G];
  auxpos = new double[2*N_G];
  pos = new double[2*N_G];
  ReparamMeshVelo = new double[2*N_G];  
  ReparamDisp = new double[2*N_G]; 
  ReparamPos = new TFEVectFunct2D(FESpaces_All[2], WString, WString, ReparamDisp, N_G, 2);  
  RefGridPos = new TFEVectFunct2D(FESpaces_All[2], WString, WString, refpos, N_G, 2);
  AuxGridPos = new TFEVectFunct2D(FESpaces_All[2], WString, WString, auxpos, N_G, 2);
  GridPos = new TFEVectFunct2D(FESpaces_All[2], WString, WString, pos, N_G, 2);

  
  RefGridPos->GridToData();
  AuxGridPos->GridToData();
  GridPos->GridToData();

  // prepare output (maxn_fespaces,  maxn_scalar,  maxn_vect, maxn_parameters, domain)
  Output = new TOutput2D(1, 3, 1, 2, Domain);
  Output->AddFEVectFunct(FEVectFuncts_All[0]);
  Output->AddFEFunction(FEFunctions_All[2]);   
#ifdef __SURFACT__      
  Output->AddFEFunction(FEFunctions_All[5]);
  
  Output->AddFEFunction(FEFunctions_All[6]);
#endif      
//   Output->AddFEVectFunct(FEVectFuncts_All[1]);  
      
  os.seekp(std::ios::beg);
  Output->AddParameter(TDatabase::TimeDB->CURRENTTIME,os.str().c_str());      
 
     
#ifdef __SURFACT__        
      // interface surfactant
      MapSurfToDomain(IFaceFeFunct[0], FEFunctions_All[6], N_List[0], N_List[1]);
     
      GetSurfactMass(FEFunctions_All[6], IFaceFeFunct[0], N_List[0], N_List[1], Surf_Mass);
 
      Get_FeFunction2DMass(FEFunctions_All[5], Params);          
      
      Initial_IFaceSurfactMass = TDatabase::ParamDB->REACTOR_P14*Surf_Mass[0];
      Initial_SurfactMass = Params[0];    
      
      OutPut( "Time, Surfactant_Mass, Da*InterfaceSurfactant_Mass, InterfaceSurfactant_Conc " <<TDatabase::TimeDB->CURRENTTIME<<
              " " <<Params[0]<< " "<<TDatabase::ParamDB->REACTOR_P14*Surf_Mass[0]<< " "<<Surf_Mass[0]/Surf_Mass[1]<< " "<<endl);
      OutPut( "Time, Surfactant_Mass_dif, InterfaceSurfactant_Mass_diff " <<TDatabase::TimeDB->CURRENTTIME<<
              " " <<Params[0] - Initial_SurfactMass<< " "<<TDatabase::ParamDB->REACTOR_P14*Surf_Mass[0]-Initial_IFaceSurfactMass<< " "<<endl);
      OutPut( "Time, Total Mass, Total Mass diff, RelativeTotal Mass diff " <<TDatabase::TimeDB->CURRENTTIME<<
          " " <<Params[0]+TDatabase::ParamDB->REACTOR_P14*Surf_Mass[0]<< " " << ((Params[0]+TDatabase::ParamDB->REACTOR_P14*Surf_Mass[0]) - 
            (Initial_IFaceSurfactMass + Initial_SurfactMass)) << " " << ((Params[0]+TDatabase::ParamDB->REACTOR_P14*Surf_Mass[0]) - 
            (Initial_IFaceSurfactMass + Initial_SurfactMass))
            /(Initial_IFaceSurfactMass + Initial_SurfactMass) <<endl);     
#endif       
   if(TDatabase::ParamDB->WRITE_VTK)
     {      
      os.seekp(std::ios::beg);
        if(img<10) os << "VTK/"<<VtkBaseName<<".0000"<<img<<".vtk" << ends;
         else if(img<100) os <<"VTK/"<< VtkBaseName<<".000"<<img<<".vtk" << ends;
         else if(img<1000) os <<"VTK/"<< VtkBaseName<<".00"<<img<<".vtk" << ends;
         else if(img<10000) os << "VTK/"<<VtkBaseName<<".0"<<img<<".vtk" << ends;
         else  os << "VTK/"<<VtkBaseName<<"."<<img<<".vtk" << ends;
      Output->WriteVtk(os.str().c_str());
      img++;
     }
       
   Get_KE(FEVectFuncts_All[0], Params);
   InitVolume = CurrVolume = Params[0];
//    MovBoundVert[0][0]->GetCoords(Lx, Ly);
//    MovBoundVert[2][0]->GetCoords(Rx, Ry);
//    OutPut(setw(20)<<"T, Wett Len d : " << TDatabase::TimeDB->CURRENTTIME
//                   <<"   "<< Rx-Lx<<endl);
//    OutPut(setw(20)<<"T, Volume, KE : " << TDatabase::TimeDB->CURRENTTIME
//                   <<"   "<< CurrVolume<<"   "<< Params[1] <<endl);
//    OutPut(setw(20)<<"T, Volume Diff : "<< TDatabase::TimeDB->CURRENTTIME
//                   <<"   "<< CurrVolume - InitVolume << endl);
   
   OutPut(setw(25)<<"T, KE, Volume, Diff, Rel. Diff : "<< TDatabase::TimeDB->CURRENTTIME<<"   "<< Params[1] <<"   "<< Params[0]
                  <<"   "<< Params[0] - InitVolume<<"   "<< (Params[0] - InitVolume)/InitVolume << endl);
   
   Getcellangle(FESpaces_All[2], Angle);   
   OutPut( "MinAngle : "<< Angle[0]<< "  MaxAngle : "<<Angle[1]<< endl);
//    exit(0);
   
   TDatabase::ParamDB->P10 = 0; // free surf reparam
   TDatabase::ParamDB->P5 = 0;  // move boundary with velo
  
 
// #ifndef __SOLUBLE__  
//   Max_It_scalar = 1;
// #endif  
  
     
//           os.seekp(std::ios::beg);
//         if(img<10) os << "VTK/"<<VtkBaseName<<".0000"<<img<<".vtk" << ends;
//          else if(img<100) os <<"VTK/"<< VtkBaseName<<".000"<<img<<".vtk" << ends;
//          else if(img<1000) os <<"VTK/"<< VtkBaseName<<".00"<<img<<".vtk" << ends;
//          else if(img<10000) os << "VTK/"<<VtkBaseName<<".0"<<img<<".vtk" << ends;
//          else  os << "VTK/"<<VtkBaseName<<"."<<img<<".vtk" << ends;
//       Output->WriteVtk(os.str().c_str());
//       img++;
//       
//       exit(0);
//      cout << "test main " <<endl;
#ifdef __SURFACT__    
     MapSurfToDomain(IFaceFeFunct[0], FEFunctions_All[6], N_List[0], N_List[1]);
#endif   
     
//      ReParam_axial3D_U(N_MovVert[0], Free_Cells,  IsoCellEdgeNos[1], IsoCellEdgeNos[0],
//                       FEVectFuncts_All[0], FEFunctions_All[5], FEFunctions_All[6], TRUE);  
//      
     
//      cout << "test main " <<endl;
//  exit(0);
 
  
//======================================================================
// start of time cycle
//======================================================================
  end_time = TDatabase::TimeDB->ENDTIME;
  N_SubSteps = GetN_SubSteps(); 
  TDatabase::TimeDB->CURRENTTIME = TDatabase::TimeDB->STARTTIME;
  m=0;
  limit = TDatabase::ParamDB->SC_NONLIN_RES_NORM_MIN_SADDLE;
  Max_It = TDatabase::ParamDB->SC_NONLIN_MAXIT_SADDLE; 
  
  surf_couple_var = (int)TDatabase::ParamDB->REACTOR_P7;

  limit_scalar = TDatabase::ParamDB->REACTOR_P8;
  Max_It_scalar = int(TDatabase::ParamDB->REACTOR_P9);

  if(surf_couple_var==1) Max_It_scalar = 1; // explicit coupling
  
  solver_time = 0.0;
  N_LinIter = 0; 
  t3 = GetTime();
  total_time = t3 - total_time;
  
  IntitialTimestep = TDatabase::TimeDB->TIMESTEPLENGTH;

  
#ifdef __SURFACT__ 
      PrintSurfSurfactant(N_MovVert[0], MovBoundVert[0], FEFunctions_All[6], N_BData); 
#else
      PrintSurfSurfactant(N_MovVert[0], MovBoundVert[0], NULL, N_BData); 
#endif   
  
  
  while(TDatabase::TimeDB->CURRENTTIME< end_time)
  {
    // time cycle
    m++;
    TDatabase::TimeDB->INTERNAL_STARTTIME = TDatabase::TimeDB->CURRENTTIME;

      for(l=0;l<N_SubSteps;l++)   // sub steps of fractional step theta
      {  
        SetTimeDiscParameters(1);

        if (m==1)
        {
          OutPut("Theta1: " << TDatabase::TimeDB->THETA1<< endl);
          OutPut("Theta2: " << TDatabase::TimeDB->THETA2<< endl);
          OutPut("Theta3: " << TDatabase::TimeDB->THETA3<< endl);
          OutPut("Theta4: " << TDatabase::TimeDB->THETA4<< endl);
        }

        tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
        TDatabase::TimeDB->CURRENTTIME += tau;

        // working rhs array for NSE
        memset(B, 0, N_Unknowns*SizeOfDouble);

        GridVelo_imping(Entries, tmp_Gsol, tmp_Gd, Rhs_All[1],
                        GridKCol, GridRowPtr,
                        GridPos, AuxGridPos,
                        FEVectFuncts_All[0], tau,
                        FEVectFuncts_All[1], MovBoundVert, N_MovVert,
                        Free_Cells, IsoCellEdgeNos,
                        reparam, RefGridPos);

// ============================================================================================
//  Assemble NSE
// ============================================================================================     
        DiscreteForm = DiscreteFormGalerkin;

        SQMATRICES[0] = SqMat_All[4];
        SQMATRICES[1] = SqMat_All[5];
        SQMATRICES[2] = SqMat_All[6];
        SQMATRICES[3] = SqMat_All[7];
        SQMATRICES[4] = SqMat_All[0];
        SQMATRICES[5] = SqMat_All[3];

        SQMATRICES[6] = SqMat_All[1];
        SQMATRICES[7] = SqMat_All[2];

        MATRICES[0] = Mat_All[0];
        MATRICES[1] = Mat_All[1];
        MATRICES[2] = Mat_All[2];
        MATRICES[3] = Mat_All[3];

        SQMATRICES[0]->Reset();
        SQMATRICES[1]->Reset();
        SQMATRICES[2]->Reset();
        SQMATRICES[3]->Reset();
        SQMATRICES[4]->Reset();
        SQMATRICES[5]->Reset();
        SQMATRICES[6]->Reset();
        SQMATRICES[7]->Reset();

        MATRICES[0]->Reset();
        MATRICES[1]->Reset();
        MATRICES[2]->Reset();
        MATRICES[3]->Reset();
	
        N_SquareMatrices = 6;
        N_RectMatrices = 4;

       // parameters which are the same for all NSTYPEs
        N_Rhs = 2;
        N_FESpaces = 3;

        fesp[0] = FESpaces_All[0];
        fesp[1] = FESpaces_All[1];
        fesp[2] = FESpaces_All[2];

        fefct[0] = FEFunctions_All[0];
        fefct[1] = FEFunctions_All[1];
        fefct[2] = FEFunctions_All[3];
        fefct[3] = FEFunctions_All[4];
 
        ferhs[0] = FESpaces_All[0];
        ferhs[1] = FESpaces_All[0];	

        RHSs[0] = Rhs_All[0];
        RHSs[1] = Rhs_All[0] + N_U;
        RHSs[2] = Rhs_All[0] + 2*N_U;

        memset(Rhs_All[0], 0, N_Unknowns*SizeOfDouble);

       // 4 parameters are needed for assembling (u1_old, u2_old)
        aux =  new TAuxParam2D(MovingTNSN_FESpaces_Axial3D, MovingTNSN_Fct_Axial3D,
                               MovingTNSN_ParamFct_Axial3D,
                               MovingTNSN_FEValues_Axial3D,
                               fesp, fefct,
                               MovingTNSFct_Axial3D,
                               MovingTNSFEFctIndex_Axial3D,
                               MovingTNSFEMultiIndex_Axial3D,
                               MovingTNSN_Params_Axial3D, MovingTNSBeginParam_Axial3D);                        
                        
       //======================================================================
       // assembling of matrices for each level
       // A_11 , (A_12), (A_21), (A_22), M_11, (M_22)
       //======================================================================
       Assemble2D(N_FESpaces, fesp,
                  N_SquareMatrices, SQMATRICES,
                  N_RectMatrices, MATRICES,
                  N_Rhs, RHSs, ferhs,
                  DiscreteForm,
                  BoundaryConditions,
                  BoundValues,
                  aux);

      delete aux;

     SqMat_All[8]->Reset(); // Matrix entries for freesurf int;
     SqMat_All[9]->Reset(); // no need to calculate in nonlinear steps
     
#ifdef __SURFACT__ 
     MapSurfToDomain(IFaceFeFunct[0], FEFunctions_All[6], N_List[0], N_List[1]);      
#endif     
//      PrintSurfSurfactant(N_MovVert[0], MovBoundVert[0], FEFunctions_All[6], N_BData); 
       
     FreeSurf_axial3D_new(SqMat_All[8], SqMat_All[9],  RHSs[0], RHSs[1], BoundCondition, tau,
                          FEFunctions_All[0]->GetValues(), Params
#ifdef __SURFACT__  
                          , FEFunctions_All[6] 
#endif                               
    );
 	
     // Adding freesurf entries to A11 and A22
     MatAdd(SqMat_All[4], SqMat_All[8], 1);
     MatAdd(SqMat_All[7], SqMat_All[9], 1);

     // set rows of Dirichlet dof in off diagonal matrix blocks to zero    
     // get row in off diagonal matrix where the Dirichlet nodes start
     RowPtr = SqMat_All[5]->GetRowPtr();
     // compute number of entries starting from this row to the end of the matrix
     j = RowPtr[N_Active];
     k = RowPtr[N_U]-j;
     // get number of active dof
     // set these entries to zero
     memset(SqMat_All[5]->GetEntries()+j, 0, SizeOfDouble*k);
     memset(SqMat_All[6]->GetEntries()+j, 0, SizeOfDouble*k);

     // slip type bc detected, manipulation of matrices is necessary
     // this is done only at the very beginning
     // the matrices A_12, A_12, M_11, M_12, M_21, M_22, B1T, B2T
     //     stay unchanged during the complete solution process
     // the matrices A_11, A_22 are manipulated after their new
     //     assembling during the nonlinear iteration

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
          N_SquareMatrices = 8;
          N_RectMatrices = 2;
          N_Rhs = 2;
          DiscreteForm = NULL;

          SQMATRICES[0] = SqMat_All[4];
          SQMATRICES[1] = SqMat_All[7];
          SQMATRICES[2] = SqMat_All[5];
          SQMATRICES[3] = SqMat_All[6];
          SQMATRICES[4] = SqMat_All[0];
          SQMATRICES[5] = SqMat_All[3];
          SQMATRICES[6] = SqMat_All[1];
          SQMATRICES[7] = SqMat_All[2];

          MATRICES[0] = Mat_All[2];
          MATRICES[1] = Mat_All[3];

          fesp[0] = FESpaces_All[0];
          ferhs[0] = FESpaces_All[0];
          ferhs[1] = FESpaces_All[0];

          RHSs[0] = Rhs_All[0];
          RHSs[1] = Rhs_All[0]+N_U;

          aux = new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);

          Assemble2DSlipBC(N_FESpaces, fesp,
                           N_SquareMatrices, SQMATRICES,
                           N_RectMatrices, MATRICES,
                           N_Rhs, RHSs, ferhs,
                           DiscreteForm,
                           BoundaryConditions,
                           BoundValues,
                           aux, FEFunctions_All[0], FEFunctions_All[1]);

       delete aux;

      } // if (TDatabase::ParamDB->INTERNA

     //    scale the pressure matrices
     Dscal(Mat_All[2]->GetN_Entries(), tau, Mat_All[2]->GetEntries());
     Dscal(Mat_All[3]->GetN_Entries(), tau, Mat_All[3]->GetEntries());
     Dscal(Mat_All[0]->GetN_Entries(), tau, Mat_All[0]->GetEntries());
     Dscal(Mat_All[1]->GetN_Entries(), tau, Mat_All[1]->GetEntries());
     
     // update rhs
     Daxpy(N_Active, tau, Rhs_All[0], B);
     Daxpy(N_Active, tau, Rhs_All[0]+N_U, B+N_U);

     // update rhs by Laplacian and convective term initialy by current time step
     // scaled by current sub time step length and theta2
     // currently : M := M + gamma A
     // M = M + (-gamma - tau*TDatabase::TimeDB->THETA2) A
     MatAdd(SqMat_All[0], SqMat_All[4], - tau*TDatabase::TimeDB->THETA2);
     MatAdd(SqMat_All[1], SqMat_All[5], - tau*TDatabase::TimeDB->THETA2);
     MatAdd(SqMat_All[2], SqMat_All[6], - tau*TDatabase::TimeDB->THETA2);
     MatAdd(SqMat_All[3], SqMat_All[7], - tau*TDatabase::TimeDB->THETA2);

     // set current factor of steady state matrix
     gamma = -tau*TDatabase::TimeDB->THETA2;		     

     // defect = M * Sol
     // B:= B + defect (rhs)     
     MatVectActive(SqMat_All[0], Sol_All[0], defect);
     Daxpy(N_Active, 1, defect, B);
     MatVectActive(SqMat_All[1], Sol_All[0]+N_U, defect);
     Daxpy(N_Active, 1, defect, B);
     MatVectActive(SqMat_All[2], Sol_All[0], defect+N_U);
     Daxpy(N_Active, 1, defect+N_U, B+N_U);
     MatVectActive(SqMat_All[3], Sol_All[0]+N_U, defect+N_U);
     Daxpy(N_Active, 1, defect+N_U, B+N_U);     
     
    // set Dirichlet values
    // RHSs[0] still available from assembling
    memcpy(B+N_Active, RHSs[0]+N_Active, (N_U-N_Active)*SizeOfDouble);
    memcpy(B+N_Active+N_U, RHSs[1]+N_Active,(N_U-N_Active)*SizeOfDouble);

    // copy Dirichlet values from rhs into Sol[0][mg_level-1]
    memcpy(Sol_All[0]+N_Active, RHSs[0]+N_Active, (N_U-N_Active)*SizeOfDouble);
    memcpy(Sol_All[0]+N_Active+N_U, RHSs[1]+N_Active, (N_U-N_Active)*SizeOfDouble);

    //=====================================================================
    // the stiffness matrix is stored on M11, (M12, M21, M22)
    // assembling of system matrix
    //========================================================================
    // M = M + (-gamma + tau*TDatabase::TimeDB->THETA1) A
     MatAdd(SqMat_All[0], SqMat_All[4], -gamma + tau*TDatabase::TimeDB->THETA1);
     MatAdd(SqMat_All[1], SqMat_All[5], -gamma + tau*TDatabase::TimeDB->THETA1);
     MatAdd(SqMat_All[2], SqMat_All[6], -gamma + tau*TDatabase::TimeDB->THETA1);
     MatAdd(SqMat_All[3], SqMat_All[7], -gamma + tau*TDatabase::TimeDB->THETA1);
          
     // set current factor of steady state matrix
     gamma = tau*TDatabase::TimeDB->THETA1;
     
     OutPut(endl << "CURRENT TIME: ");
     OutPut(TDatabase::TimeDB->CURRENTTIME << endl);

     //======================================================================
     // nonlinear loop
     //======================================================================   
     N_LinIterCurr = 0;
     solver_time_curr = 0;     
     
     for(j=0;j<Max_It;j++)
      {
       memset(defect, 0, N_Unknowns*SizeOfDouble);

       SQMATRICES[0] = SqMat_All[0];
       SQMATRICES[1] = SqMat_All[1];
       SQMATRICES[2] = SqMat_All[2];
       SQMATRICES[3] = SqMat_All[3];
       MATRICES[0] = Mat_All[0];
       MATRICES[1] = Mat_All[1];
       MATRICES[2] = Mat_All[2];
       MATRICES[3] = Mat_All[3];      
       
      // compute defect
      Defect(sqmatrices, matrices, Sol_All[0], B, defect);

      residual =  Ddot(N_Unknowns, defect, defect);
      impuls_residual = Ddot(2*N_U, defect, defect);
      OutPut("nonlinear step " << setw(3) << j);
      OutPut(setw(14) << impuls_residual);
      OutPut(setw(14) << Ddot(N_P,defect+2*N_U,defect+2*N_U));
      OutPut(setw(14) << sqrt(residual));
      
      if(j>0)
       {
        OutPut(setw(14) << sqrt(residual)/oldresidual << endl);
       }
      else
       {
        OutPut(endl);
       }
       
      oldresidual = sqrt(residual);

      if ((((sqrt(residual)<=limit)||(j==Max_It-1)))  && (j>=TDatabase::ParamDB->SC_MINIT))
       {
        if (j==Max_It-1)
        j++;
        OutPut("ITE : " << setw(3) << j);
        OutPut(" (" << setw(3) << N_LinIterCurr << "/");
        OutPut(setw(3) << N_LinIter << " LINITE)");
        OutPut("  TIME FOR SOLVER : " << solver_time_curr << "/" << solver_time << "s");
        OutPut("  RES : " <<  sqrt(residual) << endl);
        // count total running time
        t4 =  GetTime();
        total_time += t4 - t3;
        t3 = t4;
        OutPut(TDatabase::TimeDB->CURRENTTIME << " total elapsed running time "<< total_time << endl);
        break;
       }
       
       //======================================================================
       // solve linear system
       //======================================================================
        t1 = GetTime();
        DirectSolver(SQMATRICES[0], SQMATRICES[1], SQMATRICES[2], SQMATRICES[3],
                     MATRICES[2], MATRICES[3], MATRICES[0], MATRICES[1],
                     B, Sol_All[0]);
	
	       

        t2 = GetTime();
        solver_time_curr = t2-t1;
        solver_time += solver_time_curr;
 
       //======================================================================
       // end solve linear system
       //======================================================================
       // restore mass matrices by subtracting the A-matrices
       MatAdd(SqMat_All[0], SqMat_All[4], -gamma);
       MatAdd(SqMat_All[3], SqMat_All[7], -gamma);

       //======================================================================
       // assemble new matrix due to nonlinearity
       //======================================================================
        GridVelo_imping(Entries, tmp_Gsol, tmp_Gd, Rhs_All[1],
                        GridKCol, GridRowPtr,
                        GridPos, AuxGridPos,
                        FEVectFuncts_All[0], tau,
                        FEVectFuncts_All[1], MovBoundVert, N_MovVert,
                        Free_Cells, IsoCellEdgeNos,
                        reparam, RefGridPos);
 
       DiscreteForm = DiscreteFormNLGalerkin;	 
       N_RectMatrices = 0;
       N_Rhs = 0;
       N_FESpaces = 3;

       SQMATRICES[0] = SqMat_All[4];
       SQMATRICES[1] = SqMat_All[7];
       SQMATRICES[0]->Reset();
       SQMATRICES[1]->Reset();

       N_SquareMatrices = 2;
       last_sq = 1;
       
       fesp[0] = FESpaces_All[0];
       fesp[1] = FESpaces_All[1];
       fesp[2] = FESpaces_All[2];

       fefct[0] = FEFunctions_All[0];
       fefct[1] = FEFunctions_All[1];
       fefct[2] = FEFunctions_All[3];
       fefct[3] = FEFunctions_All[4];
 
       //======================================================================
       // assembling of matrices for each level due to nonlinearity
       // A_11, (A_22)
       // no assembling of rhs
       //======================================================================
        aux =  new TAuxParam2D(MovingTNSN_FESpaces_Axial3D, MovingTNSN_Fct_Axial3D,
                               MovingTNSN_ParamFct_Axial3D,
                               MovingTNSN_FEValues_Axial3D,
                               fesp, fefct,
                               MovingTNSFct_Axial3D,
                               MovingTNSFEFctIndex_Axial3D,
                               MovingTNSFEMultiIndex_Axial3D,
                               MovingTNSN_Params_Axial3D, MovingTNSBeginParam_Axial3D);

         Assemble2D(N_FESpaces, fesp,
                    N_SquareMatrices, SQMATRICES,
                    N_RectMatrices, MATRICES,
                    N_Rhs, RHSs, ferhs,
                    DiscreteForm,
                    BoundaryConditions,
                    BoundValues,
                    aux); 
    
       // Adding freesurf entries to A11 and A22
       MatAdd(SqMat_All[4], SqMat_All[8], 1.);
       MatAdd(SqMat_All[7], SqMat_All[9], 1.); 
       

       if (TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION >= 1)
        {
          // prepare everything for the assembling of slip with friction bc
          // on all levels
          N_FESpaces = 1;
          N_SquareMatrices = 2;
          N_RectMatrices = 0;
          N_Rhs = 2;
          DiscreteForm = NULL;

          SQMATRICES[0] = SqMat_All[4];
          SQMATRICES[1] = SqMat_All[7];

          fesp[0] = FESpaces_All[0];
          ferhs[0] = FESpaces_All[0];
          ferhs[1] = FESpaces_All[0];

          RHSs[0] = Rhs_All[0];
          RHSs[1] = Rhs_All[0]+N_U;

          aux = new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);

          Assemble2DSlipBC(N_FESpaces, fesp,
                           N_SquareMatrices, SQMATRICES,
                           N_RectMatrices, MATRICES,
                           N_Rhs, RHSs, ferhs,
                           DiscreteForm,
                           BoundaryConditions,
                           BoundValues,
                           aux, FEFunctions_All[0], FEFunctions_All[1]);

       delete aux;

      } // if (TDatabase::ParamDB->INTERNA
       
     //======================================================================
     // end of assemble new matrix due to nonlinearity
     //======================================================================       
     // build stiffness matrix for next nonlinear iteration step
     // stiffness matrix (left upper block) is stored on
     // M11, (M12, M21, M22)
     // M = M +  tau*TDatabase::TimeDB->THETA1 A      
       
     MatAdd(SqMat_All[0], SqMat_All[4], tau*TDatabase::TimeDB->THETA1);
     MatAdd(SqMat_All[3], SqMat_All[7], tau*TDatabase::TimeDB->THETA1);       
     } //   for(j=0;j<Max_It;j++)   

//======================================================================
/*     // adaptive time step based on no. of nonlinear iteration
     if(j>3)
      {
       ReduceTiemStep = TRUE;
//        DecrTimeStepCount++;
       IncrTimeStepCount = 0;
      }
     else if(j<2)
      {
       IncrTimeStepCount++;
      }
     else
      {
       IncrTimeStepCount = 0;
      }  */    
   //======================================================================
   // end NSE nonlinear iteration          
   // move the grid
   //======================================================================    
  
#ifdef __SURFACT__  
   // reparam includes surf surfact interpolation
   if(reparam)
    {
     MapSurfToDomain(IFaceFeFunct[0], FEFunctions_All[6], N_List[0], N_List[1]);
     SurfSurfactReparam = TRUE;
    }
#endif 

   MoveGrid_imping(Entries, tmp_Gsol, tmp_Gd, Rhs_All[1],
                  GridKCol, GridRowPtr,
                  GridPos, FEVectFuncts_All[0], tau,
                  FEVectFuncts_All[1], AuxGridPos, 
                  MovBoundVert, N_MovVert,
                  Free_Cells, IsoCellEdgeNos,  
                  FEFunctions_All[5], FEFunctions_All[6],
                  reparam, N_ReParam);   

 #ifdef __SURFACT__    
    if(SurfSurfactReparam)
     {
      MapDomainToSurf(FEFunctions_All[6], IFaceFeFunct[0], N_List[0], N_List[1]);
      SurfSurfactReparam = FALSE;      
     }
 #endif    


	
  //======================================================================   
  // Updating the Quard points on the solid surfaces
  //======================================================================         
  GridPos->GridToData();
  RefGridPos->GridToData();

  // Updating axis boundary (top)
  MovBoundVert[1][0]->GetCoords(x, y);
  MovBoundVert[0][0]->GetCoords(SLPX, SLPY);  
 

  UpdateAxialBound->SetParams(0, y, 0., SLPY-y);

  //cout << " y " << y << " SLPY " << SLPY << endl;
  
   for(k=0;k<N_MovVert[1];k++)
    if(k==N_MovVert[1]-1)
     { Bound_Joint[0][k]->UpdateParameters(MovBoundVert[1][k], MovBoundVert[0][0]); }
    else
     { Bound_Joint[0][k]->UpdateParameters(MovBoundVert[1][k], MovBoundVert[1][k+1]); }
     
  //======================================================================          
  // end Updating the Quard points on the solid surface 
  //====================================================================== 

 
#ifdef __SURFACT__       
   //======================================================================
   // Grid velocity
   //======================================================================   
//      GridVelo_imping(Entries, tmp_Gsol, tmp_Gd, Rhs_All[1],
//                      GridKCol, GridRowPtr,
//                      GridPos, AuxGridPos,
//                      FEVectFuncts_All[0], tau,
//                      FEVectFuncts_All[1], MovBoundVert, N_MovVert,
//                      Free_Cells, IsoCellEdgeNos,
//                      reparam, RefGridPos); 

   //======================================================================
   // soluble surfactants --- begin
   //======================================================================
   // save the old time step solution
#ifdef __SOLUBLE__     
   memcpy(Csol_old, Sol_All[2], N_surfactDOF*SizeOfDouble);
#endif   
   memcpy(Isol_old, Sol_All[3], N_IsurfactDOF*SizeOfDouble);     

//    cout << " Csol_old " <<  Ddot(N_surfactDOF, Sol_All[2], Sol_All[2]) << endl;
//     cout << " Isol_old " <<  Ddot(N_IsurfactDOF, Sol_All[3], Sol_All[3]) << endl;  
   
   for(j=0;j<Max_It_scalar;j++)
   {
#ifdef __SOLUBLE__
    memcpy(Csol_nonlinearstep,  Sol_All[2], N_surfactDOF*SizeOfDouble);
    
    //======================================================================
    // surfactant in bulk phase - begin
    //======================================================================    
    // assembling matrices
    SQMATRICES_SURFACT[0] = SqMat_All[15]; // A
    SQMATRICES_SURFACT[0]->Reset();
    SQMATRICES_SURFACT[1] = SqMat_All[14]; // M
    SQMATRICES_SURFACT[1]->Reset();

    fesp[0] = FESpaces_All[3]; // surfactant space
    fesp[1] = FESpaces_All[0];  // velocity space
    fesp[2] = FESpaces_All[2];  // mesh velocity space

    ferhs[0] = FESpaces_All[3]; // surfactant space for rhs

    fefct[0] = FEFunctions_All[0]; // u1
    fefct[1] = FEFunctions_All[1]; // u2
    fefct[2] = FEFunctions_All[3]; // w1
    fefct[3] = FEFunctions_All[4]; // w2
    
    CRHSs[0] =   Rhs_All[2];

    memset(CRHSs[0], 0, N_surfactDOF*SizeOfDouble);

     // (u1-w1, u2-w2)  parameters are needed for assembling
     // fesp is taken from fefct in aux
    aux =  new TAuxParam2D(MovingTNSN_FESpaces_Axial3D, MovingTNSN_Fct_Axial3D,
                           MovingTNSN_ParamFct_Axial3D,
                           MovingTNSN_FEValues_Axial3D,
                           fesp+1, fefct,
                           MovingTNSFct_Axial3D,
                           MovingTNSFEFctIndex_Axial3D,
                           MovingTNSFEMultiIndex_Axial3D,
                           MovingTNSN_Params_Axial3D, MovingTNSBeginParam_Axial3D);

     Assemble2D(3, fesp,
                2, SQMATRICES_SURFACT,
                0, NULL,
                1, CRHSs, ferhs,
                DiscreteFormSurfact,
                SurfactBoundaryConditions,
                SurfactBoundValues,
                aux);
     delete aux;    
    

     switch(surf_couple_var)
     {
      case 1: // explicit j=0 only, no iteration
       Surfact2D_InterfaceInt(SqMat_All[15], CRHSs[0], SurfactBoundCondition, FEFunctions_All[5],
                          IFaceFeFunct[0], FESpaces_All[3], N_List[0], N_List[1], GammaXmaxVal);
      break;

      case 2: // implicit but explicit in fixed pt iteration step
      case 3:
       Surfact2D_InterfaceInt(SqMat_All[15], CRHSs[0], SurfactBoundCondition, FEFunctions_All[5],
                          IFaceFeFunct[0], FESpaces_All[3], N_List[0], N_List[1], GammaXmaxVal);
      break;

      case 4: // implicit but explicit in fixed pt iterartion step
      case 5: // fully implicit
	  IFace_Coll = (IFaceFeFunct[0]->GetFESpace1D())->GetCollection();
       Surfact2D_InterfaceInt_Implicit(SqMat_All[15], CRHSs[0], SurfactBoundCondition, FEFunctions_All[5],
                                       IFaceFeFunct[0], FESpaces_All[3], N_List[0], N_List[1], GammaXmaxVal);
      break;

      default:
       OutPut("error in selecting linerizatin type for coupled surfaact eqns " << endl);
       exit(1);
     }
    
   
   // working rhs for surfactant
   memset(C_B, 0, N_surfactDOF*SizeOfDouble);  
   
   if(j==0)
    {
     Daxpy(N_surfactActive, tau*TDatabase::TimeDB->THETA3, CRHSs[0], C_B);
   
     //for given c^Gamma_{old}
     if(TDatabase::TimeDB->THETA2>0.) 
      MatAdd(SqMat_All[14], SqMat_All[15], -tau*TDatabase::TimeDB->THETA2);
     gamma = -tau*TDatabase::TimeDB->THETA2;   
   
     memset(C_defect, 0, N_surfactDOF*SizeOfDouble);
     MatVectActive(SqMat_All[14], Csol_old, C_defect);
     Daxpy(N_surfactActive, 1, C_defect, C_B);
       
     // save the old rhs for further iterations
     memset(CRhs_old, 0, N_surfactDOF*SizeOfDouble); 
     memcpy(CRhs_old, C_B, N_surfactActive*SizeOfDouble);      
    }
   else
    {
     //copy rhs from j=0
     memcpy(C_B, CRhs_old, N_surfactActive*SizeOfDouble);
     gamma = 0.;     
    }
    
    
   // copy the currtent iterative rhs,  
   Daxpy(N_surfactActive, tau*TDatabase::TimeDB->THETA4, CRHSs[0], C_B);
 
   // set Dirichlet values
   if( (N_surfactDOF - N_surfactActive)>0)
    {    
     // set Dirichlet values
     memcpy(C_B+N_surfactActive, Rhs_All[2]+N_surfactActive, (N_surfactDOF - N_surfactActive)*SizeOfDouble);
     memcpy(Sol_All[2]+N_surfactActive, Rhs_All[2]+N_surfactActive, (N_surfactDOF - N_surfactActive)*SizeOfDouble);
    }
    
    // system matrix
    MatAdd(SqMat_All[14], SqMat_All[15], -gamma + tau*TDatabase::TimeDB->THETA1);

    // check the convergence of the fixed point linerization
    if(Max_It_scalar>1)
     {
      memset(C_defect, 0, N_surfactDOF*SizeOfDouble);
      ScalarDefect(SqMat_All[14], Sol_All[2], C_B, C_defect, residual_scalar);

      OutPut("Scalar nonlinear step " << setw(3) << j);
      OutPut(setw(14) << residual_scalar); // sqrt of residual_scalar is alread done in ScalarDefect
      if (j>0)
       { 
	 if(fabs(oldresidual_scalar)>0.)
	 OutPut(setw(14) <<  residual_scalar/oldresidual_scalar << endl); 
       }
      else
       { OutPut(endl); }

   oldresidual_scalar = residual_scalar;    
   if( ((residual_scalar<=limit_scalar)||(j==Max_It-1))  && (j>=TDatabase::ParamDB->SC_MINIT))
    { OutPut(endl); break; }
   }    
   
   // solve the system  
   DirectSolver(SqMat_All[14], C_B, Sol_All[2]);

//         exit(0);
#endif // __SOLUBLE__   
 
    //======================================================================
    //  surfactant in bulk phase - end
    //  surfactant on interphase phase - begin
    //======================================================================  
 
    //  assembling surfactant matrices
    N_SquareMatrices = 2;
    N_FESpaces = 2;
    N_Rhs =1;
    N_FESpaces_low=1;  
    
    SQMATRICES_IFace[0] = SqMat_IFace[0];
    SQMATRICES_IFace[0]->Reset();
    SQMATRICES_IFace[1] = SqMat_IFace[1];
    SQMATRICES_IFace[1]->Reset();

    fesp[0] = FESpaces_All[0]; // mesh velospace space 
    fesp[1] = FESpaces_All[3];  // outer surfactant space 

    IFacefesp[0] = IFaceFeSpaces[0]; // Interface surfactant space
    IFaceferhs[0] = IFaceFeSpaces[0]; // Interface surfactant space for rhs
 
//     //mesh velocity, even though interface moves in a Lagrangian manner due to reparam w is diff from u on interfaces
//     fefct[0] = FEFunctions_All[3]; // ur
//     fefct[1] = FEFunctions_All[4]; // uz
    
    //mesh velocity, even though interface moves in a Lagrangian manner due to reparam w is diff from u on interfaces
    fefct[0] = FEFunctions_All[0]; // ur
    fefct[1] = FEFunctions_All[1]; // uz    
    
    SRHSs[0] =  Rhs_All[3];
    memset(SRHSs[0], 0, N_IsurfactDOF*SizeOfDouble);

    // working array for srhs is I_B, initialize I_B
    memset(I_B, 0,  N_IsurfactDOF*SizeOfDouble);  
     
     if(j==0)
      {
       switch(surf_couple_var)
        {
         case 1: // all source terms on rhs, using old solution
         case 2:  
         case 3:  
           AssembleSurf1D_SolubleSurfact(N_FESpaces, fesp, fefct, N_FESpaces_low,
                   IFacefesp, IFaceFeFunct[0], N_SquareMatrices, SQMATRICES_IFace, N_Rhs, SRHSs, 
                   IFaceferhs, N_List[0], N_List[1], Csol_old);
        break;

        case 4: // c^\Gamma terms in source are in LHS, but the solution  
        case 5: // c^\Gamma terms in source are in LHS, solution (c_old) 
          AssembleSurf1D_SolubleSurfact_Implicit(N_FESpaces, fesp, fefct, N_FESpaces_low,
                     IFacefesp, IFaceFeFunct[0], N_SquareMatrices, SQMATRICES_IFace, N_Rhs, SRHSs, 
                     IFaceferhs, N_List[0], N_List[1], Csol_old);
        break;
        default:
         OutPut("error in selecting linerizatin type for coupled surfaact eqns " << endl);
         exit(1);
        }     
        
      // copy the rhs, for given c^{old}
      Daxpy(N_IActive, tau*TDatabase::TimeDB->THETA3, SRHSs[0], I_B);     
     
      //for given c^{old},
      MatAdd(SqMat_IFace[1], SqMat_IFace[0], -tau*TDatabase::TimeDB->THETA2);
      memset(I_defect, 0, N_IsurfactDOF*SizeOfDouble);
      MatVectActive(SqMat_IFace[1], Isol_old, I_defect);
      Daxpy(N_IActive, 1, I_defect, I_B);  
      
      // save it for further iterative steps
      memset(IRhs_old, 0, N_IsurfactDOF*SizeOfDouble);
      memcpy(IRhs_old, I_B, N_IActive*SizeOfDouble);
      
      //reset the mat and rhs for assembling
      SQMATRICES_IFace[0]->Reset();
      SQMATRICES_IFace[1]->Reset();
      memset(SRHSs[0], 0, N_IsurfactDOF*SizeOfDouble);
     }
    else
     {
      // copy the rhs from j==0
      memcpy(I_B, IRhs_old, N_IActive*SizeOfDouble);
     }  // if( j==0 
      
     // assemble with new iterative solution
     switch(surf_couple_var)
     {
      case 1: // fimplicit but explicit in fixed pt iterartion step
       AssembleSurf1D_SolubleSurfact(N_FESpaces, fesp, fefct, N_FESpaces_low,
                   IFacefesp, IFaceFeFunct[0], N_SquareMatrices, SQMATRICES_IFace, N_Rhs, SRHSs, 
                   IFaceferhs, N_List[0], N_List[1], Csol_old);
       memcpy(IRhs_old, SRHSs[0], N_IsurfactDOF*SizeOfDouble);
      break;

      case 2: // fimplicit but explicit in fixed pt iterartion step
       AssembleSurf1D_SolubleSurfact(N_FESpaces, fesp, fefct, N_FESpaces_low,
                   IFacefesp, IFaceFeFunct[0], N_SquareMatrices, SQMATRICES_IFace, N_Rhs, SRHSs, 
                   IFaceferhs, N_List[0], N_List[1], Csol_nonlinearstep);

      break;

      case 3: // fimplicit but explicit in fixed pt iterartion step
       AssembleSurf1D_SolubleSurfact(N_FESpaces, fesp, fefct, N_FESpaces_low,
                   IFacefesp, IFaceFeFunct[0], N_SquareMatrices, SQMATRICES_IFace, N_Rhs, SRHSs, 
                   IFaceferhs, N_List[0], N_List[1], Sol_All[2]);

      break;

      case 4: // implicit but explicit in fixed pt iterartion step
       AssembleSurf1D_SolubleSurfact_Implicit(N_FESpaces, fesp, fefct, N_FESpaces_low,
                   IFacefesp, IFaceFeFunct[0], N_SquareMatrices, SQMATRICES_IFace, N_Rhs, SRHSs, 
                   IFaceferhs, N_List[0], N_List[1], Csol_nonlinearstep);
      break;

      case 5: // fully implicit
       AssembleSurf1D_SolubleSurfact_Implicit(N_FESpaces, fesp, fefct, N_FESpaces_low,
                     IFacefesp, IFaceFeFunct[0], N_SquareMatrices, SQMATRICES_IFace, N_Rhs, SRHSs, 
                     IFaceferhs, N_List[0], N_List[1], Sol_All[2]);
      break;
      default:
       OutPut("error in selecting linerizatin type for coupled surfaact eqns " << endl);
       exit(1);
     }
     
    // copy the currtent iterative rhs, for given c^{n+1}
    Daxpy(N_IActive, tau*TDatabase::TimeDB->THETA4,  SRHSs[0], I_B);

    // set Dirichlet values
    if( (N_IsurfactDOF - N_IActive)>0)
     {
      memcpy(I_B+N_IActive, Rhs_All[3]+N_IActive, (N_IsurfactDOF - N_IActive)*SizeOfDouble);
      memcpy(Sol_All[3]+N_IActive, Rhs_All[3]+N_IActive, (N_IsurfactDOF - N_IActive)*SizeOfDouble);
     }
 
   // assembling of system matrix
//    MatAdd(SqMat_IFace[1], SqMat_IFace[0], tau*TDatabase::TimeDB->THETA1);
              
   //solve the system
   DirectSolver(SqMat_IFace[1], I_B, Sol_All[3]);
   
   //======================================================================
   // surfactant on interphase phase - end
   //======================================================================
   
  } //  for(j=0;j<Max_It_scalar;j++)  


     
//  if((m % (int)(TDatabase::TimeDB->STEPS_PER_IMAGE) ) == 0   || m==1 )
//   {
//      if(TDatabase::ParamDB->WRITE_VTK)
//        { 
//        // interface surfactant
//        MapSurfToDomain(IFaceFeFunct[0], FEFunctions_All[6], N_List[0], N_List[1]);
//   
//        GetSurfactMass(FEFunctions_All[6], IFaceFeFunct[0], N_List[0], N_List[1], Surf_Mass);
// 
//        Get_FeFunction2DMass(FEFunctions_All[5], Params);        
//      
//      
//         OutPut( "Time, GammaMaxR, GammaMax " <<TDatabase::TimeDB->CURRENTTIME<<
//               " " <<GammaXmaxVal[0]<< " "<<GammaXmaxVal[1]<< " "<<endl);   
//       OutPut( "Time, Surfactant_Mass, Da*InterfaceSurfactant_Mass, InterfaceSurfactant_Conc " <<TDatabase::TimeDB->CURRENTTIME<<
//               " " <<Params[0]<< " "<< TDatabase::ParamDB->REACTOR_P14*Surf_Mass[0]<< " "<<Surf_Mass[0]/Surf_Mass[1]<< " "<<endl);     
//       OutPut( "Time, Surfactant Mass dif, InterfaceSurfactant Mass diff " <<TDatabase::TimeDB->CURRENTTIME<<
//               " " <<(Params[0] - Initial_SurfactMass) << " "<<(TDatabase::ParamDB->REACTOR_P14*Surf_Mass[0] - Initial_IFaceSurfactMass) << " "<<endl);
//       OutPut( "Time, Total Mass, Total Mass diff, RelativeTotal Mass diff " <<TDatabase::TimeDB->CURRENTTIME<<
//           " " <<Params[0]+TDatabase::ParamDB->REACTOR_P14*Surf_Mass[0]<< " " << ((Params[0]+TDatabase::ParamDB->REACTOR_P14*Surf_Mass[0]) - 
//             (Initial_IFaceSurfactMass + Initial_SurfactMass)) << " " << ((Params[0]+TDatabase::ParamDB->REACTOR_P14*Surf_Mass[0]) - 
//             (Initial_IFaceSurfactMass + Initial_SurfactMass))
//             /(Initial_IFaceSurfactMass + Initial_SurfactMass) <<endl);     
//      
//      
//         os.seekp(std::ios::beg);
//         if(img<10) os << "VTK/"<<VtkBaseName<<".0000"<<img<<".vtk" << ends;
//          else if(img<100) os <<"VTK/"<< VtkBaseName<<".000"<<img<<".vtk" << ends;
//           else if(img<1000) os <<"VTK/"<< VtkBaseName<<".00"<<img<<".vtk" << ends;
//            else if(img<10000) os << "VTK/"<<VtkBaseName<<".0"<<img<<".vtk" << ends;
//             else  os << "VTK/"<<VtkBaseName<<"."<<img<<".vtk" << ends;
//         Output->WriteVtk(os.str().c_str());
//         img++;
//        }       
//    }
  
//   exit(0);
  
  
  
  
  
/*     if(TDatabase::ParamDB->WRITE_VTK)
       { 
       // interface surfactant
      MapSurfToDomain(IFaceFeFunct[0], FEFunctions_All[6], N_List[0], N_List[1]);
 
        os.seekp(std::ios::beg);
        if(img<10) os << "VTK/"<<VtkBaseName<<".0000"<<img<<".vtk" << ends;
         else if(img<100) os <<"VTK/"<< VtkBaseName<<".000"<<img<<".vtk" << ends;
          else if(img<1000) os <<"VTK/"<< VtkBaseName<<".00"<<img<<".vtk" << ends;
           else if(img<10000) os << "VTK/"<<VtkBaseName<<".0"<<img<<".vtk" << ends;
            else  os << "VTK/"<<VtkBaseName<<"."<<img<<".vtk" << ends;
        Output->WriteVtk(os.str().c_str());
        img++;
       }    */   
// 
//   cout << " test main " << endl; 
//   exit(0);   
     
//      
#endif       
  //======================================================================       
  // Remeshing Begin 
  //======================================================================      
//    if((l==0) && ((m % 1) == 0))
    {
     Getcellangle(FESpaces_All[2], Angle);

    }
 
//  TDatabase::ParamDB->P15=1;

   if((Angle[0]<10.0) ||(Angle[1]>165.0))
    {
      OutPut( "MinAngle : "<< Angle[0]<< "  MaxAngle : "<<Angle[1]<< endl);     
//      if(TDatabase::ParamDB->P10)
//      {
//       ReParam_axial3D_U(N_MovVert[0], Free_Cells,  IsoCellEdgeNos[1], IsoCellEdgeNos[0],  FEVectFuncts_All[0]);      
//      }
//   GetMovingBoundData(coll, N_MovVert, Bound_Joint, MovBoundVert, Free_Joint,
//                      Free_Cells, IsoCellEdgeNos, S_BX[0], S_BY[0]);
// 
//       N_E = N_MovVert[0];
//       ReParam_axial3D_Data(N_E, Free_Cells,  IsoCellEdgeNos[1], IsoCellEdgeNos[0],  FEVectFuncts_All[0], 
//                            Intpol_Coord,  Intpol_VeloValues, h_interface);
//       reparam = TRUE;
//       cout<< "Old N_E " << N_MovVert[0] <<"New N_E " << N_E <<endl;     
//       N_MovVert[0] = N_E;
//       
      
//       exit(0);
     
#ifdef __SURFACT__   
        // interface surfactant
       MapSurfToDomain(IFaceFeFunct[0], FEFunctions_All[6], N_List[0], N_List[1]);
#endif    
       
      if(TDatabase::ParamDB->WRITE_VTK)
       { 
        os.seekp(std::ios::beg);
        if(RemeshImg<10) os << "VTK/"<<VtkBaseName<<"_remesh.0000"<<RemeshImg<<".vtk" << ends;
         else if(RemeshImg<100) os << "VTK/"<<VtkBaseName<<"_remesh.000"<<RemeshImg<<".vtk" << ends;
          else if(RemeshImg<1000) os << "VTK/"<<VtkBaseName<<"_remesh.00"<<RemeshImg<<".vtk" << ends;
           else if(RemeshImg<10000) os <<"VTK/"<< VtkBaseName<<"_remesh.0"<<RemeshImg<<".vtk" << ends;
            else  os << "VTK/"<<VtkBaseName<<"_remesh."<<RemeshImg<<".vtk" << ends;
        Output->WriteVtk(os.str().c_str());
        RemeshImg++;
       }  


//       t1 = GetTime();
     RemeshAxial3D(Domain, FESpaces_All, FEVectFuncts_All, FEFunctions_All,
                   N_MovVert, Bound_Joint, MovBoundVert, Free_Joint, Free_Cells, IsoCellEdgeNos,
                   Sol_All, Rhs_All, SquareStructure_All, Structure_All, SqMat_All, Mat_All
#ifdef __SURFACT__                     
                   ,IFaceDomain, IFaceFeSpaces, N_List,
                   IFaceFeFunct, SqMat_IFace,
                   IFaceStruct, FE1D_List  
#endif                   
                   );
      
#ifdef __SURFACT__     
     // copy to interface after bulk interpolation
     MapDomainToSurf(FEFunctions_All[6], IFaceFeFunct[0], N_List[0], N_List[1]);  
     N_S =  FESpaces_All[4]->GetN_DegreesOfFreedom();        
     memset(Sol_All[4], 0, N_S*SizeOfDouble);    
#endif      

     coll=FESpaces_All[0]->GetCollection();
     N_Cells = coll->GetN_Cells();
  
     Getcellangle(FESpaces_All[0], Angle);
      OutPut( "MinAngle : "<< Angle[0]<< "  MaxAngle : "<<Angle[1]<< endl);

     N_Remesh ++;
     remeshed=TRUE;

     GlobalNumbers = FESpaces_All[0]->GetGlobalNumbers();
     BeginIndex = FESpaces_All[0]->GetBeginIndex();
     N_Active =  FESpaces_All[0]->GetActiveBound();
     N_U = FESpaces_All[0]->GetN_DegreesOfFreedom();
     N_P = FESpaces_All[1]->GetN_DegreesOfFreedom();
     N_Unknowns = 2*N_U + N_P;
  
     N_G = FESpaces_All[2]->GetN_DegreesOfFreedom();
     N_GActive = FESpaces_All[2]->GetActiveBound();
     N_GBoundaryNodes = N_G - N_GActive;     
 
     delete [] B; delete [] defect; 
     B = new double[N_Unknowns];
     defect = new double[N_Unknowns];

#ifdef __SURFACT__    
     N_surfactDOF = FESpaces_All[3]->GetN_DegreesOfFreedom();
     N_surfactActive = FESpaces_All[3]->GetActiveBound();
     N_surfactNonActive = N_surfactDOF - N_surfactActive;  
     N_IsurfactDOF = IFaceFeSpaces[0]->GetN_DegreesOfFreedom();  
     N_IActive = IFaceFeSpaces[0]->GetActiveBound();
     
     delete [] Csol_old; 
     delete [] CRhs_old;
     delete [] C_B;
     delete [] C_defect;
     delete [] Csol_nonlinearstep;
     delete [] I_defect;
     delete [] Isol_old;
     delete [] I_B;
     delete [] IRhs_old;
     
   Csol_old = new double[N_surfactDOF];
   CRhs_old = new double[N_surfactDOF];
   C_B = new double[N_surfactDOF];  
   C_defect = new double[N_surfactDOF];  
   Csol_nonlinearstep = new double[N_surfactDOF];
   I_defect = new double[N_IsurfactDOF];
   Isol_old = new double[N_IsurfactDOF];
   I_B = new double[N_IsurfactDOF];
   IRhs_old = new double[N_IsurfactDOF];     
   
   memset(Csol_old, 0, N_surfactDOF*SizeOfDouble); 
   memset(CRhs_old, 0, N_surfactDOF*SizeOfDouble);   
   memset(C_defect, 0, N_surfactDOF*SizeOfDouble); 
   memset(Csol_nonlinearstep, 0, N_surfactDOF*SizeOfDouble); 
   memset(Csol_old, 0, N_surfactDOF*SizeOfDouble); 
   memset(I_defect, 0, N_IsurfactDOF*SizeOfDouble); 
   memset(Isol_old, 0, N_IsurfactDOF*SizeOfDouble); 
   memset(I_B, 0, N_IsurfactDOF*SizeOfDouble);    
   memset(IRhs_old, 0, N_IsurfactDOF*SizeOfDouble);      
#endif       
     
     
     GridKCol = SquareStructure_All[1]->GetKCol();
     GridRowPtr = SquareStructure_All[1]->GetRowPtr();
          
     delete [] refpos; delete []  auxpos;  delete [] pos;   delete [] ReparamMeshVelo;
     delete [] ReparamDisp; delete [] tmp_Gsol; delete [] tmp_Gd;

     tmp_Gd = new double[2*N_G];     
     tmp_Gsol = new double[2*N_G];
     refpos = new double[2*N_G];
     auxpos = new double[2*N_G];
     pos = new double[2*N_G];
     ReparamMeshVelo = new double[2*N_G];
     ReparamDisp = new double[2*N_G]; 
     
     delete RefGridPos; delete AuxGridPos; delete GridPos; delete ReparamPos;  

     ReparamPos = new TFEVectFunct2D(FESpaces_All[2], WString, WString, ReparamDisp, N_G, 2);       
     RefGridPos = new TFEVectFunct2D(FESpaces_All[2], WString, WString, refpos, N_G, 2);
     AuxGridPos = new TFEVectFunct2D(FESpaces_All[2], WString, WString, auxpos, N_G, 2);
     GridPos = new TFEVectFunct2D(FESpaces_All[2], WString, WString, pos, N_G, 2);

     delete Output;
     // prepare output (maxn_fespaces,  maxn_scalar,  maxn_vect, maxn_parameters, domain)
     Output = new TOutput2D(1, 3, 1, 2, Domain);
     Output->AddFEVectFunct(FEVectFuncts_All[0]);
     Output->AddFEFunction(FEFunctions_All[2]);      
       
#ifdef __SURFACT__      
    Output->AddFEFunction(FEFunctions_All[5]);
  
    Output->AddFEFunction(FEFunctions_All[6]);
#endif   
  
      if(TDatabase::ParamDB->WRITE_VTK)
       { 
#ifdef __SURFACT__ 	 
        // interface surfactant
        MapSurfToDomain(IFaceFeFunct[0], FEFunctions_All[6], N_List[0], N_List[1]);
#endif  
        os.seekp(std::ios::beg);
        if(RemeshImg<10) os << "VTK/"<<VtkBaseName<<"_remesh.0000"<<RemeshImg<<".vtk" << ends;
         else if(RemeshImg<100) os <<"VTK/"<< VtkBaseName<<"_remesh.000"<<RemeshImg<<".vtk" << ends;
          else if(RemeshImg<1000) os <<"VTK/"<< VtkBaseName<<"_remesh.00"<<RemeshImg<<".vtk" << ends;
           else if(RemeshImg<10000) os <<"VTK/"<< VtkBaseName<<"_remesh.0"<<RemeshImg<<".vtk" << ends;
            else  os << "VTK/"<<VtkBaseName<<"_remesh."<<RemeshImg<<".vtk" << ends;
        Output->WriteVtk(os.str().c_str());
        RemeshImg++;
       }  
//             exit(0);

    } // if((Angle[0]<10.0) ||(Angle[1]>165.0))
  //======================================================================  
  // end Remeshing Begin 
  // Assembeling the grid matrix - Begin
  //======================================================================  

    fesp[0] = FESpaces_All[2];
    SQMATRICES_GRID[0] = SqMat_All[10];
    SQMATRICES_GRID[0]->Reset();
    SQMATRICES_GRID[1] = SqMat_All[11];
    SQMATRICES_GRID[1]->Reset();
    SQMATRICES_GRID[2] = SqMat_All[12];
    SQMATRICES_GRID[2]->Reset();
    SQMATRICES_GRID[3] = SqMat_All[13];
    SQMATRICES_GRID[3]->Reset();
    aux = new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);  
       
    Assemble2D(1, fesp,
             4, SQMATRICES_GRID,
             0, NULL,
             0, NULL, NULL,
             DiscreteFormGrid,
             GridBoundaryConditions,
             GridBoundValues,
             aux);
    delete aux;   
   
     Entries[0] = SqMat_All[10]->GetEntries();
     Entries[1] = SqMat_All[11]->GetEntries();
     Entries[2] = SqMat_All[12]->GetEntries();
     Entries[3] = SqMat_All[13]->GetEntries();

   // for Dirichlet rows in off-diagonal matrices
   memset(Entries[1] + GridRowPtr[N_GActive], 0, (GridRowPtr[N_G] - GridRowPtr[N_GActive])*SizeOfDouble);
   memset(Entries[2] + GridRowPtr[N_GActive], 0, (GridRowPtr[N_G] - GridRowPtr[N_GActive])*SizeOfDouble);  
   
  //======================================================================  
  // end Assembeling the grid matrix
  // nonlinear loop due to remeshing
  //first update the convection matrix with the new mesh velocity        
  //======================================================================
   if(remeshed)
    {           
      // working array for rhs is B, initialize B
      memset(B, 0, N_Unknowns*SizeOfDouble);
      memset(Sol_All[1], 0, 2*N_G*SizeOfDouble);  
     
        DiscreteForm = DiscreteFormGalerkin;

        SQMATRICES[0] = SqMat_All[4];
        SQMATRICES[1] = SqMat_All[5];
        SQMATRICES[2] = SqMat_All[6];
        SQMATRICES[3] = SqMat_All[7];
        SQMATRICES[4] = SqMat_All[0];
        SQMATRICES[5] = SqMat_All[3];

        SQMATRICES[6] = SqMat_All[1];
        SQMATRICES[7] = SqMat_All[2];

        MATRICES[0] = Mat_All[0];
        MATRICES[1] = Mat_All[1];
        MATRICES[2] = Mat_All[2];
        MATRICES[3] = Mat_All[3];

        SQMATRICES[0]->Reset();
        SQMATRICES[1]->Reset();
        SQMATRICES[2]->Reset();
        SQMATRICES[3]->Reset();
        SQMATRICES[4]->Reset();
        SQMATRICES[5]->Reset();
        SQMATRICES[6]->Reset();
        SQMATRICES[7]->Reset();

        MATRICES[0]->Reset();
        MATRICES[1]->Reset();
        MATRICES[2]->Reset();
        MATRICES[3]->Reset();

        N_SquareMatrices = 6;
        N_RectMatrices = 4;

       // parameters which are the same for all NSTYPEs
        N_Rhs = 2;
        N_FESpaces = 3;

        fesp[0] = FESpaces_All[0];
        fesp[1] = FESpaces_All[1];
        fesp[2] = FESpaces_All[2];

        fefct[0] = FEFunctions_All[0];
        fefct[1] = FEFunctions_All[1];
        fefct[2] = FEFunctions_All[3];
        fefct[3] = FEFunctions_All[4];
 
        ferhs[0] = FESpaces_All[0];
        ferhs[1] = FESpaces_All[0];

        RHSs[0] = Rhs_All[0];
        RHSs[1] = Rhs_All[0] + N_U;
        RHSs[2] = Rhs_All[0] + 2*N_U;

        memset(Rhs_All[0], 0, N_Unknowns*SizeOfDouble);

       // 4 parameters are needed for assembling (u1_old, u2_old)
        aux =  new TAuxParam2D(MovingTNSN_FESpaces_Axial3D, MovingTNSN_Fct_Axial3D,
                               MovingTNSN_ParamFct_Axial3D,
                               MovingTNSN_FEValues_Axial3D,
                               fesp, fefct,
                               MovingTNSFct_Axial3D,
                               MovingTNSFEFctIndex_Axial3D,
                               MovingTNSFEMultiIndex_Axial3D,
                               MovingTNSN_Params_Axial3D, MovingTNSBeginParam_Axial3D);

       //======================================================================
       // assembling of matrices for each level
       // A_11 , (A_12), (A_21), (A_22), M_11, (M_22)
       //======================================================================
       Assemble2D(N_FESpaces, fesp,
                  N_SquareMatrices, SQMATRICES,
                  N_RectMatrices, MATRICES,
                  N_Rhs, RHSs, ferhs,
                  DiscreteForm,
                  BoundaryConditions,
                  BoundValues,
                  aux);

      delete aux;

     SqMat_All[8]->Reset(); // Matrix entries for freesurf int;
     SqMat_All[9]->Reset(); // no need to calculate in nonlinear steps
     
#ifdef __SURFACT__ 
     MapSurfToDomain(IFaceFeFunct[0], FEFunctions_All[6], N_List[0], N_List[1]);      
#endif     
     
     FreeSurf_axial3D_new(SqMat_All[8], SqMat_All[9],  RHSs[0], RHSs[1], BoundCondition, tau,
                          FEFunctions_All[0]->GetValues(), Params
#ifdef __SURFACT__  
                          , FEFunctions_All[6] 
#endif                               
    );
    
     // Adding freesurf entries to A11 and A22
     MatAdd(SqMat_All[4], SqMat_All[8], 1);
     MatAdd(SqMat_All[7], SqMat_All[9], 1);

     // set rows of Dirichlet dof in off diagonal matrix blocks
     // to zero    
     // N_Active =  FESpaces_All[0]->GetActiveBound();
     // get row in off diagonal matrix where the Dirichlet nodes start
     RowPtr = SqMat_All[5]->GetRowPtr();
     // compute number of entries starting from this row to the end
     // of the matrix
     j = RowPtr[N_Active];
     k = RowPtr[N_U]-j;
     // get number of active dof
     // set these entries to zero
     memset(SqMat_All[5]->GetEntries()+j, 0, SizeOfDouble*k);
     memset(SqMat_All[6]->GetEntries()+j, 0, SizeOfDouble*k);

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
          N_SquareMatrices = 8;
          N_RectMatrices = 2;
          N_Rhs = 2;
          DiscreteForm = NULL;

          SQMATRICES[0] = SqMat_All[4];
          SQMATRICES[1] = SqMat_All[7];
          SQMATRICES[2] = SqMat_All[5];
          SQMATRICES[3] = SqMat_All[6];
          SQMATRICES[4] = SqMat_All[0];
          SQMATRICES[5] = SqMat_All[3];
          SQMATRICES[6] = SqMat_All[1];
          SQMATRICES[7] = SqMat_All[2];

          MATRICES[0] = Mat_All[2];
          MATRICES[1] = Mat_All[3];

          fesp[0] = FESpaces_All[0];
          ferhs[0] = FESpaces_All[0];
          ferhs[1] = FESpaces_All[0];

          RHSs[0] = Rhs_All[0];
          RHSs[1] = Rhs_All[0]+N_U;

          aux = new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);

          Assemble2DSlipBC(N_FESpaces, fesp,
                           N_SquareMatrices, SQMATRICES,
                           N_RectMatrices, MATRICES,
                           N_Rhs, RHSs, ferhs,
                           DiscreteForm,
                           BoundaryConditions,
                           BoundValues,
                           aux, FEFunctions_All[0], FEFunctions_All[1]);

       delete aux;

      } // if (TDatabase::ParamDB->INTERNA

     //    scale the pressure matrices
     Dscal(Mat_All[2]->GetN_Entries(), tau, Mat_All[2]->GetEntries());
     Dscal(Mat_All[3]->GetN_Entries(), tau, Mat_All[3]->GetEntries());
     Dscal(Mat_All[0]->GetN_Entries(), tau, Mat_All[0]->GetEntries());
     Dscal(Mat_All[1]->GetN_Entries(), tau, Mat_All[1]->GetEntries());
     
     // update rhs
     Daxpy(N_Active, tau, Rhs_All[0], B);
     Daxpy(N_Active, tau, Rhs_All[0]+N_U, B+N_U);

     // update rhs by Laplacian and convective term initialy by current time step
     // scaled by current sub time step length and theta2
     // currently : M := M + gamma A
     // M = M + (-gamma - tau*TDatabase::TimeDB->THETA2) A
     MatAdd(SqMat_All[0], SqMat_All[4], - tau*TDatabase::TimeDB->THETA2);
     MatAdd(SqMat_All[1], SqMat_All[5], - tau*TDatabase::TimeDB->THETA2);
     MatAdd(SqMat_All[2], SqMat_All[6], - tau*TDatabase::TimeDB->THETA2);
     MatAdd(SqMat_All[3], SqMat_All[7], - tau*TDatabase::TimeDB->THETA2);

     // set current factor of steady state matrix
     gamma = -tau*TDatabase::TimeDB->THETA2;		     

     // defect = M * Sol
     // B:= B + defect (rhs)     
     MatVectActive(SqMat_All[0], Sol_All[0], defect);
     Daxpy(N_Active, 1, defect, B);
     MatVectActive(SqMat_All[1], Sol_All[0]+N_U, defect);
     Daxpy(N_Active, 1, defect, B);
     MatVectActive(SqMat_All[2], Sol_All[0], defect+N_U);
     Daxpy(N_Active, 1, defect+N_U, B+N_U);
     MatVectActive(SqMat_All[3], Sol_All[0]+N_U, defect+N_U);
     Daxpy(N_Active, 1, defect+N_U, B+N_U);     
     
    // set Dirichlet values
    // RHSs[0] still available from assembling
    memcpy(B+N_Active, RHSs[0]+N_Active, (N_U-N_Active)*SizeOfDouble);
    memcpy(B+N_Active+N_U, RHSs[1]+N_Active,(N_U-N_Active)*SizeOfDouble);

    // copy Dirichlet values from rhs into Sol[0][mg_level-1]
    memcpy(Sol_All[0]+N_Active, RHSs[0]+N_Active, (N_U-N_Active)*SizeOfDouble);
    memcpy(Sol_All[0]+N_Active+N_U, RHSs[1]+N_Active, (N_U-N_Active)*SizeOfDouble);

    //=====================================================================
    // the stiffness matrix is stored on M11, (M12, M21, M22)
    // assembling of system matrix
    //========================================================================
    // M = M + (-gamma + tau*TDatabase::TimeDB->THETA1) A
     MatAdd(SqMat_All[0], SqMat_All[4], -gamma + tau*TDatabase::TimeDB->THETA1);
     MatAdd(SqMat_All[1], SqMat_All[5], -gamma + tau*TDatabase::TimeDB->THETA1);
     MatAdd(SqMat_All[2], SqMat_All[6], -gamma + tau*TDatabase::TimeDB->THETA1);
     MatAdd(SqMat_All[3], SqMat_All[7], -gamma + tau*TDatabase::TimeDB->THETA1);

     // set current factor of steady state matrix
     gamma = tau*TDatabase::TimeDB->THETA1;     
  
     //======================================================================
     // nonlinear loop
     //======================================================================   
     N_LinIterCurr = 0;
     solver_time_curr = 0;      

     for(j=0;j<Max_It;j++)
      {
       memset(defect, 0, N_Unknowns*SizeOfDouble);

       SQMATRICES[0] = SqMat_All[0];
       SQMATRICES[1] = SqMat_All[1];
       SQMATRICES[2] = SqMat_All[2];
       SQMATRICES[3] = SqMat_All[3];
       MATRICES[0] = Mat_All[0];
       MATRICES[1] = Mat_All[1];
       MATRICES[2] = Mat_All[2];
       MATRICES[3] = Mat_All[3];      
       
      // compute defect
      Defect(sqmatrices, matrices, Sol_All[0], B, defect);

      residual =  Ddot(N_Unknowns, defect, defect);
      impuls_residual = Ddot(2*N_U, defect, defect);
      OutPut("nonlinear step " << setw(3) << j);
      OutPut(setw(14) << impuls_residual);
      OutPut(setw(14) << Ddot(N_P,defect+2*N_U,defect+2*N_U));
      OutPut(setw(14) << sqrt(residual));
      
      if(j>0)
       {
        OutPut(setw(14) << sqrt(residual)/oldresidual << endl);
       }
      else
       {
        OutPut(endl);
       }
       
      oldresidual = sqrt(residual);

      if ((((sqrt(residual)<=limit)||(j==Max_It-1)))  && (j>=TDatabase::ParamDB->SC_MINIT))
       {
        if (j==Max_It-1)
        j++;
        OutPut("ITE : " << setw(3) << j);
        OutPut(" (" << setw(3) << N_LinIterCurr << "/");
        OutPut(setw(3) << N_LinIter << " LINITE)");
        OutPut("  TIME FOR SOLVER : " << solver_time_curr << "/" << solver_time << "s");
        OutPut("  RES : " <<  sqrt(residual) << endl);
        // count total running time
        t4 =  GetTime();
        total_time += t4 - t3;
        t3 = t4;
        OutPut(TDatabase::TimeDB->CURRENTTIME << " total elapsed running time "<< total_time << endl);
        break;
       }

       //======================================================================
       // solve linear system
       //======================================================================
        t1 = GetTime();
        DirectSolver(SQMATRICES[0], SQMATRICES[1], SQMATRICES[2], SQMATRICES[3],
                     MATRICES[2], MATRICES[3], MATRICES[0], MATRICES[1],
                     B, Sol_All[0]);
        t2 = GetTime();
        solver_time_curr = t2-t1;
        solver_time += solver_time_curr;
 
       //======================================================================
       // end solve linear system
       //======================================================================
       // restore mass matrices by subtracting the A-matrices
       MatAdd(SqMat_All[0], SqMat_All[4], -gamma);
       MatAdd(SqMat_All[3], SqMat_All[7], -gamma);

       //======================================================================
       // assemble new matrix due to nonlinearity
       //======================================================================
       DiscreteForm = DiscreteFormNLGalerkin;	 
       N_RectMatrices = 0;
       N_Rhs = 0;
       N_FESpaces = 3;

       SQMATRICES[0] = SqMat_All[4];
       SQMATRICES[1] = SqMat_All[7];
       SQMATRICES[0]->Reset();
       SQMATRICES[1]->Reset();

       N_SquareMatrices = 2;
       last_sq = 1;
       
       fesp[0] = FESpaces_All[0];
       fesp[1] = FESpaces_All[1];
       fesp[2] = FESpaces_All[2];

       fefct[0] = FEFunctions_All[0];
       fefct[1] = FEFunctions_All[1];
       fefct[2] = FEFunctions_All[3];
       fefct[3] = FEFunctions_All[4];
 
       //======================================================================
       // assembling of matrices for each level due to nonlinearity
       // A_11, (A_22), no assembling of rhs
       //======================================================================
        aux =  new TAuxParam2D(MovingTNSN_FESpaces_Axial3D, MovingTNSN_Fct_Axial3D,
                               MovingTNSN_ParamFct_Axial3D,
                               MovingTNSN_FEValues_Axial3D,
                               fesp, fefct,
                               MovingTNSFct_Axial3D,
                               MovingTNSFEFctIndex_Axial3D,
                               MovingTNSFEMultiIndex_Axial3D,
                               MovingTNSN_Params_Axial3D, MovingTNSBeginParam_Axial3D);

         Assemble2D(N_FESpaces, fesp,
                    N_SquareMatrices, SQMATRICES,
                    N_RectMatrices, MATRICES,
                    N_Rhs, RHSs, ferhs,
                    DiscreteForm,
                    BoundaryConditions,
                    BoundValues,
                    aux); 
    
       // Adding freesurf entries to A11 and A22
       MatAdd(SqMat_All[4], SqMat_All[8], 1.);
       MatAdd(SqMat_All[7], SqMat_All[9], 1.); 
       

       if (TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION >= 1)
        {
          // prepare everything for the assembling of slip with friction bc
          // on all levels
          N_FESpaces = 1;
          N_SquareMatrices = 2;
          N_RectMatrices = 0;
          N_Rhs = 2;
          DiscreteForm = NULL;

          SQMATRICES[0] = SqMat_All[4];
          SQMATRICES[1] = SqMat_All[7];

          fesp[0] = FESpaces_All[0];
          ferhs[0] = FESpaces_All[0];
          ferhs[1] = FESpaces_All[0];

          RHSs[0] = Rhs_All[0];
          RHSs[1] = Rhs_All[0]+N_U;

          aux = new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);

          Assemble2DSlipBC(N_FESpaces, fesp,
                           N_SquareMatrices, SQMATRICES,
                           N_RectMatrices, MATRICES,
                           N_Rhs, RHSs, ferhs,
                           DiscreteForm,
                           BoundaryConditions,
                           BoundValues,
                           aux, FEFunctions_All[0], FEFunctions_All[1]);

       delete aux;

      } // if (TDatabase::ParamDB->INTERN    
       
     MatAdd(SqMat_All[0], SqMat_All[4], tau*TDatabase::TimeDB->THETA1);
     MatAdd(SqMat_All[3], SqMat_All[7], tau*TDatabase::TimeDB->THETA1);  


//         cout <<  l << "temain " <<endl; 
    } //   for(j=0;j<Max_It;j++)   
            
    remeshed = FALSE;
        
    
/*      if(TDatabase::ParamDB->WRITE_VTK)
       { 
#ifdef __SURFACT__ 	 
        // interface surfactant
        MapSurfToDomain(IFaceFeFunct[0], FEFunctions_All[6], N_List[0], N_List[1]);
#endif  
        os.seekp(std::ios::beg);
        if(img<10) os << "VTK/"<<VtkBaseName<<"_remesh.0000"<<RemeshImg<<".vtk" << ends;
         else if(img<100) os <<"VTK/"<< VtkBaseName<<"_remesh.000"<<RemeshImg<<".vtk" << ends;
          else if(img<1000) os <<"VTK/"<< VtkBaseName<<"_remesh.00"<<RemeshImg<<".vtk" << ends;
           else if(img<10000) os <<"VTK/"<< VtkBaseName<<"_remesh.0"<<RemeshImg<<".vtk" << ends;
            else  os << "VTK/"<<VtkBaseName<<"_remesh."<<RemeshImg<<".vtk" << ends;
        Output->WriteVtk(os.str().c_str());
        RemeshImg++;
       }   */   
   } //if (remeshed)  
//     cout <<  l << "testmain " <<endl;
 
  }// for(l=0;l<N_SubSteps;l++) 


 

    MovBoundVert[0][0]->GetCoords(x1, y1);
    
     tn_2 = tn_1;
     tn_1 = tn;
     tn   = TDatabase::TimeDB->CURRENTTIME;
    
     ztn_2 = ztn_1;
     ztn_1 = ztn;
     ztn   = y1;
    
     if( ((ztn_2< ztn_1)&&(ztn_1>ztn)) ||  ((ztn_2> ztn_1)&&(ztn_1<ztn))  )
      {
       OutPut(setw(25)<<"T, Bottom : "<< TDatabase::TimeDB->CURRENTTIME<<"   "<< ztn_1 <<  endl);  
   
       if( (ztn_2< ztn_1)&&(ztn_1>ztn) )
        { OutPut(setw(25)<<"T, Bott_Up : "<< TDatabase::TimeDB->CURRENTTIME<<"   "<< ztn_1 <<  endl); }
        else
         { OutPut(setw(25)<<"T, Bott_Down : "<< TDatabase::TimeDB->CURRENTTIME<<"   "<< ztn_1 <<  endl); }     
      }
      
      
 if((m % 500 ) == 0   || m==1 )
  {
//    MovBoundVert[0][N_MovVert[0]-1]->GetCoords(x1, y1);
//    MovBoundVert[2][N_MovVert[2]-1]->GetCoords(x2, y2);
//    MovBoundVert[2][N_MovVert[2]-2]->GetCoords(x3, y3);   
//    MovBoundVert[2][N_MovVert[2]-3]->GetCoords(x4, y4);
//    
//    tx = x1-Rx;
//    sx = x2-Rx;
//    ty = y1-Ry;
//    sy = y2-Ry;
//    R_Theta[0] = acos( (tx*sx+ty*sy)/(sqrt(tx*tx+ty*ty)* sqrt(sx*sx+sy*sy)) )*(180/3.141592654);
// 
//    sx = x3-Rx;
//    sy = y3-Ry;
//    R_Theta[1] = acos( (tx*sx+ty*sy)/(sqrt(tx*tx+ty*ty)* sqrt(sx*sx+sy*sy)) )*(180./3.141592654);
// 
//    sx = ((x4))-Rx;
//    sy = ((y4))-Ry;
//    R_Theta[2] = acos( (tx*sx+ty*sy)/(sqrt(tx*tx+ty*ty)* sqrt(sx*sx+sy*sy)) )*(180./3.141592654); 
//    
//    if(!remeshed)
//     OutPut(setw(25)<<"T, wd,Ucl,RAng 1,2,3: " << TDatabase::TimeDB->CURRENTTIME<<"   "<< Rx-Lx
//                    <<"   "<< Params[2]<<"   "<<R_Theta[0]<<"   "<<R_Theta[1]<<"   "<<R_Theta[2]<<endl);    
// 
 
   MovBoundVert[0][0]->GetCoords(x1, y1);
   OutPut(setw(25)<<"T, Bottom : "<< TDatabase::TimeDB->CURRENTTIME<<"   "<< y1 <<  endl);  
   
   Get_KE(FEVectFuncts_All[0], Params);  
   OutPut(setw(25)<<"T, KE, Volume, Diff, Rel. Diff : "<< TDatabase::TimeDB->CURRENTTIME<<"   "<< Params[1] <<"   "<< Params[0]
                  <<"   "<< Params[0] - InitVolume<<"   "<< (Params[0] - InitVolume)/InitVolume << endl);
  }
 
 if((m % (int)(TDatabase::TimeDB->STEPS_PER_IMAGE) ) == 0   || m==1 )
  {

#ifdef __SURFACT__ 	 
       // interface surfactant
       MapSurfToDomain(IFaceFeFunct[0], FEFunctions_All[6], N_List[0], N_List[1]);
  
       GetSurfactMass(FEFunctions_All[6], IFaceFeFunct[0], N_List[0], N_List[1], Surf_Mass);

       Get_FeFunction2DMass(FEFunctions_All[5], Params);        
     
     
        OutPut( "Time, GammaMaxR, GammaMax " <<TDatabase::TimeDB->CURRENTTIME<<
              " " <<GammaXmaxVal[0]<< " "<<GammaXmaxVal[1]<< " "<<endl);   
      OutPut( "Time, Surfactant_Mass, Da*InterfaceSurfactant_Mass, InterfaceSurfactant_Conc " <<TDatabase::TimeDB->CURRENTTIME<<
              " " <<Params[0]<< " "<< TDatabase::ParamDB->REACTOR_P14*Surf_Mass[0]<< " "<<Surf_Mass[0]/Surf_Mass[1]<< " "<<endl);     
      OutPut( "Time, Surfactant Mass dif, InterfaceSurfactant Mass diff " <<TDatabase::TimeDB->CURRENTTIME<<
              " " <<(Params[0] - Initial_SurfactMass) << " "<<(TDatabase::ParamDB->REACTOR_P14*Surf_Mass[0] - Initial_IFaceSurfactMass) << " "<<endl);
      OutPut( "Time, Total Mass, Total Mass diff, RelativeTotal Mass diff " <<TDatabase::TimeDB->CURRENTTIME<<
          " " <<Params[0]+TDatabase::ParamDB->REACTOR_P14*Surf_Mass[0]<< " " << ((Params[0]+TDatabase::ParamDB->REACTOR_P14*Surf_Mass[0]) - 
            (Initial_IFaceSurfactMass + Initial_SurfactMass)) << " " << ((Params[0]+TDatabase::ParamDB->REACTOR_P14*Surf_Mass[0]) - 
            (Initial_IFaceSurfactMass + Initial_SurfactMass))
            /(Initial_IFaceSurfactMass + Initial_SurfactMass) <<endl);     
#endif     
     if(TDatabase::ParamDB->WRITE_VTK)
       {      
        os.seekp(std::ios::beg);
        if(img<10) os << "VTK/"<<VtkBaseName<<".0000"<<img<<".vtk" << ends;
         else if(img<100) os <<"VTK/"<< VtkBaseName<<".000"<<img<<".vtk" << ends;
          else if(img<1000) os <<"VTK/"<< VtkBaseName<<".00"<<img<<".vtk" << ends;
           else if(img<10000) os << "VTK/"<<VtkBaseName<<".0"<<img<<".vtk" << ends;
            else  os << "VTK/"<<VtkBaseName<<"."<<img<<".vtk" << ends;
        Output->WriteVtk(os.str().c_str());
        img++;
       }       
   }

    if( ((m % 200)== 0) || (m==1) )
     {
#ifdef __SURFACT__ 
      PrintSurfSurfactant(N_MovVert[0], MovBoundVert[0], FEFunctions_All[6], N_BData); 
#else
      PrintSurfSurfactant(N_MovVert[0], MovBoundVert[0], NULL, N_BData); 
#endif 
     } //if( ((m % 20)== 0) || (m==1) )
   
   if( m % 200== 0 )
    {
     OutPut(setw(25)<<TDatabase::TimeDB->CURRENTTIME<<" No. ReParam : " << N_ReParam <<endl);  
     OutPut(setw(25)<<TDatabase::TimeDB->CURRENTTIME<<" No. Remeshed : " << N_Remesh <<endl);  
    }

//======================================================================
     // adaptive time step based on no. of nonlinear iteration
//======================================================================    
//      if(ReduceTiemStep)
//       {
//        ReduceTiemStep = FALSE;
//        tau = TDatabase::TimeDB->TIMESTEPLENGTH;
//        if(tau/1.5>TDatabase::TimeDB->MIN_TIMESTEPLENGTH)
//        {
//         TDatabase::TimeDB->TIMESTEPLENGTH = tau/1.5; 
//        
//         OutPut("Time, Adaptive time step reduced to " << TDatabase::TimeDB->CURRENTTIME<< " " <<TDatabase::TimeDB->TIMESTEPLENGTH <<endl);
//        }
//       }
//      else if(IncrTimeStepCount>200)
//       {
//        tau = TDatabase::TimeDB->TIMESTEPLENGTH;
//        
//        if(tau*1.15<TDatabase::TimeDB->MAX_TIMESTEPLENGTH)
//         {
//          IncrTimeStepCount=0;  
//          TDatabase::TimeDB->TIMESTEPLENGTH = tau*1.25;   
//          OutPut("Time, Adaptive time step increased to " << TDatabase::TimeDB->CURRENTTIME<< " " << TDatabase::TimeDB->TIMESTEPLENGTH <<endl);       
//         }
//       }    
// 
//       
//     //======================================================================    
//     if(IncrTimeStepCount>200)
//      {
//       if(  TDatabase::TimeDB->MAX_TIMESTEPLENGTH<2.e-3   && TDatabase::TimeDB->CURRENTTIME>100)
//          TDatabase::TimeDB->MAX_TIMESTEPLENGTH=2.e-3;
//       else if ( TDatabase::TimeDB->MAX_TIMESTEPLENGTH<3.e-3  && TDatabase::TimeDB->CURRENTTIME>250)
//           TDatabase::TimeDB->MAX_TIMESTEPLENGTH=3.e-3;
// //       else if ( TDatabase::TimeDB->MAX_TIMESTEPLENGTH< 5.e-3  && TDatabase::TimeDB->CURRENTTIME>500)
// //           TDatabase::TimeDB->MAX_TIMESTEPLENGTH=5.e-3;
//       OutPut("Time, Adaptive time MAX_TIMESTEPLENGTH increased to " << TDatabase::TimeDB->CURRENTTIME<< " " << TDatabase::TimeDB->MAX_TIMESTEPLENGTH <<endl)
//       IncrTimeStepCount=0;
//       
// //       exit(0);
//      } 
//    //======================================================================     
      
      
 //======================================================================       
    
  } // while(TDatabase::TimeDB->CURRENTTIME< end_time)

  //======================================================================
  // end of time cycle
  //======================================================================   
     if(TDatabase::ParamDB->WRITE_VTK)
       { 
#ifdef __SURFACT__  
       // interface surfactant
       MapSurfToDomain(IFaceFeFunct[0], FEFunctions_All[6], N_List[0], N_List[1]);
#endif 
     
        os.seekp(std::ios::beg);
        if(img<10) os << "VTK/"<<VtkBaseName<<".0000"<<img<<".vtk" << ends;
         else if(img<100) os << "VTK/"<<VtkBaseName<<".000"<<img<<".vtk" << ends;
          else if(img<1000) os << "VTK/"<<VtkBaseName<<".00"<<img<<".vtk" << ends;
           else if(img<10000) os << "VTK/"<<VtkBaseName<<".0"<<img<<".vtk" << ends;
            else  os << "VTK/"<<VtkBaseName<<"."<<img<<".vtk" << ends;
        Output->WriteVtk(os.str().c_str());
        img++;
       }    

    // count total running time
    OutPut(TDatabase::TimeDB->CURRENTTIME << " total elapsed running time "<< total_time << endl);

   CloseFiles();  
   return 0;  
}

