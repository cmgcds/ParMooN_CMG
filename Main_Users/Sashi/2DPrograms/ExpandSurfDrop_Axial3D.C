// =======================================================================
//
// Purpose:     main program foraxisymmetric droplet oscillating with surfactant transport
//
// Author:     S.Ganesan
//             09.05.2007
// =======================================================================

#include <Domain.h>
#include <Database.h>
#include <FEDatabase2D.h>
#include <FESpace2D.h>
#include <SquareStructure1D.h>
#include <SquareStructure2D.h>
#include <SquareMatrix1D.h>
#include <Structure2D.h>
#include <FEFunction1D.h>
#include <AuxParam2D.h>
#include <Solver.h>
#include <DirectSolver.h>
#include <Assemble2D.h>
#include <Output2D.h>
#include <DiscreteForm2D.h>
#include <LinAlg.h>
#include <TNSE2D_ParamRout.h>

#include <Collection.h>
#include <NodalFunctional2D.h>
#include <string.h>
#include <sstream>
#include <MooNMD_Io.h>
#include <stdlib.h>
#include <malloc.h>

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
#include <TimeUtilities.h>

#include <MacroCell.h>
#include <BoundEdge.h>
#include <IsoBoundEdge.h>
#include <gridgen.h>
#include <Remesh2D.h>
#include <IsoInterfaceJoint.h>
#include <BdLine.h>
#include <BdCircle.h>
#include <GridCell.h>

#include <QuadAffin.h>
#include <QuadBilinear.h>
#include <QuadIsoparametric.h>
#include <TriaAffin.h>
#include <TriaIsoparametric.h>

#define AMG 0
#define GMG 1

// #include "../Examples/TNSE_2D/Droponsolid.h"
#include "../Examples/TNSE_2D/Dropinair_axial3D.h"

extern "C"
{
  void triangulate(char*, struct triangulateio*,
		   struct triangulateio*, struct triangulateio*);
}

/** convert current grid to vector-values FE function */
void GridTOVelocity(TFEVectFunct2D *velo, double t, double dt, double &max_r)
{
 int N_U, i;
 TFESpace2D *Velospace;
 double *x, *y, *ValuesU, u, v, newx, newy, L, r;


 Velospace = velo->GetFESpace2D();

 ValuesU = velo->GetValues();
 N_U = velo->GetLength();


 x = new double[N_U];
 y = new double[N_U];

 Velospace->GetDOFPosition(x, y);

//  cout << " dt " << dt << " time " << t <<endl;
// exit(0);
max_r = -1e10;
 for(i=0;i<N_U;i++)
  {
   r = sqrt(x[i]*x[i] + y[i]*y[i]);
   if(max_r<r) max_r =r;
 }

OutPut("Time, Radius of the Sphere : " << TDatabase::TimeDB->CURRENTTIME<< " " << max_r<<endl;)

 for(i=0;i<N_U;i++)
  {
//    newx = cos(100.*t)*x[i] - sin(100.*t)*y[i];
//    newy = sin(100.*t)*x[i] + cos(100.*t)*y[i];
   u = 0.;
   v = 0.;

   r = x[i]*x[i] + y[i]*y[i] ;
   if(sqrt(r)>0.25 )
    {
     L = pow(r,3./2.);
     u = x[i]/(L);
     v = y[i]/(L);
    }


   ValuesU[i] = u;
   ValuesU[N_U+i] =  v; 
    }


  delete []x;
  delete []y;
}




void  MoveGrid(TFEVectFunct2D *GridPos, TFEVectFunct2D *AuxGridPos, 
               TFESpace2D  *Surf_space, double t, double dt)
{
 int i, j, k, l, N_G, N_Cells, N_Vertices, *GridBeginIndex, *GridGlobalNumbers, *DOF;
 int N_Edges;
 double *X, *Y, *NewX, *NewY, u, L, r, u_r, u_z, x, y, max_r=-1e10;
 TCollection *Coll;
 TBaseCell *cell;
 TFESpace2D  *GridSpace;
 TJoint *joint;
 TIsoBoundEdge *isojoint;
 TVertex **Vertices;
 int IIso, m, *VeloGlobalNumbers, *VeloBeginIndex;
 FE2D FEId;
 TFEDesc2D *FEDesc;
 int *VeloDOF, *JointDOF;

 GridSpace = GridPos->GetFESpace2D();
 Coll = GridSpace->GetCollection();
 N_Cells = Coll->GetN_Cells();
 GridBeginIndex = GridSpace->GetBeginIndex();
 GridGlobalNumbers = GridSpace->GetGlobalNumbers();

 GridPos->GridToData();
 AuxGridPos->GridToData();


 N_G=GridPos->GetLength();
 X = GridPos->GetValues();
 Y = X + N_G;

 NewX = AuxGridPos->GetValues();
 NewY = NewX + N_G;


 VeloGlobalNumbers = Surf_space->GetBeginIndex();
 VeloBeginIndex = Surf_space->GetGlobalNumbers();


//  for(i=0;i<N_G;i++)
//   {
//    r = X[i]*X[i] + Y[i]*Y[i] ;
//    if(max_r<r) max_r =r;
//  }


//  u = 10.;
 for(i=0;i<N_G;i++)
  {
   r = X[i]*X[i] + Y[i]*Y[i] ;
   if(sqrt(r)>0.25 )
    {
     L = pow(r,3./2.);
     u_r = X[i]/(L);
     u_z = Y[i]/(L);

   NewX[i] = X[i] + dt*u_r;
   NewY[i] = Y[i] + dt*u_z;

//   cout<<  " X[i] " <<X[i] << " Y[i] "<< Y[i] << " NewX[i] "  <<NewX[i] << " NewY[i] " <<NewY[i]<<endl;

   } //    if(sqrt(r)>0.5 
  } // for(i=0;i<N_G;
//   IIso = 0;
  for(i=0;i<N_Cells;i++)
   {
    cell = Coll->GetCell(i);
    DOF = GridGlobalNumbers + GridBeginIndex[i];
    N_Vertices = cell->GetN_Vertices();
    switch(N_Vertices)
    {
      case 3:
        for(j=0;j<N_Vertices;j++)
        {
          k = DOF[j];
          cell->GetVertex(j)->SetCoords(NewX[k], NewY[k]);
        }
      break;

      case 4:
        k = DOF[0];
        cell->GetVertex(0)->SetCoords(NewX[k], NewY[k]);
        k = DOF[1];
        cell->GetVertex(1)->SetCoords(NewX[k], NewY[k]);
        k = DOF[3];
        cell->GetVertex(2)->SetCoords(NewX[k], NewY[k]);
        k = DOF[2];
        cell->GetVertex(3)->SetCoords(NewX[k], NewY[k]);
      break;
    } // endswitch



    N_Edges = cell->GetN_Edges();
    for(j=0;j<N_Edges;j++)
     {
      joint = cell->GetJoint(j);
      if(joint->GetType() == IsoBoundEdge)
      {
        isojoint = (TIsoBoundEdge *)joint;
        k = isojoint->GetN_Vertices();
        Vertices = isojoint->GetVertices();
        FEId = Surf_space->GetFE2D(i, cell);
        FEDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(FEId);
        m = FEDesc->GetN_JointDOF();
        if(m == k+2)
         {
          JointDOF = FEDesc->GetJointDOF(j);
          VeloDOF =  VeloGlobalNumbers+VeloBeginIndex[i];
          for(l=0;l<k;l++)
          {
            m = VeloDOF[JointDOF[l+1]];
            Vertices[l]->GetCoords(x, y);
            r = x*x + y*y;
            if(sqrt(r)>0.25)
             {
              L = pow(r,3./2.);
              u_r = x/(L);
              u_z = y/(L);
              x += dt*u_r;
              y += dt*u_z;

             Vertices[l]->SetCoords(x, y);
            } //  if(sqrt(r)
//             IIso++;
          } // endfor l

         } // if(m == k+2)
        else
        {
          // approximation order of isoparametric boundary and velocity
          // element must be the same
          Error("No match in isoparametric case" << endl);
          exit(-1);
        }
       } // endif
     } // endfor j
   }
// cout<< " r " << r <<endl;
//  AuxGridPos->DataToGrid();
}


void PrintSurfSurfactant(TBaseCell **Cell, TVertex **Vertex, int *EdgeNo, int N, 
                         TFEFunction2D *Surfact, int &N_BData)
{
 int i, j, k, Cell_No, IJoint, N_DOF_Local;
 double  T_val[3], x1, y1, x2, y2, ArcLength=0.;
 double char_L = TDatabase::ParamDB->CHAR_L0;
 char *VtkBaseName;
 VtkBaseName = TDatabase::ParamDB->VTKBASENAME;

 std::ostringstream os;
 os << " ";

 os.seekp(std::ios::beg);
//       os << "surfact"<< i << ".dat" << ends;
  if(N_BData<10) os << VtkBaseName<<"_0000"<<N_BData<<".data" << ends;
  else if(N_BData<100) os <<VtkBaseName<<"_000"<<N_BData<<".data" << ends;
  else if(N_BData<1000) os <<VtkBaseName<<"_00"<<N_BData<<".data" << ends;
  else if(N_BData<10000) os <<VtkBaseName<<"_0"<<N_BData<<".data" << ends;
  else  os <<VtkBaseName<<"_"<<N_BData<<".data" << ends;

  std::ofstream dat(os.str().c_str());

  if (!dat)
   {
    cerr << "cannot open file for output" << endl;
    exit(0);
   }
  dat << "%% Surfactant data created by MooNMD" << endl;
  dat << "%% Current Reference Time :" << TDatabase::TimeDB->CURRENTTIME << endl;
  dat << "%% x, y, Orig_ArcLength,  Surfact" << endl;

  Vertex[0]->GetCoords(x1, y1);
  for(i=0;i<N;i++)
   {
    Vertex[i]->GetCoords(x2, y2);
//    cout<<i<< " Angle of free Vertices "<<(180./Pi)*atan2(y,x)<<endl;

    ArcLength += sqrt( (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) );

    Surfact->FindGradient( x1,  y1, T_val);
    dat << x1<< " " << y1 << " "<< ArcLength*char_L << "  " << T_val[0] <<endl;

    x1=x2; y1=y2;
   }

      dat.close();
      cout << endl;
      cout << "Surfactant data wrote into file " << endl;
 N_BData++;

}



void GetSurfactMass(TFEFunction2D *fefunction, TFEFunction1D *fefunct_low,
                   int *Cell_array, int *Joint_array, double *errors, double max_r)
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
  double  X_B[MaxN_QuadPoints_2D], Y_B[MaxN_QuadPoints_2D], Exact_Surf[5], v, Mult;
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
          ((TTriaIsoparametric *)F_K)->SetApproximationOrder(ORDER);
      break;

      case BFUnitTriangle:
          RefTrans = TriaIsoparametric;
          F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
          ((TTriaIsoparametric *)F_K)->SetCell(Me);
          ((TTriaIsoparametric *)F_K)->GetOrigBoundFromRef(IJoint,
	                            N_LinePoints, zeta, X_B, Y_B);
          ((TTriaIsoparametric *)F_K)->SetApproximationOrder(ORDER);

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
              ((TTriaIsoparametric *)F_K)->SetCell(Me); 
              ((TTriaIsoparametric *)F_K)->SetApproximationOrder(ORDER);
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

      errors[0] +=U*Mult;
      errors[1] +=Mult;
//       errors[2] +=Mult*(Exact_Surf[0]-U)*(Exact_Surf[0]-U);
     } //  for(k=0;k<N_LinePoints;k++)

//      if(h_K_min>h_K)
//          h_K_min = h_K;
//      if(h_K_max<h_K)
//          h_K_max = h_K;

  } // for(i=0;i<N
//     OutPut("h_K_min and h_K_max of free surface: "<< h_K_min << " " << h_K_max<<endl; );

   OutPut( "Time, Surfactant Mass " <<TDatabase::TimeDB->CURRENTTIME<< " "<<2.*Pi*errors[0]<< " "<<endl);
   OutPut( "Time, Surface area " <<TDatabase::TimeDB->CURRENTTIME<< " "<<2.*Pi*errors[1]<< " "<<endl);
   OutPut( "Time, Surfactant Concentration " <<TDatabase::TimeDB->CURRENTTIME<< " "<< 2.*Pi*errors[0] / (4.*Pi*max_r*max_r)<< " "<<endl);   
// exit(0);
} // void GetSurfactMass(


void GetSurfErrors(TFEFunction2D *fefunction, TFEFunction1D *fefunct_low, double *errors,
                   int *Cell_array, int *Joint_array)
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

  errors[0] =0.; // L2-norm
  errors[1] =0.; // H1-seminorm

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
          ((TTriaIsoparametric *)F_K)->SetApproximationOrder(ORDER);
      break;

      case BFUnitTriangle:
          RefTrans = TriaIsoparametric;
          F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
          ((TTriaIsoparametric *)F_K)->SetCell(Me);
          ((TTriaIsoparametric *)F_K)->GetOrigBoundFromRef(IJoint,
	                            N_LinePoints, zeta, X_B, Y_B);
          ((TTriaIsoparametric *)F_K)->SetApproximationOrder(ORDER);

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
              ((TTriaIsoparametric *)F_K)->SetCell(Me); 
              ((TTriaIsoparametric *)F_K)->SetApproximationOrder(ORDER);
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

//       get solution and its gradients
     U=0.; ux=0.;  uy=0.; 
     for(l=0;l<N_BaseFunct_low;l++)
      {
       local_dof   = JointDOF[l];

       test00 = uorig[local_dof];
//        test10 = uxorig[local_dof];
//        test01 = uyorig[local_dof];

//        ngrad = test10*n0 + test01*n1;
//        d1 = test10 - ngrad*n0;
//        d2 = test01 - ngrad*n1;

       m = DOF[local_dof];

       U  += u[m]*test00;
//        ux += u[m]*d1;
//        uy += u[m]*d2;
      }

      ExactS(X_P, Y_P, Exact_Surf);

      Mult = LineWeights[k]*normn*r_axial; // since the line integral is -1 to +1
      h_K +=normn;
      errors[0] +=Mult*(Exact_Surf[0]-U)*(Exact_Surf[0]-U);
     } //  for(k=0;k<N_LinePoints;k++)

     if(h_K_min>h_K)
         h_K_min = h_K;
     if(h_K_max<h_K)
         h_K_max = h_K;

  } // for(i=0;i<N
    OutPut("h_K_min and h_K_max of free surface: "<< h_K_min << " " << h_K_max<<endl; );
   errors[0] = sqrt(errors[0]);
//    cout << "errors[0] " << errors[0]<< " "<<endl;

// exit(0);
}


void  AssembleSurf1D(int n_fespaces, TFESpace2D **fespaces, TFEFunction2D **fefunctions, 
                     int N_FESpaces_low, TFESpace1D **fespaces_low, int N_SquareMatrices,
                     TSquareMatrix1D **sqmatrices_low, int N_Rhs, double **RHSs, 
                     TFESpace1D **ferhs_low, int *Cell_array, int *Joint_array)
{
  int i, j, k, l, m, n, N_Cells_low, N, N_LocalUsedElements, local_i, local_j;
  int N_BaseFunct, N_BaseFunct_low,  N_Points, N_JointDOF, begin, end, *N_BaseFuncts;
  TCollection *Coll, *Coll_low;
  TBaseCell *Me, *Me_low;
  FE2D FEId;
  FE1D FEId_low;
  TFE1D *Element;
  TFE2D *ele;
  double x0, y0, x1, y1, t0, t1, n0, n1, normn;
  int LocN_BF[N_BaseFuncts2D], N_LinePoints, *KCol, *RowPtr, *JointDOF;
  BaseFunct2D LocBF[N_BaseFuncts2D];
  BaseFunct2D *BaseFuncts;
  double AbsDetjk[MaxN_QuadPoints_2D], Mult;
  double *weights, *xi, *eta;
  double **uref, **uxiref, **uetaref;
  double *LineWeights, *zeta, *ValuesA, *ValuesM;
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];
  double uorig[MaxN_BaseFunctions2D], uxorig[MaxN_BaseFunctions2D], uyorig[MaxN_BaseFunctions2D];
  boolean *SecondDer;
  BF2DRefElements RefElement;
  RefTrans2D RefTrans;
  TRefTrans2D *F_K;
  QuadFormula1D LineQuadFormula;
  TQuadFormula1D *qf1;
  int *BeginIndex_low, *GlobalNumbers_low, *DOF, *DOF_LOW, TestDOF, AnsatzDOF, IJoint ;
  int *BeginIndex, *GlobalNumbers,ORDER;
  TFEDesc2D *FeDesc;
  TFEDesc1D *FeDesc_low;
  double val, theta, ngrad_ansatz, ngrad_test, TangDivU;
  double  X_B[100], Y_B[100], r_axial, d1, d2, e1, e2;
  double LocMatrixA[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double LocMatrixM[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double ansatz00, ansatz10, ansatz01 ;
  double test00, test10, test01, *u1, *u2, u1x, u2x, u1y, u2y, U1;

  double c0, Re = TDatabase::ParamDB->RE_NR, r2;
  double Pr = TDatabase::ParamDB->Pr_NR;
  
  if(Re==0 || Pr==0)
    c0=0.;
  else
    c0 =  1.0/(Re*Pr);  // Peclet number

  SecondDer = new boolean[n_fespaces];
// ########################################################################
// store information in local arrays
// ########################################################################
  BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
  N_BaseFuncts = TFEDatabase2D::GetN_BaseFunctFromFE2D();


  Coll = fespaces[0]->GetCollection(); // all spaces use same Coll
  BeginIndex = fespaces[0]->GetBeginIndex();
  GlobalNumbers = fespaces[0]->GetGlobalNumbers();
  u1 = fefunctions[0]->GetValues();
  u2 = fefunctions[1]->GetValues();
//   N_Cells = Coll->GetN_Cells();

  Coll_low = fespaces_low[0]->GetCollection(); // all low spaces use same Coll
  N_Cells_low = Coll_low->GetN_Cells();
  BeginIndex_low =  fespaces_low[0]->GetBeginIndex();
  GlobalNumbers_low =  fespaces_low[0]->GetGlobalNumbers();

  RowPtr = sqmatrices_low[0]->GetRowPtr();
  KCol = sqmatrices_low[0]->GetKCol();

  ValuesA = sqmatrices_low[0]->GetEntries();
  ValuesM = sqmatrices_low[1]->GetEntries();


  N_LocalUsedElements = n_fespaces;
  for(j=0;j<n_fespaces;j++)
    SecondDer[j]=FALSE;
// ########################################################################
// loop over all low space cells
// ########################################################################
  for(i=0;i<N_Cells_low;i++)
  {
    N = Cell_array[i];
    Me = Coll->GetCell(N);
    IJoint = Joint_array[i];

    FEId = fespaces[0]->GetFE2D(N, Me);  // velocity space in the entire domain
    ele = TFEDatabase2D::GetFE2D(FEId);
    RefElement = TFEDatabase2D::GetRefElementFromFE2D(FEId);
    FeDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(FEId);
    N_JointDOF = FeDesc->GetN_JointDOF();
    JointDOF = FeDesc->GetJointDOF(IJoint);
    N_BaseFunct = FeDesc->GetN_DOF();
    DOF = GlobalNumbers + BeginIndex[N];

    DOF_LOW = GlobalNumbers_low + BeginIndex_low[i];
    Me_low = Coll_low->GetCell(i);
    FEId_low = fespaces_low[0]->GetFE1D(i, Me_low);
    Element = TFEDatabase2D::GetFE1D(FEId_low);
    N_BaseFunct_low = Element->GetN_DOF();
    if(N_JointDOF != N_BaseFunct_low )
     {
      cout<< " N_JointDOF != N_BaseFunct_low " <<endl;
      exit(0);
    }

    memset(LocMatrixA, 0, N_BaseFunct_low*N_BaseFunct_low*SizeOfDouble);
    memset(LocMatrixM, 0, N_BaseFunct_low*N_BaseFunct_low*SizeOfDouble);


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
          ((TTriaIsoparametric *)F_K)->SetApproximationOrder(ORDER);
        break;

        case BFUnitTriangle:
          RefTrans = TriaIsoparametric;
          F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
          ((TTriaIsoparametric *)F_K)->SetCell(Me);
          ((TTriaIsoparametric *)F_K)->SetApproximationOrder(ORDER);
        break;
      } // endswitch


      TFEDatabase2D::GetBaseFunct2DFromFE2D(FEId)->MakeRefElementData(LineQuadFormula);

//     cout<< 2*l <<" LineQuadFormula " << N_LinePoints <<endl;
//       for(k=0;k<N_LinePoints;k++)
//        cout << "zeta: " << zeta[k] << " LineWeights[k] " << LineWeights[k]<<endl;
// exit(0);
//            cout << endl;
//            r2 = 0;
      for(k=0;k<N_LinePoints;k++)
        {
          F_K->GetTangent(IJoint, zeta[k], t0, t1);
          normn = sqrt(t0*t0+t1*t1);
          n0 =  t1/normn;
          n1 = -t0/normn;
//            cout << endl;
//     cout << "zeta: " << zeta[k] << " LineWeights[k] " << LineWeights[k]<<endl;
//           cout << "k= " << k << "  tangent: " << t0 << " " << t1 << endl;
          // cout << "length: " << sqrt(t0*t0+t1*t1) << endl;
          uref = TFEDatabase2D::GetJointValues2D(BaseFuncts[FEId],
                        LineQuadFormula, IJoint);
          uxiref = TFEDatabase2D::GetJointDerivatives2D(BaseFuncts[FEId],
                        LineQuadFormula, IJoint, D10);
          uetaref = TFEDatabase2D::GetJointDerivatives2D(BaseFuncts[FEId],
                        LineQuadFormula, IJoint, D01);

          switch(RefElement)
          {
            case BFUnitSquare:
              ((TQuadIsoparametric *)F_K)->GetOrigValues(IJoint,  zeta[k],
                        N_BaseFunct, uref[k], uxiref[k], uetaref[k],
                        uorig, uxorig, uyorig);
          ((TQuadIsoparametric *)F_K)->GetOrigBoundFromRef(IJoint, N_LinePoints, zeta, X_B, Y_B);
            break;

            case BFUnitTriangle:
              ((TTriaIsoparametric *)F_K)->GetOrigValues(IJoint,  zeta[k],
                        N_BaseFunct, uref[k], uxiref[k], uetaref[k],
                        uorig, uxorig, uyorig);
          ((TTriaIsoparametric *)F_K)->GetOrigBoundFromRef(IJoint, N_LinePoints, zeta, X_B, Y_B);
            break;
          } // endswitch

          r_axial = fabs(X_B[k]); // r value in the axial symmetric integral
//           cout << " x " << r_axial<< " y " << Y_B[k]<< endl;

//       get velocity gradients
          U1 = 0.;  u1x=0.; u2x=0.; u1y=0.; u2y=0.;
//           for(l=0;l<N_BaseFunct;l++)
//             {
//              m = DOF[l];
//              U1  += u1[m]*uorig[l];
//              u1x += u1[m]*uxorig[l];
//              u1y += u1[m]*uyorig[l];
//              u2x += u2[m]*uxorig[l];
//              u2y += u2[m]*uyorig[l];
//             }

          for(l=0;l<N_BaseFunct_low;l++)
            {
             local_j = JointDOF[l];
             m = DOF[local_j];
             U1  += u1[m]*uorig[local_j];
             u1x += u1[m]*uxorig[local_j];
             u1y += u1[m]*uyorig[local_j];
             u2x += u2[m]*uxorig[local_j];
             u2y += u2[m]*uyorig[local_j];
            }


          TangDivU =  u1x - (u1x*n0 + u1y*n1)*n0  + U1/r_axial
                    + u2y - (u2x*n0 + u2y*n1)*n1;

// // if div U = 0
//           TangDivU =  - (u1x*n0 + u1y*n1)*n0  - (u2x*n0 + u2y*n1)*n1;

//           cout <<" u1x:  " << u1x <<" u1y:  " << u1y <<endl;
//           cout <<" u2x:  " << u2x <<" u2y:  " << u2y <<endl;
//           cout <<" TangDivU:  " << TangDivU <<endl;

          Mult = sqrt(t0*t0+t1*t1)*(LineWeights[k]);
          for(l=0;l<N_BaseFunct_low;l++)
           {
            local_j   = JointDOF[l];

            test00  = uorig[local_j];
            test10  = uxorig[local_j];
            test01  = uyorig[local_j];

            ngrad_test= n0*test10 + n1*test01;
            d1 = test10 - ngrad_test*n0;
            d2 = test01 - ngrad_test*n1;

//            rhs

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
              val  = c0*(d1*e1 + d2*e2);
              val +=TangDivU*test00*ansatz00;
              val *= (Mult*r_axial);
              LocMatrixA[l*N_BaseFunct_low+m] += val;

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
// exit(0);
}


// ====================================================================
// modify axisymmetriv drop fluid matrices due to integrals 
// on free surface with variable surface tension w.r.t Surfact
// ====================================================================

void FreeSurf_axial3D_Surfact(TSquareMatrix2D *A11, TSquareMatrix2D *A22,
                          double *rhs1, double *rhs2,
                          BoundCondFunct2D *BoundaryCondition,
                          double dt, TFEFunction2D *Surfact)
{
  int i, j, k, l, DOF_R, DOF_L, TDOF, TDOF_R, m, TN_DOF_Local;
  TBaseCell *cell;
  TCollection *Coll;
  int N_Cells, N_Vertices, N_Edges, Semi_implicit=0;
  TJoint *joint;
  TIsoBoundEdge *isoboundedge;
  TBoundComp *BoundComp;
  int comp, N_U, test_L=1, test_R=1;
  double t0, t1, n0, n1, normn, line_wgt;
  BoundCond Cond0, Cond1;
  int JointNumbers[MAXN_JOINTS], IJoint, N_IsoJoints;
  FE2D FEId, TFEId;
  TFE2D *ele;
  RefTrans2D RefTrans;
  TRefTrans2D *F_K;
  TFESpace2D *fespace, *thermalspace;
  BF2DRefElements RefElement;
  QuadFormula1D LineQuadFormula;
  TQuadFormula1D *qf1;
  int N_LinePoints;
  double *LineWeights, *zeta;
  double x0, y0, x1, y1,tx,ty,mod_t, x, y;
  int N_BaseFunct, *N_BaseFuncts;
  double **uref, **uxiref, **uetaref;
  double uorig[MaxN_BaseFunctions2D], uxorig[MaxN_BaseFunctions2D];
  double uyorig[MaxN_BaseFunctions2D];
  BaseFunct2D *BaseFuncts;
  double r2, r;
  int *KCol, *RowPtr, *JointDOF, N_DOF, *TJointDOF;
  double *ValuesA11, *ValuesA12, *ValuesA21, *ValuesA22;
  int *BeginIndex, *GlobalNumbers, *DOF, TestDOF, AnsatzDOF;
  int index1, index2;
  double val, theta, factor1, factor2, angle, T_val[3], T_Marangoni, *Tvalues;
  int count=0, count1=0, count2=0;
  double  X_B[100], Y_B[100], T[100], r_axial, T_DOF[10], d1, d2, e1, e2, ngrad_test, ngrad_ansatz;
  int *TGlobalNumbers, *TBeginIndex;

  TFEDesc2D *FeDesc, *TFeDesc;

  BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
  N_BaseFuncts = TFEDatabase2D::GetN_BaseFunctFromFE2D();

  fespace = A11->GetFESpace();
  Coll = fespace->GetCollection();
  N_Cells = Coll->GetN_Cells();

  thermalspace = Surfact->GetFESpace2D();
  Tvalues=Surfact->GetValues();
  TGlobalNumbers = thermalspace->GetGlobalNumbers();
  TBeginIndex = thermalspace->GetBeginIndex();

  BeginIndex = fespace->GetBeginIndex();
  GlobalNumbers = fespace->GetGlobalNumbers();

  RowPtr = A11->GetRowPtr();
  KCol = A11->GetKCol();

  ValuesA11 = A11->GetEntries();
  ValuesA22 = A22->GetEntries();

//   double Re = TDatabase::ParamDB->RE_NR;
  double We = TDatabase::ParamDB->WB_NR;
  double E = TDatabase::ParamDB->P13;  // surfactant elasticity E
  int EOS = int(TDatabase::ParamDB->P14); //Equation of state, 0 linear, 1 non-linear
  double D = TDatabase::ParamDB->P15 / TDatabase::ParamDB->P10 ; //\Gamma_0/Gamma_\infty

  double ST_Angle_R = TDatabase::ParamDB->P2;
  double ST_Angle_L = TDatabase::ParamDB->P3;
    ST_Angle_R = (3.141592654/180)*ST_Angle_R;
    ST_Angle_L = (3.141592654/180)*ST_Angle_L;


//  cout << "Cell  has free surface." << endl;
 for(i=0;i<N_Cells;i++)
  {
//      cout << endl << "CELL number: " << i << endl;
    cell = Coll->GetCell(i);
    N_Edges = cell->GetN_Edges();
    IJoint = 0;
    for(j=0;j<N_Edges;j++)
     {
      joint = cell->GetJoint(j);
      if(joint->GetType() == IsoBoundEdge)
      {
        isoboundedge = (TIsoBoundEdge *)joint;
        BoundComp = isoboundedge->GetBoundComp();
        isoboundedge->GetParameters(t0, t1);
        comp=BoundComp->GetID();
        BoundaryCondition(comp, t0, Cond0);
        BoundaryCondition(comp, t1, Cond1);

        if(Cond0 == FREESURF)
        {
          JointNumbers[IJoint] = j;
          IJoint++;
        }
      } // endif
     } // endfor j


    N_IsoJoints = IJoint;
    if(N_IsoJoints > 0)
    {

      FEId = fespace->GetFE2D(i, cell);

      for(j=0;j<N_IsoJoints;j++)
      {
//      cout << "Cell " << i << " has free surface." << endl;
        IJoint = JointNumbers[j];

        cell->GetVertex(IJoint)->GetCoords(x0, y0);
        cell->GetVertex((IJoint+1) % N_Edges)->GetCoords(x1, y1);
     // entries for wetting DOF
      if(y0==0) // right wett point edge (bottom)
       {
        FeDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(FEId);
        JointDOF = FeDesc->GetJointDOF(IJoint);
        N_DOF = FeDesc->GetN_JointDOF();
        for(m=0;m<N_DOF;m++)
         {
          DOF_R =  GlobalNumbers[BeginIndex[i]+JointDOF[m]];
          fespace->GetDOFPosition(DOF_R, x, y);
          if(y==0) // right wett point
          {
//             cout<< "  x= "<< x <<"  y= "<< y <<endl;
//             RefTrans = TriaIsoparametric;
//             Surfact->FindGradientLocal(cell, i,  x, y, T_val);
//             Surfact->FindGradient( x,  y, T_val);
//             if(fabs(T_val[0])>0.5 || T_val[0]<-0.5 )
             TFEId = thermalspace->GetFE2D(i, cell);
             TFeDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(TFEId);
             TJointDOF = TFeDesc->GetJointDOF(IJoint);
//              TN_DOF_Local = TFeDesc->GetN_JointDOF();
             TDOF_R =  TGlobalNumbers[TBeginIndex[i]+TJointDOF[0]];
             T_val[0]=Tvalues[TDOF_R];
             if(fabs(T[k])>0.5 || T[k]<-0.5 )
              OutPut("x : "<<x<< " y: " << y <<"  Surfactant exceeds the reference value, T= " <<T_val[0]<<endl);
            r_axial = x;       // r value in the axial symmetric integral

            if(EOS==0)
             {
//             rhs1[DOF_R] +=  (1. - E* T_val[0] )*r_axial*((cos(ST_Angle_R))/We);
              rhs1[DOF_R] +=  (1. + E*(D - T_val[0]) )*r_axial*((cos(ST_Angle_R))/We);
             }
            else
             {
              rhs1[DOF_R] +=  (1. + E*log(1. - T_val[0]) )*r_axial*((cos(ST_Angle_R))/We);
             }
           break;
          }
        }
       }

      DOF = GlobalNumbers + BeginIndex[i];
      N_BaseFunct = N_BaseFuncts[FEId];
      ele = TFEDatabase2D::GetFE2D(FEId);
      RefElement = TFEDatabase2D::GetRefElementFromFE2D(FEId);

      l = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
      LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(2*l);
      qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
      qf1->GetFormulaData(N_LinePoints, LineWeights, zeta);
      TFEDatabase2D::GetBaseFunct2DFromFE2D(FEId)
                  ->MakeRefElementData(LineQuadFormula);

      switch(RefElement)
      {
        case BFUnitSquare:
          RefTrans = QuadIsoparametric;
          F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
          ((TQuadIsoparametric *)F_K)->SetCell(cell);
          ((TQuadIsoparametric *)F_K)->GetOrigBoundFromRef(IJoint,
	                            N_LinePoints, zeta, X_B, Y_B);
        break;

        case BFUnitTriangle:
          RefTrans = TriaIsoparametric;
          F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
          ((TTriaIsoparametric *)F_K)->SetCell(cell);
          ((TTriaIsoparametric *)F_K)->GetOrigBoundFromRef(IJoint,
	                            N_LinePoints, zeta, X_B, Y_B);

        break;
      } // endswitch


      uref = TFEDatabase2D::GetJointValues2D(BaseFuncts[FEId],
                        LineQuadFormula, IJoint);
      uxiref = TFEDatabase2D::GetJointDerivatives2D(BaseFuncts[FEId],
                        LineQuadFormula, IJoint, D10);
      uetaref = TFEDatabase2D::GetJointDerivatives2D(BaseFuncts[FEId],
                        LineQuadFormula, IJoint, D01);


     TFEId = thermalspace->GetFE2D(i, cell);
     TFeDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(TFEId);
     TJointDOF = TFeDesc->GetJointDOF(IJoint);
     TN_DOF_Local = TFeDesc->GetN_JointDOF();
      for(k=0;k<TN_DOF_Local;k++)
        {
          TDOF =  TGlobalNumbers[TBeginIndex[i]+TJointDOF[k]];
           T_DOF[k]  = Tvalues[TDOF];
//          cout<<k<< " T " << Tvalues[TDOF] <<endl;
        }

    if(TN_DOF_Local==2) // linear element
       for(k=0;k<N_LinePoints;k++)
          {
            T[k]  = zeta[k]*T_DOF[0] + (1.0 - zeta[k])*T_DOF[1];
//              cout << "zeta[k]  " << zeta[k] << " T[k]  " << T[k] <<endl;
          }
    else if(TN_DOF_Local==3) // quadratic element
       for(k=0;k<N_LinePoints;k++)
         {
          if(zeta[k]<0.)        //-1< zeta[k] < +1
           {
            t0 = (1.0 + zeta[k]);
            T[k]  = t0*T_DOF[0] + (1.0 - t0)*T_DOF[1];
//             cout << "zeta[k]  " << zeta[k] << " to " << t0<< " T[k]  " << T[k] <<endl;
           }
         else 
           {
            t0 =  zeta[k];
            T[k]  = t0*T_DOF[1] + (1.0 - t0)*T_DOF[2];
//             cout << "zeta[k]  " << zeta[k] << " to " << t0<< " T[k]  " << T[k] <<endl;
           }

         }
    else
       {
         OutPut( "surfactant space should be linear or quadratic !! Check Freesurfint " <<endl);
         exit(0);
       }

       for(k=0;k<N_LinePoints;k++)
        {
         switch(RefElement)
          {
            case BFUnitSquare:
              ((TQuadIsoparametric *)F_K)->SetCell(cell);
              ((TQuadIsoparametric *)F_K)->GetOrigValues(IJoint,  zeta[k],
                        N_BaseFunct, uref[k], uxiref[k], uetaref[k],
                        uorig, uxorig, uyorig);
            break;

            case BFUnitTriangle:
              ((TTriaIsoparametric *)F_K)->SetCell(cell); 
              ((TTriaIsoparametric *)F_K)->GetOrigValues(IJoint,  zeta[k],
                        N_BaseFunct, uref[k], uxiref[k], uetaref[k],
                        uorig, uxorig, uyorig);
            break;
          } // endswitch

          // modify matrices

         F_K->GetTangent(IJoint, zeta[k], t0, t1);  // old line
         r_axial = fabs(X_B[k]);   // r value in the axial symmetric integral
         normn = sqrt(t0*t0+t1*t1);

         t0 = t0/normn;
         t1 = t1/normn;

         n0 =  t1;
         n1 = -t0;

          // affine transformation only !!!
          Surfact->FindGradientLocal(cell, i, X_B[k], Y_B[k], T_val); 
//        cout<< "  x= "<< X_B[k] <<"  y= "<< Y_B[k] <<endl;
        if(fabs(T[k])>1. || T[k]<0. )
            OutPut(i<< "x : "<<X_B[k]<< " y: " << Y_B[k] <<"  Surfactant exceeds the reference value, T= " <<T[k]<<endl);
          // Multiply with time step dt in the main program not here
          r = normn/We;

          if(EOS==0)
            r *=(1. + E*(D - T_val[0]) );
          else
            r *=(1. + E*log(1. - T_val[0]) );


//        (c_1\sigma_sa) \tau\cdot \grad T,       norm for integral weights
          if(EOS==0)
            T_Marangoni = normn*E*(t0*T_val[1] + t1*T_val[2])/We;
          else
            T_Marangoni = normn*E*( t0*T_val[1] + t1*T_val[2]) / ( We*(1. - T_val[0]) );


          for(l=0;l<N_BaseFunct;l++)
          {
           TestDOF = DOF[l];

           // updating rhs
            ngrad_test= n0*uxorig[l] + n1*uyorig[l];
            d1 = uxorig[l] - ngrad_test*n0;
            d2 = uyorig[l] - ngrad_test*n1;

//  rhs1
            val = r_axial*( (1-n0*n0)*d1 - n0*n1*d2 );
            val +=uorig[l];
            val *= LineWeights[k]*r;
            rhs1[TestDOF] -= val;
//     Marangoni convection
            val = r_axial*t0*uorig[l];
            val *= LineWeights[k]*T_Marangoni;
            rhs1[TestDOF] -= val;

// rhs2
            val =  r_axial*( -n1*n0*d1 + (1-n1*n1)*d2 );
            val *= LineWeights[k]*r;
            rhs2[TestDOF] -= val;
//     Marangoni convection
            val = r_axial*t1*uorig[l];
            val *= LineWeights[k]*T_Marangoni;
            rhs2[TestDOF] -= val;

            index2 = RowPtr[TestDOF+1];

            for(m=0;m<N_BaseFunct;m++)
            {
              AnsatzDOF = DOF[m];
              // cout << AnsatzDOF << " -- " << TestDOF << endl;
              index1 = RowPtr[TestDOF];
              if(index1+1 == index2) continue;
              while(KCol[index1] != AnsatzDOF) index1++;

              ngrad_ansatz= n0*uxorig[m] + n1*uyorig[m];
              e1 = uxorig[m] - ngrad_ansatz*n0;
              e2 = uyorig[m] - ngrad_ansatz*n1;

              val =d1*e1 + d2*e2 + (uorig[l]*uorig[m]/(r_axial*r_axial));
              val *= dt*LineWeights[k]*r*r_axial;
              // cout << "A11: " << TestDOF << " ";
              // cout << AnsatzDOF << " " << val << endl;
//               ValuesA11[index1] += val;

              val = d1*e1 + d2*e2;
              val *= dt*LineWeights[k]*r*r_axial;

              // cout << "A22: " << TestDOF << " ";
              // cout << AnsatzDOF << " " << val << endl;
//               ValuesA22[index1] += val;

            } // endfor m
          } // endfor l
        } // endfor k
      } // endfor j

    } // end (N_IsoJoints > 0)
    else
    {
      // cout << "Cell " << i << " has NO free surface." << endl;
    }

  } // endfor i
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
//     Phase_No = cell->GetPhase_No();
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
//     TFEDatabase2D::SetCellForRefTrans(cell, RefTrans);
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
        ((TTriaIsoparametric *)F_K)->SetApproximationOrder(polydegree);
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
        ((TQuadBilinear *)F_K)->SetCell(cell);
//         locvol = ((TQuadBilinear *)rt)->GetVolume();
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
//         if(Phase_No==0)
         {
          u1_rise += values[k][l]*U1;
          u2_rise += values[k][l]*U2;
	 }
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




// ======================================================================
// utilities for main program
// ======================================================================

int main(int argc, char* argv[])
{
  TDomain *Domain = new TDomain();
  TDomain *SurfDomain = new TDomain();
  TDatabase *Database = new TDatabase();
  TFEDatabase2D *FEDatabase = new TFEDatabase2D();
  TCollection *coll, *mortarcoll = NULL,  *Surf_Coll;
  TBaseCell *Me, *cell, **Free_Cells;
  TGridCell **DelCell;
  TFESpace2D *velocity_space, *pressure_space, *streamfunction_space, *convolution_space, *fesps;
  TFESpace2D *Grid_space, *vorticity_space, *surfact_space;
  TFESpace1D *surface_space;
  TOutput2D *Output;

  double InitVolume, CurrVolume;
  double *B, *rhs, *sol, *oldsol, *refsol,  *psi, *defect, *SB;
  double *startsol, *frac_step_sol, *Itmetho_sol[2];
  double *oldrhs, *itmethod_sol, *itmethod_rhs, *vorticity, *div, *conv_vort;
  double *sol_timestep_m1;
  int i,j,k,l,m,n, l1, l2, l3, N_, Len, m2, m3, m4, i3, Initial=0;
  int N_Rows, N_Columns, N_U, N_P, N_Unknowns, N_V, N_Vort, N_SO;
  double *l2u1, *l2u2, *h1u1, *h1u2;
  double *l2p, *h1p, *sd, *l_inf, vect[3], exactvect[3];
  int which, *permutation, N_GRAPE_images=0,  N_GNU_images=0;
  double DiffL2, DiffH1;
  double dx0, dy0, dx1, dy1, KE, ST, Surface_Energy;
  char *PRM, *GEO;
  int LEVELS, BASELEVEL, ORDER, order, img=1, N_BData=1;
  int ret, pde;
  int N_Cells, DOF, *GlobalNumbers;
  double negPower;
  double x,y,tx,ty,sx,sy,max,min,sum, TX[2], TY[2];
  double RE_NR;
  double tau1, tau2;
  double Surf_Mass[3], errors[7], Terrors[2], p1, p2, T_inftyL2, T_inftyH1, T_inftyL2_time;
  double t1, t2, res, res2, oldres, solver_time, solver_time_curr, residual, oldresidual;
  double impuls_residual,limit,linredfac, total_time, t3, t4;
  int N_LinIter, N_LinIterCurr, N_LinIterCurrIte, N_SubSteps, N_Active, n_aux;
  double gamma, tau, oldtau, Surfgamma;
  int *RowPtr, N_SurfActive;
  int N_G, N_BoundaryNodes, VSP, n_change=0, N_S;
  double SavedRedFactor, wgt;
  double *vcoarse, *vfine; 
  double y_top[3], t_top[3], max_r;
  int **Cell_No_array, **Joint_No_array;

  std::ostringstream os;

  char *PsBaseName, *GrapeBaseName, *GnuBaseName, *ReadGrapeBaseName;
  char *VtkBaseName;
  double *val, cd, cl, P1, P2;
  TFEFunction2D *u1, *u2, *p, *fefct[7], *StreamFct, *Vorticity;
  TFEFunction1D *fefct_low[2];
  TFEFunction2D *Divergence, *Conv_Vort, *Approx;
  TFEFunction2D *du1Conv, *du2Conv, *du3Conv;
  TFEFunction2D *u1Conv, *u2Conv;
  TFEFunction2D **AuxFEFunctArray;
  TFEFunction2D *GridXDot, *GridYDot;
  TFEVectFunct2D *RefGridPos, *AuxGridPos, *GridPos;
  TFEVectFunct2D *GridVelocity, **GridV, *OldVelo;
  TFEVectFunct2D *u, *uconf, *duConv, **duConvArray;
  TFEVectFunct2D *uConv, **uConvArray, **AuxFEVectFunctArray;
  TFEFunction2D **du1ConvArray, **du2ConvArray, **du3ConvArray;
  TFEFunction2D **u1ConvArray, **u2ConvArray;
  TFEVectFunct2D *GL00AuxProblemSol, **GL00AuxProblemSolArray;
  TFEFunction2D *GL00AuxProblemSol11, *GL00AuxProblemSol12;
  TFEFunction2D *GL00AuxProblemSol22;
  TFEFunction2D **GL00AuxProblemSol11Array, **GL00AuxProblemSol12Array;
  TFEFunction2D **GL00AuxProblemSol22Array, *Surfactant;
  TFEFunction1D *SurfSurfactant;
  TFESpace2D *fesp[3], *ferhs[3];
  TFESpace1D *fesp_low[3], *ferhs_low[3];

  double delta, end_time, l_infty_l_2 = 0, l_infty_l_2_time=-4711.0;
  double olderror = 0, l_2_l_2Du=0, l_2_l_2u=0 , olderror_l_2_l_2u=0;
  double l_2_h_1u=0, olderror_l_2_h_1u=0;
  double L2_norm[2][2], H1_semi[2][2], H1_norm[2][2];

  TAuxParam2D *aux;

  TSquareStructure2D *sqstructureA, *sqstructureC, *SqGridStructure;
  TSquareMatrix2D *SqGridMatrix11, *SqGridMatrix12;
  TSquareMatrix2D *SqGridMatrix21, *SqGridMatrix22;
  TSquareMatrix2D **MatricesG11, **MatricesG12;
  TSquareMatrix2D **MatricesG21, **MatricesG22;
  TSquareMatrix1D **MatricesS_A, **MatricesS_M, **MatricesK;
  TStructure2D *structureB, *structureBT;
  TSquareMatrix2D *sqmatrixA, *SQMATRICES[8], *SQMATRICES_FreeSurf[2], *SQMATRICES_GRID[4];
  TSquareMatrix1D *SQMATRICES_SURF[2];
  TSquareMatrix2D *SQMATRICES_REF[4];
  TSquareMatrix2D *sqmatrixA11, *sqmatrixA12;
  TSquareMatrix2D *sqmatrixA21, *sqmatrixA22;
  TSquareMatrix2D *sqmatrixM;
  TSquareMatrix2D *sqmatrixM11, *sqmatrixM12;
  TSquareMatrix2D *sqmatrixM21, *sqmatrixM22;
  TSquareMatrix2D **MatricesM11, **MatricesM12;
  TSquareMatrix2D **MatricesM21, **MatricesM22;
  TSquareMatrix2D **MatricesA, **MatricesM;
  TSquareMatrix2D **MatricesA11, **MatricesA12;
  TSquareMatrix2D **MatricesA21, **MatricesA22;
  TSquareMatrix2D **MatricesF11, **MatricesF22; 
  TSquareMatrix2D *sqmatrixGL00AuxProblem, **MatricesGL00AuxProblem;
  TMatrix2D *matrixB1, *matrixB2, *MATRICES[4];
  TMatrix2D *matrixB1T, *matrixB2T;
  TMatrix2D **MatricesB1, **MatricesB2, **MatricesB1T, **MatricesB2T;
  TSquareMatrix **sqmatrices = (TSquareMatrix **)SQMATRICES;
  TMatrix **matrices = (TMatrix **)MATRICES;
  MatVecProc *MatVect;
  DefectProc *Defect;

  TNSE_MGLevel *MGLevel;
  TNSE_MultiGrid *MG;
  TMGLevel2D *MGLevelGL00AuxProblem;
  TMultiGrid2D *MGGL00AuxProblem;

  double *RHSs[3];
  int *N_Uarray, *N_Parray, *N_Garray;

  TDiscreteForm2D *DiscreteFormGalerkin;
  TDiscreteForm2D *DiscreteFormGalerkin_Psep;
  TDiscreteForm2D *DiscreteFormColetti;
  TDiscreteForm2D *DiscreteFormGL00Convolution;
  TDiscreteForm2D *DiscreteFormGL00AuxProblem;
  TDiscreteForm2D *DiscreteFormUpwind;
  TDiscreteForm2D *DiscreteFormUpwind_Psep;
  TDiscreteForm2D *DiscreteFormSmagorinsky;
  TDiscreteForm2D *DiscreteFormVMS_Projection;

  TDiscreteForm2D *DiscreteFormNLGalerkin;
  TDiscreteForm2D *DiscreteFormNLUpwind;
  TDiscreteForm2D *DiscreteFormNLColetti;
  TDiscreteForm2D *DiscreteFormNLGL00Convolution;
  TDiscreteForm2D *DiscreteFormNLGL00AuxProblem;
  TDiscreteForm2D *DiscreteFormNLSmagorinsky;
  TDiscreteForm2D *DiscreteFormNLVMS_Projection;

  TDiscreteForm2D *DiscreteFormRHS;
  TDiscreteForm2D *DiscreteFormRHSColetti;
  TDiscreteForm2D *DiscreteFormRHSSmagorinskyExpl;
  TDiscreteForm2D *DiscreteFormMatrixGL00AuxProblem;
  TDiscreteForm2D *DiscreteFormGL00AuxProblemRHS;
  TDiscreteForm2D *DiscreteFormRHSLESModel;
  TDiscreteForm2D *DiscreteFormRHSAuxProblemU;
  TDiscreteForm2D *DiscreteFormMatrixAuxProblemU;
  TDiscreteForm2D *DiscreteFormGrid;

  TDiscreteForm2D *DiscreteForm;

  TDiscreteForm2D *DiscreteForm_Recalc[5];

  int N_SquareMatrices, N_RectMatrices, N_Rhs, N_FESpaces, *Edge_No, N_FESpaces_low;

  BoundCondFunct2D *BoundaryConditions[2], *BoundaryConditionsAuxProblem[3];
  BoundValueFunct2D *BoundValues[2], *BoundValuesAuxProblem[3];
  double average, B_theta=1-(sqrt(2.0)/2.0);

  TItMethod *itmethod, *prec;
  int Max_It, FirstSolve;
  double omega, alpha, alpha_fine, divergence;
  int N_Parameters=1, methods, time_discs;
  double Parameters[2], hmin, hmax, fh, fhtot=0, fhmin = 1e2, fhmax= 0, fhlimit = 5e-3;
  double shtot, shmax, shmin, Params[10], *matvalues;

  int N_UConv, level_down, ii, fully_implicit = 0;
  double *u_conv, *auxConv, *du_tensor, *u_uConv, reatt_pt, vort_zero, vort_zero_conv;
  int mg_level,mg_type,CurrentDiscType, last_sq, step_length;
  int velocity_space_code, pressure_space_code;
  int zerostart, comp_vort=0, mixing_layer=0;
  int N_INNER_DOF, N_Remesh = 0, N_ReConstruct = 0, Remesh = 0;

  BoundCondFunct2D *GridBoundaryConditions[1];
  BoundValueFunct2D *GridBoundValues[1];

  // strings
  char ReadinDat[] = "readin.dat";
  char NameString[] = "name";
  char UString[] = "u";
  char PString[] = "p";
  char PSepString[] = "p_sep";
  char PsiString[] = "psi";
  char SurfString[] = "Surfact";
  char UConvString[] = "u_conv";
  char UConfString[] = "uconf";
  char AuxProbString[] = "AuxProblem";
  char VorticityString[] = "vorticity";
  char ConvVortString[] = "conv_vort";
  char DivergenceString[] = "Divergence";
  char veloString[] = "velo";
  char gridString[] = "grid";
  char refposString[] = "refpos";
  char auxposString[] = "auxpos";
  char posString[] = "pos";

   // char img1[] = "00000";
  // Variables for Triangular grid generation

   int N_Interf_Vertices, N_Interf_Cells;
   char filename[128];
   double  t, teta, dt;


  TBoundComp *BoundComp;
  TBoundPart *BoundPart;
  TBdLine *UpdateSlipBound;
  TBdCircle *UpdateFreeBound;
  TBaseCell **CellTree;

  TVertex **FreeBound_Vert, **SlipBound_Vert, **VertexDel, **NewVertices, **IsoVertices;
  TVertex *temp_Mov, **FreeGauss_Vert, *temp_Mov2;
  TIsoBoundEdge **FreeBound_Joint, *IsoJoint;
  TBoundEdge **SlipBound_Joint, *tempSlip_Joint, *Solid_Joint;
  double FREX, FREY, SLPX, SLPY;
  double R_Theta[3], L_Theta[3], LU[3], RU[3];
  double *gridsol, *gridrhs, *d;
  double *refpos, *auxpos, *pos, *velocity, *tmp, *tmpUP, *refvelocity;
  int N_Joints, N_Vertices, N_SlipJoints, N_FreeJoints;
  int N_Velo_IsoPoints, N_FreeGauss_Vert, Error_corr, very_first_time=0;
  boolean remeshed= TRUE, reparam = FALSE;
  TJoint *Joint;
  TMultiGrid2D *GridMG;
  TMGLevel2D *GridMGLevel;
  int *GridKCol, *GridRowPtr;
  double *Entries[4];
  double friction_constant;
  std::ostringstream opts;

  struct triangulateio In, Out;

  int N_RootCells;
  double left, right, top, bottom, rad1, Rx, Ry, Lx, Ly, x1, x2, y1, y2, x3, y3, x4, y4;
  int *PointNeighb, maxEpV = 0, a, b, Neighb_tmp, Neib[2];
  int CurrNeib, len1, len2, CurrComp, comp;
// auxiliary variables
  int In_Index, CurrVertex, CurrJoint, CurrInterfCell, ID, bct, bct1;
  double  T_a, T_b, temp, temp3, T, *Coordinates, Remesh_Time=0;
  int *PartMarker, *Triangles;
  int N_SlipBound_Vert, m1, N_FreeBound_Vert, temp1, temp2;
  double *Angle = new double[2];

#ifdef __BENCH__
  double Cd, Cl, dP1[3], dP2[3];
#endif
#ifdef __MIXINGLAYERSLIP__
  mixing_layer = 1;
#endif

// // variables for remeshing
  TFESpace2D **FeSpaces = new TFESpace2D*[6];
  TFESpace2D ***FeSps_Lev = new TFESpace2D**[6];
  TFESpace1D ***SurfFeSps_Lev = new TFESpace1D**[2];
  int **N_array = new int*[4];
  int ***N_List = new int**[2];
  TSquareMatrix2D ***SqMat = new TSquareMatrix2D**[16];
  TSquareMatrix1D ***SqMat_low = new TSquareMatrix1D**[3];
  TMatrix2D ***Mat = new TMatrix2D**[4];
  double **Rhs = new double*[3];
  double ***Sol = new double**[9];
//   double ***SurfRhsArray = new double**[1];
  double ***Rhsarray = new double**[4];
  double *srhs, *ssol,  *sdefect, *surfact, *surfact_velo, *ssol_old;
  TFEVectFunct2D ***VeloVect = new TFEVectFunct2D**[2];
  TFEFunction2D ***UPArrays = new TFEFunction2D**[8];
  TFEFunction1D ***SArrays = new TFEFunction1D**[1];
  TVertex ***MovBoundVert = new TVertex**[4];
  int *N_MovVert = new int[3];
  TNSE_MultiGrid **TMG = new TNSE_MultiGrid*[1];
  TNSE_MGLevel ***TMGLevel = new TNSE_MGLevel**[1];
  TBoundEdge ***Slip_Joint = new TBoundEdge**[1];
  TIsoBoundEdge ***Free_Joint = new TIsoBoundEdge**[1];
  TSquareStructure2D ***SqrStruct = new TSquareStructure2D**[2];
  TSquareStructure1D ***SqrStruct_low = new TSquareStructure1D**[1];
  TStructure2D ***Struct = new TStructure2D**[2];

  TFEVectFunct2D *Gridpos[3];

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
//======================================================================
// copy read parameters into local variables
//======================================================================
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
  GrapeBaseName = TDatabase::ParamDB->GRAPEBASENAME;
  GnuBaseName = TDatabase::ParamDB->GNUBASENAME;
  ReadGrapeBaseName = TDatabase::ParamDB->READGRAPEBASENAME;
  VtkBaseName = TDatabase::ParamDB->VTKBASENAME;

  mg_type = TDatabase::ParamDB->SC_MG_TYPE_SADDLE;
  if (mg_type)
    mg_level = 1;
  else
    mg_level = 0;

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

//   Surf_Coll = new TCollection*[LEVELS+1];
//   Cell_No_array = new int*[LEVELS+1];
//   Joint_No_array = new int*[LEVELS+1];

  duConvArray = new TFEVectFunct2D*[LEVELS+1];
  du1ConvArray = new TFEFunction2D*[LEVELS+1];
  du2ConvArray = new TFEFunction2D*[LEVELS+1];
  du3ConvArray = new TFEFunction2D*[LEVELS+1];
  uConvArray = new TFEVectFunct2D*[LEVELS+1];
  u1ConvArray = new TFEFunction2D*[LEVELS+1];
  u2ConvArray = new TFEFunction2D*[LEVELS+1];

  GL00AuxProblemSolArray = new TFEVectFunct2D*[LEVELS+1];
  GL00AuxProblemSol11Array = new TFEFunction2D*[LEVELS+1];
  GL00AuxProblemSol12Array = new TFEFunction2D*[LEVELS+1];
  GL00AuxProblemSol22Array = new TFEFunction2D*[LEVELS+1];

  N_Uarray = new int[LEVELS+1];
  N_Parray = new int[LEVELS+1];
  N_Garray = new int[LEVELS+1];

  GridV = new TFEVectFunct2D*[LEVELS+1];

  switch(TDatabase::ParamDB->NSTYPE)
  {
    case 1:
      MatricesA = new TSquareMatrix2D*[LEVELS+1];
      MatricesM = new TSquareMatrix2D*[LEVELS+1];

      MatricesB1 = new TMatrix2D*[LEVELS+1];
      MatricesB2 = new TMatrix2D*[LEVELS+1];
      MatVect = MatVect_NSE1;
      Defect = Defect_NSE1;
    break;

    case 2:
      MatricesA = new TSquareMatrix2D*[LEVELS+1];
      MatricesM = new TSquareMatrix2D*[LEVELS+1];

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
      MatricesM11 = new TSquareMatrix2D*[LEVELS+1];
      MatricesM12 = new TSquareMatrix2D*[LEVELS+1];
      MatricesM21 = new TSquareMatrix2D*[LEVELS+1];
      MatricesM22 = new TSquareMatrix2D*[LEVELS+1];

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
      MatricesM11 = new TSquareMatrix2D*[LEVELS+1];
      MatricesM12 = new TSquareMatrix2D*[LEVELS+1];
      MatricesM21 = new TSquareMatrix2D*[LEVELS+1];
      MatricesM22 = new TSquareMatrix2D*[LEVELS+1];
      MatricesF11 = new TSquareMatrix2D*[LEVELS+1];
      MatricesF22 = new TSquareMatrix2D*[LEVELS+1];

      MatricesB1 = new TMatrix2D*[LEVELS+1];
      MatricesB2 = new TMatrix2D*[LEVELS+1];
      MatricesB1T = new TMatrix2D*[LEVELS+1];
      MatricesB2T = new TMatrix2D*[LEVELS+1];


      MatricesG11 = new TSquareMatrix2D*[LEVELS+1];
      MatricesG12 = new TSquareMatrix2D*[LEVELS+1];
      MatricesG21 = new TSquareMatrix2D*[LEVELS+1];
      MatricesG22 = new TSquareMatrix2D*[LEVELS+1];

      MatricesS_A = new TSquareMatrix1D*[LEVELS+1];
      MatricesS_M = new TSquareMatrix1D*[LEVELS+1];
      MatricesK = new TSquareMatrix1D*[LEVELS+1];

      MatVect = MatVect_NSE4;
      Defect = Defect_NSE4;
    break;
  } // endswitch

  InitializeDiscreteForms_Moving_axial3D( DiscreteFormGalerkin, DiscreteFormGalerkin_Psep,
                                  DiscreteFormUpwind, DiscreteFormUpwind_Psep,
                                  DiscreteFormNLGalerkin,
                                  DiscreteFormNLUpwind, DiscreteFormRHS,
                                  DiscreteFormGrid, LinCoeffs, GridCoeffs,
                                  TDatabase::ParamDB->NSTYPE);


//======================================================================
// read boundary parameterization and initialize coarse grid
//======================================================================

Domain->Init(PRM, GEO);

//*
      // write grid into an Postscript file
      os.seekp(std::ios::beg);
      os << "Domain_old" << ".ps" << ends;
      Domain->PS(os.str().c_str(),It_Finest,0);
//*/


  boolean AllowEdgeRef = (boolean) TDatabase::ParamDB->MESHGEN_ALLOW_EDGE_REF;
 // cout<< " MESHGEN_ALLOW_EDGE_REF" <<AllowEdgeRef<<endl;
//======================================================================
// Triangular for grid generation begin
//======================================================================
//======================================================================
  BoundPart = Domain->GetBdPart(0);
  UpdateSlipBound = (TBdLine*)BoundPart->GetBdComp(0);
  UpdateFreeBound = (TBdCircle*)BoundPart->GetBdComp(1);

  Out.pointlist = NULL;
  Out.pointattributelist = NULL;
  Out.pointmarkerlist = NULL;
  Out.trianglelist = NULL;
  Out.triangleattributelist = NULL;
  Out.trianglearealist = NULL;
  Out.neighborlist = NULL;
  Out.segmentlist = NULL;
  Out.segmentmarkerlist = NULL;
  Out.holelist = NULL;
  Out.regionlist = NULL;
  Out.edgelist = NULL;
  Out.edgemarkerlist = NULL;
  Out.normlist = NULL;

  opts.seekp(std::ios::beg);

  double area = TDatabase::ParamDB->Area;
  double mode_diff = TDatabase::ParamDB->P4;
  double r, phi, deviation=TDatabase::ParamDB->P11;
  int mode = int (TDatabase::ParamDB->P4);
//OutPut("MESHGEN_REF_QUALIT " << TDatabase::ParamDB->MESHGEN_REF_QUALITY << endl);

 opts<<'p'; // Constrained Delaunay Triangulation:
           // initial values - only points defined on the boundary of the domain;
           // triangulation near boundary may variate from Delaunay criterion
 opts<<"q"<<  TDatabase::ParamDB->MESHGEN_REF_QUALITY;
              // Quality mesh generation with no angles smaller than 20 degrees;

  opts<<"a"<< area; // Imposes a maximum triangle area.
  opts<<'e'; // Outputs a list of edges of the triangulation
  opts<<'z'; // Numbers if items starting from 0
  //opts<<"VVVV"; // Gives detailed information about what Triangle is doing
  opts<<'Q'; // Supress all explanation of what Triangle is doing, unless an error occurs
if(!AllowEdgeRef)
  opts<<'Y'; // Supress adding vertices on boundary edges
  opts<<ends;


  N_FreeBound_Vert = int (TDatabase::ParamDB->P6);    //Freesurf except end point
  N_SlipBound_Vert = 50;     // Initially only three points on solid bound (except end point)
  N_Interf_Vertices = N_FreeBound_Vert+N_SlipBound_Vert;
  In.numberofpoints = N_Interf_Vertices;
  In.pointlist = new double[2*In.numberofpoints];
  In.pointmarkerlist = new int[In.numberofpoints];
  In.numberofpointattributes = 0;

  In.numberofsegments = In.numberofpoints;
  In.segmentlist = new int[2*In.numberofsegments];
  In.segmentmarkerlist = new int[In.numberofsegments];
  In.numberofholes = 0;
  In.holelist = NULL;
  In.numberofregions = 0;
  In.regionlist = NULL;

  In_Index = 0;
  CurrComp = 1;

double y_begin, y_end, dx;

    T_a = 1.; // x axis value in ellipse
    T_b = 1.; // y axis value in ellipse !!! set in modifyGausscoord funct also
    max_r = 1.;
    x = 0.0;     // Axial bound

    y = T_b;
// // spherical harmonic of order 2 
// // spherical harmonic of order 4
//     phi = atan2(y, x);
//     if(mode==2)
//       r = 1.0 + deviation*(sqrt(5./Pi)*(3.*cos(phi+Pi/2.)*cos(phi+Pi/2.) -1. )/4.);
//     else if(mode==4) 
//      {
//       temp = cos(phi+Pi/2.);
//       r = 1.0 + deviation*(sqrt(1./Pi)*(35.*temp*temp*temp*temp - 30.*temp*temp +3. )*3./16.);
//      }
//     else 
//       {
//         OutPut("No. of mode is taken as 0 check main programme"<<endl);
//        exit(0);
//        }

//     y = r*sin(phi);
    y_begin = y;



    y = -T_b;
//     phi = atan2(y, x);
// // spherical harmonic of order 2
// // spherical harmonic of order 4
//     phi = atan2(y, x);
//     if(mode==2)
//       r = 1.0 + deviation*(sqrt(5./Pi)*(3.*cos(phi+Pi/2.)*cos(phi+Pi/2.) -1. )/4.);
//     else if(mode==4) 
//      {
//       temp = cos(phi+Pi/2.);
//       r = 1.0 + deviation*(sqrt(1./Pi)*(35.*temp*temp*temp*temp - 30.*temp*temp +3. )*3./16.);
//      }
// 
//     y = r*sin(phi);

//     ModifyCoords( x,  y);
    y_end = y;

    dx = (y_end - y_begin)/N_SlipBound_Vert;
    y = y_begin;
  // points and segments on the vertical boundary (marker=1)
  for(i=0;i<N_SlipBound_Vert;i++) // without last point
   {


    if (fabs(y)<1e-10) y = 0.;
    if (fabs(x)<1e-10) x = 0.;
    In.pointlist[2*In_Index] = 0.0;
    In.pointlist[2*In_Index+1] = y;
    cout<<" x : "<< x << " y : "<< y<<endl;
    In.pointmarkerlist[In_Index] = CurrComp;
    In.segmentlist[2*In_Index] = In_Index;
    In.segmentlist[2*In_Index+1] = In_Index+1;
    In.segmentmarkerlist[In_Index] = CurrComp;
    In_Index++;
    y = y_begin + double(i+1)*dx;
   }
   CurrComp++;

cout<<endl;


    teta = -Pi/2.; // end degree value of freesurface
    dt = Pi/N_FreeBound_Vert;

 // points and segments on the interface (marker=2)
    t = teta;

   for(i=0;i<N_FreeBound_Vert;i++) // without last point
    {
    //  cout<<" teta : "<< teta <<endl;
      x = T_a*cos(teta);
      y = T_b*sin(teta);

// spherical harmonic of order 2
// spherical harmonic of order 4
//     phi = atan2(y, x);
//     if(mode==2)
//       r = 1.0 + deviation*(sqrt(5./Pi)*(3.*cos(phi+Pi/2.)*cos(phi+Pi/2.) -1. )/4.);
//     else if(mode==4) 
//      {
//       temp = cos(phi+Pi/2.);
//       r = 1.0 + deviation*(sqrt(1./Pi)*(35.*temp*temp*temp*temp - 30.*temp*temp +3. )*3./16.);
//      }

//       x = r*cos(phi);
//       y = r*sin(phi);
//       ModifyCoords( x,  y);
      if (fabs(y)<1e-10) y = 0.;
      if (fabs(x)<1e-10) x = 0.;

      In.pointlist[2*In_Index] = x;
      In.pointlist[2*In_Index+1] = y;
      cout<<" x : "<< x << " y : "<< y<<endl;
      //In.pointmarkerlist[In_Index] = CurrComp;
      In.segmentlist[2*In_Index] = In_Index;
      In.segmentlist[2*In_Index+1] = In_Index+1;
      if (AllowEdgeRef)
       {
    //    if (i)
         In.segmentmarkerlist[In_Index] = CurrComp;
       }
      else
       {
    //    if (i)
          In.segmentmarkerlist[In_Index] = 100000 + CurrComp;
       }

      In_Index++;
      teta = t + double(i+1)*dt;
//       if(teta>360) teta -= 360.;
    }
  //  CurrComp++;

  In.segmentlist[2*(In.numberofsegments-1)+1] = 0;

// exit(0);
  if(Out.pointlist!=NULL) {
    free(Out.pointlist); Out.pointlist = NULL;}
  if(Out.pointattributelist!=NULL) {
    free(Out.pointattributelist); Out.pointattributelist = NULL;}
  if(Out.pointmarkerlist!=NULL) {
    free(Out.pointmarkerlist); Out.pointmarkerlist = NULL;}
  if(Out.trianglelist!=NULL) {
    free(Out.trianglelist); Out.trianglelist = NULL;}
  if(Out.triangleattributelist!=NULL) {
    free(Out.triangleattributelist); Out.triangleattributelist = NULL;}
  if(Out.trianglearealist!=NULL) {
    free(Out.trianglearealist); Out.trianglearealist = NULL;}
  if(Out.neighborlist!=NULL) {
    free(Out.neighborlist); Out.neighborlist = NULL;}
  if(Out.segmentlist!=NULL) {
    free(Out.segmentlist); Out.segmentlist = NULL;}
  if(Out.segmentmarkerlist!=NULL) {
    free(Out.segmentmarkerlist); Out.segmentmarkerlist = NULL;}
  if(Out.holelist!=NULL) {
    free(Out.holelist); Out.holelist = NULL;}
  if(Out.regionlist!=NULL) {
    free(Out.regionlist); Out.regionlist = NULL;}
  if(Out.edgelist!=NULL) {
    free(Out.edgelist); Out.edgelist = NULL;}
  if(Out.edgemarkerlist!=NULL) {
    free(Out.edgemarkerlist); Out.edgemarkerlist = NULL;}
  if(Out.normlist!=NULL) {
    free(Out.normlist); Out.normlist = NULL;}

/*
  for(i=0;i<In.numberofpoints;i++)
    OutPut(i<<' '<<In.pointmarkerlist[i]<<' '<<
	   In.pointlist[2*i]<<' '<<In.pointlist[2*i+1]<<endl);
cout<<endl;
//exit(0);
*/
  triangulate((char*)opts.str().c_str(), &In, &Out, (struct triangulateio *)NULL);

/*
  for(i=0;i<Out.numberofpoints;i++)
     OutPut(i<<' '<<Out.pointmarkerlist[i]<<' '<<
      Out.pointlist[2*i]<<' '<<Out.pointlist[2*i+1]<<endl);
  */

  Domain->GetTreeInfo(CellTree,N_RootCells);
  coll = Domain->GetCollection(It_Finest, 0);
  N_Cells = coll->GetN_Cells();

  // remove all existing vertices and joints
  VertexDel = new TVertex*[3*N_RootCells];
  // DelCell = new TGridCell*[N_Cells];

  CurrVertex = 0;

  for(i=0;i<N_Cells;i++)
    {
      cell = coll->GetCell(i);
      N_Joints = cell->GetN_Joints();
      N_Vertices = cell->GetN_Vertices();
      for(j=0;j<N_Joints;j++)
        {
         if(CurrVertex==0)
          {
              VertexDel[CurrVertex] = cell->GetVertex(j);
              CurrVertex++;
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
              VertexDel[CurrVertex] = cell->GetVertex(j);
              CurrVertex++;
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
            VertexDel[CurrVertex] = cell->GetVertex((j+1)%N_Vertices);
            CurrVertex++;
           }
        } // for j
    } // for i
  for(i=0;i<CurrVertex;i++)
    delete VertexDel[i];

  delete []VertexDel;
  OutPut(CurrVertex<<" vertices were deleted"<<endl);

 // remove all existing cells and joints
  for(i=0;i<N_RootCells;i++)
    delete (TGridCell*)CellTree[i];
  OutPut(N_RootCells<<" cells were deleted"<<endl);
   delete CellTree;
   delete coll;

  N_RootCells = Out.numberoftriangles;

  // allocate auxillary fields
  Coordinates = Out.pointlist;
  Triangles = Out.trianglelist;
  PartMarker = new int[Out.numberofpoints];

  // generate new vertices
  N_G = Out.numberofpoints;
  NewVertices = new TVertex*[N_G];

  for (i=0;i<N_G;i++)
     NewVertices[i] = new TVertex(Coordinates[2*i], Coordinates[2*i+1]);

      // set bounding box
  left = bottom = 1e8;
  right = top = -1e8;

   for(i=0;i<In.numberofpoints;i++)
    {
      if(left>In.pointlist[2*i]) left = In.pointlist[2*i];
      if(right<In.pointlist[2*i]) right = In.pointlist[2*i];
      if(top<In.pointlist[2*i+1]) top = In.pointlist[2*i+1];
      if(bottom>In.pointlist[2*i+1]) bottom = In.pointlist[2*i+1];
    }

  //OutPut("left: "<<left<<" right: "<<right<<" top: "<<top<<" bottom: "<<bottom<<endl);

  Domain->SetBoundBox(right-left,top-bottom);
  Domain->SetBoundBoxstart(left,bottom);

 // Solid Bound startx, starty, x length and y length
 UpdateSlipBound->SetParams(In.pointlist[0],In.pointlist[1],
                          In.pointlist[2*N_SlipBound_Vert]-In.pointlist[0],
                          In.pointlist[2*N_SlipBound_Vert+1]-In.pointlist[1]);

// Free boundary xmid, ymid, radius_a, radius_b, start angle, end angle
 UpdateFreeBound ->SetParams( 0.0, 0.0, T_a, T_b, -Pi/2., Pi/2.);
 
//  cout << " In.pointlist[N_SlipBound_Vert] : " << In.pointlist[2*N_SlipBound_Vert]<< endl;
 // cout << " In.pointlist[N_SlipBound_Vert] : " << In.pointlist[0]<< endl;
 // generate cells
   CellTree = new TBaseCell*[N_RootCells];

  for (i=0;i<N_RootCells;i++)
  {
    CellTree[i] = new TMacroCell(TDatabase::RefDescDB[Triangle], 0);

    CellTree[i]->SetVertex(0, NewVertices[Out.trianglelist[3*i    ]]);
    CellTree[i]->SetVertex(1, NewVertices[Out.trianglelist[3*i + 1]]);
    CellTree[i]->SetVertex(2, NewVertices[Out.trianglelist[3*i + 2]]);

      ((TMacroCell *) CellTree[i])->SetSubGridID(0);
  }

  Domain->SetTreeInfo(CellTree, N_RootCells);
 
  // initialize iterators
  TDatabase::IteratorDB[It_EQ]->SetParam(Domain);
  TDatabase::IteratorDB[It_LE]->SetParam(Domain);
  TDatabase::IteratorDB[It_Finest]->SetParam(Domain);
  TDatabase::IteratorDB[It_Between]->SetParam(Domain);
  TDatabase::IteratorDB[It_OCAF]->SetParam(Domain);

  // search neighbours
  N_G = Out.numberofpoints;
  PointNeighb = new int[N_G];

  memset(PointNeighb, 0, N_G *SizeOfInt);

  for (i=0;i<3*N_RootCells;i++)
    PointNeighb[Triangles[i]]++;

  for (i=0;i<N_G;i++)
    if (PointNeighb[i] > maxEpV) maxEpV = PointNeighb[i];
  delete [] PointNeighb;

  PointNeighb = new int[++maxEpV * N_G];

  memset(PointNeighb, 0, maxEpV * N_G *SizeOfInt);

   // first colomn contains the number of following elements
   // for every point at first column we set the number of neighbour points
   // at further columns we set the index of corresponding cells
  for(i=0;i<3*N_RootCells;i++)
  {
    j = Triangles[i]*maxEpV;
    PointNeighb[j]++;
    PointNeighb[j + PointNeighb[j]] = i / 3;
  }

  // generate new edges
  N_G = Out.numberofedges;
  for (i=0;i<N_G;i++)
  {
    a = Out.edgelist[2*i];
    b = Out.edgelist[2*i+1];
    Neib[0] = -1;
    Neib[1] = -1;
    CurrNeib = 0;

    len1 = PointNeighb[a*maxEpV];
    len2 = PointNeighb[b*maxEpV];

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

    if (Out.edgemarkerlist[i]) // 0 for inner edges and Boundcomp+1 for Boundedge respect
    {
      CurrComp = Out.edgemarkerlist[i] - 1;
      if (CurrComp >= 100000) CurrComp -= 100000;


      if(Domain->GetBdPart(0)->GetBdComp(CurrComp)->GetTofXY(
            NewVertices[a]->GetX(), NewVertices[a]->GetY(), T_a) ||
          Domain->GetBdPart(0)->GetBdComp(CurrComp)->GetTofXY(
            NewVertices[b]->GetX(), NewVertices[b]->GetY(), T_b))
       {
          cerr<<"Error: could not set parameter values"<<endl;
          OutPut(NewVertices[a]<<endl);
          OutPut(NewVertices[b]<<endl);
        //  exit(0);
       }

      if (CurrNeib == 2)    // 2 cells contain the current edge
        if(Domain->GetBdPart(0)->GetBdComp(CurrComp)->IsFreeBoundary())
          Joint = new TIsoInterfaceJoint(Domain->GetBdPart(0)->GetBdComp(CurrComp),
                  T_a, T_b, CellTree[Neib[0]], CellTree[Neib[1]]);
        else
          Joint = new TInterfaceJoint(Domain->GetBdPart(0)->GetBdComp(CurrComp),
                  T_a, T_b, CellTree[Neib[0]], CellTree[Neib[1]]);
      else
        if(Domain->GetBdPart(0)->GetBdComp(CurrComp)->IsFreeBoundary())
          Joint = new TIsoBoundEdge(Domain->GetBdPart(0)->GetBdComp(CurrComp), T_a, T_b);
        else
          Joint = new TBoundEdge(Domain->GetBdPart(0)->GetBdComp(CurrComp), T_a, T_b);
    }
    else // inner edge
    {
    if (CurrNeib != 2)
        cerr << "Error!!!!!!!! not enough neighbours!" << endl;

    Joint = new TJointEqN(CellTree[Neib[0]], CellTree[Neib[1]]);
    }

    // find the local index for the point 'a' on the cell
    for (j=0;j<3;j++)
      if (Triangles[3*Neib[0]+j] == a) break;

    // find the local index for the point 'b' on the cell
    for (k=0;k<3;k++)
      if (Triangles[3*Neib[0]+k] == b) break;

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

  CellTree[Neib[0]]->SetJoint(j, Joint);

   if (Neib[1] != -1)
    {
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

      CellTree[Neib[1]]->SetJoint(j, Joint);
    }

  if (Joint->GetType() == InterfaceJoint ||
        Joint->GetType() == IsoInterfaceJoint)
      ((TInterfaceJoint *) Joint)->CheckOrientation();
  }


  delete [] NewVertices;
  delete [] PointNeighb;
  delete [] In.pointlist;
  delete [] In.pointmarkerlist;
  delete [] In.segmentlist;
  delete [] In.segmentmarkerlist;
  
  if(Out.pointlist!=NULL) {
    free(Out.pointlist); Out.pointlist = NULL;}
  if(Out.pointattributelist!=NULL) { 
    free(Out.pointattributelist); Out.pointattributelist = NULL;}
  if(Out.pointmarkerlist!=NULL) {
    free(Out.pointmarkerlist); Out.pointmarkerlist = NULL;}
  if(Out.trianglelist!=NULL) {
    free(Out.trianglelist); Out.trianglelist = NULL;}
  if(Out.triangleattributelist!=NULL) {
    free(Out.triangleattributelist); Out.triangleattributelist = NULL;}
  if(Out.trianglearealist!=NULL) {
    free(Out.trianglearealist); Out.trianglearealist = NULL;}
  if(Out.neighborlist!=NULL) {
    free(Out.neighborlist); Out.neighborlist = NULL;}
  if(Out.segmentlist!=NULL) {
    free(Out.segmentlist); Out.segmentlist = NULL;}
  if(Out.segmentmarkerlist!=NULL) {
    free(Out.segmentmarkerlist); Out.segmentmarkerlist = NULL;}
  if(Out.holelist!=NULL) {
    free(Out.holelist); Out.holelist = NULL;}
  if(Out.regionlist!=NULL) {
    free(Out.regionlist); Out.regionlist = NULL;}
  if(Out.edgelist!=NULL) {
    free(Out.edgelist); Out.edgelist = NULL;}
  if(Out.edgemarkerlist!=NULL) {
    free(Out.edgemarkerlist); Out.edgemarkerlist = NULL;}
  if(Out.normlist!=NULL) {
    free(Out.normlist); Out.normlist = NULL;} 
  

//======================================================================
// Triangular for grid generation end
//======================================================================

  // write grid into an Postscript file
  os.seekp(std::ios::beg);
  os << "Domain" << ".ps" << ends;
  Domain->PS(os.str().c_str(),It_Finest,0);
//  exit(0);

  BoundaryConditions[0] = BoundCondition;
  BoundaryConditions[1] = BoundCondition;


  BoundValues[0] = U1BoundValue;
  BoundValues[1] = U2BoundValue;

  GridBoundaryConditions[0] = GridBoundCondition;
  GridBoundValues[0] = GridBoundValue;

  BoundaryConditionsAuxProblem[0] = BoundConditionAuxProblem;
  BoundaryConditionsAuxProblem[1] = BoundConditionAuxProblem;
  BoundaryConditionsAuxProblem[2] = BoundConditionAuxProblem;

  BoundValuesAuxProblem[0] = BoundValueAuxProblem;
  BoundValuesAuxProblem[1] = BoundValueAuxProblem;
  BoundValuesAuxProblem[2] = BoundValueAuxProblem;

 for(i=0;i<TDatabase::ParamDB->SC_COARSEST_LEVEL_SADDLE;i++)
    Domain->RegRefineAll();

  limit = TDatabase::ParamDB->SC_NONLIN_RES_NORM_MIN_SADDLE;
  Max_It = TDatabase::ParamDB->SC_NONLIN_MAXIT_SADDLE;
  omega = TDatabase::ParamDB->SC_NONLIN_DAMP_FACTOR_SADDLE;
  alpha = TDatabase::ParamDB->SC_GMG_DAMP_FACTOR_SADDLE;
  alpha_fine = TDatabase::ParamDB->SC_GMG_DAMP_FACTOR_FINE_SADDLE;
  divergence = TDatabase::ParamDB->SC_DIV_FACTOR;

  Parameters[0] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_SADDLE;
  Parameters[1] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_FINE_SADDLE;

  TDatabase::TimeDB->CURRENTTIME = TDatabase::TimeDB->STARTTIME;

  if(TDatabase::ParamDB->SOLVER_TYPE == GMG)
  {
    i=1;
    MG = new TNSE_MultiGrid(i, N_Parameters, Parameters);
    GridMG = new TMultiGrid2D(i, N_Parameters, Parameters);
  }

  t3 = GetTime();
  total_time = t3 - total_time;
  mg_level = LEVELS+mg_level;

  for(i=0;i<5;i++) Sol[i] = new double*[mg_level];
  Sol[5] = new double*[1];
  for(i=6;i<9;i++) Sol[i] = new double*[mg_level];

  for(i=0;i<2;i++)  VeloVect[i] = new TFEVectFunct2D*[mg_level];
//       VeloVect[2] = new TFEVectFunct2D*[1];
  for(i=0;i<4;i++) N_array[i] = new int[mg_level];
  for(i=0;i<2;i++) N_List[i] = new int*[mg_level];
  for(i=0;i<2;i++) SqrStruct[i] = new TSquareStructure2D*[mg_level];
  SqrStruct_low[0] = new TSquareStructure1D*[mg_level];
  for(i=0;i<2;i++) Struct[i] = new TStructure2D*[mg_level];
  for(i=0;i<4;i++) Rhsarray[i] = new double*[mg_level];
//   for(i=0;i<1;i++) SurfRhsArray[i] = new double*[mg_level];
  for(i=0;i<8;i++) UPArrays[i] = new TFEFunction2D*[mg_level];
  for(i=0;i<1;i++) SArrays[i] = new TFEFunction1D*[mg_level];

  TMGLevel[0] = new TNSE_MGLevel*[mg_level];

  for(i=0;i<4;i++) Mat[i] = new TMatrix2D*[mg_level];
  for(i=0;i<16;i++) SqMat[i] = new TSquareMatrix2D*[mg_level];
  for(i=0;i<3;i++) SqMat_low[i] = new TSquareMatrix1D*[mg_level];

  for(i=0;i<6;i++) FeSps_Lev[i] = new TFESpace2D*[mg_level];
  for(i=0;i<2;i++) SurfFeSps_Lev[i] = new TFESpace1D*[mg_level];
  ORDER = 0;

  VSP = TDatabase::ParamDB->VELOCITY_SPACE;

if (abs(VSP) > 20)
  {ORDER = abs(VSP) - 20;}
else if ( abs(VSP) > 10)
    {ORDER = abs(VSP) - 10;}
  else ORDER = abs(VSP);

 FE1D *FE1D_List;
 int N_SCells;
//======================================================================
// loop over all levels
//======================================================================
  for(i=0;i<mg_level;i++)
  {
    if (i<LEVELS)
    {
      OutPut("*******************************************************" << endl);
      OutPut("******           GEOMETRY  LEVEL ");
      OutPut(i << "              *******" << endl);
    }
    else
    {
      OutPut("*******************************************************" << endl);
      OutPut("******           GEOMETRY  LEVEL ");
      OutPut(i-1 << "              *******" << endl);
    }
    OutPut("******           MULTIGRID LEVEL ");
    OutPut(i << "              *******" << endl);
    OutPut("*******************************************************" << endl);
    OutPut("memory before: " << setw(10) << GetMemory() << endl);

   if(i && (i<LEVELS)) Domain->RegRefineAll();

    coll=Domain->GetCollection(It_Finest, 0);
    cout << endl << endl;
    N_Cells = coll->GetN_Cells();
      
    Cell_No_array = N_List[0];
    Joint_No_array = N_List[1];
   
    Cell_No_array[i] = new int[N_Cells];
    Joint_No_array[i] = new int[N_Cells];
   
    Domain2DSurface(Domain, SurfDomain, Cell_No_array, Joint_No_array, i);
//     N_List[0] = Cell_No_array;
//     N_List[1] = Joint_No_array;

    Surf_Coll = SurfDomain->GetCollection(It_Finest, 0);
    N_SCells= Surf_Coll->GetN_Cells();
    FE1D_List = new FE1D[N_SCells];
    for(j=0;j<N_SCells;j++)
     FE1D_List[j] = FE1D(TDatabase::ParamDB->ANSATZ_ORDER);

    TDatabase::IteratorDB[It_EQ]->SetParam(Domain);
    TDatabase::IteratorDB[It_LE]->SetParam(Domain);
    TDatabase::IteratorDB[It_Finest]->SetParam(Domain);
    TDatabase::IteratorDB[It_Between]->SetParam(Domain);
    TDatabase::IteratorDB[It_OCAF]->SetParam(Domain);

//     if(TDatabase::ParamDB->WRITE_PS)
    {
      // write grid into an Postscript file
      os.seekp(std::ios::beg);
      os << PsBaseName << i << ".ps" << ends;
      Domain->PS(os.str().c_str(),It_Finest,0);
    }

    // get spaces for low order disc on finest geo grid
    if ((mg_type==1)&&(i<mg_level-1))
    {

      velocity_space = new TFESpace2D(coll,NameString, UString, BoundCondition,
                                      Non_USpace, 1, NULL);


      pressure_space = new TFESpace2D(coll,NameString, PString, BoundCondition,
                                      DiscP_PSpace, 0, NULL);

      velocity_space_code = -1;
      pressure_space_code = 0;
      order = -1;
    }
    // get spaces of high order disc on finest geo grid
    else
    {
      GetVelocityAndPressureSpace(coll,BoundCondition,
                                  mortarcoll, velocity_space,
                                  pressure_space, &pressure_space_code,
                                  TDatabase::ParamDB->VELOCITY_SPACE,
                                  TDatabase::ParamDB->PRESSURE_SPACE);

      velocity_space_code =   TDatabase::ParamDB->VELOCITY_SPACE;
      TDatabase::ParamDB->INTERNAL_PRESSURE_SPACE = pressure_space_code;

    }

    Grid_space = new TFESpace2D(coll, NameString, PsiString, GridBoundCondition,
                                1, NULL);


    surface_space = new TFESpace1D(Surf_Coll , NameString, SurfString,
                                   FE1D_List);

    surfact_space = new TFESpace2D(coll , NameString, SurfString, GridBoundCondition,
                                 TDatabase::ParamDB->ANSATZ_ORDER,   NULL);
    N_SO =  surfact_space->GetN_DegreesOfFreedom();

    // build fespace hierarchy
    // set values and pointers for fe space

    FeSps_Lev[0][i] = velocity_space;
    FeSps_Lev[1][i] = pressure_space;
    FeSps_Lev[2][i] = Grid_space;
    FeSps_Lev[5][i] = surfact_space;
    SurfFeSps_Lev[0][i] = surface_space;

    N_U = velocity_space->GetN_DegreesOfFreedom();
    N_P = pressure_space->GetN_DegreesOfFreedom();
    N_G = FeSps_Lev[2][i]->GetN_DegreesOfFreedom();
    N_BoundaryNodes = N_G - FeSps_Lev[2][i]->GetN_Inner();
    N_S = SurfFeSps_Lev[0][i]->GetN_DegreesOfFreedom();

    N_Uarray[i] = N_U;
    N_Parray[i] = N_P;
    N_Garray[i] = N_G;


    N_Active = velocity_space->GetActiveBound();
    cout << "N_ActiveBound: " << N_Active << endl;
    

 //  if (i<LEVELS)
    {
    
    streamfunction_space = new TFESpace2D(coll,NameString, PsiString, BoundCondition,
                                          1, NULL);

    vorticity_space = new TFESpace2D(coll,NameString, UString, BoundCondition,
                                      ContP_USpace,1, mortarcoll);

      FeSps_Lev[3][i] = streamfunction_space;
      N_V = streamfunction_space->GetN_DegreesOfFreedom();
      FeSps_Lev[4][i] = vorticity_space;
      N_Vort = vorticity_space->GetN_DegreesOfFreedom();
    }

    cout << "No. of Grid BoundaryNodes: " << N_BoundaryNodes << endl;
    cout << "No. of DOF of Grid cells: " << N_G << endl;

     // find on finest grid all vertices on moving boundary
    if(i==mg_level-1)
     {
      N_Vertices = 0;
      N_SlipJoints = 0;
      N_FreeJoints = 0;
      N_FreeGauss_Vert = 0;
      N_Velo_IsoPoints = 0;

      for(j=0;j<N_Cells;j++)
      {
        Me = coll->GetCell(j);
        k = Me->GetN_Edges();
        for(l=0;l<k;l++)
        {
          Joint = Me->GetJoint(l);
           if(Joint->GetType() == BoundaryEdge)
            {
               N_SlipJoints++,
               N_Velo_IsoPoints += ORDER-1;
               N_Vertices++;

             }
          if(Joint->GetType() == IsoBoundEdge)
          {
              N_FreeJoints++;
              N_FreeGauss_Vert += ORDER-1;
              N_Velo_IsoPoints += ORDER-1;
              N_Vertices++;
         }
        } // endfor l
      } // endfor j

      N_FreeBound_Vert = N_FreeGauss_Vert;
      cout << "Total N_BoundVertices: " << N_Vertices << endl;
      cout << "N_SlipJoints: " << N_SlipJoints << endl;
      cout << "N_FreeJoints: " << N_FreeJoints << endl;
      cout << "N_FreeGauss_Vert: " << N_FreeGauss_Vert << endl;
      OutPut("Number of Elements = " <<N_Cells<< endl);

      N_SlipBound_Vert = N_SlipJoints+1; //addding the right wetting point to BC2(slip)
      MovBoundVert[1] = new TVertex*[N_FreeBound_Vert];
      MovBoundVert[0] = new TVertex*[N_SlipBound_Vert+N_FreeBound_Vert];
      MovBoundVert[2] = new TVertex*[N_FreeGauss_Vert];
      MovBoundVert[3] = new TVertex*[ORDER];

      SlipBound_Vert = MovBoundVert[0];
      FreeBound_Vert = MovBoundVert[1];
      FreeGauss_Vert = MovBoundVert[2];
      IsoVertices = MovBoundVert[3];

      Free_Joint[0] = new TIsoBoundEdge*[N_FreeJoints];
      Slip_Joint[0] = new TBoundEdge*[N_SlipJoints+N_FreeJoints];

      FreeBound_Joint = Free_Joint[0];
      SlipBound_Joint = Slip_Joint[0];

      Free_Cells = new TBaseCell*[10*N_FreeJoints];
      N_array[3] = new int [N_FreeJoints];
      Edge_No = N_array[3];

      m = 0;
      m1 = 0;
      m2 = 0;
      m3 = 0;
      m4 = 0;

    for(j=0;j<N_Cells;j++)
      {
        Me = coll->GetCell(j);
        k = Me->GetN_Edges();
        for(l=0;l<k;l++)
        {
          Joint = Me->GetJoint(l);
          if(Joint->GetType() == BoundaryEdge)
          {
              SlipBound_Joint[m4++] = (TBoundEdge *)Joint;
              SlipBound_Vert[m++] = Me->GetVertex(l);
          }
         else if(Joint->GetType() == IsoBoundEdge)
          {
            Me->GetVertex(l)->GetCoords(TX[0], TY[0]);
            Me->GetVertex((l+1) % k)->GetCoords(TX[1], TY[1]);

            FreeBound_Vert[m1] = Me->GetVertex(l);
            Free_Cells[m1] = Me;
            Edge_No[m1++] = l;

            FreeBound_Joint[m2] = (TIsoBoundEdge *)Joint;
//             if(i==mg_level-1)
//              {
//               FreeBound_Joint[m2]->GenerateVertices(ORDER-1);
              FreeBound_Joint[m2]->GeneratemidVert(ORDER-1, TX, TY);
              IsoVertices = FreeBound_Joint[m2]->GetVertices();
              for(i3=0;i3<(ORDER-1);i3++)
               {
//                 IsoVertices[i3]->GetCoords(IsoX[m3], IsoY[m3]);
                FreeGauss_Vert[m3++] =  IsoVertices[i3];
               }
//               }
           m2++;
         }
       } // endfor l
      } // endfor j

   N_FreeBound_Vert = m1;
   // sort free bound vertices from -Pi/2 to +Pi/2
   SortIsoVertices_axial3D(Free_Cells, FreeBound_Vert, Edge_No, N_FreeJoints);
   // sort free bound isovertices after SortIsoVertices() only
   // according to the Free_Cells, Edge_No ordered in SortIsoVertices(()

//    SortQuardIsopts(Free_Cells, FreeGauss_Vert, Edge_No, N_FreeGauss_Vert);


  // sort bound points on axial surface in top to down order
    for(k=0;k<N_SlipBound_Vert-2;k++)
      {
      for(l=k+1;l<N_SlipBound_Vert-1;l++)
      {
        SlipBound_Vert[k]->GetCoords(SLPX, SLPY);
	SlipBound_Vert[l]->GetCoords(FREX, FREY);
	if(FREY > SLPY)
	 {
	  temp_Mov = SlipBound_Vert[k];
          SlipBound_Vert[k] = SlipBound_Vert[l];
          SlipBound_Vert[l] = temp_Mov;

	  tempSlip_Joint = SlipBound_Joint[k];
	  SlipBound_Joint[k] = SlipBound_Joint[l];
	  SlipBound_Joint[l] = tempSlip_Joint;
	 }
        }
       }
     //Adding the right wetting point
     SlipBound_Vert[N_SlipBound_Vert-1] = FreeBound_Vert[0];

     for(k=1;k<N_FreeBound_Vert;k++)
      FreeBound_Vert[k-1] = FreeBound_Vert[k];

      N_FreeBound_Vert--;



//      for(k=0;k<N_SlipBound_Vert;k++)
//       {
//        SlipBound_Vert[k]->GetCoords(SLPX, SLPY);
//        cout<< " SLPX " << SLPX<<" SLPY " << SLPY<<endl;
//        }
//
// exit(0);
   /*
     for(k=0;k<N_FreeGauss_Vert;k++)
      {
//        FreeBound_Vert[k]->GetCoords(FREX, FREY);
       FreeGauss_Vert[k]->GetCoords(SLPX, SLPY);

//        OutPut(setw(5)<< k << "   Freesurf:Iso_points: " << setw(15) << FREX);
//        OutPut(setw(15) << FREY << setw(15) << SLPX << setw(15) << SLPY << endl);
        OutPut(k << "FreeGaus " << (180/Pi)*atan2(SLPY, SLPX) <<endl);
       }
    */


//     delete [] Edge_No;
   } // endif i==mg_level-1

 // build matrices

    SqGridStructure = new TSquareStructure2D(FeSps_Lev[2][i]);
    SqGridStructure->Sort();
    structureB = new TStructure2D(pressure_space, velocity_space);
    structureBT = new TStructure2D(velocity_space, pressure_space);
    sqstructureA = new TSquareStructure2D(velocity_space);
    sqstructureA->Sort();

    SqrStruct[1][i] = SqGridStructure;
    SqrStruct[0][i] = sqstructureA;
    Struct[0][i] = structureB;
    Struct[1][i] = structureBT;

    TSquareStructure1D *sqSurfStructure;
    sqSurfStructure = new TSquareStructure1D(SurfFeSps_Lev[0][i]);
    sqSurfStructure->Sort();
    SqrStruct_low[0][i] = sqSurfStructure;

    switch(TDatabase::ParamDB->NSTYPE)
    {
      case 1:
        matrixB1 = new TMatrix2D(structureB);
        matrixB2 = new TMatrix2D(structureB);

        MatricesB1[i] = matrixB1;
        MatricesB2[i] = matrixB2;

        sqmatrixA = new TSquareMatrix2D(sqstructureA);
        MatricesA[i] = sqmatrixA;

        sqmatrixM = new TSquareMatrix2D(sqstructureA);
        MatricesM[i] = sqmatrixM;
        break;

      case 2:
        matrixB1 = new TMatrix2D(structureB);
        matrixB2 = new TMatrix2D(structureB);
        matrixB1T = new TMatrix2D(structureBT);
        matrixB2T = new TMatrix2D(structureBT);

        MatricesB1[i] = matrixB1;
        MatricesB2[i] = matrixB2;
        MatricesB1T[i] = matrixB1T;
        MatricesB2T[i] = matrixB2T;

        sqmatrixA = new TSquareMatrix2D(sqstructureA);
        MatricesA[i] = sqmatrixA;

        sqmatrixM = new TSquareMatrix2D(sqstructureA);
        MatricesM[i] = sqmatrixM;
        break;

      case 3:
        matrixB1 = new TMatrix2D(structureB);
        matrixB2 = new TMatrix2D(structureB);

        MatricesB1[i] = matrixB1;
        MatricesB2[i] = matrixB2;

        sqmatrixA11 = new TSquareMatrix2D(sqstructureA);
        sqmatrixA12 = new TSquareMatrix2D(sqstructureA);
        sqmatrixA21 = new TSquareMatrix2D(sqstructureA);
        sqmatrixA22 = new TSquareMatrix2D(sqstructureA);

        MatricesA11[i] = sqmatrixA11;
        MatricesA12[i] = sqmatrixA12;
        MatricesA21[i] = sqmatrixA21;
        MatricesA22[i] = sqmatrixA22;

        sqmatrixM11 = new TSquareMatrix2D(sqstructureA);
        sqmatrixM12 = new TSquareMatrix2D(sqstructureA);
        sqmatrixM21 = new TSquareMatrix2D(sqstructureA);
        sqmatrixM22 = new TSquareMatrix2D(sqstructureA);

        MatricesM11[i] = sqmatrixM11;
        MatricesM12[i] = sqmatrixM12;
        MatricesM21[i] = sqmatrixM21;
        MatricesM22[i] = sqmatrixM22;
        break;

      case 4:
        matrixB1 = new TMatrix2D(structureB);
        matrixB2 = new TMatrix2D(structureB);
        matrixB1T = new TMatrix2D(structureBT);
        matrixB2T = new TMatrix2D(structureBT);

        MatricesB1[i] = matrixB1;
        MatricesB2[i] = matrixB2;
        MatricesB1T[i] = matrixB1T;
        MatricesB2T[i] = matrixB2T;

        sqmatrixA11 = new TSquareMatrix2D(sqstructureA);
        sqmatrixA12 = new TSquareMatrix2D(sqstructureA);
        sqmatrixA21 = new TSquareMatrix2D(sqstructureA);
        sqmatrixA22 = new TSquareMatrix2D(sqstructureA);

        MatricesA11[i] = sqmatrixA11;
        MatricesA12[i] = sqmatrixA12;
        MatricesA21[i] = sqmatrixA21;
        MatricesA22[i] = sqmatrixA22;

        MatricesF11[i] = new TSquareMatrix2D(sqstructureA);
        MatricesF22[i] = new TSquareMatrix2D(sqstructureA);

        sqmatrixM11 = new TSquareMatrix2D(sqstructureA);
        sqmatrixM12 = new TSquareMatrix2D(sqstructureA);
        sqmatrixM21 = new TSquareMatrix2D(sqstructureA);
        sqmatrixM22 = new TSquareMatrix2D(sqstructureA);

        MatricesM11[i] = sqmatrixM11;
        MatricesM12[i] = sqmatrixM12;
        MatricesM21[i] = sqmatrixM21;
        MatricesM22[i] = sqmatrixM22;

        MatricesG11[i] = new TSquareMatrix2D(SqGridStructure);
        MatricesG12[i] = new TSquareMatrix2D(SqGridStructure);
        MatricesG21[i] = new TSquareMatrix2D(SqGridStructure);
        MatricesG22[i] = new TSquareMatrix2D(SqGridStructure);


        MatricesS_A[i] = new TSquareMatrix1D(sqSurfStructure);
        MatricesS_M[i] = new TSquareMatrix1D(sqSurfStructure);
        MatricesK[i] = new TSquareMatrix1D(sqSurfStructure);
        break;
    }

    N_Unknowns = 2*N_U + N_P;

    OutPut("dof velocity : "<< setw(10) << 2* N_U << endl);
    OutPut("dof pressure : "<< setw(10) << N_P << endl);
    OutPut("dof all      : "<<  setw(10) << N_Unknowns  << endl);
    OutPut("dof freesurf surfactant : "<< setw(10) << N_S << endl);





    gridrhs = new double[2*N_G];
    gridsol = new double[2*N_G];
    d =  new double[2*N_G];
    tmpUP = new double[N_Unknowns];
    memset(gridsol, 0, 2*N_G*SizeOfDouble);
    memset(gridrhs, 0, 2*N_G*SizeOfDouble);
    memset(d, 0, 2*N_G*SizeOfDouble);
    memset(tmpUP, 0, N_Unknowns*SizeOfDouble);

    Rhsarray[2][i] = gridrhs;
    Sol[4][i] = gridsol;

   if(TDatabase::ParamDB->SOLVER_TYPE == GMG)
    {

      GridMGLevel = new TMGLevel2D(i, MatricesG11[i], gridrhs, gridsol,
                                 4, NULL);
      GridMG->AddLevel(GridMGLevel);
    }

  t1 = GetTime();

    // assemble matrix for grid moving
    fesp[0] = FeSps_Lev[2][i];
    SQMATRICES_GRID[0] = MatricesG11[i];
    SQMATRICES_GRID[0]->Reset();
    SQMATRICES_GRID[1] = MatricesG12[i];
    SQMATRICES_GRID[1]->Reset();
    SQMATRICES_GRID[2] = MatricesG21[i];
    SQMATRICES_GRID[2]->Reset();
    SQMATRICES_GRID[3] = MatricesG22[i];
    SQMATRICES_GRID[3]->Reset();
    aux = new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL,
                        0, NULL);

   Assemble2D(1, fesp,
             4, SQMATRICES_GRID,
             0, NULL,
             0, NULL, NULL,
             DiscreteFormGrid,
             GridBoundaryConditions,
             GridBoundValues,
             aux);
    delete aux;

 t2 = GetTime();

  if(TDatabase::ParamDB->SC_VERBOSE > 1)
  {
    cout << "Grid assembling done"<< endl;
    OutPut("Time for Grid assembling: " << t2-t1 << endl);
  }
    rhs = new double[N_Unknowns];
  //  RhsArray[i] = rhs;
    Rhsarray[0][i] = rhs;
    memset(rhs, 0, N_Unknowns*SizeOfDouble);
    sol = new double[N_Unknowns];
    oldsol = new double[N_Unknowns];
    refsol = new double[N_Unknowns];
    memset(sol, 0, N_Unknowns*SizeOfDouble);
    memset(oldsol, 0, N_Unknowns*SizeOfDouble);
    memset(refsol, 0, N_Unknowns*SizeOfDouble);
    if (i==mg_level-1)
       sol_timestep_m1 =  new double[N_Unknowns];

    B = new double [N_Unknowns];
    memset(B, 0, N_Unknowns*SizeOfDouble);

    Sol[0][i] = sol;
    Sol[1][i] = oldsol;
    Sol[2][i] = refsol;
    Rhsarray[1][i] = B;

    // ( A B' )
    // ( B 0  )

   if(TDatabase::ParamDB->SOLVER_TYPE == GMG)
    {

        // determine number of auxiliary arrays
       if((TDatabase::ParamDB->SC_STEP_LENGTH_CONTROL_ALL_SADDLE)
          || (TDatabase::ParamDB->SC_STEP_LENGTH_CONTROL_FINE_SADDLE))
         n_aux=4;
       else
         n_aux=2;
       switch(TDatabase::ParamDB->NSTYPE)
      {
        case 1:
          MGLevel = new TNSE_MGLevel1(i, sqmatrixM, matrixB1, matrixB2,
                                      structureBT, B, sol, n_aux,  alpha,
                                      velocity_space_code ,
                                      pressure_space_code,NULL);
        break;

        case 2:
          MGLevel = new TNSE_MGLevel2(i, sqmatrixM, matrixB1, matrixB2,
                                      matrixB1T, matrixB2T,
                                      B, sol, n_aux, alpha,
                                      velocity_space_code ,
                                      pressure_space_code,NULL);
        break;

        case 3:
          MGLevel = new TNSE_MGLevel3(i, sqmatrixM11, sqmatrixM12,
                                      sqmatrixM21, sqmatrixM22,
                                      matrixB1, matrixB2,
                                      structureBT,
                                      B, sol, n_aux, alpha,
                                      velocity_space_code ,
                                      pressure_space_code,NULL);
        break;

        case 4:
	  MGLevel = new TNSE_MGLevel4(i, sqmatrixM11, sqmatrixM12,
                                       sqmatrixM21, sqmatrixM22,
                                       matrixB1, matrixB2,
                                       matrixB1T, matrixB2T,
                                       B, sol, n_aux, alpha,
                                       velocity_space_code ,
                                       pressure_space_code,NULL);
        break;
      } // end switch(NSTYPE)

      TMGLevel[0][i] = MGLevel;
      MG->AddLevel(MGLevel);
    }

    u = new TFEVectFunct2D(velocity_space, UString,  UString,  sol, N_U, 2);
    u1 = u->GetComponent(0);
    u2 = u->GetComponent(1);
    p = new TFEFunction2D(pressure_space, PString,  PString,  sol+2*N_U, N_P);

    UPArrays[0][i] = u1;
    UPArrays[1][i] = u2;
    UPArrays[2][i] = p;
    VeloVect[0][i] = u;

//     u1->Interpolate(InitialU1);
//     u2->Interpolate(InitialU2);
//     p->Interpolate(InitialP);

    velocity = new double[2*N_G];
    memset(velocity, 0, 2*N_G*SizeOfDouble);
    Sol[3][i] = velocity;
    GridVelocity = new TFEVectFunct2D(FeSps_Lev[2][i], veloString, gridString,
                                      velocity, N_G, 2);

    VeloVect[1][i] = GridVelocity;

    GridXDot = GridVelocity->GetComponent(0);
    GridYDot = GridVelocity->GetComponent(1);
    UPArrays[3][i] = GridXDot;
    UPArrays[4][i] = GridYDot;

//     GridXDot->Interpolate(GridU1);
//     GridYDot->Interpolate(GridU2);


// surfactant

    srhs = new double[N_S];
    memset(srhs, 0, N_S*SizeOfDouble);
    Rhsarray[3][i] = srhs;

    ssol = new double[N_S];
    Sol[6][i] = ssol;
    memset(ssol, 0, N_S*SizeOfDouble);
    SurfSurfactant = new TFEFunction1D(SurfFeSps_Lev[0][i], SurfString, SurfString,
                                       ssol, N_S);
//     SurfSurfactant->Interpolate(InitialS);
    SArrays[0][i] = SurfSurfactant;

// //  surfacetant in entire 2D domain
    surfact = new double[N_SO];
    Sol[7][i] = surfact;
    memset(surfact, 0, N_SO*SizeOfDouble);
    UPArrays[5][i] = new TFEFunction2D(FeSps_Lev[5][i], SurfString, SurfString, 
                                       surfact, N_SO);

    UPArrays[5][i]->Interpolate(InitialS);
    MapDomainToSurf(UPArrays[5][i], SArrays[0][i], N_List[0][i], N_List[1][i]);
    memset(surfact, 0, N_SO*SizeOfDouble);

    surfact_velo = new double[2*N_SO];
    Sol[8][i] = surfact_velo;
    memset(surfact_velo, 0, 2*N_SO*SizeOfDouble);
    UPArrays[6][i] = new TFEFunction2D(FeSps_Lev[5][i], PsiString, PsiString, 
                                       surfact_velo, N_SO);
    UPArrays[7][i] = new TFEFunction2D(FeSps_Lev[5][i], PsiString, PsiString,
                                       surfact_velo+N_SO, N_SO);

    } // endfor i

  TMG[0] = MG;

  tmpUP = new double[N_Unknowns];

  refpos = new double[2*N_G];
  auxpos = new double[2*N_G];
  pos = new double[2*N_G];
  tmp = new double[2*N_G];

  left = -TDatabase::ParamDB->P2;
  bottom = 0.0;
  right = TDatabase::ParamDB->P2;
  top = TDatabase::ParamDB->P2;

 // ******************************************************************

  RefGridPos = new TFEVectFunct2D(FeSps_Lev[2][mg_level-1], refposString, gridString,
                                  refpos, N_G, 2);
  AuxGridPos = new TFEVectFunct2D(FeSps_Lev[2][mg_level-1], auxposString, gridString,
                                  auxpos, N_G, 2);
  GridPos = new TFEVectFunct2D(FeSps_Lev[2][mg_level-1], posString, gridString,
                               pos, N_G, 2);

  // ******************************************************************

  RefGridPos->GridToData();

  RefGridPos->DataToGrid();
  GridPos->GridToData();
  RefGridPos->DataToGrid();

  Entries[0] = MatricesG11[mg_level-1]->GetEntries();
  Entries[1] = MatricesG12[mg_level-1]->GetEntries();
  Entries[2] = MatricesG21[mg_level-1]->GetEntries();
  Entries[3] = MatricesG22[mg_level-1]->GetEntries();

  GridKCol = SqGridStructure->GetKCol();
  GridRowPtr = SqGridStructure->GetRowPtr();

  memset(Entries[1] + GridRowPtr[N_G-N_BoundaryNodes], 0,
         N_BoundaryNodes*SizeOfDouble);
  memset(Entries[2] + GridRowPtr[N_G-N_BoundaryNodes], 0,
         N_BoundaryNodes*SizeOfDouble);

  GridPos->DataToGrid();

// assembling surfactant matrices
  N_SquareMatrices = 2;
  SQMATRICES_SURF[0] = MatricesS_A[mg_level -1];
  SQMATRICES_SURF[0]->Reset();
  SQMATRICES_SURF[1] = MatricesS_M[mg_level -1];
  SQMATRICES_SURF[1]->Reset();

  N_FESpaces = 3;
  fesp[0] = FeSps_Lev[0][mg_level -1];  // velocity space
  fesp[1] = FeSps_Lev[2][mg_level -1];  // grid vello space
  fesp[2] = FeSps_Lev[5][mg_level -1];  //  surfact_space in entire domain

  N_FESpaces_low=1;
  fesp_low[0] = SurfFeSps_Lev[0][mg_level -1];

  fefct[0] = UPArrays[0][mg_level -1];  // ur
  fefct[1] = UPArrays[1][mg_level -1];  //uz
  fefct[2] = UPArrays[3][mg_level -1];  // wr
  fefct[3] = UPArrays[4][mg_level -1];  //wz

  fefct_low[0] = SArrays[0][mg_level -1];

  N_Rhs =1;
  ferhs_low[0] = SurfFeSps_Lev[0][mg_level -1];
  rhs =  Rhsarray[3][mg_level -1];
  RHSs[0] = rhs;

  AssembleSurf1D(N_FESpaces, fesp, fefct, N_FESpaces_low, fesp_low,
                 N_SquareMatrices, SQMATRICES_SURF, N_Rhs, RHSs, 
                 ferhs_low, N_List[0][mg_level-1], 
                 N_List[1][mg_level-1]);

// 
    // prepare output (maxn_fespaces,  maxn_scalar,  maxn_vect, maxn_parameters, domain)
    // prepare output
      Output = new TOutput2D(3, 3, 2, 2, Domain);
      Output->AddFEVectFunct(u);
//       Output->AddFEVectFunct(GridVelocity);
//       Output->AddFEFunction(p);
      Output->AddFEFunction(UPArrays[5][mg_level-1]);
      os.seekp(std::ios::beg);
      Output->AddParameter(TDatabase::TimeDB->CURRENTTIME,os.str().c_str());

      vorticity = new double[N_Vort];
      memset(vorticity, 0, N_Vort*SizeOfDouble);

      div = new double[N_Vort];
      memset(div, 0, N_Vort*SizeOfDouble);
//      ComputeVorticityDivergence(velocity_space, u1, u2, vorticity_space, vorticity, div);
      Divergence = new TFEFunction2D(vorticity_space, DivergenceString, DivergenceString, div, N_Vort);
//       Output->AddFEFunction(Divergence);


  t = TDatabase::TimeDB->CURRENTTIME-TDatabase::TimeDB->STARTTIME;

  defect = new double[N_Unknowns];
  startsol = new double[N_Unknowns];
  frac_step_sol = new double[N_Unknowns];
  oldrhs =  new double[N_Unknowns];

  N_Active = velocity_space->GetActiveBound();

  SB = new double [N_S];
  memset(SB, 0, N_S*SizeOfDouble);
  sdefect = new double[N_S];
  memset(sdefect, 0, N_S*SizeOfDouble);
  ssol_old = new double[N_S];
  memcpy(ssol_old, ssol, N_S*SizeOfDouble);
  N_SurfActive = SurfFeSps_Lev[0][mg_level -1]->GetActiveBound();

  solver_time = 0.0;
  N_LinIter = 0;

  gamma = 0;
  Surfgamma = 0;
  m = 0;
  N_SubSteps = GetN_SubSteps();
  oldtau = 1;
  end_time = TDatabase::TimeDB->ENDTIME;
  if (TDatabase::TimeDB->TIMESTEPLENGTH_CONTROL)
    time_discs = 2;
  else
    time_discs = 1;
  // initialize solver
  if (TDatabase::ParamDB->SOLVER_TYPE==GMG)
  {
    switch (TDatabase::ParamDB->SC_SOLVER_SADDLE)
    {
      case 11:
        zerostart = 1;
        break;
      case 16:
        zerostart = 0;
        break;
    }

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
          }
        else
        {
          itmethod_sol = sol;
          itmethod_rhs = rhs;
        }
        break;
      case 16:
        itmethod = new TFgmresIte(MatVect, Defect, prec,
                                  0, N_Unknowns, 0);
        if (TDatabase::ParamDB->SC_PRECONDITIONER_SADDLE == 5)
        {
          itmethod_sol = new double[N_Unknowns];
          itmethod_rhs = new double[N_Unknowns];
        }
        else
        {
          itmethod_sol = sol;
          itmethod_rhs = rhs;
        }
        break;
      default:
        OutPut("Unknown solver !!!" << endl);
        exit(4711);
    }
  }

  coll->GetHminHmax(&hmin,&hmax);
  OutPut("h_min : " << hmin << " h_max : " << hmax << endl);
  if (fabs(hmin-hmax)<1e-6)
  {
    OutPut("delta " <<  TDatabase::ParamDB->FILTER_WIDTH_CONSTANT *
           pow(hmin,TDatabase::ParamDB->FILTER_WIDTH_POWER) << endl);
  }
  else
  {
    if (TDatabase::ParamDB->FILTER_WIDTH_POWER==0)
    {
      OutPut("delta " <<  TDatabase::ParamDB->FILTER_WIDTH_CONSTANT << endl);
    }
    else
    {
      OutPut("delta is non--constant" << endl);
    }
  }

  // copy sol for extrapolation after time step
  memcpy(sol_timestep_m1,sol,N_Unknowns*SizeOfDouble);

  if(TDatabase::ParamDB->WRITE_PS)
     {
        // write grid into an Postscript file
      os.seekp(std::ios::beg);
      if(img<10)  os << "final.0000" << img << ".ps" << ends;
       else if(img<100)  os << "final.000" << img << ".ps" << ends;
       else if(img<1000)  os << "final.00" << img << ".ps" << ends;
       else if(img<10000)  os << "final.00" << img << ".ps" << ends;
       else  os << "final" << img << ".ps" << ends;
      Domain->PS(os.str().c_str(),It_Finest,0);
     // if(!(TDatabase::ParamDB->WRITE_VTK)) img++;
    }


     if(TDatabase::ParamDB->WRITE_VTK)
       {
//        MapDomainToSurf(U1Array[mg_level-1],SU1Array[mg_level-1], Cell_No_array[mg_level-1], Joint_No_array[mg_level-1]);
//        MapSurfToDomain(SU1Array[mg_level-1], Surfactant, Cell_No_array[mg_level-1], Joint_No_array[mg_level-1]);
       MapSurfToDomain(SArrays[0][mg_level-1], UPArrays[5][mg_level-1], 
                       N_List[0][mg_level-1], N_List[1][mg_level-1]);
        os.seekp(std::ios::beg);
        if(img<10) os << VtkBaseName<<".0000"<<img<<".vtk" << ends;
         else if(img<100) os << VtkBaseName<<".000"<<img<<".vtk" << ends;
         else if(img<1000) os << VtkBaseName<<".00"<<img<<".vtk" << ends;
         else if(img<10000) os << VtkBaseName<<".0"<<img<<".vtk" << ends;
         else  os << VtkBaseName<<"."<<img<<".vtk" << ends;
        Output->WriteVtk(os.str().c_str());
        img++;
       }
// exit(0);
   Getboundcellangle(coll, FeSps_Lev[2][mg_level-1], Angle);
 // Getcellangle(FeSps_Lev[0][mg_level-1], Angle);
// remeshing the domain

 N_array[0] = N_Uarray;
 N_array[1] = N_Parray;
 N_array[2] = N_Garray;
//  N_array[3] = GridKCol;
//  N_array[4] = GridRowPtr;
//  N_array[5] =  Edge_No;

 SqMat[0]  = MatricesA11;
 SqMat[1]  = MatricesA12;
 SqMat[2]  = MatricesA21;
 SqMat[3]  = MatricesA22;
 SqMat[4]  = MatricesM11;
 SqMat[5]  = MatricesM12;
 SqMat[6]  = MatricesM21;
 SqMat[7]  = MatricesM22;
 SqMat[8]  = MatricesG11;
 SqMat[9]  = MatricesG12;
 SqMat[10] = MatricesG21;
 SqMat[11] = MatricesG22;
 SqMat[12] = MatricesA;
 SqMat[13] = MatricesM;
 SqMat[14] = MatricesF11;
 SqMat[15] = MatricesF22;

 Mat[0] = MatricesB1;
 Mat[1] = MatricesB2;
 Mat[2] = MatricesB1T;
 Mat[3] = MatricesB2T;

 SqMat_low[0] = MatricesS_A;
 SqMat_low[1] = MatricesS_M;
 SqMat_low[2] = MatricesK;

 Sol[5][0] = sol_timestep_m1;

 N_MovVert[0] = N_SlipBound_Vert;
 N_MovVert[1] = N_FreeBound_Vert;
 N_MovVert[2] = N_FreeGauss_Vert;

 if(TDatabase::ParamDB->P5 > 0)
       cout << "Moving domain using (U.n)n " << endl;
 else  cout << "Moving domain using (dt.U) " << endl;


//     InitVolume = CurrVolume = Volume(FeSps_Lev[0][mg_level-1]);
    Get_KE(VeloVect[0][mg_level-1], Params);
    InitVolume = CurrVolume = Params[0];
    MovBoundVert[0][0]->GetCoords(Lx, Ly);
    MovBoundVert[0][N_MovVert[0]-1]->GetCoords(Rx, Ry);
    OutPut(setw(20)<<"T, Wett Len d : " << TDatabase::TimeDB->CURRENTTIME
                   <<"   "<< Rx-Lx<<endl);
    OutPut(setw(20)<<"T, Volume : " << TDatabase::TimeDB->CURRENTTIME
                   <<"   "<< CurrVolume<<endl);
    OutPut(setw(20)<<"T, Volume Diff : "<< TDatabase::TimeDB->CURRENTTIME
                   <<"   "<< CurrVolume - InitVolume << endl);

//     min = 1e5;
//     tx  = 0.;
// 
//     for(k=0;k<N_MovVert[1];k++) // no need to set end vertices again
//      {
//       MovBoundVert[1][k]->GetCoords(x1, y1);
// //         cout << x1 << " " <<  y1<< endl;
//       if(fabs(y1)<min)
//        {
//         tx = x1;
//         min = fabs(y1);
//        }
//      }

//      MovBoundVert[0][0]->GetCoords(x1, y1);


//      OutPut(setw(25)<<"Time, Right tip, Top tip: "<< TDatabase::TimeDB->CURRENTTIME
//                    <<"   "<< tx <<"   "<< y1 << endl);

//   MovBoundVert[1][N_MovVert[1]-1]->GetCoords(x1, y1);
//   y_top[0] = y_top[1] = y_top[2] = y1;
//   t_top[0] = t_top[1] = t_top[2] = 0.;


        PrintSurfSurfactant(Free_Cells, MovBoundVert[1], N_array[3], 
                            N_MovVert[1], UPArrays[5][mg_level-1], N_BData);

Terrors[0] = 0;
Terrors[1] = 0;
T_inftyL2 = -1e10;
T_inftyH1 = -1e10;
errors[3] = 0;

     MapSurfToDomain(SArrays[0][mg_level-1], UPArrays[5][mg_level-1], 
                      N_List[0][mg_level-1], N_List[1][mg_level-1]);

      GetSurfactMass(UPArrays[5][mg_level-1], SArrays[0][mg_level-1], 
                     N_List[0][mg_level-1], N_List[1][mg_level-1], Surf_Mass, max_r);

      l_2_l_2u =  4.*Pi - 2.*Pi*Surf_Mass[0]; // not really L2

//       if(T_inftyL2<l_2_l_2u)
//         {
//          T_inftyL2=l_2_l_2u;
//          T_inftyL2_time=TDatabase::TimeDB->CURRENTTIME;
//         }

        olderror = l_2_l_2u;
        OutPut(TDatabase::TimeDB->CURRENTTIME <<  " L2(0,T;L2) " << sqrt(errors[3]) << " ");

        OutPut("time: " << TDatabase::TimeDB->CURRENTTIME);
        OutPut(" Mass Error: " <<  l_2_l_2u << endl);
        OutPut(T_inftyL2_time <<  " L2error_Max " << T_inftyL2 << " ");
	


//======================================================================
// start of time cycle
//======================================================================
  while(TDatabase::TimeDB->CURRENTTIME< end_time)
  {
    // time cycle
    m++;
    TDatabase::TimeDB->INTERNAL_STARTTIME = TDatabase::TimeDB->CURRENTTIME;

    for (methods=0;methods<time_discs;methods++)
    {
      if (time_discs==2)
      {
        if (methods==0) // fractional-step-theta-scheme
        {
          TDatabase::TimeDB->TIME_DISC = 3;
          memcpy(startsol,Sol[0][mg_level-1],N_Unknowns*SizeOfDouble); // save start Sol[0][mg_level-1] (nec. for gl00)
          memcpy(oldrhs,rhs,N_Unknowns*SizeOfDouble); // save rhs
        }
        else           // crank nicolson scheme
        {              // take solution of first scheme as initial iterate
          TDatabase::TimeDB->TIME_DISC = 2;
          TDatabase::TimeDB->CURRENTTIME = TDatabase::TimeDB->INTERNAL_STARTTIME;
          memcpy(frac_step_sol,Sol[0][mg_level-1],N_Unknowns*SizeOfDouble); // save solution of fract.step
          memcpy(Sol[0][mg_level-1],startsol,N_Unknowns*SizeOfDouble); // get former startsol
          memcpy(rhs,oldrhs,N_Unknowns*SizeOfDouble); // get old rhs
        }
        N_SubSteps = GetN_SubSteps();
      }

      for(l=0;l<N_SubSteps;l++)      // sub steps of fractional step theta
      {
        if (!very_first_time)
        {
          SetTimeDiscParameters();
        }

        if (m==1)
        {
          OutPut("Theta1: " << TDatabase::TimeDB->THETA1<< endl);
          OutPut("Theta2: " << TDatabase::TimeDB->THETA2<< endl);
          OutPut("Theta3: " << TDatabase::TimeDB->THETA3<< endl);
          OutPut("Theta4: " << TDatabase::TimeDB->THETA4<< endl);
        }

        tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
        if (!very_first_time)
          TDatabase::TimeDB->CURRENTTIME += tau;

        if (very_first_time)
            oldtau=tau;


//======================================================================
// solve the surfactant equation begin
//======================================================================

   GridTOVelocity(VeloVect[0][mg_level -1], TDatabase::TimeDB->CURRENTTIME, tau, max_r);
   MoveGrid(GridPos, AuxGridPos, FeSps_Lev[5][mg_level -1], TDatabase::TimeDB->CURRENTTIME, tau );

//      if(TDatabase::ParamDB->WRITE_VTK)
//        {
//         MapSurfToDomain(SArrays[0][mg_level-1], UPArrays[5][mg_level-1], 
//                        N_List[0][mg_level-1], N_List[1][mg_level-1]);
// 
//         os.seekp(std::ios::beg);
//         if(img<10) os << VtkBaseName<<".0000"<<img<<".vtk" << ends;
//          else if(img<100) os << VtkBaseName<<".000"<<img<<".vtk" << ends;
//          else if(img<1000) os << VtkBaseName<<".00"<<img<<".vtk" << ends;
//          else if(img<10000) os << VtkBaseName<<".0"<<img<<".vtk" << ends;
//          else  os << VtkBaseName<<"."<<img<<".vtk" << ends;
//         Output->WriteVtk(os.str().c_str());
//         img++;
//        } 
// exit(0);

 // working array for rhs is SB, initialize SB
   memset(SB, 0, N_S*SizeOfDouble);
   // old srhs multiplied with current subtime step and theta3 on SB
   //   N_SurfActive= N_S on fespace on a 1D surface
   Daxpy(N_S, tau*TDatabase::TimeDB->THETA3, Rhsarray[3][mg_level -1], SB);


// assembling surfactant matrices
  N_SquareMatrices = 2;
  SQMATRICES_SURF[0] = MatricesS_A[mg_level -1];
  SQMATRICES_SURF[0]->Reset();
  SQMATRICES_SURF[1] = MatricesS_M[mg_level -1];
  SQMATRICES_SURF[1]->Reset();

  N_FESpaces = 3;
  fesp[0] = FeSps_Lev[0][mg_level -1];  // velocity space
  fesp[1] = FeSps_Lev[2][mg_level -1];  // grid vello space
  fesp[2] = FeSps_Lev[5][mg_level -1];  //  surfact_space in entire domain

  N_FESpaces_low=1;
  fesp_low[0] = SurfFeSps_Lev[0][mg_level -1];

  fefct[0] = UPArrays[0][mg_level -1];  // ur
  fefct[1] = UPArrays[1][mg_level -1];  //uz
  fefct[2] = UPArrays[3][mg_level -1];  // wr
  fefct[3] = UPArrays[4][mg_level -1];  //wz

  fefct_low[0] = SArrays[0][mg_level -1];

  N_Rhs =1;
  ferhs_low[0] = SurfFeSps_Lev[0][mg_level -1];

  memset(Rhsarray[3][mg_level -1], 0, N_S*SizeOfDouble); 
  srhs = Rhsarray[3][mg_level -1];
  RHSs[0] =  srhs;

  AssembleSurf1D(N_FESpaces, fesp, fefct, N_FESpaces_low, fesp_low,
                 N_SquareMatrices, SQMATRICES_SURF, N_Rhs, RHSs, 
                 ferhs_low, N_List[0][mg_level-1], 
                 N_List[1][mg_level-1]);


       if (very_first_time==1)
        {
          very_first_time=0;
          l--;
          continue;
        }

  // add rhs from current sub time step to rhs array B
    Daxpy(N_S, tau*TDatabase::TimeDB->THETA4, srhs, SB);


   memset(sdefect, 0, N_S*SizeOfDouble);

   for(i=0;i<mg_level;i++)
     MatAdd(MatricesS_M[i], MatricesS_A[i], 
            -tau*TDatabase::TimeDB->THETA2);

   Surfgamma = -tau*TDatabase::TimeDB->THETA2;

   MatVectActive(MatricesS_M[mg_level-1], ssol, sdefect);
   Daxpy(N_S, 1, sdefect, SB);

    // copy Dirichlet values from rhs

   //=====================================================================
   // assembling of system matrix
   //========================================================================
   for(i=0;i<mg_level;i++)
      MatAdd(MatricesS_M[i], MatricesS_A[i],
             -Surfgamma + tau*TDatabase::TimeDB->THETA1);

   Surfgamma = tau*TDatabase::TimeDB->THETA1;

   memset(sdefect, 0, N_S*SizeOfDouble);
        // compute defect
// // //    Defect(MatricesS_M[mg_level-1],NULL,ssol,SB,sdefect);
//    residual =  Ddot(N_S, sdefect, sdefect);
//    residual = sqrt(residual);
// 
//    OutPut("initial residual ");
//    OutPut(setw(14) << residual << endl);

   DirectSolver(MatricesS_M[mg_level-1], SB, ssol);
// exit(0);
//======================================================================
// solve the surfactant equation end
//======================================================================

    } // endfor l (sub steps of fractional step theta)
   } // endfor two time discs of adaptive time step control


  if(((m % TDatabase::TimeDB->STEPS_PER_IMAGE) == 0) || m==1)
    {
   if(TDatabase::ParamDB->WRITE_PS)
     {
        // write grid into an Postscript file
      os.seekp(std::ios::beg);
       if(img<10)  os << "final.0000" << img << ".ps" << ends;
       else if(img<100)  os << "final.000" << img << ".ps" << ends;
       else if(img<1000)  os << "final.00" << img << ".ps" << ends;
       else if(img<10000)  os << "final.00" << img << ".ps" << ends;
       else  os << "final" << img << ".ps" << ends;
      Domain->PS(os.str().c_str(),It_Finest,0);
     }

    if(TDatabase::ParamDB->WRITE_GRAPE)
      {
        os.seekp(std::ios::beg);
        os << GrapeBaseName << img << ".dat" << ends;
        Output->WriteGrape(os.str().c_str());
        N_GRAPE_images++;
       }

    if(TDatabase::ParamDB->WRITE_GNU)
       {

        os.seekp(std::ios::beg);
        os << GnuBaseName << img << ".gnu" << ends;
        Output->WriteGnuplot(os.str().c_str());

        }

     if(TDatabase::ParamDB->WRITE_VTK)
       {
        MapSurfToDomain(SArrays[0][mg_level-1], UPArrays[5][mg_level-1], 
                       N_List[0][mg_level-1], N_List[1][mg_level-1]);

        os.seekp(std::ios::beg);
        if(img<10) os << VtkBaseName<<".0000"<<img<<".vtk" << ends;
         else if(img<100) os << VtkBaseName<<".000"<<img<<".vtk" << ends;
         else if(img<1000) os << VtkBaseName<<".00"<<img<<".vtk" << ends;
         else if(img<10000) os << VtkBaseName<<".0"<<img<<".vtk" << ends;
         else  os << VtkBaseName<<"."<<img<<".vtk" << ends;
        Output->WriteVtk(os.str().c_str());
        img++;
       } 
        PrintSurfSurfactant(Free_Cells, MovBoundVert[1], N_array[3], 
                            N_MovVert[1], UPArrays[5][mg_level-1], N_BData);
    } 

//   min = 1e5;
//   tx  = 0.;
//   m4   = 0;
//   if((m % 20) == 0)

/*
      UPArrays[5][mg_level-1]->Interpolate(InitialS);
      MapDomainToSurf(UPArrays[5][mg_level-1], SArrays[0][mg_level-1],
                    N_List[0][mg_level-1], N_List[1][mg_level-1]);
      memset(surfact, 0, N_SO*SizeOfDouble);
*/
//    if((m % 8) == 0)
//     {
//       MapSurfToDomain(SArrays[0][mg_level-1], UPArrays[5][mg_level-1], 
//                       N_List[0][mg_level-1], N_List[1][mg_level-1]);
// 
//       GetSurfErrors(UPArrays[5][mg_level-1], SArrays[0][mg_level-1], 
//                     errors, N_List[0][mg_level-1], N_List[1][mg_level-1]);
// 
//       OutPut( "L2: " << errors[0] << endl);
//       OutPut( "H1-semi: " << errors[1] << endl);
// 
//       Terrors[0] += (errors[0]*errors[0]);
//       Terrors[1] += (errors[1]*errors[1]);
// 
//       if(T_inftyL2<errors[0]) T_inftyL2=errors[0];
//       if(T_inftyH1<errors[1]) T_inftyH1=errors[1];
//     }

//   if((m % 5) == 0)
//    {
// 
// 
//      MovBoundVert[0][0]->GetCoords(x1, y1);
// 
// 
//      OutPut(setw(25)<<"Time, Right tip, Top tip: "<< TDatabase::TimeDB->CURRENTTIME
//                    <<"   "<< tx <<"   "<< y1 << endl);
// 
//     Get_KE(VeloVect[0][mg_level-1], Params);
// 
//     OutPut(setw(25)<<"Time, Volume : " << TDatabase::TimeDB->CURRENTTIME
//                    <<"   "<< Params[0]<<endl);
//     OutPut(setw(25)<<"Time, Volume Diff : "<< TDatabase::TimeDB->CURRENTTIME
//                    <<"   "<< Params[0] - InitVolume << endl)
//     OutPut(TDatabase::TimeDB->CURRENTTIME<<
//           " x_center_mass: "<< Params[2] <<
//           " y_center_mass: "<< Params[3] << " U1_Rise: "<< Params[4] <<
//           " U2_Rise: "<< Params[5]<< endl);
//   OutPut(TDatabase::TimeDB->CURRENTTIME<<
//          " Kinetic Energy: "<< Params[1] << endl);
// 
//   }
// exit(0);
     MapSurfToDomain(SArrays[0][mg_level-1], UPArrays[5][mg_level-1], 
                      N_List[0][mg_level-1], N_List[1][mg_level-1]);

      GetSurfactMass(UPArrays[5][mg_level-1], SArrays[0][mg_level-1], 
                     N_List[0][mg_level-1], N_List[1][mg_level-1], Surf_Mass, max_r);

//       l_2_l_2u =  2.*Pi - 2.*Pi*Surf_Mass[0]; // not really L2
      l_2_l_2u =  4.*Pi - 2.*Pi*Surf_Mass[0]; // not really L2
      
      if(T_inftyL2<l_2_l_2u)
        {
         T_inftyL2=l_2_l_2u;
         T_inftyL2_time=TDatabase::TimeDB->CURRENTTIME;
        }

        errors[3] += (l_2_l_2u*l_2_l_2u +olderror * olderror)*TDatabase::TimeDB->TIMESTEPLENGTH/2.0;
        olderror = l_2_l_2u;
        OutPut(TDatabase::TimeDB->CURRENTTIME <<  " L2(0,T;L2) " << sqrt(errors[3]) << " ");
        OutPut(" Mass Error: " <<  l_2_l_2u << " ");
        OutPut(T_inftyL2_time <<  " L2error_Max " << T_inftyL2 << endl);


// exit(0);
  } // while
//
//======================================================================
// end of time cycle
//======================================================================

   if(TDatabase::ParamDB->WRITE_GRAPE)
    {
      os.seekp(std::ios::beg);
      os << GrapeBaseName <<  "end." << m << ".dat" << ends;
      Output->WriteGrape(os.str().c_str());
    }

    if(TDatabase::ParamDB->WRITE_GNU)
    {
      os.seekp(std::ios::beg);
      os << GnuBaseName << m << ".gnu" << ends;
      Output->WriteGnuplot(os.str().c_str());
    }

     if(TDatabase::ParamDB->WRITE_VTK)
       {
        os.seekp(std::ios::beg);
        if(img<10) os << VtkBaseName<<".0000"<<img<<".vtk" << ends;
         else if(img<100) os << VtkBaseName<<".000"<<img<<".vtk" << ends;
          else if(img<1000) os << VtkBaseName<<".00"<<img<<".vtk" << ends;
           else if(img<10000) os << VtkBaseName<<".0"<<img<<".vtk" << ends;
            else  os << VtkBaseName<<"."<<img<<".vtk" << ends;
        Output->WriteVtk(os.str().c_str());
       }

/*    CurrVolume = Volume(FeSps_Lev[0][mg_level-1]);
    MovBoundVert[0][0]->GetCoords(Lx, Ly); // left wetting points
    MovBoundVert[0][N_MovVert[0]-1]->GetCoords(Rx, Ry);  // right wetting points
    OutPut(setw(25)<<"Total No of time steps : " << m <<endl);
    OutPut(setw(25)<<"No of time Remeshed : " << N_Remesh <<endl);
    OutPut(setw(25)<<"No of time ReConstructed : "<<N_ReConstruct<< endl);
    OutPut(setw(25)<<"Time, Wett Len d : " << TDatabase::TimeDB->CURRENTTIME
                   <<"   "<< Rx-Lx<<endl);
    OutPut(setw(25)<<"Time, Volume : " << TDatabase::TimeDB->CURRENTTIME
                   <<"   "<< CurrVolume<<endl);
    OutPut(setw(25)<<"Time, Volume Diff : "<< TDatabase::TimeDB->CURRENTTIME
                   <<"   "<< CurrVolume - InitVolume << endl);

  t4 =  GetTime();
  total_time += t4 - t3;
  OutPut("total time for remeshing: " << Remesh_Time <<" Sec"<< endl);
  OutPut("total running time: " << total_time<<" Sec" << endl);*/

/*
      OutPut( "T_inftyL2: " << T_inftyL2 << endl);
      OutPut( "T_L2L2: " << sqrt(Terrors[0]) << endl);
      OutPut( "T_inftyH1: " << T_inftyH1 << endl);
      OutPut( "T_L2H1: " << sqrt(Terrors[1]) << endl);
*/

  CloseFiles();
  return 0;
}

