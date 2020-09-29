// =======================================================================
//
// Purpose:     main program for axisymmetric droplet oscillating with surfactant transport
//
// Author:     S.Ganesan
//                  22.9.2007
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
//       errors[0] +=Mult*(Exact_Surf[0]-U)*(Exact_Surf[0]-U);
      errors[0] +=U*Mult;
      errors[1] +=Mult;
     } //  for(k=0;k<N_LinePoints;k++)

//      if(h_K_min>h_K)
//          h_K_min = h_K;
//      if(h_K_max<h_K)
//          h_K_max = h_K;

  } // for(i=0;i<N
//     OutPut("h_K_min and h_K_max of free surface: "<< h_K_min << " " << h_K_max<<endl; );

   OutPut( "Time, Surfactant Mass " <<TDatabase::TimeDB->CURRENTTIME<< " " <<2.*Pi*errors[0]<< " "<<endl);
   OutPut( "Time, Surface area " <<TDatabase::TimeDB->CURRENTTIME<< " " <<2.*Pi*errors[1]<< " "<<endl);

// exit(0);
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
    dat <<  x1<< "  " << y1<<"  "<<ArcLength*char_L << "  " << T_val[0] <<endl;

    x1=x2; y1=y2;
   }


Surfact->FindGradient( x1,  y1, T_val);
dat <<  x1<< "  " << y1<<"  "<<ArcLength*char_L << "  " << T_val[0] <<endl;

      dat.close();
      cout << endl;
      cout << "Surfactant data wrote into file " << endl;
 N_BData++;

}


void Remesh2D_axial3D_Surf(TDomain *Domain, TDomain *SurfDomain, TBaseCell **Free_Cells,
              TVertex ***MovBoundVert, int *N_MovVert, TFESpace2D ***FeSps_Lev,
              TFESpace1D ***SurfFeSps_Lev, int **N_array, int ***N_List, 
              TSquareMatrix2D ***SqMat, TMatrix2D ***Mat, TSquareMatrix1D ***SqMat_low, 
              double ***Sol, double ***Rharray, BoundCondFunct2D *BoundCondition, 
              BoundCondFunct2D *GridBoundCondition ,
              TNSE_MultiGrid **TMG, TFEVectFunct2D ***VeloVect,
              TFEFunction2D  ***UPArrays, TFEFunction1D  ***SArrays, 
              TBoundEdge ***Slip_Joint,  TSquareStructure2D ***SqrStruct, 
              TStructure2D ***Struct, TSquareStructure1D ***SqrStruct_low, FE1D *FE1D_List)
{
  int i, j, k, l, N_RootCells, Old_N_RootCells, N_Cells, Old_N_Cells, N_Joints;
  int N_Vertices, i3, N_G;
  int m, m1, m2, m3, m4, N_SlipJoints, N_FreeJoints, N_Velo_IsoPoints;
  int ORDER, VSP, N_OldU, N_OldP, *Edge_No;
  int LEVELS, mg_level, mg_type, order, N_Active, N_P, N_U, N_V, N_Vort;
  int N_BoundaryNodes, N_Unknowns, n_aux, mixing_layer=0;
  int velocity_space_code, pressure_space_code;
  double left, right, top, bottom, rad1, alpha, Lx =-100, Ly =-100, Rx=-100, Ry=-100;
//  double NLx, NLy, NRx, NRy;
  int *PointNeighb, maxEpV = 0, a, b, Neighb_tmp, Neib[2];
  int CurrNeib, len1, len2, CurrComp, comp, N_Interf_Vertices ;
  int In_Index, CurrVertex, CurrJoint, CurrInterfCell, ID, bct, bct1;
  double  T_a, T_b, temp, temp3, T, *Coordinates, x, y=0, tx, ty, TX[2], TY[2], x1, x2;
  int *Triangles;
  int ChangeBound1 = 0, ChangeBound2 = 0;
  double *S_BX, *S_BY, d, t, hE, y1, y2, *ValuesVX, *oldssol;
  TNSE_MGLevel *MGLevel;
  int  *BeginIndex, *JointDOF, N_DOF, *GlobalNumbers, *OldGlobalNumbers;
  FE2D FeId;
  TFEDesc2D *FeDesc;


  TIsoBoundEdge *Free_J;
  TBoundPart *BoundPart;
  TBdLine *UpdateSlipBound;
  TBaseCell **CellTree, **Old_CellTree, *cell, *Me, **del_Cell_Array, **Surf_CellTree; //**Cell_Array,
  TCollection *Coll, *Old_Coll, *SOld_Coll,*Space_Old_Coll, *mortarcoll = NULL, *Surf_Coll;
  TVertex **VertexDel, **NewVertices, *temp_Mov, *temp_Mov2;
  double *New_Sol;
  TJoint *Joint;
  TBoundEdge *tempSlip_Joint;
  TFEVectFunct2D *New_u;
  TFEFunction2D *New_p, *New_s;

  char *PsBaseName, *GrapeBaseName, *GnuBaseName, *ReadGrapeBaseName;
  char *VtkBaseName;

  TFEFunction2D *GridXDot, *GridYDot;
  int N_Parameters=1;
  double Parameters[2];

  TFESpace2D *Velo_Space, *Pres_Space, *OldSurfact_space;
//   TFESpace1D *Old_Surfspace;
  std::ostringstream opts;
  std::ostringstream os;

  os << " ";
  opts << " ";
  #define AMG 0
  #define GMG 1

  struct triangulateio In, Out;

  boolean AllowEdgeRef = (boolean) TDatabase::ParamDB->MESHGEN_ALLOW_EDGE_REF;
  alpha = TDatabase::ParamDB->SC_GMG_DAMP_FACTOR_SADDLE;
  LEVELS = TDatabase::ParamDB->LEVELS;
  mg_type = TDatabase::ParamDB->SC_MG_TYPE_SADDLE;
  if (mg_type)
    mg_level = 1;
  else
    mg_level = 0;

  mg_level = LEVELS+mg_level;

  PsBaseName = TDatabase::ParamDB->PSBASENAME;
  GrapeBaseName = TDatabase::ParamDB->GRAPEBASENAME;
  GnuBaseName = TDatabase::ParamDB->GNUBASENAME;
  ReadGrapeBaseName = TDatabase::ParamDB->READGRAPEBASENAME;
  VtkBaseName = TDatabase::ParamDB->VTKBASENAME;

  char NameString[]  = "name";
  char UString[] = "u";
  char PString[] = "p";
  char PsiString[] = "psi";
  char SurfString[] = "Surfact";
  char VorticityString[] = "vorticity";
  char ConvVortString[] = "conv_vort";
  char DivergenceString[] = "Divergence";
  char UConfString[] = "uconf";
  char OldUString[] = "oldu";
  char OldveloString[] = "Oldvelo";
  char OldPressString[] = "OldPres";
  char OldPString[] = "OldP";
  char veloString[] = "velo";
  char gridString[] = "grid";
  char U1String[] = "u1";
  char U2String[] = "u2";


// ********************************************************
   TDatabase::IteratorDB[It_EQ]->SetParam(SurfDomain);
   TDatabase::IteratorDB[It_LE]->SetParam(SurfDomain);
   TDatabase::IteratorDB[It_Finest]->SetParam(SurfDomain);
   TDatabase::IteratorDB[It_Between]->SetParam(SurfDomain);
   TDatabase::IteratorDB[It_OCAF]->SetParam(SurfDomain);

  SOld_Coll = SurfDomain->GetCollection(It_Finest, 0);
  Old_N_Cells = SOld_Coll->GetN_Cells();
  SurfDomain->GetTreeInfo(Surf_CellTree, Old_N_RootCells);
//    OutPut("Number of rootcells: "<<Old_N_RootCells<<endl);
//    if(Old_N_Cells!=Old_N_RootCells) exit(-1);
  // remove all existing vertices and joints
//   VertexDel = new TVertex*[2*Old_N_RootCells];
//  cout <<"test remesh vert"<< 3*Old_N_RootCells<<endl;
//   CurrVertex = 0;

//   for(i=0;i<Old_N_RootCells;i++)
//     {
//       cell = SOld_Coll->GetCell(i);
//       Joint = cell->GetJoint(i);
//       delete Joint;
//     }

  // remove all existing cells and joints 
  // not vertices, since all vertices are vertices of 2D domain
  for(i=0;i<Old_N_RootCells;i++)
     delete (TGridCell*)Surf_CellTree[i];

  OutPut(Old_N_RootCells<<" surface cells were deleted"<<endl);
   delete [] Surf_CellTree;

    delete SOld_Coll;

  // initialize iterators
  TDatabase::IteratorDB[It_EQ]->SetParam(Domain);
  TDatabase::IteratorDB[It_LE]->SetParam(Domain);
  TDatabase::IteratorDB[It_Finest]->SetParam(Domain);
  TDatabase::IteratorDB[It_Between]->SetParam(Domain);
  TDatabase::IteratorDB[It_OCAF]->SetParam(Domain);

// ********************************************************

     MovBoundVert[0][0]->GetCoords(x1, y1);
     MovBoundVert[0][N_MovVert[0]-1]->GetCoords(x2, y2);
     d = y2-y1;
     k = int(d/0.1); // No of intervals with step length 0.01
     if(k<2) k=2;     // minimum two intervals
     t = d/double(k);
     N_MovVert[0] = k+1;
     S_BX = new double[N_MovVert[0]];
     S_BY = new double[N_MovVert[0]];

     for(i=0;i<N_MovVert[0]-1;i++)
      {
       S_BX[i]= 0.;
 /*     if(i>0)
        {
         if(S_BX[i]<0)
          {
           T = 2*(fabs(x1) - fabs(S_BX[i]) )/d;
           S_BX[i] -= (1-T)*0.015;
          }
         else if(S_BX[i]>0)
          {
           T = 2*(fabs(x2) - fabs(S_BX[i]) )/d;
           S_BX[i] += (1-T)*0.015;
          }
       //cout<< "S_BX[i] : "<< S_BX[i]<< endl;
      }*/

     S_BY[i]= y1 + double(i)*t;
//      cout<<i<< " x :" << S_BX[i] << " -----------------y: " <<S_BY[i]<< endl;
      }
     S_BX[N_MovVert[0]-1]= 0.0;
     S_BY[N_MovVert[0]-1]= y2;

    Lx = S_BX[0];
    Ly = S_BY[0];
    Rx = S_BX[N_MovVert[0]-1];
    Ry = S_BY[N_MovVert[0]-1];
    cout<<"Lx : " << Lx << "  Ly : "<<Ly <<" Rx : " << Rx<< "  Ry : "<<Ry<< endl;

    double area = TDatabase::ParamDB->Area;

    //exit(0);
//======================================================================
// Triangular for grid generation begin
//======================================================================
  BoundPart = Domain->GetBdPart(0);
  UpdateSlipBound = (TBdLine*)BoundPart->GetBdComp(0);


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

//OutPut("MESHGEN_REF_QUALIT " << TDatabase::ParamDB->MESHGEN_REF_QUALITY << endl);

  opts<<'p'; // Constrained Delaunay Triangulation:
           // initial values - only points defined on the boundary of the domain;
           // triangulation near boundary may variate from Delaunay criterion
  opts<<"q"<<TDatabase::ParamDB->MESHGEN_REF_QUALITY;
    // Quality mesh generation with no angles smaller than given degrees;
    // an alternate minimum angle may be specified after the `q'.
  //opts<<"a0.00001";
  opts<<'e'; // Outputs a list of edges of the triangulation
  opts<<'z'; // Numbers if items starting from 0
  //opts<<"VVVV"; // Gives detailed information about what Triangle is doing
  opts<<'Q'; // Supress all explanation of what Triangle is doing, unless an error occurs
 // opts<<'Y'; // Supress adding vertices on boundary edges
//  opts<<"a0.04"; // Imposes a maximum triangle area.
  opts<<"a"<< area; // Imposes a maximum triangle area.
  opts<<ends;


  N_MovVert[0]--;     //sub right wet point as not an end point of solid surface
  N_MovVert[1]++;     //Add right wet point as begin of freesurf

  In.numberofpoints = N_MovVert[1]+N_MovVert[0];
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

  // points and segments on the solid boundary (marker=1)
  for(i=0;i<N_MovVert[0];i++) // without last point
   {
    In.pointlist[2*In_Index] = S_BX[i];
    In.pointlist[2*In_Index+1] = S_BY[i];
//     cout<<In_Index<< " x :" << S_BX[i]<< " -----------------y: " <<S_BY[i]<< endl;
    In.pointmarkerlist[In_Index] = CurrComp;
    In.segmentlist[2*In_Index] = In_Index;
    In.segmentlist[2*In_Index+1] = In_Index+1;
    In.segmentmarkerlist[In_Index] = CurrComp;
    In_Index++;
   }

   CurrComp++;
//  exit(0);

  // points and segments on the free boundary (marker=2)
  for(i=0;i<N_MovVert[1];i++) // without last point
   {
    if(i==0)
    {
     tx=Rx; ty=Ry;
     //cout<<In_Index<< " x :" << tx<< " -----------------y: " <<ty<< endl;
    }
    else
     MovBoundVert[1][i-1]->GetCoords(tx, ty);

    In.pointlist[2*In_Index] = tx;
    In.pointlist[2*In_Index+1] = ty;
   // cout<<In_Index<< " x :" << tx<< " -----------------y: " <<ty<< endl;
    In.pointmarkerlist[In_Index] = CurrComp;
    In.segmentlist[2*In_Index] = In_Index;
    In.segmentlist[2*In_Index+1] = In_Index+1;
//     if(AllowEdgeRef)
//      {
      In.segmentmarkerlist[In_Index] = CurrComp;
//      }
//     else
//      {
//       In.segmentmarkerlist[In_Index] = 100000 + CurrComp;
//      }
    In_Index++;
   }

  In.segmentlist[2*(In.numberofsegments-1)+1] = 0;

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
if(ChangeBound2)
{
  for(i=0;i<In.numberofpoints;i++)
    OutPut(i<<' '<<In.pointmarkerlist[i]<<' '<<
	   In.pointlist[2*i]<<' '<<In.pointlist[2*i+1]<<endl);
cout<<endl;
}
  */
  triangulate((char*)opts.str().c_str(), &In, &Out, (struct triangulateio *)NULL);

 /*
 if(ChangeBound)
  for(i=0;i<Out.numberofpoints;i++)
      OutPut(i<<' '<<Out.pointmarkerlist[i]<<' '<<
       Out.pointlist[2*i]<<' '<<Out.pointlist[2*i+1]<<endl);
 */

  Old_Coll = Domain->GetCollection(It_Finest, 0);
  Old_N_Cells = Old_Coll->GetN_Cells();
  Domain->GetTreeInfo(Old_CellTree, Old_N_RootCells);
//    OutPut("Number of rootcells: "<<Old_N_RootCells<<endl);
   if(Old_N_Cells!=Old_N_RootCells) exit(-1);
  // remove all existing vertices and joints
  VertexDel = new TVertex*[3*Old_N_RootCells];
//  cout <<"test remesh vert"<< 3*Old_N_RootCells<<endl;
  CurrVertex = 0;

  for(i=0;i<Old_N_Cells;i++)
    {
      cell = Old_Coll->GetCell(i);
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

   N_RootCells = Out.numberoftriangles;

  // allocate auxillary fields
  Coordinates = Out.pointlist;
  Triangles = Out.trianglelist;

  // generate all vertices
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

 // OutPut("left: "<<left<<" right: "<<right<<" top: "<<top<<" bottom: "<<bottom<<endl);
  Domain->SetBoundBox(right-left,top-bottom);
  Domain->SetBoundBoxstart(left,bottom);

 // Bound startx, starty, x length and y length

   UpdateSlipBound->SetParams(Lx, Ly, Rx-Lx, Ry-Ly);
//    cout<<"Lx : " << Lx << "  Ly : "<<Ly <<" Rx-Lx : " << Rx-Lx<< "  Ry-Rx : "<<Ry-Ly<< endl;
//    cout<<"Lx : " << Lx << "  Ly : "<<Ly <<" Rx : " << Rx<< "  Ry : "<<Ry<< endl;

 // generate cells
  CellTree = new TBaseCell*[N_RootCells];

  for (i=0;i<N_RootCells;i++)
  {
    CellTree[i] = new TMacroCell(TDatabase::RefDescDB[Triangle], LEVELS-1);

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
    PointNeighb[j + PointNeighb[j]] = i/3;
  }

  // generate edges
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

      if (Domain->GetBdPart(0)->GetBdComp(CurrComp)->GetTofXY(
            NewVertices[a]->GetX(), NewVertices[a]->GetY(), T_a) ||
          Domain->GetBdPart(0)->GetBdComp(CurrComp)->GetTofXY(
            NewVertices[b]->GetX(), NewVertices[b]->GetY(), T_b))
       {
          OutPut("Error: could not set parameter values"<<endl);
          OutPut("CurrComp : "<< CurrComp << endl);
          OutPut(NewVertices[a]<<endl);
          OutPut(NewVertices[b]<<endl);
//           exit(0);
       }

      if (Domain->GetBdPart(0)->GetBdComp(CurrComp)->GetType() != Line)
        if (ABS(T_a) < 1e-4 || ABS(T_a) > 0.9999 ||
            ABS(T_b) < 1e-4 || ABS(T_b) > 0.9999)
        {
          x = (NewVertices[a]->GetX() + NewVertices[b]->GetX()) / 2;
          y = (NewVertices[a]->GetY() + NewVertices[b]->GetY()) / 2;

          Domain->GetBdPart(0)->GetBdComp(CurrComp)->GetTofXY(x, y, T);

          if ((T_a - T)*(T - T_b) < 0)
          {
            if (ABS(T_a) < 1e-4) T_a = 1.0;
            if (ABS(T_b) < 1e-4) T_b = 1.0;
            if (ABS(T_a) > 0.9999) T_a = 0.0;
            if (ABS(T_b) > 0.9999) T_b = 0.0;
          }
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

  ORDER = 0;
//exit(0);

  VSP = TDatabase::ParamDB->VELOCITY_SPACE;

if (abs(VSP) > 20)
  {ORDER = abs(VSP) - 20;}
else if ( abs(VSP) > 10)
    {ORDER = abs(VSP) - 10;}
  else ORDER = abs(VSP);

   for(i=0;i<mg_level;i++)
    for(j=1;j<5;j++)
          delete [] Sol[j][i];

   delete [] Sol[5][0];

//    for(i=0;i<mg_level;i++)
//     for(j=7;j<9;j++)
//           delete [] Sol[j][i];


      for(i=0;i<mg_level;i++)
          delete VeloVect[1][i];

      for(i=0;i<mg_level;i++)
        for(j=0;j<2;j++)
         {
          delete Struct[j][i];
          delete SqrStruct[j][i];
         }

   for(i=0;i<mg_level;i++)
    delete SqrStruct_low[0][i];

     for(i=0;i<mg_level;i++)
      for(j=2;j<5;j++) delete FeSps_Lev[j][i];

      for(i=0;i<4;i++)
       for(j=0;j<mg_level;j++)
        delete [] Rharray[i][j];

  for(i=0;i<mg_level;i++)
   {

    switch(TDatabase::ParamDB->NSTYPE)
    {
      case 1:
        delete Mat[0][i];
        delete Mat[1][i];

        delete SqMat[12][i];
        delete SqMat[13][i];
        break;

      case 2:
        delete Mat[0][i];
        delete Mat[1][i];
        delete Mat[2][i];
        delete Mat[3][i];
        delete SqMat[12][i];
        delete SqMat[13][i];
        break;

      case 3:
        delete Mat[0][i];
        delete Mat[1][i];

        delete SqMat[0][i];
        delete SqMat[1][i];
        delete SqMat[2][i];
        delete SqMat[3][i];

        delete SqMat[4][i];
        delete SqMat[5][i];
        delete SqMat[6][i];
        delete SqMat[7][i];
        break;

      case 4:

        delete Mat[0][i];
        delete Mat[1][i];
        delete Mat[2][i];
        delete Mat[3][i];

        delete SqMat[0][i];
        delete SqMat[1][i];
        delete SqMat[2][i];
        delete SqMat[3][i];

        delete SqMat[4][i];
        delete SqMat[5][i];
        delete SqMat[6][i];
        delete SqMat[7][i];


        delete SqMat[8][i];
        delete SqMat[9][i];
        delete SqMat[10][i];
        delete SqMat[11][i];

        delete SqMat[14][i];
        delete SqMat[15][i];

        delete SqMat_low[0][i];
        delete SqMat_low[1][i];
        delete SqMat_low[2][i];
        break;
    }
  }

 int N_SCells, N_S, N_SO;
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

    Coll=Domain->GetCollection(It_Finest, 0);
    cout << endl << endl;

    Domain2DSurface(Domain, SurfDomain, N_List[0], N_List[1], i);
    Surf_Coll = SurfDomain->GetCollection(It_Finest, 0);
    N_SCells= Surf_Coll->GetN_Cells();
    FE1D_List = new FE1D[N_SCells];
    for(j=0;j<N_SCells;j++)
     FE1D_List[j] = FE1D(TDatabase::ParamDB->ANSATZ_ORDER);

    delete SurfFeSps_Lev[0][i];

    SurfFeSps_Lev[0][i] = new TFESpace1D(Surf_Coll, NameString, SurfString,
                                   FE1D_List);

    TDatabase::IteratorDB[It_EQ]->SetParam(Domain);
    TDatabase::IteratorDB[It_LE]->SetParam(Domain);
    TDatabase::IteratorDB[It_Finest]->SetParam(Domain);
    TDatabase::IteratorDB[It_Between]->SetParam(Domain);
    TDatabase::IteratorDB[It_OCAF]->SetParam(Domain);


    // get spaces for low order disc on finest geo grid
    if ((mg_type==1)&&(i<mg_level-1))
    {
      Velo_Space = new TFESpace2D(Coll,NameString, UString, BoundCondition,
                                      Non_USpace, 1, NULL);

      Pres_Space = new TFESpace2D(Coll,NameString, PString, BoundCondition,
                                      DiscP_PSpace,0, NULL);

      velocity_space_code = -1;
      pressure_space_code = 0;
      order = -1;
    }
    // get spaces of high order disc on finest geo grid
  else
    {
      GetVelocityAndPressureSpace(Coll,BoundCondition,
                                  mortarcoll, Velo_Space,
                                  Pres_Space, &pressure_space_code,
                                  TDatabase::ParamDB->VELOCITY_SPACE,
                                  TDatabase::ParamDB->PRESSURE_SPACE);

      velocity_space_code =   TDatabase::ParamDB->VELOCITY_SPACE;
      TDatabase::ParamDB->INTERNAL_PRESSURE_SPACE = pressure_space_code;
    }

    FeSps_Lev[2][i] = new TFESpace2D(Coll, NameString, PsiString, GridBoundCondition,
                                1, NULL);

    OldSurfact_space = FeSps_Lev[5][i];
    FeSps_Lev[5][i] = new TFESpace2D(Coll, NameString, SurfString, GridBoundCondition,
                                 TDatabase::ParamDB->VELOCITY_SPACE,   NULL);

    N_S = SurfFeSps_Lev[0][i]->GetN_DegreesOfFreedom();
    N_SO =  FeSps_Lev[5][i]->GetN_DegreesOfFreedom();


    N_U = Velo_Space->GetN_DegreesOfFreedom();
    N_P = Pres_Space->GetN_DegreesOfFreedom();
    N_G = FeSps_Lev[2][i]->GetN_DegreesOfFreedom();
    N_BoundaryNodes = N_G - FeSps_Lev[2][i]->GetN_Inner();

     N_array[0][i] = N_U;
     N_array[1][i] = N_P;
     N_array[2][i] = N_G;


    N_Active = Velo_Space->GetActiveBound();
    cout << "N_ActiveBound: " << N_Active << endl;

 // if (i<LEVELS)
    {

      FeSps_Lev[3][i] =  new TFESpace2D(Coll,NameString, PsiString, BoundCondition,
                                          1, NULL);

      N_V = FeSps_Lev[3][i]->GetN_DegreesOfFreedom();

      FeSps_Lev[4][i] =  new TFESpace2D(Coll,NameString, UString, BoundCondition,
                                      ContP_USpace,1, mortarcoll);

      N_Vort = FeSps_Lev[4][i]->GetN_DegreesOfFreedom();
    }


    cout << "No. of Grid BoundaryNodes: " << N_BoundaryNodes << endl;
    cout << "No. of DOF of Grid cells: " << N_G << endl;
     // find on finest grid all vertices on moving boundary

   if(i==mg_level-1)
    {

      N_Cells = Coll->GetN_Cells();
      N_SlipJoints = 0;
      N_FreeJoints = 0;
      N_MovVert[2] = 0;
      N_Velo_IsoPoints = 0;
      for(j=0;j<N_Cells;j++)
      {
        Me = Coll->GetCell(j);
        k = Me->GetN_Edges();
        for(l=0;l<k;l++)
        {
          Joint = Me->GetJoint(l);
          if(Joint->GetType() == BoundaryEdge)
            {
               N_SlipJoints++,
               N_Velo_IsoPoints += ORDER-1;
            }
          if(Joint->GetType() == IsoBoundEdge)
          {
              N_FreeJoints++;
              N_MovVert[2] += ORDER-1;
              N_Velo_IsoPoints += ORDER-1;
          }
        } // endfor l
      } // endfor j

      cout << "N_SlipJoints: " << N_SlipJoints << endl;
      cout << "N_FreeJoints: " << N_FreeJoints << endl;
      cout << "N_MovVert_QuardPts: " << N_MovVert[2] << endl;
      OutPut("Number of Elements = " <<N_Cells<< endl);

      delete [] MovBoundVert[0];
      delete [] MovBoundVert[1];
      delete [] MovBoundVert[2];
      delete [] MovBoundVert[3];
      delete [] Slip_Joint[0];
      delete [] N_array[3];
      delete [] N_array[4];

      N_MovVert[0] = N_SlipJoints+1; //addding the right wetting point to BC2(slip)
      N_MovVert[1] = N_FreeJoints;
      MovBoundVert[1] = new TVertex*[N_MovVert[1]];
      MovBoundVert[0] = new TVertex*[N_MovVert[0]+N_MovVert[1]];
      MovBoundVert[2] = new TVertex*[N_MovVert[2]];
      MovBoundVert[3] = new TVertex*[ORDER];

      Slip_Joint[0] = new TBoundEdge*[N_SlipJoints+N_FreeJoints];
      N_array[3] = new int [N_MovVert[1]];
      N_array[4] = new int [N_MovVert[1]];

      m = 0;
      m1 = 0;
      m2 = 0;
      m3 = 0;
      m4 = 0;

     for(j=0;j<N_Cells;j++)
      {
        Me = Coll->GetCell(j);
        k = Me->GetN_Edges();
        for(l=0;l<k;l++)
        {
          Joint = Me->GetJoint(l);
          if(Joint->GetType() == BoundaryEdge)
           {
              MovBoundVert[0][m] = Me->GetVertex(l);
              Slip_Joint[0][m4++] = (TBoundEdge *)Joint;
              m++;
           }
          else if(Joint->GetType() == IsoBoundEdge)
           {
            Me->GetVertex(l)->GetCoords(TX[0], TY[0]);
            Me->GetVertex((l+1) % k)->GetCoords(TX[1], TY[1]);

            MovBoundVert[1][m1] = Me->GetVertex(l);
            Free_Cells[m1] = Me;
            N_array[3][m1] = l;
            N_array[4][m1] = j;
            m1++;

//             Free_Joint[0][m2] = (TIsoBoundEdge *)Joint;
            Free_J = (TIsoBoundEdge *)Joint;

//             Free_Joint[0][m2]->GeneratemidVert(ORDER-1, TX, TY);
            Free_J->GeneratemidVert(ORDER-1, TX, TY);
//             MovBoundVert[3] = Free_Joint[0][m2]->GetVertices();
            MovBoundVert[3] = Free_J->GetVertices();
            for(i3=0;i3<(ORDER-1);i3++)
             {
//                MovBoundVert[3][i3]->GetCoords(BoundVert[0][m3], BoundVert[1][m3]);
               MovBoundVert[2][m3] =  MovBoundVert[3][i3];
               m3++;
             }
            m2++;
           }
         } // endfor l
      } // endfor j
 
 // sort free bound vertices from +x to -x
 SortIsoVertices_axial3D(Free_Cells, MovBoundVert[1], N_array[3], N_array[4], N_MovVert[1]);
  // sort free bound isovertices after SortIsoVertices() only
   // according to the Free_Cells, Edge_No ordered in SortIsoVertices(()
//  for(j=0;j<N_MovVert[1];j++)
//       {
//         Me = Free_Cells[j];
//         Free_Joint[0][j] = Me->GetJoint(j);
//       }
 

//  SortQuardIsopts(Free_Cells, MovBoundVert[2], N_array[5], N_MovVert[2]);

  // sort bound points on solid surface from -x to +x
  // first two points are wettting points
    for(k=0;k<N_MovVert[0]-2;k++)
      {
      for(l=k+1;l<N_MovVert[0]-1;l++)
      {
        MovBoundVert[0][k]->GetCoords(x, y);
	MovBoundVert[0][l]->GetCoords(tx, ty);
	if(ty > y)
	 {
	  temp_Mov = MovBoundVert[0][k];
          MovBoundVert[0][k] = MovBoundVert[0][l];
          MovBoundVert[0][l] = temp_Mov;

	  tempSlip_Joint = Slip_Joint[0][k];
	  Slip_Joint[0][k] = Slip_Joint[0][l];
	  Slip_Joint[0][l] = tempSlip_Joint;
	 }
       }
      }

    //Adding first point of freesurface i.e wetting point on MovBoundVert[0]
    MovBoundVert[0][N_MovVert[0]-1] = MovBoundVert[1][0];

   for(k=1;k<N_MovVert[1];k++)
      MovBoundVert[1][k-1] = MovBoundVert[1][k];
    N_MovVert[1]--;

       // exit(0);

  } // endif i==mg_level-1

    // build matrices
    SqrStruct[1][i] = new TSquareStructure2D(FeSps_Lev[2][i]);
    SqrStruct[1][i]->Sort();
    Struct[0][i] = new TStructure2D(Pres_Space, Velo_Space);
    Struct[1][i] = new TStructure2D(Velo_Space, Pres_Space);
    SqrStruct[0][i] = new TSquareStructure2D(Velo_Space);
    SqrStruct[0][i]->Sort();

    SqrStruct_low[0][i] = new TSquareStructure1D(SurfFeSps_Lev[0][i]);
    SqrStruct_low[0][i]->Sort();

    switch(TDatabase::ParamDB->NSTYPE)
    {
      case 1:
        Mat[0][i] = new TMatrix2D(Struct[0][i]);
        Mat[1][i] = new TMatrix2D(Struct[0][i]);

        SqMat[12][i] = new TSquareMatrix2D(SqrStruct[0][i]);
        SqMat[13][i] = new TSquareMatrix2D(SqrStruct[0][i]);
        break;

      case 2:
        Mat[0][i] = new TMatrix2D(Struct[0][i]);
        Mat[1][i] = new TMatrix2D(Struct[0][i]);
        Mat[2][i] = new TMatrix2D(Struct[1][i]);
        Mat[3][i] = new TMatrix2D(Struct[1][i]);

        SqMat[12][i] = new TSquareMatrix2D(SqrStruct[0][i]);

        SqMat[13][i] = new TSquareMatrix2D(SqrStruct[0][i]);
        break;

      case 3:
        Mat[0][i] = new TMatrix2D(Struct[0][i]);
        Mat[1][i] = new TMatrix2D(Struct[0][i]);

        SqMat[0][i] = new TSquareMatrix2D(SqrStruct[0][i]);
        SqMat[1][i] = new TSquareMatrix2D(SqrStruct[0][i]);
        SqMat[2][i] = new TSquareMatrix2D(SqrStruct[0][i]);
        SqMat[3][i] = new TSquareMatrix2D(SqrStruct[0][i]);

        SqMat[4][i] = new TSquareMatrix2D(SqrStruct[0][i]);
        SqMat[5][i] = new TSquareMatrix2D(SqrStruct[0][i]);
        SqMat[6][i] = new TSquareMatrix2D(SqrStruct[0][i]);
        SqMat[7][i] = new TSquareMatrix2D(SqrStruct[0][i]);
        break;

      case 4:

        Mat[0][i] = new TMatrix2D(Struct[0][i]);
        Mat[1][i] = new TMatrix2D(Struct[0][i]);
        Mat[2][i] = new TMatrix2D(Struct[1][i]);
        Mat[3][i] = new TMatrix2D(Struct[1][i]);

        SqMat[0][i] = new TSquareMatrix2D(SqrStruct[0][i]);
        SqMat[1][i] = new TSquareMatrix2D(SqrStruct[0][i]);
        SqMat[2][i] = new TSquareMatrix2D(SqrStruct[0][i]);
        SqMat[3][i] = new TSquareMatrix2D(SqrStruct[0][i]);

        SqMat[4][i] = new TSquareMatrix2D(SqrStruct[0][i]);
        SqMat[5][i] = new TSquareMatrix2D(SqrStruct[0][i]);
        SqMat[6][i] = new TSquareMatrix2D(SqrStruct[0][i]);
        SqMat[7][i] = new TSquareMatrix2D(SqrStruct[0][i]);

        SqMat[14][i] = new TSquareMatrix2D(SqrStruct[0][i]);
        SqMat[15][i] = new TSquareMatrix2D(SqrStruct[0][i]);

        SqMat[8][i] = new TSquareMatrix2D(SqrStruct[1][i]);
        SqMat[9][i] = new TSquareMatrix2D(SqrStruct[1][i]);
        SqMat[10][i] = new TSquareMatrix2D(SqrStruct[1][i]);
        SqMat[11][i] = new TSquareMatrix2D(SqrStruct[1][i]);

        SqMat_low[0][i] = new TSquareMatrix1D(SqrStruct_low[0][i]);
        SqMat_low[1][i] = new TSquareMatrix1D(SqrStruct_low[0][i]);
        SqMat_low[2][i] = new TSquareMatrix1D(SqrStruct_low[0][i]);

        break;
    }

    N_Unknowns = 2*N_U + N_P;

    OutPut("dof velocity : "<< setw(10) << 2* N_U << endl);
    OutPut("dof pressure : "<< setw(10) << N_P << endl);
    OutPut("dof all      : "<<  setw(10) << N_Unknowns  << endl);
    OutPut("dof freesurf surfactant : "<< setw(5) << N_S << endl);

    Rharray[0][i] = new double[N_Unknowns];
    memset(Rharray[0][i], 0, N_Unknowns*SizeOfDouble);

    New_Sol = new double[N_Unknowns];
    Sol[1][i] = new double[N_Unknowns];
    Sol[2][i] = new double[N_Unknowns];

    memset(New_Sol, 0, N_Unknowns*SizeOfDouble);
    memset(Sol[1][i], 0, N_Unknowns*SizeOfDouble);
    memset(Sol[2][i], 0, N_Unknowns*SizeOfDouble);
    if (i==mg_level-1)
    Sol[5][0] =  new double[N_Unknowns];

    Rharray[1][i] = new double[N_Unknowns];
    memset(Rharray[1][i], 0, N_Unknowns*SizeOfDouble);

    Rharray[2][i] = new double[2*N_G];
    Sol[4][i] = new double[2*N_G];
    memset(Sol[4][i], 0, 2*N_G*SizeOfDouble);
    memset(Rharray[2][i], 0, 2*N_G*SizeOfDouble);

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
          MGLevel = new TNSE_MGLevel1(i, SqMat[13][i], Mat[0][i], Mat[1][i],
                                      Struct[1][i], Rharray[1][i], New_Sol, n_aux,  alpha,
                                      velocity_space_code ,
                                      pressure_space_code,NULL);
        break;

        case 2:
          MGLevel = new TNSE_MGLevel2(i, SqMat[13][i], Mat[0][i], Mat[1][i],
                                      Mat[2][i], Mat[3][i],
                                      Rharray[1][i], New_Sol, n_aux, alpha,
                                      velocity_space_code ,
                                      pressure_space_code,NULL);
        break;

        case 3:
        MGLevel = new TNSE_MGLevel3(i, SqMat[4][i], SqMat[5][i],
                                      SqMat[6][i], SqMat[7][i],
                                      Mat[0][i], Mat[1][i],
                                      Struct[1][i],
                                      Rharray[1][i], New_Sol, n_aux, alpha,
                                      velocity_space_code ,
                                      pressure_space_code,NULL);
        break;

        case 4:
          MGLevel = new TNSE_MGLevel4(i, SqMat[4][i], SqMat[5][i],
                                       SqMat[6][i], SqMat[7][i],
                                       Mat[0][i], Mat[1][i],
                                       Mat[2][i], Mat[3][i],
                                       Rharray[1][i], New_Sol, n_aux, alpha,
                                       velocity_space_code ,
                                       pressure_space_code,NULL);
       break;
      } // end switch(NSTYPE)
//     TMG[0]->Replace(i, TMGLevel[0][i]);
    TMG[0]->Replace(i, MGLevel);
    }

    New_u  = new TFEVectFunct2D(Velo_Space, UString,  UString,  New_Sol, N_U, 2);
    New_p  = new TFEFunction2D(Pres_Space, PString,  PString,  New_Sol+2*N_U, N_P);

   if(i==mg_level-1)
    {
     New_u->GetComponent(0)->Intpol_Remesh(UPArrays[0][i], U1String);
//   U1String bcz in Intpol_Remesh function U2String and Y==0 means UPArrays[1][i]=0
     New_u->GetComponent(1)->Intpol_Remesh(UPArrays[1][i], U1String);
     New_p->Intpol_Remesh(UPArrays[2][i], PString);
    }


     ValuesVX = New_u->GetValues();
//      ValuesVY = ValuesVX + N_U;

     GlobalNumbers = Velo_Space->GetGlobalNumbers();
     BeginIndex = Velo_Space->GetBeginIndex();

     for(j=0;j<N_Cells;j++)
      {
        Me = Coll->GetCell(j);
        k = Me->GetN_Edges();
        for(l=0;l<k;l++)
        {
          Joint = Me->GetJoint(l);
          if(Joint->GetType() == BoundaryEdge)
          {
           //cout<< " cell : " << j <<"  BoundaryEdge Joint : "<< l <<endl;
           FeId = Velo_Space->GetFE2D(j, Me);
           FeDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(FeId);
           JointDOF = FeDesc->GetJointDOF(l);
           N_DOF = FeDesc->GetN_JointDOF();
           for(m=0;m<N_DOF;m++)
            {
             m1 = GlobalNumbers[BeginIndex[j]+JointDOF[m]];
             Velo_Space->GetDOFPosition(m1, x, y);
//              if(y == 0.) ValuesVY[m1] = 0.;
             if(fabs(x)<1e-12) ValuesVX[m1] = 0.;
            }
         }
        }
      }


//     Space_Old_Coll = FeSps_Lev[0][i]->GetCollection();
//     delete Space_Old_Coll;
    delete FeSps_Lev[0][i];
    delete FeSps_Lev[1][i];
    
    FeSps_Lev[0][i]=Velo_Space;
    FeSps_Lev[1][i]=Pres_Space;

    delete [] Sol[0][i];
    Sol[0][i] = New_Sol;

    delete VeloVect[0][i];

//     delete UPArrays[0][i];
//     delete UPArrays[1][i];

   VeloVect[0][i] = New_u;
   UPArrays[0][i] = VeloVect[0][i]->GetComponent(0);
   UPArrays[1][i] = VeloVect[0][i]->GetComponent(1);

   delete UPArrays[2][i];

   UPArrays[2][i] = New_p;

   Sol[3][i] = new double[2*N_G];
   memset(Sol[3][i], 0, 2*N_G*SizeOfDouble);
   VeloVect[1][i] = new TFEVectFunct2D(FeSps_Lev[2][i], veloString, gridString,
                                      Sol[3][i], N_G, 2);

//    delete UPArrays[3][i];
//    delete UPArrays[4][i];
    UPArrays[3][i] = VeloVect[1][i]->GetComponent(0);
    UPArrays[4][i] = VeloVect[1][i]->GetComponent(1);

// surfactant
   Rharray[3][i] = new double[N_S];
   memset(Rharray[3][i], 0, N_S*SizeOfDouble);

   delete  [] Sol[6][i];
   Sol[6][i] = new double[N_S];
   memset(Sol[6][i], 0, N_S*SizeOfDouble);
   delete SArrays[0][i];
   SArrays[0][i] = new TFEFunction1D(SurfFeSps_Lev[0][i], SurfString, SurfString,
                                       Sol[6][i], N_S);


// //  surfactant in entire 2D domain
   oldssol = Sol[7][i];
   Sol[7][i] = new double[N_SO];
   memset(Sol[7][i], 0, N_SO*SizeOfDouble); 
   New_s = new TFEFunction2D(FeSps_Lev[5][i], SurfString, SurfString, 
                                       Sol[7][i], N_SO);

   New_s->Intpol_Remesh(UPArrays[5][i], PString);
   MapDomainToSurf(New_s, SArrays[0][i], N_List[0][i], N_List[1][i]);
   memset(Sol[7][i], 0, N_SO*SizeOfDouble);

   delete UPArrays[5][i];
   UPArrays[5][i] = New_s;
   delete OldSurfact_space;

   delete [] Sol[8][i];
   Sol[8][i] = new double[2*N_SO];
   memset(Sol[8][i], 0, 2*N_SO*SizeOfDouble);
   delete UPArrays[6][i];
   delete UPArrays[7][i];
   UPArrays[6][i] = new TFEFunction2D(FeSps_Lev[5][i], PsiString, PsiString, 
                                       Sol[8][i], N_SO);
   UPArrays[7][i] = new TFEFunction2D(FeSps_Lev[5][i], PsiString, PsiString,
                                       Sol[8][i]+N_SO, N_SO);

 } // end for mg_level


   for(i=0;i<CurrVertex;i++)
    delete VertexDel[i];

  delete [] VertexDel;
  OutPut(CurrVertex<<" vertices were deleted"<<endl);

  // remove all existing cells and joints
  for(i=0;i<Old_N_RootCells;i++)
     delete (TGridCell*)Old_CellTree[i];
  OutPut(Old_N_RootCells<<" cells were deleted"<<endl);
   delete [] Old_CellTree;

    delete Old_Coll;

   memset(Sol[5][0], 0, N_Unknowns*SizeOfDouble);

//exit(0);

//   delete [] Bd_CngVert[0];
//   delete [] Bd_CngVert[1];
//   delete [] Bd_CngVert;
  delete [] S_BX;
  delete [] S_BY;
//   delete [] taux;
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
  int *BeginIndex, *GlobalNumbers;
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

    switch(RefElement)
      {
        case BFUnitSquare:
          RefTrans = QuadIsoparametric;
          F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
          ((TQuadIsoparametric *)F_K)->SetCell(Me);
        break;

        case BFUnitTriangle:
          RefTrans = TriaIsoparametric;
          F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
          ((TTriaIsoparametric *)F_K)->SetCell(Me);
        break;
      } // endswitch


      l = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
      LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(2*l);
      qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
      qf1->GetFormulaData(N_LinePoints, LineWeights, zeta);
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
              val  =c0*(d1*e1 + d2*e2);
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
  int i, j, k, l, DOF_R, DOF_L, TDOF_R, m, TN_DOF_Local;
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
  int N_LinePoints, local_dof;
  double *LineWeights, *zeta;
  double x0, y0, x1, y1,tx,ty,mod_t, x, y;
  int N_BaseFunct, *N_BaseFuncts, TN_BaseFunct;
  double **uref, **uxiref, **uetaref;
  double uorig[MaxN_BaseFunctions2D], uxorig[MaxN_BaseFunctions2D];
  double uyorig[MaxN_BaseFunctions2D];
  BaseFunct2D *BaseFuncts;
  double r2, r;
  int *KCol, *RowPtr, *JointDOF, N_DOF, *TJointDOF;
  double *ValuesA11, *ValuesA12, *ValuesA21, *ValuesA22;
  int *BeginIndex, *GlobalNumbers, *DOF, *TDOF, TestDOF, AnsatzDOF;
  int index1, index2;
  double val, theta, factor1, factor2, angle, T_val[3], T_Marangoni, *Tvalues;
  int count=0, count1=0, count2=0;
  double  X_B[100], Y_B[100], r_axial, T_DOF[10], d1, d2, e1, e2, ngrad_test, ngrad_ansatz;
  int *TGlobalNumbers, *TBeginIndex;
  double ngrad_Gamma, GammaE1, GammaE2;

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
  double We = TDatabase::ParamDB->WB_NR, Gamma;
  double E = TDatabase::ParamDB->P13;  // surfactant elasticity E
  int EOS = int(TDatabase::ParamDB->P14); //Equation of state, 0 linear, 1 non-linear
  double D = TDatabase::ParamDB->P15 / TDatabase::ParamDB->P10 ; //\Gamma_0/Gamma_\infty


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
     TN_BaseFunct = N_BaseFuncts[TFEId];
     TJointDOF = TFeDesc->GetJointDOF(IJoint);
     TN_DOF_Local = TFeDesc->GetN_JointDOF();
     TDOF = TGlobalNumbers + TBeginIndex[i];

     if(TN_BaseFunct != N_BaseFunct)
      {
        cout << "TN_BaseFunct != N_BaseFunct in free surface.integral" << endl;
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

       T_val[0] = 0.; T_val[1] = 0.; T_val[2] = 0.;
       for(l=0;l<TN_DOF_Local;l++)
        {
          local_dof   = TJointDOF[l];      // assumed that the velo space and 2D surfactant space are same fe space
          m = TDOF[local_dof];
          T_val[0] += Tvalues[m]*uorig[local_dof];  // Surfactant C
//           T_val[0] +=  uorig[local_dof];  // Surfactant C
          T_val[1] += Tvalues[m]*uxorig[local_dof];  // C_x
          T_val[2] += Tvalues[m]*uyorig[local_dof];  // C_y
        } // for(l=0;l<TN_

// cout << " T_val[0] " << T_val[0] <<endl;


        if(fabs(T_val[0])>1. || T_val[0]<0. )
            OutPut(i<< "x : "<<X_B[k]<< " y: " << Y_B[k] <<"  Surfactant exceeds the reference value, T_val= " <<T_val[0]<<endl);

// exit(0);
          // Multiply with time step dt in the main program not here
          r = normn/We;

          if(EOS==0)
            Gamma =(1. + E*(D - T_val[0]) );
          else
             Gamma =(1. + E*log(1. - T_val[0]) );


//        (c_1\sigma_sa) \tau\cdot \grad T,       norm for integral weights
          if(EOS==0)
            T_Marangoni = normn*E*(t0*T_val[1] + t1*T_val[2])/We;
          else
            T_Marangoni = normn*E*( t0*T_val[1] + t1*T_val[2]) / ( We*(1. - T_val[0]) );
// cout<< " r " <<r << " T_Marangoni " << T_Marangoni<<endl;

          for(l=0;l<N_BaseFunct;l++)
          {
           TestDOF = DOF[l];

           // updating rhs
            ngrad_test= n0*uxorig[l] + n1*uyorig[l];
            d1 = uxorig[l] - ngrad_test*n0;
            d2 = uyorig[l] - ngrad_test*n1;

            d1 *= Gamma;
            d2 *= Gamma;

            ngrad_Gamma = n0*T_val[1] + n1*T_val[2];
            GammaE1 = (T_val[1] - ngrad_Gamma*n0)*uorig[l];
            GammaE2 = (T_val[2] - ngrad_Gamma*n1)*uorig[l];

          if(EOS==0)
            {
             GammaE1 *=E;
             GammaE2 *=E;
            }
          else
            {
             GammaE1 *=(E/(1.- T_val[0]));
             GammaE2 *=(E/(1.- T_val[0]));
            }
//             GammaE1 = 0.;
//             GammaE2 = 0.;
//  rhs1
            val = r_axial*( (1-n0*n0)*(d1 - GammaE1) - n0*n1*(d2 -GammaE2));
            val +=(Gamma*uorig[l]);
            val *= LineWeights[k]*r;
            rhs1[TestDOF] -= val;
//     Marangoni convection
            val = r_axial*t0*uorig[l];
            val *= LineWeights[k]*T_Marangoni;
            rhs1[TestDOF] -= val;

// rhs2
            val =  r_axial*( -n1*n0*(d1 - GammaE1) + (1-n1*n1)*(d2 -GammaE2) );
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

              val =(d1 - GammaE1)*e1 + (d2 -GammaE2)*e2
                   + Gamma*(uorig[l]*uorig[m]/(r_axial*r_axial));
              val *= dt*LineWeights[k]*r*r_axial;
              // cout << "A11: " << TestDOF << " ";
              // cout << AnsatzDOF << " " << val << endl;
              ValuesA11[index1] += val;

              val = (d1 - GammaE1)*e1 + (d2 -GammaE2)*e2;
              val *= dt*LineWeights[k]*r*r_axial;

              // cout << "A22: " << TestDOF << " ";
              // cout << AnsatzDOF << " " << val << endl;
              ValuesA22[index1] += val;

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




void ReParam_axial3D_U(int N_E, TBaseCell **cell, int *EdgeNo,  int *CellNo, TFEVectFunct2D *Velocity,    
                        TFEFunction2D *Surfactant)
{
  int i, j, ISpline, N_Splines, N_V, ORDER, VSP, m, m1, k, i3, USpline, FeDof;
  double *h, *t, u0, u1, u2;
  double *a, *b, *c, *x, *y, teta;
  double *rhs, *u1rhs, *u2rhs, *srhs, *Mx, *My,*Mu1, *Mu2, *Msurf, *Params, *Param9, *FEParams;
  double phi1, phi2, phi3, phi4, X, Y, T;
  double dx0, dy0, dx1, dy1, surf, *u1_spl, *u2_spl, *surf_spl;
  TIsoBoundEdge *isojoint;
  TVertex **IsoVertices;
  TJoint *Joint;
  TBaseCell *Me;
  TFESpace2D *VelocitySpace, *SurfactSpace;
  int *VeloBeginIndex, *VeloGlobalNumbers, *JointDOF, *DOF, N_DOF_Joint, *U_DOF;
  int *SurfactBeginIndex, *SurfactGlobalNumbers, *SJointDOF, *SDOF, SN_DOF_Joint,*Surf_DOF;
  double *ValuesUX, *ValuesUY, *Surfact;
  FE2D FEId, SFEId;
  TFE2D *ele, *Sele;
  TFEDesc2D *FeDesc, *SFeDesc;
  TCollection *coll;

  ORDER = 0;
  VSP = TDatabase::ParamDB->VELOCITY_SPACE;

  if (abs(VSP) > 20)
   {ORDER = abs(VSP) - 20;}
  else if ( abs(VSP) > 10)
   {ORDER = abs(VSP) - 10;}
  else ORDER = abs(VSP);

  N_V = N_E+1 + N_E*(ORDER-1);

  N_Splines = N_V-1;
  h = new double[N_Splines+1];
  t = new double[N_Splines+1];
  a = new double[N_Splines+1];
  b = new double[N_Splines+1];
  c = new double[N_Splines+1];
  rhs = new double[N_Splines+1];
  u1rhs = new double[N_Splines+1];
  u2rhs = new double[N_Splines+1];
  srhs = new double[N_Splines+1];
  u1_spl = new double[N_Splines+1];
  u2_spl = new double[N_Splines+1];
  surf_spl = new double[N_Splines+1];
  Mu1 = new double[N_Splines+1];
  Mu2 = new double[N_Splines+1];
  Msurf = new double[N_Splines+1];
  Mx = new double[N_Splines+1];
  My = new double[N_Splines+1];
  Params = new double [10*N_Splines];
  Param9 = new double [N_Splines+1];
  FEParams = new double [3*2*N_Splines]; // 3 fe functions, u1, u2, surfact

  x = new double[N_V];
  y = new double[N_V];
//   UX = new double[N_V];
//   UY = new double[N_V];
  U_DOF = new int[N_V];
  Surf_DOF = new int[N_V];

  VelocitySpace = Velocity->GetFESpace2D();
  VeloBeginIndex = VelocitySpace->GetBeginIndex();
  VeloGlobalNumbers = VelocitySpace->GetGlobalNumbers();
  ValuesUX = Velocity->GetValues();
  ValuesUY = ValuesUX + Velocity->GetLength();


  SurfactSpace = Surfactant->GetFESpace2D();
  SurfactBeginIndex = SurfactSpace->GetBeginIndex();
  SurfactGlobalNumbers = SurfactSpace->GetGlobalNumbers();
  Surfact = Surfactant->GetValues();

  coll = VelocitySpace->GetCollection();

   m = 0;
   m1 = 0;
   for(i=0;i<N_E;i++) // i<N_E
   {
//     Me = coll->GetCell(CellNo[i]);
    Me = cell[i];
    Me->GetVertex(EdgeNo[i])->GetCoords(x[m], y[m]);
    m++;

    Joint = cell[i]->GetJoint(EdgeNo[i]);
    isojoint = (TIsoBoundEdge *)Joint;
    k = isojoint->GetN_Vertices();
    if(k==ORDER-1)
     {
      IsoVertices = isojoint->GetVertices();
      for(i3=0;i3<k;i3++)
       {   
        IsoVertices[i3]->GetCoords(x[m], y[m]);
//          cout<< i<<" FreeGaus " << (180/Pi)*atan2(y[m], x[m]) <<endl;
        m++;
       } 
     }
    else
     {
      // only second order conforming elements implimented
      cout<< " No match in isopoints per free edge "<<endl;
      exit(0);
     }

// for velocity
    FEId = VelocitySpace->GetFE2D(CellNo[i], Me);
    ele = TFEDatabase2D::GetFE2D(FEId);
    FeDesc = ele->GetFEDesc2D();   // fe descriptor
    JointDOF = FeDesc->GetJointDOF(EdgeNo[i]);
    N_DOF_Joint = FeDesc->GetN_JointDOF();
    DOF = VeloGlobalNumbers + VeloBeginIndex[CellNo[i]];

// for surfactant
    SFEId = SurfactSpace->GetFE2D(CellNo[i], Me);
    Sele = TFEDatabase2D::GetFE2D(SFEId);
    SFeDesc = Sele->GetFEDesc2D();   // fe descriptor
    SJointDOF = SFeDesc->GetJointDOF(EdgeNo[i]);
    SN_DOF_Joint = SFeDesc->GetN_JointDOF();
    SDOF = SurfactGlobalNumbers + SurfactBeginIndex[CellNo[i]];


    if((N_DOF_Joint-1)!=ORDER)
     {
      // only second order conforming elements implimented
      cout<< " No match in no. velo-nodal functunals and vertices on the free surface edge  "<<endl;
      cout<< " (N_DOF_Joint-1) " << N_DOF_Joint-1 << " ORDER " << ORDER <<endl;
      exit(0);
     }

    if(i !=N_E-1)// -1 due to end dof will be the start dof of the next edge except on last edge
     N_DOF_Joint--; // assumed that velocity and surfactant having same no. of dof on edge

//   cout << " CellNo[i] " << CellNo[i] << endl;
     for (i3=0;i3<N_DOF_Joint;i3++)
       {
         U_DOF[m1] = DOF[JointDOF[i3]]; // needed for later update
         u1_spl[m1] = ValuesUX[DOF[JointDOF[i3]]];
         u2_spl[m1] = ValuesUY[DOF[JointDOF[i3]]];
         Surf_DOF[m1] = SDOF[SJointDOF[i3]]; // needed for later update
         surf_spl[m1] = Surfact[SDOF[SJointDOF[i3]]];
//  cout << "  SJointDOf " << JointDOF[i3] << " DOF " << DOF[JointDOF[i3]] <<  endl;
         m1++;
       }

   } // for(i=0;i<N_E


// exit(0);

//   end vertex of the freeboundary
  k = cell[N_E-1]->GetN_Edges();
  cell[N_E-1]->GetVertex((EdgeNo[N_E-1]+1) % k)->GetCoords(x[m], y[m]);

  if(m+1!=m1)
   {
      // only second order conforming elements implimented
      cout<< " No match in no. velo-nodal functunals and vertices on the free surface edge  "<<endl;
      cout<< " m " << m << " m1 " << m1 <<endl;
      exit(0);
     }


//  for(i=0;i<N_V;i++)
//    OutPut("OldX: "<< i <<' '<<x[i] <<' '<< y[i] <<endl);
//  OutPut("OldX: "<< i <<' '<<x[0] <<' '<< y[0] <<endl);
//  OutPut("OldX: "<< i <<' '<<x[N_V-2] <<' '<< y[N_V-2] <<endl);
//  OutPut("OldX: "<< i <<' '<<x[N_V-1] <<' '<< y[N_V-1] <<endl);
// cout << "Surfact[Surf_DOF[0]] " <<Surfact[Surf_DOF[0]] << " Surfact[Surf_DOF[m1]] " << Surfact[Surf_DOF[m1-1]] <<endl;

  h[0] = 0.0; t[0] = 0.0;

 for(i=1;i<=N_Splines;i++)
  {
    h[i] = sqrt((x[i]-x[i-1])*(x[i]-x[i-1])+(y[i]-y[i-1])*(y[i]-y[i-1]));
    t[i] = t[i-1] + h[i];
  }

  dx0 = (x[1]-x[0])/h[1];
  dy0 = (y[1]-y[0])/h[1];

  dx1 = (x[N_Splines]-x[N_Splines-1])/h[N_Splines];
  dy1 = (y[N_Splines]-y[N_Splines-1])/h[N_Splines];


  a[0] = 2.; c[0] = 1.; rhs[0] = -6./h[1]*(dx0 - (x[1]-x[0])/h[1]);
  for(i=1;i<N_Splines;i++)
  {
    a[i] = 2.;  
    b[i] = h[i]/(h[i]+h[i+1]); // \mu_i in PhD thesis
    c[i] = h[i+1]/(h[i]+h[i+1]); // \lambda_i in PhD thesis
    rhs[i] = 6./(h[i]+h[i+1])*((x[i+1]-x[i])/h[i+1]-(x[i]-x[i-1])/h[i]);
  }
  b[N_Splines] = 1.; a[N_Splines] = 2.;
  rhs[N_Splines] = 6./h[N_Splines]*(dx1 - (x[N_Splines]-x[N_Splines-1])/h[N_Splines]);

  Solver_3dia(N_Splines, a, b, c, rhs, Mx);

  rhs[0] = -6./h[1]*(dy0 - (y[1]-y[0])/h[1]);
  for(i=1;i<N_Splines;i++)
  {
    rhs[i] = 6./(h[i]+h[i+1])*((y[i+1]-y[i])/h[i+1]-(y[i]-y[i-1])/h[i]);
  }
  rhs[N_Splines] = 6./h[N_Splines]*(dy1 - (y[N_Splines]-y[N_Splines-1])/h[N_Splines]);

  Solver_3dia(N_Splines, a, b, c, rhs, My);

  for(i=0;i<N_Splines;i++)
  {
    ISpline = i*10;
    Params[ISpline    ] = x[i]; 
    Params[ISpline + 1] = y[i];
    Params[ISpline + 2] = x[i+1]; 
    Params[ISpline + 3] = y[i+1];
    Params[ISpline + 4] = -Mx[i]*h[i+1]*h[i+1]/2. +
                          ((x[i+1]-x[i])/h[i+1]-h[i+1]/6.*(Mx[i+1]-Mx[i]))*h[i+1];

//     Params[ISpline + 4] = Mx[i]*h[i];
    Params[ISpline + 5] = -My[i]*h[i+1]*h[i+1]/2. +
                          ((y[i+1]-y[i])/h[i+1]-h[i+1]/6.*(My[i+1]-My[i]))*h[i+1];
    Params[ISpline + 6] = Mx[i+1]*h[i+1]*h[i+1]/2. +
                          ((x[i+1]-x[i])/h[i+1]-h[i+1]/6.*(Mx[i+1]-Mx[i]))*h[i+1];

//     Params[ISpline + 6] = -Mx[i+1];
    Params[ISpline + 7] = My[i+1]*h[i+1]*h[i+1]/2. +
                          ((y[i+1]-y[i])/h[i+1]-h[i+1]/6.*(My[i+1]-My[i]))*h[i+1];
    Params[ISpline + 8] = t[i+1]/t[N_Splines];
    Params[ISpline + 9] = 0.;

   //cout<<"  "<<Params[ISpline + 8]<<'\t'<<Params[ISpline + 9]<<endl;
  }

// ******************************************************************
// u1 component
  for(i=1;i<N_Splines;i++)
   {
     u0 = u1_spl[i-1];
     u1 = u1_spl[i];
     u2 = u1_spl[i+1];

     u1rhs[i] = 6./(h[i]+h[i+1])*((u2-u1)/h[i+1]-(u1-u0)/h[i]);
    }

   u1rhs[0] = u1rhs[1];
   u1rhs[N_Splines] = u1rhs[N_Splines-1];

  Solver_3dia(N_Splines, a, b, c, u1rhs, Mu1);


// u2 component
  for(i=1;i<N_Splines;i++)
   {
     u0 = u2_spl[i-1];
     u1 = u2_spl[i];
     u2 = u2_spl[i+1];

     u2rhs[i] = 6./(h[i]+h[i+1])*((u2-u1)/h[i+1]-(u1-u0)/h[i]);
    }

   u2rhs[0] = u2rhs[1];
   u2rhs[N_Splines] = u2rhs[N_Splines-1];

  Solver_3dia(N_Splines, a, b, c, u2rhs, Mu2);

// surfactant
  for(i=1;i<N_Splines;i++)
   {
     u0 = surf_spl[i-1];
     u1 = surf_spl[i];
     u2 = surf_spl[i+1];

     srhs[i] = 6./(h[i]+h[i+1])*((u2-u1)/h[i+1]-(u1-u0)/h[i]);
    }

   srhs[0] = srhs[1];
   srhs[N_Splines] = srhs[N_Splines-1];

  Solver_3dia(N_Splines, a, b, c, srhs, Msurf);



  for(i=0;i<N_Splines;i++)
  {
    ISpline = i*6;

    FEParams[ISpline  ] = -Mu1[i]*h[i+1]*h[i+1]/2. +
                          ((u1_spl[i+1]-u1_spl[i])/h[i+1]-h[i+1]/6.*(Mu1[i+1]-Mu1[i]))*h[i+1];
    FEParams[ISpline + 1] = Mu1[i+1]*h[i+1]*h[i+1]/2. +
                          ((u1_spl[i+1]-u1_spl[i])/h[i+1]-h[i+1]/6.*(Mu1[i+1]-Mu1[i]))*h[i+1];


    FEParams[ISpline + 2  ] = -Mu2[i]*h[i+1]*h[i+1]/2. +
                          ((u2_spl[i+1]-u2_spl[i])/h[i+1]-h[i+1]/6.*(Mu2[i+1]-Mu2[i]))*h[i+1];
    FEParams[ISpline + 3] = Mu2[i+1]*h[i+1]*h[i+1]/2. +
                          ((u2_spl[i+1]-u2_spl[i])/h[i+1]-h[i+1]/6.*(Mu2[i+1]-Mu2[i]))*h[i+1];


    FEParams[ISpline + 4  ] = -Msurf[i]*h[i+1]*h[i+1]/2. +
                          ((surf_spl[i+1]-surf_spl[i])/h[i+1]-h[i+1]/6.*(Msurf[i+1]-Msurf[i]))*h[i+1];
    FEParams[ISpline + 5] = Msurf[i+1]*h[i+1]*h[i+1]/2. +
                          ((surf_spl[i+1]-surf_spl[i])/h[i+1]-h[i+1]/6.*(Msurf[i+1]-Msurf[i]))*h[i+1];

  }


// ******************************************************************


   teta = 1.0/N_Splines;
   T = 0;

   Param9[0] = 0;
   for(i=1;i<=N_Splines;i++) 
    Param9[i] = Params[(i-1)*10+8];

   m = 0;
   m1 = 0;
//    wetting points fe values noneed to set
   for(j=0;j<N_E;j++)
    {
     T = double(m)*teta;
     for(i=1;i<=N_Splines;i++)
      {
       ISpline = (i-1)*10;
       USpline = (i-1)*6;
       FeDof   = i-1;
       if((T>=Param9[i-1]) && (T<=Param9[i]))
        {
      // further T must be from [0;1] on a subspline
         T = (T-Param9[i-1])/(Param9[i]-Param9[i-1]);
         break;
        }
      }

   phi1 = (2.*T*T - 3.*T)*T + 1.;
   phi2 = (-2.*T + 3.)*T*T;
   phi3 = (T*T - 2.*T + 1.)*T;
   phi4 = (T - 1)*T*T;

   X = Params[ISpline    ]*phi1 + Params[ISpline + 2]*phi2 +
       Params[ISpline + 4]*phi3 + Params[ISpline + 6]*phi4;
   Y = Params[ISpline + 1]*phi1 + Params[ISpline + 3]*phi2 +
       Params[ISpline + 5]*phi3 + Params[ISpline + 7]*phi4;

//     if(Y < 0 || fabs(Y)<1e-8) Y = 0.0; // no penetration on solid boundary
    cell[j]->GetVertex(EdgeNo[j])->SetCoords(X, Y);

//     OutPut("NewX:"<<' '<< m <<' '<<X<<' '<< Y<<endl);
    m++;

// **************************************************************************************
// for fe values

     u0 = u1_spl[FeDof]*phi1 + u1_spl[FeDof+1]*phi2 +
              FEParams[USpline]*phi3 + FEParams[USpline + 1]*phi4;
     u1 = u2_spl[FeDof]*phi1 + u2_spl[FeDof+1]*phi2 +
              FEParams[USpline+2]*phi3 + FEParams[USpline + 3]*phi4;

    surf = surf_spl[FeDof]*phi1 + surf_spl[FeDof+1]*phi2 +
              FEParams[USpline+4]*phi3 + FEParams[USpline + 5]*phi4;

// //    if(Y == 0.0) u1 = 0.;// no penetration on solid boundary
    if(j!=0) // endpoints no need to set
     {
      ValuesUX[U_DOF[m1]] = u0;
      ValuesUY[U_DOF[m1]] = u1;
      Surfact[Surf_DOF[m1]] = surf;
     }
// //    if(fabs(u1)<1e-8)
// //    OutPut("NewU:"<<' '<< X <<' '<<Y<<' '<< u0<<' '<< u1<<endl);

  m1++;
// **************************************************************************************

    Joint = cell[j]->GetJoint(EdgeNo[j]);
    isojoint = (TIsoBoundEdge *)Joint;
    k = isojoint->GetN_Vertices();
    if(k==ORDER-1)
     {
      IsoVertices = isojoint->GetVertices();
      for(i3=0;i3<k;i3++)
       {
       T = double(m)*teta;
       for(i=1;i<=N_Splines;i++)
        {
         ISpline = (i-1)*10;
         USpline = (i-1)*6;
         FeDof   = i-1;
         if((T>=Param9[i-1]) && (T<=Param9[i]))
          {
      // further T must be from [0;1] on a subspline
//          cout<< ISpline << ' ' << T;
           T = (T-Param9[i-1])/(Param9[i]-Param9[i-1]);
//          cout<< ' ' << T <<endl;
          break;
         }
       }
  
     phi1 = (2.*T*T - 3.*T)*T + 1.;
     phi2 = (-2.*T + 3.)*T*T;
     phi3 = (T*T - 2.*T + 1.)*T;
     phi4 = (T - 1)*T*T;

     X = Params[ISpline    ]*phi1 + Params[ISpline + 2]*phi2 +
         Params[ISpline + 4]*phi3 + Params[ISpline + 6]*phi4;
     Y = Params[ISpline + 1]*phi1 + Params[ISpline + 3]*phi2 +
         Params[ISpline + 5]*phi3 + Params[ISpline + 7]*phi4;     

//      if(Y < 0 || fabs(Y)<1e-8) Y = 1e-3; // no penetration on solid boundary
        IsoVertices[i3]->SetCoords(X, Y);
//         OutPut("NewX:"<<' '<< m <<' '<<X<<' '<< Y<<endl);
        m++;

// **************************************************************************************
// for fe values

     u0 = u1_spl[FeDof]*phi1 + u1_spl[FeDof+1]*phi2 +
              FEParams[USpline]*phi3 + FEParams[USpline + 1]*phi4;
     u1 = u2_spl[FeDof]*phi1 + u2_spl[FeDof+1]*phi2 +
              FEParams[USpline+2]*phi3 + FEParams[USpline + 3]*phi4;

    surf = surf_spl[FeDof]*phi1 + surf_spl[FeDof+1]*phi2 +
              FEParams[USpline+4]*phi3 + FEParams[USpline + 5]*phi4;

// //    if(Y == 0.0) u1 = 0.;// no penetration on solid boundary

   ValuesUX[U_DOF[m1]] = u0;
   ValuesUY[U_DOF[m1]] = u1;
   Surfact[Surf_DOF[m1]] = surf;

 
// //    if(fabs(u1)<1e-8)
// //    OutPut("NewU:"<<' '<< X <<' '<<Y<<' '<< u0<<' '<< u1<<endl);

  m1++;
// **************************************************************************************


       }
     }
   }  //  for(j=0;j<N_E

// cout << "Surfact[Surf_DOF[0]] " <<Surfact[Surf_DOF[0]] << " Surfact[Surf_DOF[m1]] " << Surfact[Surf_DOF[m1]] <<endl;


   delete [] h; delete [] t; delete [] a; delete [] b;
   delete [] c; delete [] rhs; delete [] Mx; delete [] My;
   delete [] Params; delete [] Param9;  delete [] x; delete [] y;
   delete [] U_DOF;  delete [] Surf_DOF;
   delete []  u1rhs ;
   delete []  u2rhs;
   delete []  srhs;
   delete []  u1_spl;
   delete []  u2_spl;
   delete []  surf_spl;
   delete []  Mu1;
   delete []  Mu2;
    delete [] Msurf;
//    exit(0);
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
  double Surf_Mass[2], errors[7], p1, p2;
  double t1, t2, res, res2, oldres, solver_time, solver_time_curr, residual, oldresidual;
  double impuls_residual,limit,linredfac, total_time, t3, t4;
  int N_LinIter, N_LinIterCurr, N_LinIterCurrIte, N_SubSteps, N_Active, n_aux;
  double gamma, tau, oldtau, Surfgamma;
  int *RowPtr, N_SurfActive;
  int N_G, N_BoundaryNodes, VSP, n_change=0, N_S;
  double SavedRedFactor, wgt;
  double *vcoarse, *vfine; 
  double y_top[3], t_top[3], r_max;
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
  double FREX, FREY, SLPX, SLPY, x_max;
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
  int **N_array = new int*[5];
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
//   opts<<'Y'; // Supress adding vertices on boundary edges
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
    x = 0.0;     // Axial bound

    y = T_b;
// spherical harmonic of order 2 
// spherical harmonic of order 4
    phi = atan2(y, x);
    if(mode==2)
      r = 1.0 + deviation*(sqrt(5./Pi)*(3.*cos(phi+Pi/2.)*cos(phi+Pi/2.) -1. )/4.);
    else if(mode==4) 
     {
      temp = cos(phi+Pi/2.);
      r = 1.0 + deviation*(sqrt(1./Pi)*(35.*temp*temp*temp*temp - 30.*temp*temp +3. )*3./16.);
     }
    else 
      {
        OutPut("No. of mode is taken as 0 check main programme"<<endl);
       exit(0);
       }

    y = r*sin(phi);
    y_begin = y;



    y = -T_b;
    phi = atan2(y, x);
// spherical harmonic of order 2
// spherical harmonic of order 4
    phi = atan2(y, x);
    if(mode==2)
      r = 1.0 + deviation*(sqrt(5./Pi)*(3.*cos(phi+Pi/2.)*cos(phi+Pi/2.) -1. )/4.);
    else if(mode==4) 
     {
      temp = cos(phi+Pi/2.);
      r = 1.0 + deviation*(sqrt(1./Pi)*(35.*temp*temp*temp*temp - 30.*temp*temp +3. )*3./16.);
     }

    y = r*sin(phi);

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
    phi = atan2(y, x);
    if(mode==2)
      r = 1.0 + deviation*(sqrt(5./Pi)*(3.*cos(phi+Pi/2.)*cos(phi+Pi/2.) -1. )/4.);
    else if(mode==4) 
     {
      temp = cos(phi+Pi/2.);
      r = 1.0 + deviation*(sqrt(1./Pi)*(35.*temp*temp*temp*temp - 30.*temp*temp +3. )*3./16.);
     }

      x = r*cos(phi);
      y = r*sin(phi);
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
 UpdateFreeBound ->SetParams( 0.0, 0.0, 1.0, 1.0, -Pi/2., Pi/2.);
 
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
  for(i=0;i<3;i++) N_array[i] = new int[mg_level];
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
    N_List[0][i] = new int[N_Cells];;
    N_List[1][i] = new int[N_Cells];;
    Domain2DSurface(Domain, SurfDomain, N_List[0], N_List[1], i);
//     N_List[0] = Cell_No_array;
//     N_List[1] = Joint_No_array;

    Surf_Coll = SurfDomain->GetCollection(It_Finest, 0);
    N_SCells= Surf_Coll->GetN_Cells();
    FE1D_List = new FE1D[N_SCells];
    for(j=0;j<N_SCells;j++)
     FE1D_List[j] = FE1D(TDatabase::ParamDB->ANSATZ_ORDER);
//      FE1D_List[j] = FE1D(TDatabase::ParamDB->ANSATZ_ORDER);

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
                                 TDatabase::ParamDB->VELOCITY_SPACE,   NULL);
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

    cout << "N_Inner in Surfact Space: " << surfact_space->GetActiveBound() << endl;
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
      N_Cells = coll->GetN_Cells();
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
      N_array[4] = new int [N_FreeJoints];
//       Edge_No = N_array[3];

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
            N_array[4][m1] = j;
            N_array[3][m1] = l;
            m1++;
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
   SortIsoVertices_axial3D(Free_Cells, FreeBound_Vert, N_array[3],  N_array[4], N_FreeJoints);
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

    u1->Interpolate(InitialU1);
    u2->Interpolate(InitialU2);
    p->Interpolate(InitialP);



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

    GridXDot->Interpolate(GridU1);
    GridYDot->Interpolate(GridU2);


// surfactant

    srhs = new double[N_S];
    memset(srhs, 0, N_S*SizeOfDouble);
    Rhsarray[3][i] = srhs;

    ssol = new double[N_S];
    Sol[6][i] = ssol;
    memset(ssol, 0, N_S*SizeOfDouble);
    SurfSurfactant = new TFEFunction1D(SurfFeSps_Lev[0][i], SurfString, SurfString,
                                       ssol, N_S);
    SurfSurfactant->Interpolate(InitialS);
    SArrays[0][i] = SurfSurfactant;

// //  surfacetant in entire 2D domain
    surfact = new double[N_SO];
    Sol[7][i] = surfact;
    memset(surfact, 0, N_SO*SizeOfDouble);
    UPArrays[5][i] = new TFEFunction2D(FeSps_Lev[5][i], SurfString, SurfString, 
                                       surfact, N_SO);

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
  rhs =  Rhsarray[3][i];
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
      Output->AddFEFunction(p);
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
    OutPut(setw(20)<<"Time, Wett Len d : " << TDatabase::TimeDB->CURRENTTIME
                   <<"   "<< Rx-Lx<<endl);
    OutPut(setw(20)<<"Time, Volume : " << TDatabase::TimeDB->CURRENTTIME
                   <<"   "<< CurrVolume<<endl);
    OutPut(setw(20)<<"Time, Volume Diff : "<< TDatabase::TimeDB->CURRENTTIME
                   <<"   "<< CurrVolume - InitVolume << endl);

    min = 1e5;
    tx  = 0.;

    for(k=0;k<N_MovVert[1];k++) // no need to set end vertices again
     {
      MovBoundVert[1][k]->GetCoords(x1, y1);
//         cout << x1 << " " <<  y1<< endl;
      if(fabs(y1)<min)
       {
        tx = x1;
        min = fabs(y1);
       }
     }

     MovBoundVert[0][0]->GetCoords(x1, y1);


     OutPut(setw(25)<<"Time, Right tip, Top tip: "<< TDatabase::TimeDB->CURRENTTIME
                   <<"   "<< tx <<"   "<< y1 << endl);

  MovBoundVert[1][N_MovVert[1]-1]->GetCoords(x1, y1);
  y_top[0] = y_top[1] = y_top[2] = y1;
  t_top[0] = t_top[1] = t_top[2] = 0.;


        PrintSurfSurfactant(Free_Cells, MovBoundVert[1], N_array[3], 
                            N_MovVert[1], UPArrays[5][mg_level-1], N_BData);

      GetSurfactMass(UPArrays[5][mg_level-1], SArrays[0][mg_level-1], 
                     N_List[0][mg_level-1], N_List[1][mg_level-1], Surf_Mass);
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

   // **********************************************************************
        // extrapolate solution to get starting value of new time step
        // from the previous time step solution Sol[0]
        // solution of last time step sol_timestep_m1

        // save current sol to the previous iteration sol
       memcpy(Sol[1][mg_level-1], Sol[0][mg_level-1], SizeOfDouble*N_Unknowns);
        // save current solution to the previous time step sol
       memcpy(Sol[5][0], Sol[1][mg_level-1], SizeOfDouble*N_Unknowns);

        //  (take solution of previous discrete time as start solution)

//        if(TDatabase::TimeDB->EXTRAPOLATE_WEIGHT!=0)
//         {
//          tau2 = TDatabase::TimeDB->EXTRAPOLATE_WEIGHT*tau/oldtau;
//          cout<< "tau2" <<tau2 <<endl;
//          tau1 = 1 + tau2;
//          // sol := tau1 *sol - tau2 * sol_timestep_m1
//          // at first time step: sol = sol_timestep_m1 -> result is sol
//          for (k=0;k<2*N_U;k++)
//           Sol[0][mg_level-1][k] = tau1*Sol[0][mg_level-1][k] - tau2*Sol[5][0][k];
//         }

       if (TDatabase::ParamDB->SOLVER_TYPE==GMG)
           TMG[0]->RestrictToAllGrids();

      // **********************************************************************
      // working array for rhs is B, initialize B

       memset(Rhsarray[1][mg_level-1], 0, N_Unknowns*SizeOfDouble);

      // =============================================================
      // determine velocity of grid movement on all levels
      // restrict grid movement to coarser grid
      // =============================================================


      GetGridVelocity_B2_axial3D(Entries, Sol[4][mg_level-1], Rhsarray[2][mg_level-1],
                         GridKCol, GridRowPtr,
                         GridPos, AuxGridPos,
                         VeloVect[0][mg_level-1], tau,
                         VeloVect[1][mg_level-1]);

      for(i=mg_level-1;i>0;i--)
      {
        vcoarse = UPArrays[3][i-1]->GetValues();
        vfine = UPArrays[3][i]->GetValues();
        RestrictFunction(FeSps_Lev[2][i-1], FeSps_Lev[2][i],
                         vcoarse, vfine, tmp);
        vcoarse = UPArrays[4][i-1]->GetValues();
        vfine = UPArrays[4][i]->GetValues();
        RestrictFunction(FeSps_Lev[2][i-1], FeSps_Lev[2][i],
                         vcoarse, vfine, tmp);
      } // endfor i


   for(i=0;i<mg_level;i++)
    {
    // set discrete forms
    if ((mg_type==1) && (i<mg_level-1))
    {
      DiscreteForm = DiscreteFormUpwind;

      CurrentDiscType =  UPWIND;
    }
    else
      switch(TDatabase::ParamDB->DISCTYPE)
      {
        case GALERKIN:
          DiscreteForm = DiscreteFormGalerkin;

          CurrentDiscType =  GALERKIN;
          break;

        case UPWIND:
          DiscreteForm = DiscreteFormUpwind;
          CurrentDiscType =  UPWIND;
          break;

      default:
        OutPut("Unknown DISCTYPE" << endl);
        exit(1);
    }

     // set matrices
    switch(TDatabase::ParamDB->NSTYPE)
    {
      case 1:
        SQMATRICES[0] = SqMat[12][i];
        SQMATRICES[1] = SqMat[13][i];
        MATRICES[0] = Mat[0][i];
        MATRICES[1] = Mat[1][i];

        SQMATRICES[0]->Reset();
        SQMATRICES[1]->Reset();
        MATRICES[0]->Reset();
        MATRICES[1]->Reset();

        N_SquareMatrices = 2;
        N_RectMatrices = 2;

      break;

      case 2:
        SQMATRICES[0] = SqMat[12][i];
        SQMATRICES[1] = SqMat[13][i];
        MATRICES[0] = Mat[0][i];
        MATRICES[1] = Mat[1][i];
        MATRICES[2] = Mat[2][i];
        MATRICES[3] = Mat[3][i];

        SQMATRICES[0]->Reset();
        SQMATRICES[1]->Reset();
        MATRICES[0]->Reset();
        MATRICES[1]->Reset();
        MATRICES[2]->Reset();
        MATRICES[3]->Reset();

        N_SquareMatrices = 2;
        N_RectMatrices = 4;

     break;

      case 3:
        SQMATRICES[0] = SqMat[0][i];
        SQMATRICES[1] = SqMat[1][i];
        SQMATRICES[2] = SqMat[2][i];
        SQMATRICES[3] = SqMat[3][i];
        SQMATRICES[4] = SqMat[4][i];
        SQMATRICES[5] = SqMat[7][i];
        MATRICES[0] = Mat[0][i];
        MATRICES[1] = Mat[1][i];

        SQMATRICES[0]->Reset();
        SQMATRICES[1]->Reset();
        SQMATRICES[2]->Reset();
        SQMATRICES[3]->Reset();
        SQMATRICES[4]->Reset();
        SQMATRICES[5]->Reset();
        MATRICES[0]->Reset();
        MATRICES[1]->Reset();

        N_SquareMatrices = 6;
        N_RectMatrices = 2;

     break;

      case 4:
        SQMATRICES[0] = SqMat[0][i];
        SQMATRICES[1] = SqMat[1][i];
        SQMATRICES[2] = SqMat[2][i];
        SQMATRICES[3] = SqMat[3][i];
        SQMATRICES[4] = SqMat[4][i];
        SQMATRICES[5] = SqMat[7][i];

        SQMATRICES[6] = SqMat[5][i];
        SQMATRICES[7] = SqMat[6][i];

        MATRICES[0] = Mat[0][i];
        MATRICES[1] = Mat[1][i];
        MATRICES[2] = Mat[2][i];
        MATRICES[3] = Mat[3][i];

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
      break;
    }

       // parameters which are the same for all NSTYPEs

        N_Rhs = 2;
        N_FESpaces = 3;

        fesp[0] = FeSps_Lev[0][i];
        fesp[1] = FeSps_Lev[1][i];
        fesp[2] = FeSps_Lev[2][i];

        fefct[0] = UPArrays[0][i];
        fefct[1] = UPArrays[1][i];
        fefct[2] = UPArrays[3][i];
        fefct[3] = UPArrays[4][i];

        ferhs[0] = FeSps_Lev[0][i];
        ferhs[1] = FeSps_Lev[0][i];

        rhs = Rhsarray[0][i];
        RHSs[0] = rhs;
        RHSs[1] = rhs + N_array[0][i];
        RHSs[2] = rhs + 2*N_array[0][i];

        memset(rhs, 0, (2*N_array[0][i]+N_array[1][i])*SizeOfDouble);

        // 4 parameters are needed for assembling (u1_old, u2_old)
//         aux =  new TAuxParam2D(MovingNSN_FESpaces4, MovingNSN_Fct4,
//                                MovingNSN_ParamFct4,
//                                MovingNSN_FEValues4,
//                                fesp, fefct,
//                                MovingNSFct4,
//                                MovingNSFEFctIndex4, MovingNSFEMultiIndex4,
//                                MovingNSN_Params4, MovingNSBeginParam4);

        aux =  new TAuxParam2D(MovingNSN_FESpaces4_axial3D, MovingNSN_Fct4_axial3D,
                               MovingNSN_ParamFct4_axial3D,
                               MovingNSN_FEValues4_axial3D,
                               fesp, fefct,
                               MovingNSFct4_axial3D,
                               MovingNSFEFctIndex4_axial3D, MovingNSFEMultiIndex4_axial3D,
                               MovingNSN_Params4_axial3D, MovingNSBeginParam4_axial3D);

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

     SqMat[14][i]->Reset(); // Matrix entries for freesurf int;
     SqMat[15][i]->Reset(); // no need to calculate in nonlinear steps

      if(i==mg_level-1)
       {
//         FreeSurfInt_BC2(SqMat[14][i], SqMat[15][i],
//                         Rhsarray[0][i], Rhsarray[0][i]+N_array[0][i],
//                         BoundCondition, tau);
//         FreeSurf_axial3D(SqMat[14][i], SqMat[15][i],
//                          Rhsarray[0][i], Rhsarray[0][i]+N_array[0][i],
//                          BoundCondition, tau);


        MapSurfToDomain(SArrays[0][mg_level-1], UPArrays[5][mg_level-1], 
                        N_List[0][mg_level-1], N_List[1][mg_level-1]);

        FreeSurf_axial3D_Surfact(SqMat[14][i], SqMat[15][i],
                         Rhsarray[0][i], Rhsarray[0][i]+N_array[0][i],
                         BoundCondition, tau, UPArrays[5][mg_level-1]);

       }
      else
       {
        cout << " FreeSurfInt " <<endl;
        exit(0);
        FreeSurfInt(SqMat[14][i], SQMATRICES[1],
                    SQMATRICES[2], SqMat[15][i],
                    Rhsarray[0][i], RHSs[1],
                    BoundCondition, tau,
                    TDatabase::ParamDB->P7);
      }

     // Adding freesurf entries to A11 and A22
     MatAdd(SqMat[0][i], SqMat[14][i], 1);
     MatAdd(SqMat[3][i], SqMat[15][i], 1);

    // add convective term in the upwind discretizations
    if(DiscreteForm == DiscreteFormUpwind)
    {
      switch(TDatabase::ParamDB->NSTYPE)
      {
        case 1:
        case 2:
          // do upwinding with one matrix
          UpwindForNavierStokes(SQMATRICES[0], UPArrays[0][i], UPArrays[1][i]);
          //cout << "UPWINDING DONE : level " << i << endl;
        break;

        case 3:
        case 4:
          // do upwinding with two matrices
          UpwindForNavierStokes(SQMATRICES[0], UPArrays[0][i], UPArrays[1][i]);
          UpwindForNavierStokes(SQMATRICES[3], UPArrays[0][i], UPArrays[1][i]);
         // cout << "UPWINDING DONE(2) : level " << i << endl;
        break;
      } // endswitch
    }

    // set rows of Dirichlet dof in off diagonal matrix blocks
    // to zero
    switch(TDatabase::ParamDB->NSTYPE)
    {
      case 3:
      case 4:
        N_Active = FeSps_Lev[0][i]->GetActiveBound();
        // get row in off diagonal matrix where the Dirichlet nodes start
        RowPtr = SqMat[1][i]->GetRowPtr();
        // compute number of entries starting from this row to the end
        // of the matrix
        j = RowPtr[N_Active];
        k = RowPtr[N_array[0][i]]-j;
        // get number of active dof
        // set these entries to zero
        memset(SqMat[1][i]->GetEntries()+j, 0, SizeOfDouble*k);
        memset(SqMat[2][i]->GetEntries()+j, 0, SizeOfDouble*k);
      break;
     }
    } // end for i

      if (very_first_time==1)
        {
          very_first_time=0;
          l--;
          continue;
        }

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
          for(i=0;i<mg_level;i++)
          {
            SQMATRICES[0] = SqMat[0][i];
            SQMATRICES[1] = SqMat[3][i];
            SQMATRICES[2] = SqMat[1][i];
            SQMATRICES[3] = SqMat[2][i];
            SQMATRICES[4] = SqMat[4][i];
            SQMATRICES[5] = SqMat[7][i];
            SQMATRICES[6] = SqMat[5][i];
            SQMATRICES[7] = SqMat[6][i];

            MATRICES[0] = Mat[2][i];
            MATRICES[1] = Mat[3][i];

            fesp[0] = FeSps_Lev[0][i];
            ferhs[0] = FeSps_Lev[0][i];
            ferhs[1] = FeSps_Lev[0][i];

            RHSs[0] = Rhsarray[0][i];
            RHSs[1] = Rhsarray[0][i]+N_array[0][i];

            aux = new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);

            Assemble2DSlipBC(N_FESpaces, fesp,
                             N_SquareMatrices, SQMATRICES,
                             N_RectMatrices, MATRICES,
                             N_Rhs, RHSs, ferhs,
                             DiscreteForm,
                             BoundaryConditions,
                             BoundValues,
                             aux, UPArrays[0][i],
                             UPArrays[1][i]);

       TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION = 1;
       delete aux;
        }
       }

        for(i=0;i<mg_level;i++)
        {
          switch(TDatabase::ParamDB->NSTYPE)
          {
            case 1:
            case 3:
              Dscal(Mat[0][i]->GetN_Entries(), tau,
                    Mat[0][i]->GetEntries());
              Dscal(Mat[1][i]->GetN_Entries(), tau,
                    Mat[1][i]->GetEntries());
              break;

            case 2:
            case 4:
              Dscal(Mat[2][i]->GetN_Entries(), tau,
                    Mat[2][i]->GetEntries());
              Dscal(Mat[3][i]->GetN_Entries(), tau,
                    Mat[3][i]->GetEntries());
              Dscal(Mat[0][i]->GetN_Entries(), tau,
                    Mat[0][i]->GetEntries());
              Dscal(Mat[1][i]->GetN_Entries(), tau,
                    Mat[1][i]->GetEntries());
              break;
          }
        } // endfor

     gamma = 0;

     if(remeshed)
      {
        // update rhs from previous rhs
//         Daxpy(N_Active, tau*TDatabase::TimeDB->THETA3, Rhsarray[0][mg_level-1],
//               Rhsarray[1][mg_level-1]);
//         Daxpy(N_Active, tau*TDatabase::TimeDB->THETA3, Rhsarray[0][mg_level-1]+N_U,
//               Rhsarray[1][mg_level-1]+N_U);
        remeshed = FALSE;
      }

         // update rhs
//         Daxpy(N_Active, tau*TDatabase::TimeDB->THETA4, Rhsarray[0][mg_level-1],
//               Rhsarray[1][mg_level-1]);
//         Daxpy(N_Active, tau*TDatabase::TimeDB->THETA4, Rhsarray[0][mg_level-1]+N_U,
//               Rhsarray[1][mg_level-1]+N_U);

        Daxpy(N_Active, tau, Rhsarray[0][mg_level-1],
              Rhsarray[1][mg_level-1]);
        Daxpy(N_Active, tau, Rhsarray[0][mg_level-1]+N_U,
              Rhsarray[1][mg_level-1]+N_U);



         // update rhs by Laplacian and convective term initialy
        // by current time step
        // scaled by current sub time step length and theta2
        // currently : M := M + gamma A
        // M = M + (-gamma - tau*TDatabase::TimeDB->THETA2) A

       for(i=0;i<mg_level;i++)
        {
          switch(TDatabase::ParamDB->NSTYPE)
          {
            case 1:
            case 2:
              MatAdd(SqMat[13][i], SqMat[12][i],
                     -gamma - tau*TDatabase::TimeDB->THETA2);
              break;

            case 3:
            case 4:
              MatAdd(SqMat[4][i], SqMat[0][i],
                     -gamma - tau*TDatabase::TimeDB->THETA2);
              MatAdd(SqMat[5][i], SqMat[1][i],
                     -gamma - tau*TDatabase::TimeDB->THETA2);
              MatAdd(SqMat[6][i], SqMat[2][i],
                     -gamma - tau*TDatabase::TimeDB->THETA2);
              MatAdd(SqMat[7][i], SqMat[3][i],
                     -gamma - tau*TDatabase::TimeDB->THETA2);
              break;
          } // endswitch
        }
        // set current factor of steady state matrix
          gamma = -tau*TDatabase::TimeDB->THETA2;

        // defect = M * Sol
        // B:= B + defect (rhs)
        switch(TDatabase::ParamDB->NSTYPE)
        {
          case 1:
          case 2:
            MatVectActive(SqMat[13][mg_level-1], Sol[5][0], defect);
            Daxpy(N_Active, 1, defect, Rhsarray[1][mg_level-1]);
            MatVectActive(SqMat[13][mg_level-1], Sol[5][0]+N_U, defect+N_U);
            Daxpy(N_Active, 1, defect+N_U, Rhsarray[1][mg_level-1]+N_U);
            break;

          case 3:
          case 4:
            MatVectActive(SqMat[4][mg_level-1], Sol[5][0], defect);
            Daxpy(N_Active, 1, defect, Rhsarray[1][mg_level-1]);
            MatVectActive(SqMat[5][mg_level-1], Sol[5][0]+N_U, defect);
            Daxpy(N_Active, 1, defect, Rhsarray[1][mg_level-1]);
            MatVectActive(SqMat[6][mg_level-1], Sol[5][0], defect+N_U);
            Daxpy(N_Active, 1, defect+N_U, Rhsarray[1][mg_level-1]+N_U);
            MatVectActive(SqMat[7][mg_level-1], Sol[5][0]+N_U, defect+N_U);
            Daxpy(N_Active, 1, defect+N_U, Rhsarray[1][mg_level-1]+N_U);
            break;
        }

    // set Dirichlet values
    // RHSs[0] still available from assembling
    memcpy(Rhsarray[1][mg_level-1]+N_Active, RHSs[0]+N_Active, (N_U-N_Active)*SizeOfDouble);
    memcpy(Rhsarray[1][mg_level-1]+N_Active+N_U, RHSs[1]+N_Active,(N_U-N_Active)*SizeOfDouble);

    // copy Dirichlet values from rhs into Sol[0][mg_level-1]
    memcpy(Sol[0][mg_level-1]+N_Active, RHSs[0]+N_Active, (N_U-N_Active)*SizeOfDouble);
    memcpy(Sol[0][mg_level-1]+N_Active+N_U, RHSs[1]+N_Active, (N_U-N_Active)*SizeOfDouble);


     //=====================================================================
     // the stiffness matrix is stored on M11, (M12, M21, M22)
     // assembling of system matrix
     //========================================================================

        // M = M + (-gamma + tau*TDatabase::TimeDB->THETA1) A
      for(i=0;i<mg_level;i++)
       {
        switch(TDatabase::ParamDB->NSTYPE)
          {
            case 1:
            case 2:
              MatAdd(SqMat[13][i], SqMat[12][i],
                     -gamma + tau*TDatabase::TimeDB->THETA1);
              break;

            case 3:
            case 4:
              MatAdd(SqMat[4][i], SqMat[0][i],
                     -gamma + tau*TDatabase::TimeDB->THETA1);
              MatAdd(SqMat[5][i], SqMat[1][i],
                     -gamma + tau*TDatabase::TimeDB->THETA1);
              MatAdd(SqMat[6][i], SqMat[2][i],
                     -gamma + tau*TDatabase::TimeDB->THETA1);
              MatAdd(SqMat[7][i], SqMat[3][i],
                     -gamma + tau*TDatabase::TimeDB->THETA1);
              break;
          } // endswitch
       }

       // set current factor of steady state matrix
        gamma = tau*TDatabase::TimeDB->THETA1;
      //======================================================================
        OutPut(endl << "CURRENT TIME: ");
        OutPut(TDatabase::TimeDB->CURRENTTIME << endl);

      //======================================================================
      // nonlinear loop
      //======================================================================
        N_LinIterCurr = 0;
        solver_time_curr = 0;

        for(j=0;j<Max_It;j++)
        {
          memcpy(Sol[1][mg_level-1], Sol[0][mg_level-1], SizeOfDouble*N_Unknowns);
          memset(defect, 0, N_Unknowns*SizeOfDouble);

          switch(TDatabase::ParamDB->NSTYPE)
          {
            case 1:
              SQMATRICES[0] = SqMat[13][mg_level-1];
              MATRICES[0] = Mat[0][mg_level-1];
              MATRICES[1] = Mat[1][mg_level-1];
             break;
            case 2:
              SQMATRICES[0] = SqMat[13][mg_level-1];
              MATRICES[0] = Mat[0][mg_level-1];
              MATRICES[1] = Mat[1][mg_level-1];
              MATRICES[2] = Mat[2][mg_level-1];
              MATRICES[3] = Mat[3][mg_level-1];
              break;
            case 3:
              SQMATRICES[0] = SqMat[4][mg_level-1];
              SQMATRICES[1] = SqMat[5][mg_level-1];
              SQMATRICES[2] = SqMat[6][mg_level-1];
              SQMATRICES[3] = SqMat[7][mg_level-1];
              MATRICES[0] = Mat[0][mg_level-1];
              MATRICES[1] = Mat[1][mg_level-1];
              break;
            case 4:
              SQMATRICES[0] = SqMat[4][mg_level-1];
              SQMATRICES[1] = SqMat[5][mg_level-1];
              SQMATRICES[2] = SqMat[6][mg_level-1];
              SQMATRICES[3] = SqMat[7][mg_level-1];
              MATRICES[0] = Mat[0][mg_level-1];
              MATRICES[1] = Mat[1][mg_level-1];
              MATRICES[2] = Mat[2][mg_level-1];
              MATRICES[3] = Mat[3][mg_level-1];
              break;
          }
          // compute defect

          if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
           IntoL20FEFunction(Sol[0][mg_level-1]+2*N_U, N_P,FeSps_Lev[1][mg_level-1],
                        velocity_space_code, pressure_space_code);

          Defect(sqmatrices,matrices,Sol[0][mg_level-1],Rhsarray[1][mg_level-1],defect);

          if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE) 
           IntoL20Vector2D(defect+2*N_U, N_P,pressure_space_code);


//         for(k=0;k<N_array[1][mg_level-1];k+=3)
//          {
//           OutPut(" defect: " << setw(5) << k << setw(20) << defect[k+2*N_U] << endl);
//          }

          residual =  Ddot(N_Unknowns, defect, defect);
          impuls_residual = Ddot(2*N_U, defect, defect);
          OutPut("nonlinear step " << setw(3) << j);
          OutPut(setw(14) << impuls_residual);
          OutPut(setw(14) << Ddot(N_P,defect+2*N_U,defect+2*N_U));
          OutPut(setw(14) << sqrt(residual));
          if (j>0)
          {
            OutPut(setw(14) << sqrt(residual)/oldresidual << endl);
          }
          else
          {
            OutPut(endl);
          }
          oldresidual = sqrt(residual);

          if ((((sqrt(residual)<=limit)||(j==Max_It-1)))
              && (j>=TDatabase::ParamDB->SC_MINIT))
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
            OutPut(TDatabase::TimeDB->CURRENTTIME << " total elapsed running time "
                                                  << total_time << endl);
            break;
          }

          //======================================================================
          // solve linear system
          //======================================================================
          switch(TDatabase::ParamDB->SOLVER_TYPE)
          {
            case AMG:
              TDatabase::ParamDB->SC_VERBOSE=0;
              t1 = GetTime();
              switch(TDatabase::ParamDB->NSTYPE)
              {
                case 1:
                  Solver(sqmatrixM, matrixB1, matrixB2, Rhsarray[1][mg_level-1],
                         Sol[0][mg_level-1]);
                  break;

                case 2:
                  Solver(sqmatrixM, matrixB1T, matrixB2T,
                         matrixB1, matrixB2, Rhsarray[1][mg_level-1], Sol[0][mg_level-1]);
              break;

                case 3:
                  Error("AMG does not work for NSTYPE = 3." << endl);
                  return -1;
                  break;

                case 4:
//                   Error("AMG does not work for NSTYPE = 4." << endl);
//                   return -1;
                  DirectSolver(SQMATRICES[0], SQMATRICES[1],
                               SQMATRICES[2], SQMATRICES[3],
                               MATRICES[2], MATRICES[3],
                               MATRICES[0], MATRICES[1],
                               Rhsarray[1][mg_level-1], Sol[0][mg_level-1]);

                break;
              }
              t2 = GetTime();
              solver_time_curr = t2-t1;
              solver_time += solver_time_curr;
              break;

            case GMG:

              t1 = GetTime();
              if (TDatabase::ParamDB->SC_PRECONDITIONER_SADDLE == 5)
              {
                memcpy(itmethod_sol, Sol[0][mg_level-1], N_Unknowns*SizeOfDouble);
                memcpy(itmethod_rhs, Rhsarray[1][mg_level-1], N_Unknowns*SizeOfDouble);
              }
              N_LinIterCurrIte = itmethod->Iterate(sqmatrices,matrices,
                                                   itmethod_sol,itmethod_rhs);
              N_LinIterCurr += N_LinIterCurrIte;
              N_LinIter += N_LinIterCurrIte;
              if (TDatabase::ParamDB->SC_PRECONDITIONER_SADDLE == 5)
              {
                memcpy(Sol[0][mg_level-1], itmethod_sol, N_Unknowns*SizeOfDouble);
                memcpy(Rhsarray[1][mg_level-1], itmethod_rhs, N_Unknowns*SizeOfDouble);
              }
              t2 = GetTime();
              solver_time_curr += t2-t1;
              solver_time += solver_time_curr;

//            if(omega!=1)
//              {
//               for(k=0;k<2*N_U;k++)
//               {
//                 p2 = Sol[0][mg_level-1][k]-Sol[1][mg_level-1][k];
//                 Sol[0][mg_level-1][k] = Sol[1][mg_level-1][k] + omega * p2;
// //                 Sol[0][mg_level-1][k] += omega * p2;
//               }
//             }
              break;
          } // endswitch SOLVER_TYPE
          //======================================================================
          // end solve linear system
          //======================================================================

         // restore mass matrices by subtracting the A-matrices
          for(i=0;i<mg_level;i++)
          {
            switch(TDatabase::ParamDB->NSTYPE)
            {
              case 1:
              case 2:
                MatAdd(SqMat[13][i], SqMat[12][i], -gamma);
                break;

              case 3:
              case 4:
                MatAdd(SqMat[4][i], SqMat[0][i], -gamma);
                MatAdd(SqMat[5][i], SqMat[1][i], -gamma);
                MatAdd(SqMat[6][i], SqMat[2][i], -gamma);
                MatAdd(SqMat[7][i], SqMat[3][i], -gamma);
                break;
            } // endswitch
          } // endfor i
          // set current factor of steady state matrix
         gamma = 0;

         if (TDatabase::ParamDB->SOLVER_TYPE==GMG)
            TMG[0]->RestrictToAllGrids();

         //======================================================================
         // assemble new matrix due to nonlinearity
         //======================================================================

         GetGridVelocity_B2_axial3D(Entries, Sol[4][mg_level-1], Rhsarray[2][mg_level-1],
                         GridKCol, GridRowPtr,
                         GridPos, AuxGridPos,
                         VeloVect[0][mg_level-1], tau,
                         VeloVect[1][mg_level-1]);

         for(i=mg_level-1;i>0;i--)
          {
           vcoarse = UPArrays[3][i-1]->GetValues();
           vfine = UPArrays[3][i]->GetValues();
           RestrictFunction(FeSps_Lev[2][i-1], FeSps_Lev[2][i],
                            vcoarse, vfine, tmp);
           vcoarse = UPArrays[4][i-1]->GetValues();
           vfine = UPArrays[4][i]->GetValues();
           RestrictFunction(FeSps_Lev[2][i-1], FeSps_Lev[2][i],
                              vcoarse, vfine, tmp);
          } // endfor i

          // for all levels
          for(i=0;i<mg_level;i++)
          {
            if ((mg_type==1) && (i<mg_level-1))
            {
              DiscreteForm = DiscreteFormNLUpwind;
              CurrentDiscType =  UPWIND;
            }
            else
              switch(TDatabase::ParamDB->DISCTYPE)
              {
                case GALERKIN:
                case SMAGORINSKY_EXPL:
                  DiscreteForm = DiscreteFormNLGalerkin;
                  CurrentDiscType =  GALERKIN;
                  break;

                case UPWIND:
                  DiscreteForm = DiscreteFormNLUpwind;
                  CurrentDiscType =  UPWIND;
                  break;

                default:
                  OutPut("Unknown DISCTYPE" << endl);
                  exit(1);
              }
            // set pointers to matrices
            switch(TDatabase::ParamDB->NSTYPE)
            {
              case 1:
              case 2:
                SQMATRICES[0] = SqMat[12][i];
                SQMATRICES[0]->Reset();

                N_SquareMatrices = 1;
                N_RectMatrices = 0;
                N_Rhs = 0;
                N_FESpaces = 1;
                break;

              case 3:
              case 4:
                N_RectMatrices = 0;

                N_Rhs = 0;
                N_FESpaces = 3;

                SQMATRICES[0] = SqMat[0][i];
                SQMATRICES[1] = SqMat[3][i];
                SQMATRICES[0]->Reset();
                SQMATRICES[1]->Reset();

                N_SquareMatrices = 2;
                last_sq = 1;

                break;
            }

            fesp[0] = FeSps_Lev[0][i];
            fesp[1] = FeSps_Lev[1][i];
            fesp[2] = FeSps_Lev[2][i];


            fefct[0] = UPArrays[0][i];
            fefct[1] = UPArrays[1][i];
            fefct[2] = UPArrays[3][i];
            fefct[3] = UPArrays[4][i];

            ferhs[0] = FeSps_Lev[0][i];
            ferhs[1] = FeSps_Lev[0][i];


            //======================================================================
            // assembling of matrices for each level due to nonlinearity
            // A_11, (A_22)
            // no assembling of rhs
            //======================================================================
//             aux =  new TAuxParam2D(MovingNSN_FESpaces4, MovingNSN_Fct4,
//                                    MovingNSN_ParamFct4,
//                                    MovingNSN_FEValues4,
//                                    fesp, fefct,
//                                    MovingNSFct4,
//                                    MovingNSFEFctIndex4, MovingNSFEMultiIndex4,
//                                    MovingNSN_Params4, MovingNSBeginParam4);


        aux =  new TAuxParam2D(MovingNSN_FESpaces4_axial3D, MovingNSN_Fct4_axial3D,
                               MovingNSN_ParamFct4_axial3D,
                               MovingNSN_FEValues4_axial3D,
                               fesp, fefct,
                               MovingNSFct4_axial3D,
                               MovingNSFEFctIndex4_axial3D, MovingNSFEMultiIndex4_axial3D,
                               MovingNSN_Params4_axial3D, MovingNSBeginParam4_axial3D);


            Assemble2D(N_FESpaces, fesp,
                       N_SquareMatrices, SQMATRICES,
                       N_RectMatrices, MATRICES,
                       N_Rhs, RHSs, ferhs,
                       DiscreteForm,
                       BoundaryConditions,
                       BoundValues,
                       aux);

           delete aux;

        // Adding freesurf entries to A11 and A22
         MatAdd(SqMat[0][i], SqMat[14][i], 1);
         MatAdd(SqMat[3][i], SqMat[15][i], 1);

            if(DiscreteForm == DiscreteFormNLUpwind)
            {
              switch(TDatabase::ParamDB->NSTYPE)
              {
                case 1:
                case 2:
                  // do upwinding with one matrix
                  UpwindForNavierStokes(SQMATRICES[0], UPArrays[0][i], UPArrays[1][i]);
                  //cout << "UPWINDING DONE : level " << i << endl;
                  break;

                case 3:
                case 4:
                  // do upwinding with two matrices
                  UpwindForNavierStokes(SQMATRICES[0], UPArrays[0][i], UPArrays[1][i]);
                  UpwindForNavierStokes(SQMATRICES[last_sq], UPArrays[0][i], UPArrays[1][i]);
                  //cout << "UPWINDING DONE(2) : level " << i << endl;
                  //cout << "check correct sqmatrix !!!! " << endl;
                  break;
              } // endswitch
            }


            // slip type bc detected, modify matrices accordingly
            if (TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION >= 1)
            {

              // prepare everything for the assembling of slip with friction bc
              // on level i
              N_FESpaces = 1;
              N_SquareMatrices = 2;
              N_RectMatrices = 0;
              N_Rhs = 2;
              DiscreteForm = NULL;

              SQMATRICES[0] = SqMat[0][i];
              SQMATRICES[1] = SqMat[3][i];

              MATRICES[0] = Mat[2][i];
              MATRICES[1] = Mat[3][i];
            
              fesp[0] = FeSps_Lev[0][i];
              ferhs[0] = FeSps_Lev[0][i];
              ferhs[1] = FeSps_Lev[0][i];

              RHSs[0] = Rhsarray[0][i];
              RHSs[1] = Rhsarray[0][i]+N_array[0][i];

              aux = new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);

              Assemble2DSlipBC(N_FESpaces, fesp,
                               N_SquareMatrices, SQMATRICES,
                               N_RectMatrices, MATRICES,
                               N_Rhs, RHSs, ferhs,
                               DiscreteForm,
                               BoundaryConditions,
                               BoundValues,
                               aux, UPArrays[0][i],
                               UPArrays[1][i]);
              delete aux;
             }

            //======================================================================
            // end of assemble new matrix due to nonlinearity
            //======================================================================

            // build stiffness matrix for next nonlinear iteration step
            // stiffness matrix (left upper block) is stored on
            // M11, (M12, M21, M22)
            // M = M +  tau*TDatabase::TimeDB->THETA1 A
            switch(TDatabase::ParamDB->NSTYPE)
            {
              case 1:
              case 2:
                MatAdd(SqMat[13][i], SqMat[12][i],
                       tau*TDatabase::TimeDB->THETA1);
                break;

              case 3:
              case 4:
                MatAdd(SqMat[4][i], SqMat[0][i],
                       tau*TDatabase::TimeDB->THETA1);
                MatAdd(SqMat[5][i], SqMat[1][i],
                       tau*TDatabase::TimeDB->THETA1);
                MatAdd(SqMat[6][i], SqMat[2][i],
                       tau*TDatabase::TimeDB->THETA1);
                MatAdd(SqMat[7][i], SqMat[3][i],
                       tau*TDatabase::TimeDB->THETA1);
                break;
            }
      } // endfor i
          // set current factor of steady state matrix
          gamma = tau*TDatabase::TimeDB->THETA1;

    } // endfor Max_It (solution of nonlinear equation)

//======================================================================
// end of nonlinear loop
//======================================================================


//======================================================================
// solve the surfactant equation begin
//======================================================================


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
  rhs =  Rhsarray[3][i];
  RHSs[0] = rhs;

  AssembleSurf1D(N_FESpaces, fesp, fefct, N_FESpaces_low, fesp_low,
                 N_SquareMatrices, SQMATRICES_SURF, N_Rhs, RHSs, 
                 ferhs_low, N_List[0][mg_level-1], 
                 N_List[1][mg_level-1]);


 // working array for rhs is SB, initialize SB
   memset(SB, 0, N_S*SizeOfDouble);
   // old rhs multiplied with current subtime step and theta3 on SB
   //   N_SurfActive= N_S on fespace on a 1D surface
//    Daxpy(N_S, tau*TDatabase::TimeDB->THETA3, srhs, SB);

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

   DirectSolver(MatricesS_M[mg_level-1], SB, ssol);

//======================================================================
// solve the surfactant equation end
//======================================================================

      MoveGrid_BC2_axial3D(Entries, Sol[4][mg_level-1], Rhsarray[2][mg_level-1],
                   GridKCol, GridRowPtr, GridPos,
                   VeloVect[0][mg_level-1], tau,
                   AuxGridPos);

//     memcpy(pos, auxpos, 2*N_G*SizeOfDouble);
       GridPos->GridToData();
       RefGridPos->GridToData();

  if((l==0) && ((m % 5) == 0))
   {
     MovBoundVert[0][0]->GetCoords(x1, y1);
     y_top[0] = y_top[1];
     y_top[1] = y_top[2];
     y_top[2] = y1;

     t_top[0] = t_top[1];
     t_top[1] = t_top[2];
     t_top[2] = TDatabase::TimeDB->CURRENTTIME;

    if(((y_top[2]-y_top[1])*(y_top[1]-y_top[0])) < 0)
      {
     if(t_top[2]>5e-2)
      {
       n_change ++;
       r_max = sqrt(x1*x1 + y1*y1);
       for(k=0;k<N_MovVert[1];k++) // no need to set end vertices again
        {
         MovBoundVert[1][k]->GetCoords(x1, y1);
//          cout << x1 << " " <<  y1<< endl;
         if(sqrt(x1*x1 + y1*y1)> r_max)
          r_max = sqrt(x1*x1 + y1*y1);
        }

       OutPut("time "<< TDatabase::TimeDB->CURRENTTIME<<
              " No.time change " << n_change << " top pos" <<
                y_top[1] << " r max " << r_max<<
		  " Y_TopOld " << y_top[0] << " Y_Top " << y_top[2] <<endl);
       }
      }
     }

    // Updating the Quard points on the solid surface
    MovBoundVert[0][0]->GetCoords(x, y);
    MovBoundVert[0][N_MovVert[0]-1]->GetCoords(SLPX, SLPY);
    // Bound startx, starty, x length and y length
//     cout << " x " <<x << " SLPX " << SLPX <<  " y " <<   y<< "  SLPY" <<   SLPY<< endl;
    UpdateSlipBound->SetParams(x, y, SLPX-x, SLPY-y);

    for(k=0;k<N_MovVert[0]-1;k++)
    {
     Slip_Joint[0][k]->UpdateParameters(MovBoundVert[0][k], MovBoundVert[0][k+1]);
    }

   if((l==0) && ((m % 1) == 0))
    {
     Getcellangle(FeSps_Lev[0][mg_level-1], Angle);
     OutPut( "MinAngle : "<< Angle[0]<< "  MaxAngle : "<<Angle[1]<< endl);
    }

   // **************************************************************************************
   // Remeshing Begin
   // **************************************************************************************
   if((Angle[0]<10.0) ||(Angle[1]>160.0) ) // || (Ly <= 1e-8) || (Ry <= 1e-8)
    {

      t1 = GetTime();
      MapSurfToDomain(SArrays[0][mg_level-1], UPArrays[5][mg_level-1], 
                       N_List[0][mg_level-1], N_List[1][mg_level-1]);
// OutPut( "Remeshing and Interpolation to be implemented " <<endl);
//     if(TDatabase::ParamDB->WRITE_VTK)
// //      {
//        MapSurfToDomain(SArrays[0][mg_level-1], UPArrays[5][mg_level-1], 
//                        N_List[0][mg_level-1], N_List[1][mg_level-1]);
//        os.seekp(std::ios::beg);
//        if(img<10) os << VtkBaseName<<".0000"<<img<<".vtk" << ends;
//        else if(img<100) os << VtkBaseName<<".000"<<img<<".vtk" << ends;
//        else if(img<1000) os << VtkBaseName<<".00"<<img<<".vtk" << ends;
//        else if(img<10000) os << VtkBaseName<<".0"<<img<<".vtk" << ends;
//        else  os << VtkBaseName<<"."<<img<<".vtk" << ends;
//        Output->WriteVtk(os.str().c_str());
//        img++;
// //      }
// exit(0);

      Remesh2D_axial3D_Surf(Domain, SurfDomain, Free_Cells, MovBoundVert, N_MovVert,
                        FeSps_Lev, SurfFeSps_Lev, N_array, N_List, SqMat, Mat, 
                        SqMat_low, Sol, Rhsarray, BoundCondition, 
                        GridBoundCondition, TMG, VeloVect, UPArrays, SArrays,
                        Slip_Joint, SqrStruct, Struct, SqrStruct_low, FE1D_List);

//        if(TDatabase::ParamDB->SC_VERBOSE > 0)
       OutPut( "Remeshing and Interpolation of velocity into new domain were done"<<endl);

       coll=FeSps_Lev[0][mg_level-1]->GetCollection();

       Getcellangle(FeSps_Lev[0][mg_level-1], Angle);
        OutPut( "MinAngle : "<< Angle[0]<< "  MaxAngle : "<<Angle[1]<< endl);

       N_Remesh ++;
       remeshed=TRUE;

       t2 = GetTime();

       if(TDatabase::ParamDB->SC_VERBOSE > 1)
        {
         cout<<endl;
         OutPut("Time for Remeshing: " << t2-t1 <<" Sec"<< endl);
        }
       Remesh_Time +=t2-t1;
   }

   // **************************************************************************************
   // end Remeshing
   // **************************************************************************************



    if(remeshed)
      {
       N_Active = FeSps_Lev[0][mg_level-1]->GetActiveBound();
       N_U = N_array[0][mg_level-1];
       N_P = N_array[1][mg_level-1];
       N_Unknowns = 2*N_U + N_P;
       N_SO = FeSps_Lev[5][mg_level-1]->GetN_DegreesOfFreedom();
       N_S = SurfFeSps_Lev[0][mg_level-1]->GetN_DegreesOfFreedom();
       N_SurfActive = SurfFeSps_Lev[0][mg_level -1]->GetActiveBound();

       N_G = N_array[2][mg_level-1];
       N_BoundaryNodes = N_G - FeSps_Lev[2][mg_level-1]->GetN_Inner();
       //  cout << "N_G " << N_G<<endl;

       delete [] refpos;
       delete [] auxpos;
       delete [] pos;
       delete [] tmp;
       delete [] d;

       refpos = new double[2*N_G];
       auxpos = new double[2*N_G];
       pos = new double[2*N_G];
       tmp = new double[2*N_G];
       d = new double[2*N_G];

       delete RefGridPos;
       delete AuxGridPos;
       delete GridPos;

       RefGridPos = new TFEVectFunct2D(FeSps_Lev[2][mg_level-1], refposString, gridString,
                                       refpos, N_G, 2);
       AuxGridPos = new TFEVectFunct2D(FeSps_Lev[2][mg_level-1], auxposString, gridString,
                                       auxpos, N_G, 2);
       GridPos = new TFEVectFunct2D(FeSps_Lev[2][mg_level-1], posString, gridString,
                                    pos, N_G, 2);

       RefGridPos->GridToData();
       GridPos->GridToData();
       AuxGridPos->GridToData();


       N_V = FeSps_Lev[3][mg_level-1]->GetN_DegreesOfFreedom();
       N_Vort = FeSps_Lev[4][mg_level-1]->GetN_DegreesOfFreedom();

     // cout << "N_U " << N_U<<endl;
     // cout << "N_Active " << N_Active<<endl;
    //  cout << "N_P " << N_P<<endl;

       delete [] frac_step_sol;
       delete [] defect;
       delete [] startsol;
       delete [] oldrhs;

       rhs = Rhsarray[0][mg_level-1];
       defect = new double[N_Unknowns];
       startsol = new double[N_Unknowns];
       frac_step_sol = new double[N_Unknowns];
       oldrhs =  new double[N_Unknowns];

       gamma = 0;
       // copy sol for extrapolation after time step
       memcpy(Sol[5][0], Sol[0][mg_level-1], N_Unknowns*SizeOfDouble);
       memcpy(frac_step_sol, Sol[0][mg_level-1],N_Unknowns*SizeOfDouble);


       delete [] SB;
       delete [] sdefect;
       delete [] ssol_old;

       SB = new double [N_S];
       memset(SB, 0, N_S*SizeOfDouble);
       sdefect = new double[N_S];
       memset(sdefect, 0, N_S*SizeOfDouble);
       ssol_old = new double[N_S];
       memcpy(ssol_old, ssol, N_S*SizeOfDouble);

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
         delete prec;
         prec = new TMultiGridIte(MatVect, Defect, NULL,
                                 0, N_Unknowns, TMG[0], zerostart);
        break;
      default:
        OutPut("Unknown preconditioner !!!" << endl);
       exit(4711);
    }
    switch (TDatabase::ParamDB->SC_SOLVER_SADDLE)
    {
      case 11:
        delete itmethod;
        itmethod = new TFixedPointIte(MatVect, Defect, prec,
                                      0, N_Unknowns, 0);
        if (TDatabase::ParamDB->SC_PRECONDITIONER_SADDLE == 5)
        {
         delete [] itmethod_sol;
         delete [] itmethod_rhs;
         itmethod_sol = new double[N_Unknowns];
         itmethod_rhs = new double[N_Unknowns];
          }
        else
        {
          itmethod_sol = Sol[0][mg_level-1];
          itmethod_rhs = rhs;
        }
        break;
      case 16:
        delete itmethod;
        itmethod = new TFgmresIte(MatVect, Defect, prec,
                                  0, N_Unknowns, 0);
        if (TDatabase::ParamDB->SC_PRECONDITIONER_SADDLE == 5)
        {
          delete [] itmethod_sol;
          delete [] itmethod_rhs;
          itmethod_sol = new double[N_Unknowns];
          itmethod_rhs = new double[N_Unknowns];
        }
        else
        {
          itmethod_sol = Sol[0][mg_level-1];
          itmethod_rhs = rhs;
        }
        break;
      default:
        OutPut("Unknown solver !!!" << endl);
        exit(4711);
      }
    }

   // prepare output(maxn_fespaces,  maxn_scalar,  maxn_vect, maxn_parameters, domain)
      delete Output;
      Output = new TOutput2D(3, 3, 2, 2, Domain);
      Output->AddFEVectFunct(VeloVect[0][mg_level-1]);
//       Output->AddFEVectFunct(VeloVect[1][mg_level-1]);
      Output->AddFEFunction(UPArrays[2][mg_level-1]);
      Output->AddFEFunction(UPArrays[5][mg_level-1]);
      os.seekp(std::ios::beg);
      Output->AddParameter(TDatabase::TimeDB->CURRENTTIME,os.str().c_str());

//   if((m % TDatabase::TimeDB->STEPS_PER_IMAGE != 0))
//    {
//      if(TDatabase::ParamDB->WRITE_VTK)
//       {
//        MapSurfToDomain(SArrays[0][mg_level-1], UPArrays[5][mg_level-1], 
//                        N_List[0][mg_level-1], N_List[1][mg_level-1]);
//        os.seekp(std::ios::beg);
//        if(img<10) os << VtkBaseName<<".0000"<<img<<".vtk" << ends;
//        else if(img<100) os << VtkBaseName<<".000"<<img<<".vtk" << ends;
//        else if(img<1000) os << VtkBaseName<<".00"<<img<<".vtk" << ends;
//        else if(img<10000) os << VtkBaseName<<".0"<<img<<".vtk" << ends;
//        else  os << VtkBaseName<<"."<<img<<".vtk" << ends;
//        Output->WriteVtk(os.str().c_str());
//        img++;
//       }
//     }
// exit(0);
  } // end manupulation due to remesh or reconstruction


   // ******************************************************************
   // Assembeling the grid matrix
   // ******************************************************************

   t1 = GetTime();
   // check in readin file

//    for(i=0;i<mg_level;i++)
    for(i=0;i<mg_level;i++)
    {
    fesp[0] = FeSps_Lev[2][i];
    SQMATRICES_GRID[0] = SqMat[8][i];
    SQMATRICES_GRID[0]->Reset();
    SQMATRICES_GRID[1] = SqMat[9][i];
    SQMATRICES_GRID[1]->Reset();
    SQMATRICES_GRID[2] = SqMat[10][i];
    SQMATRICES_GRID[2]->Reset();
    SQMATRICES_GRID[3] = SqMat[11][i];
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
     
    }

    t2 = GetTime();

    if((TDatabase::ParamDB->SC_VERBOSE > 1) )
     {
      cout << "Grid assembling done"<< endl;
      OutPut("Time for Grid assembling: " << t2-t1 << endl);
     }

     Entries[0] = SqMat[8][mg_level-1]->GetEntries();
     Entries[1] = SqMat[9][mg_level-1]->GetEntries();
     Entries[2] = SqMat[10][mg_level-1]->GetEntries();
     Entries[3] = SqMat[11][mg_level-1]->GetEntries();

     GridKCol = SqrStruct[1][mg_level-1]->GetKCol();
     GridRowPtr = SqrStruct[1][mg_level-1]->GetRowPtr();

     // for Dirichlet rows in off-diagonal matrices
     memset(Entries[1] + GridRowPtr[N_G-N_BoundaryNodes], 0,
            N_BoundaryNodes*SizeOfDouble);
     memset(Entries[2] + GridRowPtr[N_G-N_BoundaryNodes], 0,
            N_BoundaryNodes*SizeOfDouble);

 // ******************************************************************
 // End Assembeling the reference grid
 // ******************************************************************

 // ******************************************************************
 //  Reparametrization - Begin
 // *******************************************************************
// /*
   fhtot = 0.;
   fhmin = 100;
   fhmax = 0.0;
   x_max = -1e10;
   for(k=0;k<=N_MovVert[1];k++)
   {
    if(k==0)  MovBoundVert[0][N_MovVert[0]-1]->GetCoords(x1, y1);
    else MovBoundVert[1][k-1]->GetCoords(x1, y1);
    if(k==N_MovVert[1])  MovBoundVert[0][0]->GetCoords(x2, y2);
    else MovBoundVert[1][k]->GetCoords(x2, y2);

    if(x_max<x1) x_max = x1; 

    fh = sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
    fhtot +=fh;
    if (fh < fhmin) fhmin = fh;
    if (fh > fhmax) fhmax = fh;
   }

  fhtot /=N_MovVert[1];
  fhlimit = fhtot/2.0;


   if (((fhmin < fhlimit) || (fhmax > 3.5*fhtot/2.)) )
    {
      OutPut("FreeBound Edge Reparam:  "<<img<<' ' << fhlimit <<' '<< 3*fhtot/2.0
             <<' '<< fhmin <<' '<<fhmax<< endl);
//       if(TDatabase::ParamDB->WRITE_VTK)
//       {
//        os.seekp(std::ios::beg);
//        if(img<10) os << VtkBaseName<<".0000"<<img<<".vtk" << ends;
//        else if(img<100) os << VtkBaseName<<".000"<<img<<".vtk" << ends;
//        else if(img<1000) os << VtkBaseName<<".00"<<img<<".vtk" << ends;
//        else if(img<10000) os << VtkBaseName<<".0"<<img<<".vtk" << ends;
//        else  os << VtkBaseName<<"."<<img<<".vtk" << ends;
//        Output->WriteVtk(os.str().c_str());
//        img++;
//       }

    RefGridPos->GridToData();

    MapSurfToDomain(SArrays[0][mg_level-1], UPArrays[5][mg_level-1], 
                    N_List[0][mg_level-1], N_List[1][mg_level-1]);

//       GetSurfactMass(UPArrays[5][mg_level-1], SArrays[0][mg_level-1], 
//                      N_List[0][mg_level-1], N_List[1][mg_level-1], Surf_Mass);

    ReParam_axial3D_U(N_MovVert[1], Free_Cells,  N_array[3], N_array[4], 
                      VeloVect[0][mg_level-1], UPArrays[5][mg_level-1]);
    reparam = TRUE;

//       GetSurfactMass(UPArrays[5][mg_level-1], SArrays[0][mg_level-1], 
//                      N_List[0][mg_level-1], N_List[1][mg_level-1], Surf_Mass);
      MapDomainToSurf(UPArrays[5][mg_level-1], SArrays[0][mg_level-1], 
                      N_List[0][mg_level-1], N_List[1][mg_level-1]);
       memset(Sol[7][mg_level-1], 0, N_SO*SizeOfDouble);


//       if(TDatabase::ParamDB->WRITE_VTK)
//       {
//        MapSurfToDomain(SArrays[0][mg_level-1], UPArrays[5][mg_level-1], 
//                     N_List[0][mg_level-1], N_List[1][mg_level-1]);
// 
//        os.seekp(std::ios::beg);
//        if(img<10) os << VtkBaseName<<".0000"<<img<<".vtk" << ends;
//        else if(img<100) os << VtkBaseName<<".000"<<img<<".vtk" << ends;
//        else if(img<1000) os << VtkBaseName<<".00"<<img<<".vtk" << ends;
//        else if(img<10000) os << VtkBaseName<<".0"<<img<<".vtk" << ends;
//        else  os << VtkBaseName<<"."<<img<<".vtk" << ends;
//        Output->WriteVtk(os.str().c_str());
//        img++;
//       }
//       exit(0);
// */

   }
// */

// ********************************************************************************
//  Reparametrization of freesurface - end
//  Reparametrization of solid boundary begin
//   ******************************************************************************
    MovBoundVert[0][0]->GetCoords(x, y);
    MovBoundVert[0][N_MovVert[0]-1]->GetCoords(SLPX, SLPY);

    shtot = SLPY-y;
    shmax = 0.0;
    shmin = 100;
    for(k=0;k<N_MovVert[0]-1;k++)
     {
      MovBoundVert[0][k]->GetCoords(x1, y1);
      MovBoundVert[0][k+1]->GetCoords(x2, y2);
      if((y2-y1) < shmin ) shmin = y2-y1;
      if((y2-y1) > shmax ) shmax = y2-y1;
     }

    shtot /= (N_MovVert[0]-1);

//      cout<< SLPX<< ' ' << x + (N_MovVert[0]-1)*shtot  << endl;

    if(((fabs(shmin) < fabs(shtot/2.0)) || (fabs(shmax) > fabs(3.*shtot/2.0))) )
     {
         OutPut("Solid Edge Reparam: "<< img << ' ' << shtot/2.0 <<' '<< 3*shtot/2.0
                << ' ' << shmin << ' ' << shmax <<endl);

      if(!reparam) RefGridPos->GridToData();

      for(k=1;k<N_MovVert[0]-1;k++) // no need to set end vertices
       {
        x1 = 0.0;
        y1 = y + double(k)*shtot;;
//      cout<< " y1  reparam " << y1  << endl;
        MovBoundVert[0][k]->SetCoords(x1, y1);
       }
      reparam = TRUE;

    for(k=0;k<N_MovVert[0]-1;k++)
    {
     Slip_Joint[0][k]->UpdateParameters(MovBoundVert[0][k], MovBoundVert[0][k+1]);
    }
   }

 //  *****************************************************************************
 //  Reparametrization of solid boundary - end
 //  *****************************************************************************

   if(reparam)
    {
     GridPos->GridToData();
     RefGridPos->DataToGrid();

     memset(Rhsarray[2][mg_level-1], 0, 2*N_G*SizeOfDouble);
     memcpy(d, pos, 2*N_G*SizeOfDouble);
     Daxpy(2*N_G, -1, refpos, d);
     memcpy(Rhsarray[2][mg_level-1]+(N_G-N_BoundaryNodes), d+(N_G-N_BoundaryNodes),
           N_BoundaryNodes*SizeOfDouble);
     memcpy(Rhsarray[2][mg_level-1]+(2*N_G-N_BoundaryNodes), d+(2*N_G-N_BoundaryNodes),
           N_BoundaryNodes*SizeOfDouble);

     memset(Sol[4][mg_level-1], 0 , 2*N_G*SizeOfDouble);
     memcpy(Sol[4][mg_level-1]+(N_G-N_BoundaryNodes), d+(N_G-N_BoundaryNodes),
           N_BoundaryNodes*SizeOfDouble);
     memcpy(Sol[4][mg_level-1]+(2*N_G-N_BoundaryNodes), d+(2*N_G-N_BoundaryNodes),
           N_BoundaryNodes*SizeOfDouble);

     SolveGridEquation(Entries, Sol[4][mg_level-1], Rhsarray[2][mg_level-1],
                       GridKCol, GridRowPtr, N_G );

     memcpy(d, refpos, 2*N_G*SizeOfDouble);
     Daxpy(2*N_G, 1, Sol[4][mg_level-1], d);

     memcpy(pos, d , 2*N_G*SizeOfDouble);

     GridPos->DataToGrid();
     RefGridPos->GridToData();

//    if((m % TDatabase::TimeDB->STEPS_PER_IMAGE != 0))
//     {
//      if(TDatabase::ParamDB->WRITE_VTK)
//       {
// 
//         MapSurfToDomain(SArrays[0][mg_level-1], UPArrays[5][mg_level-1], 
//                     N_List[0][mg_level-1], N_List[1][mg_level-1]);
// 
//        os.seekp(std::ios::beg);
//        if(img<10) os << VtkBaseName<<".0000"<<img<<".vtk" << ends;
//        else if(img<100) os << VtkBaseName<<".000"<<img<<".vtk" << ends;
//        else if(img<1000) os << VtkBaseName<<".00"<<img<<".vtk" << ends;
//        else if(img<10000) os << VtkBaseName<<".0"<<img<<".vtk" << ends;
//        else  os << VtkBaseName<<"."<<img<<".vtk" << ends;
//        Output->WriteVtk(os.str().c_str());
//        img++;
//       }

// exit(0);
//     }
     reparam = FALSE;
  }
// *********************************************************************************************
//  Reparametrization - end
// *********************************************************************************************
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
    }

  min = 1e5;
  tx  = 0.;
//   m4   = 0;
//   if((m % TDatabase::TimeDB->STEPS_PER_IMAGE) == 0)


  if((m % TDatabase::TimeDB->STEPS_PER_IMAGE) == 0)
   {

     PrintSurfSurfactant(Free_Cells, MovBoundVert[1], N_array[3], 
                         N_MovVert[1], UPArrays[5][mg_level-1], N_BData);


     MovBoundVert[0][0]->GetCoords(x1, y1);


     OutPut(setw(25)<<"Time,  x_max, Top tip: "<< TDatabase::TimeDB->CURRENTTIME
                   <<"   "<< x_max <<"   "<< y1 << endl);

    Get_KE(VeloVect[0][mg_level-1], Params);

//     Get_Spheriity();

    OutPut(setw(25)<<"Time, Volume : " << TDatabase::TimeDB->CURRENTTIME
                   <<"   "<< Params[0]<<endl);
    OutPut(setw(25)<<"Time, Volume Diff : "<< TDatabase::TimeDB->CURRENTTIME
                   <<"   "<< Params[0] - InitVolume << endl)
    OutPut(TDatabase::TimeDB->CURRENTTIME<<
          " x_center_mass: "<< Params[2] <<
          " y_center_mass: "<< Params[3] << " U1_Rise: "<< Params[4] <<
          " U2_Rise: "<< Params[5]<< endl);
  OutPut(TDatabase::TimeDB->CURRENTTIME<<
         " Kinetic Energy: "<< Params[1] << endl);

      GetSurfactMass(UPArrays[5][mg_level-1], SArrays[0][mg_level-1], 
                     N_List[0][mg_level-1], N_List[1][mg_level-1], Surf_Mass);
  }


/*   if((m % 5) == 0)
   {
//     CurrVolume = Volume(FeSps_Lev[0][mg_level-1]);
    MovBoundVert[0][0]->GetCoords(Lx, Ly); // left wetting points
    MovBoundVert[0][N_MovVert[0]-1]->GetCoords(Rx, Ry); // right wetting points
    OutPut(setw(25)<<"Total No of time steps : " << m <<endl);
    OutPut(setw(25)<<"No of time Remeshed : " << N_Remesh <<endl);
    OutPut(setw(25)<<"No of time ReConstructed : "<<N_ReConstruct<< endl);
    OutPut(setw(25)<<"Time, Wett Len d : " << TDatabase::TimeDB->CURRENTTIME
                   <<"   "<< Ry-Ly<<endl);
    OutPut(setw(25)<<"Time, Volume : " << TDatabase::TimeDB->CURRENTTIME
                   <<"   "<< CurrVolume<<endl);
    OutPut(setw(25)<<"Time, Volume Diff : "<< TDatabase::TimeDB->CURRENTTIME
                   <<"   "<< CurrVolume - InitVolume << endl);


   MovBoundVert[0][N_MovVert[0]-2]->GetCoords(x1, y1);
   MovBoundVert[1][0]->GetCoords(x2, y2);
   tx = x1-Rx;
   sx = x2-Rx;
   ty = y1-Ry;
   sy = y2-Ry;
   R_Theta[0] = acos( (tx*sx+ty*sy)/(sqrt(tx*tx+ty*ty)* sqrt(sx*sx+sy*sy)) )*(180/3.141592654);

   UPArrays[1][mg_level-1]->FindGradient(x1, y1, RU);
   OutPut(setw(25)<<"T, 1_pt_bef R_Wet_pt : " << TDatabase::TimeDB->CURRENTTIME
                   <<"   "<< RU[0]<<endl);


   MovBoundVert[0][1]->GetCoords(x1, y1);
   MovBoundVert[1][N_MovVert[1]-1]->GetCoords(x2, y2);
   tx = x1-Lx;
   sx = x2-Lx;
   ty = y1-Ly;
   sy = y2-Ly;
   L_Theta[0] = acos( (tx*sx+ty*sy)/(sqrt(tx*tx+ty*ty)* sqrt(sx*sx+sy*sy)) )*(180/3.141592654);

   UPArrays[1][mg_level-1]->FindGradient(x1, y1, LU);
   OutPut(setw(25)<<"T, 1_pt_nex L_Wet_pt : " << TDatabase::TimeDB->CURRENTTIME
                  <<"   "<< LU[0]<<endl);

    UPArrays[1][mg_level-1]->FindGradient(Rx, Ry, RU);
    UPArrays[1][mg_level-1]->FindGradient(Lx, Ly, LU);
    OutPut(setw(25)<<"T, Right Wet_pt_velo : " << TDatabase::TimeDB->CURRENTTIME
                   <<"   "<< RU[0]<<endl);
    OutPut(setw(25)<<"T, Left Wet_pt_velo : " << TDatabase::TimeDB->CURRENTTIME
                   <<"   "<< LU[0]<<endl);


   MovBoundVert[0][N_MovVert[0]-2]->GetCoords(x1, y1);
   MovBoundVert[1][0]->GetCoords(x2, y2);
   MovBoundVert[1][1]->GetCoords(x3, y3);
   x2 = (x2+x3)/2.0;
   y2 = (y2+y3)/2.0;
   tx = x1-Rx;
   sx = x2-Rx;
   ty = y1-Ry;
   sy = y2-Ry;
   R_Theta[1] = acos( (tx*sx+ty*sy)/(sqrt(tx*tx+ty*ty)* sqrt(sx*sx+sy*sy)) )*(180/3.141592654);

   MovBoundVert[0][1]->GetCoords(x1, y1);
   MovBoundVert[1][N_MovVert[1]-1]->GetCoords(x2, y2);
   MovBoundVert[1][N_MovVert[1]-2]->GetCoords(x3, y3);
   x2 = (x2+x3)/2.0;
   y2 = (y2+y3)/2.0;
   tx = x1-Lx;
   sx = x2-Lx;
   ty = y1-Ly;
   sy = y2-Ly;
   L_Theta[1] = acos( (tx*sx+ty*sy)/(sqrt(tx*tx+ty*ty)* sqrt(sx*sx+sy*sy)) )*(180/3.141592654);

   MovBoundVert[0][N_MovVert[0]-2]->GetCoords(x1, y1);
   MovBoundVert[1][0]->GetCoords(x2, y2);
   MovBoundVert[1][1]->GetCoords(x3, y3);
   MovBoundVert[1][2]->GetCoords(x4, y4);
   x2 = (x2+x3+x4)/3.0;
   y2 = (y2+y3+y4)/3.0;
   tx = x1-Rx;
   sx = x2-Rx;
   ty = y1-Ry;
   sy = y2-Ry;
   R_Theta[2] = acos( (tx*sx+ty*sy)/(sqrt(tx*tx+ty*ty)* sqrt(sx*sx+sy*sy)) )*(180/3.141592654);

   MovBoundVert[0][1]->GetCoords(x1, y1);
   MovBoundVert[1][N_MovVert[1]-1]->GetCoords(x2, y2);
   MovBoundVert[1][N_MovVert[1]-2]->GetCoords(x3, y3);
   MovBoundVert[1][N_MovVert[1]-3]->GetCoords(x4, y4);
   x2 = (x2+x3+x4)/3.0;
   y2 = (y2+y3+y4)/3.0;
   tx = x1-Lx;
   sx = x2-Lx;
   ty = y1-Ly;
   sy = y2-Ly;
   L_Theta[2] = acos( (tx*sx+ty*sy)/(sqrt(tx*tx+ty*ty)* sqrt(sx*sx+sy*sy)) )*(180/3.141592654);


    OutPut(setw(25)<<"T, Right Angle_1 : " << TDatabase::TimeDB->CURRENTTIME
                   <<"   "<< R_Theta[0]<<endl);
    OutPut(setw(25)<<"T, Left Angle_1 : " << TDatabase::TimeDB->CURRENTTIME
                   <<"   "<< L_Theta[0]<<endl);

    OutPut(setw(25)<<"T, Right Angle_2 : " << TDatabase::TimeDB->CURRENTTIME
                   <<"   "<< R_Theta[1]<<endl);
    OutPut(setw(25)<<"T, Left Angle_2 : " << TDatabase::TimeDB->CURRENTTIME
                   <<"   "<< L_Theta[1]<<endl);

    OutPut(setw(25)<<"T, Right Angle_3 : " << TDatabase::TimeDB->CURRENTTIME
                   <<"   "<< R_Theta[2]<<endl);
    OutPut(setw(25)<<"T, Left Angle_3 : " << TDatabase::TimeDB->CURRENTTIME
                   <<"   "<< L_Theta[2]<<endl);


//   if(l==0  && (m % 5) == 0)
   {
    Get_KE(VeloVect[0][mg_level-1], Params);

//  parameters[0] = volume;
//  parameters[1] = KE;
//  parameters[2] = x_mass;
//  parameters[3] = y_mass;
//  parameters[4] = U1_Rise;
//  parameters[5] = U2_Rise;

  KE = Params[1];

  OutPut(TDatabase::TimeDB->CURRENTTIME<<
         " x_center_mass: "<< Params[2] <<
         " y_center_mass: "<< Params[3] << " U1_Rise: "<< Params[4] <<
         " U2_Rise: "<< Params[5]<< endl);

  FreeSurf_Methods(N_MovVert[1], Free_Cells,  N_array[3], Params);
  
//      Params[0] = t;
//   Params[1] = x_max;
//   Params[2] = y_max;
//   Params[3] = x[m-1]-x[0];       // wetting diameter 
//   Params[4] = y[m-1]-y[0];        // to check that 0th and (m-1)th points are wetting points
  ST=TDatabase::ParamDB->CHAR_L0;
  MovBoundVert[0][0]->GetCoords(Lx, Ly); // left wetting points
  MovBoundVert[0][N_MovVert[0]-1]->GetCoords(Rx, Ry); // right wetting points
  tx=Rx-Lx;
  sx = (L_Theta[0]+R_Theta[0])/2.0;

  Surface_Energy = tx-(Params[0]*cos((Pi/180.)*sx));
//   Surface_Energy = ST*(tx-(Params[0]*cos(sx))); 

  OutPut(TDatabase::TimeDB->CURRENTTIME<<
         " Kinetic Energy: "<< KE <<" Surface Energy: "<< Surface_Energy <<
         " FreeSurf Length: "<< Params[0] <<" x_max: "<< Params[1] << " y_max: "<< Params[2] << endl);

  }

   }

 */

//   if((m % 100) == 0)
//    {
//       os.seekp(std::ios::beg);
// //       os << "Boundary"<< i << ".dat" << ends;
//        if(N_BData<10) os << "Boundary.0000"<<N_BData<<".dat" << ends;
//          else if(N_BData<100) os << "Boundary.000"<<N_BData<<".dat" << ends;
//          else if(N_BData<1000) os << "Boundary.00"<<N_BData<<".dat" << ends;
//          else if(N_BData<10000) os << "Boundary.0"<<N_BData<<".dat" << ends;
//          else  os << "Boundary."<<N_BData<<".dat" << ends;
// 
//       std::ofstream dat(os.str().c_str());
//       if (!dat)
//        {
//         cerr << "cannot open file for output" << endl;
//         return -1;
//        }
//       dat << "# Boundary data created for droplet by MooNMD" << endl;
//       dat << "# Current Time :" << TDatabase::TimeDB->CURRENTTIME << endl;
//       for(k=0;k<N_MovVert[0];k++) // no need to set end vertices again
//        {
//         MovBoundVert[0][k]->GetCoords(x1, y1);
// 
//         dat << x1 << " " <<  y1<< endl;
//        }
//       for(k=0;k<N_MovVert[1];k++) // no need to set end vertices again
//        {
//         MovBoundVert[1][k]->GetCoords(x1, y1);
//         dat << x1 << " " <<  y1<< endl;
//        }
//       dat.close();
//       cout << endl;
//       cout << "Boundary wrote output into file " << endl;
//       N_BData++;
// 
//    }

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

    CurrVolume = Volume(FeSps_Lev[0][mg_level-1]);
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
  OutPut("total running time: " << total_time<<" Sec" << endl);
  CloseFiles();
  return 0;
}

