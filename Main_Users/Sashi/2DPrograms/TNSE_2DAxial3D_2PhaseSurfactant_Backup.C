// =======================================================================
//
// Purpose:     main program for 2phase flows with soluble surfactant
//
// Author:     Sashikumaar Ganesan
// History:
//  Strat of implementation:    01.01.2009
//  update                 :    18.04.2010
//  Latest update          :    31.05.2011
// =======================================================================
#include <mcheck.h>

#include <Domain.h>
#include <Database.h>
#include <FEDatabase2D.h>
#include <FESpace2D.h>
#include <SquareStructure1D.h>
#include <SquareStructure2D.h>
#include <SquareMatrix1D.h>
#include <FEFunction1D.h>
#include <Structure2D.h>
#include <AuxParam2D.h>
#include <QuadAffin.h>
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
// #include <Remesh2D.h>
#include <IsoInterfaceJoint.h>
#include <BdLine.h>
#include <BdCircle.h>
#include <GridCell.h>

#include <QuadAffin.h>
#include <QuadBilinear.h>
#include <QuadIsoparametric.h>
#include <TriaAffin.h>
#include <TriaIsoparametric.h>

#include <ParDirectSolver.h>

#include <sys/stat.h>
#include <sys/types.h>

#define AMG 0
#define GMG 1

#include "../Examples_All/TNSE_2D/2Phase_SurfDrop.h"
// #include "../Examples/TNSE_2D/2Phase_Cleavage.h"
// #include "../Examples/TNSE_2D/MassTransferTest.h"

extern "C"
{
  void triangulate(char*, struct triangulateio*,
                   struct triangulateio*, struct triangulateio*);
}


void Get_KE(TFEVectFunct2D *Velocity, double *parameters)
 {
  int i,j,k,l, polydegree, Phase_No, ORDER;
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
   Phase_No = cell->GetPhase_ID();
   if(Phase_No==0)
   {
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
   } // if(Phase_No==0)
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
//    cout<<i<< " Angle of free Vertices "<<(180./Pi)*atan2(y,x)<<endl;

    ArcLength += sqrt( (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) );

    Surfact->FindGradient( x1,  y1, T_val);
    dat << x1<< " " << y1 << " "<< ArcLength << "  " << T_val[0] <<endl;

    x1=x2; y1=y2;
   }

      dat.close();
      cout << endl;
      cout << "Surfactant data wrote into file " << endl;
 N_BData++;

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

}


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



// sorting curved surfrace vertices - general case
void SortIsoVert_Gen(TBaseCell **cell, TVertex **Vertex, int *CellNo, int *EdgeNo, int N,
                     double X0, double Y0, int indicator)
{
 int i, j, k, temp, test=0;
 double x, y, x1, y1;
 TVertex *temp_vert;
 TBaseCell *temp_cell;

 // finding the right(starting) vertex
  for(i=0;i<N;i++)
   {
    Vertex[i]->GetCoords(x, y);
// sort with X0 as the first vertex ordinate if indicator==0 (or) ....
    if( (indicator==0 && fabs(x-X0)<1e-10 ) || (indicator==1 && fabs(y-Y0)<1e-10 ) )
      {
//        cout << " sorting " << x << ' ' << y<<endl;
       temp_vert = Vertex[0];
       Vertex[0] = Vertex[i];
       Vertex[i] = temp_vert;

       temp_cell = cell[0];
       cell[0] = cell[i];
       cell[i] = temp_cell;
       
       temp = CellNo[0];
       CellNo[0] = CellNo[i];
       CellNo[i] = temp;       
       
       temp = EdgeNo[0];
       EdgeNo[0] = EdgeNo[i];
       EdgeNo[i] = temp;
       test++;
      }
    if(test) break;
   }

//   test = 0;
   //
  for(i=0;i<N-1;i++)
   {
    test = 0; 
    k = cell[i]->GetN_Edges();
    cell[i]->GetVertex((EdgeNo[i]+1) % k)->GetCoords(x, y);
     for(j=i;j<N;j++)
      {
       Vertex[j]->GetCoords(x1, y1);
       if((x==x1) && (y==y1))
        {
         temp_vert = Vertex[j];
         Vertex[j] = Vertex[i+1];
         Vertex[i+1] = temp_vert;
       
         temp_cell = cell[j];
         cell[j] = cell[i+1];
         cell[i+1] = temp_cell;

         temp = CellNo[j];
         CellNo[j] = CellNo[i+1];
         CellNo[i+1] = temp; 
	 
         temp = EdgeNo[j];
         EdgeNo[j] = EdgeNo[i+1];
         EdgeNo[i+1] = temp; 
         test++;
        }
      if(test) break;  
     }
   }

   //print   

//   for(i=0;i<N;i++)
//    {
//     Vertex[i]->GetCoords(x, y);
//     cout<<i<< " Angle of free Vertices "<<(180/Pi)*atan2(y-0.5,x)<<endl;
//    }
//
//   Vertex[0]->GetCoords(x, y);
//   cout << x << ' ' << y<< ' '<<endl;
//   Vertex[N-1]->GetCoords(x, y);
//   cout << x << ' ' << y<< ' '<<endl;
//   Vertex[N-2]->GetCoords(x, y);
//   cout << x << ' ' << y<< ' '<<endl;

// exit(0);

 }




void  Domain2DSurf_2Phase(TCollection *Coll, TDomain *SurfDomain, int **N_List)
  {
  int i, j, k, l, m, n,  m1, N_SurfCells, N_Cells, N_SurfVert, ID, N;
  int maxEpV, a, b, len1, len2, Neighb_tmp, PhaseID;
  int *Bd_Part, *Lines, N_G, *PointNeighb, CurrNeib,  Neib[2];
  int *Cell_No, *Joint_No;

  const int *TmpEV;
  TBaseCell *Me;
  TVertex **SurfVetrex;
  TJoint *Joint;
  TBaseCell  **SurfCellTree, *NullCell = NULL;
  boolean StartBD=FALSE, EndBD=FALSE;


  Cell_No = N_List[0];
  Joint_No = N_List[1];

  N_SurfCells = 0;
  N_Cells = Coll->GetN_Cells();
  for(i=0;i<N_Cells;i++)
  {
   Me = Coll->GetCell(i);
   PhaseID = Me->GetPhase_ID();

   if(PhaseID == 1) // outer phase cell is better for soluble surfactant coupling
   {
    k = Me->GetN_Edges();
    for(l=0;l<k;l++)
     {
       Joint = Me->GetJoint(l);
       if(Joint->GetType() == IsoBoundEdge   || Joint->GetType() == IsoInterfaceJoint
           || Joint->GetType() == InterfaceJoint  )
         N_SurfCells++;
     } // endfor l
   } //  if(PhaseID == 1)
  }// endfor i

  OutPut("N_SurfCells: " << N_SurfCells << endl);

  Bd_Part = new int[N_SurfCells];
  Lines = new int[2*N_SurfCells];
  SurfVetrex = new TVertex*[2*N_SurfCells]; // maximum possible vertex

//   m = 0;
  m1 = 0;
  N    = 0; 
  N_SurfVert = 0;
  for(i=0;i<N_Cells;i++)
  {
   Me = Coll->GetCell(i);
   PhaseID = Me->GetPhase_ID();

   if(PhaseID == 1) // outer phase cell is better for soluble surfactant coupling
   {
    k = Me->GetN_Edges();
    for(l=0;l<k;l++)
     {
       Joint = Me->GetJoint(l);
       if(Joint->GetType() == IsoBoundEdge   || Joint->GetType() == IsoInterfaceJoint
           || Joint->GetType() == InterfaceJoint  )
        {
          Cell_No[N] = i;
          Joint_No[N] = l;
          Bd_Part[N++] =(((TBoundEdge *)Joint)->GetBoundComp())->GetID();
          Me->GetShapeDesc()->GetEdgeVertex(TmpEV);

           for(n=0;n<2;n++) // In 2D each edge contains only two vertex
            {
              ID = 0;
              for(m=0;m<N_SurfVert;m++)
              if( SurfVetrex[m]==Me->GetVertex(TmpEV[2*l + n ]) )
               {
//             vertex is already in the collection
                ID = 1;
                Lines[m1] = m;
                m1++;
                break; 
              }

              if(ID==0)
               {
                 Lines[m1] = N_SurfVert;
                 m1++;
                 SurfVetrex[N_SurfVert] = Me->GetVertex(TmpEV[2*l + n ]);
                 N_SurfVert++;
               }

           } //  for(n=0
        } // if(Joint->GetType()
     } // endfor l
   } // if(PhaseID == 1)
  }// endfor i

   OutPut("N_SurfVert: " << N_SurfVert << endl);

   SurfCellTree = new TBaseCell*[N_SurfCells];

   for (i=0;i<N_SurfCells;i++)
   {
    SurfCellTree[i] = new TMacroCell(TDatabase::RefDescDB[S_Line], 0);
    SurfCellTree[i]->SetVertex(0, SurfVetrex[ Lines[ 2*i       ]]);
    SurfCellTree[i]->SetVertex(1, SurfVetrex[ Lines[ 2*i + 1]]);

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
   maxEpV = 2; // max possible degree of each vertex
   N_G = N_SurfVert;
   maxEpV++;
   PointNeighb = new int[maxEpV * N_G];
   memset(PointNeighb, 0, maxEpV * N_G *SizeOfInt);

    // every row contains "maxEpV" columns
    // for every row at first colomn contains the number of cells containing this vertex
    // at further columns we set the index of the corresponding cells containing this vertex
   for(i=0;i<2*N_SurfCells;i++)
   {
     j = Lines[i]*maxEpV;
     PointNeighb[j]++;
     PointNeighb[j + PointNeighb[j]] = i / 2;
   }

  N_G = N_SurfCells;

  for(i=0;i<N_G;i++)
   {
    StartBD=FALSE;
    EndBD=FALSE;
    a =  Lines[2*i ];
    b =  Lines[2*i + 1];
    Neib[0] = -1;
    Neib[1] = -1;
    CurrNeib = 0;
    len1 = PointNeighb[a*maxEpV];
    len2 = PointNeighb[b*maxEpV];


    if(len2==1)
     {
      Neib[0] = PointNeighb[b*maxEpV + 1];
      Joint = new TJointEqN(SurfCellTree[Neib[0]]);
      SurfCellTree[Neib[0]]->SetJoint(1, Joint);
      cout<<" End cell " << Neib[0] <<" Neib[1] " << Neib[1] <<endl;
     }

  // find indexes of cells containing the current vertex (joint)
    if(len1==1)
      Neib[0] = PointNeighb[a*maxEpV + 1];
    else
      for(j=1;j<=len1;j++)
       {
        Neighb_tmp = PointNeighb[a*maxEpV + j];
        for (k=1;k<=len2;k++)
         {
          if(Neighb_tmp == PointNeighb[b*maxEpV + k])
           {
            Neib[0] = Neighb_tmp;

            if(j==1)
             Neib[1] = PointNeighb[a*maxEpV + j + 1];
            else if(j==2)
             Neib[1] = PointNeighb[a*maxEpV + j - 1];

            break;
           }
         }
       } //  for (j=1;j<=len1


    if(len1==2)
     {
      Joint = new TJointEqN(SurfCellTree[Neib[0]], SurfCellTree[Neib[1]]);

      SurfCellTree[Neib[0]]->SetJoint(0, Joint);
      SurfCellTree[Neib[1]]->SetJoint(1, Joint);
     }
    else if(len1==1)
     {
      Joint = new TJointEqN(SurfCellTree[Neib[0]]);
      SurfCellTree[Neib[0]]->SetJoint(0, Joint);
     cout<< " Start cell " << Neib[0] <<" Neib[1] " << Neib[1] <<endl;
     }

   } //  for(i=0;i<N_G;i++)

  delete [] SurfVetrex;
  delete [] PointNeighb;
  delete []  Bd_Part;
  delete []  Lines;
 }

void ReParam_axial3D_Data(int &N_E, TBaseCell **cell, int *EdgeNo,  int *CellNo, 
                          TFEVectFunct2D *Velocity, TFEFunction2D *Surfactant,
                          double *&Intpol_Coord, double *&Intpol_VeloValues, double *&Intpol_Values,
                          double h_min, double **&FreePts)
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
        // cout<< i<<" FreeGaus " << (180/Pi)*atan2(y[m], x[m]) <<endl;
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

    if(i != N_E-1)// -1 due to end dof will be the start dof of the next edge except on last edge
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
     cout<< " m " << m << " m1 " << m1 <<" N_Splines " << N_Splines <<endl;
     exit(0);
    }
//           cout<< " m " << m << " m1 " << m1 <<" N_Splines " << N_Splines <<endl;
// exit(0);

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

// =====================================================================================
// spline construction done
// =====================================================================================
   delete [] FreePts[0];
   delete [] FreePts[1];
   
   N_E = int(t[N_Splines]/h_min);
 
   FreePts[0] = new double[N_E+1];
   FreePts[1] = new double[N_E+1];

   Intpol_Coord = new double[4*(N_E+1)];  
   Intpol_Values = new double[2*(N_E+1)];   
   Intpol_VeloValues = new double[4*(N_E+1)];
   
   FreePts[0][0] = x[0];
   FreePts[1][0] = y[0];   
   Intpol_Coord[0] = x[0];
   Intpol_Coord[1] = y[0];

   Intpol_VeloValues[0] = u1_spl[0];
   Intpol_VeloValues[1] = u2_spl[0];    
   Intpol_Values[0] = surf_spl[0];
   
   teta = 1.0/(double)(2.*N_E); 
   
   //cout<< "X " <<FreePts[0][0] <<" Y " <<FreePts[1][0] <<endl; 
   
//    m=0;
//    for(j=0;j<N_E;j++)
//     {
//      T = (double)m*teta;   
//      cout<< "T " <<T << endl;  
//      m++;
//      
//      T = (double)m*teta;   
//      cout<< "T " <<T << endl;   
//      m++;
//     }
    
    
//    cout<< "Spline " <<endl;
//    exit(0);
   
   T = 0;
   Param9[0] = 0;
   for(i=1;i<=N_Splines;i++) 
    Param9[i] = Params[(i-1)*10+8];

   m = 0;
   m1 = 0;
 
   for(j=0;j<N_E;j++)
    {
     T = (double)m*teta;
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
      } // for(i=1;i<=N_Splines;i++)

   phi1 = (2.*T*T - 3.*T)*T + 1.;
   phi2 = (-2.*T + 3.)*T*T;
   phi3 = (T*T - 2.*T + 1.)*T;
   phi4 = (T - 1)*T*T;

   X = Params[ISpline    ]*phi1 + Params[ISpline + 2]*phi2 +
       Params[ISpline + 4]*phi3 + Params[ISpline + 6]*phi4;
   Y = Params[ISpline + 1]*phi1 + Params[ISpline + 3]*phi2 +
       Params[ISpline + 5]*phi3 + Params[ISpline + 7]*phi4;

   if(fabs(X)<1.e-12)
    X = 0.;       

   if(m!=0)
    {
     FreePts[0][j] = X;
     FreePts[1][j] = Y;
     
     Intpol_Coord[2*m] = X;
     Intpol_Coord[2*m+1] = Y; 
    }
   
    m++;
// ================================================================================
// for fe values

     u0 = u1_spl[FeDof]*phi1 + u1_spl[FeDof+1]*phi2 +
              FEParams[USpline]*phi3 + FEParams[USpline + 1]*phi4;
     u1 = u2_spl[FeDof]*phi1 + u2_spl[FeDof+1]*phi2 +
              FEParams[USpline+2]*phi3 + FEParams[USpline + 3]*phi4;

     surf = surf_spl[FeDof]*phi1 + surf_spl[FeDof+1]*phi2 +
              FEParams[USpline+4]*phi3 + FEParams[USpline + 5]*phi4;

    if(j!=0) // endpoints no need to set
     {      
      Intpol_VeloValues[2*m1] = u0;
      Intpol_VeloValues[2*m1+1] = u1;
      Intpol_Values[m1] = surf;
     }
    m1++;
// ================================================================================
//     Joint = cell[j]->GetJoint(EdgeNo[j]);
//     isojoint = (TIsoBoundEdge *)Joint;
//     k = isojoint->GetN_Vertices();
       k = ORDER-1;
//     if(k==ORDER-1)
//      {
//       IsoVertices = isojoint->GetVertices();
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
           // cout<< ISpline << ' ' << T;
           T = (T-Param9[i-1])/(Param9[i]-Param9[i-1]);
           // cout<< ' ' << T <<endl;
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

     Intpol_Coord[2*m] = X;
     Intpol_Coord[2*m+1] = Y; 

     // OutPut("NewX:"<<' '<< m <<' '<<X<<' '<< Y<<endl);
     m++;

     // ====================================================================================
     // for fe values

     u0 = u1_spl[FeDof]*phi1 + u1_spl[FeDof+1]*phi2 +
              FEParams[USpline]*phi3 + FEParams[USpline + 1]*phi4;
     u1 = u2_spl[FeDof]*phi1 + u2_spl[FeDof+1]*phi2 +
              FEParams[USpline+2]*phi3 + FEParams[USpline + 3]*phi4;

     surf = surf_spl[FeDof]*phi1 + surf_spl[FeDof+1]*phi2 +
               FEParams[USpline+4]*phi3 + FEParams[USpline + 5]*phi4;
    
     Intpol_VeloValues[2*m1] = u0;
     Intpol_VeloValues[2*m1+1] = u1;   
     Intpol_Values[m1] = surf;      
     m1++;
     // ====================================================================================
     }
//     }
   }  //  for(j=0;j<N_E
   
//    T = ((double)m)*teta;   
//    cout<< "T " <<T << endl;    

   FreePts[0][N_E] = x[N_Splines];
   FreePts[1][N_E] = y[N_Splines];   
   
  Intpol_Coord[2*m] = x[N_Splines];
  Intpol_Coord[2*m+1] = y[N_Splines];
  
  Intpol_VeloValues[2*m1] = u1_spl[N_Splines];;
  Intpol_VeloValues[2*m1+1] = u2_spl[N_Splines];  
  Intpol_Values[m1] = surf_spl[N_Splines];     

//    cout<< "EndX " <<FreePts[0][N_E] <<" EndY " <<FreePts[1][N_E] <<endl;  
//    cout<< m1 << "VeloX " <<Intpol_VeloValues[2*m1] <<" VeloY " <<Intpol_VeloValues[2*m1+1] <<endl;     

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
   delete [] FEParams;   
   
//       cout<< "Spline " << Intpol_Values[m1] <<endl;
//    exit(0);

} // ReParam_axial3D_Data


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
   delete [] FEParams;
//    exit(0);
}


void SurfValuesForInterpol(int N_E, int *Cell_Nos, int *Edge_Nos, TFEFunction2D *fefunction,
                           double *Intpol_Coord, double *Intpol_Values)
{
  int i, i3, j, k, l, m, N, *GlobalNumbers, *BeginIndex, *DOF;
  int *JointDOF, N_JointDOF;
  
  double *Values;
  
  TFESpace2D *FeSpace; 
  TCollection *Coll;
  TBaseCell *Me;
  FE2D FeId;
  TFEDesc2D *FeDesc;
  TIsoInterfaceJoint *isoIntjoint;
  TJoint *Joint;
  TVertex **IsoVertices;
  
  FeSpace = fefunction->GetFESpace2D();
  GlobalNumbers = FeSpace->GetGlobalNumbers();
  BeginIndex = FeSpace->GetBeginIndex();
  
  Coll = FeSpace->GetCollection();
  Values = fefunction->GetValues(); 
  
  m=0;
  for(i=0; i<N_E; i++)
   {
    N = Cell_Nos[i];
    Me = Coll->GetCell(N);
    FeId = FeSpace->GetFE2D(N, Me);
    FeDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(FeId);
    
    l= Edge_Nos[i];
    JointDOF = FeDesc->GetJointDOF(l);
    N_JointDOF = FeDesc->GetN_JointDOF();
    DOF = GlobalNumbers + BeginIndex[N];

    //starting vertex
    (Me->GetVertex(l))->GetCoords(Intpol_Coord[2*m], Intpol_Coord[2*m+1]);
    Intpol_Values[m] = Values[DOF[JointDOF[0]]];
    m++; 

    Joint = Me->GetJoint(l);
    isoIntjoint = (TIsoInterfaceJoint *)Joint;
    k = isoIntjoint->GetN_Vertices();
    IsoVertices = isoIntjoint->GetVertices();
    
    if(k!=N_JointDOF-2)
     {
      OutPut("No match in isopoints per interface in SurfValuesForInterpol" <<endl;) 
      exit(0);
     }
     
     for(i3=0;i3<k;i3++)
      {    
       IsoVertices[i3]->GetCoords(Intpol_Coord[2*m], Intpol_Coord[2*m+1]);
       Intpol_Values[m] = Values[DOF[JointDOF[i3+1]]];
       m++;               
      }//  for(i3=0;i3<k;i3++
     
   } //  for(i=0; i<N_E; i++)
    
  //end vertex
  (Me->GetVertex((l+1)%3))->GetCoords(Intpol_Coord[2*m], Intpol_Coord[2*m+1]);
  Intpol_Values[m] = Values[DOF[JointDOF[N_JointDOF-1]]];
  m++;    
   
//    for(i=0;i<m;i++)
//     cout << i << " X[m] " << Intpol_Coord[2*i] <<" Y[m] " << Intpol_Coord[2*i+1] << " Intpol_Values[m] " << Intpol_Values[i] <<  endl;
//   exit(0);
}

void Remesh2D_2PhaseAxial3D(TDomain *Domain, TBaseCell ***Coll_Cells,
              TCollection  **Coll_Multi, TVertex ***MovBoundVert, int *N_MovVert,
              int N_FESpaces_All, TFESpace2D ***FeSpaces, int **N_DOFs, int **GlobalCell_Index, 
              TSquareMatrix2D ****SqMat, 
              TMatrix2D ***Mat, double ***Sol, double ***Rhs, TFEVectFunct2D ***VeloVect,
              TFEFunction2D ***FeFunct, TBoundEdge  ***Slip_Joint, TSquareStructure2D ***SqrStruct,
              TStructure2D ***Struct, double **FreePts,
              TDomain *SurfDomain, TFESpace1D **IFaceFeSpaces, int **N_List,
              TFEFunction1D **IFaceFeFunct, TSquareMatrix1D **SqMat_IFace,
              TSquareStructure1D **IFaceStruct, FE1D *FE1D_List,
              int ***GlobalNumbers, int ***BeginIndex, int **Bound_DOFs, double *Intpol_Coord, double *Intpol_Values, double *Intpol_VeloValues)
{
  int i, j, k, l, N, N_RootCells, Old_N_RootCells, Old_S_N_RootCells, N_Cells, Old_N_Cells, N_Joints;
  int N_Vertices, i3, N_G, N_Old_Face_Vert, N_G_P1, N_G_P2, N_BoundaryNodes_P1, N_BoundaryNodes_P2;
  int m, m1, m2, m3, m4, m5, m6, m7, m8, N_SlipJoints, N_FreeJoints, N_Velo_IsoPoints;
  int ORDER, VSP, N_OldU, N_OldP, *Edge_No, l1, jj, kk, length, *DOF_P2, *DOF;
  int order, N_Active, N_P, N_U, N_V, N_Vort;
  int N_BoundaryNodes, N_Unknowns, n_aux, mixing_layer=0;
  int velocity_space_code, pressure_space_code;
  int *PointNeighb, maxEpV = 0, a, b, Neighb_tmp, Neib[2];
  int CurrNeib, len1, len2, CurrComp, comp, N_Interf_Vertices, N_Cells_P2,  N_Cells_P1;
  int In_Index, CurrVertex, CurrJoint, CurrInterfCell, ID, bct, bct1;
  int *Triangles, *Innercell_No, N_Parameters=1;
  int ChangeBound1 = 0, ChangeBound2 = 0, temp_segment;
  int N_Hori, N_Verti, N_vert1, N_vert2, N_vert3, beg_index, end_index,bdpart;
  int sp_No, sp_no0, sp_no1, N_IFaceCells, *mm = new int[8];
  int Old_N_IsoEdges, *GlobalNumbers_S, *BeginIndex_S, *JointDOF, N_JointDOF, Cell_No;


  double x_short, y_short,  temp0, temp1, *usol, *psol, *csol, *ssol;
  double left, right, top, bottom, rad1, alpha, Lx =-100, Ly =-100, Rx=-100, Ry=-100;
  double T_a, T_b, temp, temp3, T, *Coordinates, x, y=0, tx, ty, TX[2], TY[2], x1, x2;
  double x0, y0, hi, area, x_mid, y_mid, Parameters[2];
  double d, t, hE, Pos_Indica[3], *I_FaceX,  *I_FaceY, Pos_Indicator;
  double *New_Sol, *oldssol, C_x, C_y, nx, ny, hmin, *Surfact_sol, *ValuesU2;

  char *PsBaseName, *GrapeBaseName, *GnuBaseName, *ReadGrapeBaseName;
  char *VtkBaseName;

  boolean inner;


  TIsoInterfaceJoint *Free_J;
  TBoundPart *BoundPart;
  TBoundComp *BoundComp;
  TBdLine *UpdateSlipBound, *UpdateBound[6];
  TBaseCell **CellTree, **Old_CellTree, *cell, *Me, **del_Cell_Array;
  TBaseCell **Cells_P1, **Cells_P2, **Free_Cells, **Surf_OldCellTree;
  TCollection *Coll, *Old_Coll, *Space_Old_Coll, *mortarcoll = NULL,  *Coll_P1, *Coll_P2 ;
  TCollection *SOld_Coll;
  TVertex **VertexDel, **NewVertices, *temp_Mov, *temp_Mov2, **IsoVertices;
  TJoint *Joint;
  TBoundEdge *tempSlip_Joint;
  TFEVectFunct2D *New_u;
  TFEFunction2D *New_p, *New_c, *New_s, *Surfact_FeFunct;
  TFEFunction2D *GridXDot, *GridYDot;
  TFESpace2D *Velo_Space, *Pres_Space;
  TCollection *IFace_Coll;
  TFESpace2D *velocity_space,  *velocity_space_P1, *velocity_space_P2;
  TFESpace2D *pressure_space, *pressure_space_P1, *pressure_space_P2;
  TFESpace2D *Grid_space, *Grid_space_P1, *Grid_space_P2;
  TFESpace2D *Surfact_space, *Surfact_space_P1, *Surfact_space_P2, *Surfact_Outputspace;
  TFESpace1D *IFaceSurfact_space;
  TIsoInterfaceJoint *IntFaceJoint;
  TBoundEdge *Solid_Joint;
  FE2D FeId;
  TFEDesc2D *FeDesc;
  
  std::ostringstream opts;
  std::ostringstream os;

  os << " ";
  opts << " ";

  struct triangulateio In, Out;

  boolean AllowEdgeRef = (boolean) TDatabase::ParamDB->MESHGEN_ALLOW_EDGE_REF;

  PsBaseName = TDatabase::ParamDB->PSBASENAME;
  GrapeBaseName = TDatabase::ParamDB->GRAPEBASENAME;
  GnuBaseName = TDatabase::ParamDB->GNUBASENAME;
  ReadGrapeBaseName = TDatabase::ParamDB->READGRAPEBASENAME;
  VtkBaseName = TDatabase::ParamDB->VTKBASENAME;

  char NameString[] = "name";
  char UString[] = "u";
  char PString[] = "p";
  char PsiString[] = "psi";
  char TString[] = "T";
  char SurfactString[] = "C";
  char gridString[] = "W";
  char refposString[] = "refpos";
  char auxposString[] = "auxpos";
  char posString[] = "pos";
  char IFaceSString[] = "I_C";
  
  TDatabase::IteratorDB[It_EQ]->SetParam(SurfDomain);
  TDatabase::IteratorDB[It_LE]->SetParam(SurfDomain);
  TDatabase::IteratorDB[It_Finest]->SetParam(SurfDomain);
  TDatabase::IteratorDB[It_Between]->SetParam(SurfDomain);
  TDatabase::IteratorDB[It_OCAF]->SetParam(SurfDomain);

  SOld_Coll = SurfDomain->GetCollection(It_Finest, 0);
  Old_N_Cells = SOld_Coll->GetN_Cells();
  SurfDomain->GetTreeInfo(Surf_OldCellTree, Old_S_N_RootCells);


  // initialize iterators
  TDatabase::IteratorDB[It_EQ]->SetParam(Domain);
  TDatabase::IteratorDB[It_LE]->SetParam(Domain);
  TDatabase::IteratorDB[It_Finest]->SetParam(Domain);
  TDatabase::IteratorDB[It_Between]->SetParam(Domain);
  TDatabase::IteratorDB[It_OCAF]->SetParam(Domain);
    
//======================================================================
// Triangular for grid generation begin
//======================================================================
  cout<< "Start triangulation " << endl;
  area = TDatabase::ParamDB->Area;

  BoundPart = Domain->GetBdPart(0);

  UpdateBound[0]  = (TBdLine*)BoundPart->GetBdComp(0);
  UpdateBound[1]  = (TBdLine*)BoundPart->GetBdComp(1);
  UpdateBound[2]  = (TBdLine*)BoundPart->GetBdComp(2);
  UpdateBound[3]  = (TBdLine*)BoundPart->GetBdComp(3);
  UpdateBound[4]  = (TBdLine*)BoundPart->GetBdComp(4);
  UpdateBound[5]  = (TBdLine*)BoundPart->GetBdComp(5);


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
//   opts<<"VVVV"; // Gives detailed information about what Triangle is doing
//   opts<<'Q'; // Supress all explanation of what Triangle is doing, unless an error occurs
 // opts<<'Y'; // Supress adding vertices on boundary edges
  opts<<'j'; //Jettisons(discard) vertices that are not part of the final triangulation from the output .node file (including duplicate input vertices and vertices ``eaten'' by holes).

//  opts<<"a0.04"; // Imposes a maximum triangle area.
  opts<<"a"<< area; // Imposes a maximum triangle area.

  opts << "nA" << ends;
//   opts<<ends;

  MovBoundVert[0][0]->GetCoords(x, y);
  MovBoundVert[1][0]->GetCoords(tx, ty);
  double Xi[6] = {x, tx, 0., 1., 1., 0.};
  double Yi[6] = {y, ty, 0., 0., 6., 6.};

  N_Hori  = 10;      // number of horoyontal vertices
  N_Verti = 25;

  hE = (Yi[4]-Yi[3])/double(N_Verti);
  cout<< "hE " << hE  << endl;
  N_vert1 = 3*abs(int((Yi[1]-Yi[0])/hE));
  if(N_vert1<4 )N_vert1 =4;
  N_vert2 = abs(int((Yi[2]-Yi[1])/hE));
  if(N_vert2<4 )N_vert2 =4;
  N_vert3 = abs(int((Yi[0]-Yi[5])/hE));
  if(N_vert3<4 )N_vert3 =4;

  N_Verti = N_vert1+N_vert2+N_vert3;       // number of horoyontal vertices
//   cout<< "N_vert1 " << N_vert1 << " N_vert2 " << N_vert2 << " N_vert3 " << N_vert3 << endl;

  N_MovVert[0] = 2*N_Hori + 2*N_Verti;
  I_FaceX = new double[N_MovVert[6]];
  I_FaceY = new double[N_MovVert[6]];

  N_Old_Face_Vert = N_MovVert[6];

  In.numberofpoints = N_MovVert[6]+N_MovVert[0];
  In.pointlist = new double[2*In.numberofpoints];
  In.pointmarkerlist = new int[In.numberofpoints];
  In.numberofpointattributes = 0;

//   start point already given but segment has to be given
  In.numberofsegments = In.numberofpoints;
  In.segmentlist = new int[2*In.numberofsegments];
  In.segmentmarkerlist = new int[In.numberofsegments];
  In.numberofholes = 0;
  In.holelist = NULL;
  In.numberofregions = 0;
  In.regionlist = NULL;

  In_Index = 0;
  CurrComp = 1;

  hi = (Yi[1] - Yi[0])/N_vert1;
  x0 = Xi[0];
  y0 = Yi[0];
  y  = y0;

  end_index = 0;
  // points and segments on the horizontal boundary (marker=1)
 for(i=0;i<N_vert1;i++) // without last point
   {
    In.pointlist[2*In_Index] = x0;
    In.pointlist[2*In_Index+1] = y;
// cout<<" x : "<< In.pointlist[2*In_Index] << " y : "<< In.pointlist[2*In_Index+1] <<endl;
    In.pointmarkerlist[In_Index] = CurrComp;
    In.segmentlist[2*In_Index] = In_Index;
    In.segmentlist[2*In_Index+1] = In_Index+1;
    In.segmentmarkerlist[In_Index] = CurrComp;
    In_Index++;
    y = y0 + double(i+1)*hi;
   }
  CurrComp++;
// cout<<endl;

  hi = (Yi[2] - Yi[1])/N_vert2;
  x0 = Xi[1];
  y0 = Yi[1];
  y  = y0;

 beg_index=In_Index;
  // points and segments on the horizontal boundary (marker=1)
 for(i=0;i<N_vert2;i++) // without last point
   {
    In.pointlist[2*In_Index] = x0;
    In.pointlist[2*In_Index+1] = y;
// cout<<" x : "<< In.pointlist[2*In_Index] << " y : "<< In.pointlist[2*In_Index+1] <<endl;
    In.pointmarkerlist[In_Index] = CurrComp;
    In.segmentlist[2*In_Index] = In_Index;
    In.segmentlist[2*In_Index+1] = In_Index+1;
    In.segmentmarkerlist[In_Index] = CurrComp;
    In_Index++;
    y = y0 + double(i+1)*hi;
   }
  CurrComp++;
// cout<<endl;

  hi = (Xi[3] - Xi[2])/N_Hori;
  x0 = Xi[2];
  y0 = Yi[2];
  x  = Xi[2];
  // points and segments on the horizontal boundary (marker=1)
 for(i=0;i<N_Hori;i++) // without last point
   {
    In.pointlist[2*In_Index] = x;
    In.pointlist[2*In_Index+1] = y0;
// cout<<" x : "<< In.pointlist[2*In_Index] << " y : "<< In.pointlist[2*In_Index+1] <<endl;
    In.pointmarkerlist[In_Index] = CurrComp;
    In.segmentlist[2*In_Index] = In_Index;
    In.segmentlist[2*In_Index+1] = In_Index+1;
    In.segmentmarkerlist[In_Index] = CurrComp;
    In_Index++;
    x = x0 + double(i+1)*hi;
   }
  CurrComp++;
// cout<<endl;


  hi = (Yi[4] - Yi[3])/N_Verti;
  x0 = Xi[3];
  y0 = Yi[3];
  y  = Yi[3];
  // points and segments on the horizontal boundary (marker=1)
  for(i=0;i<N_Verti;i++) // without last point
   {
    In.pointlist[2*In_Index] = x0;
    In.pointlist[2*In_Index+1] = y;
// cout<<" x : "<< In.pointlist[2*In_Index] << " y : "<< In.pointlist[2*In_Index+1] <<endl;
    In.pointmarkerlist[In_Index] = CurrComp;
    In.segmentlist[2*In_Index] = In_Index;
    In.segmentlist[2*In_Index+1] = In_Index+1;
    In.segmentmarkerlist[In_Index] = CurrComp;
    In_Index++;
    y = y0 + double(i+1)*hi;
   }
  CurrComp++;
// cout<<endl;

  hi = (Xi[5] - Xi[4])/N_Hori;
  x0 = Xi[4];
  y0 = Yi[4];
  x  = Xi[4];
  // points and segments on the horizontal boundary (marker=1)
 for(i=0;i<N_Hori;i++) // without last point
   {
    In.pointlist[2*In_Index] = x;
    In.pointlist[2*In_Index+1] = y0;
// cout<<" x : "<< In.pointlist[2*In_Index] << " y : "<< In.pointlist[2*In_Index+1] <<endl;
    In.pointmarkerlist[In_Index] = CurrComp;
    In.segmentlist[2*In_Index] = In_Index;
    In.segmentlist[2*In_Index+1] = In_Index+1;
    In.segmentmarkerlist[In_Index] = CurrComp;
    In_Index++;
    x = x0 + double(i+1)*hi;
   }
  CurrComp++;
// cout<<endl;

  hi = (Yi[0] - Yi[5])/N_vert3;
  x0 = Xi[5];
  y0 = Yi[5];
  y  = Yi[5];
  // points and segments on the horizontal boundary (marker=1)
 for(i=0;i<N_vert3;i++) // without last point
   {
    In.pointlist[2*In_Index] = x0;
    In.pointlist[2*In_Index+1] = y;
// cout<<" x : "<< In.pointlist[2*In_Index] << " y : "<< In.pointlist[2*In_Index+1] <<endl;
    In.pointmarkerlist[In_Index] = CurrComp;
    In.segmentlist[2*In_Index] = In_Index;
    In.segmentlist[2*In_Index+1] = In_Index+1;
    In.segmentmarkerlist[In_Index] = CurrComp;
    In_Index++;
    y = y0 + double(i+1)*hi;
   }
  CurrComp++;
// cout<<endl;
  In.segmentlist[2*(In_Index-1)+1] = 0;

  // points and segments on the free boundary (marker=2)
  memcpy(I_FaceX, FreePts[0], N_MovVert[6]*SizeOfDouble);
  memcpy(I_FaceY, FreePts[1], N_MovVert[6]*SizeOfDouble);

  In.segmentlist[2*In_Index] = beg_index;
  In.segmentlist[2*In_Index+1] = In_Index;
  In.segmentmarkerlist[In_Index] = CurrComp;

  for(i=1;i<N_MovVert[6];i++) // without last point
   {
//     MovBoundVert[1][i]->GetCoords(I_FaceX[i], I_FaceY[i]);
    In.pointlist[2*In_Index] = I_FaceX[i];
    In.pointlist[2*In_Index+1] = I_FaceY[i];
//     cout<<In_Index << " x :" << I_FaceX[i]<< " -----------------y: " <<I_FaceY[i]<< endl;
//   cout<< " theta "<<(180/Pi)*atan2(In.pointlist[2*In_Index+1]-0.5, In.pointlist[2*In_Index])<<endl;

    In.pointmarkerlist[In_Index] = CurrComp;
    In.segmentlist[2*(In_Index+1)] = In_Index;
    In.segmentlist[2*(In_Index+1)+1] = In_Index+1;
    In.segmentmarkerlist[(In_Index+1)] = CurrComp;

    In_Index++;
   }


  In.pointlist[2*In_Index] = -100.;
  In.pointlist[2*(In_Index)+1] = -100.;
  In.pointmarkerlist[In_Index] = CurrComp;

  In.segmentlist[2*(In_Index)+1] = end_index;


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
//   cout<< "Domain triangulated (but not yet assigned in MooNMD)" << endl;
 /*
  for(i=0;i<Out.numberofpoints;i++)
      OutPut(i<<' '<<Out.pointmarkerlist[i]<<' '<<
       Out.pointlist[2*i]<<' '<<Out.pointlist[2*i+1]<<endl);
 */

  Old_Coll = Domain->GetCollection(It_Finest, 0);
  Old_N_Cells = Old_Coll->GetN_Cells();
  Domain->GetTreeInfo(Old_CellTree, Old_N_RootCells);
   OutPut("Number of rootcells: "<<Old_N_RootCells<<endl);
   if(Old_N_Cells!=Old_N_RootCells) exit(-1);
  // remove all existing vertices and joints
  VertexDel = new TVertex*[3*Old_N_RootCells];

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


 // Bound startx, starty, x length and y length
 // Solid Bound startx, starty, x length and y length
  UpdateBound[0]->SetParams(Xi[0], Yi[0], Xi[1]-Xi[0],Yi[1]-Yi[0]);
  UpdateBound[1]->SetParams(Xi[1], Yi[1], Xi[2]-Xi[1],Yi[2]-Yi[1]);
  UpdateBound[2]->SetParams(Xi[2], Yi[2], Xi[3]-Xi[2],Yi[3]-Yi[2]);
  UpdateBound[3]->SetParams(Xi[3], Yi[3], Xi[4]-Xi[3],Yi[4]-Yi[3]);
  UpdateBound[4]->SetParams(Xi[4], Yi[4], Xi[5]-Xi[4],Yi[5]-Yi[4]);
  UpdateBound[5]->SetParams(Xi[5], Yi[5], Xi[0]-Xi[5],Yi[0]-Yi[5]);

//    UpdateSlipBound->SetParams(Lx, Ly, Rx-Lx, Ry-Ly);
//    cout<<"Lx : " << Lx << "  Ly : "<<Ly <<" Rx-Lx : " << Rx-Lx<< "  Ry-Rx : "<<Ry-Ly<< endl;
//    cout<<"Lx : " << Lx << "  Ly : "<<Ly <<" Rx : " << Rx<< "  Ry : "<<Ry<< endl;

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
//     os.seekp(std::ios::beg);
//     os << "Domain" << ".ps" << ends;
//     Domain->PS(os.str().c_str(),It_Finest,0);



  ORDER = 0;
  VSP = TDatabase::ParamDB->VELOCITY_SPACE;
if (abs(VSP) > 20)
  {ORDER = abs(VSP) - 20;}
else if ( abs(VSP) > 10)
    {ORDER = abs(VSP) - 10;}
  else ORDER = abs(VSP);

// delete subspaces all other spaces are needed for interpolation
  for(i=0;i<N_FESpaces_All;i++)
   for(j=1;j<3;j++)
     delete FeSpaces[i][j];

// delete grid space
  delete FeSpaces[2][0];

  delete IFaceFeSpaces[0];

   for(j=1;j<3;j++)
    {
     delete  [] Sol[0][j];
     delete  [] Sol[3][j];
    }

   for(j=0;j<3;j++)
     delete  [] Sol[2][j];

   for(j=0;j<3;j++)
    {
     delete  [] Rhs[0][j];
     delete  [] Rhs[2][j];
//      delete  [] Rhs[3][j];
    }

// interface surfactant subspaces
    delete  [] Sol[N_FESpaces_All][0];
    delete  [] Rhs[N_FESpaces_All][0];


//    for(j=1;j<3;j++)
//     delete VeloVect[0][j];

//    for(j=0;j<3;j++)
//     delete VeloVect[1][j];

  for(i=0;i<6;i++)
   for(j=1;j<3;j++)
    delete FeFunct[i][j];

  delete FeFunct[3][0];
  delete FeFunct[4][0];

  for(j=0;j<3;j++)
   {
//     delete SqrStruct[0][j];
//     delete (TSquareStructure*)SqrStruct[0][j];
    delete (TSquareStructure*)SqrStruct[2][j];
    delete (TSquareStructure*)SqrStruct[3][j];
   }

  for(i=0;i<2;i++)
   for(j=0;j<3;j++)
    delete (TStructure*)Struct[i][j];

  for(k=0;k<10;k++)
   delete SqMat[0][0][k];

  for(j=1;j<3;j++)
  for(k=0;k<8;k++)
   delete SqMat[0][j][k];

  for(j=0;j<3;j++)
   for(k=0;k<4;k++)
   delete Mat[j][k];

  for(j=0;j<3;j++)
  for(k=0;k<4;k++)
   delete SqMat[2][j][k];

  for(j=0;j<3;j++)
  for(k=0;k<2;k++)
   delete SqMat[3][j][k];

   delete SqMat_IFace[0];
   delete SqMat_IFace[1];
  


// *************************************************************************
//  split the doamin in to interior and exterior begin
// *************************************************************************

  Coll=Domain->GetCollection(It_Finest, 0);
  N_Cells = Coll->GetN_Cells();
  N_Cells_P1 = 0;
  N_Cells_P2 = 0;


  for(i=0;i<N_Cells;i++)
   {
    cell = Coll->GetCell(i);
    cell->SetGlobalCellNo(i);     
    
    x_mid = 0.;  y_mid = 0.;

    for(j=0;j<3;j++)
     {
      cell->GetVertex(j)->GetCoords(x, y);
      x_mid +=x;
      y_mid +=y;
     } // for j

      x_mid /=3.;
      y_mid /=3.;
      hmin = 1000.0;
// find point P on interface for the i'th cell mid point X
// such that ||P - X|| is minimum
      for(k=0;k<N_Old_Face_Vert;k++)
       {
        temp = sqrt( (I_FaceX[k]-x_mid)*(I_FaceX[k]-x_mid) + (I_FaceY[k]-y_mid)*(I_FaceY[k]-y_mid) );
        if(temp<hmin)
         {
          hmin=temp;  T_a = I_FaceX[k];  T_b = I_FaceY[k];
          x_short = T_a;  y_short = T_b;
          sp_No = k;
         }
       } // for k

// find next shortest point 
// i.e other shortest point of a spline containing (T_a,T_b)
// previous point
   if(sp_No==0)
    sp_no0 = N_Old_Face_Vert-1;
   else
    sp_no0 = sp_No-1;

    temp0 = sqrt( (I_FaceX[sp_no0]-x_mid)*(I_FaceX[sp_no0]-x_mid)
             + (I_FaceY[sp_no0]-y_mid)*(I_FaceY[sp_no0]-y_mid) );

// next point 
   if(sp_No==N_Old_Face_Vert-1)
    sp_no1 = 0;
   else
    sp_no1 = sp_No+1;

    temp1 = sqrt( (I_FaceX[sp_no1]-x_mid)*(I_FaceX[sp_no1]-x_mid)
                + (I_FaceY[sp_no1]-y_mid)*(I_FaceY[sp_no1]-y_mid) );

   if( temp0 < temp1)
     {
     T_a = I_FaceX[sp_no0];  T_b = I_FaceY[sp_no0];
     C_x = I_FaceX[sp_No];  C_y = I_FaceY[sp_No];
     }
   else 
     {
     T_a = I_FaceX[sp_No];  T_b = I_FaceY[sp_No];
     C_x = I_FaceX[sp_no1];  C_y = I_FaceY[sp_no1];
     }

     tx = x_mid - x_short; //x distance between the point and the shortest distance point on interface
     ty = y_mid - y_short; //y distance between the point and the shortest distance point on interface

     nx = -(C_y - T_b); // (-) normal at (T_a,T_b) pointing into the inner domain
     ny = -(T_a - C_x); // (-) normal at (T_a,T_b) pointing into the inner domain

     Pos_Indicator = tx*nx + ty*ny;

   if(Pos_Indicator > 0.)
    {
//  cell is in inner domain
//  cout<< "Inner cell " << i << endl;
     cell->SetPhase_ID(0);
     N_Cells_P1++;
    }
   else if(Pos_Indicator < 0.)
    {
//  cell is in outer domain
//  cout<< "outer cell " << i << endl;
     cell->SetPhase_ID(1);
     N_Cells_P2++;
    }
   else
    {
     cout<< " error in identifying phase check remesh2d.c" << endl;
     exit(0);
    }

  }// for i

  OutPut( "Number of inner cells " << N_Cells_P1 << endl);
  OutPut( "Number of outer cells " << N_Cells_P2 << endl);

  delete [] Coll_Cells[1];
  delete [] Coll_Cells[2];

  delete [] GlobalCell_Index[0];
  delete [] GlobalCell_Index[1];

  GlobalCell_Index[0] = new int[N_Cells_P1];
  GlobalCell_Index[1] = new int[N_Cells_P2];

  Coll_Cells[1] = new TBaseCell *[N_Cells_P1];
  Coll_Cells[2] = new TBaseCell *[N_Cells_P2];


  m1 = 0;
  m2 = 0;
  for(i=0;i<N_Cells;i++)
   {
    cell = Coll->GetCell(i);
    ID = cell->GetPhase_ID();

     if(ID == 0)
      {
//    inner cell
       Coll_Cells[1][m1] = cell;
       cell->SetLocalCellNo(m1);
       GlobalCell_Index[0][m1] = i;
       m1++;

      }
    else
      {
//    outer cell
       Coll_Cells[2][m2] = cell;
       cell->SetLocalCellNo(m2);
       GlobalCell_Index[1][m2] = i;
       m2++;
      }
   }

//   Coll_Multi[0] = new TCollection(N_Cells_P1, Coll_Cells[1]);
//   Coll_Multi[1] = new TCollection(N_Cells_P2, Coll_Cells[2]);
  Coll_Multi[0]->Replace_Coll(N_Cells_P1, Coll_Cells[1]);
  Coll_Multi[1]->Replace_Coll(N_Cells_P2, Coll_Cells[2]);

  Coll_P1 = Coll_Multi[0];
  Coll_P2 = Coll_Multi[1];

  delete [] I_FaceX;
  delete [] I_FaceY;

//   delete [] N_List[0];
//   delete [] N_List[1];

// *************************************************************************
//  split the doamin in to interior and exterior -- end
//*************************************************************************
// *************************************************************************
//   FESpaces memory allocation
// *************************************************************************

    Coll=Domain->GetCollection(It_Finest, 0);
    cout << endl << endl;
    TDatabase::IteratorDB[It_EQ]->SetParam(SurfDomain);
    TDatabase::IteratorDB[It_LE]->SetParam(SurfDomain);
    TDatabase::IteratorDB[It_Finest]->SetParam(SurfDomain);
    TDatabase::IteratorDB[It_Between]->SetParam(SurfDomain);
    TDatabase::IteratorDB[It_OCAF]->SetParam(SurfDomain);

//  list of outer phase cells containing interafce
//     N_List[0] = new int[N_Cells];   // Cell_No
//     N_List[1] = new int[N_Cells];   //Joint_No
    Domain2DSurf_2Phase(Coll, SurfDomain, N_List);

    IFace_Coll = SurfDomain->GetCollection(It_Finest, 0);
    N_IFaceCells= IFace_Coll->GetN_Cells();
    cout<< " N_IFaceCells " << N_IFaceCells <<endl;

    
    FE1D_List = new FE1D[N_IFaceCells];
    for(j=0;j<N_IFaceCells;j++)
     FE1D_List[j] = FE1D(TDatabase::ParamDB->ANSATZ_ORDER);

    TDatabase::IteratorDB[It_EQ]->SetParam(Domain);
    TDatabase::IteratorDB[It_LE]->SetParam(Domain);
    TDatabase::IteratorDB[It_Finest]->SetParam(Domain);
    TDatabase::IteratorDB[It_Between]->SetParam(Domain);
    TDatabase::IteratorDB[It_OCAF]->SetParam(Domain);

/* velocity space and pressure spaces */
    GetVelocityAndPressureSpace(Coll,BoundCondition,
                                mortarcoll, velocity_space,
                                pressure_space, &pressure_space_code,
                                TDatabase::ParamDB->VELOCITY_SPACE,
                                TDatabase::ParamDB->PRESSURE_SPACE);

    velocity_space_code =   TDatabase::ParamDB->VELOCITY_SPACE;
    TDatabase::ParamDB->INTERNAL_PRESSURE_SPACE = pressure_space_code;

    GetVelocityAndPressureSpace(Coll_P1,BoundCondition,
                                mortarcoll, velocity_space_P1,
                                pressure_space_P1, &pressure_space_code,
                                TDatabase::ParamDB->VELOCITY_SPACE,
                                TDatabase::ParamDB->PRESSURE_SPACE);
    GetVelocityAndPressureSpace(Coll_P2,BoundCondition,
                                mortarcoll, velocity_space_P2,
                                pressure_space_P2, &pressure_space_code,
                                TDatabase::ParamDB->VELOCITY_SPACE,
                                TDatabase::ParamDB->PRESSURE_SPACE);

//    FeSpaces[0][0] = velocity_space;
   FeSpaces[0][1] = velocity_space_P1;
   FeSpaces[0][2] = velocity_space_P2;
   N_DOFs[0][0] =  velocity_space->GetN_DegreesOfFreedom();
   N_DOFs[0][1] =  velocity_space_P1->GetN_DegreesOfFreedom();
   N_DOFs[0][2] =  velocity_space_P2->GetN_DegreesOfFreedom();

//    FeSpaces[1][0] = pressure_space;
   FeSpaces[1][1] = pressure_space_P1;
   FeSpaces[1][2] = pressure_space_P2;
   N_DOFs[1][0] =  pressure_space->GetN_DegreesOfFreedom();
   N_DOFs[1][1] =  pressure_space_P1->GetN_DegreesOfFreedom();
   N_DOFs[1][2] =  pressure_space_P2->GetN_DegreesOfFreedom();

   GlobalNumbers[0][0] = velocity_space->GetGlobalNumbers();
   BeginIndex[0][0]    = velocity_space->GetBeginIndex();
   GlobalNumbers[0][1] = FeSpaces[0][1]->GetGlobalNumbers();
   BeginIndex[0][1]    = FeSpaces[0][1]->GetBeginIndex();
   GlobalNumbers[0][2] = FeSpaces[0][2]->GetGlobalNumbers();
   BeginIndex[0][2]    = FeSpaces[0][2]->GetBeginIndex();

/* grid spaces */
   Grid_space = new TFESpace2D(Coll, NameString, PsiString, GridBoundCondition,
                                1, NULL);
   Grid_space_P1 = new TFESpace2D(Coll_P1, NameString, PsiString, GridBoundCondition,
                                1, NULL);
   Grid_space_P2 = new TFESpace2D(Coll_P2, NameString, PsiString, GridBoundCondition,
                                1, NULL);

   FeSpaces[2][0] = Grid_space;
   FeSpaces[2][1] = Grid_space_P1;
   FeSpaces[2][2] = Grid_space_P2;
   N_DOFs[2][0] =  Grid_space->GetN_DegreesOfFreedom();
   N_DOFs[2][1] =  Grid_space_P1->GetN_DegreesOfFreedom();
   N_DOFs[2][2] =  Grid_space_P2->GetN_DegreesOfFreedom();

   Bound_DOFs[2][0] = N_DOFs[2][0] - FeSpaces[2][0]->GetN_Inner();
   Bound_DOFs[2][1] = N_DOFs[2][1] - FeSpaces[2][1]->GetN_Inner();
   Bound_DOFs[2][2] = N_DOFs[2][2] - FeSpaces[2][2]->GetN_Inner();

   GlobalNumbers[2][0] = FeSpaces[2][0]->GetGlobalNumbers();
   BeginIndex[2][0]    = FeSpaces[2][0]->GetBeginIndex();
   GlobalNumbers[2][1] = FeSpaces[2][1]->GetGlobalNumbers();
   BeginIndex[2][1]    = FeSpaces[2][1]->GetBeginIndex();
   GlobalNumbers[2][2] = FeSpaces[2][2]->GetGlobalNumbers();
   BeginIndex[2][2]    = FeSpaces[2][2]->GetBeginIndex();

// soluble surfactant in the outer phase 
   Surfact_space = new TFESpace2D(Coll, NameString, SurfactString, SurfactBoundCondition,
                                   ContP_USpace, TDatabase::ParamDB->ANSATZ_ORDER, NULL);
   Surfact_space_P1 = new TFESpace2D(Coll_P1, NameString, SurfactString, SurfactBoundCondition,
                                   ContP_USpace, TDatabase::ParamDB->ANSATZ_ORDER, NULL);
   Surfact_space_P2 = new TFESpace2D(Coll_P2, NameString, SurfactString, SurfactBoundCondition,
                                   ContP_USpace, TDatabase::ParamDB->ANSATZ_ORDER, NULL);

   FeSpaces[3][1] = Surfact_space_P1;
   FeSpaces[3][2] = Surfact_space_P2;
   N_DOFs[3][0] =  Surfact_space->GetN_DegreesOfFreedom();
   N_DOFs[3][1] =  Surfact_space_P1->GetN_DegreesOfFreedom();
   N_DOFs[3][2] =  Surfact_space_P2->GetN_DegreesOfFreedom();

   GlobalNumbers[3][0] = Surfact_space->GetGlobalNumbers();
   BeginIndex[3][0]    = Surfact_space->GetBeginIndex();
   GlobalNumbers[3][1] = FeSpaces[3][1]->GetGlobalNumbers();
   BeginIndex[3][1]    = FeSpaces[3][1]->GetBeginIndex();
   GlobalNumbers[3][2] = FeSpaces[3][2]->GetGlobalNumbers();
   BeginIndex[3][2]    = FeSpaces[3][2]->GetBeginIndex();


// Surfactant on the interface
   IFaceSurfact_space = new TFESpace1D(IFace_Coll , IFaceSString, IFaceSString, FE1D_List);

   IFaceFeSpaces[0] = IFaceSurfact_space;
   N_DOFs[N_FESpaces_All][0] =  IFaceFeSpaces[0]->GetN_DegreesOfFreedom();

   Surfact_Outputspace = new TFESpace2D(Coll, NameString, SurfactString, SurfactBoundCondition,
                                   ContP_USpace, TDatabase::ParamDB->ANSATZ_ORDER, NULL);

   N_DOFs[4][0] =  Surfact_Outputspace->GetN_DegreesOfFreedom();
// *************************************************************************
// adding aditional vertices for the insinterface/isoboundedge edge
// *************************************************************************
   N_Cells = Coll->GetN_Cells();

   Old_N_IsoEdges =  N_MovVert[6];
   
   for(i=0;i<8;i++)
    N_MovVert[i] = 0;

   for(j=0;j<N_Cells;j++)
    {
     Me = Coll->GetCell(j);
     k = Me->GetN_Edges();

     for(l=0;l<k;l++)
      {
       Joint = Me->GetJoint(l);

       if(Joint->GetType() == BoundaryEdge)
        {
         Solid_Joint = (TBoundEdge *)Joint;
         BoundComp = Solid_Joint->GetBoundComp();
         comp=BoundComp->GetID();
         N_MovVert[comp]++;
        }
       else if(Joint->GetType() == IsoInterfaceJoint && Me->GetPhase_ID()==0)
        {
         IntFaceJoint = (TIsoInterfaceJoint *)Joint;
         BoundComp = IntFaceJoint->GetBoundComp();
         comp=BoundComp->GetID();
         N_MovVert[comp]++;

         Me->GetVertex(l)->GetCoords(TX[0], TY[0]);
         Me->GetVertex((l+1) % k)->GetCoords(TX[1], TY[1]);
         IntFaceJoint->GeneratemidVert(ORDER-1, TX, TY);
// cout<< "TX[0] " <<TX[0] << " TX[1] " << TX[1] << endl;
//       IntFaceJoint->GenerateVertices(ORDER-1); // cannot be used, sinece no parametrised interface available
        }
      } //  for(l=0;
     } //  for(j=0;j<N_

   for(i=0;i<7;i++)
    cout<< " Number of Boundary edge in comp " << i << " is " << N_MovVert[i]<<endl;

// store boundary vertices for later use
// assumed that we have only 7 boundary components
   for(i=0;i<7;i++)
    {
     delete [] Slip_Joint[i];
     delete [] MovBoundVert[i];

     Slip_Joint[i] = new TBoundEdge*[N_MovVert[i]];
     MovBoundVert[i] = new TVertex*[N_MovVert[i]];
     mm[i] = 0;
    }

    delete [] N_List[2];
    delete [] N_List[3];

    N_List[2] = new int[N_MovVert[6]]; // cell list
    N_List[3] = new int[N_MovVert[6]];  // edge list
//     Coll_Cells[0] = new TBaseCell*[10*N_MovVert[6]];

   for(j=0;j<N_Cells;j++)
    {
     Me = Coll->GetCell(j);
     k = Me->GetN_Edges();

     for(l=0;l<k;l++)
      {
       Joint = Me->GetJoint(l);
       if(Joint->GetType() == BoundaryEdge)
        {
         Solid_Joint = (TBoundEdge *)Joint;
         BoundComp = Solid_Joint->GetBoundComp();
         comp=BoundComp->GetID();

         MovBoundVert[comp][ mm[comp] ] = Me->GetVertex(l);
         Slip_Joint[comp][ mm[comp] ] = (TBoundEdge *)Joint;

         mm[comp]++;
        }
       else if(Joint->GetType() == IsoInterfaceJoint && Me->GetPhase_ID()==0)
        {
         IntFaceJoint = (TIsoInterfaceJoint *)Joint;
         BoundComp = IntFaceJoint->GetBoundComp();
         comp=BoundComp->GetID();
         Coll_Cells[0][mm[comp] ] = Me;
         MovBoundVert[comp][mm[comp] ] = Me->GetVertex(l);
         N_List[2][mm[comp] ] = j;
         N_List[3][mm[comp] ] = l;

         mm[comp]++;
        }

      } //  for(l=0;
     } //  for(j=0;j<N_


   // sort with X0 as the first vertex ordinate if indicator==0 (or) ....
   SortIsoVert_Gen(Coll_Cells[0], MovBoundVert[6], N_List[2], N_List[3],
                   N_MovVert[6], Xi[0], Yi[0], 0);

   OutPut(" Old_N_IsoEdges " << Old_N_IsoEdges <<" N_MovVert[6] " <<  N_MovVert[6] <<endl);
   
   if(Old_N_IsoEdges==N_MovVert[6])
    {
     // no additional vertices are added during remeshing on the interface, so copy the old isovertces
     m=0;
     comp=0;
     for(j=0;j<Old_N_IsoEdges;j++)
      {
       Me = Coll_Cells[0][j];
       l = N_List[3][j];
       k = Me->GetN_Edges();
 
       Me->GetVertex(l)->GetCoords(TX[0], TY[0]);
       Me->GetVertex((l+1) % k)->GetCoords(TX[1], TY[1]);
       
       temp0 = ((Intpol_Coord[2*m] - TX[0])*(Intpol_Coord[2*m] - TX[0]) + (Intpol_Coord[2*m+1] - TY[0])*(Intpol_Coord[2*m+1] - TY[0]) );
       comp = m+ORDER;
       temp1 = ((Intpol_Coord[2*comp] - TX[1])*(Intpol_Coord[2*comp] - TX[1]) + (Intpol_Coord[2*comp+1] - TY[1])*(Intpol_Coord[2*comp+1] - TY[1]) );
       
       if(temp0<1.e-8 && temp1<1.e-8 ) // update the isovertex from the old value
        {
         m++;
         Joint = Coll_Cells[0][j]->GetJoint(N_List[3][j]);
         IntFaceJoint = (TIsoInterfaceJoint *)Joint;
         IsoVertices = IntFaceJoint->GetVertices();
         for(k=0; k<ORDER-1; k++)
         {
	   IsoVertices[k]->SetCoords(Intpol_Coord[2*m], Intpol_Coord[2*m+1]);
	   m++;
         }
        }
       else
        {
         // only second order conforming elements implimented
         cout<< " No match in Old_N_IsoEdges  and N_MovVert[6] in Remesh "<<endl;
         
         cout <<j <<  "X: " <<  Intpol_Coord[2*m] << " " << FreePts[0][j] <<  " " << TX[0] <<" Y: " <<  Intpol_Coord[2*m+1] << " " << FreePts[1][j] <<" " << TY[0] << endl;

         cout << j << "X: " <<  Intpol_Coord[2*comp] << " " << FreePts[0][j+1] <<  " " << TX[1] <<" Y: " <<  Intpol_Coord[2*comp+1] << " " << FreePts[1][j+1] <<" " << TY[1] << endl;

         exit(0);
        }
      } // for(j=0;j<Old_N_IsoEdges-1;j++) 
    }//if(Old_N_IsoEdges==N_MovVert[6])
   else
    {
     cout<< " New vertices are added on the interface in Remesh "<<endl;
     exit(0);
    }
  
//    cout<< " SLPX0 " << Xi[0] <<" SLPY " <<  Yi[0] <<endl;
// 
//      for(k=0;k<N_MovVert[6];k++)
//       {
//        MovBoundVert[6][k]->GetCoords(x, y);
//        cout<<k << " x :" << x << " --y: " <<y;
//        cout<< " theta "<<(180./Pi)*atan2(y-0.5, x)<<endl;
//        }
//   cout<<endl;


// sort boundary vertces
    for(k=0;k<N_MovVert[0]-1;k++)
      {
      for(l=k+1;l<N_MovVert[0];l++)
      {
        MovBoundVert[0][k]->GetCoords(x, y);
	MovBoundVert[0][l]->GetCoords(tx, ty);
	if(y < ty)
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
//      for(k=0;k<N_MovVert[0];k++)
//       {
//        MovBoundVert[0][k]->GetCoords(x, y);
//        cout<< " SLPX0 " << x<<" SLPY " << y<<endl;
//        }
//   cout<<endl;

    for(k=0;k<N_MovVert[1]-1;k++)
      {
      for(l=k+1;l<N_MovVert[1];l++)
      {
        MovBoundVert[1][k]->GetCoords(x, y);
	MovBoundVert[1][l]->GetCoords(tx, ty);
	if(y < ty)
	 {
	  temp_Mov = MovBoundVert[1][k];
          MovBoundVert[1][k] = MovBoundVert[1][l];
          MovBoundVert[1][l] = temp_Mov;

	  tempSlip_Joint = Slip_Joint[1][k];
	  Slip_Joint[1][k] = Slip_Joint[1][l];
	  Slip_Joint[1][l] = tempSlip_Joint;
	 }
        }
       }
//      for(k=0;k<N_MovVert[1];k++)
//       {
//        MovBoundVert[1][k]->GetCoords(x, y);
//        cout<< " SLPX1 " << x<<" SLPY " << y<<endl;
//        }
//   cout<<endl;

    for(k=0;k<N_MovVert[2]-1;k++)
      {
      for(l=k+1;l<N_MovVert[2];l++)
      {
        MovBoundVert[2][k]->GetCoords(x, y);
	MovBoundVert[2][l]->GetCoords(tx, ty);
	if(tx < x)
	 {
	  temp_Mov = MovBoundVert[2][k];
          MovBoundVert[2][k] = MovBoundVert[2][l];
          MovBoundVert[2][l] = temp_Mov;

	  tempSlip_Joint = Slip_Joint[2][k];
	  Slip_Joint[2][k] = Slip_Joint[2][l];
	  Slip_Joint[2][l] = tempSlip_Joint;
	 }
        }
       }
//      for(k=0;k<N_MovVert[2];k++)
//       {
//        MovBoundVert[2][k]->GetCoords(x, y);
//        cout<< " SLPX " << x<<" SLPY " << y<<endl;
//        }
//   cout<<endl;
    for(k=0;k<N_MovVert[3]-1;k++)
      {
      for(l=k+1;l<N_MovVert[3];l++)
      {
        MovBoundVert[3][k]->GetCoords(x, y);
	MovBoundVert[3][l]->GetCoords(tx, ty);
	if(ty < y)
	 {
	  temp_Mov = MovBoundVert[3][k];
          MovBoundVert[3][k] = MovBoundVert[3][l];
          MovBoundVert[3][l] = temp_Mov;

	  tempSlip_Joint = Slip_Joint[3][k];
	  Slip_Joint[3][k] = Slip_Joint[3][l];
	  Slip_Joint[3][l] = tempSlip_Joint;
	 }
        }
       }
//      for(k=0;k<N_MovVert[3];k++)
//       {
//        MovBoundVert[3][k]->GetCoords(x, y);
//        cout<< " SLPX1 " << x<<" SLPY " << y<<endl;
//        }
// 
//   cout<<endl;
    for(k=0;k<N_MovVert[4]-1;k++)
      {
      for(l=k+1;l<N_MovVert[4];l++)
      {
        MovBoundVert[4][k]->GetCoords(x, y);
	MovBoundVert[4][l]->GetCoords(tx, ty);
	if(tx > x)
	 {
	  temp_Mov = MovBoundVert[4][k];
          MovBoundVert[4][k] = MovBoundVert[4][l];
          MovBoundVert[4][l] = temp_Mov;

	  tempSlip_Joint = Slip_Joint[4][k];
	  Slip_Joint[4][k] = Slip_Joint[4][l];
	  Slip_Joint[4][l] = tempSlip_Joint;
	 }
        }
       }
//      for(k=0;k<N_MovVert[4];k++)
//       {
//        MovBoundVert[4][k]->GetCoords(x, y);
//        cout<< " SLPX1 " << x<<" SLPY " << y<<endl;
//        }
//   cout<<endl;
    for(k=0;k<N_MovVert[5]-1;k++)
      {
      for(l=k+1;l<N_MovVert[5];l++)
      {
        MovBoundVert[5][k]->GetCoords(x, y);
	MovBoundVert[5][l]->GetCoords(tx, ty);
	if(y < ty)
	 {
	  temp_Mov = MovBoundVert[5][k];
          MovBoundVert[5][k] = MovBoundVert[5][l];
          MovBoundVert[5][l] = temp_Mov;

	  tempSlip_Joint = Slip_Joint[5][k];
	  Slip_Joint[5][k] = Slip_Joint[5][l];
	  Slip_Joint[5][l] = tempSlip_Joint;
	 }
        }
       }
//      for(k=0;k<N_MovVert[5];k++)
//       {
//        MovBoundVert[5][k]->GetCoords(x, y);
//        cout<< " SLPX1 " << x<<" SLPY " << y<<endl;
//        }

// *************************************************************************
//   Matrices memory allocation
// *************************************************************************
// cout<< " test 1 " <<endl;
   /* velocity */
    SqrStruct[0][0] = new TSquareStructure2D(velocity_space);
    SqrStruct[0][1] = new TSquareStructure2D(FeSpaces[0][1]);
    SqrStruct[0][2] = new TSquareStructure2D(FeSpaces[0][2]);
    for(i=0;i<3;i++)
     SqrStruct[0][i]->Sort();

   /* velo - press coupled  B */
    Struct[0][0] = new TStructure2D(pressure_space, velocity_space);
    Struct[0][1] = new TStructure2D(FeSpaces[1][1], FeSpaces[0][1]);
    Struct[0][2] = new TStructure2D(FeSpaces[1][2], FeSpaces[0][2]);

   /* velo - press coupled  BT */
    Struct[1][0] = new TStructure2D(velocity_space, pressure_space);
    Struct[1][1] = new TStructure2D(FeSpaces[0][1], FeSpaces[1][1]);
    Struct[1][2] = new TStructure2D(FeSpaces[0][2], FeSpaces[1][2]);

 /* grid matrices */
    SqrStruct[2][0] = new TSquareStructure2D(FeSpaces[2][0]);
    SqrStruct[2][1] = new TSquareStructure2D(FeSpaces[2][1]);
    SqrStruct[2][2] = new TSquareStructure2D(FeSpaces[2][2]);


 /* surfactant matrices */
    SqrStruct[3][0] = new TSquareStructure2D(Surfact_space);
    SqrStruct[3][1] = new TSquareStructure2D(FeSpaces[3][1]);
    SqrStruct[3][2] = new TSquareStructure2D(FeSpaces[3][2]);
    for(i=0;i<3;i++)
     SqrStruct[3][i]->Sort();
    
 /* interface surfactant matrices */
    IFaceStruct[0] = new TSquareStructure1D(IFaceFeSpaces[0]);
    IFaceStruct[0]->Sort();

 /* velocity */
    SqMat[0][0][0] = new TSquareMatrix2D(SqrStruct[0][0]); // A11
    SqMat[0][0][1] = new TSquareMatrix2D(SqrStruct[0][0]); // A12
    SqMat[0][0][2] = new TSquareMatrix2D(SqrStruct[0][0]); // A21
    SqMat[0][0][3] = new TSquareMatrix2D(SqrStruct[0][0]); // A22
    SqMat[0][0][4] = new TSquareMatrix2D(SqrStruct[0][0]); // M11
    SqMat[0][0][5] = new TSquareMatrix2D(SqrStruct[0][0]); // M12
    SqMat[0][0][6] = new TSquareMatrix2D(SqrStruct[0][0]); // M21
    SqMat[0][0][7] = new TSquareMatrix2D(SqrStruct[0][0]); // M22
    SqMat[0][0][8] = new TSquareMatrix2D(SqrStruct[0][0]); // F11 interfaceint
    SqMat[0][0][9] = new TSquareMatrix2D(SqrStruct[0][0]); // F22 interfaceint

//   Defect = Defect_NSE4;
/* velocity phase_1*/
    SqMat[0][1][0] = new TSquareMatrix2D(SqrStruct[0][1]); // A11
    SqMat[0][1][1] = new TSquareMatrix2D(SqrStruct[0][1]); // A12
    SqMat[0][1][2] = new TSquareMatrix2D(SqrStruct[0][1]); // A21
    SqMat[0][1][3] = new TSquareMatrix2D(SqrStruct[0][1]); // A22
    SqMat[0][1][4] = new TSquareMatrix2D(SqrStruct[0][1]); // M11
    SqMat[0][1][5] = new TSquareMatrix2D(SqrStruct[0][1]); // M12
    SqMat[0][1][6] = new TSquareMatrix2D(SqrStruct[0][1]); // M21
    SqMat[0][1][7] = new TSquareMatrix2D(SqrStruct[0][1]); // M22

/* velocity phase_2*/
    SqMat[0][2][0] = new TSquareMatrix2D(SqrStruct[0][2]); // A11
    SqMat[0][2][1] = new TSquareMatrix2D(SqrStruct[0][2]); // A12
    SqMat[0][2][2] = new TSquareMatrix2D(SqrStruct[0][2]); // A21
    SqMat[0][2][3] = new TSquareMatrix2D(SqrStruct[0][2]); // A22
    SqMat[0][2][4] = new TSquareMatrix2D(SqrStruct[0][2]); // M11
    SqMat[0][2][5] = new TSquareMatrix2D(SqrStruct[0][2]); // M12
    SqMat[0][2][6] = new TSquareMatrix2D(SqrStruct[0][2]); // M21
    SqMat[0][2][7] = new TSquareMatrix2D(SqrStruct[0][2]); // M22

    Mat[0][0]  = new TMatrix2D(Struct[0][0]); // B1
    Mat[0][1]  = new TMatrix2D(Struct[0][0]); // B2
    Mat[0][2]  = new TMatrix2D(Struct[1][0]); // B1T
    Mat[0][3]  = new TMatrix2D(Struct[1][0]); // B2T

 /* phase 1 */
    Mat[1][0]  = new TMatrix2D(Struct[0][1]); // B1
    Mat[1][1]  = new TMatrix2D(Struct[0][1]); // B2
    Mat[1][2]  = new TMatrix2D(Struct[1][1]); // B1T
    Mat[1][3]  = new TMatrix2D(Struct[1][1]); // B2T

 /* phase 2 */
    Mat[2][0]  = new TMatrix2D(Struct[0][2]); // B1
    Mat[2][1]  = new TMatrix2D(Struct[0][2]); // B2
    Mat[2][2]  = new TMatrix2D(Struct[1][2]); // B1T
    Mat[2][3]  = new TMatrix2D(Struct[1][2]); // B2T

 /* grid */
  SqMat[2][0][0] = new TSquareMatrix2D(SqrStruct[2][0]); // G11
  SqMat[2][0][1] = new TSquareMatrix2D(SqrStruct[2][0]); // G12
  SqMat[2][0][2] = new TSquareMatrix2D(SqrStruct[2][0]); // G21
  SqMat[2][0][3] = new TSquareMatrix2D(SqrStruct[2][0]); // G22

 /* grid phase_1* */
  SqMat[2][1][0] = new TSquareMatrix2D(SqrStruct[2][1]); // G11
  SqMat[2][1][1] = new TSquareMatrix2D(SqrStruct[2][1]); // G12
  SqMat[2][1][2] = new TSquareMatrix2D(SqrStruct[2][1]); // G21
  SqMat[2][1][3] = new TSquareMatrix2D(SqrStruct[2][1]); // G22

 /* grid phase_2* */
  SqMat[2][2][0] = new TSquareMatrix2D(SqrStruct[2][2]); // G11
  SqMat[2][2][1] = new TSquareMatrix2D(SqrStruct[2][2]); // G12
  SqMat[2][2][2] = new TSquareMatrix2D(SqrStruct[2][2]); // G21
  SqMat[2][2][3] = new TSquareMatrix2D(SqrStruct[2][2]); // G22

/* surfactant */
  SqMat[3][0][0] = new TSquareMatrix2D(SqrStruct[3][0]); // C_A
  SqMat[3][0][1] = new TSquareMatrix2D(SqrStruct[3][0]); // C_M

  SqMat[3][1][0] = new TSquareMatrix2D(SqrStruct[3][1]); // C_A
  SqMat[3][1][1] = new TSquareMatrix2D(SqrStruct[3][1]); // C_M

  SqMat[3][2][0] = new TSquareMatrix2D(SqrStruct[3][2]); // C_A
  SqMat[3][2][1] = new TSquareMatrix2D(SqrStruct[3][2]); // C_M

 /* interface surfactant matrices */
   SqMat_IFace[0] = new TSquareMatrix1D(IFaceStruct[0]); // Gamma_A
   SqMat_IFace[1] = new TSquareMatrix1D(IFaceStruct[0]); // Gamma_M


   N_U = N_DOFs[0][0];
   N_P = N_DOFs[1][0];
   N_Unknowns = 2*N_U+N_P;

   OutPut("dof velocity : "<< setw(10) << 2*N_U << endl);
   OutPut("dof pressure : "<< setw(10) << N_P << endl);
   OutPut("dof all      : "<<  setw(10) << N_Unknowns  << endl);

   OutPut("dof surfactant (outer phase) : "<< setw(5) << N_DOFs[3][2] << endl);

// *************************************************************************
//   memory allocation for solution and rhs
// *************************************************************************

 /* velocity */
   usol = new double[N_Unknowns];
   Sol[0][1] = new double[2*N_DOFs[0][1]+N_DOFs[1][1]];
   Sol[0][2] = new double[2*N_DOFs[0][2]+N_DOFs[1][2]];

   Rhs[0][0] = new double[N_Unknowns];
   Rhs[0][1] = new double[2*N_DOFs[0][1]+N_DOFs[1][1]];
   Rhs[0][2] = new double[2*N_DOFs[0][2]+N_DOFs[1][2]];

 /* pressure */
   psol = usol+2*N_U;
   Sol[1][1] = Sol[0][1]+2*N_DOFs[0][1];
   Sol[1][2] = Sol[0][2]+2*N_DOFs[0][2];
   Rhs[1][0] = Rhs[0][0]+2*N_U;
   Rhs[1][1] = Rhs[0][1]+2*N_DOFs[0][1];
   Rhs[1][2] = Rhs[0][2]+2*N_DOFs[0][2];

   memset(usol, 0, N_Unknowns*SizeOfDouble);
   memset(Sol[0][1], 0, (2*N_DOFs[0][1]+N_DOFs[1][1])*SizeOfDouble);
   memset(Sol[0][2], 0, (2*N_DOFs[0][2]+N_DOFs[1][2])*SizeOfDouble);
   memset(Rhs[0][0], 0, N_Unknowns*SizeOfDouble);
   memset(Rhs[0][1], 0, (2*N_DOFs[0][1]+N_DOFs[1][1])*SizeOfDouble);
   memset(Rhs[0][2], 0, (2*N_DOFs[0][2]+N_DOFs[1][2])*SizeOfDouble);


 /* grid velocity */
   Sol[2][0] = new double[2*N_DOFs[2][0]];
   Sol[2][1] = new double[2*N_DOFs[2][1]];
   Sol[2][2] = new double[2*N_DOFs[2][2]];
   Rhs[2][0] = new double[2*N_DOFs[2][0]];
   Rhs[2][1] = new double[2*N_DOFs[2][1]];
   Rhs[2][2] = new double[2*N_DOFs[2][2]];

   memset(Sol[2][0], 0, 2*N_DOFs[2][0]*SizeOfDouble);
   memset(Sol[2][1], 0, 2*N_DOFs[2][1]*SizeOfDouble);
   memset(Sol[2][2], 0, 2*N_DOFs[2][2]*SizeOfDouble);
   memset(Rhs[2][0], 0, 2*N_DOFs[2][0]*SizeOfDouble);
   memset(Rhs[2][1], 0, 2*N_DOFs[2][1]*SizeOfDouble);
   memset(Rhs[2][2], 0, 2*N_DOFs[2][2]*SizeOfDouble);

/* surfactant */
   csol = new double[N_DOFs[3][0]];
   Sol[3][1] = new double[N_DOFs[3][1]];
   Sol[3][2] = new double[N_DOFs[3][2]];
   Rhs[3][0] = new double[N_DOFs[3][0]];
   Rhs[3][1] = new double[N_DOFs[3][1]];
   Rhs[3][2] = new double[N_DOFs[3][2]];

   memset(csol, 0, N_DOFs[3][0]*SizeOfDouble);
   memset(Sol[3][1], 0, N_DOFs[3][1]*SizeOfDouble);
   memset(Sol[3][2], 0, N_DOFs[3][2]*SizeOfDouble);
   memset(Rhs[3][0], 0, N_DOFs[3][0]*SizeOfDouble);
   memset(Rhs[3][1], 0, N_DOFs[3][1]*SizeOfDouble);
   memset(Rhs[3][2], 0, N_DOFs[3][2]*SizeOfDouble);

 /* interface surfactant */
   Sol[N_FESpaces_All][0] = new double[N_DOFs[N_FESpaces_All][0]];
   Rhs[N_FESpaces_All][0] = new double[N_DOFs[N_FESpaces_All][0]];

    memset(Sol[N_FESpaces_All][0], 0, N_DOFs[N_FESpaces_All][0]*SizeOfDouble);
    memset(Rhs[N_FESpaces_All][0], 0, N_DOFs[N_FESpaces_All][0]*SizeOfDouble);

   Surfact_sol = new double[N_DOFs[4][0]];
   memset(Surfact_sol, 0, N_DOFs[4][0]*SizeOfDouble);

  /* velocity and pressure */
   New_u = new TFEVectFunct2D(velocity_space, UString,  UString,  usol, N_U, 2);
   New_p = new TFEFunction2D(pressure_space, PString,  PString, usol+2*N_U, N_P);

   New_u->Interpolate(VeloVect[0][0]);
   
   
   //no need to interpolate Pressure for direct solver
//    New_p->Interpolate(FeFunct[2][0]);

//copy freeboundary velo values
   ValuesU2 = usol + N_U;
   m=0;
   if(Old_N_IsoEdges==N_MovVert[6])
    {
     // no additional vertices are added during remeshing on the interface, so copy the old isovertces
     m=0;
     comp=0;
     for(j=0;j<Old_N_IsoEdges;j++)
      {
       Me = Coll_Cells[0][j];
       l = N_List[3][j];
       k = Me->GetN_Edges();
       
       Cell_No = Me->GetGlobalCellNo();
       FeId = velocity_space->GetFE2D(Cell_No, Me);  
       FeDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(FeId);
       JointDOF = FeDesc->GetJointDOF(l); 
       N_JointDOF = FeDesc->GetN_JointDOF();
       DOF = GlobalNumbers[0][0] + BeginIndex[0][0][Cell_No];       

       Me->GetVertex(l)->GetCoords(TX[0], TY[0]);
       Me->GetVertex((l+1) % k)->GetCoords(TX[1], TY[1]);
       
       temp0 = ((Intpol_Coord[2*m] - TX[0])*(Intpol_Coord[2*m] - TX[0]) + (Intpol_Coord[2*m+1] - TY[0])*(Intpol_Coord[2*m+1] - TY[0]) );
       comp = m+ORDER;
       temp1 = ((Intpol_Coord[2*comp] - TX[1])*(Intpol_Coord[2*comp] - TX[1]) + (Intpol_Coord[2*comp+1] - TY[1])*(Intpol_Coord[2*comp+1] - TY[1]) );
       
       if(temp0<1.e-8 && temp1<1.e-8 ) // update the isovertex from the old value
        {         
         for(jj=0;jj<N_JointDOF-1;jj++)
         {
          kk = DOF[JointDOF[jj]];

          usol[kk] = Intpol_VeloValues[2*m];
          ValuesU2[kk] = Intpol_VeloValues[2*m+1];
          m++;
         }
        }
       else
        {
         // only second order conforming elements implimented
         cout<< " No match in Old_N_IsoEdges  and N_MovVert[6] in Remesh "<<endl;
         
         cout <<j <<  "X: " <<  Intpol_Coord[2*m] << " " << FreePts[0][j] <<  " " << TX[0] <<" Y: " <<  Intpol_Coord[2*m+1] << " " << FreePts[1][j] <<" " << TY[0] << endl;

         cout << j << "X: " <<  Intpol_Coord[2*comp] << " " << FreePts[0][j+1] <<  " " << TX[1] <<" Y: " <<  Intpol_Coord[2*comp+1] << " " << FreePts[1][j+1] <<" " << TY[1] << endl;

         exit(0);
        }
      } // for(j=0;j<Old_N_IsoEdges-1;j++) 
    }//if(Old_N_IsoEdges==N_MovVert[6])
   else
    {
     cout<< " New vertices are added on the interface in Remesh "<<endl;
     exit(0);
    }

   delete FeSpaces[0][0];
   delete FeSpaces[1][0];
   FeSpaces[0][0] = velocity_space;
   FeSpaces[1][0] = pressure_space;

   delete [] Sol[0][0];
   Sol[0][0] = usol;

//    delete VeloVect[0][0];
   VeloVect[0][0] = New_u;
   FeFunct[0][0] = VeloVect[0][0]->GetComponent(0);
   FeFunct[1][0] = VeloVect[0][0]->GetComponent(1);

   delete FeFunct[2][0];
   FeFunct[2][0] = New_p;

  VeloVect[0][1] = new TFEVectFunct2D(FeSpaces[0][1], UString,  UString,  Sol[0][1], N_DOFs[0][1], 2);
  FeFunct[0][1] = VeloVect[0][1]->GetComponent(0);
  FeFunct[1][1] = VeloVect[0][1]->GetComponent(1);
  FeFunct[2][1] = new TFEFunction2D(FeSpaces[1][1], PString,  PString,
                                    Sol[0][1]+2*N_DOFs[0][1], N_DOFs[1][1]);


  VeloVect[0][2] = new TFEVectFunct2D(FeSpaces[0][2], UString,  UString,  Sol[0][2], N_DOFs[0][2], 2);
  FeFunct[0][2] = VeloVect[0][2]->GetComponent(0);
  FeFunct[1][2] = VeloVect[0][2]->GetComponent(1);
  FeFunct[2][2] = new TFEFunction2D(FeSpaces[1][2], PString,  PString,
                                    Sol[0][2]+2*N_DOFs[0][2], N_DOFs[1][2]);


 /* grid velocity */
   VeloVect[1][0] = new TFEVectFunct2D(FeSpaces[2][0], gridString,  gridString,  Sol[2][0], N_DOFs[2][0], 2);
   FeFunct[3][0] = VeloVect[1][0]->GetComponent(0);
   FeFunct[4][0] = VeloVect[1][0]->GetComponent(1);

   VeloVect[1][1] = new TFEVectFunct2D(FeSpaces[2][1], gridString,  gridString,  Sol[2][1], N_DOFs[2][1], 2);
   FeFunct[3][1] = VeloVect[1][1]->GetComponent(0);
   FeFunct[4][1] = VeloVect[1][1]->GetComponent(1);

   VeloVect[1][2] = new TFEVectFunct2D(FeSpaces[2][2], gridString,  gridString,  Sol[2][2], N_DOFs[2][2], 2);
   FeFunct[3][2] = VeloVect[1][2]->GetComponent(0);
   FeFunct[4][2] = VeloVect[1][2]->GetComponent(1);

/* Surfactant */
   New_c = new TFEFunction2D(Surfact_space, SurfactString,  SurfactString, csol, N_DOFs[3][0]);
   FeFunct[5][1] = new TFEFunction2D(FeSpaces[3][1], SurfactString,  SurfactString, Sol[3][1], N_DOFs[3][1]);
   FeFunct[5][2] = new TFEFunction2D(FeSpaces[3][2], SurfactString,  SurfactString, Sol[3][2], N_DOFs[3][2]);

   New_c->Interpolate(FeFunct[5][0]);

   delete FeSpaces[3][0];
   delete [] Sol[3][0];
   delete FeFunct[5][0];

   FeSpaces[3][0] = Surfact_space;
   Sol[3][0] = csol;
   FeFunct[5][0] = New_c;

   length = BeginIndex[3][0][1] -  BeginIndex[3][0][0];
   for(i=0;i<N_Cells_P2;i++)
    {
     DOF = GlobalNumbers[3][0] + BeginIndex[3][0][GlobalCell_Index[1][i]];
     DOF_P2 = GlobalNumbers[3][2] + BeginIndex[3][2][i];

     for(jj=0;jj<length;jj++)
      {
       k = DOF[jj];
       l1 = DOF_P2[jj];
       Sol[3][2][l1] = Sol[3][0][k];
      }
     }
   memset(Sol[3][0], 0, N_DOFs[3][0]*SizeOfDouble);

 /* interface surfactant */
   delete IFaceFeFunct[0];
   IFaceFeFunct[0] = new TFEFunction1D(IFaceFeSpaces[0], IFaceSString, IFaceSString,
                                       Sol[N_FESpaces_All][0], N_DOFs[N_FESpaces_All][0]);

   Surfact_FeFunct  = new TFEFunction2D(Surfact_Outputspace, IFaceSString,  IFaceSString,
                                        Surfact_sol, N_DOFs[4][0]);

//    Surfact_FeFunct->Intpol_2Phase(FeFunct[6][0]);
   // copy the values from the previous interface (no interpolation)
   
   GlobalNumbers_S = Surfact_Outputspace->GetGlobalNumbers();
   BeginIndex_S = Surfact_Outputspace->GetBeginIndex();
   
   if(Old_N_IsoEdges==N_MovVert[6])
    {
     // no additional vertices are added during remeshing on the interface, so the No. DOF in old and new interface are same
     m=0;
     
     for(j=0;j<Old_N_IsoEdges;j++)
      {
       Me = Coll_Cells[0][j];
       N = N_List[2][j]; //local cell number
       FeId = Surfact_Outputspace->GetFE2D(N, Me);
       FeDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(FeId);
       
       l = N_List[3][j]; //local edge number
       JointDOF = FeDesc->GetJointDOF(l);
       N_JointDOF = FeDesc->GetN_JointDOF();
       DOF = GlobalNumbers_S + BeginIndex_S[N];

       // dof at the starting vertex
       Surfact_sol[DOF[JointDOF[0]]] = Intpol_Values[m]; 
       m++;

       Joint = Me->GetJoint(l);
       IntFaceJoint = (TIsoInterfaceJoint *)Joint;
       k = IntFaceJoint->GetN_Vertices();
       IsoVertices = IntFaceJoint->GetVertices();
    
       if(k!=N_JointDOF-2)
        {
         OutPut("No match in isopoints per interface in  Remesh" <<endl;) 
         exit(0);
        }
     
       for(i3=0;i3<k;i3++)
        {
         Surfact_sol[DOF[JointDOF[i3+1]]] = Intpol_Values[m];
         m++;               
        }//  for(i3=0;i3<k;i3++     
     
      } // for(j=0;j<Old_N_IsoEdges;j++)
      
     //end vertex dof
     Surfact_sol[DOF[JointDOF[N_JointDOF-1]]] = Intpol_Values[m];
     m++;        
    } //  if(Old_N_IsoEdges==N_MovVert[6])
   else
    {
     cout<< " New vertices are added on the interface in Remesh "<<endl;
     exit(0);
    }

   delete FeSpaces[4][0];
   delete [] Sol[4][0];
   delete FeFunct[6][0];

   FeSpaces[4][0] = Surfact_Outputspace;
   Sol[4][0] = Surfact_sol;
   FeFunct[6][0] = Surfact_FeFunct;

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

  // remove all existing cells and joints 
  // not vertices, since all vertices are vertices of 2D domain
  for(i=0;i<Old_S_N_RootCells;i++)
     delete (TGridCell*)Surf_OldCellTree[i];

  OutPut(Old_S_N_RootCells<<" surface cells were deleted"<<endl);
   delete [] Surf_OldCellTree;

    delete SOld_Coll;

  delete [] mm;
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

    N_Outer = Me->GetLocalCellNo();
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
      cout<< " N_JointDOF != N_BaseFunct_low " <<endl;
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


      uref_Outer = TFEDatabase2D::GetJointValues2D(BaseFuncts[FEId_Outer],
                     LineQuadFormula, IJoint);
      uxiref_Outer = TFEDatabase2D::GetJointDerivatives2D(BaseFuncts[FEId_Outer],
                        LineQuadFormula, IJoint, D10);
      uetaref_Outer = TFEDatabase2D::GetJointDerivatives2D(BaseFuncts[FEId_Outer],
                        LineQuadFormula, IJoint, D01);

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
              ((TQuadIsoparametric *)F_K)->GetOrigValues(IJoint,  zeta[k],
                        N_BaseFunct_Outer, uref_Outer[k], uxiref_Outer[k], uetaref_Outer[k],
                        uorig_Outer, uxorig_Outer, uyorig_Outer);
            break;

            case BFUnitTriangle:
              ((TTriaIsoparametric *)F_K)->GetOrigValues(IJoint,  zeta[k],
                        N_BaseFunct, uref[k], uxiref[k], uetaref[k],
                        uorig, uxorig, uyorig);
              ((TTriaIsoparametric *)F_K)->GetOrigValues(IJoint,  zeta[k],
                        N_BaseFunct_Outer, uref_Outer[k], uxiref_Outer[k], uetaref_Outer[k],
                        uorig_Outer, uxorig_Outer, uyorig_Outer);


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
         for(l=0;l<N_BaseFunct_Outer;l++)
           C += C_Outer[DOF_Outer[l]]*uorig_Outer[l];

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

    N_Outer = Me->GetLocalCellNo();
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
      cout<< " N_JointDOF != N_BaseFunct_low " <<endl;
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


      uref_Outer = TFEDatabase2D::GetJointValues2D(BaseFuncts[FEId_Outer],
                     LineQuadFormula, IJoint);
      uxiref_Outer = TFEDatabase2D::GetJointDerivatives2D(BaseFuncts[FEId_Outer],
                        LineQuadFormula, IJoint, D10);
      uetaref_Outer = TFEDatabase2D::GetJointDerivatives2D(BaseFuncts[FEId_Outer],
                        LineQuadFormula, IJoint, D01);

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
              ((TQuadIsoparametric *)F_K)->GetOrigValues(IJoint,  zeta[k],
                        N_BaseFunct_Outer, uref_Outer[k], uxiref_Outer[k], uetaref_Outer[k],
                        uorig_Outer, uxorig_Outer, uyorig_Outer);
            break;

            case BFUnitTriangle:
              ((TTriaIsoparametric *)F_K)->GetOrigValues(IJoint,  zeta[k],
                        N_BaseFunct, uref[k], uxiref[k], uetaref[k],
                        uorig, uxorig, uyorig);
              ((TTriaIsoparametric *)F_K)->GetOrigValues(IJoint,  zeta[k],
                        N_BaseFunct_Outer, uref_Outer[k], uxiref_Outer[k], uetaref_Outer[k],
                        uorig_Outer, uxorig_Outer, uyorig_Outer);
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
         for(l=0;l<N_BaseFunct_Outer;l++)
           C += C_Outer[DOF_Outer[l]]*uorig_Outer[l];

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
//           cout <<" TangDivU:  " << TangDivU <<endl;

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

 // ====================================================================
// modify matrices due to integrals on free surface
// axialsymmetric case
// ====================================================================
void FreeSurf_2PhaseSurfAxial3D(TSquareMatrix2D *A11, TSquareMatrix2D *A22,
                     double *rhs1, double *rhs2,
                     BoundCondFunct2D *BoundaryCondition,
                     double dt, int PhaseNo, TFEFunction2D *Surfact)
{
  int i, j, k, l, DOF_R, DOF_L, m;
  TBaseCell *cell;
  TCollection *Coll;
  int N_Cells, N_Vertices, N_Edges, Semi_implicit=0;
  TJoint *joint;
  TIsoBoundEdge *isoboundedge;
  TBoundComp *BoundComp;
  int comp, N_U, test_L=1, test_R=1, ORDER;
  double t0, t1, n0, n1, normn, line_wgt;
  BoundCond Cond0, Cond1;
  int JointNumbers[MAXN_JOINTS], IJoint, N_IsoJoints;
  FE2D FEId, TFEId;
  TFE2D *ele;
  RefTrans2D RefTrans;
  TRefTrans2D *F_K;
  TFESpace2D *fespace, *surfactantspace;
  BF2DRefElements RefElement;
  QuadFormula1D LineQuadFormula;
  TQuadFormula1D *qf1;
  int N_LinePoints;
  double *LineWeights, *zeta;
  double x0, y0, x1, y1,tx,ty,mod_t, x, y;
  int N_BaseFunct, *N_BaseFuncts, TN_BaseFunct, *TJointDOF, TN_DOF_Local, *TDOF;
  double **uref, **uxiref, **uetaref;
  double **Turef, **Tuxiref, **Tuetaref;
  double uorig[MaxN_BaseFunctions2D], uxorig[MaxN_BaseFunctions2D];
  double uyorig[MaxN_BaseFunctions2D];
  double Tuorig[MaxN_BaseFunctions2D], Tuxorig[MaxN_BaseFunctions2D];
  double Tuyorig[MaxN_BaseFunctions2D];
  BaseFunct2D *BaseFuncts;
  double r2, r, T_val[3], T_Marangoni, *S_Values;
  int *KCol, *RowPtr, *JointDOF, N_DOF, N_Surf;
  double *ValuesA11, *ValuesA12, *ValuesA21, *ValuesA22;
  int *BeginIndex, *GlobalNumbers, *DOF, TestDOF, AnsatzDOF;
  int index1, index2, Phase_No;
  double val, theta, factor1, factor2, angle, Gamma;
  int count=0, count1=0, count2=0;
  double  X_B[100], Y_B[100], r_axial, d1, d2, e1, e2, ngrad_test, ngrad_ansatz;
  double ngrad_Gamma, GammaE1, GammaE2;
  int *TGlobalNumbers, *TBeginIndex, local_dof;
  double Gamma_Max, Gamma_Infty;

  TFEDesc2D *FeDesc, *TFeDesc;

  BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
  N_BaseFuncts = TFEDatabase2D::GetN_BaseFunctFromFE2D();

  fespace = A11->GetFESpace();
  Coll = fespace->GetCollection();
  N_Cells = Coll->GetN_Cells();

  BeginIndex = fespace->GetBeginIndex();
  GlobalNumbers = fespace->GetGlobalNumbers();

  RowPtr = A11->GetRowPtr();
  KCol = A11->GetKCol();

  ValuesA11 = A11->GetEntries();
  ValuesA22 = A22->GetEntries();

  double We = TDatabase::ParamDB->WB_NR;

  surfactantspace = Surfact->GetFESpace2D();
  S_Values=Surfact->GetValues();
  TGlobalNumbers = surfactantspace->GetGlobalNumbers();
  TBeginIndex = surfactantspace->GetBeginIndex();
  
#ifdef __WITHSURFACTANT__  
#ifndef __MARANGONISTRESSTEST__  
  Gamma_Max = 0;
  N_Surf = Surfact->GetLength();

  for(i=0;i<N_Surf;i++)
   if(Gamma_Max<S_Values[i]) Gamma_Max=S_Values[i];

   OutPut("Gamma_Max " << Gamma_Max<<endl);
#endif
#endif

// surfactant elasticity E
  double E = TDatabase::ParamDB->REACTOR_P10;
//Equation of state, 0 linear, 1 non-linear
  int EOS = int(TDatabase::ParamDB->REACTOR_P11);
//\Gamma_1/Gamma_\infty
  double D = TDatabase::ParamDB->REACTOR_P12 / TDatabase::ParamDB->REACTOR_P16 ;
  double CHAR_L = TDatabase::ParamDB->CHAR_L0;
//  if(EOS==1 && Gamma_Max>1.e-12) // non-linear equation of state
//   {
// // Assumed that REACTOR_P13 = 1
// // Gamma_Infty is 10% more than Gamma_Max
//    Gamma_Infty = 1.1*Gamma_Max;
// 
//    E = E *Gamma_Infty;
//    Dscal(N_Surf, 1./Gamma_Infty, SurfactValues);
//   }

 for(i=0;i<N_Cells;i++)
  {
//      cout << endl << "CELL number: " << i << endl;
   cell = Coll->GetCell(i);
   Phase_No = cell->GetPhase_ID();
   if(Phase_No==0)
    {
    N_Edges = cell->GetN_Edges();
    IJoint = 0;
    for(j=0;j<N_Edges;j++)
     {
      joint = cell->GetJoint(j);
      if(joint->GetType() == IsoInterfaceJoint || joint->GetType() == InterfaceJoint)
      {
       FEId = fespace->GetFE2D(i, cell);
       DOF = GlobalNumbers + BeginIndex[i];
       N_BaseFunct = N_BaseFuncts[FEId];
       ele = TFEDatabase2D::GetFE2D(FEId);
       RefElement = TFEDatabase2D::GetRefElementFromFE2D(FEId);
       ORDER = TFEDatabase2D::GetAccuracyFromFE2D(FEId);
       l = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
//        l = 3;
       LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(2*l);
       qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
       qf1->GetFormulaData(N_LinePoints, LineWeights, zeta);
       TFEDatabase2D::GetBaseFunct2DFromFE2D(FEId)->MakeRefElementData(LineQuadFormula);
       
#ifdef __WITHSURFACTANT__   
#ifndef __MARANGONISTRESSTEST__  
       TFEId = surfactantspace->GetFE2D(i, cell);
       TFEDatabase2D::GetBaseFunct2DFromFE2D(TFEId)->MakeRefElementData(LineQuadFormula);
       TFeDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(TFEId);
       TN_BaseFunct = N_BaseFuncts[TFEId];
       TJointDOF = TFeDesc->GetJointDOF(j);
       TN_DOF_Local = TFeDesc->GetN_JointDOF();
       TDOF = TGlobalNumbers + TBeginIndex[i];
#endif
#endif      
       switch(RefElement)
       {
        case BFUnitSquare:
          RefTrans = QuadIsoparametric;
          F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
          ((TQuadIsoparametric *)F_K)->SetCell(cell);
          ((TQuadIsoparametric *)F_K)->GetOrigBoundFromRef(j, N_LinePoints, zeta, X_B, Y_B);
        break;

        case BFUnitTriangle:
          RefTrans = TriaIsoparametric;
          F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
          ((TTriaIsoparametric *)F_K)->SetApproximationOrder(ORDER);
          ((TTriaIsoparametric *)F_K)->SetCell(cell);
          ((TTriaIsoparametric *)F_K)->GetOrigBoundFromRef(j, N_LinePoints, zeta, X_B, Y_B);
        break;
       } // endswitch


     uref = TFEDatabase2D::GetJointValues2D(BaseFuncts[FEId],
                        LineQuadFormula, j);
     uxiref = TFEDatabase2D::GetJointDerivatives2D(BaseFuncts[FEId],
                        LineQuadFormula, j, D10);
     uetaref = TFEDatabase2D::GetJointDerivatives2D(BaseFuncts[FEId],
                        LineQuadFormula, j, D01);
#ifdef __WITHSURFACTANT__  
#ifndef __MARANGONISTRESSTEST__  
     Turef = TFEDatabase2D::GetJointValues2D(BaseFuncts[TFEId],
                        LineQuadFormula, j);
     Tuxiref = TFEDatabase2D::GetJointDerivatives2D(BaseFuncts[TFEId],
                        LineQuadFormula, j, D10);
     Tuetaref = TFEDatabase2D::GetJointDerivatives2D(BaseFuncts[TFEId],
                        LineQuadFormula, j, D01);
#endif
#endif
       for(k=0;k<N_LinePoints;k++)
        {
         switch(RefElement)
          {
            case BFUnitSquare:
              ((TQuadIsoparametric *)F_K)->GetOrigValues(j,  zeta[k],
                        N_BaseFunct, uref[k], uxiref[k], uetaref[k],
                        uorig, uxorig, uyorig);
#ifdef __WITHSURFACTANT__  
#ifndef __MARANGONISTRESSTEST__  
              ((TQuadIsoparametric *)F_K)->GetOrigValues(j,  zeta[k],
                        TN_BaseFunct, Turef[k], Tuxiref[k], Tuetaref[k],
                        Tuorig, Tuxorig, Tuyorig);
            break;
#endif
#endif
            case BFUnitTriangle:

              ((TTriaIsoparametric *)F_K)->GetOrigValues(j,  zeta[k],
                        N_BaseFunct, uref[k], uxiref[k], uetaref[k],
                        uorig, uxorig, uyorig);
#ifdef __WITHSURFACTANT__ 
#ifndef __MARANGONISTRESSTEST__  
              ((TTriaIsoparametric *)F_K)->GetOrigValues(j,  zeta[k],
                        TN_BaseFunct, Turef[k], Tuxiref[k], Tuetaref[k],
                        Tuorig, Tuxorig, Tuyorig);
            break;
#endif 
#endif 
          } // endswitch

          // modify matrices
         F_K->GetTangent(j, zeta[k], t0, t1);  // old line
         r_axial = fabs(X_B[k]);   // r value in the axial symmetric integral
         if(X_B[k]<=0)
          {
           cell->GetVertex(j)->GetCoords(x0, y0);       
           cell->GetVertex((j+1) % N_Edges)->GetCoords(x1, y1);    
           cout <<"X_B[k] negative in freesurface int, change Quad rule " <<  X_B[k] <<endl;
           cout << j << " " << k <<" X  " << x0<<" Y  " << y0<<" X  " << x1<<" Y  " << y1 <<endl;
           cout << "If the quad formula is correct, then the time step might be large !!!! " <<endl; 
           exit(0);
          }

         normn = sqrt(t0*t0+t1*t1);
	 n0 =  t1/normn;
         n1 = -t0/normn;

       T_val[0] = 0.; T_val[1] = 0.; T_val[2] = 0.;
       
#ifdef __WITHSURFACTANT__     

#ifdef __MARANGONISTRESSTEST__  

      T_val[0] = Y_B[k];  // see Trygvasion L=15a, a=0.5
      T_val[1] = 0.;   
      T_val[2] = 1.;
#else
       for(l=0;l<TN_DOF_Local;l++)
        {
          // assumed that the velo space and 2D surfactant space are same fe space
          local_dof   = TJointDOF[l];
          m = TDOF[local_dof];
          val = S_Values[m];


          T_val[0] += val*Tuorig[local_dof];  // Surfactant C
          T_val[1] += val*Tuxorig[local_dof];  // C_x
          T_val[2] += val*Tuyorig[local_dof];  // C_y
        } // for(l=0;l<TN_

        if(T_val[0]<0. )
	 {
          OutPut(i<< "x : "<<X_B[k]<< " y: " << Y_B[k] <<"  Surfactant exceeds the reference value, T_val= " <<T_val[0]<<endl);
          //numerical correction
          T_val[0]=0.; 
         }         
#endif         
          if(EOS==0)
          {
           Gamma =(1. + E*(D - T_val[0]) ); 
           T_Marangoni = normn*E*(t0*T_val[1] + t1*T_val[2])/We;           
           }
          else
          {
           Gamma =(1. + E*log(1. - T_val[0]) );
   
          if(T_val[0]!=1)
           {T_Marangoni = normn*E*( t0*T_val[1] + t1*T_val[2]) / ( We*(1. - T_val[0]) );}
          else
           {T_Marangoni = 0.; }
     
           //see SolubleSurf JCP paper
           if(Gamma<0.1)
            {
             Gamma = 0.1;
             T_Marangoni = 0.;
            }
          }                 
#else
    D = 0.;
    Gamma = 1;
    T_Marangoni = 0.;    
#endif           
         // cout<< " Gamma " << Gamma<<" T_Marangoni " << T_Marangoni<<endl;  
          // Multiply with time step dt in the main program not here
          r = normn/We;           
          for(l=0;l<N_BaseFunct;l++)
          {
           TestDOF = DOF[l];
           
           // updating rhs
            ngrad_test= n0*uxorig[l] + n1*uyorig[l];
            d1 = uxorig[l] - ngrad_test*n0;
            d2 = uyorig[l] - ngrad_test*n1;
            
	    
	    //surface tension is inside the Dive, so no additional terms due tointegration by parts
#ifdef __WITHSURFACTANT__   
            d1 *= Gamma;
            d2 *= Gamma;
// 
//             ngrad_Gamma = n0*T_val[1] + n1*T_val[2];
//             GammaE1 = (T_val[1] - ngrad_Gamma*n0)*uorig[l];
//             GammaE2 = (T_val[2] - ngrad_Gamma*n1)*uorig[l];
//             
//             if(EOS==0)
//              {
//               GammaE1 *=E;
//               GammaE2 *=E;
//              }
//             else
//              {
//               GammaE1 *=(E/(1.- T_val[0]));
//               GammaE2 *=(E/(1.- T_val[0]));
//              }
// #else
//               GammaE1 = 0.;  
//               GammaE2 = 0.;  
#endif
             
// rhs1
//             val = r_axial*( (1.-n0*n0)*(d1 - GammaE1) - n0*n1*(d2 -GammaE2) );
            val = r_axial*( (1.-n0*n0)*d1 - n0*n1*d2);    
            val +=(Gamma*uorig[l]); // due to axialsymmetric
            val *= LineWeights[k]*r;
            rhs1[TestDOF] -= val;
//     Marangoni convection
#ifdef __WITHSURFACTANT__  
            val = r_axial*t0*uorig[l];
            val *= LineWeights[k]*T_Marangoni;
            rhs1[TestDOF] -= val;
#endif
// rhs2
//             val =  r_axial*( -n1*n0*(d1 - GammaE1) + (1.-n1*n1)*(d2 - GammaE2) );
            val =  r_axial*( -n1*n0*d1 + (1.-n1*n1)*d2 );    
            val *= LineWeights[k]*r;
            rhs2[TestDOF] -= val;
//     Marangoni convection
#ifdef __WITHSURFACTANT__  
            val = r_axial*t1*uorig[l];
            val *= LineWeights[k]*T_Marangoni;
            rhs2[TestDOF] -= val;
#endif
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

//               val =(d1 - GammaE1)*e1 + (d2 -GammaE2)*e2 
              val =d1*e1 + d2*e2 + Gamma*(uorig[l]*uorig[m]/(r_axial*r_axial));
              val *= dt*LineWeights[k]*r*r_axial;
              // cout << "A11: " << TestDOF << " ";
              // cout << AnsatzDOF << " " << val << endl;
              ValuesA11[index1] += val;

//               val =(d1 - GammaE1)*e1 + (d2 -GammaE2)*e2;
              val =d1*e1 + d2*e2;
              val *= dt*LineWeights[k]*r*r_axial;

              // cout << "A22: " << TestDOF << " ";
              // cout << AnsatzDOF << " " << val << endl;
              ValuesA22[index1] += val;

            } // endfor m
          } // endfor l
        } // endfor k
       } // 
     } // endfor j
     }    // if(Phase_No==0)
   } // endfor i
//  delete [] SurfactValues;
 }

void Surfact2D_InterfaceInt(TSquareMatrix2D *A, double *rhs, BoundCondFunct2D SurfactBoundCondition, 
                           TFEFunction2D *FeFunct, TFEFunction1D *SurfaceFeFunct, 
                           TFESpace2D *GlobalFE_Space, int *Cell_array, int *Joint_array, double *XmaxVal)
{
  int i, j, k, l, N_Cells, Phase_ID, N_Vertices, N_Edges, ORDER;
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
//   surf_flux =  0.;

 for(i=0;i<N_Cells_Surface;i++)
  {
    Me_Surface = Coll_Surface->GetCell(i);
    FEId_Surface = fespace_surface->GetFE1D(i, Me_Surface);
    Element = TFEDatabase2D::GetFE1D(FEId_Surface);
    N_BaseFunct_Surface = Element->GetN_DOF();
    DOF_Surface = GlobalNumbers_Surface + BeginIndex_Surface[i];


    Global_N = Cell_array[i];     // outerphase
    Me = GlobalColl->GetCell(Global_N);
    N = Me->GetLocalCellNo();
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

//         // update rhs for all test functions
     for(l=0;l<N_BaseFunct;l++)
       {
         if((index1 = DOF[l])<ActiveBound)
          {
           rhs[index1] += rhsval*uorig[l];
//            surf_flux += rhsval*uorig[l];
//            cout <<"l " << l << " rhs "<<  val*uorig[l] <<endl;
          }
       }//  for(l=0;l<N_BaseFunct;l++)
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
  int i, j, k, l, N_Cells, Phase_ID, N_Vertices, N_Edges, ORDER;
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
  Coll_Surface = fespace_surface->GetCollection();
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
    N = Me->GetLocalCellNo();
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

int main(int argc, char* argv[])
{
  TDomain *Domain = new TDomain();
  TDomain *IFaceDomain = new TDomain();
  TDatabase *Database = new TDatabase();
  TFEDatabase2D *FEDatabase = new TFEDatabase2D();
  TCollection *coll, *mortarcoll = NULL, *Coll_P1, *Coll_P2;
  TCollection  *Coll_Multi[2];
  TBaseCell *Me, *cell;
  TBaseCell **Coll_Cells[3];
  TFESpace2D *velocity_space,  *velocity_space_P1, *velocity_space_P2;
  TFESpace2D *pressure_space, *pressure_space_P1, *pressure_space_P2;
  TFESpace2D *Grid_space, *Grid_space_P1, *Grid_space_P2;
  TFESpace2D *Surfact_space, *Surfact_space_P1, *Surfact_space_P2;
  TFESpace2D ***FeSpaces;

// variables for surfactant
  FE1D *FE1D_List;
  TCollection *IFace_Coll;
  TFESpace1D *IFaceSurfact_space;
  TFESpace1D **IFaceFeSpaces = new TFESpace1D*[2];

  BoundCondFunct2D *BoundaryConditions[2];
  BoundCondFunct2D *SurfactBoundaryConditions[1];
  BoundValueFunct2D *BoundValues[2], *BoundValuesAuxProblem[3];
  BoundValueFunct2D *SurfactBoundValues[1];
  BoundCondFunct2D *GridBoundaryConditions[1];
  BoundValueFunct2D *GridBoundValues[1];


  TDiscreteForm2D *DiscreteFormGalerkin;
  TDiscreteForm2D *DiscreteFormNLGalerkin;
  TDiscreteForm2D *DiscreteFormRHS;
  TDiscreteForm2D *DiscreteFormGrid, *DiscreteFormSurfact_OutPhase, *DiscreteFormSurfact_SUPG;

  TDiscreteForm2D *DiscreteForm;

  TSquareStructure2D ***SqrStruct;
  TStructure2D ***Struct = new TStructure2D**[3];
  TSquareStructure1D **IFaceStruct = new TSquareStructure1D*[1];

  TSquareMatrix2D ****SqMat, *SQMATRICES_GRID[4];
  TMatrix2D ***Mat;
  TSquareMatrix1D **SqMat_IFace,*SQMATRICES_IFace[2] ;

  TFEVectFunct2D ***VeloVect;
  TFEVectFunct2D *RefGridPos, *AuxGridPos, *GridPos;
  TFEVectFunct2D *RefGridPos_P1, *AuxGridPos_P1, *GridPos_P1;
  TFEVectFunct2D *RefGridPos_P2, *AuxGridPos_P2, *GridPos_P2;
  TFEFunction2D ***FeFunct, *fefct[7], *Surfact_FeFunct, *OutputFeunction;
  TFEFunction1D **IFaceFeFunct = new TFEFunction1D*[1];

  TOutput2D *Output;

  TSquareMatrix2D *SQMATRICES_SURFACT[2], *SQMATRICES[8];
  TMatrix2D *MATRICES[4];
  TSquareMatrix **sqmatrices = (TSquareMatrix **)SQMATRICES;
  TMatrix **matrices = (TMatrix **)MATRICES;
  TFESpace2D *fesp[3], *ferhs[3];
  TFESpace1D *IFacefesp[1], *IFaceferhs[1];
  TAuxParam2D *aux;
  DefectProc *Defect;
  TJoint *Joint; 
  TIsoBoundEdge *FreeBound_Joint;
  TIsoInterfaceJoint *IntfaceBound_Joint;
  TBoundComp *BoundComp;
  TBoundPart *BoundPart[2];
  TBdLine *UpdateBound[6];
  TBdCircle *UpdateIntface;
  TBaseCell **CellTree;
  TVertex **InterFace_Vert, **CollVertex, **SlipBound_Vert, **VertexDel, **NewVertices;
  std::ostringstream os, opts;
  TIsoInterfaceJoint *IntFaceJoint;
  TBoundEdge  ***Slip_Joint = new TBoundEdge**[7], *Solid_Joint, *tempSlip_Joint;
  TVertex ***MovBoundVert = new TVertex**[9], *temp_Mov, **IsoVertices;

  struct triangulateio In, Out;

  int i, j, jj, k, l, m, m1, m2, ID;
  int N_Cells, N_Cells_P1, N_Cells_P2, N_IFaceCells;
  int ret, N_Interface_Vert, N_Old_Face_Vert;
  int sp_No, sp_no0, sp_no1, N_FESpaces_All, N_FESpaces, N_FESpaces_low, ORDER;
  int **GlobalCell_Index = new int*[2], methods, N_SubSteps,N_Active, *RowPtr;
  int pressure_space_code, velocity_space_code, N_Rhs, N_LinIterCurr;
  int **N_DOFs, **Bound_DOFs, N_Unknowns, N_U, N_P, img=1, ReParam_img=1, N_RectMatrices, Max_It;
  int time_discs, very_first_time=0, N_CActive, N_IActive, N_SquareMatrices, N_LinIter;
  int ***GlobalNumbers, ***BeginIndex, length, *DOF, *DOF_P2, *DOF_P1, l1, VSP;
  int N_Hori, N_Verti, N_vert1, N_vert2, N_vert3, beg_index, end_index,bdpart;
  int N_Boundary_Vert, N_Interf_Vertices, CurrComp;
  int N_RootCells, N_E, N_BoundaryNodes;
  int *PointNeighb, maxEpV = 0, a, b, Neighb_tmp, Neib[2];
  int CurrNeib, len1, len2, comp, Currregion, temp_segment;
  int In_Index, CurrVertex, N_CollVertex, CurrJoint, CurrInterfCell, bct, bct1;
  int *PartMarker, *Triangles;
  int N_Joints, N_Vertices, N_SlipJoints, N_FreeJoints, N_G;
  int *N_MovVert = new int[8], N_Remesh=0, iso_update, *mm = new int[8];
  int *GridKCol, *GridRowPtr,  *GridKCol_P1, *GridRowPtr_P1,  *GridKCol_P2, *GridRowPtr_P2;
  int **N_List = new int*[4], N_BData=1, Max_It_scalar;
  int surf_couple_var, N_AllIntVertices;


  double RE_NR, total_time, *I_FaceX, *I_FaceY, t3, *CRHSs[1], *SRHSs[1], GammaXmaxVal[2];
  double x_short, y_short,  temp0, temp1, Pos_Indicator, gamma, Tgamma;
  double x_mid, y_mid, x, y, hmin, temp, T_a, T_b, C_x, C_y, solver_time_curr;
  double tx, ty, nx, ny, ***Sol, ***Rhs, t, end_time, oldtau, tau;
  double **gridsol, **gridrhs;
  double *C_defect, *Csol_old, *oldsol, *B, *oldrhs, *RHSs[3], *defect, *C_B, *CRhs_old;
  double *I_defect, *Isol_old, *I_B, *IRhs_old, *Csol_nonlinearstep;
  double *S_defect, *Ssol_old, TX[2], TY[2], area;
  double residual, impuls_residual, oldresidual, solver_time, limit, t1, t2, t4;
  double residual_scalar, oldresidual_scalar;
  double x0, y0, hi, phi1, phi2, theta, h_interface;
  double T, *Coordinates, Remesh_Time=0;
  double *Angle = new double[2], phi, r, deviation;
  double left, right, top, bottom, rad1, Rx, Ry, Lx, Ly, x1, x2, y1, y2, x3, y3, x4, y4;
  double *Entries[4], *Entries_P1[4], *Entries_P2[4];
  double *refpos, *auxpos, *pos, *tmp, limit_scalar;
  double *refpos_P1, *auxpos_P1, *pos_P1, *tmp_P1, *refpos_P2, *auxpos_P2, *pos_P2;
  double Surf_Mass[2], Params[10], Initial_IFaceSurfactMass, Initial_SurfactMass;
  double h, c_x[50], c_xVal[50], Xvalues[5], **FreePts, InitVolume, CurrVolume;
  double fh, fhtot=0, fhmin = 1e2, fhmax= 0, fhlimit = 5e-3;
  double sphericity, radius, *outputsol, *Intpol_Coord,  *Intpol_Values, *Intpol_VeloValues;
  
  FreePts = new double *[2];

  boolean DirichletBC, remeshed= FALSE, reparam=FALSE;
  boolean FluidFlow=TRUE, MovingMesh=TRUE;

  if(!FluidFlow) MovingMesh=FALSE;


  char *PRM, *GEO;
  char *PsBaseName, *GrapeBaseName, *GnuBaseName, *ReadGrapeBaseName;
  char *VtkBaseName;
  // strings
  char ReadinDat[] = "readin.dat";
  char NameString[] = "name";
  char UString[] = "u";
  char PString[] = "p";
  char PsiString[] = "psi";
  char TString[] = "T";
  char SurfactString[] = "C";
  char gridString[] = "W";
  char refposString[] = "refpos";
  char auxposString[] = "auxpos";
  char posString[] = "pos";
  char IFaceSString[] = "I_C";
  const char vtkdir[] = "VTK";
  const char BDdir[] = "BDData";
  
  mtrace();
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
  Database->WriteParamDB(argv[0]);
  Database->WriteTimeDB();
  ExampleFile();

  mkdir(vtkdir, 0777);
  mkdir(BDdir, 0777);  
 
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

  limit = TDatabase::ParamDB->SC_NONLIN_RES_NORM_MIN_SADDLE;
  Max_It = TDatabase::ParamDB->SC_NONLIN_MAXIT_SADDLE;
  surf_couple_var = int(TDatabase::ParamDB->REACTOR_P7);

  limit_scalar = TDatabase::ParamDB->REACTOR_P8;
  Max_It_scalar = int(TDatabase::ParamDB->REACTOR_P9);

  if(surf_couple_var==1) Max_It_scalar = 1; // explicit coupling

  VSP = TDatabase::ParamDB->VELOCITY_SPACE;

  if (abs(VSP) > 20)
   {ORDER = abs(VSP) - 20;}
  else if ( abs(VSP) > 10)
    {ORDER = abs(VSP) - 10;}
  else ORDER = abs(VSP);

/* No multigrid */
  if (TDatabase::ParamDB->SC_MG_TYPE_SADDLE)
    TDatabase::ParamDB->SC_MG_TYPE_SADDLE = 0;

/* velocity and pressure, mesh, surfactant */
  N_FESpaces_All = 5;

  FeSpaces = new TFESpace2D**[N_FESpaces_All];
  N_DOFs = new int*[N_FESpaces_All+1]; // interface surfactant
  Bound_DOFs = new int*[N_FESpaces_All];
  SqrStruct = new TSquareStructure2D**[N_FESpaces_All];
  SqMat = new TSquareMatrix2D***[N_FESpaces_All];
  SqMat_IFace = new TSquareMatrix1D*[2];

  Mat = new TMatrix2D**[3];
  Sol = new double**[N_FESpaces_All+1]; // interface surfactant
  Rhs = new double**[N_FESpaces_All+1]; // interface surfactant
  VeloVect = new TFEVectFunct2D**[2];
  FeFunct = new TFEFunction2D**[7];
  GlobalNumbers = new int**[N_FESpaces_All];
  BeginIndex = new int**[N_FESpaces_All];
  gridsol = new double*[3];
  gridrhs = new double*[3];

 /* Two-Phase flows */
  for(i=0;i<N_FESpaces_All;i++)
   {
    FeSpaces[i] = new TFESpace2D*[3];
    N_DOFs[i] = new int[3];
    Bound_DOFs[i] = new int[3];
    SqrStruct[i] = new TSquareStructure2D*[3];
    SqMat[i] = new TSquareMatrix2D**[3];

    Sol[i] = new double*[3];
    Rhs[i] = new double*[3];

    GlobalNumbers[i] = new int*[3];
    BeginIndex[i] = new int*[3];
   }
    N_DOFs[N_FESpaces_All] = new int[1]; // interface surfactant
    Sol[N_FESpaces_All] = new double*[1]; // interface surfactant
    Rhs[N_FESpaces_All] = new double*[1]; // interface surfactant

   Struct[0] = new TStructure2D*[3]; // BT
   Struct[1] = new TStructure2D*[3]; // B

   SqMat[0][0] = new TSquareMatrix2D*[10]; // Aij and Mij, F11, F22
   SqMat[0][1] = new TSquareMatrix2D*[10];
   SqMat[0][2] = new TSquareMatrix2D*[10];

   SqMat[2][0] = new TSquareMatrix2D*[4]; // grid Gij
   SqMat[2][1] = new TSquareMatrix2D*[4];
   SqMat[2][2] = new TSquareMatrix2D*[4];

   SqMat[3][0] = new TSquareMatrix2D*[2]; // Surfactant Tij
   SqMat[3][1] = new TSquareMatrix2D*[2];
   SqMat[3][2] = new TSquareMatrix2D*[2];

   Mat[0] = new TMatrix2D*[4];// BT, B
   Mat[1] = new TMatrix2D*[4];// BT, B Phase_1
   Mat[2] = new TMatrix2D*[4];// BT, B Phase_2

   VeloVect[0] = new TFEVectFunct2D*[3];
   VeloVect[1] = new TFEVectFunct2D*[3];

  for(i=0;i<7;i++)
   FeFunct[i] = new TFEFunction2D*[3];

  InitializeDiscreteForms_2PhaseAxial3D( DiscreteFormGalerkin,
                                  DiscreteFormNLGalerkin,
                                  DiscreteFormRHS,
                                  DiscreteFormGrid, LinCoeffs, GridCoeffs,
                                  TDatabase::ParamDB->NSTYPE);

  InitializeDiscreteForms_Moving(DiscreteFormSurfact_OutPhase, DiscreteFormSurfact_SUPG, SurfactCoeffs);
 
  {  
  //======================================================================
// read boundary parameterization and generate the mesh
//======================================================================

  Domain->Init(PRM, GEO);
// exit(0);
//======================================================================
// Triangular for grid generation begin
//======================================================================
  boolean AllowEdgeRef = (boolean) TDatabase::ParamDB->MESHGEN_ALLOW_EDGE_REF;
  BoundPart[0] = Domain->GetBdPart(0);
  BoundPart[1] = Domain->GetBdPart(0);

  UpdateBound[0]  = (TBdLine*)BoundPart[1]->GetBdComp(0);
  UpdateBound[1]  = (TBdLine*)BoundPart[1]->GetBdComp(1);
  UpdateBound[2]  = (TBdLine*)BoundPart[1]->GetBdComp(2);
  UpdateBound[3]  = (TBdLine*)BoundPart[1]->GetBdComp(3);
  UpdateBound[4]  = (TBdLine*)BoundPart[1]->GetBdComp(4);
  UpdateBound[5]  = (TBdLine*)BoundPart[1]->GetBdComp(5);
  UpdateIntface = (TBdCircle*)BoundPart[0]->GetBdComp(6);

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

  area = TDatabase::ParamDB->Area;

  opts<<'p'; // Constrained Delaunay Triangulation:
           // initial values - only points defined on the boundary of the domain;
           // triangulation near boundary may variate from Delaunay criterion
  opts<<"q"<<  TDatabase::ParamDB->MESHGEN_REF_QUALITY;
              // Quality mesh generation with no angles smaller than 20 degrees;

  opts<<"a"<< area; // Imposes a maximum triangle area.
  opts<<'e'; // Outputs a list of edges of the triangulation
  opts<<'z'; // Numbers if items starting from 0
//   opts<<"VVVV"; // Gives detailed information about what Triangle is doing
//   opts<<'Q'; // Supress all explanation of what Triangle is doing, unless an error occurs
//   opts<<'Y'; // Supress adding vertices on boundary edges
  opts<<'j'; //Jettisons(discard) vertices that are not part of the final triangulation from the output .node file (including duplicate input vertices and vertices ``eaten'' by holes).
  opts << "nA";
  opts<<ends;

  N_Interface_Vert = int (TDatabase::ParamDB->P6);    //interface points
#ifdef __MASSTRANSTEST__
  double Xi[6] = {0., 0., 0., 5., 5., 0.};
  double Yi[6] = {5., 3., 0., 0., 8., 8.};  
  T_a = 1.;
  T_b = 1.; 
  C_x = 0.0;  // center of the inner phase
  C_y = 4., // center of the inner phase 
  FluidFlow=FALSE, MovingMesh=FALSE;
#else
  double Xi[6] = {0., 0., 0., 1., 1., 0.};
  double Yi[6] = {0.75, .25, 0., 0., 6., 6.};
  T_a = .25;
  T_b = .25;
  C_x = 0.0;  // center of the inner phase
  C_y = 0.5, // center of the inner phase
#endif
  
  deviation=TDatabase::ParamDB->P3;
  int mode = int (TDatabase::ParamDB->P4);

// prolate
//  T_a = 0.7925; // x axis value in ellipse
//  T_b = 1.592; // y axis value in ellipse !!! set in modifyGausscoord funct also
//    T_a = 1.8; // x axis value in ellipse
//    T_b = 4.; // y axis value in ellipse !!! set in modifyGausscoord funct also

// oblate
//  T_a = 1.592; // x axis value in ellipse
//  T_b = 0.7925; // y axis value in ellipse !!! set in modifyGausscoord funct also


  phi1 = -Pi/2.0; // initial degree value of interface
  phi2 = Pi/2.0; // end degree value of interface

  x = Xi[0];
  y = Yi[0];

  phi = atan2(y-C_y, x-C_x);
    if(mode==2)
     {
      r = 1.0 + deviation*(sqrt(5./Pi)*(3.*cos(phi+Pi/2.)*cos(phi+Pi/2.) -1. )/4.);
      T_a = r;
      T_b = r;
     }
    else if(mode==4) 
     {
      temp = cos(phi+Pi/2.);
      r = 1.0 + deviation*(sqrt(1./Pi)*(35.*temp*temp*temp*temp - 30.*temp*temp +3. )*3./16.);
      T_a = r;
      T_b = r;
     }

//    Xi[0] = C_x + T_a*cos(phi);
   Yi[0] = C_y + T_b*sin(phi);

  x = Xi[1];
  y = Yi[1];

  phi = atan2(y-C_y, x-C_x);
    if(mode==2)
     {
      r = 1.0 + deviation*(sqrt(5./Pi)*(3.*cos(phi+Pi/2.)*cos(phi+Pi/2.) -1. )/4.);
      T_a = r;
      T_b = r;
     }
    else if(mode==4) 
     {
      temp = cos(phi+Pi/2.);
      r = 1.0 + deviation*(sqrt(1./Pi)*(35.*temp*temp*temp*temp - 30.*temp*temp +3. )*3./16.);
      T_a = r;
      T_b = r;
     }

   Xi[1] = C_x + T_a*cos(phi);
   Yi[1] = C_y + T_b*sin(phi);

      if(fabs(Xi[1])<1e-10) Xi[1] = 0;

   N_Hori  = 10;      // number of horoyontal vertices
   N_Verti = 25;

   hi = (Yi[4]-Yi[3])/N_Verti;

   N_vert1 = 3*abs(int((Yi[1]-Yi[0])/hi));
    if(N_vert1<4 )N_vert1 =4;
   N_vert2 = abs(int((Yi[2]-Yi[1])/hi));
    if(N_vert2<4 )N_vert2 =4;
   N_vert3 = abs(int((Yi[0]-Yi[5])/hi));
    if(N_vert3<4 )N_vert3 =4;

  N_Verti = N_vert1+N_vert2+N_vert3;       // number of horoyontal vertices

  N_Boundary_Vert = 2*N_Hori + 2*N_Verti;

  N_Interf_Vertices = N_Interface_Vert+N_Boundary_Vert;
  In.numberofpoints = N_Interf_Vertices;
//   cout<< "N_Hori " <<N_Hori << " N_Verti " << N_Verti<< " In.numberofpoints " <<In.numberofpoints << endl;
  In.pointlist = new double[2*In.numberofpoints];
  In.pointmarkerlist = new int[In.numberofpoints];
  In.numberofpointattributes = 0;

//   start vertex already given but segment has to be given
  In.numberofsegments = In.numberofpoints;
  In.segmentlist = new int[2*In.numberofsegments];
  In.segmentmarkerlist = new int[In.numberofsegments];
  In.numberofholes = 0;
  In.holelist = NULL;
  In.numberofregions = 0;
  In.regionlist = NULL;

  In_Index = 0;
  CurrComp = 1;
  Currregion = 0;

  hi = (Yi[1] - Yi[0])/N_vert1;
  x0 = Xi[0];
  y0 = Yi[0];
  y  = y0;

  end_index = 0;
  if(x0<1.e-8) x0 = 0.;
  // points and segments on the horizontal boundary (marker=1)
  for(i=0;i<N_vert1;i++) // without last point
   {
    In.pointlist[2*In_Index] = x0;
    In.pointlist[2*In_Index+1] = y;
// cout<<" x : "<< In.pointlist[2*In_Index] << " y : "<< In.pointlist[2*In_Index+1] <<endl;
    In.pointmarkerlist[In_Index] = CurrComp;
    In.segmentlist[2*In_Index] = In_Index;
    In.segmentlist[2*In_Index+1] = In_Index+1;
    In.segmentmarkerlist[In_Index] = CurrComp;
    In_Index++;
    y = y0 + double(i+1)*hi;
   }
  CurrComp++;
//   cout<<endl;

  hi = (Yi[2] - Yi[1])/N_vert2;
  x0 = Xi[1];
  y0 = Yi[1];
  y  = y0;

  beg_index=In_Index;
  // points and segments on the horizontal boundary (marker=1)
  for(i=0;i<N_vert2;i++) // without last point
   {
    In.pointlist[2*In_Index] = x0;
    In.pointlist[2*In_Index+1] = y;
// cout<<" x : "<< In.pointlist[2*In_Index] << " y : "<< In.pointlist[2*In_Index+1] <<endl;
    In.pointmarkerlist[In_Index] = CurrComp;
    In.segmentlist[2*In_Index] = In_Index;
    In.segmentlist[2*In_Index+1] = In_Index+1;
    In.segmentmarkerlist[In_Index] = CurrComp;
    In_Index++;
    y = y0 + double(i+1)*hi;
   }
  CurrComp++;
//   cout<<endl;

  hi = (Xi[3] - Xi[2])/N_Hori;
  x0 = Xi[2];
  y0 = Yi[2];
  x  = Xi[2];
  // points and segments on the horizontal boundary (marker=1)
  for(i=0;i<N_Hori;i++) // without last point
   {
    In.pointlist[2*In_Index] = x;
    In.pointlist[2*In_Index+1] = y0;
// cout<<" x : "<< In.pointlist[2*In_Index] << " y : "<< In.pointlist[2*In_Index+1] <<endl;
    In.pointmarkerlist[In_Index] = CurrComp;
    In.segmentlist[2*In_Index] = In_Index;
    In.segmentlist[2*In_Index+1] = In_Index+1;
    In.segmentmarkerlist[In_Index] = CurrComp;
    In_Index++;
    x = x0 + double(i+1)*hi;
   }
  CurrComp++;

//   cout<<endl;
  hi = (Yi[4] - Yi[3])/N_Verti;
  x0 = Xi[3];
  y0 = Yi[3];
  y  = Yi[3];
  // points and segments on the horizontal boundary (marker=1)
  for(i=0;i<N_Verti;i++) // without last point
   {
    In.pointlist[2*In_Index] = x0;
    In.pointlist[2*In_Index+1] = y;
// cout<<" x : "<< In.pointlist[2*In_Index] << " y : "<< In.pointlist[2*In_Index+1] <<endl;
    In.pointmarkerlist[In_Index] = CurrComp;
    In.segmentlist[2*In_Index] = In_Index;
    In.segmentlist[2*In_Index+1] = In_Index+1;
    In.segmentmarkerlist[In_Index] = CurrComp;
    In_Index++;
    y = y0 + double(i+1)*hi;
   }
  CurrComp++;
  cout<<endl;
  hi = (Xi[5] - Xi[4])/N_Hori;
  x0 = Xi[4];
  y0 = Yi[4];
  x  = Xi[4];
  // points and segments on the horizontal boundary (marker=1)
  for(i=0;i<N_Hori;i++) // without last point
   {
    In.pointlist[2*In_Index] = x;
    In.pointlist[2*In_Index+1] = y0;
// cout<<" x : "<< In.pointlist[2*In_Index] << " y : "<< In.pointlist[2*In_Index+1] <<endl;
    In.pointmarkerlist[In_Index] = CurrComp;
    In.segmentlist[2*In_Index] = In_Index;
    In.segmentlist[2*In_Index+1] = In_Index+1;
    In.segmentmarkerlist[In_Index] = CurrComp;
    In_Index++;
    x = x0 + double(i+1)*hi;
   }
  CurrComp++;

  cout<<endl;
  hi = (Yi[0] - Yi[5])/N_vert3;
  x0 = Xi[5];
  y0 = Yi[5];
  y  = Yi[5];
  // points and segments on the horizontal boundary (marker=1)
 for(i=0;i<N_vert3;i++) // without last point
   {
    In.pointlist[2*In_Index] = x0;
    In.pointlist[2*In_Index+1] = y;
// cout<<" x : "<< In.pointlist[2*In_Index] << " y : "<< In.pointlist[2*In_Index+1] <<endl;
    In.pointmarkerlist[In_Index] = CurrComp;
    In.segmentlist[2*In_Index] = In_Index;
    In.segmentlist[2*In_Index+1] = In_Index+1;
    In.segmentmarkerlist[In_Index] = CurrComp;
    In_Index++;
    y = y0 + double(i+1)*hi;
   }
  CurrComp++;

  In.segmentlist[2*(In_Index-1)+1] = 0;

//   temp_segment=In_Index;
// exit(0);
//   cout<<endl;

  t = (phi2-phi1)/(N_Interface_Vert);

  N_Old_Face_Vert = N_Interface_Vert;
  I_FaceX = new double[N_Interface_Vert]; // for spline construction
  I_FaceY = new double[N_Interface_Vert]; // for spline construction

  I_FaceX[0] = In.pointlist[2*beg_index];
  I_FaceY[0] = In.pointlist[2*beg_index+1];

  In.segmentlist[2*In_Index] = beg_index;
  In.segmentlist[2*In_Index+1] = In_Index;
  In.segmentmarkerlist[In_Index] = CurrComp;

//   cout<< " x : "<< I_FaceX[0] << " y : "<< I_FaceY[0] <<endl;

 // points and segments on the interface (marker=2)
  theta = phi1;
  double t0 = theta;

  theta = t0 + double(0+1)*t;
  for(i=0;i<N_Interface_Vert-1;i++)
    {
//      cout<<" theta : "<< theta <<endl;
      x = C_x + T_a*cos(theta);
      y = C_y + T_b*sin(theta);

// spherical harmonic of order 2
// spherical harmonic of order 4
    phi = atan2(y-C_y, x-C_x);
    if(mode==2)
     {
      r = 1.0 + deviation*(sqrt(5./Pi)*(3.*cos(phi+Pi/2.)*cos(phi+Pi/2.) -1. )/4.);
      T_a = r;
      T_b = r;
     }
    else if(mode==4) 
     {
      temp = cos(phi+Pi/2.);
      r = 1.0 + deviation*(sqrt(1./Pi)*(35.*temp*temp*temp*temp - 30.*temp*temp +3. )*3./16.);
      T_a = r;
      T_b = r;
     }

      x = C_x + T_a*cos(phi);
      y = C_y + T_b*sin(phi);

      if(fabs(x)<1e-10) x = 0;
      
      
      I_FaceX[i+1] = x;
      I_FaceY[i+1] = y;

      In.pointlist[2*In_Index] = I_FaceX[i+1];
      In.pointlist[2*In_Index+1] = I_FaceY[i+1];
      In.pointmarkerlist[In_Index] = CurrComp;

//       cout<< In_Index<<" x : "<< I_FaceX[i+1] << " y : "<< I_FaceY[i+1] <<endl;
      if(i==1)
       {
        h_interface = sqrt((I_FaceX[i-1]-I_FaceX[i])*(I_FaceX[i-1]-I_FaceX[i]) +
                            (I_FaceY[i-1]-I_FaceY[i])*(I_FaceY[i-1]-I_FaceY[i]) );
      cout << "h_interface " <<h_interface << endl;
       }
//   cout<<(180./Pi)*theta<< " x : "<< In.pointlist[2*In_Index] << " y : "<<In.pointlist[2*In_Index+1] <<endl;

      In.segmentlist[2*(In_Index+1)] = In_Index;
      In.segmentlist[2*(In_Index+1)+1] = In_Index+1;
      In.segmentmarkerlist[(In_Index+1)] = CurrComp;
      In_Index++;

      theta = t0 + double(i+2)*t;
    }

  In.pointlist[2*In_Index] = -100.;
  In.pointlist[2*(In_Index)+1] = -100.;
  In.pointmarkerlist[In_Index] = CurrComp;


  In.segmentlist[2*(In_Index)+1] = end_index;
  cout<< In_Index<<endl;


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
    OutPut(i<<' '<<In.segmentlist[2*i]<< ' '<<' '<<In.segmentlist[2*i+1]<<' '<<
	   In.pointlist[2*i]<<' '<<In.pointlist[2*i+1]<<endl);
cout<<endl;
  exit(0);
*/
  triangulate((char*)opts.str().c_str(), &In, &Out, (struct triangulateio *)NULL);

/*  for(i=0;i<Out.numberofpoints;i++)
     OutPut(i<<' '<<Out.pointmarkerlist[i]<<' '<<
      Out.pointlist[2*i]<<' '<<Out.pointlist[2*i+1]<<endl);
  exit(0)*/
// cout<< " test triangulate " << endl;
// exit(0);
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

  delete [] VertexDel;
  OutPut(CurrVertex<<" vertices were deleted"<<endl);

 // remove all existing cells and joints
  for(i=0;i<N_RootCells;i++)
    delete (TGridCell*)CellTree[i];
  OutPut(N_RootCells<<" cells were deleted"<<endl);
   delete [] CellTree;
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

//   OutPut("left: "<<left<<" right: "<<right<<" top: "<<top<<" bottom: "<<bottom<<endl);

  Domain->SetBoundBox(right-left,top-bottom);
  Domain->SetBoundBoxstart(left,bottom);

 // Solid Bound startx, starty, x length and y length
  UpdateBound[0]->SetParams(Xi[0], Yi[0], Xi[1]-Xi[0],Yi[1]-Yi[0]);
  UpdateBound[1]->SetParams(Xi[1], Yi[1], Xi[2]-Xi[1],Yi[2]-Yi[1]);
  UpdateBound[2]->SetParams(Xi[2], Yi[2], Xi[3]-Xi[2],Yi[3]-Yi[2]);
  UpdateBound[3]->SetParams(Xi[3], Yi[3], Xi[4]-Xi[3],Yi[4]-Yi[3]);
  UpdateBound[4]->SetParams(Xi[4], Yi[4], Xi[5]-Xi[4],Yi[5]-Yi[4]);
  UpdateBound[5]->SetParams(Xi[5], Yi[5], Xi[0]-Xi[5],Yi[0]-Yi[5]);
// Free boundary xmid, ymid, radius_a, radius_b, start angle, end angle
  UpdateIntface->SetParams(C_x, C_y, T_a, T_b, phi1, phi2);

//   cout << " Xi[5] : " << Xi[5]<< endl;
 // cout << " In.pointlist[N_Boundary_Vert] : " << In.pointlist[0]<< endl;
 // generate cells

  CellTree = new TBaseCell*[N_RootCells];

  for (i=0;i<N_RootCells;i++)
  {
    CellTree[i] = new TMacroCell(TDatabase::RefDescDB[Triangle], 0);

    CellTree[i]->SetVertex(0, NewVertices[Out.trianglelist[3*i    ]]);
    CellTree[i]->SetVertex(1, NewVertices[Out.trianglelist[3*i + 1]]);
    CellTree[i]->SetVertex(2, NewVertices[Out.trianglelist[3*i + 2]]);

//       if(N_Regions)
//       ((TMacroCell *) CellTree[i])->SetSubGridID(
//           (int)(outt.triangleattributelist[i]+0.05));
//      else
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
// cout<< " test main N_G "<< N_G <<endl;
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
     bdpart = 0;
     CurrComp = Out.edgemarkerlist[i] - 1;
//         cout << " CurrComp "<< CurrComp<<endl;
     if (CurrComp >= 100)
      {
       CurrComp -= 100;
       bdpart = 0;
      }


      if(Domain->GetBdPart(bdpart)->GetBdComp(CurrComp)->GetTofXY(
            NewVertices[a]->GetX(), NewVertices[a]->GetY(), T_a) ||
          Domain->GetBdPart(bdpart)->GetBdComp(CurrComp)->GetTofXY(
            NewVertices[b]->GetX(), NewVertices[b]->GetY(), T_b))
       {
          cerr<<"Error: could not set parameter values"<<endl;
          OutPut(NewVertices[a]<<endl);
          OutPut(NewVertices[b]<<endl);
        //  exit(0);
       }

      if (CurrNeib == 2)    // 2 cells contain the current edge
        if(Domain->GetBdPart(bdpart)->GetBdComp(CurrComp)->IsFreeBoundary())
          Joint = new TIsoInterfaceJoint(Domain->GetBdPart(bdpart)->GetBdComp(CurrComp),
                  T_a, T_b, CellTree[Neib[0]], CellTree[Neib[1]]);
        else
          Joint = new TInterfaceJoint(Domain->GetBdPart(bdpart)->GetBdComp(CurrComp),
                  T_a, T_b, CellTree[Neib[0]], CellTree[Neib[1]]);
      else
        if(Domain->GetBdPart(bdpart)->GetBdComp(CurrComp)->IsFreeBoundary())
          Joint = new TIsoBoundEdge(Domain->GetBdPart(bdpart)->GetBdComp(CurrComp), T_a, T_b);
        else
          Joint = new TBoundEdge(Domain->GetBdPart(bdpart)->GetBdComp(CurrComp), T_a, T_b);
      

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

  if ( (Joint->GetType()) == InterfaceJoint ||
        Joint->GetType() == IsoInterfaceJoint )
      ((TInterfaceJoint *) Joint)->CheckOrientation();
  }


//   delete [] NewVertices;
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
  }
   for(i=0;i<TDatabase::ParamDB->REFINEMENT; i++)
     Domain->RegRefineAll();


    os.seekp(std::ios::beg);
    os << "Domain" << ".ps" << ends;
    Domain->PS(os.str().c_str(),It_Finest,0);


  BoundaryConditions[0] = BoundCondition;
  BoundaryConditions[1] = BoundCondition;

  SurfactBoundaryConditions[0] = SurfactBoundCondition;

  BoundValues[0] = U1BoundValue;
  BoundValues[1] = U2BoundValue;

  GridBoundaryConditions[0] = GridBoundCondition;
  GridBoundValues[0] = GridBoundValue;

  SurfactBoundValues[0] = SurfactBoundValue;

//   BoundaryConditionsAuxProblem[0] = BoundConditionAuxProblem;
//   BoundaryConditionsAuxProblem[1] = BoundConditionAuxProblem;
//   BoundaryConditionsAuxProblem[2] = BoundConditionAuxProblem;

//   BoundValuesAuxProblem[0] = BoundValueAuxProblem;
//   BoundValuesAuxProblem[1] = BoundValueAuxProblem;
//   BoundValuesAuxProblem[2] = BoundValueAuxProblem;

  TDatabase::TimeDB->CURRENTTIME = TDatabase::TimeDB->STARTTIME;

  t3 = GetTime();
  total_time = t3 - total_time;

// // **********************************************************************
//  split the doamin in to interior and exterior -- begin
// *************************************************************************

  coll=Domain->GetCollection(It_Finest, 0);
  N_Cells = coll->GetN_Cells();
  N_Cells_P1 = 0;
  N_Cells_P2 = 0;

  for(i=0;i<N_Cells;i++)
   {
    cell = coll->GetCell(i);
    cell->SetGlobalCellNo(i);         
    
    x_mid = 0.;  y_mid = 0.;

    for(j=0;j<3;j++)
     {
      cell->GetVertex(j)->GetCoords(x, y);
      x_mid +=x;
      y_mid +=y;
     } // for j

      x_mid /=3.;
      y_mid /=3.;
      hmin = 1000.0;
// find point P on interface for the i'th cell mid point X
// such that ||P - X|| is minimum
      for(k=0;k<N_Old_Face_Vert;k++)
       {
        temp = sqrt( (I_FaceX[k]-x_mid)*(I_FaceX[k]-x_mid) + (I_FaceY[k]-y_mid)*(I_FaceY[k]-y_mid) );
        if(temp<hmin)
         {
          hmin=temp;  T_a = I_FaceX[k];  T_b = I_FaceY[k];
          x_short = T_a;  y_short = T_b;
          sp_No = k;
         }
       } // for k

// find next shortest point 
// i.e other shortest point of a spline containing (T_a,T_b)
// previous point
   if(sp_No==0)
    sp_no0 = N_Old_Face_Vert-1;
   else
    sp_no0 = sp_No-1;

    temp0 = sqrt( (I_FaceX[sp_no0]-x_mid)*(I_FaceX[sp_no0]-x_mid)
             + (I_FaceY[sp_no0]-y_mid)*(I_FaceY[sp_no0]-y_mid) );

// next point 
   if(sp_No==N_Old_Face_Vert-1)
    sp_no1 = 0;
   else
    sp_no1 = sp_No+1;

    temp1 = sqrt( (I_FaceX[sp_no1]-x_mid)*(I_FaceX[sp_no1]-x_mid)
                + (I_FaceY[sp_no1]-y_mid)*(I_FaceY[sp_no1]-y_mid) );

   if( temp0 < temp1)
     {
     T_a = I_FaceX[sp_no0];  T_b = I_FaceY[sp_no0];
     C_x = I_FaceX[sp_No];  C_y = I_FaceY[sp_No];
     }
   else 
     {
     T_a = I_FaceX[sp_No];  T_b = I_FaceY[sp_No];
     C_x = I_FaceX[sp_no1];  C_y = I_FaceY[sp_no1];
     }

     tx = x_mid - x_short; //x distance between the point and the shortest distance point on interface
     ty = y_mid - y_short; //y distance between the point and the shortest distance point on interface

     nx = -(C_y - T_b); // (-) normal at (T_a,T_b) pointing into the inner domain
     ny = -(T_a - C_x); // (-) normal at (T_a,T_b) pointing into the inner domain

     Pos_Indicator = tx*nx + ty*ny;

   if(Pos_Indicator > 0.)
    {
//  cell is in inner domain
//  cout<< "Inner cell " << i << endl;
     cell->SetPhase_ID(0);
     N_Cells_P1++;
    }
   else if(Pos_Indicator < 0.)
    {
//  cell is in outer domain
//  cout<< "outer cell " << i << endl;
     cell->SetPhase_ID(1);
     N_Cells_P2++;
    }
   else
    {
     cout<< " error in identifying phase check remesh2d.c" << endl;
     exit(0);
    }

  }// for i

  OutPut( "Number of inner cells " << N_Cells_P1 << endl);
  OutPut( "Number of outer cells " << N_Cells_P2 << endl);

  GlobalCell_Index[0] = new int[N_Cells_P1];
  GlobalCell_Index[1] = new int[N_Cells_P2];

  Coll_Cells[1] = new TBaseCell *[10*N_Cells_P1];
  Coll_Cells[2] = new TBaseCell *[10*N_Cells_P2];

  m1 = 0;
  m2 = 0;
  for(i=0;i<N_Cells;i++)
   {
    cell = coll->GetCell(i);
    ID = cell->GetPhase_ID();

     if(ID == 0)
      {
//    inner cell
       Coll_Cells[1][m1] = cell;
       cell->SetLocalCellNo(m1);
       GlobalCell_Index[0][m1++] = i;
      }
    else
      {
//    outer cell
       Coll_Cells[2][m2] = cell;
       cell->SetLocalCellNo(m2);
       GlobalCell_Index[1][m2++] = i;
      }
   }

  Coll_Multi[0] = new TCollection(N_Cells_P1, Coll_Cells[1]);
  Coll_Multi[1] = new TCollection(N_Cells_P2, Coll_Cells[2]);

  Coll_P1 = Coll_Multi[0];
  Coll_P2 = Coll_Multi[1];

  delete [] I_FaceX;
  delete [] I_FaceY;
}
// *************************************************************************
//  split the doamin in to interior and exterior -- end
// ************************************************************************
// *************************************************************************
//   FESpaces memory allocation
// *************************************************************************

    coll=Domain->GetCollection(It_Finest, 0);
    cout << endl << endl;

//  list of outer phase cells containing interafce
    N_List[0] = new int[N_Cells];   // Cell_No
    N_List[1] = new int[N_Cells];   //Joint_No
    Domain2DSurf_2Phase(coll, IFaceDomain, N_List);

    IFace_Coll = IFaceDomain->GetCollection(It_Finest, 0);
    N_IFaceCells= IFace_Coll->GetN_Cells();
    cout<< " N_IFaceCells " << N_IFaceCells <<endl;
    FE1D_List = new FE1D[N_Cells];
    for(j=0;j<N_Cells;j++)
     FE1D_List[j] = FE1D(TDatabase::ParamDB->ANSATZ_ORDER);

    TDatabase::IteratorDB[It_EQ]->SetParam(Domain);
    TDatabase::IteratorDB[It_LE]->SetParam(Domain);
    TDatabase::IteratorDB[It_Finest]->SetParam(Domain);
    TDatabase::IteratorDB[It_Between]->SetParam(Domain);
    TDatabase::IteratorDB[It_OCAF]->SetParam(Domain);


/* velocity space and pressure spaces */
    GetVelocityAndPressureSpace(coll,BoundCondition,
                                mortarcoll, velocity_space,
                                pressure_space, &pressure_space_code,
                                TDatabase::ParamDB->VELOCITY_SPACE,
                                TDatabase::ParamDB->PRESSURE_SPACE);

    velocity_space_code =   TDatabase::ParamDB->VELOCITY_SPACE;
    TDatabase::ParamDB->INTERNAL_PRESSURE_SPACE = pressure_space_code;

    GetVelocityAndPressureSpace(Coll_P1,BoundCondition,
                                mortarcoll, velocity_space_P1,
                                pressure_space_P1, &pressure_space_code,
                                TDatabase::ParamDB->VELOCITY_SPACE,
                                TDatabase::ParamDB->PRESSURE_SPACE);
    GetVelocityAndPressureSpace(Coll_P2,BoundCondition,
                                mortarcoll, velocity_space_P2,
                                pressure_space_P2, &pressure_space_code,
                                TDatabase::ParamDB->VELOCITY_SPACE,
                                TDatabase::ParamDB->PRESSURE_SPACE);

   FeSpaces[0][0] = velocity_space;
   FeSpaces[0][1] = velocity_space_P1;
   FeSpaces[0][2] = velocity_space_P2;
   N_DOFs[0][0] =  velocity_space->GetN_DegreesOfFreedom();
   N_DOFs[0][1] =  velocity_space_P1->GetN_DegreesOfFreedom();
   N_DOFs[0][2] =  velocity_space_P2->GetN_DegreesOfFreedom();

   FeSpaces[1][0] = pressure_space;
   FeSpaces[1][1] = pressure_space_P1;
   FeSpaces[1][2] = pressure_space_P2;
   N_DOFs[1][0] =  pressure_space->GetN_DegreesOfFreedom();
   N_DOFs[1][1] =  pressure_space_P1->GetN_DegreesOfFreedom();
   N_DOFs[1][2] =  pressure_space_P2->GetN_DegreesOfFreedom();

   GlobalNumbers[0][0] = FeSpaces[0][0]->GetGlobalNumbers();
   BeginIndex[0][0]    = FeSpaces[0][0]->GetBeginIndex();
   GlobalNumbers[0][1] = FeSpaces[0][1]->GetGlobalNumbers();
   BeginIndex[0][1]    = FeSpaces[0][1]->GetBeginIndex();
   GlobalNumbers[0][2] = FeSpaces[0][2]->GetGlobalNumbers();
   BeginIndex[0][2]    = FeSpaces[0][2]->GetBeginIndex();

//  if(TDatabase::ParamDB->Moving_Mesh)
  {
/* grid spaces */
    Grid_space = new TFESpace2D(coll, NameString, PsiString, GridBoundCondition,
                                1, NULL);
    Grid_space_P1 = new TFESpace2D(Coll_P1, NameString, PsiString, GridBoundCondition,
                                1, NULL);
    Grid_space_P2 = new TFESpace2D(Coll_P2, NameString, PsiString, GridBoundCondition,
                                1, NULL);

   FeSpaces[2][0] = Grid_space;
   FeSpaces[2][1] = Grid_space_P1;
   FeSpaces[2][2] = Grid_space_P2;
   N_DOFs[2][0] =  Grid_space->GetN_DegreesOfFreedom();
   N_DOFs[2][1] =  Grid_space_P1->GetN_DegreesOfFreedom();
   N_DOFs[2][2] =  Grid_space_P2->GetN_DegreesOfFreedom();

   Bound_DOFs[2][0] = N_DOFs[2][0] - FeSpaces[2][0]->GetN_Inner();
   Bound_DOFs[2][1] = N_DOFs[2][1] - FeSpaces[2][1]->GetN_Inner();
   Bound_DOFs[2][2] = N_DOFs[2][2] - FeSpaces[2][2]->GetN_Inner();


   GlobalNumbers[2][0] = FeSpaces[2][0]->GetGlobalNumbers();
   BeginIndex[2][0]    = FeSpaces[2][0]->GetBeginIndex();
   GlobalNumbers[2][1] = FeSpaces[2][1]->GetGlobalNumbers();
   BeginIndex[2][1]    = FeSpaces[2][1]->GetBeginIndex();
   GlobalNumbers[2][2] = FeSpaces[2][2]->GetGlobalNumbers();
   BeginIndex[2][2]    = FeSpaces[2][2]->GetBeginIndex();
 }

// soluble surfactant in the outer phase 
//  if(TDatabase::ParamDB->Compute_Energy)
  {
/* energy spaces */
   Surfact_space = new TFESpace2D(coll, NameString, SurfactString, SurfactBoundCondition,
                                   ContP_USpace, TDatabase::ParamDB->ANSATZ_ORDER, NULL);
   Surfact_space_P1 = new TFESpace2D(Coll_P1, NameString, SurfactString, SurfactBoundCondition,
                                   ContP_USpace, TDatabase::ParamDB->ANSATZ_ORDER, NULL);
   Surfact_space_P2 = new TFESpace2D(Coll_P2, NameString, SurfactString, SurfactBoundCondition,
                                   ContP_USpace, TDatabase::ParamDB->ANSATZ_ORDER, NULL);

   FeSpaces[3][0] = Surfact_space;
   FeSpaces[3][1] = Surfact_space_P1;
   FeSpaces[3][2] = Surfact_space_P2;
   N_DOFs[3][0] =  Surfact_space->GetN_DegreesOfFreedom();
   N_DOFs[3][1] =  Surfact_space_P1->GetN_DegreesOfFreedom();
   N_DOFs[3][2] =  Surfact_space_P2->GetN_DegreesOfFreedom();

   GlobalNumbers[3][0] = FeSpaces[3][0]->GetGlobalNumbers();
   BeginIndex[3][0]    = FeSpaces[3][0]->GetBeginIndex();
   GlobalNumbers[3][1] = FeSpaces[3][1]->GetGlobalNumbers();
   BeginIndex[3][1]    = FeSpaces[3][1]->GetBeginIndex();
   GlobalNumbers[3][2] = FeSpaces[3][2]->GetGlobalNumbers();
   BeginIndex[3][2]    = FeSpaces[3][2]->GetBeginIndex();
  }

// Surfactant on the interface
  {
   IFaceSurfact_space = new TFESpace1D(IFace_Coll , IFaceSString, IFaceSString, FE1D_List);

   IFaceFeSpaces[0] = IFaceSurfact_space;
   N_DOFs[N_FESpaces_All][0] =  IFaceFeSpaces[0]->GetN_DegreesOfFreedom();


//    Surfact_Outputspace = new TFESpace2D(coll, NameString, SurfactString, SurfactBoundCondition,
//                                    ContP_USpace, TDatabase::ParamDB->ANSATZ_ORDER, NULL);
   FeSpaces[4][0] =  new TFESpace2D(coll, NameString, SurfactString, SurfactBoundCondition,
                                   ContP_USpace, TDatabase::ParamDB->ANSATZ_ORDER, NULL);

//    N_S =  Surfact_Outputspace->GetN_DegreesOfFreedom();
   N_DOFs[4][0] =  FeSpaces[4][0]->GetN_DegreesOfFreedom();
  }

// *************************************************************************
// adding aditional vertices for the insinterface/isoboundedge edge
// *************************************************************************
   N_Cells = coll->GetN_Cells();

   for(i=0;i<8;i++)
    N_MovVert[i] = 0;

   for(j=0;j<N_Cells;j++)
    {
     Me = coll->GetCell(j);
     k = Me->GetN_Edges();

     for(l=0;l<k;l++)
      {
       Joint = Me->GetJoint(l);

       if(Joint->GetType() == BoundaryEdge)
        {
         Solid_Joint = (TBoundEdge *)Joint;
         BoundComp = Solid_Joint->GetBoundComp();
         comp=BoundComp->GetID();
         N_MovVert[comp]++;
        }
       else if(Joint->GetType() == IsoInterfaceJoint && Me->GetPhase_ID()==0) // w.r.t inner phase
        {
         IntFaceJoint = (TIsoInterfaceJoint *)Joint;
         BoundComp = IntFaceJoint->GetBoundComp();
         comp=BoundComp->GetID();
         N_MovVert[comp]++;

         Me->GetVertex(l)->GetCoords(TX[0], TY[0]);
         Me->GetVertex((l+1) % k)->GetCoords(TX[1], TY[1]);
//          IntFaceJoint->GeneratemidVert(ORDER-1, TX, TY);
         IntFaceJoint->GenerateVertices(ORDER-1);
        }
      } //  for(l=0;
     } //  for(j=0;j<N_

   for(i=0;i<7;i++)
    cout<< " Number of Boundary edge in comp " << i << " is " << N_MovVert[i]<<endl;

// exit(0);
// store boundary vertices for later use
// assumed that we have only 7 boundary components
   for(i=0;i<7;i++)
    {
     Slip_Joint[i] = new TBoundEdge*[N_MovVert[i]];
     MovBoundVert[i] = new TVertex*[N_MovVert[i]];
     mm[i] = 0;
    }

    N_List[2] = new int[N_MovVert[6]]; // cell list
    N_List[3] = new int[N_MovVert[6]];  // edge list
    Coll_Cells[0] = new TBaseCell*[10*N_MovVert[6]];

   for(j=0;j<N_Cells;j++)
    {
     Me = coll->GetCell(j);
     k = Me->GetN_Edges();

     for(l=0;l<k;l++)
      {
       Joint = Me->GetJoint(l);
       if(Joint->GetType() == BoundaryEdge)
        {
         Solid_Joint = (TBoundEdge *)Joint;
         BoundComp = Solid_Joint->GetBoundComp();
         comp=BoundComp->GetID();

         MovBoundVert[comp][ mm[comp] ] = Me->GetVertex(l);
         Slip_Joint[comp][ mm[comp] ] = (TBoundEdge *)Joint;

         mm[comp]++;
        }
       else if(Joint->GetType() == IsoInterfaceJoint && Me->GetPhase_ID()==0)
        {
         IntFaceJoint = (TIsoInterfaceJoint *)Joint;
         BoundComp = IntFaceJoint->GetBoundComp();
         comp=BoundComp->GetID();

         Coll_Cells[0][mm[comp] ] = Me;
         MovBoundVert[comp][mm[comp] ] = Me->GetVertex(l);
         N_List[2][mm[comp] ] = j;
         N_List[3][mm[comp] ] = l;

         mm[comp]++;
        }

      } //  for(l=0;
     } //  for(j=0;j<N_

//    for(k=0;k<N_MovVert[6];k++)
//    {
//     MovBoundVert[6][k]->GetCoords(x2, y2);
//     cout<< "New x " << x2 <<  "New y " << y2 <<  endl;   
//    }   
// exit(0);

// sort with X0 as the first vertex ordinate if indicator==0 (or) ....
   SortIsoVert_Gen(Coll_Cells[0], MovBoundVert[6], N_List[2], N_List[3],
                   N_MovVert[6], Xi[0], Yi[0], 0);

// sort boundary vertces
    for(k=0;k<N_MovVert[0]-1;k++)
      {
      for(l=k+1;l<N_MovVert[0];l++)
      {
        MovBoundVert[0][k]->GetCoords(x, y);
	MovBoundVert[0][l]->GetCoords(tx, ty);
	if(y < ty)
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
//      for(k=0;k<N_MovVert[0];k++)
//       {
//        MovBoundVert[0][k]->GetCoords(x, y);
//        cout<< " SLPX0 " << x<<" SLPY " << y<<endl;
//        }
//   cout<<endl;


    for(k=0;k<N_MovVert[1]-1;k++)
      {
      for(l=k+1;l<N_MovVert[1];l++)
      {
        MovBoundVert[1][k]->GetCoords(x, y);
	MovBoundVert[1][l]->GetCoords(tx, ty);
	if(y < ty)
	 {
	  temp_Mov = MovBoundVert[1][k];
          MovBoundVert[1][k] = MovBoundVert[1][l];
          MovBoundVert[1][l] = temp_Mov;

	  tempSlip_Joint = Slip_Joint[1][k];
	  Slip_Joint[1][k] = Slip_Joint[1][l];
	  Slip_Joint[1][l] = tempSlip_Joint;
	 }
        }
       }
//      for(k=0;k<N_MovVert[1];k++)
//       {
//        MovBoundVert[1][k]->GetCoords(x, y);
//        cout<< " SLPX1 " << x<<" SLPY " << y<<endl;
//        }
//   cout<<endl;

    for(k=0;k<N_MovVert[2]-1;k++)
      {
      for(l=k+1;l<N_MovVert[2];l++)
      {
        MovBoundVert[2][k]->GetCoords(x, y);
	MovBoundVert[2][l]->GetCoords(tx, ty);
	if(tx < x)
	 {
	  temp_Mov = MovBoundVert[2][k];
          MovBoundVert[2][k] = MovBoundVert[2][l];
          MovBoundVert[2][l] = temp_Mov;

	  tempSlip_Joint = Slip_Joint[2][k];
	  Slip_Joint[2][k] = Slip_Joint[2][l];
	  Slip_Joint[2][l] = tempSlip_Joint;
	 }
        }
       }
//      for(k=0;k<N_MovVert[2];k++)
//       {
//        MovBoundVert[2][k]->GetCoords(x, y);
//        cout<< " SLPX " << x<<" SLPY " << y<<endl;
//        }
//   cout<<endl;
    for(k=0;k<N_MovVert[3]-1;k++)
      {
      for(l=k+1;l<N_MovVert[3];l++)
      {
        MovBoundVert[3][k]->GetCoords(x, y);
	MovBoundVert[3][l]->GetCoords(tx, ty);
	if(ty < y)
	 {
	  temp_Mov = MovBoundVert[3][k];
          MovBoundVert[3][k] = MovBoundVert[3][l];
          MovBoundVert[3][l] = temp_Mov;

	  tempSlip_Joint = Slip_Joint[3][k];
	  Slip_Joint[3][k] = Slip_Joint[3][l];
	  Slip_Joint[3][l] = tempSlip_Joint;
	 }
        }
       }
//      for(k=0;k<N_MovVert[3];k++)
//       {
//        MovBoundVert[3][k]->GetCoords(x, y);
//        cout<< " SLPX1 " << x<<" SLPY " << y<<endl;
//        }
// 
//   cout<<endl;
    for(k=0;k<N_MovVert[4]-1;k++)
      {
      for(l=k+1;l<N_MovVert[4];l++)
      {
        MovBoundVert[4][k]->GetCoords(x, y);
	MovBoundVert[4][l]->GetCoords(tx, ty);
	if(tx > x)
	 {
	  temp_Mov = MovBoundVert[4][k];
          MovBoundVert[4][k] = MovBoundVert[4][l];
          MovBoundVert[4][l] = temp_Mov;

	  tempSlip_Joint = Slip_Joint[4][k];
	  Slip_Joint[4][k] = Slip_Joint[4][l];
	  Slip_Joint[4][l] = tempSlip_Joint;
	 }
        }
       }
//      for(k=0;k<N_MovVert[4];k++)
//       {
//        MovBoundVert[4][k]->GetCoords(x, y);
//        cout<< " SLPX1 " << x<<" SLPY " << y<<endl;
//        }
//   cout<<endl;
    for(k=0;k<N_MovVert[5]-1;k++)
      {
      for(l=k+1;l<N_MovVert[5];l++)
      {
        MovBoundVert[5][k]->GetCoords(x, y);
	MovBoundVert[5][l]->GetCoords(tx, ty);
	if(y < ty)
	 {
	  temp_Mov = MovBoundVert[5][k];
          MovBoundVert[5][k] = MovBoundVert[5][l];
          MovBoundVert[5][l] = temp_Mov;

	  tempSlip_Joint = Slip_Joint[5][k];
	  Slip_Joint[5][k] = Slip_Joint[5][l];
	  Slip_Joint[5][l] = tempSlip_Joint;
	 }
        }
       }
//      for(k=0;k<N_MovVert[5];k++)
//       {
//        MovBoundVert[5][k]->GetCoords(x, y);
//        cout<< " SLPX1 " << x<<" SLPY " << y<<endl;
//        }

// cout << " test main " << endl;
// exit(0);

// *************************************************************************
//   Matrices memory allocation
// *************************************************************************
   /* velocity */
    SqrStruct[0][0] = new TSquareStructure2D(FeSpaces[0][0]);
    SqrStruct[0][1] = new TSquareStructure2D(FeSpaces[0][1]);
    SqrStruct[0][2] = new TSquareStructure2D(FeSpaces[0][2]);
    for(i=0;i<3;i++)
     SqrStruct[0][i]->Sort();

   /* velo - press coupled  B */
    Struct[0][0] = new TStructure2D(FeSpaces[1][0], FeSpaces[0][0]);
    Struct[0][1] = new TStructure2D(FeSpaces[1][1], FeSpaces[0][1]);
    Struct[0][2] = new TStructure2D(FeSpaces[1][2], FeSpaces[0][2]);

   /* velo - press coupled  BT */
    Struct[1][0] = new TStructure2D(FeSpaces[0][0], FeSpaces[1][0]);
    Struct[1][1] = new TStructure2D(FeSpaces[0][1], FeSpaces[1][1]);
    Struct[1][2] = new TStructure2D(FeSpaces[0][2], FeSpaces[1][2]);

 /* grid matrices */
//  if(TDatabase::ParamDB->Moving_Mesh)
  {
    SqrStruct[2][0] = new TSquareStructure2D(FeSpaces[2][0]);
    SqrStruct[2][1] = new TSquareStructure2D(FeSpaces[2][1]);
    SqrStruct[2][2] = new TSquareStructure2D(FeSpaces[2][2]);
  }

 /* surfactant matrices */
//  if(TDatabase::ParamDB->Compute_Energy)
  {
    SqrStruct[3][0] = new TSquareStructure2D(FeSpaces[3][0]);
    SqrStruct[3][1] = new TSquareStructure2D(FeSpaces[3][1]);
    SqrStruct[3][2] = new TSquareStructure2D(FeSpaces[3][2]);
    for(i=0;i<3;i++)
     SqrStruct[3][i]->Sort();
  }

 /* interface surfactant matrices */
  {
    IFaceStruct[0] = new TSquareStructure1D(IFaceFeSpaces[0]);
    IFaceStruct[0]->Sort();
  }

 /* velocity */
  SqMat[0][0][0] = new TSquareMatrix2D(SqrStruct[0][0]); // A11
  SqMat[0][0][1] = new TSquareMatrix2D(SqrStruct[0][0]); // A12
  SqMat[0][0][2] = new TSquareMatrix2D(SqrStruct[0][0]); // A21
  SqMat[0][0][3] = new TSquareMatrix2D(SqrStruct[0][0]); // A22
  SqMat[0][0][4] = new TSquareMatrix2D(SqrStruct[0][0]); // M11
  SqMat[0][0][5] = new TSquareMatrix2D(SqrStruct[0][0]); // M12
  SqMat[0][0][6] = new TSquareMatrix2D(SqrStruct[0][0]); // M21
  SqMat[0][0][7] = new TSquareMatrix2D(SqrStruct[0][0]); // M22
  SqMat[0][0][8] = new TSquareMatrix2D(SqrStruct[0][0]); // F11 interfaceint
  SqMat[0][0][9] = new TSquareMatrix2D(SqrStruct[0][0]); // F22 interfaceint

  Defect = Defect_NSE4;
/* velocity phase_1*/
  SqMat[0][1][0] = new TSquareMatrix2D(SqrStruct[0][1]); // A11
  SqMat[0][1][1] = new TSquareMatrix2D(SqrStruct[0][1]); // A12
  SqMat[0][1][2] = new TSquareMatrix2D(SqrStruct[0][1]); // A21
  SqMat[0][1][3] = new TSquareMatrix2D(SqrStruct[0][1]); // A22
  SqMat[0][1][4] = new TSquareMatrix2D(SqrStruct[0][1]); // M11
  SqMat[0][1][5] = new TSquareMatrix2D(SqrStruct[0][1]); // M12
  SqMat[0][1][6] = new TSquareMatrix2D(SqrStruct[0][1]); // M21
  SqMat[0][1][7] = new TSquareMatrix2D(SqrStruct[0][1]); // M22

/* velocity phase_2*/
  SqMat[0][2][0] = new TSquareMatrix2D(SqrStruct[0][2]); // A11
  SqMat[0][2][1] = new TSquareMatrix2D(SqrStruct[0][2]); // A12
  SqMat[0][2][2] = new TSquareMatrix2D(SqrStruct[0][2]); // A21
  SqMat[0][2][3] = new TSquareMatrix2D(SqrStruct[0][2]); // A22
  SqMat[0][2][4] = new TSquareMatrix2D(SqrStruct[0][2]); // M11
  SqMat[0][2][5] = new TSquareMatrix2D(SqrStruct[0][2]); // M12
  SqMat[0][2][6] = new TSquareMatrix2D(SqrStruct[0][2]); // M21
  SqMat[0][2][7] = new TSquareMatrix2D(SqrStruct[0][2]); // M22

  Mat[0][0]  = new TMatrix2D(Struct[0][0]); // B1
  Mat[0][1]  = new TMatrix2D(Struct[0][0]); // B2
  Mat[0][2]  = new TMatrix2D(Struct[1][0]); // B1T
  Mat[0][3]  = new TMatrix2D(Struct[1][0]); // B2T

 /* phase 1 */
  Mat[1][0]  = new TMatrix2D(Struct[0][1]); // B1
  Mat[1][1]  = new TMatrix2D(Struct[0][1]); // B2
  Mat[1][2]  = new TMatrix2D(Struct[1][1]); // B1T
  Mat[1][3]  = new TMatrix2D(Struct[1][1]); // B2T

 /* phase 2 */
  Mat[2][0]  = new TMatrix2D(Struct[0][2]); // B1
  Mat[2][1]  = new TMatrix2D(Struct[0][2]); // B2
  Mat[2][2]  = new TMatrix2D(Struct[1][2]); // B1T
  Mat[2][3]  = new TMatrix2D(Struct[1][2]); // B2T


//  if(TDatabase::ParamDB->Moving_Mesh)
  {
 /* grid */
  SqMat[2][0][0] = new TSquareMatrix2D(SqrStruct[2][0]); // G11
  SqMat[2][0][1] = new TSquareMatrix2D(SqrStruct[2][0]); // G12
  SqMat[2][0][2] = new TSquareMatrix2D(SqrStruct[2][0]); // G21
  SqMat[2][0][3] = new TSquareMatrix2D(SqrStruct[2][0]); // G22

 /* grid phase_1* */
  SqMat[2][1][0] = new TSquareMatrix2D(SqrStruct[2][1]); // G11
  SqMat[2][1][1] = new TSquareMatrix2D(SqrStruct[2][1]); // G12
  SqMat[2][1][2] = new TSquareMatrix2D(SqrStruct[2][1]); // G21
  SqMat[2][1][3] = new TSquareMatrix2D(SqrStruct[2][1]); // G22

 /* grid phase_2* */
  SqMat[2][2][0] = new TSquareMatrix2D(SqrStruct[2][2]); // G11
  SqMat[2][2][1] = new TSquareMatrix2D(SqrStruct[2][2]); // G12
  SqMat[2][2][2] = new TSquareMatrix2D(SqrStruct[2][2]); // G21
  SqMat[2][2][3] = new TSquareMatrix2D(SqrStruct[2][2]); // G22
 }

/* surfactant */
//  if(TDatabase::ParamDB->Compute_Energy)
  {
  SqMat[3][0][0] = new TSquareMatrix2D(SqrStruct[3][0]); // C_A
  SqMat[3][0][1] = new TSquareMatrix2D(SqrStruct[3][0]); // C_M

  SqMat[3][1][0] = new TSquareMatrix2D(SqrStruct[3][1]); // C_A
  SqMat[3][1][1] = new TSquareMatrix2D(SqrStruct[3][1]); // C_M

  SqMat[3][2][0] = new TSquareMatrix2D(SqrStruct[3][2]); // C_A
  SqMat[3][2][1] = new TSquareMatrix2D(SqrStruct[3][2]); // C_M
  }

 /* interface surfactant matrices */
  {
   SqMat_IFace[0] = new TSquareMatrix1D(IFaceStruct[0]); // Gamma_A
   SqMat_IFace[1] = new TSquareMatrix1D(IFaceStruct[0]); // Gamma_M
  }


   N_U = N_DOFs[0][0];
   N_P = N_DOFs[1][0];
   N_Unknowns = 2*N_U+N_P;

   OutPut("dof velocity : "<< setw(10) << 2*N_U << endl);
   OutPut("dof pressure : "<< setw(10) << N_P << endl);
   OutPut("dof all      : "<<  setw(10) << N_Unknowns  << endl);

//  if(TDatabase::ParamDB->Compute_Energy)
   OutPut("dof surfactant (outer phase) : "<< setw(5) << N_DOFs[3][2] << endl);
// *************************************************************************
//   memory allocation for solution and rhs
// *************************************************************************

 /* velocity */
   Sol[0][0] = new double[N_Unknowns];
   Sol[0][1] = new double[2*N_DOFs[0][1]+N_DOFs[1][1]];
   Sol[0][2] = new double[2*N_DOFs[0][2]+N_DOFs[1][2]];

   oldsol =  new double[N_Unknowns];
   oldrhs =  new double[N_Unknowns];
   B =  new double[N_Unknowns];
   defect =  new double[N_Unknowns];

   Rhs[0][0] = new double[N_Unknowns];
   Rhs[0][1] = new double[2*N_DOFs[0][1]+N_DOFs[1][1]];
   Rhs[0][2] = new double[2*N_DOFs[0][2]+N_DOFs[1][2]];


 /* pressure */
   Sol[1][0] = Sol[0][0]+2*N_U;
   Sol[1][1] = Sol[0][1]+2*N_DOFs[0][1];
   Sol[1][2] = Sol[0][2]+2*N_DOFs[0][2];
   Rhs[1][0] = Rhs[0][0]+2*N_U;
   Rhs[1][1] = Rhs[0][1]+2*N_DOFs[0][1];
   Rhs[1][2] = Rhs[0][2]+2*N_DOFs[0][2];

   memset(Sol[0][0], 0, N_Unknowns*SizeOfDouble);
   memset(Sol[0][1], 0, (2*N_DOFs[0][1]+N_DOFs[1][1])*SizeOfDouble);
   memset(Sol[0][2], 0, (2*N_DOFs[0][2]+N_DOFs[1][2])*SizeOfDouble);
   memset(Rhs[0][0], 0, N_Unknowns*SizeOfDouble);
   memset(Rhs[0][1], 0, (2*N_DOFs[0][1]+N_DOFs[1][1])*SizeOfDouble);
   memset(Rhs[0][2], 0, (2*N_DOFs[0][2]+N_DOFs[1][2])*SizeOfDouble);

//  if(TDatabase::ParamDB->Moving_Mesh)
  {
 /* grid velocity */
   Sol[2][0] = new double[2*N_DOFs[2][0]];
   Sol[2][1] = new double[2*N_DOFs[2][1]];
   Sol[2][2] = new double[2*N_DOFs[2][2]];
   Rhs[2][0] = new double[2*N_DOFs[2][0]];
   Rhs[2][1] = new double[2*N_DOFs[2][1]];
   Rhs[2][2] = new double[2*N_DOFs[2][2]];

   gridsol[0] =new double[2*N_DOFs[2][0]];
   gridsol[1] =new double[2*N_DOFs[2][1]];
   gridsol[2] =new double[2*N_DOFs[2][2]];
   gridrhs[0] =new double[2*N_DOFs[2][0]];
   gridrhs[1] =new double[2*N_DOFs[2][1]];
   gridrhs[2] =new double[2*N_DOFs[2][2]];

   memset(Sol[2][0], 0, 2*N_DOFs[2][0]*SizeOfDouble);
   memset(Sol[2][1], 0, 2*N_DOFs[2][1]*SizeOfDouble);
   memset(Sol[2][2], 0, 2*N_DOFs[2][2]*SizeOfDouble);
   memset(Rhs[2][0], 0, 2*N_DOFs[2][0]*SizeOfDouble);
   memset(Rhs[2][1], 0, 2*N_DOFs[2][1]*SizeOfDouble);
   memset(Rhs[2][2], 0, 2*N_DOFs[2][2]*SizeOfDouble);
   memset(gridsol[0], 0, 2*N_DOFs[2][0]*SizeOfDouble);
   memset(gridsol[1], 0, 2*N_DOFs[2][1]*SizeOfDouble);
   memset(gridsol[2], 0, 2*N_DOFs[2][2]*SizeOfDouble);
   memset(gridrhs[0], 0, 2*N_DOFs[2][0]*SizeOfDouble);
   memset(gridrhs[1], 0, 2*N_DOFs[2][1]*SizeOfDouble);
   memset(gridrhs[2], 0, 2*N_DOFs[2][2]*SizeOfDouble);
  }


/* surfactant */
//  if(TDatabase::ParamDB->Compute_Energy)
  {
   Sol[3][0] = new double[N_DOFs[3][0]];
   Sol[3][1] = new double[N_DOFs[3][1]];
   Sol[3][2] = new double[N_DOFs[3][2]];
   Rhs[3][0] = new double[N_DOFs[3][0]];
   Rhs[3][1] = new double[N_DOFs[3][1]];
   Rhs[3][2] = new double[N_DOFs[3][2]];

   C_defect = new double[N_DOFs[3][2]];
   Csol_old = new double[N_DOFs[3][2]];
   C_B = new double[N_DOFs[3][2]];
   CRhs_old = new double[N_DOFs[3][2]];
   Csol_nonlinearstep = new double[N_DOFs[3][2]];

   memset(Sol[3][0], 0, N_DOFs[3][0]*SizeOfDouble);
   memset(Sol[3][1], 0, N_DOFs[3][1]*SizeOfDouble);
   memset(Sol[3][2], 0, N_DOFs[3][2]*SizeOfDouble);
   memset(Rhs[3][0], 0, N_DOFs[3][0]*SizeOfDouble);
   memset(Rhs[3][1], 0, N_DOFs[3][1]*SizeOfDouble);
   memset(Rhs[3][2], 0, N_DOFs[3][2]*SizeOfDouble);

   memset(C_defect, 0, N_DOFs[3][2]*SizeOfDouble);
   memset(Csol_old, 0, N_DOFs[3][2]*SizeOfDouble);
   memset(C_B, 0, N_DOFs[3][2]*SizeOfDouble);
   
   outputsol = new double[N_DOFs[2][0]];
   memset(outputsol, 0, N_DOFs[2][0]*SizeOfDouble);
  }


 /* interface surfactant */
  {
   Sol[N_FESpaces_All][0] = new double[N_DOFs[N_FESpaces_All][0]];
   Rhs[N_FESpaces_All][0] = new double[N_DOFs[N_FESpaces_All][0]];

   I_defect = new double[N_DOFs[N_FESpaces_All][0]];
   Isol_old = new double[N_DOFs[N_FESpaces_All][0]];
   I_B = new double[N_DOFs[N_FESpaces_All][0]];
   IRhs_old = new double[N_DOFs[N_FESpaces_All][0]];

   memset(Sol[N_FESpaces_All][0], 0, N_DOFs[N_FESpaces_All][0]*SizeOfDouble);
   memset(Rhs[N_FESpaces_All][0], 0, N_DOFs[N_FESpaces_All][0]*SizeOfDouble);

   memset(I_defect, 0, N_DOFs[N_FESpaces_All][0]*SizeOfDouble);
   memset(Isol_old, 0, N_DOFs[N_FESpaces_All][0]*SizeOfDouble);
   memset(I_B, 0, N_DOFs[N_FESpaces_All][0]*SizeOfDouble);

   Sol[4][0] = new double[ N_DOFs[4][0]];
   memset(Sol[4][0], 0,  N_DOFs[4][0]*SizeOfDouble);
  }


  /* velocity and pressure */
  VeloVect[0][0] = new TFEVectFunct2D(FeSpaces[0][0], UString,  UString,  Sol[0][0], N_U, 2);
  FeFunct[0][0] = VeloVect[0][0]->GetComponent(0);
  FeFunct[1][0] = VeloVect[0][0]->GetComponent(1);
  FeFunct[2][0] = new TFEFunction2D(FeSpaces[1][0], PString,  PString, Sol[0][0]+2*N_U, N_P);

  FeFunct[0][0]->Interpolate(InitialU1);
  FeFunct[1][0]->Interpolate(InitialU2);

  VeloVect[0][1] = new TFEVectFunct2D(FeSpaces[0][1], UString,  UString,  Sol[0][1], N_DOFs[0][1], 2);
  FeFunct[0][1] = VeloVect[0][1]->GetComponent(0);
  FeFunct[1][1] = VeloVect[0][1]->GetComponent(1);
  FeFunct[2][1] = new TFEFunction2D(FeSpaces[1][1], PString,  PString,
                                    Sol[0][1]+2*N_DOFs[0][1], N_DOFs[1][1]);

  VeloVect[0][2] = new TFEVectFunct2D(FeSpaces[0][2], UString,  UString,  Sol[0][2], N_DOFs[0][2], 2);
  FeFunct[0][2] = VeloVect[0][2]->GetComponent(0);
  FeFunct[1][2] = VeloVect[0][2]->GetComponent(1);
  FeFunct[2][2] = new TFEFunction2D(FeSpaces[1][2], PString,  PString,
                                    Sol[0][2]+2*N_DOFs[0][2], N_DOFs[1][2]);

//  if(TDatabase::ParamDB->Moving_Mesh)
  {
 /* grid velocity */
   VeloVect[1][0] = new TFEVectFunct2D(FeSpaces[2][0], gridString,  gridString,  Sol[2][0], N_DOFs[2][0], 2);
   FeFunct[3][0] = VeloVect[1][0]->GetComponent(0);
   FeFunct[4][0] = VeloVect[1][0]->GetComponent(1);

   VeloVect[1][1] = new TFEVectFunct2D(FeSpaces[2][1], gridString,  gridString,  Sol[2][1], N_DOFs[2][1], 2);
   FeFunct[3][1] = VeloVect[1][1]->GetComponent(0);
   FeFunct[4][1] = VeloVect[1][1]->GetComponent(1);

   VeloVect[1][2] = new TFEVectFunct2D(FeSpaces[2][2], gridString,  gridString,  Sol[2][2], N_DOFs[2][2], 2);
   FeFunct[3][2] = VeloVect[1][2]->GetComponent(0);
   FeFunct[4][2] = VeloVect[1][2]->GetComponent(1);
  }

/* Surfactant */
//  if(TDatabase::ParamDB->Compute_Energy)
  {
   FeFunct[5][0] = new TFEFunction2D(FeSpaces[3][0], SurfactString,  SurfactString, Sol[3][0], N_DOFs[3][0]);
   FeFunct[5][1] = new TFEFunction2D(FeSpaces[3][1], SurfactString,  SurfactString, Sol[3][1], N_DOFs[3][1]);
   FeFunct[5][2] = new TFEFunction2D(FeSpaces[3][2], SurfactString,  SurfactString, Sol[3][2], N_DOFs[3][2]);

   FeFunct[5][2]->Interpolate(InitialS_Outer);
  }

 /* interface surfactant */
  {
   IFaceFeFunct[0] = new TFEFunction1D(IFaceFeSpaces[0], IFaceSString, IFaceSString, Sol[N_FESpaces_All][0],
                                       N_DOFs[N_FESpaces_All][0]);


   IFaceFeFunct[0]->Interpolate(InitialS);


   FeFunct[6][0]  = new TFEFunction2D(FeSpaces[4][0], IFaceSString,  IFaceSString, Sol[4][0],  N_DOFs[4][0]);
  }

  OutputFeunction = new TFEFunction2D(FeSpaces[2][0], gridString,  gridString, outputsol, N_DOFs[2][0]);

  /*assemble matrix for grid moving - begin*/
  t1 = GetTime();

  fesp[0] = FeSpaces[2][0];
  SQMATRICES_GRID[0] = SqMat[2][0][0];
  SQMATRICES_GRID[0]->Reset();
  SQMATRICES_GRID[1] = SqMat[2][0][1];
  SQMATRICES_GRID[1]->Reset();
  SQMATRICES_GRID[2] = SqMat[2][0][2];
  SQMATRICES_GRID[2]->Reset();
  SQMATRICES_GRID[3] = SqMat[2][0][3];
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

  Entries[0] = SqMat[2][0][0]->GetEntries();
  Entries[1] = SqMat[2][0][1]->GetEntries();
  Entries[2] = SqMat[2][0][2]->GetEntries();
  Entries[3] = SqMat[2][0][3]->GetEntries();
  GridKCol = SqrStruct[2][0]->GetKCol();
  GridRowPtr = SqrStruct[2][0]->GetRowPtr();
  memset(Entries[1] + GridRowPtr[N_DOFs[2][0]-Bound_DOFs[2][0]], 0,
         Bound_DOFs[2][0]*SizeOfDouble);
  memset(Entries[2] + GridRowPtr[N_DOFs[2][0]-Bound_DOFs[2][0]], 0,
         Bound_DOFs[2][0]*SizeOfDouble);

 t2 = GetTime();

  if(TDatabase::ParamDB->SC_VERBOSE > 1)
  {
    cout << "Grid assembling done"<< endl;
    OutPut("Time for Grid assembling: " << t2-t1 << endl);
  }

// inner Phase
  t1 = GetTime();

  fesp[0] = FeSpaces[2][1];
  SQMATRICES_GRID[0] = SqMat[2][1][0];
  SQMATRICES_GRID[0]->Reset();
  SQMATRICES_GRID[1] = SqMat[2][1][1];
  SQMATRICES_GRID[1]->Reset();
  SQMATRICES_GRID[2] = SqMat[2][1][2];
  SQMATRICES_GRID[2]->Reset();
  SQMATRICES_GRID[3] = SqMat[2][1][3];
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

  Entries_P1[0] = SqMat[2][1][0]->GetEntries();
  Entries_P1[1] = SqMat[2][1][1]->GetEntries();
  Entries_P1[2] = SqMat[2][1][2]->GetEntries();
  Entries_P1[3] = SqMat[2][1][3]->GetEntries();
  GridKCol_P1 = SqrStruct[2][1]->GetKCol();
  GridRowPtr_P1 = SqrStruct[2][1]->GetRowPtr();
  memset(Entries_P1[1] + GridRowPtr_P1[N_DOFs[2][1]-Bound_DOFs[2][1]], 0,
         Bound_DOFs[2][1]*SizeOfDouble);
  memset(Entries_P1[2] + GridRowPtr_P1[N_DOFs[2][1]-Bound_DOFs[2][1]], 0,
         Bound_DOFs[2][1]*SizeOfDouble);

 t2 = GetTime();

  if(TDatabase::ParamDB->SC_VERBOSE > 1)
  {
    cout << "Internal Phase Grid assembling done"<< endl;
    OutPut("Time for Internal Grid assembling: " << t2-t1 << endl);
  }

/* outer Phase */
  t1 = GetTime();

  fesp[0] = FeSpaces[2][2];
  SQMATRICES_GRID[0] = SqMat[2][2][0];
  SQMATRICES_GRID[0]->Reset();
  SQMATRICES_GRID[1] = SqMat[2][2][1];
  SQMATRICES_GRID[1]->Reset();
  SQMATRICES_GRID[2] = SqMat[2][2][2];
  SQMATRICES_GRID[2]->Reset();
  SQMATRICES_GRID[3] = SqMat[2][2][3];
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

  Entries_P2[0] = SqMat[2][2][0]->GetEntries();
  Entries_P2[1] = SqMat[2][2][1]->GetEntries();
  Entries_P2[2] = SqMat[2][2][2]->GetEntries();
  Entries_P2[3] = SqMat[2][2][3]->GetEntries();
  GridKCol_P2 = SqrStruct[2][2]->GetKCol();
  GridRowPtr_P2 = SqrStruct[2][2]->GetRowPtr();
  memset(Entries_P2[1] + GridRowPtr_P2[N_DOFs[2][2]-Bound_DOFs[2][2]], 0,
         Bound_DOFs[2][2]*SizeOfDouble);
  memset(Entries_P2[2] + GridRowPtr_P2[N_DOFs[2][2]-Bound_DOFs[2][2]], 0,
         Bound_DOFs[2][2]*SizeOfDouble);

 t2 = GetTime();

  if(TDatabase::ParamDB->SC_VERBOSE > 1)
  {
    cout << "Internal Phase Grid assembling done"<< endl;
    OutPut("Time for Internal Grid assembling: " << t2-t1 << endl);
  }

  /*assemble matrix for grid moving - end*/

  refpos = new double[2*N_DOFs[2][0]];
  auxpos = new double[2*N_DOFs[2][0]];
  pos = new double[2*N_DOFs[2][0]];
  tmp = new double[2*N_DOFs[2][0]];

  memset(tmp, 0, 2*N_DOFs[2][0]*SizeOfDouble);

  refpos_P1 = new double[2*N_DOFs[2][1]];
  auxpos_P1 = new double[2*N_DOFs[2][1]];
  pos_P1 = new double[2*N_DOFs[2][1]];

  refpos_P2 = new double[2*N_DOFs[2][2]];
  auxpos_P2 = new double[2*N_DOFs[2][2]];
  pos_P2 = new double[2*N_DOFs[2][2]];


  RefGridPos = new TFEVectFunct2D(FeSpaces[2][0], refposString, gridString,
                                  refpos, N_DOFs[2][0], 2);

  AuxGridPos = new TFEVectFunct2D(FeSpaces[2][0], auxposString, gridString,
                                  auxpos, N_DOFs[2][0], 2);
  GridPos = new TFEVectFunct2D(FeSpaces[2][0], posString, gridString,
                                 pos, N_DOFs[2][0], 2);


  RefGridPos_P1 = new TFEVectFunct2D(FeSpaces[2][1], refposString, refposString,
                                  refpos_P1, N_DOFs[2][1], 2);
  AuxGridPos_P1 = new TFEVectFunct2D(FeSpaces[2][1], auxposString, auxposString,
                                  auxpos_P1, N_DOFs[2][1], 2);
  GridPos_P1 = new TFEVectFunct2D(FeSpaces[2][1], posString, posString,
                                  pos_P1, N_DOFs[2][1], 2);

  RefGridPos_P2 = new TFEVectFunct2D(FeSpaces[2][2], refposString, gridString,
                                  refpos_P2, N_DOFs[2][2], 2);
  AuxGridPos_P2 = new TFEVectFunct2D(FeSpaces[2][2], auxposString, gridString,
                                  auxpos_P2, N_DOFs[2][2], 2);
  GridPos_P2 = new TFEVectFunct2D(FeSpaces[2][2], posString, gridString,
                                  pos_P2, N_DOFs[2][2], 2);

   GridPos->GridToData();
   memcpy(refpos, pos, 2*N_DOFs[2][0]*SizeOfDouble); 
   memcpy(auxpos, pos, 2*N_DOFs[2][0]*SizeOfDouble); 
   
   GridPos_P1->GridToData();
   memcpy(refpos_P1, pos_P1, 2*N_DOFs[2][1]*SizeOfDouble); 
   memcpy(auxpos_P1, pos_P1, 2*N_DOFs[2][1]*SizeOfDouble); 
   
   GridPos_P2->GridToData();
   memcpy(refpos_P2, pos_P2, 2*N_DOFs[2][2]*SizeOfDouble); 
   memcpy(auxpos_P2, pos_P2, 2*N_DOFs[2][2]*SizeOfDouble);    

  
 // prepare output (maxn_fespaces,  maxn_scalar,  maxn_vect, maxn_parameters, domain)
   Output = new TOutput2D(4, 5, 1, 2, Domain);
   
#ifdef __FLUIDFLOW__  
   Output->AddFEVectFunct(VeloVect[0][0]);
   Output->AddFEFunction(FeFunct[2][0]);   
#endif
 
// mesh velocity
//    Output->AddFEVectFunct(VeloVect[1][0]);

 /* surfactant  */
#ifdef __WITHSURFACTANT__ 
//    if(TDatabase::ParamDB->Compute_Energy)
   Output->AddFEFunction(FeFunct[5][0]);

 /* interface surfactant  */
   Output->AddFEFunction(FeFunct[6][0]);
#endif
 /* rise velo in bubble and gridvelo in outer  */
   Output->AddFEFunction(OutputFeunction);
   
   os.seekp(std::ios::beg);
   Output->AddParameter(TDatabase::TimeDB->CURRENTTIME,os.str().c_str());
/********************************************************************************************************/

  t = TDatabase::TimeDB->CURRENTTIME-TDatabase::TimeDB->STARTTIME;

  if(TDatabase::ParamDB->WRITE_VTK)
   {
       // bulk phase surfactants
       length = BeginIndex[3][0][1] -  BeginIndex[3][0][0];
       for(i=0;i<N_Cells_P2;i++)
        {
         DOF = GlobalNumbers[3][0] + BeginIndex[3][0][GlobalCell_Index[1][i]];
         DOF_P2 = GlobalNumbers[3][2] + BeginIndex[3][2][i];

         for(jj=0;jj<length;jj++)
          {
           k = DOF[jj];
           l1 = DOF_P2[jj];
           Sol[3][0][k] = Sol[3][2][l1];
          }
         }
    // interface surfactant
    MapSurfToDomain(IFaceFeFunct[0], FeFunct[6][0], N_List[0], N_List[1]);

    os.seekp(std::ios::beg);
    if(img<10) os << "VTK/"<<VtkBaseName<<".0000"<<img<<".vtk" << ends;
    else if(img<100) os << "VTK/"<<VtkBaseName<<".000"<<img<<".vtk" << ends;
    else if(img<1000) os << "VTK/"<<VtkBaseName<<".00"<<img<<".vtk" << ends;
    else if(img<10000) os << "VTK/"<<VtkBaseName<<".0"<<img<<".vtk" << ends;
    else  os << "VTK/"<<VtkBaseName<<"."<<img<<".vtk" << ends;
    Output->WriteVtk(os.str().c_str());
    img++;
    
    PrintSurfSurfactant(N_MovVert[6], MovBoundVert[6], FeFunct[6][0], N_BData);    
   }
//    cout << " test " << endl;
//    exit(0);
   
   
  m = 0;
  end_time = TDatabase::TimeDB->ENDTIME;
  if (TDatabase::TimeDB->TIMESTEPLENGTH_CONTROL)
    time_discs = 2;
  else
    time_discs = 1;

  N_SubSteps = GetN_SubSteps();
  solver_time = 0.0;
  N_LinIter = 0;

//  if(TDatabase::ParamDB->Compute_Energy)
#ifdef __WITHSURFACTANT__ 
  memcpy(Csol_old, Sol[3][2], N_DOFs[3][2]*SizeOfDouble);
  memcpy(Isol_old, Sol[N_FESpaces_All][0], N_DOFs[N_FESpaces_All][0]*SizeOfDouble);


  GetSurfactMass(FeFunct[6][0], IFaceFeFunct[0], N_List[0], N_List[1], Surf_Mass);

  Get_FeFunction2DMass(FeFunct[5][2], Params);

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
#else
  GetSurfactMass(FeFunct[6][0], IFaceFeFunct[0], N_List[0], N_List[1], Surf_Mass); // for sphericity
#endif
    Get_KE(VeloVect[0][0], Params);
    InitVolume = CurrVolume = 2.*Pi*Params[0];
    MovBoundVert[0][0]->GetCoords(Lx, Ly);
    MovBoundVert[1][0]->GetCoords(Rx, Ry);

    radius = pow((3.*CurrVolume/(4.*Pi)), (1./3.));
    //cout<< "CurrVolume : " << CurrVolume << " radius : " << radius <<endl;
    sphericity =  4.*Pi*radius*radius/(Surf_Mass[1]);
    
    OutPut(setw(20)<<"Time, Volume , Volume Diff, Top, Bottom: " << TDatabase::TimeDB->CURRENTTIME
           <<"   "<< CurrVolume<<"   "<< CurrVolume - InitVolume<<"   "<< Ly<<"   "<<Ry<< endl);
    OutPut(setw(20)<<"Time, KE , x_mass, y_mass, U1_Rise, U2_Rise, SurfArea, sphericity: " 
           << TDatabase::TimeDB->CURRENTTIME
           <<"   "<< Params[1]<<"   "<< Params[2] <<"   "<< Params[3]<<"   "<< Params[4]
           <<"   "<<Params[5]<<"   "<< Surf_Mass[1]<<"   "<< sphericity<<endl);
 
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
         cout<< " very_first_time : "<<very_first_time<<endl;
        }

        tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
        if (!very_first_time)
          TDatabase::TimeDB->CURRENTTIME += tau;

        if (m==1 && l==0)
            oldtau=tau;

       // cout << " tau " << tau << endl;
#ifdef __FLUIDFLOW__  
        // save current solution to the previous time step sol
       memcpy(oldsol, Sol[0][0], SizeOfDouble*N_Unknowns);
      // working array for rhs is B, initialize B
       memset(B, 0, N_Unknowns*SizeOfDouble);   
       
      //=======================================================================
      // get mesh velo
      // inner domain velocity
      GetGridVelocity(Entries_P1, gridsol[1], gridrhs[1],
                      GridKCol_P1, GridRowPtr_P1,
                      GridPos_P1, AuxGridPos_P1,
                      VeloVect[0][0], tau,
                      VeloVect[1][1], GlobalCell_Index[0]);

      //  outer domain velocity
      GetGridVelo_outer(Entries_P2, gridsol[2], gridrhs[2],
                        GridKCol_P2, GridRowPtr_P2,
                        GridPos_P2, AuxGridPos_P2,
                        VeloVect[0][0], tau,
                        VeloVect[1][2], GlobalCell_Index[1]);

      length = BeginIndex[2][0][1] -  BeginIndex[2][0][0];
      for(i=0;i<N_Cells_P2;i++)
       {
        DOF = GlobalNumbers[2][0] + BeginIndex[2][0][GlobalCell_Index[1][i]];
        DOF_P2 = GlobalNumbers[2][2] + BeginIndex[2][2][i];

        for(j=0;j<length;j++)
         {
          k = DOF[j];
          l1 = DOF_P2[j];
          Sol[2][0][k] = Sol[2][2][l1];
          Sol[2][0][k+N_DOFs[2][0]] = Sol[2][2][l1+N_DOFs[2][2]];
         }
       }
      for(i=0;i<N_Cells_P1;i++)
       {
        DOF = GlobalNumbers[2][0] + BeginIndex[2][0][GlobalCell_Index[0][i]];
        DOF_P1 = GlobalNumbers[2][1] + BeginIndex[2][1][i];
        for(j=0;j<length;j++)
         {
          k = DOF[j];
          l1 = DOF_P1[j];
          Sol[2][0][k] = Sol[2][1][l1];
          Sol[2][0][k+N_DOFs[2][0]] = Sol[2][1][l1+N_DOFs[2][1]];
//           cout << "val1: " << Sol[2][1][l1] << " va2:" <<  Sol[2][1][l1+N_DOFs[2][1]] << endl;
         }
       }     

    //=======================================================================
    // Assemble matrices
    //=======================================================================

       DiscreteForm = DiscreteFormGalerkin;

       SQMATRICES[0] = SqMat[0][0][0]; // A11
       SQMATRICES[1] = SqMat[0][0][1]; // A12
       SQMATRICES[2] = SqMat[0][0][2]; // A21
       SQMATRICES[3] = SqMat[0][0][3]; // A22
       SQMATRICES[4] = SqMat[0][0][4]; // M11
       SQMATRICES[5] = SqMat[0][0][7]; // M22

       SqMat[0][0][5]->Reset();  // M12
       SqMat[0][0][6]->Reset();  // M21

       MATRICES[0] = Mat[0][0]; // B1
       MATRICES[1] = Mat[0][1]; // B2
       MATRICES[2] = Mat[0][2]; // B1T
       MATRICES[3] = Mat[0][3]; // B2T

       SQMATRICES[0]->Reset();
       SQMATRICES[1]->Reset();
       SQMATRICES[2]->Reset();
       SQMATRICES[3]->Reset();
       SQMATRICES[4]->Reset();
       SQMATRICES[5]->Reset();

       MATRICES[0]->Reset();
       MATRICES[1]->Reset();
       MATRICES[2]->Reset();
       MATRICES[3]->Reset();

       N_SquareMatrices = 6;
       N_RectMatrices = 4;

       N_Rhs = 2;
       N_FESpaces = 3;

       fesp[0] = FeSpaces[0][0];
       fesp[1] = FeSpaces[1][0];
       fesp[2] = FeSpaces[2][0];

       fefct[0] = FeFunct[0][0];
       fefct[1] = FeFunct[1][0];
       fefct[2] = FeFunct[3][0];
       fefct[3] = FeFunct[4][0];

       ferhs[0] = FeSpaces[0][0];
       ferhs[1] = FeSpaces[0][0];

       RHSs[0] = Rhs[0][0];
       RHSs[1] = RHSs[0] + N_U;
       RHSs[2] = RHSs[0] + 2*N_U;

       memset(RHSs[0], 0, N_Unknowns*SizeOfDouble);

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


      SqMat[0][0][8]->Reset(); // F11
      SqMat[0][0][9]->Reset(); // F22   
      
#ifdef __WITHSURFACTANT__ 
      MapSurfToDomain(IFaceFeFunct[0], FeFunct[6][0], N_List[0], N_List[1]);
#endif      
      FreeSurf_2PhaseSurfAxial3D(SqMat[0][0][8], SqMat[0][0][9],
                                 RHSs[0], RHSs[0]+N_U,
                                 BoundCondition, tau, 0, FeFunct[6][0]);
     
      // Adding freesurf entries to A11 and A22
      MatAdd(SqMat[0][0][0], SqMat[0][0][8], 1.);
      MatAdd(SqMat[0][0][3], SqMat[0][0][9], 1.);

      N_Active = FeSpaces[0][0]->GetActiveBound();
        // get row in off diagonal matrix where the Dirichlet nodes start
      RowPtr = SqMat[0][0][0]->GetRowPtr();
        // compute number of entries starting from this row to the end
        // of the matrix
      j = RowPtr[N_Active];
      k = RowPtr[N_U]-j;
        // get number of active dof
        // set these entries to zero
      memset(SqMat[0][0][1]->GetEntries()+j, 0, SizeOfDouble*k);
      memset(SqMat[0][0][2]->GetEntries()+j, 0, SizeOfDouble*k);

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

          SQMATRICES[0] = SqMat[0][0][0]; // A11
          SQMATRICES[1] = SqMat[0][0][3]; // A22
          SQMATRICES[2] = SqMat[0][0][1]; // A12
          SQMATRICES[3] = SqMat[0][0][2]; // A21
          SQMATRICES[4] = SqMat[0][0][4]; // M11
          SQMATRICES[5] = SqMat[0][0][7]; // M22
          SQMATRICES[6] = SqMat[0][0][5]; // M12
          SQMATRICES[7] = SqMat[0][0][6]; // M21

          MATRICES[0] =Mat[0][2]; // B1T
          MATRICES[1] = Mat[0][3]; // B2T

          fesp[0] = FeSpaces[0][0];
          ferhs[0] = FeSpaces[0][0];
          ferhs[1] = FeSpaces[0][0];

          RHSs[0] = Rhs[0][0];
          RHSs[1] = RHSs[0] + N_U;

          aux = new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);

          Assemble2DSlipBC(N_FESpaces, fesp,
                           N_SquareMatrices, SQMATRICES,
                           N_RectMatrices, MATRICES,
                           N_Rhs, RHSs, ferhs,
                           DiscreteForm,
                           BoundaryConditions,
                           BoundValues,
                           aux, FeFunct[0][0],
                           FeFunct[1][0]);

//           TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION = 1;
       delete aux;

        } // if (TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION >= 1)


        Dscal(Mat[0][0]->GetN_Entries(), tau, Mat[0][0]->GetEntries() );  // B1
        Dscal(Mat[0][1]->GetN_Entries(), tau, Mat[0][1]->GetEntries());  // B2
        Dscal(Mat[0][2]->GetN_Entries(), tau, Mat[0][2]->GetEntries());  // B1T
        Dscal(Mat[0][3]->GetN_Entries(), tau, Mat[0][3]->GetEntries());  // B2T

        gamma = 0.;

        // since rhs depends on moving grid so its better to use current geo
        Daxpy(N_Active, tau, Rhs[0][0], B);
        Daxpy(N_Active, tau, Rhs[0][0]+N_U, B+N_U);

        // update rhs by Laplacian and convective term initialy
        // by current time step
        // scaled by current sub time step length and theta2
        // currently : M := M + gamma A
        // M = M + (-gamma - tau*TDatabase::TimeDB->THETA2) A
        MatAdd(SqMat[0][0][4], SqMat[0][0][0],
                     -gamma - tau*TDatabase::TimeDB->THETA2);
        MatAdd(SqMat[0][0][5], SqMat[0][0][1],
                     -gamma - tau*TDatabase::TimeDB->THETA2);
        MatAdd(SqMat[0][0][6], SqMat[0][0][2],
                     -gamma - tau*TDatabase::TimeDB->THETA2);
        MatAdd(SqMat[0][0][7], SqMat[0][0][3],
                     -gamma - tau*TDatabase::TimeDB->THETA2);

       // set current factor of steady state matrix
        gamma = -tau*TDatabase::TimeDB->THETA2;

        memset(defect, 0, N_Unknowns*SizeOfDouble);
        MatVectActive(SqMat[0][0][4], oldsol, defect);
        Daxpy(N_Active, 1, defect, B);
        memset(defect, 0, N_Unknowns*SizeOfDouble);
        MatVectActive(SqMat[0][0][5], oldsol+N_U, defect);
        Daxpy(N_Active, 1, defect, B);
        memset(defect, 0, N_Unknowns*SizeOfDouble);
        MatVectActive(SqMat[0][0][6], oldsol, defect+N_U);
        Daxpy(N_Active, 1, defect+N_U, B+N_U);
        memset(defect, 0, N_Unknowns*SizeOfDouble);
        MatVectActive(SqMat[0][0][7], oldsol+N_U, defect+N_U);
        Daxpy(N_Active, 1, defect+N_U, B+N_U);

//     cout << "sol " << Ddot(N_Unknowns, B, B)<<endl;
//     exit(0);
       //      for(i=0; i<N_U; i++)
       //       cout << i <<" rhs 1 " << B[i] << " rhs 2 " << B[i+N_U] << endl;
       // exit(0);
       // set Dirichlet values
       // RHSs[0] still available from assembling
       memcpy(B+N_Active, RHSs[0]+N_Active, (N_U-N_Active)*SizeOfDouble);
       memcpy(B+N_Active+N_U, RHSs[1]+N_Active,(N_U-N_Active)*SizeOfDouble);
  
       // copy Dirichlet values from rhs into Sol[0][mg_level-1]
       memcpy(Sol[0][0]+N_Active, RHSs[0]+N_Active, (N_U-N_Active)*SizeOfDouble);
       memcpy(Sol[0][0]+N_Active+N_U, RHSs[1]+N_Active, (N_U-N_Active)*SizeOfDouble);


   //=====================================================================
   // the stiffness matrix is stored on M11, (M12, M21, M22)
   // assembling of system matrix
   //========================================================================
    MatAdd(SqMat[0][0][4], SqMat[0][0][0],
                     -gamma + tau*TDatabase::TimeDB->THETA1);
    MatAdd(SqMat[0][0][5], SqMat[0][0][1],
                     -gamma + tau*TDatabase::TimeDB->THETA1);
    MatAdd(SqMat[0][0][6], SqMat[0][0][2],
                     -gamma + tau*TDatabase::TimeDB->THETA1);
    MatAdd(SqMat[0][0][7], SqMat[0][0][3],
                     -gamma + tau*TDatabase::TimeDB->THETA1);

       // set current factor of steady state matrix
    gamma = tau*TDatabase::TimeDB->THETA1;   
        
    //======================================================================
    OutPut(endl << "CURRENT TIME: ");
    OutPut(TDatabase::TimeDB->CURRENTTIME << endl);
   //=====================================================================   
    
  //======================================================================
  // nonlinear loop
  //======================================================================
    N_LinIterCurr = 0;
    solver_time_curr = 0;

    for(j=0;j<Max_It;j++)
     {
      memcpy(oldsol, Sol[0][0], SizeOfDouble*N_Unknowns);
      memset(defect, 0, N_Unknowns*SizeOfDouble);

      SQMATRICES[0] = SqMat[0][0][4];
      SQMATRICES[1] = SqMat[0][0][5];
      SQMATRICES[2] = SqMat[0][0][6];
      SQMATRICES[3] = SqMat[0][0][7];
      MATRICES[0] = Mat[0][0];
      MATRICES[1] = Mat[0][1];
      MATRICES[2] = Mat[0][2];
      MATRICES[3] = Mat[0][3];

      if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE) 
        IntoL20FEFunction(Sol[0][0]+2*N_U, N_P, FeSpaces[1][0],
                        velocity_space_code, pressure_space_code);

      Defect(sqmatrices,matrices,Sol[0][0],B,defect);

//      for(i=0; i<N_U; i++)
//       cout << i<< " defect 1 " << Sol[0][0][i] << " defect 2 " << Sol[0][0][i+N_U] << endl;
//        cout << i<< "defect 1 " << defect[i] << " defect 2 " << defect[i+N_U] << endl;
//       cout << i<< " rhs 1 " << B[i] << " rhs 2 " << B[i+N_U] << endl;
// // exit(0);

      if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE) 
          IntoL20Vector2D(defect+2*N_U, N_P,pressure_space_code);

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


      if((((sqrt(residual)<=limit)||(j==Max_It-1))) && (j>=TDatabase::ParamDB->SC_MINIT))
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
       t1 = GetTime();
       DirectSolver(SQMATRICES[0], SQMATRICES[1],
                    SQMATRICES[2], SQMATRICES[3],
                    MATRICES[2], MATRICES[3],
                    MATRICES[0], MATRICES[1],
                    B, Sol[0][0]);

       t2 = GetTime();
       solver_time_curr = t2-t1;
       solver_time += solver_time_curr;

     //======================================================================
     // end solve linear system
     //======================================================================

     // restore mass matrices by subtracting the A-matrices
        MatAdd(SqMat[0][0][4], SqMat[0][0][0], -gamma);
        MatAdd(SqMat[0][0][5], SqMat[0][0][1], -gamma);
        MatAdd(SqMat[0][0][6], SqMat[0][0][2], -gamma);
        MatAdd(SqMat[0][0][7], SqMat[0][0][3], -gamma);

        gamma = 0;
     //======================================================================
     // Grid velocity
     //======================================================================
      // inner domain velocity
      GetGridVelocity(Entries_P1, gridsol[1], gridrhs[1],
                      GridKCol_P1, GridRowPtr_P1,
                      GridPos_P1, AuxGridPos_P1,
                      VeloVect[0][0], tau,
                      VeloVect[1][1], GlobalCell_Index[0]);

//  outer domain velocity
      GetGridVelo_outer(Entries_P2, gridsol[2], gridrhs[2],
                        GridKCol_P2, GridRowPtr_P2,
                        GridPos_P2, AuxGridPos_P2,
                        VeloVect[0][0], tau,
                        VeloVect[1][2], GlobalCell_Index[1]);

      length = BeginIndex[2][0][1] -  BeginIndex[2][0][0];
      for(i=0;i<N_Cells_P2;i++)
       {
        DOF = GlobalNumbers[2][0] + BeginIndex[2][0][GlobalCell_Index[1][i]];
        DOF_P2 = GlobalNumbers[2][2] + BeginIndex[2][2][i];

        for(jj=0;jj<length;jj++)
         {
          k = DOF[jj];
          l1 = DOF_P2[jj];
          Sol[2][0][k] = Sol[2][2][l1];
          Sol[2][0][k+N_DOFs[2][0]] = Sol[2][2][l1+N_DOFs[2][2]];
         }
       }

      for(i=0;i<N_Cells_P1;i++)
       {
        DOF = GlobalNumbers[2][0] + BeginIndex[2][0][GlobalCell_Index[0][i]];
        DOF_P1 = GlobalNumbers[2][1] + BeginIndex[2][1][i];

        for(jj=0;jj<length;jj++)
         {
          k = DOF[jj];
          l1 = DOF_P1[jj];
          Sol[2][0][k] = Sol[2][1][l1];
          Sol[2][0][k+N_DOFs[2][0]] = Sol[2][1][l1+N_DOFs[2][1]];
         }
       }

     //======================================================================
     // assemble new matrix due to nonlinearity
     //======================================================================
        DiscreteForm = DiscreteFormNLGalerkin;

        N_RectMatrices = 0;
        N_Rhs = 0;
        N_FESpaces = 3;
        N_SquareMatrices = 2;

        SQMATRICES[0] = SqMat[0][0][0];
        SQMATRICES[1] = SqMat[0][0][3];
        SQMATRICES[0]->Reset();
        SQMATRICES[1]->Reset();

        fesp[0] = FeSpaces[0][0];
        fesp[1] = FeSpaces[1][0];
        fesp[2] = FeSpaces[2][0];

        fefct[0] = FeFunct[0][0]; // u1
        fefct[1] = FeFunct[1][0]; // u2
        fefct[2] = FeFunct[3][0]; // w1
        fefct[3] = FeFunct[4][0]; // w1

      //======================================================================
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
        MatAdd(SqMat[0][0][0], SqMat[0][0][8], 1.);
        MatAdd(SqMat[0][0][3], SqMat[0][0][9], 1.);

        delete aux;

        // slip type bc detected, modify matrices accordingly
       if (TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION >= 1)
        {
          N_FESpaces = 1;
          N_SquareMatrices = 2;
          N_RectMatrices = 0;
          N_Rhs = 2;
          DiscreteForm = NULL;

          SQMATRICES[0] = SqMat[0][0][0]; // A11
          SQMATRICES[1] = SqMat[0][0][3]; // A22

          MATRICES[0] =Mat[0][2]; // B1T
          MATRICES[1] = Mat[0][3]; // B2T

          fesp[0] = FeSpaces[0][0];
          ferhs[0] = FeSpaces[0][0];
          ferhs[1] = FeSpaces[0][0];

          RHSs[0] = Rhs[0][0];
          RHSs[1] = RHSs[0] + N_U;

          aux = new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);

          Assemble2DSlipBC(N_FESpaces, fesp,
                           N_SquareMatrices, SQMATRICES,
                           N_RectMatrices, MATRICES,
                           N_Rhs, RHSs, ferhs,
                           DiscreteForm,
                           BoundaryConditions,
                           BoundValues,
                           aux, FeFunct[0][0],
                           FeFunct[1][0]);
          delete aux;

        } // if (TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION >= 1)


       //======================================================================
       // end of assemble new matrix due to nonlinearity
       //======================================================================

       // build stiffness matrix for next nonlinear iteration step
       // stiffness matrix (left upper block) is stored on
       // M11, (M12, M21, M22)
       // M = M +  tau*TDatabase::TimeDB->THETA1 A
          MatAdd(SqMat[0][0][4], SqMat[0][0][0],
                 tau*TDatabase::TimeDB->THETA1);
          MatAdd(SqMat[0][0][5], SqMat[0][0][1],
                 tau*TDatabase::TimeDB->THETA1);
          MatAdd(SqMat[0][0][6], SqMat[0][0][2],
                 tau*TDatabase::TimeDB->THETA1);
          MatAdd(SqMat[0][0][7], SqMat[0][0][3],
                 tau*TDatabase::TimeDB->THETA1);

          // set current factor of steady state matrix
          gamma = tau*TDatabase::TimeDB->THETA1;

  } //  for(j=0;j<Max_It;j++)
  
#endif //  __FLUIDFLOW__  


#ifdef __WITHSURFACTANT__    
#ifndef __MARANGONISTRESSTEST__  
#ifdef __MOVINGMESH__
     //======================================================================
     // Grid velocity
     //======================================================================
      // inner domain velocity
      GetGridVelocity(Entries_P1, gridsol[1], gridrhs[1],
                      GridKCol_P1, GridRowPtr_P1,
                      GridPos_P1, AuxGridPos_P1,
                      VeloVect[0][0], tau,
                      VeloVect[1][1], GlobalCell_Index[0]);

//  outer domain velocity
      GetGridVelo_outer(Entries_P2, gridsol[2], gridrhs[2],
                        GridKCol_P2, GridRowPtr_P2,
                        GridPos_P2, AuxGridPos_P2,
                        VeloVect[0][0], tau,
                        VeloVect[1][2], GlobalCell_Index[1]);



      length = BeginIndex[2][0][1] -  BeginIndex[2][0][0];
      for(i=0;i<N_Cells_P2;i++)
       {
        DOF = GlobalNumbers[2][0] + BeginIndex[2][0][GlobalCell_Index[1][i]];
        DOF_P2 = GlobalNumbers[2][2] + BeginIndex[2][2][i];

        for(jj=0;jj<length;jj++)
         {
          k = DOF[jj];
          l1 = DOF_P2[jj];
          Sol[2][0][k] = Sol[2][2][l1];
          Sol[2][0][k+N_DOFs[2][0]] = Sol[2][2][l1+N_DOFs[2][2]];
         }
        }

      for(i=0;i<N_Cells_P1;i++)
       {
        DOF = GlobalNumbers[2][0] + BeginIndex[2][0][GlobalCell_Index[0][i]];
        DOF_P1 = GlobalNumbers[2][1] + BeginIndex[2][1][i];

        for(jj=0;jj<length;jj++)
         {
          k = DOF[jj];
          l1 = DOF_P1[jj];
          Sol[2][0][k] = Sol[2][1][l1];
          Sol[2][0][k+N_DOFs[2][0]] = Sol[2][1][l1+N_DOFs[2][1]];
         }
        }
#endif

// update the velocity in outer phase
#ifdef __FLUIDFLOW__
      length = BeginIndex[0][2][1] -  BeginIndex[0][2][0];
      for(i=0;i<N_Cells_P2;i++)
       {
        DOF = GlobalNumbers[0][0] + BeginIndex[0][0][GlobalCell_Index[1][i]];
        DOF_P2 = GlobalNumbers[0][2] + BeginIndex[0][2][i];

        for(jj=0;jj<length;jj++)
         {
          k = DOF[jj];
          l1 = DOF_P2[jj];
          Sol[0][2][l1] = Sol[0][0][k];
          Sol[0][2][l1+N_DOFs[0][2]] = Sol[0][0][k+N_DOFs[0][0]];
         }
        }
#else
     //=====================================================================
      OutPut(endl << "CURRENT TIME: ");
      OutPut(TDatabase::TimeDB->CURRENTTIME << endl);
     //=====================================================================
#endif
   
   //======================================================================
   // soluble surfactants --- begin
   //======================================================================
   // save the old time step solution
   memcpy(Csol_old, Sol[3][2], N_DOFs[3][2]*SizeOfDouble);
   memcpy(Isol_old, Sol[N_FESpaces_All][0], N_DOFs[N_FESpaces_All][0]*SizeOfDouble);
   
    
   for(j=0;j<Max_It_scalar;j++)
   {
    memcpy(Csol_nonlinearstep, Sol[3][2], N_DOFs[3][2]*SizeOfDouble);

    //======================================================================
    // surfactant in outer phase - begin
    //======================================================================

    // assembling matrices
    SQMATRICES_SURFACT[0] = SqMat[3][2][0]; // A
    SQMATRICES_SURFACT[0]->Reset();
    SQMATRICES_SURFACT[1] = SqMat[3][2][1]; // M
    SQMATRICES_SURFACT[1]->Reset();

    fesp[0] = FeSpaces[3][2];  // surfactant space
    fesp[1] = FeSpaces[0][2];  // velocity space
    fesp[2] = FeSpaces[2][2];  // mesh velocity space

    ferhs[0] = FeSpaces[3][2]; // surfactant space for rhs

    fefct[0] = FeFunct[0][2]; // u1
    fefct[1] = FeFunct[1][2]; // u2
    fefct[2] = FeFunct[3][2]; // w1
    fefct[3] = FeFunct[4][2]; // w2

    CRHSs[0] =  Rhs[3][2];

    memset(CRHSs[0], 0, N_DOFs[3][2]*SizeOfDouble);

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
                DiscreteFormSurfact_OutPhase,
                SurfactBoundaryConditions,
                SurfactBoundValues,
                aux);
   delete aux;

     switch(surf_couple_var)
     {
      case 1: // explicit j=0 only, no iteration
       Surfact2D_InterfaceInt(SqMat[3][2][0], CRHSs[0], SurfactBoundCondition, FeFunct[5][2],
                          IFaceFeFunct[0], FeSpaces[3][0], N_List[0], N_List[1], GammaXmaxVal);
       memcpy(CRhs_old, CRHSs[0], N_DOFs[3][2]*SizeOfDouble);
      break;

      case 2: // implicit but explicit in fixed pt iteration step
      case 3:
       Surfact2D_InterfaceInt(SqMat[3][2][0], CRHSs[0], SurfactBoundCondition, FeFunct[5][2],
                          IFaceFeFunct[0], FeSpaces[3][0], N_List[0], N_List[1], GammaXmaxVal);
       if(j==0)
        memcpy(CRhs_old, CRHSs[0], N_DOFs[3][2]*SizeOfDouble);
      break;

      case 4: // implicit but explicit in fixed pt iterartion step
      case 5: // fully implicit
       Surfact2D_InterfaceInt_Implicit(SqMat[3][2][0], CRHSs[0], SurfactBoundCondition, FeFunct[5][2],
                          IFaceFeFunct[0], FeSpaces[3][0], N_List[0], N_List[1], GammaXmaxVal);
      if(j==0)
       memcpy(CRhs_old, CRHSs[0], N_DOFs[3][2]*SizeOfDouble);
      break;

      default:
       OutPut("error in selecting linerizatin type for coupled surfaact eqns " << endl);
       exit(1);
     }


// working rhs for surfactant
   memset(C_B, 0, N_DOFs[3][2]*SizeOfDouble);

   N_CActive = FeSpaces[3][2]->GetActiveBound();

   Daxpy(N_CActive, tau, CRHSs[0], C_B);
//    Daxpy(N_CActive, tau*TDatabase::TimeDB->THETA3, CRhs_old, C_B);
//    Daxpy(N_CActive, tau*TDatabase::TimeDB->THETA4, CRHSs[0], C_B);


   MatAdd(SqMat[3][2][1], SqMat[3][2][0], -tau*TDatabase::TimeDB->THETA2);
   gamma = -tau*TDatabase::TimeDB->THETA2;

   memset(C_defect, 0, N_DOFs[3][2]*SizeOfDouble);
   MatVectActive(SqMat[3][2][1], Csol_old, C_defect);
   Daxpy(N_CActive, 1, C_defect, C_B);

   // set Dirichlet values
   memcpy(C_B+N_CActive, Rhs[3][2]+N_CActive,
          (N_DOFs[3][2] - N_CActive)*SizeOfDouble);

   // copy Dirichlet values from rhs into Sol[6][mg_level-1]
   memcpy(Sol[3][2]+N_CActive,  Rhs[3][2]+N_CActive,
          (N_DOFs[3][2] - N_CActive)*SizeOfDouble);

   MatAdd(SqMat[3][2][1], SqMat[3][2][0], -gamma + tau*TDatabase::TimeDB->THETA1);


   // check the convergence of the fixed point linerization
   memset(C_defect, 0, N_DOFs[3][2]*SizeOfDouble);

   ScalarDefect(SqMat[3][2][1], Sol[3][2], C_B, C_defect, residual_scalar);

   OutPut("Scalar nonlinear step " << setw(3) << j);
   OutPut(setw(14) << residual_scalar); // sqrt of residual_scalar is alread done in ScalarDefect

   if (j>0)
    { OutPut(setw(14) <<  residual_scalar/oldresidual_scalar << endl); }
   else
    { OutPut(endl); }

   oldresidual_scalar = residual_scalar;

   if((((residual_scalar<=limit_scalar)||(j==Max_It-1))) && (j>=TDatabase::ParamDB->SC_MINIT))
    {  break; }

   DirectSolver(SqMat[3][2][1], C_B, Sol[3][2]);

    //======================================================================
    //  surfactant in outer phase - end
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

    fesp[0] = FeSpaces[0][0];  // velospace space 
    fesp[1] = FeSpaces[3][2];  // outer surfactant space 

    IFacefesp[0] = IFaceFeSpaces[0]; // Interface surfactant space
    IFaceferhs[0] = IFaceFeSpaces[0]; // Interface surfactant space for rhs

    // no mesh velocity since interface is moving in an Lagrangian manner
    fefct[0] = FeFunct[0][0]; // ur
    fefct[1] = FeFunct[1][0]; // uz

    fefct[2] = FeFunct[5][2]; // // surfactant in the outer phase

    SRHSs[0] =  Rhs[N_FESpaces_All][0];

    memset(SRHSs[0], 0, N_DOFs[N_FESpaces_All][0]*SizeOfDouble);


    if( (surf_couple_var!=1) && j==0)
     {
      AssembleSurf1D_SolubleSurfact_Implicit(N_FESpaces, fesp, fefct, N_FESpaces_low,
                   IFacefesp, IFaceFeFunct[0], N_SquareMatrices, SQMATRICES_IFace, N_Rhs, SRHSs, 
                   IFaceferhs, N_List[0], N_List[1], Csol_old);

      memcpy(IRhs_old, SRHSs[0], N_DOFs[N_FESpaces_All][0]*SizeOfDouble);

      SQMATRICES_IFace[0]->Reset();
      SQMATRICES_IFace[1]->Reset();
      memset(SRHSs[0], 0, N_DOFs[N_FESpaces_All][0]*SizeOfDouble);
     }
   

     switch(surf_couple_var)
     {
      case 1: // fimplicit but explicit in fixed pt iterartion step
       AssembleSurf1D_SolubleSurfact(N_FESpaces, fesp, fefct, N_FESpaces_low,
                   IFacefesp, IFaceFeFunct[0], N_SquareMatrices, SQMATRICES_IFace, N_Rhs, SRHSs, 
                   IFaceferhs, N_List[0], N_List[1], Csol_old);

       memcpy(IRhs_old, SRHSs[0], N_DOFs[N_FESpaces_All][0]*SizeOfDouble);
      break;

      case 2: // fimplicit but explicit in fixed pt iterartion step
       AssembleSurf1D_SolubleSurfact(N_FESpaces, fesp, fefct, N_FESpaces_low,
                   IFacefesp, IFaceFeFunct[0], N_SquareMatrices, SQMATRICES_IFace, N_Rhs, SRHSs, 
                   IFaceferhs, N_List[0], N_List[1], Csol_nonlinearstep);

      break;

      case 3: // fimplicit but explicit in fixed pt iterartion step
       AssembleSurf1D_SolubleSurfact(N_FESpaces, fesp, fefct, N_FESpaces_low,
                   IFacefesp, IFaceFeFunct[0], N_SquareMatrices, SQMATRICES_IFace, N_Rhs, SRHSs, 
                   IFaceferhs, N_List[0], N_List[1], Sol[3][2]);

      break;

      case 4: // implicit but explicit in fixed pt iterartion step
       AssembleSurf1D_SolubleSurfact_Implicit(N_FESpaces, fesp, fefct, N_FESpaces_low,
                   IFacefesp, IFaceFeFunct[0], N_SquareMatrices, SQMATRICES_IFace, N_Rhs, SRHSs, 
                   IFaceferhs, N_List[0], N_List[1], Csol_nonlinearstep);
      break;

      case 5: // fully implicit
       AssembleSurf1D_SolubleSurfact_Implicit(N_FESpaces, fesp, fefct, N_FESpaces_low,
                     IFacefesp, IFaceFeFunct[0], N_SquareMatrices, SQMATRICES_IFace, N_Rhs, SRHSs, 
                     IFaceferhs, N_List[0], N_List[1], Sol[3][2]);
      break;
      default:
       OutPut("error in selecting linerizatin type for coupled surfaact eqns " << endl);
       exit(1);
     }
        
     
    // working array for srhs is I_B, initialize I_B
    memset(I_B, 0,  N_DOFs[N_FESpaces_All][0]*SizeOfDouble);
    memset(I_defect, 0, N_DOFs[N_FESpaces_All][0]*SizeOfDouble);

    N_IActive = IFaceFeSpaces[0]->GetActiveBound();
    
    Daxpy(N_IActive, tau, Rhs[N_FESpaces_All][0], I_B);
//     Daxpy(N_IActive, tau*TDatabase::TimeDB->THETA3, IRhs_old, I_B);
//     Daxpy(N_IActive, tau*TDatabase::TimeDB->THETA4, Rhs[N_FESpaces_All][0], I_B);

    MatAdd(SqMat_IFace[1], SqMat_IFace[0], -tau*TDatabase::TimeDB->THETA2);
    gamma = -tau*TDatabase::TimeDB->THETA2;

    memset(I_defect, 0, N_DOFs[N_FESpaces_All][0]*SizeOfDouble);
    MatVectActive(SqMat_IFace[1], Isol_old, I_defect);
    Daxpy(N_IActive, 1, I_defect, I_B);

   // set Dirichlet values
   // Rhs[N_FESpaces4][0] still available from assembling
   memcpy(I_B+N_IActive, Rhs[N_FESpaces_All][0]+N_IActive, 
          (N_DOFs[N_FESpaces_All][0] - N_IActive)*SizeOfDouble);

    // copy Dirichlet values from rhs into Sol
   memcpy(Sol[N_FESpaces_All][0]+N_IActive, Rhs[N_FESpaces_All][0]+N_IActive, 
          (N_DOFs[N_FESpaces_All][0] - N_IActive)*SizeOfDouble);

   //=====================================================================
   // assembling of system matrix
   //========================================================================
   MatAdd(SqMat_IFace[1], SqMat_IFace[0], -gamma + tau*TDatabase::TimeDB->THETA1);
    
   //solve the system
   DirectSolver(SqMat_IFace[1], I_B, Sol[N_FESpaces_All][0]);

   
   //======================================================================
   // surfactant on interphase phase - end
   //======================================================================
  } //  for(j=0;j<Max_It;j++)

#endif // Solvesurfactant
#endif //withsurfactant

  //======================================================================
  // soluble surfactants --- end
  //======================================================================
#ifdef __MOVINGMESH__  
    // move grid by known boundary velocity
    RefGridPos_P1->GridToData();
    iso_update = 0;
    MoveGrid_2Phase(Entries_P2, gridsol[2], gridrhs[2],
                    GridKCol_P2, GridRowPtr_P2,
                    GridPos_P2, AuxGridPos_P2,
                    VeloVect[0][0], tau,
                    GlobalCell_Index[1], iso_update);

    RefGridPos_P1->DataToGrid();

    iso_update = 1;    
    MoveGrid_2Phase(Entries_P1, gridsol[1], gridrhs[1],
                    GridKCol_P1, GridRowPtr_P1,
                    GridPos_P1, AuxGridPos_P1,
                    VeloVect[0][0], tau,
                    GlobalCell_Index[0], iso_update);  

   GridPos_P1->GridToData();
   memcpy(refpos_P1, pos_P1, 2*N_DOFs[2][1]*SizeOfDouble); 

   GridPos_P2->GridToData();
   memcpy(refpos_P2, pos_P2, 2*N_DOFs[2][2]*SizeOfDouble); 
   
   GridPos->GridToData();
   memcpy(refpos, pos, 2*N_DOFs[2][0]*SizeOfDouble); 
   
// updating the boundary parameter
// only axial boundary components will change
   MovBoundVert[0][0]->GetCoords(x, y);
   MovBoundVert[1][0]->GetCoords(tx, ty);
   UpdateBound[0]->SetParams(Xi[0], y, Xi[1]-Xi[0],ty-y);

   MovBoundVert[1][0]->GetCoords(x, y);
   MovBoundVert[2][0]->GetCoords(tx, ty);
   UpdateBound[1]->SetParams(Xi[1], y, Xi[2]-Xi[1],ty-y);

   MovBoundVert[5][0]->GetCoords(x, y);
   MovBoundVert[0][0]->GetCoords(tx, ty);
   UpdateBound[5]->SetParams(Xi[5], y, Xi[0]-Xi[5],ty-y);

   // updating all moving boundary vertices
    for(k=0;k<N_MovVert[0];k++)
      {
       MovBoundVert[0][k]->GetCoords(x, y);
       //        cout<< " SLPX0 " << Xi[0]<<" SLPY " << y<<endl;
       MovBoundVert[0][k]->SetCoords(Xi[0], y);
       }

    for(k=0;k<N_MovVert[1];k++)
      {
       MovBoundVert[1][k]->GetCoords(x, y);
       //        cout<< " SLPX1 " << Xi[1]<<" SLPY " << y<<endl;
       MovBoundVert[1][k]->SetCoords(Xi[1], y);
      }

    // BD comp 2 and 4 are Diriclet, i.e., no movement
    for(k=0;k<N_MovVert[3];k++)
      {
       MovBoundVert[3][k]->GetCoords(x, y);
       //        cout<< " SLPX " << Xi[3] <<" SLPY " << y<<endl;
       MovBoundVert[3][k]->SetCoords(Xi[3], y );
       }

    for(k=0;k<N_MovVert[5];k++)
      {
       MovBoundVert[5][k]->GetCoords(x, y);
       //        cout<< " SLPX1 " << Xi[5]<<" SLPY " << y<<endl;
       MovBoundVert[5][k]->SetCoords(Xi[5], y);
       }

    for(k=0;k<N_MovVert[0];k++)
    {
     if(k==N_MovVert[0]-1)
      Slip_Joint[0][k]->UpdateParameters(MovBoundVert[0][k], MovBoundVert[1][0]);
     else
      Slip_Joint[0][k]->UpdateParameters(MovBoundVert[0][k], MovBoundVert[0][k+1]);
    }

    for(k=0;k<N_MovVert[1];k++)
    {
     if(k==N_MovVert[1]-1)
      Slip_Joint[1][k]->UpdateParameters(MovBoundVert[1][k], MovBoundVert[2][0]);
     else
      Slip_Joint[1][k]->UpdateParameters(MovBoundVert[1][k], MovBoundVert[1][k+1]);
    }


    for(k=0;k<N_MovVert[3];k++)
    {
     if(k==N_MovVert[3]-1)
      Slip_Joint[3][k]->UpdateParameters(MovBoundVert[3][k], MovBoundVert[4][0]);
     else
      Slip_Joint[3][k]->UpdateParameters(MovBoundVert[3][k], MovBoundVert[3][k+1]);
    }

    for(k=0;k<N_MovVert[5];k++)
    {
     if(k==N_MovVert[5]-1)
      Slip_Joint[5][k]->UpdateParameters(MovBoundVert[5][k], MovBoundVert[0][0]);
     else
      Slip_Joint[5][k]->UpdateParameters(MovBoundVert[5][k], MovBoundVert[5][k+1]);
    }

   if((l==0) && ((m % 1) == 0))
    {
     Getcellangle(FeSpaces[0][0], Angle);
     OutPut( "MinAngle : "<< Angle[0]<< "  MaxAngle : "<<Angle[1]<< endl);
    }  

   if((Angle[0]<10.0) ||(Angle[1]>165.0))
    {
      
     FreePts[0] = new double[N_MovVert[6]];
     FreePts[1] = new double[N_MovVert[6]];

// =============================================================================================
// reparametrize the interface points before remeshing

     fhtot = 0.;
     fhmin = 100;
     fhmax = 0.0;

     for(k=0;k<N_MovVert[6];k++)
     {
      MovBoundVert[6][k]->GetCoords(x1, y1);
      MovBoundVert[6][k]->GetCoords(FreePts[0][k], FreePts[1][k]);
    
      if(k==N_MovVert[6]-1)
       MovBoundVert[0][0]->GetCoords(x2, y2);
      else
       MovBoundVert[6][k+1]->GetCoords(x2, y2);

      fh = sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));

      fhtot +=fh;
      if (fh < fhmin) fhmin = fh;
      if (fh > fhmax) fhmax = fh;
     }

    fhtot /= (double)N_MovVert[6];
    fhlimit = fhtot*fhtot;
    
//     if( ((fhmin < (fhtot - fhlimit) ) || (fhmax > (fhtot + fhlimit)))  && (!remeshed)  )   
//      {
      OutPut("Interface Reparam:  "<<img<<' ' << fhlimit <<' '<< 3*fhtot/2.0 <<' '<< fhmin <<' '<<fhmax<< endl);
      
      N_E = N_MovVert[6];
      ReParam_axial3D_Data(N_E, Coll_Cells[0],  N_List[3], N_List[2], VeloVect[0][0], 
                           FeFunct[6][0], Intpol_Coord,  Intpol_VeloValues, Intpol_Values, h_interface, FreePts);
      reparam = TRUE;
      cout<< "Old N_E " << N_MovVert[6] <<"New N_E " << N_E <<endl;     
      N_MovVert[6] = N_E;
//      }//if( ((fhmin < (fhtot - fhlimit) ) || (fhmax > (
//     else
//     {
//      N_E = N_MovVert[6];
//     
//      N_AllIntVertices = N_E + 1 + (ORDER-1) *N_E;
//      cout << " N_AllIntVertices " << N_AllIntVertices << endl;
// 
//      Intpol_Coord = new double[2*N_AllIntVertices];
//      Intpol_Values = new double[N_AllIntVertices];
//      
//      SurfValuesForInterpol(N_E, N_List[2], N_List[3], FeFunct[6][0], Intpol_Coord,  Intpol_Values);
//     }
// ============================================================================================
#ifdef __WITHSURFACTANT__        
    // bulk phase surfactants
    length = BeginIndex[3][0][1] -  BeginIndex[3][0][0];
    for(i=0;i<N_Cells_P2;i++)
     {
      DOF = GlobalNumbers[3][0] + BeginIndex[3][0][GlobalCell_Index[1][i]];
      DOF_P2 = GlobalNumbers[3][2] + BeginIndex[3][2][i];

      for(jj=0;jj<length;jj++)
       {
        k = DOF[jj];
        l1 = DOF_P2[jj];
        Sol[3][0][k] = Sol[3][2][l1];
       }
      }
    // interface surfactant  
    MapSurfToDomain(IFaceFeFunct[0],  FeFunct[6][0], N_List[0], N_List[1]);
#endif

     os.seekp(std::ios::beg);
     if(N_Remesh<10) os << "VTK/"<<VtkBaseName<<"_remesh.0000"<<2*N_Remesh + 1<<".vtk" << ends;
     else if(N_Remesh<100) os << "VTK/"<<VtkBaseName<<"_remesh.000"<<2*N_Remesh+ 1<<".vtk" << ends;
     else if(N_Remesh<1000) os << "VTK/"<<VtkBaseName<<"_remesh.00"<<2*N_Remesh+ 1<<".vtk" << ends;
     else if(N_Remesh<10000) os << "VTK/"<<VtkBaseName<<"_remesh.0"<<2*N_Remesh+ 1<<".vtk" << ends;
     else  os << "VTK/"<<VtkBaseName<<"_remesh."<<2*N_Remesh+ 1<<".vtk" << ends;
     Output->WriteVtk(os.str().c_str());

// ===========================================================================================

    t1 = GetTime();

    Remesh2D_2PhaseAxial3D(Domain, Coll_Cells, Coll_Multi, MovBoundVert, N_MovVert,
                     N_FESpaces_All, FeSpaces, N_DOFs, GlobalCell_Index, SqMat, Mat, Sol,
                     Rhs, VeloVect, FeFunct, Slip_Joint, SqrStruct, Struct, FreePts,
                     IFaceDomain, IFaceFeSpaces, N_List, IFaceFeFunct, SqMat_IFace,
                     IFaceStruct, FE1D_List, GlobalNumbers, BeginIndex, Bound_DOFs,
                     Intpol_Coord,  Intpol_Values, Intpol_VeloValues);
    
    Getcellangle(FeSpaces[0][0], Angle);
    remeshed = TRUE;
    
    if(reparam == TRUE)
     {
      delete [] Intpol_Coord;
      delete [] Intpol_Values;
      delete [] Intpol_VeloValues;
    
      reparam =FALSE;
     }
     
    delete [] FreePts[0];
    delete [] FreePts[1];    
    
    MapDomainToSurf(FeFunct[6][0], IFaceFeFunct[0], N_List[0], N_List[1]);
    memset(Sol[4][0], 0,  N_DOFs[4][0]*SizeOfDouble);

    OutPut( "Remeshing and Interpolation were done"<<endl);

    N_U = N_DOFs[0][0];
    N_P = N_DOFs[1][0];
    N_Unknowns = 2*N_U+N_P;
    
    Coll_P1 = FeSpaces[0][1]->GetCollection();
    N_Cells_P1 = Coll_P1->GetN_Cells();
    Coll_P2 = FeSpaces[0][2]->GetCollection();
    N_Cells_P2 = Coll_P2->GetN_Cells();

    delete [] oldsol;
    delete [] oldrhs;
    delete [] B;
    delete [] defect;
    delete [] outputsol;
   
    oldsol =  new double[N_Unknowns];
    oldrhs =  new double[N_Unknowns];
    B =  new double[N_Unknowns];
    defect =  new double[N_Unknowns];
    outputsol = new double[N_DOFs[2][0]];

    delete [] gridsol[0];
    delete [] gridsol[1];
    delete [] gridsol[2];
    delete [] gridrhs[0];
    delete [] gridrhs[1];
    delete [] gridrhs[2];
    gridsol[0] =new double[2*N_DOFs[2][0]];
    gridsol[1] =new double[2*N_DOFs[2][1]];
    gridsol[2] =new double[2*N_DOFs[2][2]];
    gridrhs[0] =new double[2*N_DOFs[2][0]];
    gridrhs[1] =new double[2*N_DOFs[2][1]];
    gridrhs[2] =new double[2*N_DOFs[2][2]];
    memset(gridsol[0], 0, 2*N_DOFs[2][0]*SizeOfDouble);
    memset(gridsol[1], 0, 2*N_DOFs[2][1]*SizeOfDouble);
    memset(gridsol[2], 0, 2*N_DOFs[2][2]*SizeOfDouble);
    memset(gridrhs[0], 0, 2*N_DOFs[2][0]*SizeOfDouble);
    memset(gridrhs[1], 0, 2*N_DOFs[2][1]*SizeOfDouble);
    memset(gridrhs[2], 0, 2*N_DOFs[2][2]*SizeOfDouble);
    memset(outputsol, 0, N_DOFs[2][0]*SizeOfDouble);

    delete [] C_defect;
    delete [] Csol_old;
    delete [] C_B;
    delete [] CRhs_old;
    delete [] Csol_nonlinearstep;   

    C_defect = new double[N_DOFs[3][2]];
    Csol_old = new double[N_DOFs[3][2]];
    C_B = new double[N_DOFs[3][2]];
    CRhs_old = new double[N_DOFs[3][2]];
    Csol_nonlinearstep = new double[N_DOFs[3][2]];

    delete [] I_defect;
    delete [] Isol_old;
    delete [] I_B;
    delete [] IRhs_old;

    I_defect = new double[N_DOFs[N_FESpaces_All][0]];
    Isol_old = new double[N_DOFs[N_FESpaces_All][0]];
    I_B = new double[N_DOFs[N_FESpaces_All][0]];
    IRhs_old = new double[N_DOFs[N_FESpaces_All][0]];
 
    delete [] refpos;
    delete [] auxpos;
    delete [] pos;
    delete [] tmp;  
 
    delete [] refpos_P1;
    delete [] auxpos_P1;
    delete [] pos_P1;
     
    delete [] refpos_P2;
    delete [] auxpos_P2;
    delete [] pos_P2;
 
    refpos = new double[2*N_DOFs[2][0]];
    auxpos = new double[2*N_DOFs[2][0]];
    pos = new double[2*N_DOFs[2][0]];
    tmp = new double[2*N_DOFs[2][0]];

    memset(tmp, 0, 2*N_DOFs[2][0]*SizeOfDouble);
 
    refpos_P1 = new double[2*N_DOFs[2][1]];
    auxpos_P1 = new double[2*N_DOFs[2][1]];
    pos_P1 = new double[2*N_DOFs[2][1]];

    refpos_P2 = new double[2*N_DOFs[2][2]];
    auxpos_P2 = new double[2*N_DOFs[2][2]];
    pos_P2 = new double[2*N_DOFs[2][2]];

    RefGridPos = new TFEVectFunct2D(FeSpaces[2][0], refposString, gridString,
                                  refpos, N_DOFs[2][0], 2);

    AuxGridPos = new TFEVectFunct2D(FeSpaces[2][0], auxposString, gridString,
                                  auxpos, N_DOFs[2][0], 2);
    GridPos = new TFEVectFunct2D(FeSpaces[2][0], posString, gridString,
                                 pos, N_DOFs[2][0], 2);


    RefGridPos_P1 = new TFEVectFunct2D(FeSpaces[2][1], refposString, refposString,
                                  refpos_P1, N_DOFs[2][1], 2);
    AuxGridPos_P1 = new TFEVectFunct2D(FeSpaces[2][1], auxposString, auxposString,
                                  auxpos_P1, N_DOFs[2][1], 2);
    GridPos_P1 = new TFEVectFunct2D(FeSpaces[2][1], posString, posString,
                                  pos_P1, N_DOFs[2][1], 2);


    RefGridPos_P2 = new TFEVectFunct2D(FeSpaces[2][2], refposString, gridString,
                                  refpos_P2, N_DOFs[2][2], 2);
    AuxGridPos_P2 = new TFEVectFunct2D(FeSpaces[2][2], auxposString, gridString,
                                  auxpos_P2, N_DOFs[2][2], 2);
    GridPos_P2 = new TFEVectFunct2D(FeSpaces[2][2], posString, gridString,
                                  pos_P2, N_DOFs[2][2], 2);
 
    OutputFeunction = new TFEFunction2D(FeSpaces[2][0], gridString,  gridString, outputsol, N_DOFs[2][0]);
  
    RefGridPos->GridToData();
    GridPos->GridToData();
    AuxGridPos->GridToData();

    RefGridPos_P1->GridToData();
    GridPos_P1->GridToData();
    AuxGridPos_P1->GridToData();

    RefGridPos_P2->GridToData();
    GridPos_P2->GridToData();
    AuxGridPos_P2->GridToData();

    //==================================================================================
    delete Output;

 // prepare output (maxn_fespaces,  maxn_scalar,  maxn_vect, maxn_parameters, domain)
    Output = new TOutput2D(4, 5, 1, 2, Domain);
    if(FluidFlow)
     {
      Output->AddFEVectFunct(VeloVect[0][0]);
      Output->AddFEFunction(FeFunct[2][0]);
     }

// mesh velocity
//    Output->AddFEVectFunct(VeloVect[1][0]);

 /* surfactant  */
 #ifdef __WITHSURFACTANT__ 
//    if(TDatabase::ParamDB->Compute_Energy)
     Output->AddFEFunction(FeFunct[5][0]);

 /* interface surfactant  */
     Output->AddFEFunction(FeFunct[6][0]);
#endif     
 /* rise velo in bubble and gridvelo in outer  */
     Output->AddFEFunction(OutputFeunction);
   
     os.seekp(std::ios::beg);
     Output->AddParameter(TDatabase::TimeDB->CURRENTTIME,os.str().c_str());
   //==================================================================================
   length = BeginIndex[3][0][1] -  BeginIndex[3][0][0];
   for(i=0;i<N_Cells_P2;i++)
    {
     DOF = GlobalNumbers[3][0] + BeginIndex[3][0][GlobalCell_Index[1][i]];
     DOF_P2 = GlobalNumbers[3][2] + BeginIndex[3][2][i];

     for(jj=0;jj<length;jj++)
      {
       k = DOF[jj];
       l1 = DOF_P2[jj];
       Sol[3][0][k] = Sol[3][2][l1];
      }
     }
     MapSurfToDomain(IFaceFeFunct[0],  FeFunct[6][0], N_List[0], N_List[1]);
       os.seekp(std::ios::beg);
       if(N_Remesh<10) os << "VTK/"<<VtkBaseName<<"_remesh.0000"<<2*N_Remesh+2<<".vtk" << ends;
       else if(N_Remesh<100) os << "VTK/"<<VtkBaseName<<"_remesh.000"<<2*N_Remesh+2<<".vtk" << ends;
       else if(N_Remesh<1000) os << "VTK/"<<VtkBaseName<<"_remesh.00"<<2*N_Remesh+2<<".vtk" << ends;
       else if(N_Remesh<10000) os << "VTK/"<<VtkBaseName<<"_remesh.0"<<2*N_Remesh+2<<".vtk" << ends;
       else  os << "VTK/"<<VtkBaseName<<"_remesh."<<2*N_Remesh+2<<".vtk" << ends;
       Output->WriteVtk(os.str().c_str());
       N_Remesh++;
       
//        exit(0);
    }
 
  //==================================================================================
  // Assembeling the grid matrix
  //==================================================================================
  t1 = GetTime();

  fesp[0] = FeSpaces[2][0];
  SQMATRICES_GRID[0] = SqMat[2][0][0];
  SQMATRICES_GRID[0]->Reset();
  SQMATRICES_GRID[1] = SqMat[2][0][1];
  SQMATRICES_GRID[1]->Reset();
  SQMATRICES_GRID[2] = SqMat[2][0][2];
  SQMATRICES_GRID[2]->Reset();
  SQMATRICES_GRID[3] = SqMat[2][0][3];
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

  Entries[0] = SqMat[2][0][0]->GetEntries();
  Entries[1] = SqMat[2][0][1]->GetEntries();
  Entries[2] = SqMat[2][0][2]->GetEntries();
  Entries[3] = SqMat[2][0][3]->GetEntries();
  GridKCol = SqrStruct[2][0]->GetKCol();
  GridRowPtr = SqrStruct[2][0]->GetRowPtr();
  memset(Entries[1] + GridRowPtr[N_DOFs[2][0]-Bound_DOFs[2][0]], 0,
         Bound_DOFs[2][0]*SizeOfDouble);
  memset(Entries[2] + GridRowPtr[N_DOFs[2][0]-Bound_DOFs[2][0]], 0,
         Bound_DOFs[2][0]*SizeOfDouble);

 t2 = GetTime();

  if(TDatabase::ParamDB->SC_VERBOSE > 1)
  {
    cout << "Grid assembling done"<< endl;
    OutPut("Time for Grid assembling: " << t2-t1 << endl);
  }

  // inner Phase
  t1 = GetTime();

  fesp[0] = FeSpaces[2][1];
  SQMATRICES_GRID[0] = SqMat[2][1][0];
  SQMATRICES_GRID[0]->Reset();
  SQMATRICES_GRID[1] = SqMat[2][1][1];
  SQMATRICES_GRID[1]->Reset();
  SQMATRICES_GRID[2] = SqMat[2][1][2];
  SQMATRICES_GRID[2]->Reset();
  SQMATRICES_GRID[3] = SqMat[2][1][3];
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

  Entries_P1[0] = SqMat[2][1][0]->GetEntries();
  Entries_P1[1] = SqMat[2][1][1]->GetEntries();
  Entries_P1[2] = SqMat[2][1][2]->GetEntries();
  Entries_P1[3] = SqMat[2][1][3]->GetEntries();
  GridKCol_P1 = SqrStruct[2][1]->GetKCol();
  GridRowPtr_P1 = SqrStruct[2][1]->GetRowPtr();
  memset(Entries_P1[1] + GridRowPtr_P1[N_DOFs[2][1]-Bound_DOFs[2][1]], 0,
         Bound_DOFs[2][1]*SizeOfDouble);
  memset(Entries_P1[2] + GridRowPtr_P1[N_DOFs[2][1]-Bound_DOFs[2][1]], 0,
         Bound_DOFs[2][1]*SizeOfDouble);

 t2 = GetTime();

  if(TDatabase::ParamDB->SC_VERBOSE > 1)
  {
    cout << "Internal Phase Grid assembling done"<< endl;
    OutPut("Time for Internal Grid assembling: " << t2-t1 << endl);
  }

/* outer Phase */
  t1 = GetTime();

  fesp[0] = FeSpaces[2][2];
  SQMATRICES_GRID[0] = SqMat[2][2][0];
  SQMATRICES_GRID[0]->Reset();
  SQMATRICES_GRID[1] = SqMat[2][2][1];
  SQMATRICES_GRID[1]->Reset();
  SQMATRICES_GRID[2] = SqMat[2][2][2];
  SQMATRICES_GRID[2]->Reset();
  SQMATRICES_GRID[3] = SqMat[2][2][3];
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

  Entries_P2[0] = SqMat[2][2][0]->GetEntries();
  Entries_P2[1] = SqMat[2][2][1]->GetEntries();
  Entries_P2[2] = SqMat[2][2][2]->GetEntries();
  Entries_P2[3] = SqMat[2][2][3]->GetEntries();
  GridKCol_P2 = SqrStruct[2][2]->GetKCol();
  GridRowPtr_P2 = SqrStruct[2][2]->GetRowPtr();
  memset(Entries_P2[1] + GridRowPtr_P2[N_DOFs[2][2]-Bound_DOFs[2][2]], 0,
         Bound_DOFs[2][2]*SizeOfDouble);
  memset(Entries_P2[2] + GridRowPtr_P2[N_DOFs[2][2]-Bound_DOFs[2][2]], 0,
         Bound_DOFs[2][2]*SizeOfDouble);

  t2 = GetTime();

  if(TDatabase::ParamDB->SC_VERBOSE > 1)
  {
    cout << "Internal Phase Grid assembling done"<< endl;
    OutPut("Time for Internal Grid assembling: " << t2-t1 << endl);
  }

 //==================================================================================
 // End Assembeling the reference grid
 //  Reparametrization - Begin
 //==================================================================================
//   if(!remeshed)
//    {
//      fhtot = 0.;
//      fhmin = 100;
//      fhmax = 0.0;
// 
//      for(k=0;k<N_MovVert[6];k++)
//      {
//       MovBoundVert[6][k]->GetCoords(x1, y1);
// 
//       if(k==N_MovVert[6]-1)
// 	MovBoundVert[0][0]->GetCoords(x2, y2);
//       else
// 	MovBoundVert[6][k+1]->GetCoords(x2, y2);
// 
//       fh = sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
// 
//       fhtot +=fh;
//       if (fh < fhmin) fhmin = fh;
//       if (fh > fhmax) fhmax = fh;
//      }
// 
//     fhtot /= (double)N_MovVert[6];
//     fhlimit = fhtot*fhtot;
//     
//     if( ((fhmin < (fhtot - fhlimit) ) || (fhmax > (fhtot + fhlimit)))  && (!remeshed)  )   
//      {
//       OutPut("Interface Reparam:  "<<img<<' ' << fhlimit <<' '<< 3*fhtot/2.0 <<' '<< fhmin <<' '<<fhmax<< endl);
//       
// // =====================================================================================================      
//        length = BeginIndex[3][0][1] -  BeginIndex[3][0][0];
//        for(i=0;i<N_Cells_P2;i++)
//         {
//          DOF = GlobalNumbers[3][0] + BeginIndex[3][0][GlobalCell_Index[1][i]];
//          DOF_P2 = GlobalNumbers[3][2] + BeginIndex[3][2][i];
// 
//          for(jj=0;jj<length;jj++)
//           {
//            k = DOF[jj];
//            l1 = DOF_P2[jj];
//            Sol[3][0][k] = Sol[3][2][l1];
//           }
//          }
// 
//        MapSurfToDomain(IFaceFeFunct[0], FeFunct[6][0], N_List[0], N_List[1]);
//        
//        os.seekp(std::ios::beg);
//        if(img<10) os <<"VTK/"<< VtkBaseName<<"ReParam.0000"<<ReParam_img<<".vtk" << ends;
//        else if(img<100) os << "VTK/"<<VtkBaseName<<"ReParam.000"<<ReParam_img<<".vtk" << ends;
//        else if(img<1000) os <<"VTK/"<<VtkBaseName<<"ReParam.00"<<ReParam_img<<".vtk" << ends;
//        else if(img<10000) os <<"VTK/"<< VtkBaseName<<"ReParam.0"<<ReParam_img<<".vtk" << ends;
//        else  os << "VTK/"<<VtkBaseName<<"."<<ReParam_img<<"ReParam.vtk" << ends;
//        Output->WriteVtk(os.str().c_str());
//        ReParam_img++;
// // =====================================================================================================
//        
//       MapSurfToDomain(IFaceFeFunct[0], FeFunct[6][0], N_List[0], N_List[1]);
//       
//       RefGridPos_P1->GridToData();
//       RefGridPos_P2->GridToData();
//       ReParam_axial3D_U(N_MovVert[6], Coll_Cells[0],  N_List[3], N_List[2], VeloVect[0][0], FeFunct[6][0]);
//       reparam = TRUE;
//       
//       MapDomainToSurf(FeFunct[6][0], IFaceFeFunct[0], N_List[0], N_List[1]);
//       memset(Sol[4][0], 0,  N_DOFs[4][0]*SizeOfDouble);
//       
// // =====================================================================================================      
//        length = BeginIndex[3][0][1] -  BeginIndex[3][0][0];
//        for(i=0;i<N_Cells_P2;i++)
//         {
//          DOF = GlobalNumbers[3][0] + BeginIndex[3][0][GlobalCell_Index[1][i]];
//          DOF_P2 = GlobalNumbers[3][2] + BeginIndex[3][2][i];
// 
//          for(jj=0;jj<length;jj++)
//           {
//            k = DOF[jj];
//            l1 = DOF_P2[jj];
//            Sol[3][0][k] = Sol[3][2][l1];
//           }
//          }
// 
//        MapSurfToDomain(IFaceFeFunct[0], FeFunct[6][0], N_List[0], N_List[1]);
//        
//        os.seekp(std::ios::beg);
//        if(img<10) os <<"VTK/"<< VtkBaseName<<"ReParam.0000"<<ReParam_img<<".vtk" << ends;
//        else if(img<100) os <<"VTK/"<< VtkBaseName<<"ReParam.000"<<ReParam_img<<".vtk" << ends;
//        else if(img<1000) os << "VTK/"<<VtkBaseName<<"ReParam.00"<<ReParam_img<<".vtk" << ends;
//        else if(img<10000) os << "VTK/"<<VtkBaseName<<"ReParam.0"<<ReParam_img<<".vtk" << ends;
//        else  os << "VTK/"<<VtkBaseName<<"."<<ReParam_img<<"ReParam.vtk" << ends;
//        Output->WriteVtk(os.str().c_str());
//        ReParam_img++;       
// // =====================================================================================================
//     }
 
//    } // if(!remeshed)
 //==================================================================================
 // End Reparametrization 
 //================================================================================== 
 
 
   if(remeshed || reparam)
    { 
 
     if(remeshed)
      {
       // mesh velo is zero
       memset(Sol[2][0], 0, 2*N_DOFs[2][0]*SizeOfDouble); 
       memset(Sol[2][1], 0, 2*N_DOFs[2][1]*SizeOfDouble);      
       memset(Sol[2][2], 0, 2*N_DOFs[2][2]*SizeOfDouble);        
       memcpy(oldsol, Sol[0][0], SizeOfDouble*N_Unknowns); 
      }
      else
      {
       GridPos_P1->GridToData();
       GridPos_P2->GridToData();

       memset(Sol[2][0], 0, 2*N_DOFs[2][0]*SizeOfDouble);
     
       // Phase 1
       N_G = N_DOFs[2][1];
       N_BoundaryNodes = Bound_DOFs[2][1];
         
       memcpy(tmp, pos_P1, 2*N_G*SizeOfDouble);
       Daxpy(2*N_G, -1, refpos_P1, tmp);
     
       memset(Rhs[2][1], 0, 2*N_G*SizeOfDouble);
       memset(Sol[2][1], 0, 2*N_G*SizeOfDouble); 
     
       memcpy(Rhs[2][1]+(N_G-N_BoundaryNodes), tmp+(N_G-N_BoundaryNodes), N_BoundaryNodes*SizeOfDouble);
       memcpy(Sol[2][1]+(N_G-N_BoundaryNodes), tmp+(N_G-N_BoundaryNodes), N_BoundaryNodes*SizeOfDouble);
     
       memcpy(Rhs[2][1]+(2*N_G-N_BoundaryNodes), tmp+(2*N_G-N_BoundaryNodes), N_BoundaryNodes*SizeOfDouble); 
       memcpy(Sol[2][1]+(2*N_G-N_BoundaryNodes), tmp+(2*N_G-N_BoundaryNodes), N_BoundaryNodes*SizeOfDouble); 
   
       SolveGridEquation(Entries_P1, Sol[2][1], Rhs[2][1],  GridKCol_P1, GridRowPtr_P1, N_G );
     
       Daxpy(2*N_G, 1, Sol[2][1], refpos_P1);
       Dscal(2*N_G, 1./tau, Sol[2][1]);
     
       //no grid velo on Interface, since we interpolated already, so (u-w)=0
       memset(Sol[2][1]+(N_G-N_BoundaryNodes), 0, N_BoundaryNodes*SizeOfDouble);
       memset(Sol[2][1]+(2*N_G-N_BoundaryNodes), 0, N_BoundaryNodes*SizeOfDouble);
     
//      RefGridPos_P1->DataToGrid(); //nonlinear iteration on old mesh, move to new mesh after iteration
     
       // Phase 2
       N_G = N_DOFs[2][2];
       N_BoundaryNodes = Bound_DOFs[2][2];    
          
       memcpy(tmp, pos_P2, 2*N_G*SizeOfDouble);
       Daxpy(2*N_G, -1, refpos_P2, tmp);     

       memset(Rhs[2][2], 0, 2*N_G*SizeOfDouble);
       memset(Sol[2][2], 0, 2*N_G*SizeOfDouble); 
     
       memcpy(Rhs[2][2]+(N_G-N_BoundaryNodes), tmp+(N_G-N_BoundaryNodes), N_BoundaryNodes*SizeOfDouble);
       memcpy(Sol[2][2]+(N_G-N_BoundaryNodes), tmp+(N_G-N_BoundaryNodes), N_BoundaryNodes*SizeOfDouble);
     
       memcpy(Rhs[2][2]+(2*N_G-N_BoundaryNodes), tmp+(2*N_G-N_BoundaryNodes), N_BoundaryNodes*SizeOfDouble); 
       memcpy(Sol[2][2]+(2*N_G-N_BoundaryNodes), tmp+(2*N_G-N_BoundaryNodes), N_BoundaryNodes*SizeOfDouble); 
     
       SolveGridEquation(Entries_P2, Sol[2][2], Rhs[2][2],  GridKCol_P2, GridRowPtr_P2, N_G );
     
       Daxpy(2*N_G, 1, Sol[2][2], refpos_P2);     
       Dscal(2*N_G, 1./tau, Sol[2][2]);
     
       //no grid velo on Interface, since we interpolated already, so (u-w)=0
       memset(Sol[2][2]+(N_G-N_BoundaryNodes), 0, N_BoundaryNodes*SizeOfDouble);
       memset(Sol[2][2]+(2*N_G-N_BoundaryNodes), 0, N_BoundaryNodes*SizeOfDouble);
     
//      RefGridPos_P2->DataToGrid(); //nonlinear iteration on old mesh, move to new mesh after iteration

      // copy the new grid position and velocity to the full domain
       length = BeginIndex[2][0][1] -  BeginIndex[2][0][0];
       for(i=0;i<N_Cells_P2;i++)
       {
        DOF = GlobalNumbers[2][0] + BeginIndex[2][0][GlobalCell_Index[1][i]];
        DOF_P2 = GlobalNumbers[2][2] + BeginIndex[2][2][i];

        for(jj=0;jj<length;jj++)
         {
          k = DOF[jj];
          l1 = DOF_P2[jj];
          Sol[2][0][k] = Sol[2][2][l1];
          Sol[2][0][k+N_DOFs[2][0]] = Sol[2][2][l1+N_DOFs[2][2]];

          refpos[k] = refpos_P2[l1];
          refpos[k+N_DOFs[2][0]] = refpos_P2[l1+N_DOFs[2][2]];
	 }
        }

      for(i=0;i<N_Cells_P1;i++)
       {
        DOF = GlobalNumbers[2][0] + BeginIndex[2][0][GlobalCell_Index[0][i]];
        DOF_P1 = GlobalNumbers[2][1] + BeginIndex[2][1][i];

        for(jj=0;jj<length;jj++)
         {
          k = DOF[jj];
          l1 = DOF_P1[jj];
          Sol[2][0][k] = Sol[2][1][l1];
          Sol[2][0][k+N_DOFs[2][0]] = Sol[2][1][l1+N_DOFs[2][1]];
  
          refpos[k] = refpos_P1[l1];
          refpos[k+N_DOFs[2][0]] = refpos_P1[l1+N_DOFs[2][1]];  
         }
        }

        OutPut("ReParam non-linear CURRENT TIME: ");
        OutPut(TDatabase::TimeDB->CURRENTTIME << endl); 
      } //  else reparam)
       
       
      // working array for rhs is B, initialize B
      memset(B, 0, N_Unknowns*SizeOfDouble);

       DiscreteForm = DiscreteFormGalerkin;

       SQMATRICES[0] = SqMat[0][0][0]; // A11
       SQMATRICES[1] = SqMat[0][0][1]; // A12
       SQMATRICES[2] = SqMat[0][0][2]; // A21
       SQMATRICES[3] = SqMat[0][0][3]; // A22
       SQMATRICES[4] = SqMat[0][0][4]; // M11
       SQMATRICES[5] = SqMat[0][0][7]; // M22

       SqMat[0][0][5]->Reset();  // M12
       SqMat[0][0][6]->Reset();  // M21

       MATRICES[0] = Mat[0][0]; // B1
       MATRICES[1] = Mat[0][1]; // B2
       MATRICES[2] = Mat[0][2]; // B1T
       MATRICES[3] = Mat[0][3]; // B2T

       SQMATRICES[0]->Reset();
       SQMATRICES[1]->Reset();
       SQMATRICES[2]->Reset();
       SQMATRICES[3]->Reset();
       SQMATRICES[4]->Reset();
       SQMATRICES[5]->Reset();

       MATRICES[0]->Reset();
       MATRICES[1]->Reset();
       MATRICES[2]->Reset();
       MATRICES[3]->Reset();

       N_SquareMatrices = 6;
       N_RectMatrices = 4;

       N_Rhs = 2;
       N_FESpaces = 3;

       fesp[0] = FeSpaces[0][0];
       fesp[1] = FeSpaces[1][0];
       fesp[2] = FeSpaces[2][0];

       fefct[0] = FeFunct[0][0];
       fefct[1] = FeFunct[1][0];
       fefct[2] = FeFunct[3][0];
       fefct[3] = FeFunct[4][0];

       ferhs[0] = FeSpaces[0][0];
       ferhs[1] = FeSpaces[0][0];

       RHSs[0] = Rhs[0][0];
       RHSs[1] = RHSs[0] + N_U;
       RHSs[2] = RHSs[0] + 2*N_U;

       memset(RHSs[0], 0, N_Unknowns*SizeOfDouble);

       aux =  new TAuxParam2D(MovingTNSN_FESpaces_Axial3D, MovingTNSN_Fct_Axial3D,
                           MovingTNSN_ParamFct_Axial3D,
                           MovingTNSN_FEValues_Axial3D,
                           fesp+1, fefct,
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

      SqMat[0][0][8]->Reset(); // F11
      SqMat[0][0][9]->Reset(); // F22

      MapSurfToDomain(IFaceFeFunct[0], FeFunct[6][0], N_List[0], N_List[1]);
      
      FreeSurf_2PhaseSurfAxial3D(SqMat[0][0][8], SqMat[0][0][9],
                                 RHSs[0], RHSs[0]+N_U,
                                 BoundCondition, tau, 0, FeFunct[6][0]);

      // Adding freesurf entries to A11 and A22
      MatAdd(SqMat[0][0][0], SqMat[0][0][8], 1.);
      MatAdd(SqMat[0][0][3], SqMat[0][0][9], 1.);

      N_Active = FeSpaces[0][0]->GetActiveBound();
        // get row in off diagonal matrix where the Dirichlet nodes start
      RowPtr = SqMat[0][0][0]->GetRowPtr();
        // compute number of entries starting from this row to the end
        // of the matrix
      j = RowPtr[N_Active];
      k = RowPtr[N_U]-j;
      // get number of active dof
      // set these entries to zero
      memset(SqMat[0][0][1]->GetEntries()+j, 0, SizeOfDouble*k);
      memset(SqMat[0][0][2]->GetEntries()+j, 0, SizeOfDouble*k);

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

          SQMATRICES[0] = SqMat[0][0][0]; // A11
          SQMATRICES[1] = SqMat[0][0][3]; // A22
          SQMATRICES[2] = SqMat[0][0][1]; // A12
          SQMATRICES[3] = SqMat[0][0][2]; // A21
          SQMATRICES[4] = SqMat[0][0][4]; // M11
          SQMATRICES[5] = SqMat[0][0][7]; // M22
          SQMATRICES[6] = SqMat[0][0][5]; // M12
          SQMATRICES[7] = SqMat[0][0][6]; // M21

          MATRICES[0] =Mat[0][2]; // B1T
          MATRICES[1] = Mat[0][3]; // B2T

          fesp[0] = FeSpaces[0][0];
          ferhs[0] = FeSpaces[0][0];
          ferhs[1] = FeSpaces[0][0];

          RHSs[0] = Rhs[0][0];
          RHSs[1] = RHSs[0] + N_U;

          aux = new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);

          Assemble2DSlipBC(N_FESpaces, fesp,
                           N_SquareMatrices, SQMATRICES,
                           N_RectMatrices, MATRICES,
                           N_Rhs, RHSs, ferhs,
                           DiscreteForm,
                           BoundaryConditions,
                           BoundValues,
                           aux, FeFunct[0][0],
                           FeFunct[1][0]);

//           TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION = 1;
         delete aux;

        } // if (TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION >= 1)


        Dscal(Mat[0][0]->GetN_Entries(), tau, Mat[0][0]->GetEntries() );  // B1
        Dscal(Mat[0][1]->GetN_Entries(), tau, Mat[0][1]->GetEntries());  // B2
        Dscal(Mat[0][2]->GetN_Entries(), tau, Mat[0][2]->GetEntries());  // B1T
        Dscal(Mat[0][3]->GetN_Entries(), tau, Mat[0][3]->GetEntries());  // B2T

        gamma = 0.;

        // since rhs depends on moving grid so its better to use current geo
        Daxpy(N_Active, tau, Rhs[0][0], B);
        Daxpy(N_Active, tau, Rhs[0][0]+N_U, B+N_U);
 
        // update rhs by Laplacian and convective term initialy by current time step
        // scaled by current sub time step length and theta2
        // currently : M := M + gamma A
        // M = M + (-gamma - tau*TDatabase::TimeDB->THETA2) A
        MatAdd(SqMat[0][0][4], SqMat[0][0][0],
                     -gamma - tau*TDatabase::TimeDB->THETA2);
        MatAdd(SqMat[0][0][5], SqMat[0][0][1],
                     -gamma - tau*TDatabase::TimeDB->THETA2);
        MatAdd(SqMat[0][0][6], SqMat[0][0][2],
                     -gamma - tau*TDatabase::TimeDB->THETA2);
        MatAdd(SqMat[0][0][7], SqMat[0][0][3],
                     -gamma - tau*TDatabase::TimeDB->THETA2);

       // set current factor of steady state matrix
        gamma = -tau*TDatabase::TimeDB->THETA2;
    
        memset(defect, 0, N_Unknowns*SizeOfDouble);
        MatVectActive(SqMat[0][0][4], oldsol, defect);
        Daxpy(N_Active, 1, defect, B);
        memset(defect, 0, N_Unknowns*SizeOfDouble);
        MatVectActive(SqMat[0][0][5], oldsol+N_U, defect);
        Daxpy(N_Active, 1, defect, B);
        memset(defect, 0, N_Unknowns*SizeOfDouble);
        MatVectActive(SqMat[0][0][6], oldsol, defect+N_U);
        Daxpy(N_Active, 1, defect+N_U, B+N_U);
        memset(defect, 0, N_Unknowns*SizeOfDouble);
        MatVectActive(SqMat[0][0][7], oldsol+N_U, defect+N_U);
        Daxpy(N_Active, 1, defect+N_U, B+N_U);

 
    // set Dirichlet values
    // RHSs[0] still available from assembling
    memcpy(B+N_Active, RHSs[0]+N_Active, (N_U-N_Active)*SizeOfDouble);
    memcpy(B+N_Active+N_U, RHSs[1]+N_Active,(N_U-N_Active)*SizeOfDouble);

    // copy Dirichlet values from rhs into Sol[0][mg_level-1]
    memcpy(Sol[0][0]+N_Active, RHSs[0]+N_Active, (N_U-N_Active)*SizeOfDouble);
    memcpy(Sol[0][0]+N_Active+N_U, RHSs[1]+N_Active, (N_U-N_Active)*SizeOfDouble);

   //=====================================================================
   // the stiffness matrix is stored on M11, (M12, M21, M22)
   // assembling of system matrix
   //========================================================================
    MatAdd(SqMat[0][0][4], SqMat[0][0][0],
                     -gamma + tau*TDatabase::TimeDB->THETA1);
    MatAdd(SqMat[0][0][5], SqMat[0][0][1],
                     -gamma + tau*TDatabase::TimeDB->THETA1);
    MatAdd(SqMat[0][0][6], SqMat[0][0][2],
                     -gamma + tau*TDatabase::TimeDB->THETA1);
    MatAdd(SqMat[0][0][7], SqMat[0][0][3],
                     -gamma + tau*TDatabase::TimeDB->THETA1);

       // set current factor of steady state matrix
    gamma = tau*TDatabase::TimeDB->THETA1;    
    
//     cout << "Remesh sol 2 " << Ddot(N_Unknowns, Sol[0][0], Sol[0][0])<<endl;    
  //======================================================================
  // nonlinear loop
  //======================================================================
    N_LinIterCurr = 0;
    solver_time_curr = 0;

    for(j=0;j<Max_It;j++)
     {
      memcpy(oldsol, Sol[0][0], SizeOfDouble*N_Unknowns);
      memset(defect, 0, N_Unknowns*SizeOfDouble);

      SQMATRICES[0] = SqMat[0][0][4];
      SQMATRICES[1] = SqMat[0][0][5];
      SQMATRICES[2] = SqMat[0][0][6];
      SQMATRICES[3] = SqMat[0][0][7];
      MATRICES[0] = Mat[0][0];
      MATRICES[1] = Mat[0][1];
      MATRICES[2] = Mat[0][2];
      MATRICES[3] = Mat[0][3];
    
      if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE) 
        IntoL20FEFunction(Sol[0][0]+2*N_U, N_P, FeSpaces[1][0],
                        velocity_space_code, pressure_space_code);

      Defect(sqmatrices,matrices,Sol[0][0],B,defect);

//      for(i=0; i<200; i++)
//              cout << i<< " rhs 1 " << B[i] << " rhs 2 " << B[i+N_U] << endl;
//               cout << i<< "defect 1 " << defect[i] << " defect 2 " << defect[i+N_U] << endl;
     
//       cout << i<< " defect 1 " << Sol[0][0][i] << " defect 2 " << Sol[0][0][i+N_U] << endl;


// // exit(0);

      if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE) 
          IntoL20Vector2D(defect+2*N_U, N_P,pressure_space_code);

       residual =  Ddot(N_Unknowns, defect, defect);
       impuls_residual = Ddot(2*N_U, defect, defect);
       OutPut("ReMesh nonlinear step " << setw(3) << j);
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


      if((((sqrt(residual)<=limit)||(j==Max_It-1))) && (j>=TDatabase::ParamDB->SC_MINIT))
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
       t1 = GetTime();
       DirectSolver(SQMATRICES[0], SQMATRICES[1],
                    SQMATRICES[2], SQMATRICES[3],
                    MATRICES[2], MATRICES[3],
                    MATRICES[0], MATRICES[1],
                    B, Sol[0][0]);
    
       t2 = GetTime();
       solver_time_curr = t2-t1;
       solver_time += solver_time_curr;

     //======================================================================
     // end solve linear system
     //======================================================================

     // restore mass matrices by subtracting the A-matrices
        MatAdd(SqMat[0][0][4], SqMat[0][0][0], -gamma);
        MatAdd(SqMat[0][0][5], SqMat[0][0][1], -gamma);
        MatAdd(SqMat[0][0][6], SqMat[0][0][2], -gamma);
        MatAdd(SqMat[0][0][7], SqMat[0][0][3], -gamma);

        gamma = 0;
       
     //======================================================================
     // assemble new matrix due to nonlinearity
     //======================================================================
        DiscreteForm = DiscreteFormNLGalerkin;

        N_RectMatrices = 0;
        N_Rhs = 0;
        N_FESpaces = 3;
        N_SquareMatrices = 2;

        SQMATRICES[0] = SqMat[0][0][0];
        SQMATRICES[1] = SqMat[0][0][3];
        SQMATRICES[0]->Reset();
        SQMATRICES[1]->Reset();

        fesp[0] = FeSpaces[0][0];
        fesp[1] = FeSpaces[1][0];
        fesp[2] = FeSpaces[2][0];

        fefct[0] = FeFunct[0][0]; // u1
        fefct[1] = FeFunct[1][0]; // u2
        fefct[2] = FeFunct[3][0]; // w1
        fefct[3] = FeFunct[4][0]; // w1

      //======================================================================
      // A_11, (A_22)
      // no assembling of rhs
      //======================================================================
       aux =  new TAuxParam2D(MovingTNSN_FESpaces_Axial3D, MovingTNSN_Fct_Axial3D,
                           MovingTNSN_ParamFct_Axial3D,
                           MovingTNSN_FEValues_Axial3D,
                           fesp+1, fefct,
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
        MatAdd(SqMat[0][0][0], SqMat[0][0][8], 1.);
        MatAdd(SqMat[0][0][3], SqMat[0][0][9], 1.);

        delete aux;

        // slip type bc detected, modify matrices accordingly
       if (TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION >= 1)
        {
          N_FESpaces = 1;
          N_SquareMatrices = 2;
          N_RectMatrices = 0;
          N_Rhs = 2;
          DiscreteForm = NULL;

          SQMATRICES[0] = SqMat[0][0][0]; // A11
          SQMATRICES[1] = SqMat[0][0][3]; // A22

          MATRICES[0] =Mat[0][2]; // B1T
          MATRICES[1] = Mat[0][3]; // B2T

          fesp[0] = FeSpaces[0][0];
          ferhs[0] = FeSpaces[0][0];
          ferhs[1] = FeSpaces[0][0];

          RHSs[0] = Rhs[0][0];
          RHSs[1] = RHSs[0] + N_U;

          aux = new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);

          Assemble2DSlipBC(N_FESpaces, fesp,
                           N_SquareMatrices, SQMATRICES,
                           N_RectMatrices, MATRICES,
                           N_Rhs, RHSs, ferhs,
                           DiscreteForm,
                           BoundaryConditions,
                           BoundValues,
                           aux, FeFunct[0][0],
                           FeFunct[1][0]);
          delete aux;

        } // if (TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION >= 1)


       //======================================================================
       // end of assemble new matrix due to nonlinearity
       //======================================================================

       // build stiffness matrix for next nonlinear iteration step
       // stiffness matrix (left upper block) is stored on
       // M11, (M12, M21, M22)
       // M = M +  tau*TDatabase::TimeDB->THETA1 A
          MatAdd(SqMat[0][0][4], SqMat[0][0][0],
                 tau*TDatabase::TimeDB->THETA1);
          MatAdd(SqMat[0][0][5], SqMat[0][0][1],
                 tau*TDatabase::TimeDB->THETA1);
          MatAdd(SqMat[0][0][6], SqMat[0][0][2],
                 tau*TDatabase::TimeDB->THETA1);
          MatAdd(SqMat[0][0][7], SqMat[0][0][3],
                 tau*TDatabase::TimeDB->THETA1);

          // set current factor of steady state matrix
          gamma = tau*TDatabase::TimeDB->THETA1;
     } // for(j=0;j<Max_It;j++)

    if(reparam)
    {
     RefGridPos->DataToGrid();
     reparam = FALSE;
     
/*// ============================================================================================      
       length = BeginIndex[3][0][1] -  BeginIndex[3][0][0];
       for(i=0;i<N_Cells_P2;i++)
        {
         DOF = GlobalNumbers[3][0] + BeginIndex[3][0][GlobalCell_Index[1][i]];
         DOF_P2 = GlobalNumbers[3][2] + BeginIndex[3][2][i];

         for(jj=0;jj<length;jj++)
          {
           k = DOF[jj];
           l1 = DOF_P2[jj];
           Sol[3][0][k] = Sol[3][2][l1];
          }
         }

       MapSurfToDomain(IFaceFeFunct[0], FeFunct[6][0], N_List[0], N_List[1]);
       
       os.seekp(std::ios::beg);
       if(img<10) os <<"VTK/"<< VtkBaseName<<"ReParam.0000"<<ReParam_img<<".vtk" << ends;
       else if(img<100) os <<"VTK/"<< VtkBaseName<<"ReParam.000"<<ReParam_img<<".vtk" << ends;
       else if(img<1000) os << "VTK/"<<VtkBaseName<<"ReParam.00"<<ReParam_img<<".vtk" << ends;
       else if(img<10000) os << "VTK/"<<VtkBaseName<<"ReParam.0"<<ReParam_img<<".vtk" << ends;
       else  os << "VTK/"<<VtkBaseName<<"."<<ReParam_img<<"ReParam.vtk" << ends;
       Output->WriteVtk(os.str().c_str());
       ReParam_img++;
// ============================================================================================   */ 
    }
    else
    {remeshed = FALSE;}

    }//  if(remeshed)

#endif  // MovingMesh

      }// for(l=0;l<N_SubSteps;l++
    }//  for (methods=0;methods<time_discs;metho
    
  
#ifdef __MASSTRANSTEST__
   if( fabs(TDatabase::TimeDB->CURRENTTIME  - 0.1) < 1.e-8  ||
       fabs(TDatabase::TimeDB->CURRENTTIME  - 1.0) < 1.e-8  ||
       fabs(TDatabase::TimeDB->CURRENTTIME  - 2.0) < 1.e-8  ||
       fabs(TDatabase::TimeDB->CURRENTTIME  - 4.0) < 1.e-8  ||
       fabs(TDatabase::TimeDB->CURRENTTIME  - 7.0) < 1.e-8  ||
       fabs(TDatabase::TimeDB->CURRENTTIME  - 10.0) < 1.e-8  )
    {
       length = BeginIndex[3][0][1] -  BeginIndex[3][0][0];
       for(i=0;i<N_Cells_P2;i++)
        {
         DOF = GlobalNumbers[3][0] + BeginIndex[3][0][GlobalCell_Index[1][i]];
         DOF_P2 = GlobalNumbers[3][2] + BeginIndex[3][2][i];

         for(jj=0;jj<length;jj++)
          {
           k = DOF[jj];
           l1 = DOF_P2[jj];
           Sol[3][0][k] = Sol[3][2][l1];
          }
         }

      os.seekp(std::ios::beg);
     if(N_BData<10) os << "BDData/Gamma000"<<N_BData<<".data" << ends;
     else if(N_BData<10) os << "BDData/Gamma00"<<N_BData<<".data" << ends;
     std::ofstream dat(os.str().c_str());

     if (!dat)
      {
       cerr << "cannot open file for output" << endl;
       exit(0);
      }
     dat << "%% Surfactant data at Y= 4. created by MooNMD" << endl;
     dat << "%% Current Reference Time :" << TDatabase::TimeDB->CURRENTTIME << endl;
     dat << "%% x, y, Gamma" << endl;


     y=4.;
     c_x[0] = 1.;
     h = (5. - 1.)/40.;
     for(i=0;i<41;i++)
      {
       c_x[i] = c_x[0] + h*double(i);
       FeFunct[5][0]->FindGradient(c_x[i], y, Xvalues);
       c_xVal[i] = Xvalues[0];
//        cout <<" c_xVal[i] " << c_x[i]<<" c_xVal[i] " <<  c_xVal[i] << endl;
       dat <<  c_x[i]<< "  " << y<<"  " << Xvalues[0] <<endl;
      }

     dat.close();
      cout << endl;
      cout << "Surfactant data wrote into file " << endl;
     N_BData++;
     }
#endif

  if(((m % TDatabase::TimeDB->STEPS_PER_IMAGE) == 0) || m==1)
    {
 
     if(TDatabase::ParamDB->WRITE_VTK)
      {
#ifdef __WITHSURFACTANT__ 
       length = BeginIndex[3][0][1] -  BeginIndex[3][0][0];
       for(i=0;i<N_Cells_P2;i++)
        {
         DOF = GlobalNumbers[3][0] + BeginIndex[3][0][GlobalCell_Index[1][i]];
         DOF_P2 = GlobalNumbers[3][2] + BeginIndex[3][2][i];

         for(jj=0;jj<length;jj++)
          {
           k = DOF[jj];
           l1 = DOF_P2[jj];
           Sol[3][0][k] = Sol[3][2][l1];
          }
         }

      MapSurfToDomain(IFaceFeFunct[0], FeFunct[6][0], N_List[0], N_List[1]);

      GetSurfactMass(FeFunct[6][0], IFaceFeFunct[0], 
                     N_List[0], N_List[1], Surf_Mass);

      Get_FeFunction2DMass(FeFunct[5][2], Params);

//       Surf_Mass[0] *= double(TDatabase::ParamDB->REACTOR_P16);
//       Params[0] *=double(TDatabase::ParamDB->REACTOR_P20);

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
#else
      GetSurfactMass(FeFunct[6][0], IFaceFeFunct[0], 
                     N_List[0], N_List[1], Surf_Mass);  // calculation of sphericity  
#endif
      if(FluidFlow)
       {
        Get_KE(VeloVect[0][0], Params);
        CurrVolume = 2.*Pi*Params[0];
        MovBoundVert[0][0]->GetCoords(Lx, Ly);
        MovBoundVert[1][0]->GetCoords(Rx, Ry);
    
        radius = pow((3.*CurrVolume/(4.*Pi)), (1./3.));
        sphericity =  4.*Pi*radius*radius/(Surf_Mass[1]);
    
        OutPut(setw(20)<<"Time, Volume , Volume Diff, Top, Bottom: " << TDatabase::TimeDB->CURRENTTIME
               <<"   "<< CurrVolume<<"   "<< CurrVolume - InitVolume<<"   "<< Ly<<"   "<<Ry<< endl);
        OutPut(setw(20)<<"Time, KE , x_mass, y_mass, U1_Rise, U2_Rise, SurfArea, sphericity: " 
               << TDatabase::TimeDB->CURRENTTIME
               <<"   "<< Params[1]<<"   "<< Params[2] <<"   "<< Params[3]<<"   "<< Params[4]
               <<"   "<<Params[5]<<"   "<< Surf_Mass[1]<<"   "<< sphericity<<endl);
        }
        
       // output rise velo inside the bubble
      length = BeginIndex[2][0][1] -  BeginIndex[2][0][0];
      for(i=0;i<N_Cells_P1;i++)
       {
        DOF = GlobalNumbers[2][0] + BeginIndex[2][0][GlobalCell_Index[0][i]];
        for(j=0;j<length;j++)
         {
          k = DOF[j];
          outputsol[k] = Params[5];
         }
        }

       os.seekp(std::ios::beg);
       if(img<10) os << "VTK/"<<VtkBaseName<<".0000"<<img<<".vtk" << ends;
       else if(img<100) os << "VTK/"<<VtkBaseName<<".000"<<img<<".vtk" << ends;
       else if(img<1000) os << "VTK/"<<VtkBaseName<<".00"<<img<<".vtk" << ends;
       else if(img<10000) os << "VTK/"<<VtkBaseName<<".0"<<img<<".vtk" << ends;
       else  os << "VTK/"<<VtkBaseName<<"."<<img<<".vtk" << ends;
       Output->WriteVtk(os.str().c_str());
       
//        os.seekp(std::ios::beg);
//        if(img<10) os << "VTK/"<<VtkBaseName<<".0000"<<img<<".plt" << ends;
//        else if(img<100) os << "VTK/"<<VtkBaseName<<".000"<<img<<".plt" << ends;
//        else if(img<1000) os << "VTK/"<<VtkBaseName<<".00"<<img<<".plt" << ends;
//        else if(img<10000) os << "VTK/"<<VtkBaseName<<".0"<<img<<".plt" << ends;
//        else  os << "VTK/"<<VtkBaseName<<"."<<img<<".plt" << ends; 
//        Output->WriteBinaryPlt(os.str().c_str());
       img++;
       
     if(((m % (10*TDatabase::TimeDB->STEPS_PER_IMAGE)) == 0) || m==1)       
       PrintSurfSurfactant(N_MovVert[6], MovBoundVert[6], FeFunct[6][0], N_BData);         
      }
      
//    if(((m % (20*TDatabase::TimeDB->STEPS_PER_IMAGE)) == 0) || m==1)
//     { 
//      os.seekp(std::ios::beg);
//      if(N_BData<10) os << "Gamma000"<<N_BData<<".data" << ends;
//      else if(N_BData<100) os << "Gamma00"<<N_BData<<".data" << ends;
//      else if(N_BData<1000) os << "Gamma0"<<N_BData<<".data" << ends;
//      else os << "Gamma"<<N_BData<<".data" << ends;
//      std::ofstream dat(os.str().c_str());
// 
//      if (!dat)
//       {
//        cerr << "cannot open file for output" << endl;
//        exit(0);
//       }
// 
//       dat << "%% Boundary data created for droplet by MooNMD" << endl;
//       dat << "%% Current Time :" << TDatabase::TimeDB->CURRENTTIME << endl;
// 
//       for(k=0;k<N_MovVert[6];k++)
//        {
//         MovBoundVert[6][k]->GetCoords(x1, y1);
//         dat << x1 << " " <<  y1<< endl;
//        }
//        
//        // set end vertices again for closed curve
//         MovBoundVert[0][0]->GetCoords(x1, y1);
//         dat << x1 << " " <<  y1<< endl;
// 
//       dat.close();
//       cout << endl;
//       OutPut( "Boundary wrote output into file " << N_BData <<endl);
//       N_BData++;
//     }

    } // if(((m % TDatabase::TimeDB->STEPS_PER_IMAGE) == 0) || m==1)

  } // while(TDatabase::TimeDB->
 
 
  t4 =  GetTime();
  total_time += t4 - t3;
  OutPut("total running time: " << total_time << endl);
  CloseFiles();
  return 0;
}


