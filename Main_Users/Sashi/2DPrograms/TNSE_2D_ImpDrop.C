// =======================================================================
// 
// Purpose:     Main program for impinging droplet
//
// Author:     Sangeetha Rajasekaran, Sashikumaar Ganesan
// modified    10.06.2010 
// ======================================================================= 

#include <Domain.h>
#include <Database.h>
#include <DiscreteForm2D.h>
#include <FEDatabase2D.h>
#include <FESpace2D.h>
#include <SquareStructure2D.h>
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
#include "../Examples/TNSE_2D/Drop_Imping_Axial3D.h"
// #include "../Examples/TNSE_2D/DropHeat_imping_axial3D.h"

extern "C"
{
  void triangulate(char*, struct triangulateio*,
		   struct triangulateio*, struct triangulateio*);
}
// sorting curved surfrace vertices - general case
void Sort_Imping(TBaseCell **cell, TVertex **Vertex, int *CellNo, int *EdgeNo, int N,
                     double X0, double Y0 )
{
 int i, j, k, temp, test=0;
 
 double x, y, x1, y1;
 
 TVertex *temp_vert;
 TBaseCell *temp_cell;

 // finding the right(starting) vertex
  for(i=0;i<N;i++)
   {
    Vertex[i]->GetCoords(x, y);
    if(  sqrt((x-X0)*(x-X0) +(y-Y0)* (y-Y0))<1e-5  )
      {
//        cout << " sorting " << x << ' ' << y<<endl;
       temp_vert = Vertex[0];
       Vertex[0] = Vertex[i];
       Vertex[i] = temp_vert;

       temp_cell = cell[0];
       cell[0] = cell[i];
       cell[i] = temp_cell;
       
       temp = EdgeNo[0];
       EdgeNo[0] = EdgeNo[i];
       EdgeNo[i] = temp;
       
       temp = CellNo[0];
       CellNo[0] = CellNo[i];
       CellNo[i] = temp;        

       test++;
      }
    if(test) break;
   }

  if(i==N)
  {
   cout<< "Error in finding start vert " << X0 << endl;   
   exit(0);
  }

  for(i=0;i<N-1;i++)
   {
    test = 0; 
    k = cell[i]->GetN_Edges();
    cell[i]->GetVertex((EdgeNo[i]+1) % k)->GetCoords(x, y);

     for(j=i+1;j<N;j++)
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
       
         temp = EdgeNo[j];
         EdgeNo[j] = EdgeNo[i+1];
         EdgeNo[i+1] = temp; 
         
         temp = CellNo[j];
         CellNo[j] = CellNo[i+1];
         CellNo[i+1] = temp;          
         
         test++;
        }
      if(test) break;  
     }
   }

//    print   
//   for(i=0;i<N;i++)
//    {
//     Vertex[i]->GetCoords(x, y);
//     cout<<i<<  " x : "<< x << " y : "<< y<<  " Angle of free Vertices "<<(180/Pi)*atan2(y,x)<<endl;
//    }
//   exit(0); 
}

// ====================================================================
// Get the inner angles of the cells in whole domain
// ====================================================================
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


void SolveGridEquation(double **Entries, double *sol, double *rhs,
                       int *KCol, int *RowPtr, int N_DOF)
{
  int i,j,k, col, Diognal;
  double *Entries11, *Entries12, *Entries21, *Entries22;
  double sum1, sum2, max_sum1, max_sum2;
  int start, end;

  double max_error, error=1.e-12;
  int iter;

  Entries11 = Entries[0];
  Entries12 = Entries[1];
  Entries21 = Entries[2];
  Entries22 = Entries[3];
  
  max_error = 1.; iter = 0;
  while(max_error>error)
  {
    max_error = 0.0; iter++;
    for(i=0;i<N_DOF;i++)
    {
      start = RowPtr[i];
      end = RowPtr[i+1];
      sum1 = rhs[i];
      sum2 = rhs[i+N_DOF];
      for(k=start;k<end;k++)
      {
        col = KCol[k];
        if (col==i) Diognal = k;
        sum1 -= Entries11[k] * sol[col]
              +Entries12[k] * sol[col+N_DOF];
        sum2 -= Entries21[k] * sol[col]
              +Entries22[k] * sol[col+N_DOF];
      } // endfor k
     // sol[i] += sum1/Entries11[start];
     // sol[i+N_DOF] += sum2/Entries22[start];
        sol[i] += sum1/Entries11[Diognal];
        sol[i+N_DOF] += sum2/Entries22[Diognal];
      if(max_error<fabs(sum1/Entries11[Diognal])) max_error = fabs(sum1/Entries11[Diognal]);
      if(max_error<fabs(sum2/Entries22[Diognal])) max_error = fabs(sum2/Entries22[Diognal]);
    } // endfor i
    if(iter == 1000) break;
  } // end while
//OutPut("Grid Solver: Number iteration "<<iter<<endl);
}



void GridVelo_imping(double **Entries, double *Sol, double *d, double *Rhs,
                     int *KCol, int *RowPtr,
                     TFEVectFunct2D *GridPos,
                     TFEVectFunct2D *AuxGridPos,
                     TFEVectFunct2D *Velocity, double dt,
                     TFEVectFunct2D *GridVelocity, 
                     TVertex ***MovBoundVert, int *N_MovVert,
                     TBaseCell **Free_Cells, int **IsoCellEdgeNos,
                     bool &reparam, TFEVectFunct2D *RefGridPos)
{
  int i,j,k,l,m,comp, N;
  int *VeloBeginIndex, *VeloGlobalNumbers;
  int *GridBeginIndex, *GridGlobalNumbers;
  int N_Cells, N_Vertices, N_Edges, N_LocalDOFs;
  int N_Levels, *DOF, *JointDOF, GridLength;
  int N_BoundaryNodes, N_LinePoints, IIso, N_Inner, N_;
  
  double xi[4], eta[4], X[4], Y[4], VX[4], VY[4];
  double FunctValues[4][MaxN_BaseFunctions2D];
  double FEValuesX[MaxN_BaseFunctions2D];
  double FEValuesY[MaxN_BaseFunctions2D];
  double *ValuesX, *ValuesY, *RefValueX, *RefValueY;
  double *ValuesVX, *ValuesVY, *NewValuesX, *NewValuesY;
  double s, t, x, y, h_tot, x0, x1, y0, y1;
  double res, oldres, *gridvelo, *Nx, *Ny, Ay;
  double *LineWeights, *zeta;
  double normalx, normaly, tangenx, tangeny, nx, ny, tx, ty;
  double un, hE, t0,t1, temp2, eps=1e-6;
  double h, hmin, hmax, hlimit, *IsoX, *IsoY;
   
  BF2DRefElements RefElement;
  TRefTrans2D *F_K;
  RefTrans2D RefTrans;
  TIsoBoundEdge *isojoint;  
  TMGLevel2D *Level;
  TFESpace2D *VelocitySpace, *GridSpace;
  TCollection *Coll;
  TBaseCell *cell;
  BoundTypes bdtype;
  TBoundEdge *BoundEdge;
  TBoundComp2D *BoundComp;
  FE2D FEId;
  TFE2D *Element;
  TFEDesc2D *FEDesc;
  BaseFunct2D BF;
  TBaseFunct2D *bf;
  boolean OnBoundary;
  TJoint *joint;
  TVertex **Vertices;
  TQuadFormula1D *qf1;
  QuadFormula1D LineQuadFormula;
  BoundCond Cond0, Cond1;
    
  VelocitySpace = Velocity->GetFESpace2D();
  VeloBeginIndex = VelocitySpace->GetBeginIndex();
  VeloGlobalNumbers = VelocitySpace->GetGlobalNumbers();
  ValuesVX = Velocity->GetValues();
  ValuesVY = ValuesVX + Velocity->GetLength();

  GridPos->GridToData();
  GridSpace = GridPos->GetFESpace2D();
  GridBeginIndex = GridSpace->GetBeginIndex();
  GridGlobalNumbers = GridSpace->GetGlobalNumbers();
  GridLength = GridPos->GetLength();
  ValuesX = GridPos->GetValues();
  ValuesY = ValuesX + GridLength;

  N_Inner = GridSpace->GetN_Inner();
  N_BoundaryNodes = GridLength - GridSpace->GetN_Inner();
  
  if(TDatabase::ParamDB->P5 > 0)
  {
   Nx = new double[N_BoundaryNodes];
   Ny = new double[N_BoundaryNodes];
   memset(Nx, 0, N_BoundaryNodes*SizeOfDouble);
   memset(Ny, 0, N_BoundaryNodes*SizeOfDouble);
  }

 // cout << "N_BoundaryNodes: " << N_BoundaryNodes << endl;
//   cout << GridLength << " " << N_DOF << endl;

  NewValuesX = AuxGridPos->GetValues();
  NewValuesY = NewValuesX + GridLength;

  memcpy(NewValuesX, ValuesX, 2*GridLength*SizeOfDouble);

  RefValueX = RefGridPos->GetValues();
  RefValueY = RefValueX + GridLength;   
  
  Coll = VelocitySpace->GetCollection();
  N_Cells = Coll->GetN_Cells();
//  cout << "N_Cells: " <<N_Cells<< endl;
  // determine outer normal vectors

  IIso = N_BoundaryNodes;
  // Outward normal no need if we move boundary with velocity
 if(TDatabase::ParamDB->P5 > 0)
  {
  // determine outer normal vectors
   for(i=0;i<N_Cells;i++)
    {
    // cout << "cell: " << i << endl;
    cell = Coll->GetCell(i);
    cell->SetClipBoard(i);
    N_Edges = cell->GetN_Edges();

    for(j=0;j<N_Edges;j++)
    {
      if( !(cell->GetJoint(j)->InnerJoint()) )
      {
        // cout << "joint: " << j << endl;
        cell->GetVertex(j)->GetCoords(x0, y0);
        cell->GetVertex((j+1)%N_Edges)->GetCoords(x1, y1);
        t = x1-x0;
        s = y1-y0;
        FEId = VelocitySpace->GetFE2D(i, cell);
        l = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
        LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(2*l);
        qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
        qf1->GetFormulaData(N_LinePoints, LineWeights, zeta);

        RefElement = TFEDatabase2D::GetRefElementFromFE2D(FEId);
        switch(RefElement)
        {
          case BFUnitTriangle:
            RefTrans = TriaIsoparametric;
            F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
            ((TTriaIsoparametric *)F_K)->SetCell(cell);
          break;

          case BFUnitSquare:
            RefTrans = QuadIsoparametric;
            F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
            ((TQuadIsoparametric *)F_K)->SetCell(cell);
          break;

          default:
            Error("only triangles and quadrilaterals are allowes" << endl);
            Error("file: " << __FILE__ << " line " << __LINE__ << endl);
            exit(-1);
        } // endswitch

        normalx = 0;
        normaly = 0;
        hE = 0;
        for(k=0;k<N_LinePoints;k++)
        {
          F_K->GetOuterNormal(j, zeta[k], nx, ny);
          F_K->GetTangent(j, zeta[k], tx, ty);
          t = sqrt(tx*tx+ty*ty);
          normalx += t*LineWeights[k]*nx;
          normaly += t*LineWeights[k]*ny;
          hE += t*LineWeights[k];
          // cout << "k= " << k << " " << nx << " " << ny << endl;
        } // endfor k

        DOF = GridGlobalNumbers + GridBeginIndex[i];

        switch(N_Edges)
        {
          case 3:
            switch(j)
            {
              case 0:
                l = DOF[0] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[1] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;

              case 1:
                l = DOF[1] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[2] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;

              case 2:
                l = DOF[2] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[0] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;
            } // endswitch j
          break;

          case 4:
            switch(j)
            {
              case 0:
                l = DOF[0] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[1] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;

              case 1:
                l = DOF[1] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[3] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;

              case 2:
                l = DOF[3] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[2] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;

              case 3:
                l = DOF[2] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[0] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;
            } // endswitch j
          break;
        } // endswitch N_Edges
      } // !InnerJoint
    } // endfor j
  } // endfor i

  N_ = IIso;
  // normalize normal vector
  for(i=0;i<N_;i++)
  {
    x = Nx[i];
    y = Ny[i];
    t = sqrt(x*x+y*y);
    Nx[i] /= t;
    Ny[i] /= t;
  }
 }
  // determine new position of boundary vertices
  for(i=0;i<N_Cells;i++)
  {
    cell  = Coll->GetCell(i);
    N_Edges = cell->GetN_Edges();
    OnBoundary = FALSE;
    for(j=0;j<N_Edges;j++)
    {
      if( !(cell->GetJoint(j)->InnerJoint()) || 
           (cell->GetJoint(j)->GetType() == InterfaceJoint)  || 
           (cell->GetJoint(j)->GetType() == IsoInterfaceJoint) )
       {
        OnBoundary = TRUE;
        joint = cell->GetJoint(j);
       }
    } // endfor j

    if(OnBoundary)
    {

      FEId = VelocitySpace->GetFE2D(i, cell);
      Element = TFEDatabase2D::GetFE2D(FEId);
      bf = Element->GetBaseFunct2D();
      N_LocalDOFs = Element->GetN_DOF();

      switch(N_Edges)
      {
        case 3:
          xi[0]  = 0; xi[1]  = 1; xi[2]  = 0;
          eta[0] = 0; eta[1] = 0; eta[2] = 1;
          VX[0] = VX[1] = VX[2] = 0;
          VY[0] = VY[1] = VY[2] = 0;
          N_Vertices = N_Edges;
        break;

        case 4:
          // NOTE: sorted due to number of DOFs !!!
          xi[0]  = -1; xi[1]  =  1; xi[2]  = -1; xi[3]  =  1;
          eta[0] = -1; eta[1] = -1; eta[2] =  1; eta[3] =  1;
          VX[0] = VX[1] = VX[2] = VX[3] = 0;
          VY[0] = VY[1] = VY[2] = VY[3] = 0;
          N_Vertices = N_Edges;
        break;

        default:
          Error("only triangles and quadrilaterals are allowed!" << endl);
          exit(-1);
      } // endswitch

      for(j=0;j<N_Vertices;j++)
        bf->GetDerivatives(D00, xi[j], eta[j], FunctValues[j]);

      DOF = VeloGlobalNumbers + VeloBeginIndex[i];

      for(j=0;j<N_LocalDOFs;j++)
      {
        k = DOF[j];
        s = ValuesVX[k];
        t = ValuesVY[k];

        for(l=0;l<N_Vertices;l++)
        {
          VX[l] += FunctValues[l][j]*s;
          VY[l] += FunctValues[l][j]*t;
        } // endfor l
      } // endfor j

      FEId = GridSpace->GetFE2D(i, cell);
      Element = TFEDatabase2D::GetFE2D(FEId);
      BF = Element->GetBaseFunct2D_ID();
      if( (BF != BF_C_T_P1_2D) && (BF != BF_C_Q_Q1_2D) )
      {
        Error("Grid Space must be conforming and of first order!" << endl);
        exit(-1);
      }  // endif
      bf = Element->GetBaseFunct2D();
      N_LocalDOFs = Element->GetN_DOF();

      switch(N_Edges)
      {
        case 3:
          xi[0]  = 0; xi[1]  = 1; xi[2]  = 0;
          eta[0] = 0; eta[1] = 0; eta[2] = 1;
          X[0] = X[1] = X[2] = 0;
          Y[0] = Y[1] = Y[2] = 0;
          N_Vertices = N_Edges;
        break;

        case 4:
          // NOTE: sorted due to number of DOF !!!
          xi[0]  = -1; xi[1]  =  1; xi[2]  = -1; xi[3]  =  1;
          eta[0] = -1; eta[1] = -1; eta[2] =  1; eta[3] =  1;
          X[0] = X[1] = X[2] = X[3] = 0;
          Y[0] = Y[1] = Y[2] = Y[3] = 0;
          N_Vertices = N_Edges;
        break;

        default:
          Error("only triangles and quadrilaterals are allowed!" << endl);
          exit(-1);
      } // endswitch

      DOF = GridGlobalNumbers + GridBeginIndex[i];

      for(j=0;j<N_Vertices;j++)
      {
        l = DOF[j];
        k = l - N_Inner;
        if(k>=0)
        {
          if((TDatabase::ParamDB->P5 > 0) && (ValuesY[l] != 0) )
          {
            un = VX[j]*Nx[k] + VY[j]*Ny[k];
            NewValuesX[l] = ValuesX[l] + dt*un*Nx[k];
            NewValuesY[l] = ValuesY[l] + dt*un*Ny[k];
          }
          else
          {
	    if(ValuesX[l] == 0 )
	    { NewValuesX[l] = ValuesX[l]; }
	    else
	    { NewValuesX[l] = ValuesX[l] + dt*VX[j]; }

            if(ValuesY[l] == 0) 
	     { NewValuesY[l] = ValuesY[l];  }
            else    
	     { NewValuesY[l] = ValuesY[l] + dt*VY[j];  }
	    
	  }
       } //  if(k>=0)
 //    Due to spline approximation solid boundary end vertices may take negative y value
        if(NewValuesY[l]<0.0 ) NewValuesY[l] = 0.0;
        if( fabs(NewValuesX[l]) < 1e-10 ) NewValuesX[l] = 0.0;
      } // endfor j
    } // endif
  } // endfor i

// cout << " dt " << dt <<endl;
/*
   for(i=0;i<GridLength;i++)
 cout << i <<"  ---  " <<ValuesX[i] << "  ---  " << NewValuesX[i] << endl;
*/
// exit(0);
  
   MovBoundVert[1][0]->GetCoords(x, Ay);   

   AuxGridPos->DataToGrid();
   
  //======================================================================  
  //  Reparametrization of free surface - Begin
  //======================================================================  
  if(!reparam)
  {
   h_tot = 0;
   hmin = 100;
   hmax = 0.0;
   
   for(k=0;k<N_MovVert[2];k++)
    {
     MovBoundVert[2][k]->GetCoords(x1, y1);

     if(k==N_MovVert[2]-1)
     { MovBoundVert[1][0]->GetCoords(x, y);}
     else
     { MovBoundVert[2][k+1]->GetCoords(x, y); }

     h = sqrt((x-x1)*(x-x1)+(y-y1)*(y-y1));
     h_tot +=h;
     if (h < hmin) hmin = h;
     if (h > hmax) hmax = h;
    } // for(k=0;k<N_MovVert[2];k++) 
   
    h_tot /= (double)N_MovVert[2];
    hlimit =  0.8*h_tot;   
   }
   
   if ( ((hmin < hlimit) || (hmax > 3.*h_tot/2.)) ||  reparam )
   { 
 
    //before reparam update iso points, since we use interpolated cubic spline
    //which pass through iso points also
    IsoX = new double[IIso];
    IsoY = new double[IIso];
    
    IIso = 0;
    for(i=0;i<N_MovVert[2];i++)
    {
     cell = Free_Cells[i];
      j = IsoCellEdgeNos[1][i];
      joint = cell->GetJoint(j);
      isojoint = (TIsoBoundEdge *)joint;
      k = isojoint->GetN_Vertices();
      Vertices = isojoint->GetVertices();
      FEId = VelocitySpace->GetFE2D(IsoCellEdgeNos[0][i], cell);
      FEDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(FEId);
      m = FEDesc->GetN_JointDOF();
      if(m == k+2)
        {
          JointDOF = FEDesc->GetJointDOF(j);
          DOF =  VeloGlobalNumbers+VeloBeginIndex[IsoCellEdgeNos[0][i]];
          for(l=0;l<k;l++)
          {
            m = DOF[JointDOF[l+1]];
            Vertices[l]->GetCoords(IsoX[IIso], IsoY[IIso]);
            if(TDatabase::ParamDB->P5 > 0)
            {
              un = ValuesVX[m]*Nx[IIso+N_BoundaryNodes]
                  + ValuesVY[m]*Ny[IIso+N_BoundaryNodes];
              x  = IsoX[IIso] + dt*un*Nx[IIso+N_BoundaryNodes];
              y  = IsoY[IIso] + dt*un*Ny[IIso+N_BoundaryNodes];
            }
            else
            {
             x  = IsoX[IIso] + dt*ValuesVX[m];
             y  = IsoY[IIso] + dt*ValuesVY[m];
            }
            
           if(y<=0) y = 1e-5;
//            if(fabs(x)<1e-12) x = 0.;
           
           Vertices[l]->SetCoords(x, y);
           IIso++;
          } // endfor l
        }
       else
        {
          // approximation order of isoparametric boundary and velocity
          // element must be the same
          Error("No match in isoparametric case" << endl);
          exit(-1);
        } //  if(m == k+2)     
    }// for(i=0;i<N_MovVert[2];i++)

    ReParam_axial3D_U(N_MovVert[2], Free_Cells,  IsoCellEdgeNos[1], IsoCellEdgeNos[0],
                      Velocity, NULL, NULL, FALSE);  
     
    reparam = TRUE;   
    RefGridPos->GridToData();    
    Daxpy(2*GridLength, -1, ValuesX, RefValueX); // now reparamdisp in RefValueX

    //back to orig mesh (no reparam movement in calculation of free surf w)
    AuxGridPos->DataToGrid(); 
    
   //restore iso points
    IIso = 0;
    for(i=0;i<N_MovVert[2];i++)
    {
     cell = Free_Cells[i];
     j = IsoCellEdgeNos[1][i];      
     joint = cell->GetJoint(j);
     isojoint = (TIsoBoundEdge *)joint;
     k = isojoint->GetN_Vertices();
     Vertices = isojoint->GetVertices();
     FEId = VelocitySpace->GetFE2D(IsoCellEdgeNos[0][i], cell);
     FEDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(FEId);
     m = FEDesc->GetN_JointDOF();
     for(l=0;l<k;l++)
      {
        Vertices[l]->SetCoords(IsoX[IIso], IsoY[IIso]);
        IIso++;
      }
    }// for(i=0;i<N_Cells;i++)    
    
    delete [] IsoX; delete [] IsoY; 
   } // if ( ((hmin < hlimit) || (hmax > 3.*h_tot/2.)) ||  reparam )
   //======================================================================       

   MovBoundVert[2][0]->GetCoords(x, y); // right wetting point

   y=0.;
   h_tot = x;
   h_tot /= (double)N_MovVert[0];
   for(i=1;i<N_MovVert[0];i++)
     MovBoundVert[0][i]->SetCoords(h_tot*(double)i, y);
 

// axial boundary
   MovBoundVert[1][0]->GetCoords(x, y);
   
   N=N_MovVert[1];      
   h_tot = (y-Ay)/(double)N;   
   N--;
   
    for(i=1;i<N;i++)
    {
     MovBoundVert[1][i]->GetCoords(x, y);
     // cout<< " y " << y <<" new y " << y +((double)(N-i))*h_tot<<endl;      
     y += ((double)(N-i))*h_tot;
     MovBoundVert[1][i]->SetCoords(x, y);   
    }      
   
//    exit(0);
   
//    x=0.;
//    h_tot = -y;
//    h_tot /= (double)N_MovVert[1];
//    for(i=1;i<N_MovVert[1];i++)
//     MovBoundVert[1][i]->SetCoords(x,  y +  h_tot*(double)i );
//      
   AuxGridPos->GridToData();   
   GridPos->DataToGrid();
   
   
//       MovBoundVert[1][0]->GetCoords(x, y);
//     cout << " x " << x <<" y " << y <<endl;  
//     exit(0); 
    
   memset(Rhs, 0, 2*GridLength*SizeOfDouble);

   memcpy(d, NewValuesX, 2*GridLength*SizeOfDouble);
   Daxpy(2*GridLength, -1, ValuesX, d);
   memcpy(Rhs + (GridLength-N_BoundaryNodes), d+(GridLength-N_BoundaryNodes),
          N_BoundaryNodes*SizeOfDouble);
   memcpy(Rhs + (2*GridLength-N_BoundaryNodes), d+(2*GridLength-N_BoundaryNodes),
          N_BoundaryNodes*SizeOfDouble);

  memset(Sol, 0 , 2*GridLength*SizeOfDouble);
  memcpy(Sol + (GridLength-N_BoundaryNodes), d+(GridLength-N_BoundaryNodes),
         N_BoundaryNodes*SizeOfDouble);
  memcpy(Sol + (2*GridLength-N_BoundaryNodes), d+(2*GridLength-N_BoundaryNodes),
         N_BoundaryNodes*SizeOfDouble);
	 
//     for(i=0;i<GridLength;i++)
//      cout<< i <<"  ---  "<< Rhs[i] << "  ---  " << Rhs[i+GridLength] << endl;
//     
  SolveGridEquation(Entries, Sol, Rhs, KCol, RowPtr, GridLength);

  gridvelo = GridVelocity->GetValues();
  memcpy(gridvelo, Sol, 2*GridLength*SizeOfDouble);
  Dscal(2*GridLength, 1./dt, gridvelo);

  //======================================================================  
  //  Reparametrization of free surface 
  //======================================================================      
  if(reparam)
  {
   memset(Rhs, 0, 2*GridLength*SizeOfDouble);    
   memcpy(Rhs + (GridLength-N_BoundaryNodes), RefValueX+(GridLength-N_BoundaryNodes),
          N_BoundaryNodes*SizeOfDouble);
   memcpy(Rhs + (2*GridLength-N_BoundaryNodes), RefValueX+(2*GridLength-N_BoundaryNodes),
          N_BoundaryNodes*SizeOfDouble);
 
  memset(Sol, 0 , 2*GridLength*SizeOfDouble);
  memcpy(Sol + (GridLength-N_BoundaryNodes), RefValueX+(GridLength-N_BoundaryNodes),
         N_BoundaryNodes*SizeOfDouble);
  memcpy(Sol + (2*GridLength-N_BoundaryNodes), RefValueX+(2*GridLength-N_BoundaryNodes),
         N_BoundaryNodes*SizeOfDouble); 
 
  SolveGridEquation(Entries, Sol, Rhs, KCol, RowPtr, GridLength);
  
  Dscal(2*GridLength, 1./dt, Sol);
   
  //only inner mesh velo, since U will be interpolated in free surf reparm     
  Daxpy(GridLength-N_BoundaryNodes, 1., Sol, gridvelo);
  Daxpy(GridLength-N_BoundaryNodes, 1., Sol+GridLength, gridvelo+GridLength);  
  
//     cout<< "reparam : "<<endl;
//     exit(0);
  }  
  
//     for(i=0;i<GridLength;i++)
//      cout<< i <<"  ---  "<< gridvelo[i] << "  ---  " << gridvelo[i+GridLength] << endl;

//   delete [] d;
   
  if(TDatabase::ParamDB->P5 > 0)
  { delete [] Nx;  delete [] Ny; }
  
  
} // GridVelo_imping


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



void MoveGrid_imping(double **Entries, double *Sol, double *d, double *Rhs,
                  int *KCol, int *RowPtr,
                  TFEVectFunct2D *GridPos,
                  TFEVectFunct2D *Velocity, double dt,
                  TFEVectFunct2D *NewGridPos, 
                  TVertex ***MovBoundVert, int *N_MovVert,
                  TBaseCell **Free_Cells, int **IsoCellEdgeNos,
                  bool &reparam, int &N_ReParam)
{
  int i,j,k,l,m,comp, N;
  int *VeloBeginIndex, *VeloGlobalNumbers;
  int *GridBeginIndex, *GridGlobalNumbers;
  int N_Cells, N_Vertices, N_Edges, N_LocalDOFs;
  int N_Levels, *DOF, *JointDOF, GridLength, polydegree;
  int N_BoundaryNodes, N_LinePoints, IIso, N_Inner, N_;

  double xi[4], eta[4], X[4], Y[4], VX[4], VY[4];
  double FunctValues[4][MaxN_BaseFunctions2D];
  double FEValuesX[MaxN_BaseFunctions2D];
  double FEValuesY[MaxN_BaseFunctions2D];
  double *ValuesX, *ValuesY;
  double *ValuesVX, *ValuesVY;
  double *NewValuesX, *NewValuesY;
  double s, t, x, y, IsoX, IsoY, Ay;
  double x0, x1, y0, y1, h_tot, res, oldres;
  double *gridvelo, *Nx, *Ny, *LineWeights, *zeta;
  double normalx, normaly, tangenx, tangeny, nx, ny, tx, ty;
  double un, hE, t0,t1, temp2, x_max, x_min, temp, eps=1e-6;  
  
  TIsoBoundEdge *isojoint;
  TMGLevel2D *Level;  
  TFESpace2D *VelocitySpace, *GridSpace;
  TCollection *Coll;
  TBaseCell *cell;
  BoundTypes bdtype;
  TBoundEdge *BoundEdge;
  TBoundComp2D *BoundComp;
  FE2D FEId;
  TFE2D *Element;
  TFEDesc2D *FEDesc;
  BaseFunct2D BF;
  TBaseFunct2D *bf;
  boolean OnBoundary;
  TJoint *joint;
  TVertex **Vertices;
  TQuadFormula1D *qf1;
  QuadFormula1D LineQuadFormula;
  BF2DRefElements RefElement;
  TRefTrans2D *F_K;
  RefTrans2D RefTrans;  
  QuadFormula2D QuadFormula;
  BoundCond Cond0, Cond1;  
 
  VelocitySpace = Velocity->GetFESpace2D();
  VeloBeginIndex = VelocitySpace->GetBeginIndex();
  VeloGlobalNumbers = VelocitySpace->GetGlobalNumbers();
  ValuesVX = Velocity->GetValues();
  ValuesVY = ValuesVX + Velocity->GetLength();

  GridPos->GridToData();
  GridSpace = GridPos->GetFESpace2D();
  GridBeginIndex = GridSpace->GetBeginIndex();
  GridGlobalNumbers = GridSpace->GetGlobalNumbers();
  GridLength = GridPos->GetLength();
  ValuesX = GridPos->GetValues();
  ValuesY = ValuesX + GridLength;

  N_Inner = GridSpace->GetN_Inner();
  N_BoundaryNodes = GridLength - GridSpace->GetN_Inner();
  // cout << "N_BoundaryNodes: " << N_BoundaryNodes << endl;

//   d = new double[ 2*GridLength];

  if(TDatabase::ParamDB->P5 > 0)
  {  
   Nx = new double[2*N_BoundaryNodes]; // additional values for edge midpoints
   Ny = new double[2*N_BoundaryNodes]; // additional values for edge midpoints
   memset(Nx, 0, 2*N_BoundaryNodes*SizeOfDouble);
   memset(Ny, 0, 2*N_BoundaryNodes*SizeOfDouble);
  }
  //cout << "N_BoundaryNodes: " << N_BoundaryNodes << endl;
 // cout << GridLength << " " << N_DOF << endl;

  NewValuesX = NewGridPos->GetValues();
  NewValuesY = NewValuesX + GridLength;

  memcpy(NewValuesX, ValuesX, 2*GridLength*SizeOfDouble);

  Coll = VelocitySpace->GetCollection();
  N_Cells = Coll->GetN_Cells();

  // determine outer normal vectors
  IIso = N_BoundaryNodes;

    // Outward normal no need if we move boundary with velocity
  if(TDatabase::ParamDB->P5 > 0)
  {
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    N_Edges = cell->GetN_Edges();

    for(j=0;j<N_Edges;j++)
    {
      if( !(cell->GetJoint(j)->InnerJoint()) )
      {
        FEId = VelocitySpace->GetFE2D(i, cell);
        l = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
        LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(2*l);
        qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
        qf1->GetFormulaData(N_LinePoints, LineWeights, zeta);

        RefElement = TFEDatabase2D::GetRefElementFromFE2D(FEId);
        switch(RefElement)
        {
          case BFUnitTriangle:
            RefTrans = TriaIsoparametric;
            F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
            polydegree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
            QuadFormula = TFEDatabase2D::GetQFTriaFromDegree(2*polydegree-1);
//             ((TTriaIsoparametric *)F_K)->SetApproximationOrder(polydegree);
            ((TTriaIsoparametric *)F_K)->SetQuadFormula(QuadFormula);
            ((TTriaIsoparametric *)F_K)->SetCell(cell);
          break;

          case BFUnitSquare:
            RefTrans = QuadIsoparametric;
            F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
            polydegree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
            QuadFormula = TFEDatabase2D::GetQFQuadFromDegree(3*polydegree);
//             ((TQuadIsoparametric *)F_K)->SetApproximationOrder(polydegree);
            ((TQuadIsoparametric *)F_K)->SetQuadFormula(QuadFormula);
            ((TQuadIsoparametric *)F_K)->SetCell(cell);
          break;

          default:
            Error("only triangles and quadrilaterals are allowes" << endl);
            Error("file: " << __FILE__ << " line " << __LINE__ << endl);
            exit(-1);
        } // endswitch

        normalx = 0;
        normaly = 0;
        hE = 0;
        for(k=0;k<N_LinePoints;k++)
        {
          F_K->GetOuterNormal(j, zeta[k], nx, ny);
          F_K->GetTangent(j, zeta[k], tx, ty);
          t = sqrt(tx*tx+ty*ty);
          normalx += t * LineWeights[k] * nx;
          normaly += t * LineWeights[k] * ny;
          hE += t * LineWeights[k];
          // cout << "k= " << k << " " << nx << " " << ny << endl;
        } // endfor k

        DOF = GridGlobalNumbers + GridBeginIndex[i];
/*
        switch(N_Edges)
        {
          case 3:
            switch(j)
            {
              case 0:
                l = DOF[0] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
                l = DOF[1] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
              break;

              case 1:
                l = DOF[1] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
                l = DOF[2] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
              break;

              case 2:
                l = DOF[2] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
                l = DOF[0] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
              break;
            } // endswitch j
          break;

          case 4:
            switch(j)
            {
              case 0:
                l = DOF[0] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
                l = DOF[1] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
              break;

              case 1:
                l = DOF[1] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
                l = DOF[3] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
              break;

              case 2:
                l = DOF[3] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
                l = DOF[2] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
              break;

              case 3:
                l = DOF[2] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
                l = DOF[0] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
              break;
            } // endswitch j
          break;
        } // endswitch N_Edges
*/

// /*
        switch(N_Edges)
        {
          case 3:
            switch(j)
            {
              case 0:
                l = DOF[0] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[1] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;

              case 1:
                l = DOF[1] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[2] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;

              case 2:
                l = DOF[2] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[0] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;
            } // endswitch j
          break;

          case 4:
            switch(j)
            {
              case 0:
                l = DOF[0] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[1] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;

              case 1:
                l = DOF[1] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[3] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;

              case 2:
                l = DOF[3] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[2] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;

              case 3:
                l = DOF[2] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[0] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;
            } // endswitch j
          break;
        } // endswitch N_Edges
// */

        if(cell->GetJoint(j)->GetType() == IsoBoundEdge)
        {
          FEId = VelocitySpace->GetFE2D(i, cell);
          FEDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(FEId);
          N_LocalDOFs = FEDesc->GetN_JointDOF();
          t = 2.0/(N_LocalDOFs-1);
          for(k=1;k<N_LocalDOFs-1;k++)
          {
            /*
            Nx[IIso] += normalx;
            Ny[IIso] += normaly;
            */
            // /*
            s = -1.0 + k*t;
            F_K->GetOuterNormal(j, s, nx, ny);
            Nx[IIso] += nx;
            Ny[IIso] += ny;
            // */
            IIso++;
          } // endfor
        } // endif
      } // !InnerJoint
    } // endfor j
  } // endfor i

  N_ = IIso;
  // normalize normal vector
  for(i=0;i<N_;i++)
  {
    x = Nx[i];
    y = Ny[i];
    t = sqrt(x*x+y*y);
    Nx[i] /= t;
    Ny[i] /= t;

    // cout << setw(5) << i << "n = (" << Nx[i] << ", " << Ny[i] << ")" << endl;
  }
}

  // determine new position of boundary vertices
  for(i=0;i<N_Cells;i++)
  {
    cell  = Coll->GetCell(i);
    N_Edges = cell->GetN_Edges();
    OnBoundary = FALSE;
    for(j=0;j<N_Edges;j++)
    {
      if( !(cell->GetJoint(j)->InnerJoint()) )
       {
        OnBoundary = TRUE;
        joint = cell->GetJoint(j);
       }
    } // endfor j

    if(OnBoundary)
    {
      FEId = VelocitySpace->GetFE2D(i, cell);
      Element = TFEDatabase2D::GetFE2D(FEId);
      bf = Element->GetBaseFunct2D();
      N_LocalDOFs = Element->GetN_DOF();

      switch(N_Edges)
      {
        case 3:
          xi[0]  = 0; xi[1]  = 1; xi[2]  = 0;
          eta[0] = 0; eta[1] = 0; eta[2] = 1;
          VX[0] = VX[1] = VX[2] = 0;
          VY[0] = VY[1] = VY[2] = 0;
          N_Vertices = N_Edges;
        break;

        case 4:
          // NOTE: sorted due to number of DOFs !!!
          xi[0]  = -1; xi[1]  =  1; xi[2]  = -1; xi[3]  =  1;
          eta[0] = -1; eta[1] = -1; eta[2] =  1; eta[3] =  1;
          VX[0] = VX[1] = VX[2] = VX[3] = 0;
          VY[0] = VY[1] = VY[2] = VY[3] = 0;
          N_Vertices = N_Edges;
        break;

        default:
          Error("only triangles and quadrilaterals are allowed!" << endl);
          exit(-1);
      } // endswitch

      for(j=0;j<N_Vertices;j++)
        bf->GetDerivatives(D00, xi[j], eta[j], FunctValues[j]);

      DOF = VeloGlobalNumbers + VeloBeginIndex[i];

      for(j=0;j<N_LocalDOFs;j++)
      {
        k = DOF[j];
        s = ValuesVX[k];
        t = ValuesVY[k];

        for(l=0;l<N_Vertices;l++)
        {
          VX[l] += FunctValues[l][j]*s;
          VY[l] += FunctValues[l][j]*t;
        } // endfor l
      } // endfor j

      FEId = GridSpace->GetFE2D(i, cell);
      Element = TFEDatabase2D::GetFE2D(FEId);
      BF = Element->GetBaseFunct2D_ID();
      if( (BF != BF_C_T_P1_2D) && (BF != BF_C_Q_Q1_2D) )
      {
        Error("Grid Space must be conforming and of first order!" << endl);
        exit(-1);
      }  // endif
      bf = Element->GetBaseFunct2D();
      N_LocalDOFs = Element->GetN_DOF();


      DOF = GridGlobalNumbers + GridBeginIndex[i];

      for(j=0;j<N_Vertices;j++)
      {
        l = DOF[j];
        k = l - N_Inner;
        if(k>=0)
        {
          if((TDatabase::ParamDB->P5 > 0) && (ValuesY[l] != 0) )
          {
            un = VX[j]*Nx[k] + VY[j]*Ny[k];
            NewValuesX[l] = ValuesX[l] + dt*un*Nx[k];
            NewValuesY[l] = ValuesY[l] + dt*un*Ny[k];
          }
          else
          {
 
	    if(ValuesX[l] == 0 )
	    { NewValuesX[l] = ValuesX[l]; }
	    else
	    { NewValuesX[l] = ValuesX[l] + dt*VX[j]; }

            if(ValuesY[l] == 0) 
	     { NewValuesY[l] = ValuesY[l];  }
            else    
	     { NewValuesY[l] = ValuesY[l] + dt*VY[j];  }                
          }
        }
 //    Due to spline approximation solid boundary end vertices may take negative y value
         if(NewValuesY[l]<0.0 ) NewValuesY[l] = 0.0;
        if( fabs(NewValuesX[l]) < 1e-12 ) NewValuesX[l] = 0.0;
      } // endfor j
     } // endif
   } // endfor i

/*

   for(i=GridLength-N_BoundaryNodes;i<GridLength;i++)
 cout << i <<"  ---  " <<NewValuesX[i] << "  ---  " << NewValuesY[i]<<endl;

exit(0);
*/

   MovBoundVert[1][0]->GetCoords(x, Ay);   
    
   NewGridPos->DataToGrid();
    
   MovBoundVert[2][0]->GetCoords(x, y); // right wetting point
   y=0.;
   h_tot = x;
   h_tot /= (double)N_MovVert[0];
   for(i=1;i<N_MovVert[0];i++)
     MovBoundVert[0][i]->SetCoords(h_tot*(double)i, y);
 
  
   // axial boundary
   MovBoundVert[1][0]->GetCoords(x, y);
   
   N=N_MovVert[1];      
   h_tot = (y-Ay)/(double)N;   
   N--;
   
   for(i=1;i<N;i++)
    {
     MovBoundVert[1][i]->GetCoords(x, y);
//      cout<< " y " << y <<" new y " << y +((double)(N-i))*h_tot<<endl;      
     y += ((double)(N-i))*h_tot;
     MovBoundVert[1][i]->SetCoords(x, y);   
    }       
 

  //======================================================================  
  //  Reparametrization of free surface - Begin
  //======================================================================     
   if(reparam)  
   {   
    // update the iso points and then reparmetrize the free surf
    IIso = 0;
    for(i=0;i<N_MovVert[2];i++)
    {
     cell = Free_Cells[i];
     j = IsoCellEdgeNos[1][i];
     joint = cell->GetJoint(j);
     isojoint = (TIsoBoundEdge *)joint;
     k = isojoint->GetN_Vertices();
     Vertices = isojoint->GetVertices();
     FEId = VelocitySpace->GetFE2D(IsoCellEdgeNos[0][i], cell);
     FEDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(FEId);
     m = FEDesc->GetN_JointDOF();
     if(m == k+2)
      {
       JointDOF = FEDesc->GetJointDOF(j);
       DOF =  VeloGlobalNumbers+VeloBeginIndex[IsoCellEdgeNos[0][i]];
       for(l=0;l<k;l++)
        {
         m = DOF[JointDOF[l+1]];
         Vertices[l]->GetCoords(IsoX, IsoY);
         if(TDatabase::ParamDB->P5 > 0)
          {
           un = ValuesVX[m]*Nx[IIso+N_BoundaryNodes]
                + ValuesVY[m]*Ny[IIso+N_BoundaryNodes];
              IsoX += dt*un*Nx[IIso+N_BoundaryNodes];
              IsoY += dt*un*Ny[IIso+N_BoundaryNodes];
            }
            else
            {
             IsoX += dt*ValuesVX[m];
             IsoY += dt*ValuesVY[m];
            }
            
           if(IsoY<=0) IsoY = 1e-5;
//            if(fabs(x)<1e-12) x = 0.;
           
           Vertices[l]->SetCoords(IsoX, IsoY);
           IIso++;
          } // endfor l
        }
       else
        {
          // approximation order of isoparametric boundary and velocity
          // element must be the same
          Error("No match in isoparametric case" << endl);
          exit(-1);
        } //  if(m == k+2)     
    }// for(i=0;i<N_MovVert[2];i++)     
//   if(fabs(TDatabase::TimeDB->CURRENTTIME - 0.00220711)<1e-8)
// {
//    cout << "test ReParam_axial3D_U start " << reparam << endl;
// //   exit(0);
// }  
    ReParam_axial3D_U(N_MovVert[2], Free_Cells,  IsoCellEdgeNos[1], IsoCellEdgeNos[0],
                      Velocity, NULL, NULL, TRUE);   
//    if(fabs(TDatabase::TimeDB->CURRENTTIME - 0.00220711)<1e-8)
// {
//    cout << "test ReParam_axial3D_U end " << reparam << endl;
//   exit(0);
// }  
    OutPut("ReParam CURRENT TIME: ");
    OutPut(TDatabase::TimeDB->CURRENTTIME << endl);   
   }  //if (reparam)  
   
   NewGridPos->GridToData();
   GridPos->DataToGrid();

   memset(Rhs, 0, 2*GridLength*SizeOfDouble);
   memcpy(d, NewValuesX, 2*GridLength*SizeOfDouble);
   Daxpy(2*GridLength, -1, ValuesX, d);
   memcpy(Rhs + (GridLength-N_BoundaryNodes),
          d+(GridLength-N_BoundaryNodes),
          N_BoundaryNodes*SizeOfDouble);
   memcpy(Rhs + (2*GridLength-N_BoundaryNodes),
          d+(2*GridLength-N_BoundaryNodes),
          N_BoundaryNodes*SizeOfDouble);

  memset(Sol, 0 , 2*GridLength*SizeOfDouble);
  memcpy(Sol + (GridLength-N_BoundaryNodes),
         d+(GridLength-N_BoundaryNodes),
         N_BoundaryNodes*SizeOfDouble);
  memcpy(Sol + (2*GridLength-N_BoundaryNodes),
         d+(2*GridLength-N_BoundaryNodes),
         N_BoundaryNodes*SizeOfDouble);


  SolveGridEquation(Entries, Sol, Rhs, KCol, RowPtr, GridLength);

  memcpy(d, ValuesX, 2*GridLength*SizeOfDouble);
  Daxpy(2*GridLength, 1, Sol, d);
  memcpy(NewValuesX, d, (GridLength-N_BoundaryNodes)*SizeOfDouble);
  memcpy(NewValuesY, d+GridLength, (GridLength-N_BoundaryNodes)*SizeOfDouble);

// for(i=GridLength-N_BoundaryNodes;i<GridLength;i++)
//  cout << i << "  ---  "<<NewValuesX[i] << "  ---  " << NewValuesX[i+GridLength] << endl;
//cout << "test " << endl;
  // put solution into grid position
  IIso = 0;
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
          cell->GetVertex(j)->SetCoords(NewValuesX[k], NewValuesY[k]); 
        }
      break;

      case 4:
        k = DOF[0];
        cell->GetVertex(0)->SetCoords(NewValuesX[k], NewValuesY[k]);
        k = DOF[1];
        cell->GetVertex(1)->SetCoords(NewValuesX[k], NewValuesY[k]);
        k = DOF[3];
        cell->GetVertex(2)->SetCoords(NewValuesX[k], NewValuesY[k]);
        k = DOF[2];
        cell->GetVertex(3)->SetCoords(NewValuesX[k], NewValuesY[k]);
      break;
    } // endswitch
 
  if (!reparam)  
   {  
    N_Edges = cell->GetN_Edges();

    for(j=0;j<N_Edges;j++)
    {
      joint = cell->GetJoint(j);
      if(joint->GetType() == IsoBoundEdge)
      {
        isojoint = (TIsoBoundEdge *)joint;
        k = isojoint->GetN_Vertices();
        Vertices = isojoint->GetVertices();
        FEId = VelocitySpace->GetFE2D(i, cell);
        FEDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(FEId);
        m = FEDesc->GetN_JointDOF();
        if(m == k+2)
        {
          JointDOF = FEDesc->GetJointDOF(j);
          DOF =  VeloGlobalNumbers+VeloBeginIndex[i];
          for(l=0;l<k;l++)
          {
            m = DOF[JointDOF[l+1]];
            Vertices[l]->GetCoords(IsoX, IsoY);
            if(TDatabase::ParamDB->P5 > 0)
            {
              un = ValuesVX[m]*Nx[IIso+N_BoundaryNodes]
                  + ValuesVY[m]*Ny[IIso+N_BoundaryNodes];
              IsoX += dt*un*Nx[IIso+N_BoundaryNodes];
              IsoY += dt*un*Ny[IIso+N_BoundaryNodes];
              // cout << "U:   " << ValuesVX[m] << " " << ValuesVY[m] << endl;
              // cout << "N:   " << Nx[IIso+N_BoundaryNodes] << " "
              //                 << Ny[IIso+N_BoundaryNodes] << endl;
              // cout << "UNN: " << un*Nx[IIso+N_BoundaryNodes] << " "
              //                 << un*Ny[IIso+N_BoundaryNodes] << endl;
            }
            else
            {  
              IsoX += dt*ValuesVX[m];
              IsoY += dt*ValuesVY[m];
            }
           if(IsoY<=0) IsoY = 1e-5;
//            if(fabs(IsoX)<1e-12) IsoX = 0.;
           Vertices[l]->SetCoords(IsoX, IsoY);
            IIso++;
          } // endfor l
        }
       else
        {
          // approximation order of isoparametric boundary and velocity
          // element must be the same
          Error("No match in isoparametric case" << endl);
          exit(-1);
        }
      } // endif
    } // endfor j  
   } // if(reparam)
  } // endfor i

   if(reparam)
    {
     N_ReParam++; 
     reparam = FALSE;
    }    
    
//   delete [] d;
  
  if(TDatabase::ParamDB->P5 > 0)
  { delete [] Nx;  delete [] Ny;}
  
} // MoveGrid_imping



void MapNSESol(TFESpace2D *velocity_space_nse, double *sol_nse, int N_U_nse, TFESpace2D *pressure_space_nse, int N_P_nse,
	       TFESpace2D *velocity_space,  double *sol, int N_U, TFESpace2D *pressure_space, int N_P, int *NSE_GlobalCllNo)
{
  int i, j, k, l, N, N_Cells_NSE, n1, n2;
  int N_DOF, *UGlobalNumbers_NSE, *UBeginIndex_NSE, *PGlobalNumbers_NSE, *PBeginIndex_NSE;
  int *UGlobalNumbers, *UBeginIndex, *PGlobalNumbers, *PBeginIndex;
  int *DOF_NSE, *DOF;
  
  TCollection *Coll_NSE, *Coll;  
  TBaseCell *Me;
  FE2D FeId;
  TFEDesc2D *FeDesc;  
  
  
  Coll_NSE = velocity_space_nse->GetCollection();
  N_Cells_NSE = Coll_NSE->GetN_Cells();
  UGlobalNumbers_NSE = velocity_space_nse->GetGlobalNumbers();
  UBeginIndex_NSE = velocity_space_nse->GetBeginIndex(); 
  PGlobalNumbers_NSE = pressure_space_nse->GetGlobalNumbers();
  PBeginIndex_NSE = pressure_space_nse->GetBeginIndex(); 
  
//   Coll = velocity_space->GetCollection();
  UGlobalNumbers = velocity_space->GetGlobalNumbers();
  UBeginIndex = velocity_space->GetBeginIndex(); 
  PGlobalNumbers = pressure_space->GetGlobalNumbers();
  PBeginIndex = pressure_space->GetBeginIndex(); 
  
  memset(sol, 0, 2*N_U + N_P*SizeOfDouble);    
  for(i=0;i<N_Cells_NSE;i++)
   {
    N = NSE_GlobalCllNo[i];
    Me = Coll_NSE->GetCell(i); 
    
    //velocity
    FeId = velocity_space_nse->GetFE2D(i, Me);
    FeDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(FeId);
    N_DOF = FeDesc->GetN_DOF(); 
    
    DOF_NSE = UGlobalNumbers_NSE + UBeginIndex_NSE[i];
    DOF = UGlobalNumbers + UBeginIndex[N];  
    
     for(j=0;j<N_DOF;j++)
      {
       n1 = DOF_NSE[j];
       n2 = DOF[j];
       sol[n2] = sol_nse[n1];
       sol[N_U + n2] = sol_nse[N_U_nse + n1];      
      }
  
    //pressure
    FeId = pressure_space_nse->GetFE2D(i, Me);
    FeDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(FeId);
    N_DOF = FeDesc->GetN_DOF(); 
    DOF_NSE = PGlobalNumbers_NSE + PBeginIndex_NSE[i];
    DOF = PGlobalNumbers + PBeginIndex[N]; 
     for(j=0;j<N_DOF;j++)
      {
       n1 = DOF_NSE[j];
       n2 = DOF[j];
       sol[2*N_U + n2] = sol_nse[2*N_U_nse + n1];   
      }  
      
//      exit(0);  
   } 
//    exit(0);
}

void  GetMovingBoundData(TCollection *coll, int *N_MovVert, TBoundEdge *** &Bound_Joint, TVertex *** &MovBoundVert,
                         TIsoBoundEdge ** &Free_Joint, TBaseCell ** &Free_Cells, int ** &IsoCellEdgeNos, double x, double y)
{
 int i, j, k, l, m0, m1, m2, N_Cells, comp;
 int ORDER, VSP;
 
 double  TX[2], TY[2], x0, y0, x1, y1;
  
 TBaseCell *Me;
 TJoint *Joint;
 TBoundEdge *Solid_Joint;
 TBoundComp *BoundComp;  
 TVertex *temp_Mov;
 TBoundEdge *tempSlip_Joint;

  ORDER = 0;
  VSP = TDatabase::ParamDB->VELOCITY_SPACE;

  if (abs(VSP) > 20)
   {ORDER = abs(VSP) - 20;}
  else if ( abs(VSP) > 10)
    {ORDER = abs(VSP) - 10;}
  else ORDER = abs(VSP);
  
  
  N_Cells = coll->GetN_Cells();
  N_MovVert[0] = 0;
  N_MovVert[1] = 0;
  N_MovVert[2] = 0;
  
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
             
             if(comp==0)
              { N_MovVert[0]++; }  
             else if(comp==2)
              {N_MovVert[1]++; }
             else
              {
               cout<<"Error comp " <<endl;
               exit(0);      
              }                            
            }
           else if(Joint->GetType() == IsoBoundEdge)
            { N_MovVert[2]++;  }
            
          }// endfor l
        }// endfor j
  
//      for(i=0; i<3; i++)
//       cout<<"BDComp " << i << " N_Vert " << N_MovVert[i] << endl; 
     

     // solid bound
     Bound_Joint[0] = new TBoundEdge* [N_MovVert[0]];
     Bound_Joint[1] = new TBoundEdge* [N_MovVert[1]];

     MovBoundVert[0] = new TVertex* [N_MovVert[0]];
     MovBoundVert[1] = new TVertex*[N_MovVert[1]];

     // free bound
     Free_Joint = new TIsoBoundEdge*[N_MovVert[2]];
     MovBoundVert[2] = new TVertex*[N_MovVert[2]];
     Free_Cells = new TBaseCell*[N_MovVert[2]];
     IsoCellEdgeNos[0]  = new int [N_MovVert[2]];
     IsoCellEdgeNos[1]  = new int [N_MovVert[2]];
     
      m0 = 0;
      m1 = 0;
      m2 = 0;  
  
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
           if(comp==0)
           {
             Bound_Joint[0][m0] = (TBoundEdge *)Joint;
             MovBoundVert[0][m0] = Me->GetVertex(l);
             m0++;
            }
         else if(comp==2)
            {
             Bound_Joint[1][m1] = (TBoundEdge *)Joint;
             MovBoundVert[1][m1] = Me->GetVertex(l);
             m1++;
            }
         
           }
         else if(Joint->GetType() == IsoBoundEdge)
           {
            Free_Joint[m2] = (TIsoBoundEdge *)Joint;
            MovBoundVert[2][m2] = Me->GetVertex(l);
            Free_Cells[m2] = Me;
            IsoCellEdgeNos[0][m2] = j;
            IsoCellEdgeNos[1][m2] = l;
	    
            Me->GetVertex(l)->GetCoords(TX[0], TY[0]);
            Me->GetVertex((l+1) % k)->GetCoords(TX[1], TY[1]);
            Free_Joint[m2]->GeneratemidVert(ORDER-1, TX, TY);
            m2++;
           }
          }// endfor l
         }// endfor j    
         
         
  // sort         
    for(k=0;k<N_MovVert[0]-1;k++)
      {
      for(l=k+1;l<N_MovVert[0];l++)
       {MovBoundVert[0][k]->GetCoords(x0, y0);
	MovBoundVert[0][l]->GetCoords(x1, y1);
	if(x0 > x1)
	{
	  temp_Mov = MovBoundVert[0][k];
          MovBoundVert[0][k] = MovBoundVert[0][l];
          MovBoundVert[0][l] = temp_Mov;

	  tempSlip_Joint = Bound_Joint[0][k];
	  Bound_Joint[0][k] = Bound_Joint[0][l];
	  Bound_Joint[0][l] = tempSlip_Joint;
	 }
        }
       }         
//  for (k=0;k<N_MovVert[0];k++)
//       {
//        MovBoundVert[0][k]->GetCoords(x0, y0);
//        cout<< " x0 " << x0<<" y0 " << y0<<endl;
//        }

    for(k=0;k<N_MovVert[1]-1;k++)
     {
      for(l=k+1;l<N_MovVert[1];l++)
       {MovBoundVert[1][k]->GetCoords(x0, y0);
	MovBoundVert[1][l]->GetCoords(x1, y1);
	if(y1 > y0)
	{
	  temp_Mov = MovBoundVert[1][k];
          MovBoundVert[1][k] = MovBoundVert[1][l];
          MovBoundVert[1][l] = temp_Mov;

	  tempSlip_Joint = Bound_Joint[1][k];
	  Bound_Joint[1][k] = Bound_Joint[1][l];
	  Bound_Joint[1][l] = tempSlip_Joint;
	 }
        }
       }
//  for (k=0;k<N_MovVert[1];k++)
//       {
//        MovBoundVert[1][k]->GetCoords(x0, y0);
//        cout<< " x " << x0<<" y " << y0<<endl;
//        }

   Sort_Imping(Free_Cells, MovBoundVert[2], IsoCellEdgeNos[0], IsoCellEdgeNos[1], N_MovVert[2], x, y);

}// GetMovingBoundData


void RemeshAxial3D_ImpDrop(TDomain * &Domain, TFESpace2D ** &FESpaces_All,
                           TFEVectFunct2D ** &FEVectFuncts_All, TFEFunction2D ** &FEFunctions_All, int *N_MovVert, 
                           TBoundEdge *** &Bound_Joint, TVertex *** &MovBoundVert, TIsoBoundEdge ** &Free_Joint,
                           TBaseCell ** &Free_Cells, int ** &IsoCellEdgeNos, double ** &Sol_All,  double ** &Rhs_All, 
                           TSquareStructure2D ** &SquareStructure_All, TStructure2D ** &Structure_All,
                           TSquareMatrix2D ** &SqMat_All, TMatrix2D ** &Mat_All)
{
  int i, j, k, l, N_G, N_Cells, ORDER, VSP, N_DOF, N_ThermalDOF;
  int In_Index, CurrComp, Old_N_Cells, Old_N_RootCells, CurrVertex, N_Joints, N_Vertices, ID;
  int N_RootCells, *Triangles, *PointNeighb, maxEpV = 0, a, b, Neighb_tmp, Neib[2];
  int CurrNeib, len1, len2, pressure_space_code, N_U, N_P, N_Unknowns, comp;
  int *JointDOF, *DOF, *GlobalNumbers, *BeginIndex, N_refX, N;
  
  double d, h, t, tx, ty, x1, x2, y1, y2, Lx, Ly, Rx, Ry, *S_BX, *S_BY, *A_BX, *A_BY, refX;
  double area, *Coordinates, left, right, bottom, top, T_a, T_b, *sol, *ValuesU2, *Sx, *Sy;
  
  TBoundPart *BoundPart;
  TBoundComp *BoundComp;
  TBdLine *UpdateSlipBound, *UpdateAxialBound;
  TCollection *coll, *Old_Coll, *mortarcoll = NULL;
  TBaseCell *cell, **CellTree, **Old_CellTree;
  TVertex **VertexDel, **NewVertices;
  TJoint *Joint;
  TFESpace2D *velocity_space, *pressure_space, *thermal_space;
  TFEVectFunct2D *u; 
  TFEFunction2D *p; 
  TBoundEdge *Solid_Joint;
  FE2D FeId;
  TFEDesc2D *FeDesc; 
    
  // strings
  char ReadinDat[] = "readin.dat";
  char TString[] = "T";
  char NameString[]  = "name";
  char UString[] = "u";
  char PString[] = "p";
  char WString[] = "w";
  
  std::ostringstream opts;
  std::ostringstream os;
  os << " ";
  opts << " ";
// 
  struct triangulateio In, Out;
  
//   // free surface vertices
//   d = 0;
//   for(i=0;i<N_MovVert[2]-1;i++) // without last point
//    {
//     MovBoundVert[2][i]->GetCoords(x1, y1);
//     MovBoundVert[2][i+1]->GetCoords(x2, y2);
//      d +=sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2)  );
//    }
//    h = d/((double)i-1.);
  h = TDatabase::ParamDB->P9;
  
  MovBoundVert[0][0]->GetCoords(Lx, Ly);
  MovBoundVert[2][0]->GetCoords(Rx, Ry);
  
  
  d = Rx-Lx;
  k = (int)(d/h); // No of intervals with step length h
  if(k<2) k=2;     // minimum two intervals
  t = d/(double)k;
  N_MovVert[0] = k;
  S_BX = new double[N_MovVert[0]];
  S_BY = new double[N_MovVert[0]];
     
  for(i=0;i<N_MovVert[0];i++)
   {
    S_BX[i]= Lx + (double)i*t;
    S_BY[i]=0.0;
//      cout<<i<< " x :" << S_BX[i] << " -----------------y: " <<S_BY[i]<< endl;
   }
 
  //axial bound
  MovBoundVert[1][0]->GetCoords(x1, y1);
  MovBoundVert[0][0]->GetCoords(x2, y2);
  
//   if(y1>0.2)
//    {
//     refX = 0.2;  
//     N_MovVert[1] = 10;    
//    }
//   else
//    {
//     refX = y1/2.0;   
//     N_MovVert[1] = 3;   
//    }   
//   
//   N_refX = (int)(refX/(1.5*h));
//   N_MovVert[1] += N_refX;
//   A_BY = new double[N_MovVert[1]];
// 
//    N=0;
//    t= -refX / (double) N_refX ;
//    for(i=0;i< N_refX;i++)  
//     A_BY[N++]= y1 + (double)i*t;
//    
//    y1 -= refX;
//    t= -y1/(double)(N_MovVert[1] - N_refX);
//     
//    j = N_MovVert[1] - N_refX;
//    for(i=0;i< j;i++)  
//     A_BY[N++]= y1 + (double)i*t;   

   d  = y2 - y1; 
   k = (int)(fabs(d)/h);
   if(k<2) k=2;     // minimum two intervals     
   t = d/(double)k;
   N_MovVert[1] = k;
   A_BY = new double[N_MovVert[1]]; 

  for(i=0;i<N_MovVert[1];i++)
  {
   A_BY[i] = y1 + t*(double)i;
   // cout<<i<< " x :" << 0 << " -----------------y: " <<A_BY[i]<< endl;
  }
//   exit(0);
  
  area = TDatabase::ParamDB->Area; 
//======================================================================
// Triangular for grid generation begin
//======================================================================
  BoundPart = Domain->GetBdPart(0);
  UpdateSlipBound = (TBdLine*)BoundPart->GetBdComp(0);
  UpdateAxialBound = (TBdLine*)BoundPart->GetBdComp(2);

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
  
  In.numberofpoints = N_MovVert[0]+N_MovVert[1]+N_MovVert[2];
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

//   Sx = new double[N_MovVert[2]+1];
//   Sy = new double[N_MovVert[2]+1];
 
//   for(i=0;i<N_MovVert[2];i++) // without last point  
//    MovBoundVert[2][i]->GetCoords(Sx[i], Sy[i]);   
  
//   if(TDatabase::ParamDB->P10==0)
//    {
//     Sx[N_MovVert[2]] = A_BX[0];
//     Sy[N_MovVert[2]] = A_BY[0]; 
//     EqDist_Pts(N_MovVert[2]+1, Sx, Sy);         
//    }

  // points and segments on the free boundary (marker=2)
  for(i=0;i<N_MovVert[2];i++) // without last point
   {
    MovBoundVert[2][i]->GetCoords(tx, ty);
    In.pointlist[2*In_Index] = tx;
    In.pointlist[2*In_Index+1] = ty;
//    cout<<In_Index<< " x :" << tx<< " -----------------y: " <<ty<< endl;
    In.pointmarkerlist[In_Index] = CurrComp;
    In.segmentlist[2*In_Index] = In_Index;
    In.segmentlist[2*In_Index+1] = In_Index+1;
    In.segmentmarkerlist[In_Index] = CurrComp;
    In_Index++;
   }
//    delete [] Sx; delete [] Sy; 

   CurrComp++;
  // points and segments on the solid boundary (marker=1)
  for(i=0;i<N_MovVert[1];i++) // without last point
   {
    In.pointlist[2*In_Index] = 0.;
    In.pointlist[2*In_Index+1] = A_BY[i];
//     cout<<In_Index<< " x :" << A_BX[i]<< " -----------------y: " <<A_BY[i]<< endl;
    In.pointmarkerlist[In_Index] = CurrComp;
    In.segmentlist[2*In_Index] = In_Index;
    In.segmentlist[2*In_Index+1] = In_Index+1;
    In.segmentmarkerlist[In_Index] = CurrComp;
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
    
  triangulate((char*)opts.str().c_str(), &In, &Out, (struct triangulateio *)NULL);    
 
  
  Old_Coll = Domain->GetCollection(It_Finest, 0);
  Old_N_Cells = Old_Coll->GetN_Cells();
  Domain->GetTreeInfo(Old_CellTree, Old_N_RootCells);  
  if(Old_N_Cells != Old_N_RootCells) 
   exit(-1);
  
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
  Coordinates = Out.pointlist;
  Triangles = Out.trianglelist;
  
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
  
  
  Domain->SetBoundBox(right-left,top-bottom);
  Domain->SetBoundBoxstart(left,bottom); 
  
  UpdateSlipBound->SetParams(Lx, Ly, Rx-Lx, Ry-Ly);
  UpdateAxialBound->SetParams(0., A_BY[0], 0., Ly-A_BY[0]);  

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
  
  for(i=0;i<3*N_RootCells;i++)
  {
    j = Triangles[i]*maxEpV;
    PointNeighb[j]++;
    PointNeighb[j + PointNeighb[j]] = i/3;
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

   if(Out.edgemarkerlist[i]) // 0 for inner edges and Boundcomp+1 for Boundedge respect
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
          cout << " CurrComp " << CurrComp <<endl;
        //  exit(0);
       }

      if (CurrNeib == 2)    // 2 cells contain the current edge
        if(Domain->GetBdPart(0)->GetBdComp(CurrComp)->IsFreeBoundary())
          Joint = new TIsoInterfaceJoint(Domain->GetBdPart(0)->GetBdComp(CurrComp), T_a, T_b, 
                                         CellTree[Neib[0]], CellTree[Neib[1]]);
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
  
  } //  for (i=0;i<N_G;i++)

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
  
// ======================================================================
// Triangular for grid generation end
// ======================================================================
 
 
   // refine grid up to the coarsest level
  for(i=0;i<TDatabase::ParamDB->UNIFORM_STEPS;i++)
   Domain->RegRefineAll(); 
 
 
//       // write grid into an Postscript file
//       os.seekp(std::ios::beg);
//       os << "Domain" << ".ps" << ends;
//       Domain->PS(os.str().c_str(),It_Finest,0);
// //      exit(0);
  
//======================================================================
// construct all finite element spaces
//======================================================================
  coll=Domain->GetCollection(It_Finest, 0);
  N_Cells = coll->GetN_Cells();
 
  ORDER = 0;
  VSP = TDatabase::ParamDB->VELOCITY_SPACE;

  if (abs(VSP) > 20)
   {ORDER = abs(VSP) - 20;}
  else if ( abs(VSP) > 10)
    {ORDER = abs(VSP) - 10;}
  else ORDER = abs(VSP);

//======================================================================
// construct all finite element spaces
//======================================================================
  // get velocity and pressure spacess
  GetVelocityAndPressureSpace(coll,BoundCondition,
                              mortarcoll, velocity_space,
                              pressure_space, &pressure_space_code,
                              TDatabase::ParamDB->VELOCITY_SPACE,
                              TDatabase::ParamDB->PRESSURE_SPACE);
  N_U = velocity_space->GetN_DegreesOfFreedom();
  N_P = pressure_space->GetN_DegreesOfFreedom();
  GlobalNumbers = velocity_space->GetGlobalNumbers();
  BeginIndex = velocity_space->GetBeginIndex();

  // mesh velocity space 
  delete FESpaces_All[2];
  FESpaces_All[2] = new TFESpace2D(coll, NameString, WString, GridBoundCondition, 1, NULL);
  N_G = FESpaces_All[2]->GetN_DegreesOfFreedom();
  
  // thermal space
#ifdef __ENERGY__    
  delete FESpaces_All[3];
  FESpaces_All[3] = new TFESpace2D(coll, NameString, TString, HeatBoundCondition,
                                  TDatabase::ParamDB->ANSATZ_ORDER, NULL); 
  N_ThermalDOF = FESpaces_All[3]->GetN_DegreesOfFreedom();
#endif  
//===========================================================================
  delete [] Bound_Joint[0];
  delete [] Bound_Joint[1];                      

  delete [] MovBoundVert[0];
  delete [] MovBoundVert[1];   
  delete [] MovBoundVert[2];  

  delete [] Free_Joint;    
  delete [] Free_Cells;  
  delete [] IsoCellEdgeNos[0];  
  delete [] IsoCellEdgeNos[1];  
    
  Ry = 0.;
  GetMovingBoundData(coll, N_MovVert, Bound_Joint, MovBoundVert, Free_Joint,
                     Free_Cells, IsoCellEdgeNos, Rx, Ry);
     
//======================================================================
// construct all finite element functions
//======================================================================
  N_Unknowns = 2*N_U + N_P;
  
  OutPut("dof velocity : "<< setw(10) << 2* N_U << endl);
  OutPut("dof pressure : "<< setw(10) << N_P << endl);
  OutPut("dof all      : "<<  setw(10) << N_Unknowns  << endl);  
  
  delete [] Rhs_All[0];
  Rhs_All[0] = new double[N_Unknowns];
  sol = new double[N_Unknowns];
  memset(sol, 0, N_Unknowns*SizeOfDouble);
   
  u =  new TFEVectFunct2D(velocity_space, UString, UString, sol, N_U, 2);
  p = new TFEFunction2D(pressure_space, PString,  PString,  sol+2*N_U, N_P);

  u->Interpolate(FEVectFuncts_All[0]);
  
//====================================================================== 
// impose no-penetration condition 
//====================================================================== 
  ValuesU2 = sol + N_U;
  
  for(i=0;i<N_Cells;i++)
   {
    cell = coll->GetCell(i);
    k = cell->GetN_Edges();
    
    for(l=0;l<k;l++)
     {
      Joint = cell->GetJoint(l);    
      
      if(Joint->GetType() == BoundaryEdge)
        {
         FeId = velocity_space->GetFE2D(i, cell);  
         FeDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(FeId);
         JointDOF = FeDesc->GetJointDOF(l); 
         N_DOF = FeDesc->GetN_JointDOF();
         DOF = GlobalNumbers+BeginIndex[i];
    
         Solid_Joint = (TBoundEdge *)Joint;
         BoundComp = Solid_Joint->GetBoundComp();
         comp=BoundComp->GetID();

         if(comp==0) 
          {           
           for(j=0;j<N_DOF;j++)
            ValuesU2[DOF[JointDOF[j]]] = 0.;
          }  
         else if(comp==2)
          {
           for(j=0;j<N_DOF;j++)
            sol[DOF[JointDOF[j]]] = 0.;  
          }
 
        }// if(Joint->GetType() == Bo
     } // for(l=0;l<k;l++)
   } //  for(i=0;i<N_Cells;i++)
//======================================================================   
// no need to interpolate pressure for DirectSolvers !!!
//======================================================================    
  delete FESpaces_All[0];
  delete FESpaces_All[1]; 
 
  FESpaces_All[0] = velocity_space;
  FESpaces_All[1] = pressure_space;    

  delete [] Sol_All[0];
  Sol_All[0] = sol;
 
  delete FEVectFuncts_All[0]; 
  FEVectFuncts_All[0] = u;
 
  FEFunctions_All[0] = FEVectFuncts_All[0]->GetComponent(0);
  FEFunctions_All[1] = FEVectFuncts_All[0]->GetComponent(1); 

  delete FEFunctions_All[2];
  FEFunctions_All[2] = p;
  
//======================================================================
// grid space finite element functions
//====================================================================== 
  delete [] Sol_All[1];
  Sol_All[1] = new double[2*N_G];
  delete [] Rhs_All[1]; 
  Rhs_All[1] = new double[2*N_G];  
  
  memset(Sol_All[1], 0, 2*N_G*SizeOfDouble); 
  delete FEVectFuncts_All[1];
  FEVectFuncts_All[1]  = new TFEVectFunct2D(FESpaces_All[2], WString, WString, Sol_All[1], N_G, 2);
  
  FEFunctions_All[3] = FEVectFuncts_All[1]->GetComponent(0);
  FEFunctions_All[4] = FEVectFuncts_All[1]->GetComponent(1); 

    
//======================================================================
// thermal space finite element functions
//======================================================================    
#ifdef __ENERGY__  
  delete [] Rhs_All[2];
  Rhs_All[2] = new double[N_ThermalDOF]; 
  memset(Rhs_All[2], 0, N_ThermalDOF*SizeOfDouble);
    
  delete [] Sol_All[2];
  Sol_All[2] = new double[N_ThermalDOF];

  memset(Sol_All[2], 0, N_ThermalDOF*SizeOfDouble);
  
  // thermal fefunction
  delete FEFunctions_All[5];
  FEFunctions_All[5] = new TFEFunction2D(FESpaces_All[3], TString, TString, Sol_All[2], N_ThermalDOF);
#endif
// ======================================================================
// allocate memory for all matrices
// ======================================================================  
  delete Structure_All[0];  delete Structure_All[1];
  
  Structure_All[0] = new TStructure2D(FESpaces_All[1], FESpaces_All[0]); // B
  Structure_All[1] = new TStructure2D(FESpaces_All[0], FESpaces_All[1]); // BT

  delete SquareStructure_All[0];  delete SquareStructure_All[1]; 

//thermal
#ifdef __ENERGY_   
  delete SquareStructure_All[2];
#endif  
  
  //velo 
  SquareStructure_All[0] = new TSquareStructure2D(FESpaces_All[0]);  
  SquareStructure_All[0]->Sort();  
  
  // grid 
  SquareStructure_All[1] = new TSquareStructure2D(FESpaces_All[2]); 
  SquareStructure_All[1]->Sort();  
   
  //thermal
#ifdef __ENERGY__    
  SquareStructure_All[2] = new TSquareStructure2D(FESpaces_All[3]);
  SquareStructure_All[2]->Sort(); 
#endif
  // u
  for(i=0; i<10; i++)
   {
    delete SqMat_All[i];
    SqMat_All[i] = new TSquareMatrix2D(SquareStructure_All[0]);
   }

  // B
  for(i=0; i<2; i++)
   {
    delete Mat_All[i];
    Mat_All[i] = new TMatrix2D(Structure_All[0]);   
   }  

  // BT
  for(i=2; i<4; i++)
   {
    delete Mat_All[i];
    Mat_All[i] = new TMatrix2D(Structure_All[1]);   
   }  

  // for mesh
  for(i=10; i<14; i++)
   {
    delete SqMat_All[i];
    SqMat_All[i] = new TSquareMatrix2D(SquareStructure_All[1]);  
   }  
  
  // for heat
#ifdef __ENERGY__    
  for(i=14; i<16; i++)
   {
    delete SqMat_All[i];
    SqMat_All[i] = new TSquareMatrix2D(SquareStructure_All[2]);  
   }  
#endif   
   
// ======================================================================
// delete old mesh
// ======================================================================  
    
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

  delete [] S_BX;
  delete [] S_BY;
//   delete [] A_BX;
  delete [] A_BY; 
}// RemeshAxial3D_ImpDrop



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
  TFESpace2D *Grid_space, *vorticity_space, *thermal_space,*grid_space;

  TOutput2D *Output;
  TFEVectFunct2D *RefGridPos, *AuxGridPos, *GridPos, *ReparamPos;
  TFEFunction2D *fefct[4];
  TFESpace2D *fesp[3], *ferhs_T[3], *ferhs[2];
  TAuxParam2D *aux;
  TDiscreteForm2D *DiscreteFormMatrixT_MRhs;
  TDiscreteForm2D *DiscreteForm;
  TMatrix2D *MATRICES[4];
  TSquareMatrix2D *SQMATRICES[8], *SQMATRICES_GRID[4];
  TDiscreteForm2D *DiscreteFormGalerkin;
  TDiscreteForm2D *DiscreteFormNLGalerkin;
  TDiscreteForm2D *DiscreteFormGrid, *DiscreteFormHeat;
  TSquareMatrix **sqmatrices = (TSquareMatrix **)SQMATRICES;
  TMatrix **matrices = (TMatrix **)MATRICES;
  MatVecProc *MatVect;
  DefectProc *Defect;

  TFESpace2D **FESpaces_All = new TFESpace2D *[4];      
  TFEFunction2D **FEFunctions_All = new TFEFunction2D *[6];    
  TFEVectFunct2D **FEVectFuncts_All = new TFEVectFunct2D*[2];
  TStructure2D **Structure_All = new TStructure2D *[2];
  TSquareStructure2D **SquareStructure_All = new TSquareStructure2D *[3];
  TSquareMatrix2D **SqMat_All = new TSquareMatrix2D *[16];
  TMatrix2D **Mat_All = new TMatrix2D *[4];
  
  double total_time,*Coordinates;
  double  t, teta, dt,x,y,gamma, tx,ty,sx,sy, R_Theta[3];
  double left, right, top, bottom,T_a, T_b;
  double x0, y0,x1,y1,hi, residual, impuls_residual, oldresidual, solver_time;
  double *oldsol_T;
  double end_time, t1, t2, t4, t3;
  double *B, *defect;
  double *RHSs[3], *refpos, *auxpos, *pos, *ReparamDisp, *ReparamMeshVelo;
  double  TX[2], TY[2], solver_time_curr;
  double SLPX, SLPY, *Entries[4], tau, oldtau, limit, *sol_output, Params[10], InitVolume, CurrVolume;  
  double Lx, Ly, Rx, Ry,  x2, y2, x3, y3, x4, y4, fh, fhlimit, fhtot, fhmin, fhmax;
  double *Angle = new double[2], **FreePts = new double *[2];  
  double **Sol_All = new double *[3], *tmp_Gsol, *tmp_Gd;
  double **Rhs_All = new double *[3], MaxWetD=0., T_MaxWetD=0.;
  
  int ret,N_RootCells,N_Cells,N_Joints, N_Vertices,N_G, N_NSE_Cells, N_NonNSE_Cells;
  int N_SlipBound_Vert,  N_FreeBound_Vert,  N_AxialBound_Vert,N_Interf_Vertices;
  int In_Index,CurrComp,CurrVertex, img=1, RemeshImg=1, N_BData=1;
  int ID,CurrNeib,Neib[2], N_SquareMatrices, N_RectMatrices;
  int a,b,i,X,Y,j,k,l,len1, len2,Neighb_tmp,Indextemp;
  int *PartMarker, *Triangles,comp, *NSE_GlobalCllNo;
  int *PointNeighb,maxEpV = 0, Max_It, N_U_output, N_P_output;
  int  N_U, N_P,N_Unknowns,N_pressureDOF,N_Rhs,N_FESpaces;
  int  N_Hori1, N_Hori2,N_Verti,N_Boundary_Vert,N_thermalDOF,N_thermalActive,N_thermalNonActive;
  int velocity_space_code, pressure_space_code;
  int m,m0,m1,m2,m3,m4,m5,m6, N_Active, N_LinIterCurr, N_LinIter;
  int ORDER,VSP, N_GActive, N_GBoundaryNodes, N_SubSteps, very_first_time=1 ;
  int **IsoCellEdgeNos, *GridKCol, *GridRowPtr, *RowPtr, last_sq;
  int  *JointDOF, N_DOF, dof, *DOF, *GlobalNumbers, *BeginIndex, N_ReParam=0, N_Remesh=0;
  
  char *PRM, *GEO, *PsBaseName, *VtkBaseName;
  char ReadinDat[] = "readin.dat";
  char TString[] = "T";
  char NameString[]  = "name";
  char UString[] = "u";
  char PString[] = "p";
  char WString[] = "w";
  const char vtkdir[] = "VTK";
  const char BDdir[] = "BDData";
  
  bool remeshed=FALSE, reparam = FALSE;

  TBoundPart *BoundPart;
  TBoundComp *BoundComp;
  TBdLine *UpdateSlipBound, *UpdateAxialBound;
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

#ifdef __ENERGY__
  BoundCondFunct2D *HeatBoundaryConditions[1];
  BoundValueFunct2D *HeatBoundValues[1];
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

  mkdir(vtkdir, 0777);
  mkdir(BDdir, 0777);
#define __AXIAL3D__ 

//======================================================================
// read boundary parameterization and initialize coarse grid
//======================================================================

   Domain->Init(PRM, GEO);

   // write grid into an Postscript file
//    os.seekp(std::ios::beg);
//    os << "Domain_old" << ".ps" << ends;
//    Domain->PS(os.str().c_str(),It_Finest,0);

   boolean AllowEdgeRef = (boolean) TDatabase::ParamDB->MESHGEN_ALLOW_EDGE_REF;

//======================================================================
// Triangular for grid generation begin
//======================================================================
  BoundPart = Domain->GetBdPart(0);
  UpdateSlipBound = (TBdLine*)BoundPart->GetBdComp(0);
  UpdateFreeBound = (TBdCircle*)BoundPart->GetBdComp(1);
  UpdateAxialBound = (TBdLine*)BoundPart->GetBdComp(2);

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

  double *Sx, *Sy, area = TDatabase::ParamDB->Area;
  double  phi, r = TDatabase::ParamDB->P4, mod_r;
  double hE, refX;
  int N_refX;

   int *N_MovVert = new int[3];
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
  N_AxialBound_Vert = 50;
  N_SlipBound_Vert = 2;     // Initially only two points on solid bound (except end point)
  
  N_Interf_Vertices = N_FreeBound_Vert+N_SlipBound_Vert+N_AxialBound_Vert;
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
  
          a=0;
          b=r;
  teta = -Pi/2.; // end degree value of freesurface
  dt = (Pi/2. - teta)/(N_SlipBound_Vert+N_FreeBound_Vert);
  t = teta;
  for(i=0;i<N_SlipBound_Vert;i++) // without last point
   {
#ifdef __EXPERIMENTAL__ 
    x = r*cos(t);
    y = r*sin(t); 

    phi = atan2(y, x);
    mod_r = 1.0 + 0.29*(sqrt(5./Pi)*(3.*cos(phi+Pi/2.)*cos(phi+Pi/2.) -1. )/4.); 
    x = a+mod_r*cos(t);
    y = b+mod_r*sin(t);
    
    if(i==0)   ty = mod_r;      
#else
   mod_r = r;
   ty = r;    
#endif

    x = a + mod_r*cos(t);
    y=0;
  
 
    if(i==0) x = 0.;    
    
    In.pointlist[2*In_Index] = x;
    In.pointlist[2*In_Index+1] = y;

    In.pointmarkerlist[In_Index] = CurrComp;
    In.segmentlist[2*In_Index] = In_Index;
    In.segmentlist[2*In_Index+1] = In_Index+1;
    In.segmentmarkerlist[In_Index] = CurrComp;
    In_Index++;
    t = teta + double(i+1.)*dt;
   }
  
   CurrComp++;
//     cout<<endl; 
//      exit(0);
     
   Sx = new double[N_FreeBound_Vert+1];
   Sy = new double[N_FreeBound_Vert+1];
     
   for(i=0;i<N_FreeBound_Vert;i++) // without last point
    {
#ifdef __EXPERIMENTAL__   
     x = r*cos(t);
     y = r*sin(t);  

     phi = atan2(y, x); 
     mod_r = 1.0 + 0.29*(sqrt(5./Pi)*(3.*cos(phi+Pi/2.)*cos(phi+Pi/2.) -1. )/4.);           
#else
   mod_r = r;    
#endif   
   
      x = a + mod_r*cos(t);
      y = ty + mod_r*sin(t);       
      
     if (fabs(x)<1.e-10) x = 0;
     if(i==0) y =0.;  

      Sx[i] = x;
      Sy[i] = y;      
     
      t = teta + double(N_SlipBound_Vert+i+1)*dt;     
    }// for(i=0;i<N_FreeBound_Ver
     
   //endpoint 
   Sx[N_FreeBound_Vert] = 0;
   Sy[N_FreeBound_Vert] = 2*ty;       
          
#ifdef __EXPERIMENTAL__       
   EqDist_Pts(N_FreeBound_Vert+1, Sx, Sy); 
#endif       
     
  for(i=0;i<N_FreeBound_Vert;i++) // without last point
    {
      x = Sx[i];
      y = Sy[i];
      
     if (fabs(x)<1.e-10) x = 0;
     if(i==0) 
       {
        y =0.;   SLPX = x;
        Indextemp = In_Index;    
       }
      else if(i==1) 
      {
       hE = sqrt( (In.pointlist[2*(In_Index-1)] - x)*(In.pointlist[2*(In_Index-1)] - x) 
               + (In.pointlist[2*(In_Index-1) +1] - y)*(In.pointlist[2*(In_Index-1)+1] - y));
      }
      
      In.pointlist[2*In_Index] = x;
      In.pointlist[2*In_Index+1] = y;   
//        cout<<i<< " x : "<< x << " y : "<< y<<endl;
      In.pointmarkerlist[In_Index] = CurrComp;
      In.segmentlist[2*In_Index] = In_Index;
      In.segmentlist[2*In_Index+1] = In_Index+1;
      In.segmentmarkerlist[In_Index] = CurrComp;
      In_Index++;
    }
//    exit(0);
//   cout<<endl; 
   CurrComp++;
   dt=-2.*ty/N_AxialBound_Vert;
   t = -teta;
   
#ifdef __EXPERIMENTAL__  
    phi = Pi/2;     
    mod_r = 1.0 + 0.29*(sqrt(5./Pi)*(3.*cos(phi+Pi/2.)*cos(phi+Pi/2.) -1. )/4.); 
#else
    mod_r = r;       
#endif    

    //needed during remeshing
    TDatabase::ParamDB->P9 = hE;

//    refX = 0.2;   
//    N_refX = (int)(refX/(1.5*hE));
//    if( (N_refX+10)>N_AxialBound_Vert)
//    {
//     cout << "Increase  N_AxialBound_Vert points, N_refX:  "<<  N_refX <<endl;
//     exit(0);
//    }
// 
//    dt= -0.2 / (double)N_refX;
//    x = 0.;
//    y = ty+mod_r*sin(t);   
//    t = y;
//    
//    
//   UpdateAxialBound->SetParams(x, y, 0, -y);
//  
//   for(i=0;i<N_refX;i++) // without last point
//    {
//       if (fabs(y)<1e-12) y = 0.;
//       In.pointlist[2*In_Index] = x;
//       In.pointlist[2*In_Index+1] = y;
// //       cout<<" x : "<< x << " y : "<< y<<endl;
//       In.pointmarkerlist[In_Index] = CurrComp;
//       In.segmentlist[2*In_Index] = In_Index;
//       In.segmentlist[2*In_Index+1] = In_Index+1;
//       In.segmentmarkerlist[In_Index] = CurrComp;
//       In_Index++;
//       y = t + double(i+1)*dt;
//     } 
//     
//     
// //    cout<<endl;
//    
//    dt= -((ty+mod_r*sin(t)) -  refX) / (double)(N_AxialBound_Vert - N_refX);
//    x = 0.;
//    y = (ty+mod_r*sin(t)) -  refX;
//    t = y;   
//    j = N_AxialBound_Vert - N_refX;       
// 
//     for(i=0;i<j;i++) // without last point
//      {
//       if (fabs(y)<1e-12) y = 0.;
//       In.pointlist[2*In_Index] = x;
//       In.pointlist[2*In_Index+1] = y;
// //       cout<<" x : "<< x << " y : "<< y<<endl;
//       In.pointmarkerlist[In_Index] = CurrComp;
//       In.segmentlist[2*In_Index] = In_Index;
//       In.segmentlist[2*In_Index+1] = In_Index+1;
//       In.segmentmarkerlist[In_Index] = CurrComp;
//       In_Index++;
//       y = t + double(i+1)*dt;
//     }

   y = ty+mod_r*sin(t);
   dt= -y / (double)(N_AxialBound_Vert);
   x = 0.;
   t = y;   

   
   UpdateAxialBound->SetParams(x, y, 0, -y);  
    for(i=0;i<N_AxialBound_Vert;i++) // without last point
     {
      if (fabs(y)<1e-12) y = 0.;
      In.pointlist[2*In_Index] = x;
      In.pointlist[2*In_Index+1] = y;
//       cout<<" x : "<< x << " y : "<< y<<endl;
      In.pointmarkerlist[In_Index] = CurrComp;
      In.segmentlist[2*In_Index] = In_Index;
      In.segmentlist[2*In_Index+1] = In_Index+1;
      In.segmentmarkerlist[In_Index] = CurrComp;
      In_Index++;
      y = t + double(i+1)*dt;
    }   

  In.segmentlist[2*(In_Index-1)+1] = 0;
 
 delete [] Sx;
 delete [] Sy; 
//   exit(0);

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

//  OutPut("left: "<<left<<" right: "<<right<<" top: "<<top<<" bottom: "<<bottom<<endl);

  Domain->SetBoundBox(right-left,top-bottom);
  Domain->SetBoundBoxstart(left,bottom);

  cout << "N_RootCells :" << N_RootCells << endl;
 
  
  
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
          cout << " CurrComp " << CurrComp <<endl;
        //  exit(0);
       }

      if (CurrNeib == 2)    // 2 cells contain the current edge
        if(Domain->GetBdPart(0)->GetBdComp(CurrComp)->IsFreeBoundary())
          Joint = new TIsoInterfaceJoint(Domain->GetBdPart(0)->GetBdComp(CurrComp), T_a, T_b, 
                                         CellTree[Neib[0]], CellTree[Neib[1]]);
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

//      exit(0);
      // write grid into an Postscript file
//       os.seekp(std::ios::beg);
//       os << "Domain" << ".ps" << ends;
//       Domain->PS(os.str().c_str(),It_Finest,0);
//      exit(0);

   //Initialize DiscreteForms         
  InitializeDiscreteForms_Moving(DiscreteFormGalerkin, DiscreteFormNLGalerkin,
                                 DiscreteFormGrid, LinCoeffs, GridCoeffs);
 
  BoundaryConditions[0] = BoundCondition;
  BoundaryConditions[1] = BoundCondition;

  BoundValues[0] = U1BoundValue;
  BoundValues[1] = U2BoundValue;

  GridBoundaryConditions[0] = GridBoundCondition;
  GridBoundValues[0] = GridBoundValue;

#ifdef __ENERGY__
  HeatBoundaryConditions[0] = HeatBoundCondition;
  HeatBoundValues[0] = TBoundValue;
#endif

//======================================================================
// construct all finite element spaces
//======================================================================
  coll=Domain->GetCollection(It_Finest, 0);
  N_Cells = coll->GetN_Cells();
 
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
   grid_space = new TFESpace2D(coll, NameString, TString, GridBoundCondition, 1, NULL);
   FESpaces_All[2] =  grid_space;     
   N_G = FESpaces_All[2]->GetN_DegreesOfFreedom();
   N_GActive = FESpaces_All[2]->GetActiveBound();
   N_GBoundaryNodes = N_G - N_GActive;
   
   OutPut("N_G    : "<< setw(10) << N_G  << endl);
   OutPut("N_GActive    : "<< setw(10) << N_GActive  << endl);
   
// thermal space
#ifdef __ENERGY__
   thermal_space = new TFESpace2D(coll, NameString, TString, HeatBoundCondition,
                                  TDatabase::ParamDB->ANSATZ_ORDER, NULL);
   FESpaces_All[3] =  thermal_space;  
   N_thermalDOF = FESpaces_All[3]->GetN_DegreesOfFreedom();
   N_thermalActive = FESpaces_All[3]->GetActiveBound();
   N_thermalNonActive = N_thermalDOF - N_thermalActive;
   OutPut("N_thermalDOF    : "<< setw(10) << N_thermalDOF  << endl);
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
// thermal space finite element functions
//======================================================================
#ifdef __ENERGY__
  Sol_All[2] = new double[N_thermalDOF];
  Rhs_All[2] = new double[N_thermalDOF];
  oldsol_T = new double[N_thermalDOF];
  
  memset(Sol_All[2], 0, N_thermalDOF*SizeOfDouble);
  memset(oldsol_T, 0, N_thermalDOF*SizeOfDouble);
  memset(Rhs_All[2], 0, N_thermalDOF*SizeOfDouble);

  // thermal fefunction
  FEFunctions_All[5] = new TFEFunction2D(FESpaces_All[3], TString, TString, Sol_All[2], N_thermalDOF);
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
  
  //thermal
#ifdef __ENERGY__
  SquareStructure_All[2] = new TSquareStructure2D(FESpaces_All[3]);
  SquareStructure_All[2]->Sort();
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
#ifdef __ENERGY__  
  SqMat_All[14]  = new TSquareMatrix2D(SquareStructure_All[2]); // T_M
  SqMat_All[15] = new TSquareMatrix2D(SquareStructure_All[2]); // T_A
#endif     
  IsoCellEdgeNos = new int *[2];
 
 //  ====================================================================================   
  y = 0.;
  GetMovingBoundData(coll, N_MovVert, Bound_Joint, MovBoundVert, Free_Joint,
                     Free_Cells, IsoCellEdgeNos, SLPX, y);
                                          
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
#ifdef __ENERGY__      
  Output->AddFEFunction(FEFunctions_All[5]);
#endif      
//   Output->AddFEVectFunct(FEVectFuncts_All[1]);  
      
  os.seekp(std::ios::beg);
  Output->AddParameter(TDatabase::TimeDB->CURRENTTIME,os.str().c_str());      
 
     
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
   MovBoundVert[0][0]->GetCoords(Lx, Ly);
   MovBoundVert[2][0]->GetCoords(Rx, Ry);
   OutPut(setw(20)<<"T, Wett Len d : " << TDatabase::TimeDB->CURRENTTIME
                  <<"   "<< Rx-Lx<<endl);
   OutPut(setw(20)<<"T, Volume : " << TDatabase::TimeDB->CURRENTTIME
                  <<"   "<< CurrVolume<<endl);
   OutPut(setw(20)<<"T, Volume Diff : "<< TDatabase::TimeDB->CURRENTTIME
                  <<"   "<< CurrVolume - InitVolume << endl);
   
   Getcellangle(FESpaces_All[2], Angle);   
   OutPut( "MinAngle : "<< Angle[0]<< "  MaxAngle : "<<Angle[1]<< endl);
//    exit(0);
   
   TDatabase::ParamDB->P10 = 1; // free surf reparam
   TDatabase::ParamDB->P5 = 0;  // move boundary with velo
    
//======================================================================
// start of time cycle
//======================================================================
  end_time = TDatabase::TimeDB->ENDTIME;
  N_SubSteps = GetN_SubSteps(); 
  TDatabase::TimeDB->CURRENTTIME = TDatabase::TimeDB->STARTTIME;
  m=0;
  limit = TDatabase::ParamDB->SC_NONLIN_RES_NORM_MIN_SADDLE;
  Max_It = TDatabase::ParamDB->SC_NONLIN_MAXIT_SADDLE; 
  solver_time = 0.0;
  N_LinIter = 0; 
  t3 = GetTime();
  total_time = t3 - total_time;
  SetTimeDiscParameters(); 
  
  while(TDatabase::TimeDB->CURRENTTIME< end_time)
  {
    // time cycle
    m++;
    TDatabase::TimeDB->INTERNAL_STARTTIME = TDatabase::TimeDB->CURRENTTIME;

      for(l=0;l<N_SubSteps;l++)   // sub steps of fractional step theta
      {    
        SetTimeDiscParameters();

        if (m==1)
        {
          OutPut("Theta1: " << TDatabase::TimeDB->THETA1<< endl);
          OutPut("Theta2: " << TDatabase::TimeDB->THETA2<< endl);
          OutPut("Theta3: " << TDatabase::TimeDB->THETA3<< endl);
          OutPut("Theta4: " << TDatabase::TimeDB->THETA4<< endl);
        }

        tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
	
	cout << " tau " << tau <<  endl;
	
        TDatabase::TimeDB->CURRENTTIME += tau;

        if (very_first_time)
            oldtau=tau;

        // working rhs array for NSE
        memset(B, 0, N_Unknowns*SizeOfDouble);

        GridVelo_imping(Entries, tmp_Gsol, tmp_Gd, Rhs_All[1],
                        GridKCol, GridRowPtr,
                        GridPos, AuxGridPos,
                        FEVectFuncts_All[0], tau,
                        FEVectFuncts_All[1], MovBoundVert, N_MovVert,
                        Free_Cells, IsoCellEdgeNos, reparam, RefGridPos);
	


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
      
//            cout << " Ddot " <<  Ddot(N_Unknowns, Rhs_All[0], Rhs_All[0]) << endl;
//        exit(0);   
       
       
     SqMat_All[8]->Reset(); // Matrix entries for freesurf int;
     SqMat_All[9]->Reset(); // no need to calculate in nonlinear steps			       

     FreeSurf_axial3D_new(SqMat_All[8], SqMat_All[9],  RHSs[0], RHSs[1], BoundCondition, tau,
                          FEFunctions_All[0]->GetValues(), NULL, Params);
 

     
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
                        Free_Cells, IsoCellEdgeNos, reparam, RefGridPos);
 
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
   // end NSE nonlinear iteration         
   // move the grid
   //======================================================================    
  
   MoveGrid_imping(Entries, tmp_Gsol, tmp_Gd, Rhs_All[1],
                  GridKCol, GridRowPtr,
                  GridPos, FEVectFuncts_All[0], tau,
                  AuxGridPos, 
                  MovBoundVert, N_MovVert,
                  Free_Cells, IsoCellEdgeNos, reparam,
                  N_ReParam);
//   if(fabs(TDatabase::TimeDB->CURRENTTIME - 0.00220711)<1e-8)
// {
//    cout << "test " << endl;
//   exit(0);
// }
  //======================================================================   
  // check freesurf point on the solid surface
  // if so, change bounddary description
  //======================================================================                
  MovBoundVert[2][1]->GetCoords(SLPX, SLPY);
    
  if( SLPY <= 1e-8  )
   {    
    SLPY = 0.;
    MovBoundVert[0][0]->GetCoords(x, y);    
    MovBoundVert[2][1]->SetCoords(SLPX, SLPY);
    UpdateSlipBound->SetParams(x, y, SLPX-x, SLPY);

    Me = Free_Cells[0];
    Joint = Me->GetJoint(IsoCellEdgeNos[1][0]);
      
    IsoJoint = (TIsoBoundEdge *)Joint;
    IsoJoint->DeleteVertices();
    IsoJoint->ChangeEdgeID(BoundaryEdge);      

    Solid_Joint = (TBoundEdge *)Joint;
    Solid_Joint->ChangeBoundComp((Domain->GetBdPart(0))->GetBdComp(0));
    Solid_Joint->UpdateParameters(MovBoundVert[2][0], MovBoundVert[2][1]);
    Joint->ChangeType(BoundaryEdge);
    BoundComp = Solid_Joint->GetBoundComp();
    BoundComp->ChangeType(Line); 
    OutPut("Free surface first vert changed as solid vert"<< " x : "<< SLPX<< "  y : "<<SLPY<< endl); 

    // impose no-penetration condition
    FeId =  FESpaces_All[0]->GetFE2D(IsoCellEdgeNos[0][0], Me);
    FeDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(FeId);
    JointDOF = FeDesc->GetJointDOF(IsoCellEdgeNos[1][0]); 
    N_DOF = FeDesc->GetN_JointDOF();
    DOF = GlobalNumbers+BeginIndex[IsoCellEdgeNos[0][0]];
     
    for(k=0;k<N_DOF;k++)
      Sol_All[0][N_U + DOF[JointDOF[k]] ] = 0.;       
    
    delete [] Bound_Joint[0];
    delete [] Bound_Joint[1];                      

    delete [] MovBoundVert[0];
    delete [] MovBoundVert[1];   
    delete [] MovBoundVert[2];  

    delete [] Free_Joint;    
    delete [] Free_Cells;  
    delete [] IsoCellEdgeNos[0];  
    delete [] IsoCellEdgeNos[1];  
    
    GetMovingBoundData(coll, N_MovVert, Bound_Joint, MovBoundVert, Free_Joint,
                       Free_Cells, IsoCellEdgeNos, SLPX, y);
    
   }  // if( SLPY <= 1e-8  )
       
  //======================================================================          
  // end change boundary description    
  // Updating the Quard points on the solid surfaces
  //======================================================================         
  GridPos->GridToData();
  RefGridPos->GridToData();
  
  // Updating solid boundary 
  MovBoundVert[0][0]->GetCoords(x, y);
  MovBoundVert[2][0]->GetCoords(SLPX, SLPY);  
  UpdateSlipBound->SetParams(x, y, SLPX-x, SLPY-y);
  for(k=0;k<N_MovVert[0];k++)
   if(k==N_MovVert[0]-1)
    { Bound_Joint[0][k]->UpdateParameters(MovBoundVert[0][k], MovBoundVert[2][0]); }
   else
    { Bound_Joint[0][k]->UpdateParameters(MovBoundVert[0][k], MovBoundVert[0][k+1]); }
         
  // Updating axis boundary
   MovBoundVert[1][0]->GetCoords(x, y);
   MovBoundVert[0][0]->GetCoords(SLPX, SLPY);
   UpdateAxialBound->SetParams(x, y, SLPX-x, SLPY-y);
   for(k=0;k<N_MovVert[1];k++)
    if(k==N_MovVert[1]-1)
     { Bound_Joint[1][k]->UpdateParameters(MovBoundVert[1][k], MovBoundVert[0][0]); }
    else
     { Bound_Joint[1][k]->UpdateParameters(MovBoundVert[1][k], MovBoundVert[1][k+1]); }
     
  //======================================================================          
  // end Updating the Quard points on the solid surface
  // Remeshing Begin 
  //======================================================================      
   if((l==0) && ((m % 1) == 0))
    {
     Getcellangle(FESpaces_All[2], Angle);
     OutPut( "MinAngle : "<< Angle[0]<< "  MaxAngle : "<<Angle[1]<< endl);
    }
    
   if((Angle[0]<10.0) ||(Angle[1]>165.0))
    {
      
//      if(TDatabase::ParamDB->P10)
//      {
//       ReParam_axial3D_U(N_MovVert[2], Free_Cells,  IsoCellEdgeNos[1], IsoCellEdgeNos[0],  FEVectFuncts_All[0]);      
//      }

      
      if(TDatabase::ParamDB->WRITE_VTK)
       { 
        os.seekp(std::ios::beg);
        if(img<10) os << "VTK/"<<VtkBaseName<<"_remesh.0000"<<RemeshImg<<".vtk" << ends;
         else if(img<100) os << "VTK/"<<VtkBaseName<<"_remesh.000"<<RemeshImg<<".vtk" << ends;
          else if(img<1000) os << "VTK/"<<VtkBaseName<<"_remesh.00"<<RemeshImg<<".vtk" << ends;
           else if(img<10000) os <<"VTK/"<< VtkBaseName<<"_remesh.0"<<RemeshImg<<".vtk" << ends;
            else  os << "VTK/"<<VtkBaseName<<"_remesh."<<RemeshImg<<".vtk" << ends;
        Output->WriteVtk(os.str().c_str());
        RemeshImg++;
       }  
   
//       t1 = GetTime();
     RemeshAxial3D_ImpDrop(Domain, FESpaces_All, FEVectFuncts_All, FEFunctions_All,
                            N_MovVert, Bound_Joint, MovBoundVert, Free_Joint, Free_Cells, IsoCellEdgeNos,
                            Sol_All, Rhs_All, SquareStructure_All, Structure_All, SqMat_All, Mat_All);
   
   
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

#ifdef __ENERGY__    
     N_thermalDOF = FESpaces_All[3]->GetN_DegreesOfFreedom();
     N_thermalActive = FESpaces_All[3]->GetActiveBound();
     N_thermalNonActive = N_thermalDOF - N_thermalActive;  
     
     delete [] oldsol_T;
     oldsol_T = new double[N_thermalDOF];       
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
       
     
      if(TDatabase::ParamDB->WRITE_VTK)
       { 
        os.seekp(std::ios::beg);
        if(img<10) os << "VTK/"<<VtkBaseName<<"_remesh.0000"<<RemeshImg<<".vtk" << ends;
         else if(img<100) os <<"VTK/"<< VtkBaseName<<"_remesh.000"<<RemeshImg<<".vtk" << ends;
          else if(img<1000) os <<"VTK/"<< VtkBaseName<<"_remesh.00"<<RemeshImg<<".vtk" << ends;
           else if(img<10000) os <<"VTK/"<< VtkBaseName<<"_remesh.0"<<RemeshImg<<".vtk" << ends;
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
  //  Reparametrization of free surface - Begin
  //======================================================================  
//   if(!remeshed && TDatabase::ParamDB->P10!=0)
//    {
//    fhtot = 0;
//    fhmin = 100;
//    fhmax = 0.0;
//    
//    for(k=0;k<N_MovVert[2];k++)
//     {
//      MovBoundVert[2][k]->GetCoords(x1, y1);
// 
//      if(k==N_MovVert[2]-1)
//      { MovBoundVert[1][0]->GetCoords(x2, y2);}
//      else
//      { MovBoundVert[2][k+1]->GetCoords(x2, y2); }
// 
//      fh = sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
//      fhtot +=fh;
//      if (fh < fhmin) fhmin = fh;
//      if (fh > fhmax) fhmax = fh;
//     } // for(k=0;k<N_MovVert[2];k++)
// 
//    fhtot /= (double)N_MovVert[2];
//    fhlimit =  0.8*fhtot; 
//    
//    
//    if ( ((fhmin < fhlimit) || (fhmax > 3.*fhtot/2.)) ) // && (m<10000)
//    {     
//       
//     OutPut("FreeBound Edge Reparam:  "<<img<<' ' << fhlimit <<' '<< 3*fhtot/2.0
//              <<' '<< fhmin <<' '<<fhmax<< endl);    
//          
// //     RefGridPos->GridToData();
// 
//     ReParam_axial3D_U(N_MovVert[2], Free_Cells,  IsoCellEdgeNos[1], IsoCellEdgeNos[0],  FEVectFuncts_All[0]);      
// 
//        
//     N_ReParam++;
    
//     AuxGridPos->GridToData();
//     
// //        if(TDatabase::ParamDB->WRITE_VTK)
// //        { 
// //         os.seekp(std::ios::beg);
// //         if(img<10) os << "VTK/"<<VtkBaseName<<".0000"<<img<<".vtk" << ends;
// //          else if(img<100) os <<"VTK/"<< VtkBaseName<<".000"<<img<<".vtk" << ends;
// //           else if(img<1000) os <<"VTK/"<< VtkBaseName<<".00"<<img<<".vtk" << ends;
// //            else if(img<10000) os << "VTK/"<<VtkBaseName<<".0"<<img<<".vtk" << ends;
// //             else  os << "VTK/"<<VtkBaseName<<"."<<img<<".vtk" << ends;
// //         Output->WriteVtk(os.str().c_str());
// //         img++;
// //        }  
//        
//     memset(Rhs_All[1], 0, 2*N_G*SizeOfDouble);
//     Daxpy(2*N_G, -1, refpos, auxpos);
//     
//     memcpy(Rhs_All[1]+N_GActive, auxpos+N_GActive, N_GBoundaryNodes*SizeOfDouble);
//     memcpy(Rhs_All[1]+(N_G+N_GActive), auxpos+(N_G+N_GActive), N_GBoundaryNodes*SizeOfDouble);
//     
//     memset(Sol_All[1], 0 , 2*N_G*SizeOfDouble);
//     memcpy(Sol_All[1]+N_GActive, auxpos+N_GActive, N_GBoundaryNodes*SizeOfDouble);
//     memcpy(Sol_All[1]+(N_G+N_GActive), auxpos+(N_G+N_GActive), N_GBoundaryNodes*SizeOfDouble);    
//     
//     SolveGridEquation(Entries, Sol_All[1], Rhs_All[1], GridKCol, GridRowPtr, N_G );
// 
//     //no grid velo on Interface, since we interpolated already, so (u-w)=0 
//     if(TDatabase::ParamDB->P10==2)    
//     {      
//      memset(Sol_All[1]+N_GActive, 0, N_GBoundaryNodes*SizeOfDouble);
//      memset(Sol_All[1]+(N_G+N_GActive), 0, N_GBoundaryNodes*SizeOfDouble);
//     }
//     
//     Daxpy(2*N_G, 1, Sol_All[1], ReparamDisp);        
//     Dscal(2*N_G, 1./tau, Sol_All[1]);    
//         
//     memcpy(ReparamMeshVelo, Sol_All[1], 2*N_G*SizeOfDouble);   

//     OutPut("ReParam non-linear CURRENT TIME: ");
//     OutPut(TDatabase::TimeDB->CURRENTTIME << endl);       
//     } //if ( ((fhmin < fhlimit) || (fhm
//    } // if (!remeshed && TDatabase::ParamDB->P10!=0)
  //======================================================================  
  // end  Reparametrization of free surface
  //======================================================================  
  
  
    //======================================================================
    // nonlinear loop due to remeshing
    //first update the convection matrix with the new mesh velocity        
    //======================================================================
   if(remeshed)
    {           
     // working array for rhs is B, initialize B
     memset(B, 0, N_Unknowns*SizeOfDouble);

     if(remeshed)
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

  
     FreeSurf_axial3D_new(SqMat_All[8], SqMat_All[9],  RHSs[0], RHSs[1], BoundCondition, tau,
                          FEFunctions_All[0]->GetValues(), NULL, Params);
     
     
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
    } //   for(j=0;j<Max_It;j++)   
            
    remeshed = FALSE;
         
   } //if (remeshed)  
  
  
  }// for(l=0;l<N_SubSteps;l++) 


   MovBoundVert[0][0]->GetCoords(Lx, Ly);
   MovBoundVert[2][0]->GetCoords(Rx, Ry);
   
   if(MaxWetD<Rx-Lx)
   {
    MaxWetD = Rx-Lx;
    T_MaxWetD = TDatabase::TimeDB->CURRENTTIME;
   }
   
   
   
 if((m % 25 ) == 0   || m==1 )
  {
   MovBoundVert[0][N_MovVert[0]-1]->GetCoords(x1, y1);
   MovBoundVert[2][1]->GetCoords(x2, y2);
   MovBoundVert[2][2]->GetCoords(x3, y3);   
   MovBoundVert[2][3]->GetCoords(x4, y4);
   
   tx = x1-Rx;
   sx = x2-Rx;
   ty = y1-Ry;
   sy = y2-Ry;
   R_Theta[0] = acos( (tx*sx+ty*sy)/(sqrt(tx*tx+ty*ty)* sqrt(sx*sx+sy*sy)) )*(180/3.141592654);

   sx = x3-Rx;
   sy = y3-Ry;
   R_Theta[1] = acos( (tx*sx+ty*sy)/(sqrt(tx*tx+ty*ty)* sqrt(sx*sx+sy*sy)) )*(180./3.141592654);

   sx = ((x4))-Rx;
   sy = ((y4))-Ry;
   R_Theta[2] = acos( (tx*sx+ty*sy)/(sqrt(tx*tx+ty*ty)* sqrt(sx*sx+sy*sy)) )*(180./3.141592654); 
   
   MovBoundVert[1][0]->GetCoords(x1, y1);
   
   if(!remeshed)
    OutPut(setw(25)<<"T, wd,AxialZ,Ucl,RAng 1,2,3: " << TDatabase::TimeDB->CURRENTTIME<<"   "<< Rx-Lx
                   <<"   "<< y1<<"   "<< Params[2]<<"   "<<R_Theta[0]<<"   "<<R_Theta[1]<<"   "<<R_Theta[2]<<endl);    

   Get_KE(FEVectFuncts_All[0], Params);  
   OutPut(setw(25)<<"T, Volume, Diff, Rel. Diff : "<< TDatabase::TimeDB->CURRENTTIME<<"   "<< Params[0]
                  <<"   "<< Params[0] - InitVolume<<"   "<< (Params[0] - InitVolume)/InitVolume << endl);
  }

 if((m % (int)(TDatabase::TimeDB->STEPS_PER_IMAGE) ) == 0   || m==1 )
  {
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


//     if( ((m % 200)== 0) || (m==1) )
//      {
//       os.seekp(std::ios::beg);
//       if(N_BData<10) os << "BDData/Boundary.0000"<<N_BData<<".data" << ends;
//       else if(N_BData<100) os << "BDData/Boundary.000"<<N_BData<<".data" << ends;
//       else if(N_BData<1000) os << "BDData/Boundary.00"<<N_BData<<".data" << ends;
//       else if(N_BData<10000) os << "BDData/Boundary.0"<<N_BData<<".data" << ends;
//       else  os << "BDData/Boundary."<<N_BData<<".data" << ends;
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
//       for(k=0;k<N_MovVert[2];k++) // no need to set end vertices again
//        {
//         MovBoundVert[2][k]->GetCoords(x1, y1);
//         dat << x1 << " " <<  y1<< endl;
//        }
//         MovBoundVert[1][0]->GetCoords(x1, y1);
//         dat << x1 << " " <<  y1<< endl;
//         
//         MovBoundVert[0][0]->GetCoords(x1, y1);
//         dat << x1 << " " <<  y1<< endl;
//         
//       dat.close();
//       cout << endl;
//       cout << "Boundary wrote output into file " << endl;
//       N_BData++;  
//       
//    } //if( ((m % 20)== 0) || (m==1) )
   
   if( m % 200== 0 )
    {
     OutPut(setw(25)<<TDatabase::TimeDB->CURRENTTIME<<" No. ReParam : " << N_ReParam <<endl);  
     OutPut(setw(25)<<TDatabase::TimeDB->CURRENTTIME<<" No. Remeshed : " << N_Remesh <<endl);  
    }
    
//     exit(0);
  } // while(TDatabase::TimeDB->CURRENTTIME< end_time)

  //======================================================================
  // end of time cycle
  //======================================================================   
     if(TDatabase::ParamDB->WRITE_VTK)
       { 
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

