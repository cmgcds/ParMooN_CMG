 
// =======================================================================
//
// Purpose:     main program (no multigrid solver)
//
// Author:      Sashikumaar Ganesan
//
// History:     Implementation started on 24.07.2009
//              OpenMP implementation - 24.07.2009
// =======================================================================
#include <Domain.h>
#include <Database.h>
#include <DiscreteForm2D.h>
#include <FEDatabase2D.h>
#include <FESpace2D.h>
#include <FESpace1D.h>
#include <SquareStructure2D.h>
#include <Structure2D.h>
#include <AuxParam2D.h>
#include <QuadAffin.h>
#include <DirectSolver.h>
#include <Assemble2D.h>
#include <Output2D.h>
#include <LinAlg.h>
#include <CD2DErrorEstimator.h>
#include <FEFunction1D.h>
#include <QuadBilinear.h>
#include <TriaAffin.h>
#include <LineAffin.h>

#include <string.h>
#include <sstream>
#include <MooNMD_Io.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>

#include <MainUtilities.h>
#include <BoundEdge.h>
#include <MacroCell.h>
// #include <BoundEdge.h>
// #include <IsoBoundEdge.h>
// #include <gridgen.h>
//  
// #include <IsoInterfaceJoint.h>
// #include <BdLine.h>
// #include <BdCircle.h>
// #include <GridCell.h>
// =======================================================================
// include current example
// =======================================================================
//#include "../Examples/CD_2D/SineLaplace.h"
//#include "../Examples/CD_2D/TwoBoundaryLayers.h"
//#include "../Examples/CD_2D/TwoBoundaryLayersLaplace.h"
// #include "../Examples/CD_2D/Sphere.h"
// #include "../Examples/CD_2D/Constant1.h"
// #include "../Examples/CD_2D/Cylinder.h"
// #include "../Examples/CD_2D/Plane.h"
// #include "../Examples/CD_2D/MF1.h"
// #include "../Examples/CD_2D/MF2.h"
// #include "../Examples/CD_2D/MF3.h"
// #include "../Examples/CD_2D/BspLutz.h"
//#include "../Examples/CD_2D/Hexagon.h"
// #include "../Examples/CD_2D/Bw_Facing_Step_1_3.h"
// #include "../Examples/CD_2D/LaplaceComputerPraktikum.h"
//#include "../Examples/CD_2D/Robin.h"
// #include "../Examples/CD_2D/Robin_indian_stat.h"
//#include "../Examples/CD_2D/TwoInteriorLayers.h"
//#include "../Examples/CD_2D/PetrSkew.h"
//#include "../Examples/CD_2D/BurmanErn2002.h"
//#include "../Examples/CD_2D/BurmanHansbo2004.h"
//#include "../Examples/CD_2D/HMM1986.h"
//#include "../Examples/CD_2D/JohnMaubachTobiska1997.h"
// #include "../Examples/CD_2D/Hemker1996.h"
//#include "../Examples/CD_2D/ParabolicLayers.h"
//#include "../Examples/CD_2D/MG_Vorlesung_rhs1.h"
//#include "../Examples/CD_2D/pw_linear_rhs.h"
//#include "../Examples/CD_2D/TestReaction.h"
//#include "../Examples/CD_2D/SinSin.h"
//#include "../Examples/CD_2D/Smooth.h"
//#include "../Examples/CD_2D/Poisson_Const.h"
//#include "../Examples/CD_2D/Burman2004.h"
//#include "../Examples/CD_2D/Leveque.h"
#include "../Examples/CD_2D/Mortar.h"

extern "C"
{
  #include "umfpack.h"
}


// insertion sort,
int Sort(int *Array, int length)
{
 int i, j, a, tmp;
 
 for(i=0; i<length-1; i++)
  {
   for(j=i+1; j<length; j++)
    if(Array[i]>Array[j])
     {
      tmp = Array[i];
      Array[i]=Array[j];
      Array[j]=tmp;
     }
  }

}

static int GetIndex(int *Array, int Length, int Element)
{
  int l=0, r=Length, m=(r+l)/2;
  int Mid;

  Mid=Array[m];
  while(Mid!=Element)
  {
    if(Mid<Element)
    { l=m; }
    else
    { r=m; }

    m=(r+l)/2;
    Mid=Array[m];
  }
  return m;
}



void GenMortarInterface(TDomain *Domain, int BdID, TDomain *Domain_MotIntFace,
                        int *&Cell_No, int *&Joint_No, double startX, double startY)
{
 int i, j, k, l, m, n, ID, N_MortarCells, N_Cells, m1, N, N_MortarVert;
 int maxEpV, a, b, len1, len2, Neighb_tmp;
 int *Bd_Part, *Lines, N_G, *PointNeighb, CurrNeib,  Neib[2];
 const int *TmpEV;

 double x, y, x1, y1;

 TCollection *Coll;
 TBaseCell *Me, *Me2;
 TJoint *Joint;
 TVertex **MortarVetrex;
 TBaseCell  **MortarCellTree, *NullCell = NULL;
 boolean StartBD=FALSE, EndBD=FALSE;
 TBoundEdge *BD_Joint;
 TBoundComp *BoundComp;


  Coll = Domain->GetCollection(It_Finest, 0);
  N_MortarCells = 0;

  N_Cells = Coll->GetN_Cells();
  for(i=0;i<N_Cells;i++)
   {
    Me = Coll->GetCell(i);
    k = Me->GetN_Edges();
    for(l=0;l<k;l++)
     {
       Joint = Me->GetJoint(l);
       if(Joint->GetType() == BoundaryEdge)
        {
         BD_Joint = (TBoundEdge *)Joint;
         BoundComp = BD_Joint->GetBoundComp();

         if(BdID==BoundComp->GetID())
          N_MortarCells++;
        }
     } // endfor l
   }// endfor i

  OutPut("N_MortarCells: " << N_MortarCells << endl);

  Cell_No = new int[N_MortarCells];
  Joint_No = new int[N_MortarCells];

  Bd_Part = new int[N_MortarCells];
  Lines = new int[2*N_MortarCells];

  MortarVetrex = new TVertex*[2*N_MortarCells]; // maximum possible vertex

  m = 0;
  N = 0; 
  N_MortarVert = 0;
  for(i=0;i<N_Cells;i++)
   {
    Me = Coll->GetCell(i);
    k = Me->GetN_Edges();
    for(l=0;l<k;l++)
     {
       Joint = Me->GetJoint(l);
       if(Joint->GetType() == BoundaryEdge)
        {
         BD_Joint = (TBoundEdge *)Joint;
         BoundComp = BD_Joint->GetBoundComp();

         if(BdID==BoundComp->GetID())
          {
           Cell_No[N] = i;
           Joint_No[N] = l;
           Bd_Part[N] = BdID;
           N++;
          } // if(BdID==BoundComp->GetID())
        } // if(Joint->GetType()
     } // endfor l
    }// endfor i

  //find start vertex
  for(i=0;i<N_MortarCells;i++)
   {
    Me = Coll->GetCell(Cell_No[i]);
    Me->GetVertex(Joint_No[i])->GetCoords(x, y);

    if(sqrt((x-startX)*(x-startX) + (y-startY)*(y-startY))<1.e-8)
     {
      if(i!=0)
       {
        a = Cell_No[i];
        b = Joint_No[i];
        Cell_No[i] = Cell_No[0];
        Joint_No[i] = Joint_No[0];
        Cell_No[0] = a;
        Joint_No[0] = b;
       }
      break;
     }
   }

  if(i==N_MortarCells)
   {
    OutPut(i << " N_MortarCells start vert not found: " << N_MortarCells << endl);
    exit(0);
   }

  for(i=0;i<N_MortarCells-1;i++)
   {
    Me = Coll->GetCell(Cell_No[i]);
    k = Me->GetN_Edges();
    Me->GetVertex((Joint_No[i]+1) % k)->GetCoords(x, y);

    for(l=i+1;l<N_MortarCells;l++)
     {
      Me2 = Coll->GetCell(Cell_No[l]);
      Me2->GetVertex(Joint_No[l])->GetCoords(x1, y1);

      if(sqrt((x-x1)*(x-x1) + (y-y1)*(y-y1))<1.e-8)
       {
        if(l==i+1)
         { break; }
        else
         {
          a = Cell_No[l];
          b = Joint_No[l];
          Cell_No[l] = Cell_No[i+1];
          Joint_No[l] = Joint_No[i+1];
          Cell_No[i+1] = a;
          Joint_No[i+1] = b;
          break; 
         }
       }
     }// for(l=0;l<N_MortarCells;l++)
   }

  for(i=0;i<N_MortarCells;i++)
   {
    Me = Coll->GetCell(Cell_No[i]);
    MortarVetrex[i] = Me->GetVertex(Joint_No[i]);
    Lines[2*i] = i;
    Lines[2*i+1] = i+1;
   }

  //end vertex
  Me = Coll->GetCell(Cell_No[N_MortarCells-1]);
  k = Me->GetN_Edges();
  MortarVetrex[N_MortarCells] =  Me->GetVertex((Joint_No[N_MortarCells-1]+1) % k);
  N_MortarVert=N_MortarCells+1;

//   for(i=0;i<N_MortarCells+1;i++)
//    {
//     MortarVetrex[i]->GetCoords(x, y);
//     OutPut(i << " N_MortarCell : " << x << endl);
//    }
//     OutPut(i << " N_MortarCells start vert not found: " << N_MortarCells << endl);
//     exit(0);

  MortarCellTree = new TBaseCell*[N_MortarCells];

  for (i=0;i<N_MortarCells;i++)
   {
    MortarCellTree[i] = new TMacroCell(TDatabase::RefDescDB[S_Line], 0);
    MortarCellTree[i]->SetVertex(0, MortarVetrex[ Lines[ 2*i       ]]);
    MortarCellTree[i]->SetVertex(1, MortarVetrex[ Lines[ 2*i + 1]]);

    ((TMacroCell *) MortarCellTree[i])->SetSubGridID(0);
    ((TBaseCell *) MortarCellTree[i])->SetBd_Part(Bd_Part[i] );
   }

  OutPut("N_MortarVert: " << N_MortarVert << endl);

   Domain_MotIntFace->SetTreeInfo(MortarCellTree, N_MortarCells);

   TDatabase::IteratorDB[It_EQ]->SetParam(Domain_MotIntFace);
   TDatabase::IteratorDB[It_LE]->SetParam(Domain_MotIntFace);
   TDatabase::IteratorDB[It_Finest]->SetParam(Domain_MotIntFace);
   TDatabase::IteratorDB[It_Between]->SetParam(Domain_MotIntFace);
   TDatabase::IteratorDB[It_OCAF]->SetParam(Domain_MotIntFace);

  // search neighbours
   maxEpV = 2; // max possible degree of each vertex
   N_G = N_MortarVert;
   maxEpV++;
   PointNeighb = new int[maxEpV * N_G];
   memset(PointNeighb, 0, maxEpV * N_G *SizeOfInt);

    // every row contains "maxEpV" columns
    // for every row at first colomn contains the number of cells containing this vertex
    // at further columns we set the index of the corresponding cells containing this vertex
   for(i=0;i<2*N_MortarCells;i++)
   {
     j = Lines[i]*maxEpV;
     PointNeighb[j]++;
     PointNeighb[j + PointNeighb[j]] = i / 2;
   }

  N_G = N_MortarCells;


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
      Joint = new TJointEqN(MortarCellTree[Neib[0]]);
      MortarCellTree[Neib[0]]->SetJoint(1, Joint);
      cout<<"End cell " << Neib[0] <<" Neib[1] " << Neib[1] <<endl;
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
      Joint = new TJointEqN(MortarCellTree[Neib[0]], MortarCellTree[Neib[1]]);

      MortarCellTree[Neib[0]]->SetJoint(0, Joint);
      MortarCellTree[Neib[1]]->SetJoint(1, Joint);
     }
    else if(len1==1)
     {
      Joint = new TJointEqN(MortarCellTree[Neib[0]]);
      MortarCellTree[Neib[0]]->SetJoint(0, Joint);
     cout<< "Start cell " << Neib[0] <<" Neib[1] " << Neib[1] <<endl;
     }

   } //  for(i=0;i<N_G;i++)

  delete [] PointNeighb;
  delete []  Bd_Part;
  delete []  Lines;

}

void Map2DNonMot_1DMotIntFace(TCollection *Coll_2D, int BdID, TCollection *Coll_1D,
                              int *&Cell_No, int *&Joint_No, double startX, double startY,
                              int *&N_NonMotNeibs, int *&NonMotNeibs, int &N_NonMortarEdges)
{
 int i, j, k, l, N_Cells_1D, N_Cells_2D;
 int a, b;

 double x, y, x0, x1, y0, y1, m0, m1;

 TBaseCell *Me, *Me2;
 TJoint *Joint;
 TBoundEdge *BD_Joint;
 TBoundComp *BoundComp;
 TVertex *Vetrex;

 bool PtInCell;

  N_Cells_2D = Coll_2D->GetN_Cells();

  N_NonMortarEdges = 0;
  for(i=0;i<N_Cells_2D;i++)
   {
    Me = Coll_2D->GetCell(i);
    k = Me->GetN_Edges();
    for(l=0;l<k;l++)
     {
       Joint = Me->GetJoint(l);
       if(Joint->GetType() == BoundaryEdge)
        {
         BD_Joint = (TBoundEdge *)Joint;
         BoundComp = BD_Joint->GetBoundComp();

         if(BdID==BoundComp->GetID())
          N_NonMortarEdges++;
        }
     } // endfor l
   }// endfor i


  Cell_No = new int[N_NonMortarEdges];
  Joint_No = new int[N_NonMortarEdges];

  N_NonMortarEdges = 0;
  for(i=0;i<N_Cells_2D;i++)
   {
    Me = Coll_2D->GetCell(i);
    k = Me->GetN_Edges();
    for(l=0;l<k;l++)
     {
       Joint = Me->GetJoint(l);
       if(Joint->GetType() == BoundaryEdge)
        {
         BD_Joint = (TBoundEdge *)Joint;
         BoundComp = BD_Joint->GetBoundComp();

         if(BdID==BoundComp->GetID())
          {
           Cell_No[N_NonMortarEdges] = i;
           Joint_No[N_NonMortarEdges] = l;
           N_NonMortarEdges++;
          }
        }
     } // endfor l
   }// endfor i


  //find the end vertex
  for(i=0;i<N_NonMortarEdges;i++)
   {
    Me = Coll_2D->GetCell(Cell_No[i]);
    k = Me->GetN_Edges();
    Me->GetVertex((Joint_No[i]+1) %k)->GetCoords(x, y);

    if(sqrt((x-startX)*(x-startX) + (y-startY)*(y-startY))<1.e-8)
     {
      if(i!=N_NonMortarEdges-1)
       {
        a = Cell_No[i];
        b = Joint_No[i];
        Cell_No[i] = Cell_No[N_NonMortarEdges-1];
        Joint_No[i] = Joint_No[N_NonMortarEdges-1];
        Cell_No[N_NonMortarEdges-1] = a;
        Joint_No[N_NonMortarEdges-1] = b;
       }
       break;

     }
   }

  if(i==N_NonMortarEdges)
   {
    OutPut(i << " N_NonMortarEdges start vert not found: " << N_NonMortarEdges << endl);
    exit(0);
   }

  for(i=N_NonMortarEdges-1;i>0;i--)
   {
    Me = Coll_2D->GetCell(Cell_No[i]);
    Me->GetVertex(Joint_No[i])->GetCoords(x, y);

    for(l=i-1;l>=0;l--)
     {
      Me2 = Coll_2D->GetCell(Cell_No[l]);
      k = Me->GetN_Edges();
      Me2->GetVertex((Joint_No[l]+1)%k)->GetCoords(x1, y1);

      if(sqrt((x-x1)*(x-x1) + (y-y1)*(y-y1))<1.e-8)
       {
        if(l==i-1)
         { break; }
        else
         {
          a = Cell_No[l];
          b = Joint_No[l];
          Cell_No[l] = Cell_No[i-1];
          Joint_No[l] = Joint_No[i-1];
          Cell_No[i-1] = a;
          Joint_No[i-1] = b;
          break; 
         }
       }
     }// for(l=i+1;l>=0;l++)
   }

 // for each 1D cell find neibs nonmortar edges
 N_Cells_1D = Coll_1D->GetN_Cells();
 N_NonMotNeibs = new int[N_Cells_1D];
 NonMotNeibs = new int[10*N_Cells_1D];

  for(i=0;i<N_Cells_1D;i++)
    N_NonMotNeibs[i]=0;

  for(i=0;i<N_Cells_1D;i++)
   {
    Me = Coll_1D->GetCell(i);

    Vetrex = Me->GetVertex(0);
    Vetrex->GetCoords(x0, y);
    Vetrex = Me->GetVertex(1);
    Vetrex->GetCoords(x1, y);

    // find all nonmortar cells which contain this point
    // nonmot cells are ordered in opposite direction
    for(j=N_NonMortarEdges-1;j>=0;j--)
     {
      Me2 = Coll_2D->GetCell(Cell_No[j]);
      k = Me2->GetN_Edges();

      Me2->GetVertex(Joint_No[j])->GetCoords(m1, y);
      Me2->GetVertex((Joint_No[j]+1) %k)->GetCoords(m0, y);

      if((m0>=x0 && m1<=x1) || (m0<=x0 && m1>=x1) ||
          (m0>=x0 && m0<=x1) || (m1>=x0 && m1<=x1) )
        {
         k = N_NonMotNeibs[i];

         if(k==5)
          {
           cout<< "Increase the array size, only 5 NonMotEdges are allowesd now !! " <<endl;
           //change in TStructure2D also
           exit(0);
          }

         NonMotNeibs[10*i + k ] = Cell_No[j];
         NonMotNeibs[10*i + 5 + k ] = Joint_No[j];
         N_NonMotNeibs[i]++;
        // OutPut(i << " " << j << " N_MortarCell find : " << x0 <<  " " <<x1<< endl)
        }
     }
   }
} //Map2DNonMot_1DMotIntFace

void GetN_NonMortaDof(TFESpace2D *FESpace_NonMot, int N_NonMortarEdges, int **NonMortarCellEdge, 
                      int &N_NonMortaDofs, int *&GlobalDof_NonMortaEdgeDof)
{
  int i, j, k,l, IJoint, I, dof, N_DOF_NonMot;
  int *Cells_No, *Joint_No, *DOF, *BeginIndex2D, *GlobalNumbers2D;
  int *JointDOF, N_JointDOF, *AddedDof;

  bool Added;

  TBaseCell *Me;
  TCollection *Coll;
  FE2D FEId;
  TFEDesc2D *FeDesc;

  Cells_No = NonMortarCellEdge[0];
  Joint_No = NonMortarCellEdge[1];

  Coll = FESpace_NonMot->GetCollection();
  GlobalNumbers2D=FESpace_NonMot->GetGlobalNumbers();
  BeginIndex2D=FESpace_NonMot->GetBeginIndex();

  N_DOF_NonMot = FESpace_NonMot->GetN_DegreesOfFreedom();
  AddedDof = new int[N_DOF_NonMot];

  N_NonMortaDofs = 0;
  for(i=0; i<N_NonMortarEdges; i++)
   {
    I = Cells_No[i];
    IJoint = Joint_No[i];
    Me = Coll->GetCell(I);
    FEId = FESpace_NonMot->GetFE2D(I, Me);
    FeDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(FEId);
    JointDOF = FeDesc->GetJointDOF(IJoint);
    N_JointDOF = FeDesc->GetN_JointDOF();
    DOF = GlobalNumbers2D + BeginIndex2D[I];

    for(j=0; j<N_JointDOF; j++)
     {
      dof = DOF[JointDOF[j]];
      // cout << i << "   " << j << " Cell2D_No  " <<I <<" IJoint  " <<IJoint << " dof  " <<dof << endl;
      Added = FALSE;
      for(k=0; k<N_NonMortaDofs; k++)
       if(AddedDof[k] == dof)
        {
         Added = TRUE;
         break;
        }

      if(!Added)
       {
        AddedDof[N_NonMortaDofs] = dof;
        N_NonMortaDofs++;
       }

     }

   }

  GlobalDof_NonMortaEdgeDof = new int[N_NonMortaDofs];

  for(i=0; i<N_NonMortaDofs; i++)
   GlobalDof_NonMortaEdgeDof[i] = AddedDof[i];

  Sort(GlobalDof_NonMortaEdgeDof, N_NonMortaDofs);

//   for(i=0; i<N_NonMortaDofs; i++)
//   cout << i << " GlobalDof_NonMortaEdgeDof " << GlobalDof_NonMortaEdgeDof[i] <<endl;

 delete [] AddedDof;

cout << "N_NonMortarEdges " << N_NonMortarEdges <<" N_NonMortaDofs " << N_NonMortaDofs << endl;
// exit(0);

} // GetN_NonMortaDof

void GetN_NonMortaDof(TFESpace1D *FeSpace_MortIntFace, TFESpace2D *FESpace_NonMot, 
        int N_NonMortarEdges, int **NonMortarCellEdge, int *N_NonMotNeibs, int *NonMotNeibs,
        TNonMortarData *NonMortarFEData)
{
  int i, j, k,l, N_Cells, IJoint, I, dof, N_DOF_NonMot;
  int *Cells_No, *Joint_No, *DOF, *BeginIndex2D, *GlobalNumbers2D;
  int *JointDOF, N_JointDOF, *AddedDof, N_NonMortaDofs, *GlobalDof;
  int N_NonMotEdges, N_ColAnsatzs, *N_DofPerEdge, total, LocIndex;
  int *NonMotGlobalNumber, *NonMotLocGlobalNumber, *NonMotBeginIndex;
  int *AuxLocDof, N_LocDof;

  bool Added, UPDATE;

  TBaseCell *Me;
  TCollection *Coll;
  FE2D FEId;
  TFEDesc2D *FeDesc;

  Cells_No = NonMortarCellEdge[0];
  Joint_No = NonMortarCellEdge[1];

  Coll = FESpace_NonMot->GetCollection();
  GlobalNumbers2D=FESpace_NonMot->GetGlobalNumbers();
  BeginIndex2D=FESpace_NonMot->GetBeginIndex();

  N_DOF_NonMot = FESpace_NonMot->GetN_DegreesOfFreedom();
  AddedDof = new int[N_DOF_NonMot];

  N_NonMortaDofs = 0;
  for(i=0; i<N_NonMortarEdges; i++)
   {
    I = Cells_No[i];
    IJoint = Joint_No[i];
    Me = Coll->GetCell(I);
    FEId = FESpace_NonMot->GetFE2D(I, Me);
    FeDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(FEId);
    JointDOF = FeDesc->GetJointDOF(IJoint);
    N_JointDOF = FeDesc->GetN_JointDOF();
    DOF = GlobalNumbers2D + BeginIndex2D[I];

    for(j=0; j<N_JointDOF; j++)
     {
      dof = DOF[JointDOF[j]];
      // cout << i << "   " << j << " Cell2D_No  " <<I <<" IJoint  " <<IJoint << " dof  " <<dof << endl;
      Added = FALSE;
      for(k=0; k<N_NonMortaDofs; k++)
       if(AddedDof[k] == dof)
        {
         Added = TRUE;
         break;
        }

      if(!Added)
       {
        AddedDof[N_NonMortaDofs] = dof;
        N_NonMortaDofs++;
       }
     }
   }

  // fill the structure 
  NonMortarFEData->N_NonMortaDofs = N_NonMortaDofs;

  NonMortarFEData->GlobalDof_NonMortaEdgeDof =  new int[N_NonMortaDofs];
  GlobalDof =  NonMortarFEData->GlobalDof_NonMortaEdgeDof;

  for(i=0; i<N_NonMortaDofs; i++)
   GlobalDof[i] = AddedDof[i];

  Sort(GlobalDof, N_NonMortaDofs);

//   for(i=0; i<N_NonMortaDofs; i++)
//   cout << i << " GlobalDof " << GlobalDof[i] <<endl;

 delete [] AddedDof;

  N_Cells = FeSpace_MortIntFace->GetN_Cells();
  NonMortarFEData->N_DofPerEdge = new int[N_Cells];
  NonMortarFEData->EdgeNonMotBeginIndex = new int[N_Cells];;

  N_DofPerEdge = NonMortarFEData->N_DofPerEdge;
  NonMotBeginIndex = NonMortarFEData->EdgeNonMotBeginIndex;


  AuxLocDof = new int[N_NonMortaDofs];

  total = 0;
  NonMotBeginIndex[0] = 0;
  for(i=0;i<N_Cells;i++)
   {
    N_NonMotEdges = N_NonMotNeibs[i];

    N_LocDof = 0;
    for(j=0;j<N_NonMotEdges;j++)
     {
      I = NonMotNeibs[10*i + j];
      IJoint = NonMotNeibs[10*i + 5 + j];

      Me = Coll->GetCell(I);
      FEId = FESpace_NonMot->GetFE2D(I, Me);
      FeDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(FEId);

      JointDOF = FeDesc->GetJointDOF(IJoint);
      N_JointDOF = FeDesc->GetN_JointDOF();
      DOF = GlobalNumbers2D + BeginIndex2D[I];

      for(k=0;k<N_JointDOF;k++)
       {
        dof =  DOF[JointDOF[k]];

        UPDATE = TRUE;
        for(l=0;l<N_LocDof;l++)
         if(AuxLocDof[l] == dof)
          {
           UPDATE=FALSE;
           break;
          }

       if(UPDATE)
         {
          AuxLocDof[N_LocDof] = dof;
          N_LocDof++;

         }
       }
      }

     N_DofPerEdge[i] = N_LocDof;
     total += N_LocDof;
     NonMotBeginIndex[i+1] = total;
    }


  NonMortarFEData->EdgeNonMotGlobalNo = new int[total];
  NonMortarFEData->EdgeNonMotLocGlobalNo = new int[total];
  NonMotGlobalNumber = NonMortarFEData->EdgeNonMotGlobalNo;
  NonMotLocGlobalNumber = NonMortarFEData->EdgeNonMotLocGlobalNo;


  for(i=0;i<total;i++)
   NonMotLocGlobalNumber[i] = -1;

  total = 0;
  for(i=0;i<N_Cells;i++)
   {
    N_NonMotEdges = N_NonMotNeibs[i];
    N_LocDof = 0;

    for(j=0;j<N_NonMotEdges;j++)
     {
      I = NonMotNeibs[10*i + j];
      IJoint = NonMotNeibs[10*i + 5 + j];

      Me = Coll->GetCell(I);
      FEId = FESpace_NonMot->GetFE2D(I, Me);
      FeDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(FEId);

      JointDOF = FeDesc->GetJointDOF(IJoint);
      N_JointDOF = FeDesc->GetN_JointDOF();
      DOF = GlobalNumbers2D + BeginIndex2D[I];

      for(k=0;k<N_JointDOF;k++)
       {
        dof =  DOF[JointDOF[k]];

        UPDATE = TRUE;
        for(l=0;l<N_LocDof;l++)
         if(AuxLocDof[l] == dof)
          {
           UPDATE=FALSE;
           break;
          }

        if(UPDATE)
         {
          AuxLocDof[N_LocDof] = dof;
          N_LocDof++;

          LocIndex = GetIndex(GlobalDof, N_NonMortaDofs, dof);

          NonMotGlobalNumber[total] = dof;
          NonMotLocGlobalNumber[total] = LocIndex;
          total++;
         }
       }
     }
    }

 delete [] AuxLocDof;


//   //check
//   for(i=0;i<N_Cells;i++)
//    {
//     I = NonMotBeginIndex[i];
//     for(j=0;j<N_DofPerEdge[i];j++)
//      cout << "GlobalNumber " << NonMotGlobalNumber[I+j] <<" LocGlobalNumber " 
//           << NonMotLocGlobalNumber[I+j] << endl;
//     cout<<endl;
//    }

//     cout <<  "GetN_NonMortaDof   " << (*NonMortarFEData).N_NonMortaDofs << endl;

// cout << "N_NonMortarEdges " << N_NonMortarEdges <<" N_NonMortaDofs " << N_NonMortaDofs << endl;
// exit(0);

} //GetN_NonMortaDof


void AssembleB_Mat(TMatrix2D *B_Mot, TFEFunction1D *fefunct_low, TFESpace2D *fespace, 
                   int **MortarCellEdge, int *MotLocalGlobalDof, int *GlobalDof_MortaEdgeDof)
{
 int i, j, k, l, M, N, *N_BaseFuncts, N_Cells_low;
 int *Cell_array, *Joint_array, N_BaseFunct, N_BaseFunct_low, N_LinePoints;
 int *GlobalNumbers, *BeginIndex, *GlobalNumbers_low, *BeginIndex_low, *DOF, *DOF_low;
 int *JointDOF, N_JointDOF, N_Sets=1, locdof;
 int *KCol, *RowPtr, begin, end;
 
 double *LineWeights, *zeta, Mult, val, test0, test1, ansatz0, ansatz1;
 double **uref, *uorig, *ValuesB;
 double psi[MaxN_QuadPoints_1D], AbsDetjk[MaxN_QuadPoints_1D];
 double **origvaluesD0, *orgD0;
 double LocMatrixB[MaxN_BaseFunctions1D*MaxN_BaseFunctions1D];
 
 bool Needs2ndDer[1];
 
  TBaseCell *cell, *Cell_low;
  TCollection *Coll, *Coll_low;
  TFESpace1D *fespace_low;
  TJoint *joint;
  FE2D FEId;
  FE1D FEId_low;
  BF2DRefElements RefElement;
  RefTrans2D RefTrans;
  TRefTrans2D *F_K;
  QuadFormula1D LineQuadFormula;
  TQuadFormula1D *qf1;
  BaseFunct2D *BaseFuncts;
  TFEDesc2D *FeDesc;
  TRefTrans1D *F_K_low;
    
  TFE1D *Element_low;
  TBaseFunct1D *bf_low;
  BaseFunct1D BaseFunct_ID_low, BaseFunct_low[1];
 
  BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
  N_BaseFuncts = TFEDatabase2D::GetN_BaseFunctFromFE2D();

  GlobalNumbers = fespace->GetGlobalNumbers();
  BeginIndex = fespace->GetBeginIndex();  
  Coll = fespace->GetCollection();

  fespace_low = fefunct_low->GetFESpace1D();
  GlobalNumbers_low = fespace_low->GetGlobalNumbers();
  BeginIndex_low = fespace_low->GetBeginIndex();   
  Coll_low = fespace_low->GetCollection();
  N_Cells_low = Coll_low->GetN_Cells();
  
  Cell_array = MortarCellEdge[0];
  Joint_array = MortarCellEdge[1];
  
  Needs2ndDer[0] = FALSE;

  RowPtr = B_Mot->GetRowPtr();
  KCol = B_Mot->GetKCol();
  ValuesB = B_Mot->GetEntries();

  
  for(i=0;i<N_Cells_low;i++)
   {
    Cell_low = Coll_low->GetCell(i);
    FEId_low =  fespace_low->GetFE1D(i, Cell_low);
    Element_low = TFEDatabase2D::GetFE1D(FEId_low);
    bf_low = Element_low->GetBaseFunct1D();
    N_BaseFunct_low = Element_low->GetN_DOF();
    BaseFunct_ID_low = Element_low->GetBaseFunct1D_ID();     
     
     
    M = Cell_array[i];
    N = Joint_array[i];
    // ansatz
    cell = Coll->GetCell(M);
    FEId = fespace->GetFE2D(M, cell);
    RefElement = TFEDatabase2D::GetRefElementFromFE2D(FEId);
    FeDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(FEId);

    // set reference element
    switch(RefElement)
     {
      case BFUnitSquare:
          RefTrans = QuadBilinear;
          F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
          ((TQuadBilinear *)F_K)->SetCell(cell);
      break;

      case BFUnitTriangle:
          RefTrans = TriaAffin;
          F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
          ((TTriaAffin *)F_K)->SetCell(cell);
          //cout << "Cell RefElement" << RefElement  << endl;
      break;
    } // endswitch    

    // select quadrature formulae   
    l = bf_low->GetPolynomialDegree();
    if(l<TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId))
      l = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
    
    LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(2*l);
    qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
    qf1->GetFormulaData(N_LinePoints, LineWeights, zeta); 
    TFEDatabase2D::GetBaseFunct2DFromFE2D(FEId)->MakeRefElementData(LineQuadFormula);     
     
    // get the 2D cell dof
    DOF = GlobalNumbers + BeginIndex[M];
    N_BaseFunct = N_BaseFuncts[FEId];
    // get the isoedge joint dof
    JointDOF = FeDesc->GetJointDOF(N);
    N_JointDOF = FeDesc->GetN_JointDOF();

    if(N_BaseFunct_low!=N_JointDOF)
     {
       cout << " N_BaseFunct_low!=N_JointDOF " << endl;
       exit(0);
     }


    uref = TFEDatabase2D::GetJointValues2D(BaseFuncts[FEId], LineQuadFormula, N);
//     uxiref = TFEDatabase2D::GetJointDerivatives2D(BaseFuncts[FEId], LineQuadFormula, N, D10);
//     uetaref = TFEDatabase2D::GetJointDerivatives2D(BaseFuncts[FEId], LineQuadFormula, N, D01);

    // get 1D cell dof
    DOF_low = GlobalNumbers_low +  BeginIndex_low[i]; 


    for(j=0;j<N_JointDOF;j++)
     {
      MotLocalGlobalDof[ DOF[JointDOF[j]] ] = DOF_low[j];
      GlobalDof_MortaEdgeDof[DOF_low[j]] = DOF[JointDOF[j]];
     }

    F_K_low = TFEDatabase2D::GetRefTrans1D(LineAffin);
    ((TLineAffin *)F_K_low)->SetCell(Cell_low);
    
    ((TLineAffin *)F_K_low)->GetOrigFromRef(N_LinePoints, zeta, psi, AbsDetjk);
//     for(j=0;j<N_LinePoints;j++)
//      cout << " psi " << psi[j]<< endl;
    
    
    BaseFunct_low[0] = BaseFunct_ID_low;
    ((TLineAffin *)F_K_low)->GetOrigValues(N_Sets, BaseFunct_low, N_LinePoints, zeta,  LineQuadFormula,  Needs2ndDer);

    origvaluesD0=TFEDatabase2D::GetOrigElementValues(BaseFunct_ID_low, D0);

    
    memset(LocMatrixB, 0, N_BaseFunct_low*N_BaseFunct_low*SizeOfDouble);

    for(j=0;j<N_LinePoints;j++)
     {    
      //cout<< " c  " << c0 << " g " <<g0 << " c1 " << c1 <<endl;   
      Mult = LineWeights[j]*AbsDetjk[j];
      orgD0 = origvaluesD0[j];
      uorig = uref[j];  
      
      for(k=0;k<N_BaseFunct_low;k++) // test
       {
        test0  = orgD0[k];

        for(l=0;l<N_JointDOF;l++) // ansatz
         {
          locdof = JointDOF[l];
          ansatz0  = uorig[locdof];

          val  = ansatz0*test0;
          //cout<< "A ( " << k <<" , " <<l << " ) " << Mult*val << endl;
          LocMatrixB[k*N_BaseFunct_low + l] += (Mult*val);
        } 
       }         
     } //N_LinePoints  
     
     // add to global matrices
    for(j=0;j<N_BaseFunct_low;j++)
     {
      locdof = DOF_low[j];

      begin = RowPtr[locdof];
      end = RowPtr[locdof+1];
      for(k=begin;k<end;k++)
       {
       for(l=0;l<N_BaseFunct_low;l++)
        {
         if(KCol[k] == DOF_low[l])
          {
           ValuesB[k] +=LocMatrixB[j*N_BaseFunct_low + l];  
           break;
          }
        } // for(l=0;l<N_BaseFunct_low;l++)
      } // for(k=begin;k<end;k++)
     } // for(j=0;j<N_BaseFunct_low;j++)
     
   } //  for(i=0;i<N_Cells_low;i++)
  
  
//   //print matrix
//   N = fespace_low->GetN_DegreesOfFreedom(); 
//   for(j=0;j<N;j++)
//    {
//     begin = RowPtr[j];
//     end = RowPtr[j+1];
//     for(k=begin;k<end;k++)
//      {
//       cout << "B(" << j << ", "<< KCol[k] << ") = " << ValuesB[k] <<endl;
//      }
//     cout<<endl;
//    }



}  // AssembleB_Mat


void AssembleB_NonMat(TMatrix2D *B_NonMot, TMatrix2D *B_NonMotT, TFEFunction1D *fefunct_low, TFESpace2D *fespace,
		      int N_NonMortarEdges, int **NonMortarCellEdge, int *N_NonMotNeibs, int *NonMotNeibs, TNonMortarData *NonMortarFEData, int *AssembleB_NonMat)
{
 int i, j, k, l, m, M, N, *N_BaseFuncts, N_Cells_low;
 int *Cell_array, *Joint_array, N_BaseFunct, N_BaseFunct_low, N_LinePoints;
 int *GlobalNumbers, *BeginIndex, *GlobalNumbers_low, *BeginIndex_low, *DOF_1D, *DOF_low;
 int *JointDOF, N_JointDOF, N_Sets=1, locdof;
 int *KCol, *RowPtr, *KColT, *RowPtrT, begin, end, kk, ll;
 int N_NonMotNeibsInCell;
 int *N_NeededQuadPtValues, *NonMotCells, *NonMotEdges, IJoint;
 int **NonMotCellLocalIndex, disp, *N_DofPerEdge;
 int DOF_NonMot[MaxN_BaseFunctions1D*10], dof, N_NonMortaDofs, *GlobalDof, *DOF;

 double *LineWeights, *zeta, Mult, val, test0, test1, ansatz0, ansatz1;
 double **uref, *uorig, *ValuesB;
 double psi[MaxN_QuadPoints_1D], AbsDetjk[MaxN_QuadPoints_1D];
 double **origvaluesD0, *orgD0;
 double LocMatrixB[MaxN_BaseFunctions1D*MaxN_BaseFunctions1D*10];
 double LocMatrixBT[MaxN_BaseFunctions1D*MaxN_BaseFunctions1D*10];
 double **QuadPts, startX, endX, x, y, *ValuesBT;
 double DummyWeights[MaxN_QuadPoints_1D], DummyZeta[MaxN_QuadPoints_1D];
 double **NonMotBasisVal;

 bool Needs2ndDer[1], UPDATED;
 
  TBaseCell *cell, *Cell_low;
  TCollection *Coll, *Coll_low;
  TFESpace1D *fespace_low;
  TJoint *joint;
  FE2D FEId;
  FE1D FEId_low;
  BF2DRefElements RefElement;
  RefTrans2D RefTrans;
  TRefTrans2D *F_K;
  QuadFormula1D LineQuadFormula;
  TQuadFormula1D *qf1, *Dummyqf1;
  BaseFunct2D *BaseFuncts;
  TFEDesc2D *FeDesc;
  TRefTrans1D *F_K_low;

  TFE1D *Element_low;
  TBaseFunct1D *bf_low;
  BaseFunct1D BaseFunct_ID_low, BaseFunct_low[1]; 

  BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
  N_BaseFuncts = TFEDatabase2D::GetN_BaseFunctFromFE2D();

  Coll = fespace->GetCollection();
  GlobalNumbers = fespace->GetGlobalNumbers();
  BeginIndex = fespace->GetBeginIndex();

  N_DofPerEdge = NonMortarFEData->N_DofPerEdge;
  N_NonMortaDofs = NonMortarFEData->N_NonMortaDofs;
  GlobalDof =  NonMortarFEData->GlobalDof_NonMortaEdgeDof;

  fespace_low = fefunct_low->GetFESpace1D();
  GlobalNumbers_low = fespace_low->GetGlobalNumbers();
  BeginIndex_low = fespace_low->GetBeginIndex();   
  Coll_low = fespace_low->GetCollection();
  N_Cells_low = Coll_low->GetN_Cells();

  NonMotCells = NonMortarCellEdge[0];
  NonMotEdges = NonMortarCellEdge[1];

  N_NeededQuadPtValues = new int[N_NonMortarEdges];
  QuadPts = new double*[N_NonMortarEdges];
  NonMotBasisVal = new double*[N_NonMortarEdges]; 

  NonMotCellLocalIndex = new int*[N_Cells_low];
  for(i=0;i<N_Cells_low;i++)
    NonMotCellLocalIndex[i] = new int[N_NonMotNeibs[i]];

  memset(N_NeededQuadPtValues, 0, N_NonMortarEdges*SizeOfInt);  

  Needs2ndDer[0] = FALSE;

  RowPtr = B_NonMot->GetRowPtr();
  KCol = B_NonMot->GetKCol();
  ValuesB = B_NonMot->GetEntries();

  RowPtrT = B_NonMotT->GetRowPtr();
  KColT = B_NonMotT->GetKCol();
  ValuesBT = B_NonMotT->GetEntries();

  for(i=0;i<N_Cells_low;i++)
   {
    Cell_low = Coll_low->GetCell(i);
    FEId_low =  fespace_low->GetFE1D(i, Cell_low);
    Element_low = TFEDatabase2D::GetFE1D(FEId_low);
    bf_low = Element_low->GetBaseFunct1D();

    l = bf_low->GetPolynomialDegree();
    LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(2*l);
    qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
    qf1->GetFormulaData(N_LinePoints, LineWeights, zeta);

    for(j=0;j<N_NonMotNeibs[i];j++)
     {
      N = NonMotNeibs[10*i + j];
      //find local index
      k = 0;
      while(NonMotCells[k] != N ) k++;

      NonMotCellLocalIndex[i][j] = k;
      N_NeededQuadPtValues[k] +=N_LinePoints;
     }
   }

//        for(i=0;i<N_NonMortarEdges;i++)
//     cout << "N_Needed  " << N_NeededQuadPtValues[i] << endl;


  for(i=0;i<N_NonMortarEdges;i++)
   {
    QuadPts[i]  = new double[N_NeededQuadPtValues[i]];
    NonMotBasisVal[i]  = new double[N_NeededQuadPtValues[i]*MaxN_BaseFunctions2D];
    memset(NonMotBasisVal[i], 0, N_NeededQuadPtValues[i]*MaxN_BaseFunctions2D*SizeOfInt);
   }

  memset(N_NeededQuadPtValues, 0, N_NonMortarEdges*SizeOfInt);

  for(i=0;i<N_Cells_low;i++)
   {
    Cell_low = Coll_low->GetCell(i);
    FEId_low =  fespace_low->GetFE1D(i, Cell_low);
    Element_low = TFEDatabase2D::GetFE1D(FEId_low);
    bf_low = Element_low->GetBaseFunct1D(); 

    l = bf_low->GetPolynomialDegree();
    LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(2*l);
    qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
    qf1->GetFormulaData(N_LinePoints, LineWeights, zeta);

     F_K_low = TFEDatabase2D::GetRefTrans1D(LineAffin);
    ((TLineAffin *)F_K_low)->SetCell(Cell_low);

    ((TLineAffin *)F_K_low)->GetOrigFromRef(N_LinePoints, zeta, psi, AbsDetjk);

    //for(j=0;j<N_LinePoints;j++)
     //cout << " psi " << psi[j]<< endl;

    for(j=0;j<N_NonMotNeibs[i];j++)
     {
      k = NonMotCellLocalIndex[i][j];

      for(l=0;l<N_LinePoints;l++)
      {
       N = N_NeededQuadPtValues[k];
       QuadPts[k][N] = psi[l];
       N_NeededQuadPtValues[k]++;
      }
     }
   }

//     for(i=0;i<N_NonMortarEdges;i++)
//     { for(l=0;l<N_NeededQuadPtValues[i]; l++)
//         cout << "QuadPts " << QuadPts[i][l] << endl;
//       cout <<  endl;
//     }

  for(i=0;i<N_NonMortarEdges;i++)
   {
    M = NonMotCells[i];
    IJoint = NonMotEdges[i];
    cell = Coll->GetCell(M);     // ansatz

// ==========================================================================================
    // manupulate data for line quad formulae
    k=cell->GetN_Edges();
    cell->GetVertex(IJoint)->GetCoords(startX, y);
    cell->GetVertex( (IJoint+1) %k)->GetCoords(endX, y);

    // cout << "start " << startX  << " endX " << endX  << endl;
     N_LinePoints = 0;
     for(l=0;l<N_NeededQuadPtValues[i]; l++)
      {
       x=QuadPts[i][l];
       if( x>=startX && x<=endX  )
        {
         DummyWeights[N_LinePoints] = 0;
         DummyZeta[N_LinePoints] = (2.*x - startX - endX)/(endX - startX);
         // cout << "start " << startX  << " endX " << endX  <<  " X " << x 
         // << " DummyZeta[N_LinePoints] " << DummyZeta[N_LinePoints]  << endl;
         N_LinePoints++;
        }
       else if( x<=startX && x>=endX )
        {
         DummyWeights[N_LinePoints] = 0;
         DummyZeta[N_LinePoints] = (2.*x - startX - endX)/(startX - endX);
         cout << "start " << startX  << " endX " << endX  <<  "L " << l 
         << " DummyZeta[N_LinePoints] " << DummyZeta[N_LinePoints]  << endl;
         N_LinePoints++;
        }
       else // mortar point is outside of the cell
        { NonMotBasisVal[i][l*MaxN_BaseFunctions2D] = -1e8;  } 

      } // for(l=0;l<N_NeededQuadPtValues[i]; l++)


    if(i) delete Dummyqf1;

    Dummyqf1 = new TQuadFormula1D(N_LinePoints, DummyWeights,  DummyZeta, 0);
    TFEDatabase2D::RegisterQuadFormula1D(Dummy, Dummyqf1); 
// =============================================================================================

    FEId = fespace->GetFE2D(M, cell);
    RefElement = TFEDatabase2D::GetRefElementFromFE2D(FEId);
    FeDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(FEId);



    // set reference element
    switch(RefElement)
     {
      case BFUnitSquare:
          RefTrans = QuadBilinear;
          F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
          ((TQuadBilinear *)F_K)->SetCell(cell);
      break;

      case BFUnitTriangle:
          RefTrans = TriaAffin;
          F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
          ((TTriaAffin *)F_K)->SetCell(cell);
          //cout << "Cell RefElement" << RefElement  << endl;
      break;
    } // endswitch

    cout << " Mortar   Test " << endl;
    TFEDatabase2D::GetBaseFunct2DFromFE2D(FEId)->MakeRefElementData(Dummy); 


    // get the isoedge joint dof
    JointDOF = FeDesc->GetJointDOF(IJoint);
    N_JointDOF = FeDesc->GetN_JointDOF(); 

    uref = TFEDatabase2D::GetJointValues2D(BaseFuncts[FEId], Dummy, IJoint);   


   k=0; 
   for(j=0;j<N_LinePoints;j++)
    {
     while(NonMotBasisVal[i][k*MaxN_BaseFunctions2D] == -1e8 ) k++;

     uorig = uref[j];
     for(l=0;l<N_JointDOF;l++) // ansatz
      {
       locdof = JointDOF[l];
        NonMotBasisVal[i][k*MaxN_BaseFunctions2D + l] = uorig[locdof];
        // cout<< j << " k " << k << " ansatz0 " << uorig[locdof] << endl;
      }
     k++;
    }

   } //  for(i=0;i<N_NonMortarEdges;i++)



  // assume same FE is used in all cells of NonMot domain
  // so that N_JointDOF is same in all cells
  cell = Coll->GetCell(0);
  FEId = fespace->GetFE2D(0, cell);
  FeDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(FEId);
  N_JointDOF = FeDesc->GetN_JointDOF();

  //collected all info, now assemble the matrix
  for(i=0;i<N_Cells_low;i++)
   {
    Cell_low = Coll_low->GetCell(i);
    FEId_low =  fespace_low->GetFE1D(i, Cell_low);
    Element_low = TFEDatabase2D::GetFE1D(FEId_low);
    bf_low = Element_low->GetBaseFunct1D(); 
    N_BaseFunct_low = Element_low->GetN_DOF();
    BaseFunct_ID_low = Element_low->GetBaseFunct1D_ID();


    l = bf_low->GetPolynomialDegree();
    LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(2*l);
    qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
    qf1->GetFormulaData(N_LinePoints, LineWeights, zeta);


     F_K_low = TFEDatabase2D::GetRefTrans1D(LineAffin);
    ((TLineAffin *)F_K_low)->SetCell(Cell_low);

    ((TLineAffin *)F_K_low)->GetOrigFromRef(N_LinePoints, zeta, psi, AbsDetjk);

    //for(j=0;j<N_LinePoints;j++)
     //cout << " psi " << psi[j]<< endl;

    BaseFunct_low[0] = BaseFunct_ID_low;
    ((TLineAffin *)F_K_low)->GetOrigValues(N_Sets, BaseFunct_low, N_LinePoints, zeta,
                                           LineQuadFormula,  Needs2ndDer);
    origvaluesD0=TFEDatabase2D::GetOrigElementValues(BaseFunct_ID_low, D0);

    // get the 2D cell dof on 1D interface
    N_NonMotNeibsInCell = N_NonMotNeibs[i];

    memset(LocMatrixB, 0, N_BaseFunct_low*N_JointDOF*N_NonMotNeibsInCell*SizeOfDouble);
    memset(LocMatrixBT, 0, N_BaseFunct_low*N_JointDOF*N_NonMotNeibsInCell*SizeOfDouble);

    // assemble  Mat
    for(j=0;j<N_LinePoints;j++)
     {
      x = psi[j];
      //cout<<" i  " << i << " x  " << x <<endl;
      // find this point "x" in the corresponding NonMot edges associated with this 1D cell
      for(k=0;k<N_NonMotNeibsInCell ;k++)
       {
        l = NonMotCellLocalIndex[i][k];
        N = N_NeededQuadPtValues[l];

        for(m=0;m<N;m++)
         if(QuadPts[l][m] == x && NonMotBasisVal[l][m*MaxN_BaseFunctions2D] != -1e8)
          {
           //cout<< " i  " << i <<" l  " << l  << " x " << x << " m " << m <<endl;
           Mult = LineWeights[j]*AbsDetjk[j];
           orgD0 = origvaluesD0[j];
           uorig = NonMotBasisVal[l] + m*MaxN_BaseFunctions2D;
           disp = k*N_JointDOF*N_BaseFunct_low;

           // B Mat
           for(kk=0;kk<N_BaseFunct_low;kk++) // test
            {
             test0  = orgD0[kk];
             for(ll=0;ll<N_JointDOF;ll++) // ansatz
              {
               ansatz0  = uorig[ll];
               val  = ansatz0*test0;
               LocMatrixB[disp + kk*N_BaseFunct_low + ll] += (Mult*val);
               //cout<< "B ( " << kk <<" , " <<ll << " ) " << Mult*val << endl;
              } // 
            } // for(k=0;k<N_BaseFunct_lo

           // BMatT
           for(kk=0;kk<N_JointDOF;kk++) // ansatz
            {
             test0  = uorig[kk];
             for(ll=0;ll<N_BaseFunct_low;ll++) // ansatz
              {
               ansatz0  = orgD0[ll];
               val  = ansatz0*test0;
               LocMatrixBT[disp + kk*N_JointDOF + ll] += (Mult*val);
               //cout<< "BT ( " << kk <<" , " <<ll << " ) " << Mult*val << endl;
              } // 
             }// for(ll=0;ll<N_JointDOF;ll++) // ansatz

          }
       } //    for(k=0;k<N_NonMotNeibs[i];k++)
     } //N_LinePoints

   // get 1D cell dof
   DOF_low = GlobalNumbers_low +  BeginIndex_low[i];

   for(j=0;j<N_NonMotNeibsInCell;j++)
    {
     N = NonMotNeibs[10*i + j];
     IJoint = NonMotNeibs[10*i + 5 + j];
     cell = Coll->GetCell(N);
     FEId = fespace->GetFE2D(N, cell);
     FeDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(FEId);
     JointDOF = FeDesc->GetJointDOF(IJoint);
     N_JointDOF = FeDesc->GetN_JointDOF();
     DOF = GlobalNumbers + BeginIndex[N];

      for(k=0;k<N_JointDOF;k++)
       {
        dof = DOF[JointDOF[k]];

        ll = GetIndex(GlobalDof, N_NonMortaDofs, dof);
        DOF_NonMot[j*N_JointDOF + k ] = ll;
        AssembleB_NonMat[dof] = ll;
       }
    }

  // add to global matrices B
   for(j=0;j<N_BaseFunct_low;j++)
     {
      locdof = DOF_low[j];
      begin = RowPtr[locdof];
      end = RowPtr[locdof+1];
//       cout<< "locdof " << locdof;
      for(k=begin;k<end;k++)
       {
//         cout<< " KCol[k] " << KCol[k] <<endl;
        for(kk=0;kk<N_NonMotNeibsInCell;kk++)
         {
          UPDATED = FALSE;
          disp = kk*N_JointDOF*N_BaseFunct_low;
          for(l=0;l<N_JointDOF;l++)
           {
            if(KCol[k] == DOF_NonMot[kk*N_JointDOF + l])
             {
              ValuesB[k] +=LocMatrixB[disp + j*N_BaseFunct_low + l];
              UPDATED = TRUE;
//               cout<< "B ( " << locdof <<" , " << DOF_NonMot[l] << " ) " << ValuesB[k] << endl;
               break;
              }
            } // for(l=0;l<N_JointDOF;l++)
//            if(UPDATED) break;
         } // for(kk=0;kk<N_NonMotNeibsInCell;kk++)
      } // for(k=begin;k<end;k++)
     }

  // add to global matrices BT
  for(kk=0;kk<N_NonMotNeibsInCell;kk++)
   {
    UPDATED = FALSE;
    disp = kk*N_JointDOF*N_BaseFunct_low;

    for(j=0;j<N_JointDOF;j++)
     {
      locdof = DOF_NonMot[kk*N_JointDOF + j];
      begin = RowPtrT[locdof];
      end = RowPtrT[locdof+1];
      for(k=begin;k<end;k++)
       {
        for(l=0;l<N_BaseFunct_low;l++)
         {
          if(KColT[k] == DOF_low[l])
           {
            ValuesBT[k] +=LocMatrixBT[disp + j*N_JointDOF + l];
            //out<< "BT ( " << locdof <<" , " << DOF_low[l] << " ) " << ValuesBT[k] << endl;
            break;
           }
         } // for(l=0;l<N_BaseFunct_low;l++)
        } // for(k=begin;k<end;k++)
      } //for(j=0;j<N_JointDOF;j++)
     } //for(kk=0;kk<N_NonMotNeibsInCell;kk++)
  } //for(i=0;i<N_Cells_low;i++)


//   //print matrix
//   N = fespace_low->GetN_DegreesOfFreedom(); 
//   for(j=0;j<N;j++)
//    {
//     begin = RowPtr[j];
//     end = RowPtr[j+1];
//     for(k=begin;k<end;k++)
//      {
//       cout << "B(" << j << ", "<< KCol[k] << ") = " << ValuesB[k] <<endl;
//      }
//     cout<<endl;
//    }

//   //print matrix
//   N = NonMortarFEData->N_NonMortaDofs;
//   for(j=0;j<N;j++)
//    {
//     begin = RowPtrT[j];
//     end = RowPtrT[j+1];
//     for(k=begin;k<end;k++)
//      {
//       cout << "BT(" << j << ", "<< KColT[k] << ") = " << ValuesBT[k] <<endl;
//      }
//     cout<<endl;
//    }
// 

  
} // AssembleB_NonMat


void SolveMortarSystem(TSquareMatrix2D *sqmatrixA_Mot, TSquareMatrix2D *sqmatrixA_NonMot,
        TMatrix2D *B_Mot, TMatrix2D *B_NonMot, TMatrix2D *B_NonMotT,
        int *MotLocalGlobalDof, int *NonMotLocalGlobalDof,
        TNonMortarData *NonMortarFEData, int *GlobalDof_MortaEdgeDof,
        double *sol_all, double *rhs_all)
{
 int i, j, l, k, N_Mot, N_NonMot, N_L, N_, ret;
 int *KColA_Mot, *RowPtrA_Mot, *KColA_NonMot, *RowPtrA_NonMot, *KColB_Mot, *RowPtrB_Mot;
 int *KColB_NonMot, *RowPtrB_NonMot, *KColB_NonMotT, *RowPtrB_NonMotT;
 int pos, begin, end, N_Entries, *KCol, *RowPtr, row, *GlobalDof_NonMortaEdgeDof;

 double t1, t2, t3, t4;
 double *Values;
 double *EntriesA_Mot, *EntriesA_NonMot, *EntriesB_Mot, *EntriesB_NonMot, *EntriesB_NonMotT;
 double *Entries, value;

  void *Symbolic, *Numeric;

  N_Mot = sqmatrixA_Mot->GetN_Rows();
  N_NonMot = sqmatrixA_NonMot->GetN_Rows();
  N_L = B_Mot->GetN_Rows();
  N_ = N_Mot + N_NonMot + N_L;

  cout<< "Total N_Unkowns " << N_ << endl;

  KColA_Mot = sqmatrixA_Mot->GetKCol();
  RowPtrA_Mot = sqmatrixA_Mot->GetRowPtr();

  KColA_NonMot = sqmatrixA_NonMot->GetKCol();
  RowPtrA_NonMot = sqmatrixA_NonMot->GetRowPtr();

  KColB_Mot = B_Mot->GetKCol();
  RowPtrB_Mot = B_Mot->GetRowPtr();

  KColB_NonMot = B_NonMot->GetKCol();
  RowPtrB_NonMot = B_NonMot->GetRowPtr();

  KColB_NonMotT = B_NonMotT->GetKCol();
  RowPtrB_NonMotT = B_NonMotT->GetRowPtr();

  EntriesA_Mot = sqmatrixA_Mot->GetEntries();
  EntriesA_NonMot = sqmatrixA_NonMot->GetEntries();

  EntriesB_Mot = B_Mot->GetEntries();
  EntriesB_NonMot = B_NonMot->GetEntries();
  EntriesB_NonMotT = B_NonMotT->GetEntries();

  // allocate arrays for structure of combined matrix
  // total number of entries
  N_Entries = RowPtrA_Mot[N_Mot] + RowPtrA_NonMot[N_NonMot] + 2*RowPtrB_Mot[N_L]
                + 2*RowPtrB_NonMot[N_L];
  Entries = new double[N_Entries];
  KCol = new int[N_Entries];
  RowPtr = new int[N_+1];
  RowPtr[0] = 0;

  pos = 0;

  // fill combined matrix
  for(i=0;i<N_Mot;i++)
  {
    begin = RowPtrA_Mot[i];
    end = RowPtrA_Mot[i+1];

    for(j=begin;j<end;j++)
    {
      Entries[pos] = EntriesA_Mot[j];
      KCol[pos] = KColA_Mot[j];
      pos++;
    }

   if(TDatabase::ParamDB->P0==1)
    {
    // lagrangian part
    row = MotLocalGlobalDof[i];
    if(row>=0)
     {
      //cout << "row " << row << endl;
      begin = RowPtrB_Mot[row];
      end = RowPtrB_Mot[row+1];

      for(j=begin;j<end;j++)
       {
        Entries[pos] = EntriesB_Mot[j];
        KCol[pos] = N_Mot + N_NonMot + KColB_Mot[j];
        pos++;
       }
     }
    }
    RowPtr[i+1] = pos;
  } //  for(i=0;i<N_Mot;i++)

 //nonmortar part
  for(i=0;i<N_NonMot;i++)
   {
    begin = RowPtrA_NonMot[i];
    end = RowPtrA_NonMot[i+1];

    for(j=begin;j<end;j++)
     {
      Entries[pos] = EntriesA_NonMot[j];
      KCol[pos] = N_Mot + KColA_NonMot[j];
      pos++;
     }

   if(TDatabase::ParamDB->P0==1)
    {
    // lagrangian part
    row = NonMotLocalGlobalDof[i];
    if(row>=0)
     {
      //cout << "row " << row << endl;
      begin = RowPtrB_NonMotT[row];
      end = RowPtrB_NonMotT[row+1];

      for(j=begin;j<end;j++)
       {
        Entries[pos] = -EntriesB_NonMotT[j];
        KCol[pos] = N_Mot + N_NonMot + KColB_NonMotT[j];
        pos++;
       }
     }
    }// if(TDatabase::ParamDB->P0==1)
    RowPtr[N_Mot + i+1] = pos;
   } // for(i=0;i<N_NonMot;i++)

  GlobalDof_NonMortaEdgeDof =  NonMortarFEData->GlobalDof_NonMortaEdgeDof;

 if(TDatabase::ParamDB->P0==1)
  {
  // lagrangian part
  for(i=0;i<N_L;i++)
   {
    begin = RowPtrB_Mot[i];
    end = RowPtrB_Mot[i+1];

    for(j=begin;j<end;j++)
     {
      Entries[pos] = EntriesB_Mot[j];
      KCol[pos] = GlobalDof_MortaEdgeDof[ KColB_Mot[j] ];
      //cout << "col " << KCol[pos]  << endl;
      pos++;
     }

    begin = RowPtrB_NonMot[i];
    end = RowPtrB_NonMot[i+1];

    for(j=begin;j<end;j++)
     {
      Entries[pos] = -EntriesB_NonMot[j];
      KCol[pos] = N_Mot + GlobalDof_NonMortaEdgeDof[ KColB_NonMot[j] ];
      //cout << "col " << KCol[pos]  << endl;
      pos++;
     }
    RowPtr[N_Mot + N_NonMot + i+1] = pos;
   }
  } //if(TDatabase::ParamDB->P0==1)
//   N_ =  N_Mot + N_NonMot  ;

   if(TDatabase::ParamDB->P0!=1)
     N_ = N_Mot + N_NonMot;



  // sort matrix
  for(i=0;i<N_;i++)
  {
    begin=RowPtr[i];
    end=RowPtr[i+1];

    for(j=begin;j<end;j++)
    {
      for(k=j+1;k<end;k++)
      {
        if(KCol[j] > KCol[k])
        {
          l = KCol[j];      value = Entries[j];
          KCol[j] = KCol[k]; Entries[j] = Entries[k];
          KCol[k] = l;       Entries[k] = value;
        }                        // endif
      }                          // endfor k
    }                            // endfor j
  }                              // endfor i

  /*
  for(i=0;i<N_;i++)
  {
    for(j=RowPtr[i];j<RowPtr[i+1];j++)
      cout << i << " " << KCol[j] << " " << Entries[j] << endl;
  }
  */


  t1 = GetTime();
  ret = umfpack_di_symbolic(N_, N_, RowPtr, KCol, Entries,
    &Symbolic, NULL, NULL);
  t2 = GetTime();

  // error occured
  if (ret!=0)
  {
    OutPut("error in umfpack_di_symbolic " << ret << endl);
    exit(4711);
  }

  ret = umfpack_di_numeric(RowPtr, KCol, Entries, Symbolic,
    &Numeric, NULL, NULL);
  umfpack_di_free_symbolic(&Symbolic);
  t3 = GetTime();

  // error occured
  if (ret!=0)
  {
    OutPut("error in umfpack_di_numeric " << ret << endl);
    exit(4711);
  }

  ret = umfpack_di_solve(UMFPACK_At, RowPtr, KCol, Entries,
    sol_all, rhs_all, Numeric, NULL, NULL);
  umfpack_di_free_numeric(&Numeric);
  t4 = GetTime();
  if (ret!=0)
  {
    OutPut("error in umfpack_di_solve " << ret << endl);
    exit(4711);
  }
//  OutPut("umfpack: " << ret << " " << t4-t1 << " sec." << endl);

} // SolveMortarSystem


int main(int argc, char* argv[])
{
  TDomain *Domain = new TDomain();
  TDomain *Domain_NonMot = new TDomain();
  TDomain *Domain_Mot = new TDomain();
  TDomain *Domain_MotIntFace = new TDomain();

  TDatabase *Database = new TDatabase();
  TFEDatabase2D *FEDatabase = new TFEDatabase2D();
  TCollection *Coll,  *Coll_NonMot, *Coll_Mot, *MotIntFace_Coll;
  TFESpace2D *FESpace, *FESpace_Mot, *FESpace_NonMot;
  TFESpace1D *FeSpace_MortIntFace;
  TFESpace2D *fesp[2], *ferhs[3];
  TOutput2D *Output;
  TAuxParam2D *aux;
  TFEFunction2D *FEFunction_Mot, *FEFunction_NonMot, *FEFunction;
  TFEFunction1D *FEFunct_MortIntFace;

  TSquareStructure2D *sqstructureA_Mot, *sqstructureA_NonMot;
  TSquareMatrix2D *sqmatrixA_Mot, *SQMATRICES[1], *sqmatrixA_NonMot;
  TSquareMatrix2D **MatricesA_Mot, **MatricesA_NonMot; 
  TSquareMatrix **sqmatrices = (TSquareMatrix **)SQMATRICES;
  TStructure2D *Structure_Mot, *Structure_NonMot, *Structure_NonMotT;
  TMatrix2D *B_Mot, *B_NonMot, *B_NonMotT;
  
  TBaseCell **Cells;

  TNonMortarData NonMortarFEData;

  int i, j, ret, ORDER, N_DOF_Mot, N_DOF_NonMot, N_DOF;
  int N_Cells, N_Cells_Mot, N_Cells_NonMot, img=1, *MortarCellEdge[2], *NonMortarCellEdge[2];
  int N_MotIntFaceCells, N_DOF_MortIntFace, N_Unknowns, *GlobalDof_NonMortaEdgeDof;
  int *N_NonMotNeibs, *NonMotNeibs, N_NonMortarEdges, N_NonMortaDofs;
  int N_Active_Mot, N_Active_NonMot;
  int *GlobalNo, *BeginIndex, *GlobalNo_Mot, *BeginIndex_Mot, *GlobalNo_NonMot, *BeginIndex_NonMot;

  double *sol_Mot, *oldsol_Mot, *rhs_Mot, *defect_Mot;
  double *sol_NonMot, *oldsol_NonMot, *rhs_NonMot, *defect_NonMot;
  double *sol, *RHSs[1], *sol_MortIntFace, *rhs_MortIntFace;

  char *PRM_Mot, *GEO_Mot, *PRM_NonMot, *GEO_NonMot;
  char *PsBaseName, *VtkBaseName;
  char CdString[] = "Conv-Diff";
  char CMotString[] = "C_Mot";
  char CNonMotString[] = "C_NonMot";
  char CString[] = "C";
  char LString[] = "Lambda"; 
  char GalString[] = "Galerkin";
  char ReadinDat [] = "readin.dat";
  char Name[] = "name";
  char Description[] = "description";

  std::ostringstream os, opts;
  os << " ";
  opts << " ";
//======================================================================
//read parameter file Readin.dat
//======================================================================
  if(argc>=2)
    ret=Domain_Mot->ReadParam(argv[1]);
  else
    ret=Domain_Mot->ReadParam(ReadinDat);

 if(argc>=2)
    ret=Domain_NonMot->ReadParam(argv[1]);
  else
    ret=Domain_NonMot->ReadParam(ReadinDat);

  if(ret==-1)
    exit(-1);

  OpenFiles();
  OutFile.setf(std::ios::scientific);

  Database->WriteParamDB(argv[0]);
  ExampleFile();

//======================================================================
// include the boundary condition and boundary values from the example file
//======================================================================
  BoundCondFunct2D *BoundaryConditions_Mot[1] = { BoundCondition_Mortar };
  BoundValueFunct2D *BoundaryValues_Mot[1] = { BoundValue_Mortar };
  BoundCondFunct2D *BoundaryConditions_NonMot[1] = { BoundCondition_NonMortar };
  BoundValueFunct2D *BoundaryValues_NonMot[1] = { BoundValue_NonMortar }; 
//======================================================================
// copy read parameters into local variables
//======================================================================
  PRM_Mot = TDatabase::ParamDB->BNDFILE;
  GEO_Mot = TDatabase::ParamDB->GEOFILE;

  PRM_NonMot = TDatabase::ParamDB->GRAPEBASENAME;
  GEO_NonMot = TDatabase::ParamDB->GNUBASENAME;

  PsBaseName = TDatabase::ParamDB->PSBASENAME;
  VtkBaseName = TDatabase::ParamDB->VTKBASENAME;

//======================================================================
// define discrete form
//======================================================================
  TDiscreteForm2D *DiscreteForm, *DiscreteForms[10];
  TDiscreteForm2D *DiscreteFormGalerkin = new TDiscreteForm2D
    (CdString, GalString, N_Terms, Derivatives, SpacesNumbers,
    N_Matrices, N_Rhs, RowSpace, ColumnSpace, RhsSpace,
    BilinearAssemble, BilinearCoeffs, NULL);

//======================================================================
// generate mesh  
//======================================================================
  Domain_Mot->Init(PRM_Mot, GEO_Mot);
  Domain_NonMot->Init(PRM_NonMot, GEO_NonMot);
//======================================================================
// collect collection of cells
//======================================================================
  /** Get Coll_Mot */
  TDatabase::IteratorDB[It_EQ]->SetParam(Domain_Mot);
  TDatabase::IteratorDB[It_LE]->SetParam(Domain_Mot);
  TDatabase::IteratorDB[It_Finest]->SetParam(Domain_Mot);
  TDatabase::IteratorDB[It_Between]->SetParam(Domain_Mot);
  TDatabase::IteratorDB[It_OCAF]->SetParam(Domain_Mot);

  for(i=0;i<TDatabase::ParamDB->UNIFORM_STEPS;i++)
   Domain_Mot->RegRefineAll();

  Coll_Mot = Domain_Mot->GetCollection(It_Finest, 0);
  N_Cells_Mot = Coll_Mot->GetN_Cells();

  // write grid into an Postscript file
  os.seekp(std::ios::beg);
  os << "Domain_Mot" << ".ps" << ends;
  Domain_Mot->PS(os.str().c_str(),It_Finest,0);

  /**generate 1D mortar interface domain */
  GenMortarInterface(Domain_Mot, 0, Domain_MotIntFace, MortarCellEdge[0], MortarCellEdge[1], 0., 0.);
  MotIntFace_Coll = Domain_MotIntFace->GetCollection(It_Finest, 0);
  N_MotIntFaceCells= MotIntFace_Coll->GetN_Cells();

  // for(i=0;i<N_MotIntFaceCells;i++)
  // cout << i << " Cell_No " <<  MortarCellEdge[0][i] <<" Joint_Nr " << MortarCellEdge[1][i] << endl;

  /** Get Coll_NonMot */
  TDatabase::IteratorDB[It_EQ]->SetParam(Domain_NonMot);
  TDatabase::IteratorDB[It_LE]->SetParam(Domain_NonMot);
  TDatabase::IteratorDB[It_Finest]->SetParam(Domain_NonMot);
  TDatabase::IteratorDB[It_Between]->SetParam(Domain_NonMot);
  TDatabase::IteratorDB[It_OCAF]->SetParam(Domain_NonMot);

  for(i=0;i<TDatabase::ParamDB->P9;i++)
   Domain_NonMot->RegRefineAll();

  Coll_NonMot = Domain_NonMot->GetCollection(It_Finest, 0);
  N_Cells_NonMot = Coll_NonMot->GetN_Cells();

  os.seekp(std::ios::beg);
  os << "Domain_NonMot" << ".ps" << ends;
  Domain_NonMot->PS(os.str().c_str(),It_Finest,0);

  cout << "N_Cells_Mot " << N_Cells_Mot << " N_Cells_NonMot " << N_Cells_NonMot<<
          " N_MotIntFaceCells " << N_MotIntFaceCells << endl;

  /** find match for Domain_NonMot and  Domain_MotIntFace */
  Map2DNonMot_1DMotIntFace(Coll_NonMot, 2, MotIntFace_Coll, NonMortarCellEdge[0],
                           NonMortarCellEdge[1], 0., 0., N_NonMotNeibs, NonMotNeibs,
                           N_NonMortarEdges);

   /** put both coll together */
  N_Cells = N_Cells_Mot+N_Cells_NonMot;
  Cells = new TBaseCell*[N_Cells];

  j=0;
  for (i=0;i<N_Cells_Mot;i++)
   Cells[j++] = Coll_Mot->GetCell(i);

  for (i=0;i<N_Cells_NonMot;i++)
   Cells[j++] = Coll_NonMot->GetCell(i);

  Domain->SetTreeInfo(Cells, N_Cells);

  TDatabase::IteratorDB[It_EQ]->SetParam(Domain);
  TDatabase::IteratorDB[It_LE]->SetParam(Domain);
  TDatabase::IteratorDB[It_Finest]->SetParam(Domain);
  TDatabase::IteratorDB[It_Between]->SetParam(Domain);
  TDatabase::IteratorDB[It_OCAF]->SetParam(Domain);

  Coll = Domain->GetCollection(It_Finest, 0);
//========================================================================================
// Construct all FESpaces
//========================================================================================
  ORDER = TDatabase::ParamDB->ANSATZ_ORDER;

  FESpace_Mot =  new TFESpace2D(Coll_Mot, Name, Description, BoundCondition_Mortar, ORDER, NULL);

  N_DOF_Mot =  FESpace_Mot->GetN_DegreesOfFreedom();
  N_Active_Mot = FESpace_Mot->GetActiveBound();
  GlobalNo_Mot = FESpace_Mot->GetGlobalNumbers();
  BeginIndex_Mot = FESpace_Mot->GetBeginIndex();


  OutPut(setw(20) <<"Mortar   : "  << N_DOF_Mot  << endl);

  FESpace_NonMot = new TFESpace2D(Coll_NonMot, Name, Description, BoundCondition_NonMortar, ORDER, NULL);

  N_DOF_NonMot = FESpace_NonMot->GetN_DegreesOfFreedom();
  N_Active_NonMot = FESpace_NonMot->GetActiveBound();
  GlobalNo_NonMot = FESpace_NonMot->GetGlobalNumbers();
  BeginIndex_NonMot = FESpace_NonMot->GetBeginIndex();

  OutPut(setw(20) <<"NonMortar   : "  <<N_DOF_NonMot  << endl);

  FeSpace_MortIntFace = new TFESpace1D(MotIntFace_Coll, LString, LString, ORDER);
  N_DOF_MortIntFace =  FeSpace_MortIntFace->GetN_DegreesOfFreedom();
  OutPut(setw(20) <<"MortIntFace   : "  <<N_DOF_MortIntFace  << endl);

  N_Unknowns = N_DOF_Mot + N_DOF_NonMot + N_DOF_MortIntFace;
  OutPut(setw(20) <<"N_Unknowns   : "  <<N_Unknowns  << endl);

  /** for output */
  FESpace = new TFESpace2D(Coll, Name, Description, BoundCondition_Mortar, ORDER, NULL);
  N_DOF = FESpace->GetN_DegreesOfFreedom();
  GlobalNo = FESpace->GetGlobalNumbers();
  BeginIndex = FESpace->GetBeginIndex();
// ======================================================================
// construct all finite element functions
// ======================================================================
    sol_Mot = new double[N_DOF_Mot];
    oldsol_Mot = new double[N_DOF_Mot];
    rhs_Mot = new double[N_DOF_Mot];
    defect_Mot = new double[N_DOF_Mot];

    memset(sol_Mot, 0, N_DOF_Mot*SizeOfDouble);
    memset(oldsol_Mot, 0, N_DOF_Mot*SizeOfDouble);
    memset(rhs_Mot, 0, N_DOF_Mot*SizeOfDouble);

    FEFunction_Mot = new TFEFunction2D(FESpace_Mot, CMotString, CMotString, sol_Mot, N_DOF_Mot);
//     FEFunction_Mot->Interpolate(Exact)

    sol_NonMot = new double[N_DOF_NonMot];
    oldsol_NonMot = new double[N_DOF_NonMot];
    rhs_NonMot = new double[N_DOF_NonMot];
    defect_NonMot = new double[N_DOF_NonMot];

    memset(sol_NonMot, 0, N_DOF_NonMot*SizeOfDouble);
    memset(oldsol_NonMot, 0, N_DOF_NonMot*SizeOfDouble);
    memset(rhs_NonMot, 0, N_DOF_NonMot*SizeOfDouble);

    FEFunction_NonMot = new TFEFunction2D(FESpace_NonMot, CNonMotString, CNonMotString, 
                                          sol_NonMot, N_DOF_NonMot);

    sol = new double[N_DOF];
    memset(sol, 0, N_DOF*SizeOfDouble);

    FEFunction = new TFEFunction2D(FESpace, CString, CString, sol, N_DOF);


    sol_MortIntFace = new double[N_DOF_MortIntFace];
    memset(sol_MortIntFace, 0, N_DOF_MortIntFace*SizeOfDouble);
    
    // fefunction for lagrangian multipliers
    FEFunct_MortIntFace = new TFEFunction1D(FeSpace_MortIntFace, LString, LString, sol_MortIntFace, N_DOF_MortIntFace);
        
//======================================================================
// allocate memory for all matrices
//======================================================================

    sqstructureA_Mot = new TSquareStructure2D(FESpace_Mot);
    sqstructureA_Mot->Sort();
    sqmatrixA_Mot = new TSquareMatrix2D(sqstructureA_Mot);   

    sqstructureA_NonMot = new TSquareStructure2D(FESpace_NonMot);
    sqstructureA_NonMot->Sort();
    sqmatrixA_NonMot = new TSquareMatrix2D(sqstructureA_NonMot);        
    

   // B_Mot matrix and B_Mot transpose are same in the present implementation
    Structure_Mot = new TStructure2D(FeSpace_MortIntFace, FESpace_Mot, MortarCellEdge);
    B_Mot = new TMatrix2D(Structure_Mot);
    

    GetN_NonMortaDof(FeSpace_MortIntFace, FESpace_NonMot, N_NonMortarEdges, NonMortarCellEdge,
		     N_NonMotNeibs, NonMotNeibs, &NonMortarFEData);
    // B_NonMot block
    Structure_NonMot = new TStructure2D(FeSpace_MortIntFace, FESpace_NonMot, &NonMortarFEData);
    B_NonMot = new TMatrix2D(Structure_NonMot);
    
    // B_NonMotT block
    Structure_NonMotT = new TStructure2D(FESpace_NonMot, FeSpace_MortIntFace, &NonMortarFEData);
    B_NonMotT = new TMatrix2D(Structure_NonMotT);   
//======================================================================
// // assemble 2D M matrix - begin
//======================================================================

    // first mortar matrices
    RHSs[0] = rhs_Mot;
    memset(rhs_Mot, 0, N_DOF_Mot*SizeOfDouble);

    fesp[0] = FESpace_Mot;
    ferhs[0] = FESpace_Mot; 

    DiscreteForm = DiscreteFormGalerkin;

      // initialize matrices
    SQMATRICES[0] = sqmatrixA_Mot;
    SQMATRICES[0]->Reset();
    aux =  new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);

      // assemble
    Assemble2D(1, fesp,
               1, SQMATRICES,
               0, NULL,
               1, RHSs, ferhs,
               DiscreteForm,
               BoundaryConditions_Mot,
               BoundaryValues_Mot,
               aux);

     delete aux;

   memcpy(sol_Mot+N_Active_Mot,  rhs_Mot+N_Active_Mot, (N_DOF_Mot - N_Active_Mot)*SizeOfDouble);


  //assemble B_mot matrix
   int *MotLocalGlobalDof, *GlobalDof_MortaEdgeDof;
   MotLocalGlobalDof = new int[N_DOF_Mot];
   GlobalDof_MortaEdgeDof = new int[N_DOF_MortIntFace];

   for (i=0;i<N_DOF_Mot;i++)
    MotLocalGlobalDof[i] = -1;


   B_Mot->Reset();
   AssembleB_Mat(B_Mot, FEFunct_MortIntFace, FESpace_Mot, MortarCellEdge, MotLocalGlobalDof, 
                 GlobalDof_MortaEdgeDof);


   //assemble B_NonMot matrix
   int *NonMotLocalGlobalDof;
   NonMotLocalGlobalDof = new int[N_DOF_NonMot];

   for (i=0;i<N_DOF_NonMot;i++)
    NonMotLocalGlobalDof[i] = -1;

   B_NonMot->Reset();
   B_NonMotT->Reset(); 

   AssembleB_NonMat(B_NonMot, B_NonMotT, FEFunct_MortIntFace, FESpace_NonMot, N_NonMortarEdges, NonMortarCellEdge, 
		    N_NonMotNeibs, NonMotNeibs, &NonMortarFEData, NonMotLocalGlobalDof);

   cout << "Mortar Matries assembled " << endl;

    //nonmortar matrices
    RHSs[0] = rhs_NonMot;
    memset(rhs_NonMot, 0, N_DOF_NonMot*SizeOfDouble);

    fesp[0] = FESpace_NonMot;
    ferhs[0] = FESpace_NonMot; 

    DiscreteForm = DiscreteFormGalerkin;

      // initialize matrices
    SQMATRICES[0] = sqmatrixA_NonMot;
    SQMATRICES[0]->Reset();
    aux =  new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);

      // assemble
    Assemble2D(1, fesp,
               1, SQMATRICES,
               0, NULL,
               1, RHSs, ferhs,
               DiscreteForm,
               BoundaryConditions_NonMot,
               BoundaryValues_NonMot,
               aux);

     delete aux;

   memcpy(sol_NonMot+N_Active_NonMot,  rhs_NonMot+N_Active_NonMot, 
                       (N_DOF_NonMot - N_Active_NonMot)*SizeOfDouble);



    cout << "NonMortar Matries assembled " << endl;
    double *sol_all = new double[N_Unknowns];
    double *rhs_all = new double[N_Unknowns];
    memset(sol_all, 0, N_Unknowns*SizeOfDouble);
    memset(rhs_all, 0, N_Unknowns*SizeOfDouble);

    //mot part
    int N, *DOF_Low, *DOF;
    N = BeginIndex[1] - BeginIndex[0];
    for (i=0;i<N_Cells_Mot;i++)
     {
      DOF_Low= GlobalNo_Mot + BeginIndex_Mot[i];
      DOF= GlobalNo + BeginIndex[i];
      for (j=0;j<N;j++)
       sol[DOF[j]] = sol_Mot[DOF_Low[j]];
     }

    for (i=0;i<N_Cells_NonMot;i++)
     {
      DOF_Low= GlobalNo_NonMot + BeginIndex_NonMot[i];
      DOF = GlobalNo + BeginIndex[N_Cells_Mot + i];
      for (j=0;j<N;j++)
       sol[DOF[j]] = sol_NonMot[DOF_Low[j]];
     }

    memcpy(rhs_all, rhs_Mot,   N_DOF_Mot*SizeOfDouble);
    memcpy(rhs_all+N_DOF_Mot, rhs_NonMot,   N_DOF_NonMot*SizeOfDouble);


   SolveMortarSystem(sqmatrixA_Mot, sqmatrixA_NonMot, B_Mot, B_NonMot, B_NonMotT, 
                     MotLocalGlobalDof, NonMotLocalGlobalDof,
                     &NonMortarFEData, GlobalDof_MortaEdgeDof, sol_all, rhs_all);


    N = BeginIndex[1] - BeginIndex[0];
    for (i=0;i<N_Cells_Mot;i++)
     {
      DOF_Low= GlobalNo_Mot + BeginIndex_Mot[i];
      DOF= GlobalNo + BeginIndex[i];
      for (j=0;j<N;j++)
       sol[DOF[j]] = sol_all[DOF_Low[j]];
     }

    for (i=0;i<N_Cells_NonMot;i++)
     {
      DOF_Low= GlobalNo_NonMot + BeginIndex_NonMot[i];
      DOF = GlobalNo + BeginIndex[N_Cells_Mot + i];
      for (j=0;j<N;j++)
       sol[DOF[j]] = sol_all[N_DOF_Mot + DOF_Low[j]];
     }



//======================================================================
// construct output
//======================================================================
   Output = new TOutput2D(0, 1, 1, 1, Domain);
   Output->AddFEFunction(FEFunction);


    if(TDatabase::ParamDB->WRITE_VTK)
    {
      os.seekp(std::ios::beg);
       if(img<10) os << VtkBaseName<<".0000"<<img<<".vtk" << ends;
         else if(img<100) os << VtkBaseName<<".000"<<img<<".vtk" << ends;
          else if(img<1000) os << VtkBaseName<<".00"<<img<<".vtk" << ends;
           else if(img<10000) os << VtkBaseName<<".0"<<img<<".vtk" << ends;
            else  os << VtkBaseName<<"."<<img<<".vtk" << ends;
      Output->WriteVtk(os.str().c_str());
      img++;
     }
//======================================================================

  CloseFiles();
  OutPut("used time: " << GetTime() << endl);
  OutPut("used bytes: " << GetMemory() << endl);
  return 0;
}






 












