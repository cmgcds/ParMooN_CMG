// =======================================================================
// @(#)Edge.h  
// 
// Class:       TEdge
// Purpose:     superclass for edges in 3D
//
// Author:      Sashikumaar Ganesan  03.09.2010
//
// History:
//
// ======================================================================= 
#ifdef _MPI
#  include "mpi.h"
#endif

#include <Edge.h>
#include <string.h>
#include <BaseCell.h>
#include <Database.h>

/** constructor */

TEdge::TEdge(int n_Neibs, TBaseCell **neighbs)
{
 int i;

  N_Neibs = n_Neibs;
  Neighbs = new TBaseCell*[N_Neibs];

  for(i=0;i<N_Neibs;i++)
   Neighbs[i] = neighbs[i];

#ifdef  _MPI  
  SubDomainEdge = FALSE;
  SubDomainCrossEdge = FALSE;
  N_SubDomainNeibs = 0;
  SubDomainNeibsRank = NULL;
  CrossNeibsGlobalNo = NULL;
  N_CrossNeibs = 0;
  CrossNeibsLocalEdgeNo = NULL;
  CrossEdgeMaptype = NULL;
  CrossNeibsRank = NULL;
#endif 
 
}

#ifdef  _MPI
void TEdge::InitSubDomainInfo(int rank)
{
 int i, j, ID;
 bool UPDATE;

  if(N_Neibs<=0)
   {
    cout<<" TEdge::InitSubDomainInfo should be called after  Neib info added !!" <<endl;
//     MPI_Abort(MPI_COMM_WORLD, 0);
   }

  for(i=0; i<N_Neibs; i++)
   if(Neighbs[i]->GetSubDomainNo() != rank )
    {
     SubDomainEdge = TRUE;
     break;
    }

  N_SubDomainNeibs = 0;

  if(SubDomainEdge)
   {
    SubDomainNeibsRank = new int[N_Neibs];

    for(i=0; i<N_Neibs; i++)
     SubDomainNeibsRank[i] = -1;

    for(i=0; i<N_Neibs; i++)
     {
      ID = Neighbs[i]->GetSubDomainNo();
      if(ID==rank)
       continue;

     UPDATE = TRUE;
     for(j=0; j<N_SubDomainNeibs; j++)
      if(ID==SubDomainNeibsRank[j])
       {
        UPDATE = FALSE;
        break;
       }

     if(UPDATE)
      {
       SubDomainNeibsRank[N_SubDomainNeibs] = ID;
       N_SubDomainNeibs++;
      }
     } //  for(i=0; i<N_Neibs

    // now init the cross edge info (may not be cross edge but it will bw identified later)
    N_CrossNeibs = N_SubDomainNeibs;
    CrossNeibsRank = new int[N_CrossNeibs];

    for(j=0; j<N_CrossNeibs; j++)
     CrossNeibsRank[j] = SubDomainNeibsRank[j];

   } // if(SubDomainEd

}


void TEdge::SetAsNotCrossEdgeFor(int rank, int Neib_ID)
{
 int i;

  if(N_CrossNeibs>0)
   for(i=0; i<N_CrossNeibs; i++)
    if(CrossNeibsRank[i]==Neib_ID)
     {
      CrossNeibsRank[i] = CrossNeibsRank[N_CrossNeibs-1];
      N_CrossNeibs--;
      break;
     }

}

void TEdge::SetCrossNeibInfo(TVertex *Vert_a)
{
 int i, j, k, ID, Max, Max_Index, rank;
 const int *EdgeVertex;

 TBaseCell *NeibCell;

//  MPI_Comm Comm = TDatabase::ParamDB->Comm;
//  MPI_Comm_rank(Comm, &rank);

 if(N_CrossNeibs>0)
  {
   SubDomainCrossEdge = TRUE;
   CrossNeibsGlobalNo = new int[N_CrossNeibs];
   CrossNeibsLocalEdgeNo = new int[N_CrossNeibs];
   CrossEdgeMaptype = new int[N_CrossNeibs];

   for(i=0; i<N_CrossNeibs; i++)
    {
     ID=CrossNeibsRank[i];
     Max = -1;

     // if two neibs from same rank ID contain this edge, then the edge will be assigned to 
     // neib with max. global cell no
     for(j=0; j<N_Neibs; j++)
      if((ID==Neighbs[j]->GetSubDomainNo()) && (Max < Neighbs[j]->GetGlobalCellNo()) )
       {
        Max=Neighbs[j]->GetGlobalCellNo();
        Max_Index = j;
       }

      NeibCell = Neighbs[Max_Index];
      // find neibs local edge number
      k=0;
      while(this != NeibCell->GetEdge(k)) k++;

      CrossNeibsGlobalNo[i] = Max;
      CrossNeibsLocalEdgeNo[i] = k;

      (NeibCell->GetShapeDesc())->GetEdgeVertex(EdgeVertex);

      if(NeibCell->GetVertex(EdgeVertex[2*k])==Vert_a)
       { CrossEdgeMaptype[i] =  1; }
      else if(NeibCell->GetVertex(EdgeVertex[2*k+1])==Vert_a)
       { CrossEdgeMaptype[i] =  -1; }
      else
       {
        printf("NeibRank %d NeibCell No %d Error in finding cross edge maptype Edge %d\n", 
             NeibCell->GetSubDomainNo(), Max, k);
        MPI_Abort(MPI_COMM_WORLD, 0);
       }

    } //  for(i=0; i<N_CrossNeibs; i++)
   } // if
}

void TEdge::GetCrossEdgeNeibs(int &n_CrossNeibs, int *&crossNeibsRank,
                              int *&crossNeibsGlobalNo, int *&crossNeibsLocalEdgeNo,
                              int *&crossEdgeMaptype)
{
    n_CrossNeibs = N_CrossNeibs;
    crossNeibsRank = CrossNeibsRank;
    crossNeibsGlobalNo = CrossNeibsGlobalNo;
    crossNeibsLocalEdgeNo = CrossNeibsLocalEdgeNo;
    crossEdgeMaptype = CrossEdgeMaptype;
}

#endif 


// Destructor
TEdge::~TEdge()
{
  delete [] Neighbs;

#ifdef  _MPI 
  if(N_CrossNeibs)
   {
    delete [] CrossNeibsRank;

   if(SubDomainCrossEdge)
    {
     delete [] CrossNeibsGlobalNo;
     delete [] CrossNeibsLocalEdgeNo;
     delete [] CrossEdgeMaptype;
    }
   }

  if(SubDomainEdge)
    delete [] SubDomainNeibsRank;

#endif 

}



