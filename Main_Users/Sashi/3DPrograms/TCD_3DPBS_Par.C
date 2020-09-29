// =======================================================================
//
// Purpose:     Parallel main program with operator-splitting 
//                     for population balance system in 3D+1D
//
// Author:      Sashikumaar Ganesan
//
// Date:        15.07.2010 (start of implementation)
// ======================================================================= 

#ifdef _MPI
#  include "mpi.h"
#endif

#ifdef _OMP
#include <omp.h>
#endif

#include <Assemble3D.h>
#include <AuxParam3D.h>
#include <Bulk_3d4d.h>
#include <Collection.h>
#include <Convolution.h>
#include <Database.h>
#include <DirectSolver.h>
#include <DiscreteForm3D.h>
#include <Domain.h>
#include <FESpace1D.h>
#include <SquareStructure1D.h>
#include <SquareMatrix1D.h>
#include <FEFunction1D.h>
#include <FEDatabase2D.h>
#include <FEDatabase3D.h>
#include <FESpace3D.h>
#include <LinAlg.h>
#include <math.h>
#include <TCD3D.h>
#include <HexaAffin.h>
#include <TetraAffin.h>
#include <TetraIsoparametric.h>
#include <HexaIsoparametric.h>
#include <HexaTrilinear.h>
#include <ADISystem.h>
#include <ADISystem1D.h>
#include <NodalFunctional1D.h>
#include <LineAffin.h>

#include <Output3D.h>
#include <MainUtilities.h>
#include <TimeUtilities.h>


#include <string.h>
#include <sstream>
#include <MooNMD_Io.h>
#include <stdlib.h>
#include <malloc.h>

#ifdef _MPI
#include <MeshPartition.h>
#include <ParFECommunicator3D.h>
#include <MumpsSolver.h>
#include <ParVector3D.h>
#include <ParVectorNSE3D.h>
#include <Scalar_ParSolver.h>
#endif


#include <Urea_3d4d.h>
#include <LocalProjection.h>


#include <dirent.h> 
#include <unistd.h>

#include <fstream>

#include <sys/stat.h>
#include <sys/types.h>

double bound = 0;
// =======================================================================
// include current example
// =======================================================================
#include "../Examples/TCD_3D/PBS.h"
  
// #include "../Examples/TCD_3D/PBS_UnitCubeSmoothT.h"
 

void  GetVeloGradForPhysPts(TFEVectFunct3D *u, TFESpace3D *PBE_Spaces, double *velo, double *grad_velo,
                            int N_PhySpacePts)
{
 int i, j, k, l, m, N_Cells, N_Points, N_U;
 int N_LocalDOFs, *DOF, *BeginIndex, *GlobalNumbers; 
 int *U_BeginIndex, *U_GlobalNumbers, *IncidentArray;
 int U_N_LocalDOFs, *U_DOF;
 
 double s, *xi, *eta, *zeta, *U1, *U2, *U3, val1, val2, val3; 
 double X[MaxN_PointsForNodal3D], Y[MaxN_PointsForNodal3D];
 double Z[MaxN_PointsForNodal3D], AbsDetjk[MaxN_PointsForNodal3D];
 double FunctionalValuesX[MaxN_PointsForNodal3D];
 double FunctionalValuesY[MaxN_PointsForNodal3D];
 double FunctionalValuesZ[MaxN_PointsForNodal3D];
 double BasisValues[MaxN_BaseFunctions3D];
 double BasisValues_x[MaxN_BaseFunctions3D];
 double BasisValues_y[MaxN_BaseFunctions3D];
 double BasisValues_z[MaxN_BaseFunctions3D];
 
 TBaseCell *cell;
 TCollection *Coll; 
 TFESpace3D *U_FESpace;
 FE3D FEId, U_FEId;
 TFE3D *Element, *U_Element;
 TBaseFunct3D *bf, *U_bf;
 TNodalFunctional3D *nf;
 RefTrans3D F_K;  
 TRefTrans3D *rt;   
 
 
  // assume that all fespace2Ds and PBE_Spaces use same coll
  Coll = PBE_Spaces->GetCollection();
  N_Cells = Coll->GetN_Cells();  
  
  BeginIndex = PBE_Spaces->GetBeginIndex();
  GlobalNumbers = PBE_Spaces->GetGlobalNumbers();


  U_FESpace = u->GetFESpace3D();
  U_BeginIndex = U_FESpace->GetBeginIndex();
  U_GlobalNumbers = U_FESpace->GetGlobalNumbers();  
  N_U = U_FESpace->GetN_DegreesOfFreedom();
  U1 = u->GetValues();
  U2 = U1 + N_U;
  U3 = U2 + N_U; 
  
  IncidentArray = new int[N_PhySpacePts];
  
  memset(velo, 0, 3*N_PhySpacePts*SizeOfDouble);
  memset(grad_velo, 0, 9*N_PhySpacePts*SizeOfDouble);  
  memset(IncidentArray, 0, N_PhySpacePts*SizeOfInt); 
  
  
  for(i=0;i<N_Cells;i++)
  {
   cell = Coll->GetCell(i);
   FEId = PBE_Spaces->GetFE3D(i, cell);
   Element = TFEDatabase3D::GetFE3D(FEId);
   nf = Element->GetNodalFunctional3D();
   nf->GetPointsForAll(N_Points, xi, eta, zeta);
   N_LocalDOFs = Element->GetN_DOF();    
   DOF = GlobalNumbers + BeginIndex[i];
     
   F_K = Element->GetRefTransID(); 
    
   switch(F_K)
    {
      case HexaAffin:
        rt = TFEDatabase3D::GetRefTrans3D(HexaAffin);
        ((THexaAffin *)rt)->SetCell(cell);
        break;
      case HexaTrilinear:
        rt = TFEDatabase3D::GetRefTrans3D(HexaTrilinear);
        ((THexaTrilinear *)rt)->SetCell(cell);
        break;
      case HexaIsoparametric:
        rt = TFEDatabase3D::GetRefTrans3D(HexaIsoparametric);
        ((THexaIsoparametric *)rt)->SetCell(cell);
        break;
      case TetraAffin:
        rt = TFEDatabase3D::GetRefTrans3D(TetraAffin);
        ((TTetraAffin *)rt)->SetCell(cell);
        break;
      case TetraIsoparametric:
        rt = TFEDatabase3D::GetRefTrans3D(TetraIsoparametric);
        ((TTetraIsoparametric *)rt)->SetCell(cell);
        break;
     }   
     
    TFEDatabase3D::GetOrigFromRef(F_K, N_Points, xi, eta, zeta,
                                X, Y, Z, AbsDetjk);  
                                
    nf->GetAllFunctionals(Coll, cell, X, FunctionalValuesX);
    nf->GetAllFunctionals(Coll, cell, Y, FunctionalValuesY);
    nf->GetAllFunctionals(Coll, cell, Z, FunctionalValuesZ);
    
    U_FEId = U_FESpace->GetFE3D(i, cell);
    U_Element = TFEDatabase3D::GetFE3D(U_FEId);
    U_bf = U_Element->GetBaseFunct3D();
    U_N_LocalDOFs = U_Element->GetN_DOF();
    U_DOF = U_GlobalNumbers + U_BeginIndex[i];                               
    
     for(j=0;j<N_LocalDOFs;j++)
     {  
      k = DOF[j]; 
      U_bf->GetDerivatives(D000, FunctionalValuesX[j], FunctionalValuesY[j], FunctionalValuesZ[j], BasisValues);        
      U_bf->GetDerivatives(D100, FunctionalValuesX[j], FunctionalValuesY[j], FunctionalValuesZ[j], BasisValues_x);         
      U_bf->GetDerivatives(D010, FunctionalValuesX[j], FunctionalValuesY[j], FunctionalValuesZ[j], BasisValues_y);        
      U_bf->GetDerivatives(D001, FunctionalValuesX[j], FunctionalValuesY[j], FunctionalValuesZ[j], BasisValues_z);        
       
      for(l=0;l<U_N_LocalDOFs;l++)
      {
       m = U_DOF[l];
       val1 = U1[m];
       val2 = U2[m];       
       val3 = U3[m];
       
       velo[3*k    ] += BasisValues[l]*val1;
       velo[3*k + 1] += BasisValues[l]*val2;
       velo[3*k + 2] += BasisValues[l]*val3;
       
       grad_velo[9*k    ] += BasisValues_x[l]*val1;
       grad_velo[9*k + 1] += BasisValues_y[l]*val1;
       grad_velo[9*k + 2] += BasisValues_z[l]*val1;    
       
       grad_velo[9*k + 3] += BasisValues_x[l]*val2;
       grad_velo[9*k + 4] += BasisValues_y[l]*val2;
       grad_velo[9*k + 5] += BasisValues_z[l]*val2;         
       
       grad_velo[9*k + 6] += BasisValues_x[l]*val3;
       grad_velo[9*k + 7] += BasisValues_y[l]*val3;
       grad_velo[9*k + 8] += BasisValues_z[l]*val3; 
      } // for(l=0;l<U_N_LocalDOFs;l++)      
     IncidentArray[k]++;      
    } // for(j=0;j<U_N_LocalDOFs;j++)    
  } //for(i=0;i<N_Cells;i++)
  
  for(i=0;i<N_PhySpacePts;i++)
  {
    if(IncidentArray[i]==0)
    {
      cout << i<<  " IncidentArray[i] " << IncidentArray[i] <<endl;
#ifdef _MPI      
      MPI_Finalize();
#endif      
      exit(0);
    }
    
    for(j=0;j<3;j++)  
     velo[3*i+j] /=(double)IncidentArray[i];

    for(j=0;j<9;j++)  
     grad_velo[9*i+j] /=(double)IncidentArray[i];
  }    

  
//   printf("GetVeloGradForPhysPts \n" );
//   MPI_Finalize();
//   exit(0); 
  
  delete [] IncidentArray;
  
} // GetVeloGradForPhysPts


 
void SetDirDofFromNeu(double *rhs, double *sol, int N_neum_to_diri,
                      int *neum_to_diri, double *neum_to_diri_x,
                      double *neum_to_diri_y, double *neum_to_diri_z,
                      BoundValueFunct3D *BoundaryValue)
{
 
    double value;
    int i, index;
    
    //OutPut(N_neum_to_diri << endl);
    if (N_neum_to_diri == 0)
     return;

    // loop over dof to change
    for (i=0;i<N_neum_to_diri;i++)
     {
      index = neum_to_diri[i];
      
       // set boundary condition
       BoundaryValue(neum_to_diri_x[i], neum_to_diri_y[i],  neum_to_diri_z[i], rhs[index]);
       sol[index] = rhs[index];
    }
}


void GetPhyNodalFuctional(int N_U, int N_Intl_Levels, double *ValuesT, double *Sol)
                              
{
 int i, j, k;
 
 double *LocSol;
 
 for(i=0; i<N_Intl_Levels; i++)
  {
   LocSol = Sol+ (i*N_U);
   
   for(j=0; j<N_U; j++)   
    LocSol[j] = ValuesT[j*N_Intl_Levels + i ];
    
  } // for(i=0; i<N_PhySpace
  
} // GetNodalPtsValuesForIntl



void  Sol2IntlNodal(TFESpace1D *FESpace1D, double *Sol_AllL, int N_Levels,  
                    double *Sol_NodalPts, int N_XPoints)
{
  int i,j,k,l, ii;
  int N_Cells;
  int N_DOFs, N_LocalDOFs;
  int *BeginIndex, *GlobalNumbers;
  int N_Points, disp;
  int *DOF, *IndexArray, *NodalPtIndex, *LocNodalPtIndex;

  double *xi, *eta, *sol, *sol_Nodal, val;
  double Z[MaxN_PointsForNodal1D];
  double AbsDetjk[MaxN_PointsForNodal1D];
  double BasisValues[MaxN_PointsForNodal1D][MaxN_BaseFunctions1D];

  TBaseCell *cell;
  TCollection *Coll;
  FE1D FEId;
  TFE1D *Element;
  TFE1D *FE_Obj;
  TNodalFunctional1D *nf;
  TRefTrans1D *rt;
  TBaseFunct1D *bf;

  Coll = FESpace1D->GetCollection();
  N_Cells = Coll->GetN_Cells();
  BeginIndex = FESpace1D->GetBeginIndex();
  GlobalNumbers = FESpace1D->GetGlobalNumbers();
  N_DOFs = FESpace1D->GetN_DegreesOfFreedom();
  NodalPtIndex = FESpace1D->GetIntlPtIndexOfPts();

  IndexArray = new int[N_Levels];
  memset(IndexArray, 0, SizeOfInt*N_Levels);
  memset(Sol_NodalPts , 0, SizeOfDouble*N_XPoints*N_Levels);

  disp = 0;
  for(i=0; i<N_Cells; i++)
   {
    cell = Coll->GetCell(i);
    FEId = FESpace1D->GetFE1D(i, cell);
    Element = TFEDatabase2D::GetFE1D(FEId);
    FE_Obj = TFEDatabase2D::GetFE1D(FEId);
    bf = FE_Obj->GetBaseFunct1D();
    nf = Element->GetNodalFunctional1D();
    nf->GetPointsForAll(N_Points, xi, eta);
    N_LocalDOFs = Element->GetN_DOF();
    rt = TFEDatabase2D::GetRefTrans1D(LineAffin);
    ((TLineAffin *)rt)->SetCell(cell);
    ((TLineAffin *)rt)->GetOrigFromRef(N_Points, xi, Z, AbsDetjk);
    DOF = GlobalNumbers + BeginIndex[i];

     for(j=0;j<N_Points;j++)
       bf->GetDerivatives(D0, xi[j], BasisValues[j]);

     LocNodalPtIndex = NodalPtIndex+disp;
     
     for(ii=0;ii<N_XPoints;ii++)
      {
       sol = Sol_AllL + ii*N_DOFs;
       sol_Nodal =  Sol_NodalPts + ii*N_Levels; 

//        for(l=0;l<N_LocalDOFs;l++)
//         cout << " val " << sol[DOF[l]] << endl;

       for(j=0;j<N_Points;j++)
        {
         k = LocNodalPtIndex[j]; // find the right L level
         val = 0.;

          for(l=0;l<N_LocalDOFs;l++)
           val += sol[DOF[l]]*BasisValues[j][l];

         sol_Nodal[k] += val;
         if(ii==0)
          IndexArray[k]++;
        } // for(j=0;j<N_Points
      } // for(ii=0

     disp +=N_Points;
   } // for(i=0; i<N_Cells; i


   for(ii=0;ii<N_XPoints;ii++)
    {
     sol_Nodal =  Sol_NodalPts + ii*N_Levels; 
     for(i=0;i<N_Levels;i++)
      {
       if(ii==0)
        if(IndexArray[i] == 0)
         {
          cout << "Error in Sol2IntlNodal : "<< IndexArray[i] << endl;
          exit(0);
         }

       sol_Nodal[i] /= (double)IndexArray[i];
      }
    }
}

void GetSolFromNodalPtVales(TFEFunction3D *ScalarFunction, double *Values, int N_PhySpacePts, 
                            TFESpace3D *PBE_Space)
{  
 int i, j, k, l, m, N_Dof, N_Cells, N_Points, *DOF, N_LocalDOFs, *PBE_DOF, PBE_N_LocalDOFs;
 int *BeginIndex, *GlobalNumbers, *PBE_BeginIndex, *PBE_GlobalNumbers, *IncidentArray; 
 
 double s, *xi, *eta, *zeta, *sol, PBE_LocValues[MaxN_PointsForNodal3D];
 double X[MaxN_PointsForNodal3D], Y[MaxN_PointsForNodal3D];
 double Z[MaxN_PointsForNodal3D];
 double AbsDetjk[MaxN_PointsForNodal3D];
 double BasisValues[MaxN_BaseFunctions3D];
 double FunctionalValues[MaxN_PointsForNodal3D];
 
 TBaseCell *cell;
 TCollection *Coll;  
 TFESpace3D *FE_Space;
 TFE3D *Element, *PBE_Element;
 TBaseFunct3D *bf, *PBE_bf;
 TNodalFunctional3D *nf;
 RefTrans3D F_K;  
 TRefTrans3D *rt;  
 FE3D FEId, PBE_FEId; 

  PBE_BeginIndex = PBE_Space->GetBeginIndex();
  PBE_GlobalNumbers = PBE_Space->GetGlobalNumbers(); 
  
  if(PBE_Space->GetN_DegreesOfFreedom()!=N_PhySpacePts )
   {
    // Values should be the nodal functionals in this currtent implementation
    printf("N_Dof!=N_PhySpacePts in GetSolFromNodalPtVales\n" );
#ifdef _MPI      
      MPI_Finalize();
#endif  
    exit(0); 
   }
   
  sol = ScalarFunction->GetValues();
  FE_Space = ScalarFunction->GetFESpace3D(); 
  N_Dof = FE_Space->GetN_DegreesOfFreedom();    
  BeginIndex = FE_Space->GetBeginIndex();
  GlobalNumbers = FE_Space->GetGlobalNumbers();  
  Coll = FE_Space->GetCollection();// assume that all fespaces use same coll
  N_Cells = Coll->GetN_Cells();  
  
  IncidentArray = new int[N_Dof];
  memset(IncidentArray, 0, N_Dof*SizeOfInt);   
  memset(sol, 0, N_Dof*SizeOfDouble);
  
  for(i=0;i<N_Cells;i++)
  {
   cell = Coll->GetCell(i);
   FEId = FE_Space->GetFE3D(i, cell);
   Element = TFEDatabase3D::GetFE3D(FEId);
   nf = Element->GetNodalFunctional3D();
   nf->GetPointsForAll(N_Points, xi, eta, zeta);
   N_LocalDOFs = Element->GetN_DOF();    
   DOF = GlobalNumbers + BeginIndex[i];   

   F_K = Element->GetRefTransID(); 
    
   switch(F_K)
    {
      case HexaAffin:
        rt = TFEDatabase3D::GetRefTrans3D(HexaAffin);
        ((THexaAffin *)rt)->SetCell(cell);
        break;
      case HexaTrilinear:
        rt = TFEDatabase3D::GetRefTrans3D(HexaTrilinear);
        ((THexaTrilinear *)rt)->SetCell(cell);
        break;
      case HexaIsoparametric:
        rt = TFEDatabase3D::GetRefTrans3D(HexaIsoparametric);
        ((THexaIsoparametric *)rt)->SetCell(cell);
        break;
      case TetraAffin:
        rt = TFEDatabase3D::GetRefTrans3D(TetraAffin);
        ((TTetraAffin *)rt)->SetCell(cell);
        break;
      case TetraIsoparametric:
        rt = TFEDatabase3D::GetRefTrans3D(TetraIsoparametric);
        ((TTetraIsoparametric *)rt)->SetCell(cell);
        break;
     }
     
    TFEDatabase3D::GetOrigFromRef(F_K, N_Points, xi, eta, zeta,
                                X, Y, Z, AbsDetjk);    


    PBE_FEId = PBE_Space->GetFE3D(i, cell);
    PBE_Element = TFEDatabase3D::GetFE3D(PBE_FEId);
    PBE_bf = PBE_Element->GetBaseFunct3D();
    PBE_N_LocalDOFs = PBE_Element->GetN_DOF();
    PBE_DOF = PBE_GlobalNumbers + PBE_BeginIndex[i];

    memset(PBE_LocValues, 0, MaxN_PointsForNodal3D*SizeOfDouble);

    for(j=0;j<N_Points;j++)
     {
      // first find T value at this point
      PBE_bf->GetDerivatives(D000, xi[j], eta[j], zeta[j], BasisValues);      
      for(l=0;l<PBE_N_LocalDOFs;l++)
      {
        m = PBE_DOF[l];
        PBE_LocValues[j] += BasisValues[l]*Values[m];
      }    
     } // for(j=0;j<N_Points;j++)    

    nf->GetAllFunctionals(Coll, cell,  PBE_LocValues, FunctionalValues); 
 
    for(j=0;j<N_LocalDOFs;j++)
    {
      k = DOF[j];
      sol[k] += FunctionalValues[j];
      IncidentArray[k]++;
    }   

  }//  for(i=0;i<N_Cells;i++
  
  for(i=0;i<N_Dof;i++)
  {
    if(IncidentArray[i]==0)
    {
      cout << i<<  " IncidentArray[i] " << IncidentArray[i] <<endl;
#ifdef _MPI      
      MPI_Finalize();
#endif  
      exit(0);
    }
    
    sol[i] /=(double)IncidentArray[i];
  }  
 
  delete [] IncidentArray;
  
//   printf("GetSolFromNodalPtVales not yet \n" );
//   MPI_Finalize();
//   exit(0); 
 
} // GetSolFromNodalPtVales


void  IntlNodal2Sol(TFESpace1D *FESpace1D, double *Sol_QuadIntl, int N_PhySpacePts, 
                    int N_Levels, double *Sol_AllL)
{
  int i,j,k,l;
  int N_Cells;
  int N_DOFs, N_LocalDOFs, *IncidentArray;
  int *BeginIndex, *GlobalNumbers;
  int N_Points, ii;
  int *DOF, *Index, *IndexOfNodalPts, disp;

  double *xi, *eta, *Values_Level, *sol;
  double X[MaxN_PointsForNodal1D], Y[MaxN_PointsForNodal1D];
  double AbsDetjk[MaxN_PointsForNodal1D];
  double *PointValues, *PtVal;
  double FunctionalValues[MaxN_PointsForNodal1D];
  double FctVal[4];

  TBaseCell *cell;
  TCollection *Coll;
  FE1D FEId;
  TFE1D *Element;
  TNodalFunctional1D *nf;

  Coll = FESpace1D->GetCollection();
  N_Cells = Coll->GetN_Cells();
  IndexOfNodalPts = FESpace1D->GetIntlPtIndexOfPts();

  BeginIndex = FESpace1D->GetBeginIndex();
  GlobalNumbers = FESpace1D->GetGlobalNumbers();
  N_DOFs = FESpace1D->GetN_DegreesOfFreedom();

  PointValues = new double [MaxN_PointsForNodal1D*N_PhySpacePts];
  IncidentArray = new int [N_DOFs];
  memset(Sol_AllL, 0, SizeOfDouble*N_DOFs*N_PhySpacePts);
  memset(IncidentArray, 0, SizeOfInt*N_DOFs);

  disp = 0; 

  for(i=0;i<N_Cells;i++)
   {
    cell = Coll->GetCell(i);
    FEId = FESpace1D->GetFE1D(i, cell);
    Element = TFEDatabase2D::GetFE1D(FEId);
    nf = Element->GetNodalFunctional1D();
    nf->GetPointsForAll(N_Points, xi, eta);
    N_LocalDOFs = Element->GetN_DOF();
    DOF = GlobalNumbers+BeginIndex[i];

    Index = IndexOfNodalPts+disp;

    // get the Sol_QuadIntl values for all pts this Internal cell
    // corresponding to all XLocPoints
    for(ii=0;ii<N_PhySpacePts;ii++)
     {
      Values_Level  =  Sol_QuadIntl + ii*N_Levels;

      for(j=0;j<N_Points;j++)
      {
       k = Index[j]; // corresponding L level
       PointValues[ii*N_Points + j] = Values_Level[k];
      }
     }//  for(ii=0;ii<N_PhySpacePts;
     
    for(ii=0;ii<N_PhySpacePts;ii++)
     {
      PtVal = PointValues + ii*N_Points;
      nf->GetAllFunctionals(PtVal, FunctionalValues);

      sol = Sol_AllL + ii*N_DOFs;

      for(j=0;j<N_LocalDOFs;j++)
       {
        sol[DOF[j]] += FunctionalValues[j];

        if(ii==0)
         IncidentArray[DOF[j]] ++;
       }
     } //  for(ii=0;ii<N_PhySpacePts;ii++)

    disp +=N_Points;     
     
   }// for(i=0;i<N_Cells;i++)
  
  
  for(ii=0;ii<N_PhySpacePts;ii++)
    {      
     sol = Sol_AllL + ii*N_DOFs;      
     
     for(i=0;i<N_DOFs;i++)
      {

       if(ii==0)
        if(IncidentArray[i] == 0)
         {
          cout << "Error in IntlNodal2Sol : "<< IncidentArray[i] << endl;
          exit(0);
         }
         
       sol[i] /= (double)IncidentArray[i];
      } // for(i=0;i<N_DOFs;i
    } // for(ii=0;ii<N_PhySpacePts;i

  delete [] PointValues;
  delete [] IncidentArray;  

//   printf("IntlNodal2Sol not yet \n" );
//   MPI_Finalize();
//   exit(0); 
//   
  
} // IntlNodal2Sol


void GetNodalPtsValuesForIntl(int N_Levels, TFEFunction3D **ScalarFunctions, double *Sol,
                              int N_U, double *ValuesT, int N_IntlPts)
{
  int i, ii, j, k, l, m, N_Cells, N_Points, N_LocalDOFs;
  int *BeginIndex, *GlobalNumbers, *DOF;
  int PolynomialDegree, ApproxOrder, *IntlPtIndexOfPts, *PtIndexLoc, disp;
  int *IncidentArray, *Incident;

  double *xi, *eta, *zeta, *NodalValues, s, *Values_Level;
  double BasisValues[MaxN_BaseFunctions3D], maxval=0.;

  TFESpace3D *PBE_Spaces;
  TBaseCell *cell;
  TCollection *Coll;
  FE3D FEId_PBE;
  TFE3D *Element_PBE;
  TBaseFunct3D *bf_PBE;
  TNodalFunctional3D *nf_PBE;
  
  // assume all  ScalarFunctions use same fespace3D
  PBE_Spaces = ScalarFunctions[0]->GetFESpace3D();
  BeginIndex = PBE_Spaces->GetBeginIndex();
  GlobalNumbers = PBE_Spaces->GetGlobalNumbers();

  IncidentArray = new int[N_Levels*N_IntlPts];

  memset(ValuesT, 0, N_Levels*N_IntlPts*SizeOfDouble);
  memset(IncidentArray, 0, N_Levels*N_IntlPts*SizeOfInt);
  
  // assume that both fespace2D and PBE_Spaces use same coll
  Coll = PBE_Spaces->GetCollection();
  N_Cells = Coll->GetN_Cells();

  for(i=0;i<N_Cells;i++)
   {
    cell = Coll->GetCell(i);
    FEId_PBE = PBE_Spaces->GetFE3D(i, cell);
    Element_PBE = TFEDatabase3D::GetFE3D(FEId_PBE);
    nf_PBE = Element_PBE->GetNodalFunctional3D();
    nf_PBE->GetPointsForAll(N_Points, xi, eta, zeta);
    bf_PBE = Element_PBE->GetBaseFunct3D();
    N_LocalDOFs = Element_PBE->GetN_DOF();
    
    DOF = GlobalNumbers + BeginIndex[i];
    PtIndexLoc = IntlPtIndexOfPts + disp;

    for(j=0;j<N_Points;j++)
     {
      bf_PBE->GetDerivatives(D000, xi[j], eta[j], zeta[j], BasisValues);      

//       k = PtIndexLoc[j];
// 
//       Incident = IncidentArray + k*N_Levels;
//       Values_Level  =  ValuesT + k*N_Levels;

      
     }//  for(j=0;j<N_Points;j++) 
    
   } //  for(i=0;i<N_Cells;i++)
  
  
  delete []  IncidentArray;
  
  printf("GetNodalPtsValuesForIntl not yet \n" );
#ifdef _MPI      
      MPI_Finalize();
#endif  
  exit(0); 
  
  
  
} // GetNodalPtsValuesForIntl
                               
                               
void GetNodalPtsValuesForIntl(int N_U, int N_Intl_Levels, double *Sol, double *ValuesT)
                              
{
 int i, j, k;
 
 double *LocSol;
 
 for(i=0; i<N_Intl_Levels; i++)
  {
   LocSol = Sol+ (i*N_U);
   
   for(j=0; j<N_U; j++)   
    ValuesT[j*N_Intl_Levels + i ] = LocSol[j];
    
  } // for(i=0; i<N_PhySpace
  
} // GetNodalPtsValuesForIntl


void GetNodalPtsValuesForIntl(TFEFunction3D **ScalarFunctions, TFESpace3D *PBE_Spaces,
                              double *T_Values, double *C_Values,
                              double *C_Sat_Values, int N_PhySpacePts)
{
 int i, j, k, l, m, N_Cells, N_Points;
 int *BeginIndex, *GlobalNumbers; 
 int *T_BeginIndex, *T_GlobalNumbers, *C_BeginIndex, *C_GlobalNumbers;
 int *IncidentArray, N_LocalDOFs, *DOF, T_N_LocalDOFs, *T_DOF, C_N_LocalDOFs, *C_DOF;
 
 double *T_NodalValues, *C_NodalValues, T_LocValues[MaxN_PointsForNodal3D], C_LocValues[MaxN_PointsForNodal3D];
 double s, *xi, *eta, *zeta;
 double X[MaxN_PointsForNodal3D], Y[MaxN_PointsForNodal3D];
 double Z[MaxN_PointsForNodal3D];
 double AbsDetjk[MaxN_PointsForNodal3D];
 double BasisValues[MaxN_BaseFunctions3D];
 double FunctionalValuesT[MaxN_PointsForNodal3D];
 double FunctionalValuesC[MaxN_PointsForNodal3D];
 double T_ref = TDatabase::ParamDB->REACTOR_P23;
 double C_ref = TDatabase::ParamDB->REACTOR_P25;  
  
 TBaseCell *cell;
 TCollection *Coll;  
 TFEFunction3D *Heat, *Concentration;
 TFESpace3D *T_FESpace, *C_FESpace;
 FE3D FEId, T_FEId, C_FEId;
 TFE3D *Element, *T_Element, *C_Element;
 TBaseFunct3D *bf, *T_bf, *C_bf;
 TNodalFunctional3D *nf;
 RefTrans3D F_K;  
 TRefTrans3D *rt;  

 
  BeginIndex = PBE_Spaces->GetBeginIndex();
  GlobalNumbers = PBE_Spaces->GetGlobalNumbers();
 
  Heat = ScalarFunctions[0];
  Concentration   = ScalarFunctions[1];
  
  T_NodalValues =  Heat->GetValues();
  T_FESpace = Heat->GetFESpace3D();
  T_BeginIndex = T_FESpace->GetBeginIndex();
  T_GlobalNumbers = T_FESpace->GetGlobalNumbers();

  C_NodalValues =  Concentration->GetValues();
  C_FESpace = Concentration->GetFESpace3D();
  C_BeginIndex = C_FESpace->GetBeginIndex();
  C_GlobalNumbers = C_FESpace->GetGlobalNumbers();  

  IncidentArray = new int[N_PhySpacePts];
  
  // assume that all fespace2Ds and PBE_Spaces use same coll
  Coll = T_FESpace->GetCollection();
  N_Cells = Coll->GetN_Cells();
  
  memset(T_Values, 0, N_PhySpacePts*SizeOfDouble);
  memset(C_Values, 0, N_PhySpacePts*SizeOfDouble);
  memset(C_Sat_Values, 0, N_PhySpacePts*SizeOfDouble);
  memset(IncidentArray, 0, N_PhySpacePts*SizeOfInt);  
  
  for(i=0;i<N_Cells;i++)
  {
   cell = Coll->GetCell(i);
   FEId = PBE_Spaces->GetFE3D(i, cell);
   Element = TFEDatabase3D::GetFE3D(FEId);
   nf = Element->GetNodalFunctional3D();
   nf->GetPointsForAll(N_Points, xi, eta, zeta);
   N_LocalDOFs = Element->GetN_DOF();    
   DOF = GlobalNumbers + BeginIndex[i];
     
   F_K = Element->GetRefTransID(); 
    
   switch(F_K)
    {
      case HexaAffin:
        rt = TFEDatabase3D::GetRefTrans3D(HexaAffin);
        ((THexaAffin *)rt)->SetCell(cell);
        break;
      case HexaTrilinear:
        rt = TFEDatabase3D::GetRefTrans3D(HexaTrilinear);
        ((THexaTrilinear *)rt)->SetCell(cell);
        break;
      case HexaIsoparametric:
        rt = TFEDatabase3D::GetRefTrans3D(HexaIsoparametric);
        ((THexaIsoparametric *)rt)->SetCell(cell);
        break;
      case TetraAffin:
        rt = TFEDatabase3D::GetRefTrans3D(TetraAffin);
        ((TTetraAffin *)rt)->SetCell(cell);
        break;
      case TetraIsoparametric:
        rt = TFEDatabase3D::GetRefTrans3D(TetraIsoparametric);
        ((TTetraIsoparametric *)rt)->SetCell(cell);
        break;
     }
     
    TFEDatabase3D::GetOrigFromRef(F_K, N_Points, xi, eta, zeta,
                                X, Y, Z, AbsDetjk);    

    T_FEId = T_FESpace->GetFE3D(i, cell);
    T_Element = TFEDatabase3D::GetFE3D(T_FEId);
    T_bf = T_Element->GetBaseFunct3D();
    T_N_LocalDOFs = T_Element->GetN_DOF();
    T_DOF = T_GlobalNumbers + T_BeginIndex[i];

    C_FEId = C_FESpace->GetFE3D(i, cell);
    C_Element = TFEDatabase3D::GetFE3D(C_FEId);
    C_bf = C_Element->GetBaseFunct3D();
    C_N_LocalDOFs = C_Element->GetN_DOF();
    C_DOF = C_GlobalNumbers + C_BeginIndex[i];
    
    memset(T_LocValues, 0, MaxN_PointsForNodal3D*SizeOfDouble);
    memset(C_LocValues, 0, MaxN_PointsForNodal3D*SizeOfDouble);   
    
    for(j=0;j<N_Points;j++)
     {
      // first find T value at this point
      T_bf->GetDerivatives(D000, xi[j], eta[j], zeta[j], BasisValues);      
      for(l=0;l<T_N_LocalDOFs;l++)
      {
        m = T_DOF[l];
        T_LocValues[j] += BasisValues[l]*T_NodalValues[m];
      }
      
      // find C value at this point 
      C_bf->GetDerivatives(D000, xi[j], eta[j], zeta[j], BasisValues);      
      for(l=0;l<C_N_LocalDOFs;l++)
       {
        m = C_DOF[l];
        C_LocValues[j] += BasisValues[l]*C_NodalValues[m];
       }
     } // for(j=0;j<N_Points;j++)
    
    nf->GetAllFunctionals(Coll, cell,  T_LocValues, FunctionalValuesT);  
    nf->GetAllFunctionals(Coll, cell,  C_LocValues, FunctionalValuesC);    

    DOF = GlobalNumbers+BeginIndex[i];

    for(j=0;j<N_LocalDOFs;j++)
    {
      k = DOF[j];
      T_Values[k] += FunctionalValuesT[j];
      C_Values[k] += FunctionalValuesC[j];
      IncidentArray[k]++;
    }   
  } //  for(i=0;i<N_Cells;i++)
  
  for(i=0;i<N_PhySpacePts;i++)
  {
    if(IncidentArray[i]==0)
    {
      cout << i<<  " IncidentArray[i] " << IncidentArray[i] <<endl;
#ifdef _MPI      
      MPI_Finalize();
#endif  
      exit(0);
    }
    
    T_Values[i] /=(double)IncidentArray[i];
    C_Values[i] /=(double)IncidentArray[i];

    // now compute C_Sat from T
    C_Sat_Values[i] = (1.3045*(T_ref*T_Values[i] - 273.15) + 35.3642)/(C_ref*TDatabase::ParamDB->UREA_m_mol);
  }  
     
//   printf("GetNodalPtsValuesForIntl \n" );
//   MPI_Finalize();
//   exit(0); 
  
  delete [] IncidentArray;
} // GetNodalPtsValuesForIntl


void GetPhysicalSpaceQuadPts(int &N_PhySpacePts, TFESpace3D *FESpace3D, double *&IntlX, double *&IntlY, double *&IntlZ)
{
  int i,j,k,l, m;
  int N_Cells, PolynomialDegree, N_Points;

  double X[MaxN_QuadPoints_3D], Y[MaxN_QuadPoints_3D], Z[MaxN_QuadPoints_3D], AbsDetjk[MaxN_QuadPoints_3D];
  double *xi, *eta, *zeta, *weights;

  TBaseCell *cell;
  TCollection *Coll;
  FE3D FEId;
  TFE3D *Element;
  BF3DRefElements RefElement;
  QuadFormula3D QuadFormula;
  TQuadFormula3D *qf3;
  TRefTrans3D *F_K;

  Coll = FESpace3D->GetCollection();
  N_Cells = Coll->GetN_Cells();

   m = 0;
  for(i=0; i<N_Cells; i++)
  {
    cell = Coll->GetCell(i);
    FEId = FESpace3D->GetFE3D(i, cell);
    Element = TFEDatabase3D::GetFE3D(FEId);
    RefElement = TFEDatabase3D::GetRefElementFromFE3D(FEId);
    PolynomialDegree = TFEDatabase3D::GetPolynomialDegreeFromFE3D(FEId);

     switch(RefElement)
    {
      case BFUnitHexahedron:
        QuadFormula = TFEDatabase3D::GetQFHexaFromDegree(3*PolynomialDegree);
        F_K = TFEDatabase3D::GetRefTrans3D(HexaAffin);
        ((THexaAffin *)F_K)->SetCell(cell);
        break;

      case BFUnitTetrahedron:
        QuadFormula = TFEDatabase3D::GetQFTetraFromDegree(3*PolynomialDegree-1);
        F_K = TFEDatabase3D::GetRefTrans3D(TetraAffin);
        ((TTetraAffin *)F_K)->SetCell(cell);
        break;
    }                                             // endswitch

//     cout << "QuadFormula: " << QuadFormula << endl;
    qf3 = TFEDatabase3D::GetQuadFormula3D(QuadFormula);
    qf3->GetFormulaData(N_Points, weights, xi, eta, zeta);

    if(i==0)
    {
      N_PhySpacePts = N_Points*N_Cells;

      IntlX = new double[N_PhySpacePts];
      IntlY = new double[N_PhySpacePts];
      IntlZ = new double[N_PhySpacePts];      
    }

    switch(RefElement)
    {
      case BFUnitHexahedron:
        ((THexaAffin *)F_K)->GetOrigFromRef(N_Points, xi, eta, zeta, X, Y, Z, AbsDetjk);
        break;

      case BFUnitTetrahedron:
        ((TTetraAffin *)F_K)->GetOrigFromRef(N_Points, xi, eta, zeta, X, Y, Z, AbsDetjk);
        break;
    }     

    for(j=0; j<N_Points; j++)
    {
      IntlX[m] = X[j];
      IntlY[m] = Y[j];
      IntlZ[m] = Z[j];
      m++;
    }
       
  } //   for(i=0; i<N_Cells; i++)
  
//      cout<< N_PhySpacePts << " N_PhySpacePts " << m <<endl;
//     for(i=0; i<N_PhySpacePts; i++)
//      cout<< i << " IntlX " << IntlX[i] << " IntlY " << IntlY[i]  << " IntlZ " << IntlZ[i]  <<  endl;
//    exit(0);
//   cout << " GetPhysicalSpaceQuadPts " << endl;
//   exit(0);
} // GetPhysicalSpaceQuadPts


void GetInternalNodalPts(TFESpace1D *FeSpace_Intl, int &N_Intl_Levels, double *&IntlPosL)
{
  int i, j, k, l, m, r, N_Cells, N_RootPts, *RootPtIndex;
  int N_DOFs, N_LocalDOFs, N_Points;
  int N_AllLocalPoints;

  double L, L0, *xi, *eta, *L_loc, *L_loc_origOrder;
  double X[MaxN_PointsForNodal1D], Y[MaxN_PointsForNodal1D];
  double AbsDetjk[MaxN_PointsForNodal1D];

  TBaseCell *cell;
  TCollection *Coll;
  FE1D FEId;
  TFE1D *Element;
  TFE1D *FE_Obj;
  TNodalFunctional1D *nf;
  TRefTrans1D *rt;
  TBaseFunct1D *bf;

  Coll = FeSpace_Intl->GetCollection();
  N_Cells = Coll->GetN_Cells();

  N_AllLocalPoints = 0;
  for(i=0; i<N_Cells; i++)
  {
    cell = Coll->GetCell(i);
    FEId = FeSpace_Intl->GetFE1D(i, cell);
    Element = TFEDatabase2D::GetFE1D(FEId);
    FE_Obj = TFEDatabase2D::GetFE1D(FEId);
    bf = FE_Obj->GetBaseFunct1D();
    nf = Element->GetNodalFunctional1D();
    nf->GetPointsForAll(N_Points, xi, eta);
    N_LocalDOFs = Element->GetN_DOF();
    N_AllLocalPoints +=N_Points;
  }

  L_loc = new double [N_AllLocalPoints];
  L_loc_origOrder = new double [N_AllLocalPoints];
  N_AllLocalPoints = 0;
  for(i=0; i<N_Cells; i++)
  {
    cell = Coll->GetCell(i);
    FEId = FeSpace_Intl->GetFE1D(i, cell);
    Element = TFEDatabase2D::GetFE1D(FEId);
    FE_Obj = TFEDatabase2D::GetFE1D(FEId);
    bf = FE_Obj->GetBaseFunct1D();
    nf = Element->GetNodalFunctional1D();
    nf->GetPointsForAll(N_Points, xi, eta);
    N_LocalDOFs = Element->GetN_DOF();

    rt = TFEDatabase2D::GetRefTrans1D(LineAffin);
    ((TLineAffin *)rt)->SetCell(cell);

    ((TLineAffin *)rt)->GetOrigFromRef(N_Points, xi, X, Y, AbsDetjk);

    for(j=0; j<N_Points; j++)
    {
      L_loc[N_AllLocalPoints] = X[j];
      N_AllLocalPoints++;
    }
  } // for(i=0; i<N_Cells; i++)

  memcpy(L_loc_origOrder, L_loc,  N_AllLocalPoints*SizeOfDouble);

  for(i=0; i<N_AllLocalPoints-1; i++)
   for(j=i; j<N_AllLocalPoints; j++)
    if(L_loc[i]> L_loc[j])
     {
      L= L_loc[i];
      L_loc[i]= L_loc[j];
      L_loc[j]= L;
     }

  L  = L_loc[0];
  N_RootPts = 1;

  for(i=1; i<N_AllLocalPoints; i++)
   {
    if( fabs(L_loc[i]-L)>1e-5 )
     {
      N_RootPts++;
      L = L_loc[i];
     }
   }

  IntlPosL= new double[N_AllLocalPoints];
  IntlPosL[0] = L_loc[0];
  N_RootPts = 1;
  L  = L_loc[0];

  for(i=1; i<N_AllLocalPoints; i++)
   {
    if( fabs(L_loc[i]-L)>1e-5 )
     {
      IntlPosL[N_RootPts] = L_loc[i];
      N_RootPts++;
      L = L_loc[i];
     }
   }

  delete [] L_loc;

  RootPtIndex = new int[N_AllLocalPoints];

  // find the index for the local points in the root points
  for(i=0; i<N_AllLocalPoints; i++)
   {
    L = L_loc_origOrder[i];
    l=0;
    r=N_RootPts;

    m = N_RootPts/2;
    L0 = IntlPosL[m];

    while(fabs(L-L0) > 1.e-8 )
     {
      if(L < L0)  //poin lies in left side
       {
        r = m;
       }
      else
       {
        l=m;
       }

      m= (l+r)/2;
      L0 = IntlPosL[m];
     } //  while ( 

    RootPtIndex[i] = m;
   }

  FeSpace_Intl->SetIntlPtIndexOfPts(RootPtIndex);
  FeSpace_Intl->SetN_RootNodalPts(N_RootPts);

  N_Intl_Levels = N_RootPts;

//  cout << N_AllLocalPoints << " N_RootPts  "  << N_RootPts << endl;
//   for(i=0; i<N_RootPts; i++)
//   cout << i << " L: "  << IntlPosL[i] << endl;
// exit(0);
}

#ifndef __PBSConstT__   
void B_nuc(TFEFunction3D *f, TFEFunction3D *Heat, TFEFunction3D *Concentration)
{
  int i, j, k, l, m, N_Cells, N_Points, C_N_LocalDOFs, T_N_LocalDOFs;
  int N_DOF_PBE, *DOF_PBE, *BeginIndex_PBE, *GlobalNumbers_PBE, N_LocalDOFs_PBE;
  int *C_BeginIndex, *C_GlobalNumbers, *C_DOF;
  int *T_BeginIndex, *T_GlobalNumbers, *T_DOF;
  int PolynomialDegree, ApproxOrder, *IntlPtIndexOfPts, *PtIndexLoc, disp;
  int *IncidentArray;

  double *xi, *eta, *zeta, *C_NodalValues, s, *T_NodalValues, *Values_PBE;
  double BasisValues[MaxN_BaseFunctions3D], T_val, C_val, C_sat, val;
  double G, b_nuc[MaxN_PointsForNodal3D];
  double FunctionalValues[MaxN_PointsForNodal3D], Inn[2], Out[2];
  
  TFESpace3D *C_fespace3D, *T_fespace3D, *fespace3D_PBE;
  TBaseCell *cell;
  TCollection *Coll;
  FE3D C_FEId, T_FEId, FEId_PBE;
  TFE3D *C_Element, *T_Element, *Element_PBE;
  TBaseFunct3D *C_bf, *T_bf;
  TNodalFunctional3D *nf_PBE;
    
  T_NodalValues =  Heat->GetValues();
  T_fespace3D = Heat->GetFESpace3D();
  T_BeginIndex = T_fespace3D->GetBeginIndex();
  T_GlobalNumbers = T_fespace3D->GetGlobalNumbers();

  C_NodalValues =  Concentration->GetValues();
  C_fespace3D = Concentration->GetFESpace3D();
  C_BeginIndex = C_fespace3D->GetBeginIndex();
  C_GlobalNumbers = C_fespace3D->GetGlobalNumbers();  

  Values_PBE = f->GetValues();
  fespace3D_PBE = f->GetFESpace3D();
  N_DOF_PBE = fespace3D_PBE->GetN_DegreesOfFreedom(); 
  BeginIndex_PBE = fespace3D_PBE->GetBeginIndex();
  GlobalNumbers_PBE = fespace3D_PBE->GetGlobalNumbers();  
  
  IncidentArray = new int[N_DOF_PBE];   
  
  memset(Values_PBE, 0, N_DOF_PBE*SizeOfDouble);
  memset(IncidentArray, 0, N_DOF_PBE*SizeOfInt);

  // assume that all fespace3Ds and PBE_Spaces use same coll
  Coll = fespace3D_PBE->GetCollection();
  N_Cells = Coll->GetN_Cells();
 
  for(i=0;i<N_Cells;i++)
   {
    cell = Coll->GetCell(i);
    FEId_PBE = fespace3D_PBE->GetFE3D(i, cell);
    Element_PBE = TFEDatabase3D::GetFE3D(FEId_PBE);
    nf_PBE = Element_PBE->GetNodalFunctional3D();
    nf_PBE->GetPointsForAll(N_Points, xi, eta, zeta);
    N_LocalDOFs_PBE = Element_PBE->GetN_DOF();
    DOF_PBE = GlobalNumbers_PBE+BeginIndex_PBE[i];    
    
    C_FEId = C_fespace3D->GetFE3D(i, cell);
    C_Element = TFEDatabase3D::GetFE3D(C_FEId);
    C_bf = C_Element->GetBaseFunct3D();
    C_N_LocalDOFs = C_Element->GetN_DOF();
    C_DOF = C_GlobalNumbers + C_BeginIndex[i];

    T_FEId = T_fespace3D->GetFE3D(i, cell);
    T_Element = TFEDatabase3D::GetFE3D(T_FEId);
    T_bf = T_Element->GetBaseFunct3D();
    T_N_LocalDOFs = T_Element->GetN_DOF();
    T_DOF = T_GlobalNumbers + T_BeginIndex[i];

    for(j=0;j<N_Points;j++)
     {
      // C value at this point
      C_val = 0.;
      C_bf->GetDerivatives(D000, xi[j], eta[j], zeta[j], BasisValues);
      for(l=0;l<C_N_LocalDOFs;l++)
       {
        m = C_DOF[l];
        s = C_NodalValues[m];
        C_val += BasisValues[l]*s;
       }

      //find T value at this point,
      T_val = 0.;
      T_bf->GetDerivatives(D000, xi[j], eta[j], zeta[j], BasisValues);
      for(l=0;l<T_N_LocalDOFs;l++)
       {
        m = T_DOF[l];
        s = T_NodalValues[m];
        T_val += BasisValues[l]*s;
       }

      // see example file
      Inn[0] = C_val;
      Inn[1] = T_val;   
      
      Get_GrowthAndB_Nuc(2, Inn, Out);
      b_nuc[j] = Out[1]; 
     
     } // for(j=0;j<N_Points;j++)    
    
    nf_PBE->GetAllFunctionals(Coll, (TGridCell *)cell, b_nuc, FunctionalValues);

    for(j=0;j<N_LocalDOFs_PBE;j++)
     {
      k = DOF_PBE[j];
      Values_PBE[k] += FunctionalValues[j];
      IncidentArray[k]++;
     }        
   }// for(i=0;i<N_Cells;i++)
  
  
  for(i=0;i<N_DOF_PBE;i++)
   {
    if(IncidentArray[i]==0)
    {
     cout << i<<  " IncidentArray[i] " << IncidentArray[i] <<endl;
#ifdef _MPI
     MPI_Finalize();
#endif
     exit(0);
    }
    Values_PBE[i] /=(double)IncidentArray[i];
   }     
 
 
 delete [] IncidentArray;
 
//     printf("B_nuc temp and conc \n");   
// #ifdef _MPI      
//           MPI_Finalize(); 
// #endif    
//           exit(0); 
} // void B_nuc(TFEFunction3D *f
#endif

void SetLMinLMaxBoundValue(int N_Intl_Levels, TFEFunction3D **PbeFuncts, int N_U,  TFEFunction3D **ScalarFunctions,
                           double *RhsArray_Pbe)
{
  double *sol;
  BoundCond cond_Lmin, cond_Lmax;

  BoundCondition_LminLMax(cond_Lmin, cond_Lmax);
 
  //set the boundary condition at the L_Min
  if(cond_Lmin==NEUMANN)
   {
    // do nothing, it should be zero Neumann
   }
  else if(cond_Lmin==DIRICHLET)
  {
#ifdef __PBSConstT__      
    PbeFuncts[0]->Interpolate(BoundValue_LMin);  
#else    
   //nucleation 
   B_nuc(PbeFuncts[0], ScalarFunctions[0], ScalarFunctions[1]);
#endif 

//     sol = PbeFuncts[0]->GetValues();
//     memset(sol, 0, N_U*SizeOfDouble);     // no nucleation !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
//     memcpy(RhsArray_Pbe,  sol,  N_U*SizeOfDouble);
  }
  else
  {
    cout<< " Only DIRICHLET and zero NEUMANN are allowed " <<endl;
    exit(0);
  }

  //set the boundary condition at the L_Min
  if(cond_Lmax==NEUMANN)
  {
    // do nothing, it should be zero Neumann
  }
  else if(cond_Lmax==DIRICHLET)
  {
#ifdef __PBSConstT__     
    PbeFuncts[N_Intl_Levels-1]->Interpolate(BoundValue_LMax);
    sol = PbeFuncts[N_Intl_Levels-1]->GetValues();
    memcpy(RhsArray_Pbe+((N_Intl_Levels-1)*N_U),  sol,  N_U*SizeOfDouble);
#else   
    cout<< " Only Zero NEUMANN are allowed " <<endl;
    exit(0);
#endif 
  }
  else
  {
    cout<< " Only DIRICHLET and zero NEUMANN are allowed " <<endl;
    exit(0);
  }
}


void GetPhySolFromIntlVal_Nodal(int N_Intl_Levels, int N_PhySpacePts, double *PbeValues_IntlPhys,
                                double *SolPbe)
{
 int i, j, k, N_Cells;
 double *sol;
 
 for(i=0; i<N_Intl_Levels ; i++)
  {
   sol = SolPbe + (i*N_PhySpacePts);  
   for(j=0; j<N_PhySpacePts ; j++)
     sol[j] =  PbeValues_IntlPhys[j*N_Intl_Levels  + i];       
  } //for(i=0; i<N_Intl_Levels
    
} // GetPhySolFromIntlVal_Nodal



void Assemble1DInternal(TFESpace1D *FeSpace, TSquareMatrix1D *M)
{
  int i, j, k, l, N_Cells, N_BaseFunct, N_U, dGDisc;
  int N_Points, N_Sets=1, *GlobalNumbers, *BeginIndex, *DOF;
  int TestDOF, begin, end, *RowPtr, *KCol;

  double *Weights, *zeta, X[20], AbsDetjk[20];
  double LocMatrixM[MaxN_BaseFunctions1D*MaxN_BaseFunctions1D];
  double **origvaluesD0, **origvaluesD1, Mult;
  double *orgD0, *orgD1, test0, test1, ansatz0, ansatz1, *ValuesM;
  double x=0., len =0.;

  bool Needs2ndDer[1];

  TBaseCell *Cell;
  FE1D FEId;
  TFE1D *Element;
  TBaseFunct1D *bf;
  QuadFormula1D LineQuadFormula;
  TQuadFormula1D *qf1;
  TRefTrans1D *F_K;
  BaseFunct1D BaseFunct_ID, BaseFunct[1];
  TCollection *Coll;
  BoundCond BDType;
  BoundCond cond_Lmin, cond_Lmax;

  BoundCondition_LminLMax(cond_Lmin, cond_Lmax);

  dGDisc=FeSpace->IsDGSpace();
  Coll = FeSpace->GetCollection();
  N_U = FeSpace->GetN_DegreesOfFreedom();
  GlobalNumbers = FeSpace->GetGlobalNumbers();
  BeginIndex = FeSpace->GetBeginIndex();
  RowPtr = M->GetRowPtr();
  KCol = M->GetKCol();
  ValuesM = M->GetEntries();

  // all QuadPts in a cell use same FEspace in internal direction
  N_Cells = Coll->GetN_Cells();
  Needs2ndDer[0] = FALSE;

  for(i=0; i<N_Cells; i++)
  {
    Cell = Coll->GetCell(i);
    FEId = FeSpace->GetFE1D(i, Cell);
    Element = TFEDatabase2D::GetFE1D(FEId);
    bf = Element->GetBaseFunct1D();
    N_BaseFunct = Element->GetN_DOF();
    BaseFunct_ID = Element->GetBaseFunct1D_ID();

    l = bf->GetPolynomialDegree();
    LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(2*l);
    qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
    qf1->GetFormulaData(N_Points, Weights, zeta);

    F_K = TFEDatabase2D::GetRefTrans1D(LineAffin);
    ((TLineAffin *)F_K)->SetCell(Cell);
    ((TLineAffin *)F_K)->GetOrigFromRef(N_Points, zeta, X, AbsDetjk);

    BaseFunct[0] = BaseFunct_ID;
    ((TLineAffin *)F_K)->GetOrigValues(N_Sets, BaseFunct, N_Points, zeta,  LineQuadFormula,  Needs2ndDer);

    origvaluesD0=TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D0);

    memset(LocMatrixM, 0, N_BaseFunct*N_BaseFunct*SizeOfDouble);

    DOF = GlobalNumbers + BeginIndex[i];

    for(j=0;j<N_Points;j++)
    {
      Mult = Weights[j]*AbsDetjk[j];
      orgD0 = origvaluesD0[j];

      // cout<< " zeta[j] " << zeta[j]  <<endl;
      //len +=Mult;
      for(k=0;k<N_BaseFunct;k++)
      {
        test0  = orgD0[k];
        // cout<< " uref " << test0  <<endl;
        for(l=0;l<N_BaseFunct;l++)
        {
          ansatz0  = orgD0[l];
          LocMatrixM[k*N_BaseFunct + l] += (Mult*ansatz0*test0);
        }
      }
    }

    //   add to global matrices
    for(j=0;j<N_BaseFunct;j++)
    {
      TestDOF = DOF[j];

      begin = RowPtr[TestDOF];
      end = RowPtr[TestDOF+1];
      for(k=begin;k<end;k++)
      {
        for(l=0;l<N_BaseFunct;l++)
        {
          if(KCol[k] == DOF[l])
          {
            ValuesM[k] +=LocMatrixM[j*N_BaseFunct + l];
            break;
          }
        }    // for(m=0;m<N_BaseFunct_low
      }    // for(n=begin;n<end;n++)
    } // for(l=0;l<N_BaseFunct_low
  }  // for(i=0; i<N_Cells; i++)

  //update boundary data
  // starting point:  in PBS
  // end point: in PBS

  if(cond_Lmin==DIRICHLET && !dGDisc)
   {
    begin = RowPtr[0];
    end = RowPtr[1];

    for(k=begin;k<end;k++)
    {
      if(KCol[k] == 0 )
        { ValuesM[k] = 1.; }
        else
          { ValuesM[k] = 0.; }
    }
   }  //if(cond_Lmin==DIRICHLET)

  if(cond_Lmax==DIRICHLET && !dGDisc)
   {
    begin = RowPtr[N_U-1];
    end = RowPtr[N_U];

    for(k=begin;k<end;k++)
    {
      if(KCol[k] == N_U-1 )
        { ValuesM[k] = 1.; }
        else
          { ValuesM[k] = 0.; }
    }
   }  //if(cond_Lmin==DIRICHLET)

  // cout<< " len " << len << endl;
  // exit(0);
  //print matrix
//   int rank;
//    MPI_Comm_rank(TDatabase::ParamDB->Comm, &rank); 
//    
//   if(rank==0)
//     for(j=0;j<N_U;j++)
//      {
//       begin = RowPtr[j];
//       end = RowPtr[j+1];
//       for(k=begin;k<end;k++)
//        {
//         cout << "M(" << j << ", "<< KCol[k] << ") = " << ValuesM[k] <<endl;
//        }
//       cout<<endl;
//      }
// //  

    
} // Assemble1DInternal




void GetOSError(int N_Intl_Levels, TFEFunction3D **ScalarFunctions,  TFESpace1D *FESpace1D,
                double *SolPbe, int N_U, double *errors)
{
 int i, j, k, l, m, P, N_Cells, N_Cells_Intl, N_LDof, N_LocalUsedElements, *N_BaseFunct;
 int *BeginIndex, *GlobalNumbers, *BeginIndex_Intl, *GlobalNumbers_Intl;
 int N_Points, N_LocalDOFs, *DOF, N_BaseFunct_Intl, *DOF_Intl, N_LinePoints, N_Sets=1;
  
 double *weights, *xi, *eta, *zeta, *Sol_AllL;
 double X[MaxN_QuadPoints_3D], Y[MaxN_QuadPoints_3D], Z[MaxN_QuadPoints_3D];
 double AbsDetjk[MaxN_QuadPoints_3D];
 double *Sol_QuadIntl, *sol;
 double **OrigFEValues, *Orig, value, *InternalError;
 double *LineWeights, *beta, L[MaxN_QuadPoints_1D], LineAbsDetjk[MaxN_QuadPoints_1D];;
 double **origvaluesD0, **origvaluesD1, *orgD0, *orgD1, Mult;
 double x, y, z, ell, Exact[6], *sol_QI, valuegrad;

 
 TBaseCell *cell, *Cell_Intl;
 TCollection *Coll, *Coll_Intl;
 TFESpace3D *FESpace;  
 BaseFunct3D BaseFunct, *BaseFuncts;
 TFE3D *Element;
 FE3D LocalUsedElements[1], CurrentElement;
 FE1D FEId_Intl;
 TFE1D *Element_Intl;
 TBaseFunct1D *bf_Intl;
 BaseFunct1D BaseFunct_ID_Intl, BaseFunct_Intl[1];
 QuadFormula1D LineQuadFormula;
 TQuadFormula1D *qf1;
 TRefTrans1D *F_K;
  
 bool *SecondDer;
 bool Needs2ndDer[1];
  
  FESpace = ScalarFunctions[0]->GetFESpace3D();
  BeginIndex = FESpace->GetBeginIndex();
  GlobalNumbers = FESpace->GetGlobalNumbers();
  Coll = FESpace->GetCollection();
  N_Cells = Coll->GetN_Cells();
    
  Coll_Intl = FESpace1D->GetCollection();
  BeginIndex_Intl = FESpace1D->GetBeginIndex();
  GlobalNumbers_Intl = FESpace1D->GetGlobalNumbers();
  N_Cells_Intl = Coll_Intl->GetN_Cells();
  N_LDof = FESpace1D->GetN_DegreesOfFreedom();
  N_LocalUsedElements = 1;
  SecondDer = new bool[1];
  SecondDer[0] = FALSE;
  Needs2ndDer[0] = FALSE;
  
  BaseFuncts = TFEDatabase3D::GetBaseFunct3D_IDFromFE3D();
  N_BaseFunct = TFEDatabase3D::GetN_BaseFunctFromFE3D();
  
  //first check how many quad pts  in the cell
  cell = Coll->GetCell(0);
  LocalUsedElements[0] = FESpace->GetFE3D(0, cell);
  Element = TFEDatabase3D::GetFE3D(LocalUsedElements[0]);
  TFEDatabase3D::GetOrig(N_LocalUsedElements, LocalUsedElements, Coll, cell, (boolean *)SecondDer,
                         N_Points, xi, eta, zeta, weights, X, Y, Z, AbsDetjk);
  
  Sol_QuadIntl = new double[N_Points*N_Intl_Levels];
  Sol_AllL = new double[N_Points*N_LDof];
  InternalError = new double[N_Points];
//   InternalErrorgrad = new double[N_Points];

  errors[0] = 0.;
  errors[1] = 0.;
  
  for(i=0; i<N_Cells; i++)
   {
    cell = Coll->GetCell(i);
    LocalUsedElements[0] = FESpace->GetFE3D(i, cell);
    Element = TFEDatabase3D::GetFE3D(LocalUsedElements[0]);
    
    // ====================================================================
    // calculate values on original element
    // ====================================================================
    TFEDatabase3D::GetOrig(N_LocalUsedElements, LocalUsedElements, Coll, cell, (boolean *)SecondDer,
                           N_Points, xi, eta, zeta, weights, X, Y, Z, AbsDetjk);
    // calculate all needed derivatives of this FE function
    CurrentElement = FESpace->GetFE3D(i, cell);
    BaseFunct = BaseFuncts[CurrentElement];
    N_LocalDOFs = N_BaseFunct[CurrentElement];
    DOF = GlobalNumbers+BeginIndex[i];

    OrigFEValues = TFEDatabase3D::GetOrigElementValues(BaseFunct, D000);

    // find values at all quad points
    memset(Sol_QuadIntl, 0, N_Intl_Levels*N_Points*SizeOfDouble);
    for(j=0; j<N_Intl_Levels; j++)
    {
      sol = SolPbe+j*N_U;

      for(k=0; k<N_Points; k++)
      {
        Orig = OrigFEValues[k];
        value = 0.;

        for(l=0; l<N_LocalDOFs; l++)
          value += sol[DOF[l]]*Orig[l];

        Sol_QuadIntl[k*N_Intl_Levels  +  j] = value;
      }
    }      //for(j=0; j<N_Intl_Levels; j++)    
    
    
    IntlNodal2Sol(FESpace1D, Sol_QuadIntl, N_Points, N_Intl_Levels, Sol_AllL);

    memset(InternalError, 0, N_Points*SizeOfDouble);
// //     memset(InternalErrorgrad, 0, N_Points*SizeOfDouble);
    
    for(j=0; j<N_Points; j++)
    {
      x = X[j];
      y = Y[j];
      z = Z[j];
      
      sol_QI = Sol_AllL + j*N_LDof;
      
      for(k=0; k<N_Cells_Intl; k++)
      {
        Cell_Intl = Coll_Intl->GetCell(k);
        FEId_Intl = FESpace1D->GetFE1D(k, Cell_Intl);
        Element_Intl = TFEDatabase2D::GetFE1D(FEId_Intl);
        bf_Intl = Element_Intl->GetBaseFunct1D();
        N_BaseFunct_Intl = Element_Intl->GetN_DOF();
        BaseFunct_ID_Intl = Element_Intl->GetBaseFunct1D_ID();
        DOF_Intl = GlobalNumbers_Intl+BeginIndex_Intl[k];

        P = bf_Intl->GetPolynomialDegree();
        LineQuadFormula = TFEDatabase2D::GetQFLineFromDegree(2*P);
        qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
        qf1->GetFormulaData(N_LinePoints, LineWeights, beta);

        F_K = TFEDatabase2D::GetRefTrans1D(LineAffin);
        ((TLineAffin *)F_K)->SetCell(Cell_Intl);
        ((TLineAffin *)F_K)->GetOrigFromRef(N_LinePoints, beta, L, LineAbsDetjk);

        BaseFunct_Intl[0] = BaseFunct_ID_Intl;
        ((TLineAffin *)F_K)->GetOrigValues(N_Sets, BaseFunct_Intl, N_LinePoints, zeta,  LineQuadFormula,  Needs2ndDer);

        origvaluesD0=TFEDatabase2D::GetOrigElementValues(BaseFunct_ID_Intl, D0);
        origvaluesD1=TFEDatabase2D::GetOrigElementValues(BaseFunct_ID_Intl, D1);
        DOF_Intl = GlobalNumbers_Intl + BeginIndex_Intl[k];

        for(l=0; l<N_LinePoints; l++)
        {
          Mult = 0.5*LineWeights[l]*LineAbsDetjk[l];
          orgD0 = origvaluesD0[l];
          orgD1 = origvaluesD1[l];
          ell = L[l];

          //find sol at this point
          value = 0.; valuegrad= 0.;
          for(m=0; m<N_BaseFunct_Intl; m++)
          {
            value += sol_QI[DOF_Intl[m]]*orgD0[m];
            valuegrad  += sol_QI[DOF_Intl[m]]*orgD1[m];
          }

          Exact_Psd_Intl(x,  y,  z, ell, Exact);
          InternalError[j] += Mult*(Exact[0]-value)*(Exact[0]-value);
// //           InternalErrorgrad[j] += Mult*(Exact[1]-valuegrad)*(Exact[1]-valuegrad);

//           if(j==0 && l==0)
//            cout <<  j << " " << k << " " << Exact[0] << " " << value << endl;

        }  //for(l=0; l
      } // for(k=0; k<N_Cells_Intl; k++)
      
      Mult = weights[j]*AbsDetjk[j];
      errors[0] += Mult*InternalError[j];
// //       errors[1] += Mult*InternalErrorgrad[j];      
    }// for(j=0; j<N_Points; j++)
   
   } // for(i=0; i<N_Cells; i++)


  delete [] Sol_QuadIntl;
  delete [] InternalError;
  delete [] Sol_AllL;
  delete [] SecondDer;

#ifndef _MPI // sqrt(errors[j]) in the main programm after collecting error from all subdomains
    errors[0] = sqrt(errors[0]);
#endif
}//GetOSError


int main(int argc, char* argv[])
{
  // ======================================================================
  // variable declaration
  // ======================================================================
  double t_par1, t_par2, time1, start_time, t_dt;  
  
#ifdef _MPI
  const int root = 0;
  int rank, size;
  int MaxCpV, MaxSubDomainPerDof;
  int  out_rank, N_Cells_loc;

  TParVector3D  **ParSolVect, **ParRhsVect, *ParSolPbe, *ParRhsPbe;
  TParFECommunicator3D **ParComm, *ParVeloComm;

  MPI_Request request001, request002, request003, request004, request005;
  MPI_Status status;

  MPI_Comm Comm = MPI_COMM_WORLD;
  MPI_Init(&argc, &argv);

  time1 =  MPI_Wtime();
  MPI_Allreduce(&time1, &start_time, 1, MPI_DOUBLE, MPI_MIN, Comm);

  MPI_Comm_rank(Comm, &rank);
  MPI_Comm_size(Comm, &size);
#endif

  TDomain *Domain = new TDomain();
  TDomain *Domain_Intl = new TDomain();
  TDatabase *Database = new TDatabase();
#ifdef _MPI
  TDatabase::ParamDB->Comm = Comm;
  TParVectorNSE3D *ParVeloVect;
  TScalar_ParSolver **Par_Solver; 
#endif
  TFEDatabase3D *FEDatabase = new TFEDatabase3D();
  TFEDatabase2D *FEDatabase2D1D = new TFEDatabase2D();

  int i, ii, j, k, l, m, N_SubSteps, ret, ORDER, VELOCITYORDER, PBE_INDEX, img=1, N_VeloDof;
  int N_PBEqns, N_ScalarEqns, N_IndepntScalarEqns, N_U, PhyDofTimesIntlLevels;
  int Disctypes[4], N_Cells, *N_Active, *N_DofArray, VeloFunction, N_PBEAct;
  int N, N_Cells_Intl, N_Intl_Levels, N_PhySpacePts, N_PhySpaceIntlPts_All, Out_Level, start_pbe, end_pbe, dGDisc=0;
  int N_FESpaces, N_Rhs, N_SquareMatrices, N_PhyTimesIntlPts, Max_N_Unknowns=0;
  int time_discs, very_first_time=0, temp_var, N_Intl_Dof;
  int **neum_to_diri, *N_neum_to_diri, N_OwnDofs, *OwnDofs, N_f_Integ_NodalPt, *f_Integ_NodalPt, N_Q3Data=1;
  int OrderDiff;
  
  double t1, t2, t3, total_time, L0, L1, oldtau, end_time, gamma, tau;
  double **B, **RhsArray, **SolArray, *defect, **OldSolArray;
  double *IntlX, *IntlY, *IntlZ, *Sol_IntlLoc, *OldSol_IntlLoc, *B_IntlLoc, *defect_IntlLoc;
  double *Scalar_Entries, l2, H1, *PbeValues_IntlPhys, *RhsPbe_Intl_Loc, *OldPBE_IntlPtValuesT, *C_PhyPtValues;
  double *IntlSol_All_SUM, *IntlSol_All_SUM_T, *IntlSol_tmp2, q3_max;
  double *PBE_IntlPtValues, *PBE_IntlPtRhsValues, *T_PhyPtValues, *velo_PhyPtValues, *grad_velo_PhyPtValues;
  double len, *C_Sat_PhyPtValues, *SolPbe_OneLevel, *PBE_IntlPtRhsValuesT, *Rhs_Intl, *AggrRhs_Intl;
  double *SolPbe_Lmin, *SolPbe_Lmax, *L_Cube, *SolPbe, *B_Pbe, *velocity;
  double *Sol_Loc, *MatValues, *IntlPosL, *RhsArray_Pbe, *OldRhsArray_Pbe, *layers_mass;
  double *XNodalPts, *RHSs[3], *SolPbe_Intl_Loc, *Sol_Intl;
  double T, C, C_Sat, *velo_PhysPts, *grad_velo_PhysPts, params[8], hmin, hmax, hmin_all;    
  double errors[7],  *olderror, *olderror1, L2error_Max_t, L2error_Max=-1.e-8, *L2T, *H1T; 
  double L2error_Max_Pbe=-1.e-8, L2error_Max_t_Pbe, L2T_Pbe=0., H1T_Pbe=0., olderror_Pbe;
  double **neum_to_diri_x, **neum_to_diri_y, **neum_to_diri_z;
  double **lump_mass, **matrix_D_Entries, **tilde_u, **oldrhs_fem_fct0, **oldrhs_fem_fct1;
  double lpcoeff, lpexponent, integral, val_left, val_right, L;
  
  BoundCondFunct3D *BoundaryConditions[4];
  BoundCondFunct3D *BDCond[4];
  DoubleFunct3D *InitiaValues[4];
  BoundValueFunct3D *BoundValues[4];
  BoundValueFunct3D *BDValue[4];
  CoeffFct3D *Coefficients[4];
  MatVecProc *MatVect;
  DefectProc *Defect;
  TDiscreteForm3D *DiscreteForm;
  TDiscreteForm3D *DiscreteFormMatrixMRhs[4], *DiscreteFormMatrixMRhs_SUPG[4];
  TDiscreteForm3D *DiscreteFormMatrixARhs[4], *DiscreteFormMatrixARhs_SUPG[4];
  TFESpace3D **Scalar_Spaces, *velocity_space, *fesp[2], *ferhs[1];
  TFESpace1D *FeSpace_Intl;
  TSquareStructure3D **sqstructureA;
  TSquareMatrix3D *sqmatrixM, *sqmatrixA, **MatricesA, **MatricesM, *SQMATRICES[3];
  TSquareMatrix3D *sqmatrixK, **MatricesK;
  TSquareMatrix **sqmatrices = (TSquareMatrix **)SQMATRICES;  
  TSquareStructure1D *SqStruct1D_Intl;
  TSquareMatrix1D *M_Intl, *A_Intl, *S_Intl, *K_Intl;
  TFEFunction3D **ScalarFunctions, **ScalarPbeFunctions, *u1, *u2, *u3, *fefct[3];
  TCollection *coll, *Coll_Intl;
  TADISystem1D **ADI_System;
  TFEFunction1D *FeFunction_Intl, *GridFEFunction1D_Intl;
  BoundCond cond_Lmin, cond_Lmax;
  TOutput3D *Output;
  TFEVectFunct3D *u, *PhySpaceNodalPts; 
  TAuxParam3D *aux; 

  char *PRM, *GEO, *VtkBaseName, RankBaseName[1000], *ReadGrapeBaseName, test;
  char ReadinDat[] = "readin.dat"; 
  char MassMatrix[] = "Mass matrix";
  char Mass[] = "Mass";
  char Name[] = "name";
  char UString[] = "U";  
  char PosString[] = "pos";    
  char NameStrings[5][10]= {"T", "C1", "PSD", "C3", "C4"};
  char VString[] = "V";
  char GString[] = "IntlGrid";
  char LevelString[] = "PSD_L_";
  char ZeroString[] = "0";  
  char LevelString_Orig[] = "PSD_L_";    
  char SubID[] = "";

  bool UpdateConvection=FALSE, UpdateRhs=FALSE, ConvectionFirstTime=TRUE;

  bool UpdatePBEConvection=FALSE, UpdatePBERhs=FALSE, PBEConvectionFirstTime=TRUE;  
  bool *DirichletBDPt, FACTORIZE=TRUE;
  
  const char psddir[] = "PSD_DATA";
  time_t rawtime;
  struct tm * timeinfo;  
  
  std::ostringstream os;
  os << " ";
  //======================================================================
  // read parameter file
  //======================================================================

  total_time = GetTime();
  if(argc>=2)
   { ret=Domain->ReadParam(argv[1]);  }
  else
   { ret=Domain->ReadParam(ReadinDat); }

  if(ret==-1)
  {
    exit(-1);
  }

  OpenFiles();
  OutFile.setf(std::ios::scientific);

  

#ifdef _MPI
  out_rank=TDatabase::ParamDB->Par_P0;

  if(rank==out_rank)
#endif
   {
    Database->WriteParamDB(argv[0]);
    Database->WriteTimeDB();
   }
  ExampleFile();

  VtkBaseName = TDatabase::ParamDB->VTKBASENAME;
  ReadGrapeBaseName = TDatabase::ParamDB->READGRAPEBASENAME;
  
  GetExampleFileData(BoundaryConditions, BoundValues, InitiaValues, Coefficients,
                     N_PBEqns, N_IndepntScalarEqns, Disctypes);

  N_ScalarEqns=N_PBEqns+N_IndepntScalarEqns;
  
  if(TDatabase::ParamDB->MEASURE_ERRORS)
   {
    olderror = new double[N_ScalarEqns];
    olderror1 = new double[N_ScalarEqns];
    H1T = new double[N_ScalarEqns];
    L2T = new double[N_ScalarEqns];
    
    for(i=0; i<N_ScalarEqns; i++)
     {
      H1T[i] = 0;
      L2T[i] = 0;      
     }
     
   }
   
  // pointers to the routines which compute matrix-vector
  // products and the defect
  MatVect = MatVect_Scalar;
  Defect = Defect_Scalar;

  //======================================================================
  // initialize discrete forms
  //======================================================================
  for(i=0; i<N_ScalarEqns; i++)
  {
    // discrete form for assembling mass matrix and rhs (Galerkin FEM)
    DiscreteFormMatrixMRhs[i] = new TDiscreteForm3D
             (MassMatrix, Mass, N_Terms_MatrixMRhs, Derivatives_MatrixMRhs,
              SpacesNumbers_MatrixMRhs, N_Matrices_MatrixMRhs, N_Rhs_MatrixMRhs,
              RowSpace_MatrixMRhs, ColumnSpace_MatrixMRhs, RhsSpace_MatrixMRhs,
              MatrixMRhsAssemble, Coefficients[i], NULL);

    // discrete form for assembling stiffness matrix and rhs (Galerkin FEM)
    DiscreteFormMatrixARhs[i] = new TDiscreteForm3D
            (MassMatrix, Mass, N_Terms_MatrixARhs, Derivatives_MatrixARhs,
             SpacesNumbers_MatrixARhs, N_Matrices_MatrixARhs, N_Rhs_MatrixARhs,
             RowSpace_MatrixARhs, ColumnSpace_MatrixARhs, RhsSpace_MatrixARhs,
             MatrixARhsAssemble, Coefficients[i], NULL);
  }

  //======================================================================
  // read boundary parameterization and initialize coarse grid
  //======================================================================
  BoundCondition_LminLMax(cond_Lmin, cond_Lmax);
  PRM = TDatabase::ParamDB->BNDFILE;
  GEO = TDatabase::ParamDB->GEOFILE; 
  Domain->Init(PRM, GEO);

  for(i=0;i<TDatabase::ParamDB->UNIFORM_STEPS;i++)
    Domain->RegRefineAll();
 
  //======================================================================
  // Partition grid using Metis
  //======================================================================
#ifdef _MPI
  Domain->GenerateEdgeInfo();

  t_par1 = MPI_Wtime();
  Partition_Mesh3D(Comm, Domain, MaxCpV);
  t_par2 = MPI_Wtime();

  if(rank==0)
    printf("Time taken for Domain Decomposition is %e\n", (t_par2-t_par1));

  MaxSubDomainPerDof = MIN(MaxCpV, size);
#endif
  
  t3 = GetTime();
  total_time = t3 - total_time; 
  SetPolynomialDegree();

  // check the example file, to activate
  VeloFunction=TDatabase::ParamDB->INTERNAL_CONVECTION_EQ_VELOFIELD;

  Scalar_Spaces = new TFESpace3D*[N_ScalarEqns];
  sqstructureA = new TSquareStructure3D*[N_ScalarEqns];

  N_DofArray = new int[N_ScalarEqns];
  N_Active = new int[N_ScalarEqns];
  MatricesA = new TSquareMatrix3D*[N_ScalarEqns];
  MatricesM = new TSquareMatrix3D*[N_ScalarEqns];
  MatricesK = new TSquareMatrix3D*[N_ScalarEqns];
  SolArray = new double*[N_ScalarEqns];
  OldSolArray = new double*[N_ScalarEqns];  
  RhsArray = new double*[N_ScalarEqns]; 

  B = new double*[N_ScalarEqns];

  ScalarFunctions = new TFEFunction3D*[N_ScalarEqns];

  
#ifdef __FEMFCT__  
  N_neum_to_diri = new int [N_ScalarEqns];
  neum_to_diri = new int*[N_ScalarEqns];
  neum_to_diri_x = new double*[N_ScalarEqns];
  neum_to_diri_y = new double*[N_ScalarEqns];
  neum_to_diri_z = new double*[N_ScalarEqns];  
  
  lump_mass = new double*[N_ScalarEqns];
  matrix_D_Entries = new double*[N_ScalarEqns];   
  tilde_u = new double*[N_ScalarEqns]; 
  oldrhs_fem_fct0 = new double*[N_ScalarEqns];   
  oldrhs_fem_fct1 = new double*[N_ScalarEqns];   
#endif   
  
#ifdef _MPI
  ParSolVect = new TParVector3D*[N_ScalarEqns];
  ParRhsVect = new TParVector3D*[N_ScalarEqns];

  ParComm = new TParFECommunicator3D*[N_ScalarEqns];
  Par_Solver = new TScalar_ParSolver*[N_ScalarEqns];
  
  coll=Domain->GetOwnCollection(It_Finest, 0, rank);
#else
  coll=Domain->GetCollection(It_Finest, 0);
#endif  

  N_Cells = coll->GetN_Cells();  
//   printf("Rank %d: N_Cells  : %d \n",rank, N_Cells); 
   

  // ORDER should be same for all scalar and Physical PBE spaces 
  // due to coupling
  ORDER = TDatabase::ParamDB->ANSATZ_ORDER;
  VELOCITYORDER = TDatabase::ParamDB->VELOCITY_SPACE;
    
 
#ifdef _MPI
    TDatabase::ParamDB->SOLVER_TYPE = 101;
#else
    TDatabase::ParamDB->SOLVER_TYPE = 2;
#endif
  TDatabase::TimeDB->CURRENTTIME=0.;
 
  
  //=========================================================================
  // construct all finite element spaces
  //=========================================================================
  for(i=0;i<N_IndepntScalarEqns;i++)
  {
           
    TDatabase::ParamDB->PBE_P9=i;
     
    // fespaces for scalar equations
    Scalar_Spaces[i] =  new TFESpace3D(coll, Name, NameStrings[i], BoundaryConditions[i], ORDER); // ORDER
    N_DofArray[i] = Scalar_Spaces[i]->GetN_DegreesOfFreedom();
    N_Active[i] = Scalar_Spaces[i]->GetActiveBound();  
        
 
//    if(i>0) 
//    {
// #ifdef _MPI      
//    if(rank==0)
//     printf("Main Programm  N_Intl_Levels %d \n", N_Intl_Levels );  
//     MPI_Finalize(); 
// #endif    
//   exit(0);    
//    }    
    
   if(Max_N_Unknowns<N_DofArray[i])  Max_N_Unknowns=N_DofArray[i];   
#ifdef _MPI

    t_par1 = MPI_Wtime();

    Scalar_Spaces[i]->SetMaxSubDomainPerDof(MaxSubDomainPerDof);
    ParComm[i] = new TParFECommunicator3D(Comm, Scalar_Spaces[i]);
  
    t_par2 = MPI_Wtime();
    if(rank==out_rank)
     {
      printf("Time taken for FeSpace SubDomain dof mapping %e\n", (t_par2-t_par1));
      printf("DOF of FeSpace  space %d : %d \n", i, ParComm[i]->GetN_GlobalDegreesOfFreedom());      
     }
#else
    OutPut(" Rank " <<  " DOF Scalar : " << i << setw(10) << N_DofArray[i] << endl);
#endif

    //=========================================================================
    // memory allocate all vectors and construction of all fefunction
    //=========================================================================    
    SolArray[i] =  new double[N_DofArray[i]];
    OldSolArray[i] =  new double[N_DofArray[i]];    
    B[i] = new double [N_DofArray[i]];    
    
#ifdef _MPI
    ParSolVect[i] =  new TParVector3D(Comm, SolArray[i], N_DofArray[i], 1, ParComm[i]);
    ParRhsVect[i] =  new TParVector3D(Comm, B[i], N_DofArray[i], 1, ParComm[i]);
#endif

    memset(SolArray[i], 0, N_DofArray[i]*SizeOfDouble);    
    memset(B[i], 0, N_DofArray[i]*SizeOfDouble);

    ScalarFunctions[i] = new TFEFunction3D(Scalar_Spaces[i], NameStrings[i], NameStrings[i], SolArray[i], N_DofArray[i]);

    ScalarFunctions[i]->Interpolate(InitiaValues[i]);
    memcpy(OldSolArray[i], SolArray[i], N_DofArray[i]*SizeOfDouble);      
    
    RhsArray[i] = new double[N_DofArray[i]];
        
    //=========================================================================
    // allocate memory for all matrices
    //=========================================================================
     // first build matrix structure
     sqstructureA[i] = new TSquareStructure3D(Scalar_Spaces[i]);
     sqstructureA[i]->Sort();  // sort column numbers: numbers are increasing

     // two matrices used
     // M is the mass matrix, also a system matrix
     // the iterative solver uses M
     sqmatrixM = new TSquareMatrix3D(sqstructureA[i]);
     MatricesM[i] = sqmatrixM;

     // A contains the non time dependent part of the discretization
     sqmatrixA = new TSquareMatrix3D(sqstructureA[i]);
     MatricesA[i] = sqmatrixA;
   
 
# ifdef _MPI
     // initialize the parallel solver
     Par_Solver[i] = new TScalar_ParSolver(ParComm[i], sqstructureA[i], 1);
#endif     
           
  } // for(i=0;i<N_IndepntScalarEqns;i++)

 
// #ifdef _MPI      
//    if(rank==0)
//     printf("Main Programm  N_IndepntScalarEqns %d \n", N_IndepntScalarEqns );  
//     MPI_Finalize(); 
// #endif    
//   exit(0);

   
 //=========================================================================
 // NSE space
 //  each scalar equation including PBE in space will use same velo function
 //=========================================================================
  if(VeloFunction)
   {
    velocity_space =  new TFESpace3D(coll, Name, UString, BoundCondition_NSE, 
                                     TDatabase::ParamDB->VELOCITY_SPACE);
                                     
#ifdef _MPI
    i = TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE;     
    MPI_Allreduce(&i, &TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE,
                    1, MPI_INT, MPI_MIN, Comm);                         

    velocity_space->SetMaxSubDomainPerDof(MaxSubDomainPerDof);
    ParVeloComm = new TParFECommunicator3D(Comm, velocity_space);
#endif
  
    N_VeloDof = velocity_space->GetN_DegreesOfFreedom();
    velocity = new double [3*N_VeloDof];                                  
    memset(velocity, 0, 3*N_VeloDof*SizeOfDouble);                                     
#ifdef _MPI
    //  velocity only
    ParVeloVect =  new TParVectorNSE3D(Comm, velocity, N_VeloDof, 0, 3, ParVeloComm, NULL);
# endif                                     
                                     
    u = new TFEVectFunct3D(velocity_space, UString, UString, velocity, N_VeloDof, 3);
    u1 = u->GetComponent(0);
    u2 = u->GetComponent(1);
    u3 = u->GetComponent(2);
    
#ifdef _MPI    
    sprintf (RankBaseName, "%d", rank);
    strcat(ReadGrapeBaseName, RankBaseName);
    strcat(ReadGrapeBaseName, ".Sol");     
#else    
    strcat(ReadGrapeBaseName, "0.Sol");   
#endif

    u->ReadSol(ReadGrapeBaseName);  
   } // if(VeloFunction)
   
   
   
#ifndef __PBS__ 
   // allocate arrays for solver
   defect = new double[Max_N_Unknowns];   
   
  //=========================================================================
  // PBS setting begin
  // construct FESpace and matrices for population balance equation
  //=========================================================================
#else
  i = N_IndepntScalarEqns;
  PBE_INDEX = N_IndepntScalarEqns;
  N_PBEqns = 1;
  N_ScalarEqns = N_IndepntScalarEqns + N_PBEqns;
  
  //=========================================================================
  // Fespaces for population balance equation in physical space
  //=========================================================================
  Scalar_Spaces[PBE_INDEX] = new TFESpace3D(coll, Name, NameStrings[PBE_INDEX], BoundaryConditions[PBE_INDEX], ORDER);

  N_DofArray[PBE_INDEX] = Scalar_Spaces[PBE_INDEX]->GetN_DegreesOfFreedom();
  N_Active[PBE_INDEX] = Scalar_Spaces[PBE_INDEX]->GetActiveBound();  

  if(Max_N_Unknowns<N_DofArray[PBE_INDEX])  Max_N_Unknowns=N_DofArray[PBE_INDEX];   
     
#ifdef _MPI
    t_par1 = MPI_Wtime();
    Scalar_Spaces[PBE_INDEX]->SetMaxSubDomainPerDof(MaxSubDomainPerDof);
    ParComm[PBE_INDEX] = new TParFECommunicator3D(Comm, Scalar_Spaces[PBE_INDEX]);
    t_par2 = MPI_Wtime();
    
    if(rank==out_rank)
     {
      printf("Time taken for PBE physical space dof mapping %e\n", (t_par2-t_par1));
      printf("DOF of PBE physical space : %d \n", ParComm[PBE_INDEX]->GetN_GlobalDegreesOfFreedom());
     }
#else
   OutPut("DOF of PBE physical space : "<< setw(10) << N_DofArray[PBE_INDEX] << endl);
#endif 
      
   //=========================================================================
   // Memory for all matrices
   // Assume that the convection and the reaction coefficient terms are 
   //  independent of internal coordinates, so that lhs matrices are same
   //  for all internal lelvels of internal coordinate
   //=========================================================================
    sqstructureA[PBE_INDEX] = new TSquareStructure3D(Scalar_Spaces[PBE_INDEX]);
    sqstructureA[PBE_INDEX]->Sort();  

    sqmatrixM = new TSquareMatrix3D(sqstructureA[PBE_INDEX]);
    MatricesM[PBE_INDEX] = sqmatrixM;

    sqmatrixA = new TSquareMatrix3D(sqstructureA[PBE_INDEX]);
    MatricesA[PBE_INDEX] = sqmatrixA;

    //=========================================================================
    // Construct FeSpace and all data for the internal domain
    // Internal FESpace is same for all  QuadPts/NodalPts of physical space
    //=========================================================================
    N = (int)TDatabase::ParamDB->REACTOR_P11;
    L0 = TDatabase::ParamDB->REACTOR_P12;
    L1 = TDatabase::ParamDB->REACTOR_P13;

    Generate1DMesh(Domain_Intl, L0, L1, N);

    Coll_Intl = Domain_Intl->GetCollection(It_Finest, 0);
    N_Cells_Intl= Coll_Intl->GetN_Cells();

#ifdef _MPI
    if(rank==out_rank)
#endif
     printf("N_Cells_Internal : %d \n", N_Cells_Intl); 

    FeSpace_Intl = new TFESpace1D(Coll_Intl, VString, VString, TDatabase::ParamDB->TEST_ORDER);

    if(TDatabase::ParamDB->TEST_ORDER<-9)
     {
      FeSpace_Intl->SetAsDGSpace();
      dGDisc = 1;
      TDatabase::ParamDB->P10=0;  // no supg method
     }
  
    N_Intl_Dof =  FeSpace_Intl->GetN_DegreesOfFreedom();
    N_Intl_Levels = N_Intl_Dof;
    
    // N_Intl_Levels will be changed accordingly to 
    // the no of points which are needed to evaluate internal nodal functionals
    GetInternalNodalPts(FeSpace_Intl, N_Intl_Levels, IntlPosL);

    layers_mass = new double[N_Intl_Levels];
    for (ii=0;ii<N_Intl_Levels;ii++)
     {
      layers_mass[ii] = pow(IntlPosL[ii],3.0);
        /* OutPut(layers_mass[ii] );
         if (ii >0)
            OutPut(" " <<  layers_mass[ii] -  layers_mass[ii-1] << endl)
         else
            OutPut(endl);*/
      }
    
    
# ifdef _MPI
    // initialize the parallel solver
    Par_Solver[PBE_INDEX] = new TScalar_ParSolver(ParComm[PBE_INDEX], sqstructureA[PBE_INDEX],
                                                  N_Intl_Levels);
 
    if(rank==out_rank)
#endif
     printf("Dof of PBE Internal space : %d \n", N_Intl_Dof); 
     printf("N_Intl_Levels of PBE Internal space : %d \n", N_Intl_Levels);      

    SqStruct1D_Intl = new TSquareStructure1D(FeSpace_Intl);
    SqStruct1D_Intl->Sort();
    M_Intl = new TSquareMatrix1D(SqStruct1D_Intl);
    A_Intl = new TSquareMatrix1D(SqStruct1D_Intl);

    if(TDatabase::ParamDB->REACTOR_P15)
    {
      S_Intl = new TSquareMatrix1D(SqStruct1D_Intl);
      K_Intl = new TSquareMatrix1D(SqStruct1D_Intl);
    }
  
    if(TDatabase::ParamDB->REACTOR_P5==0)
     {
      N_PhySpacePts = N_DofArray[PBE_INDEX];
      XNodalPts = new double[3*N_PhySpacePts];
      PhySpaceNodalPts = new TFEVectFunct3D(Scalar_Spaces[PBE_INDEX], PosString, PosString, XNodalPts, N_PhySpacePts, 3);

      PhySpaceNodalPts->GridToData();
      IntlX = XNodalPts;
      IntlY = IntlX+N_PhySpacePts;
      IntlZ = IntlY+N_PhySpacePts;           
     }
    else
     {
      GetPhysicalSpaceQuadPts(N_PhySpacePts, Scalar_Spaces[PBE_INDEX], IntlX, IntlY, IntlZ);
     
     }
#ifdef _MPI    
      MPI_Allreduce(&N_PhySpacePts, &N_PhySpaceIntlPts_All, 1, MPI_INT, MPI_SUM, Comm);
#else
      N_PhySpaceIntlPts_All = N_PhySpacePts;
#endif  

 
  
#ifdef _MPI
    if(rank==out_rank)
#endif
     printf("N_PhySpaceIntlPts_All (including overlaping Pts on SubDomainJoints): %d \n", N_PhySpaceIntlPts_All); 

//  info for output psd sol at outlet 
#ifdef _MPI
    ParComm[PBE_INDEX]->GetOwnDofs(N_OwnDofs, OwnDofs); 
#else
    N_OwnDofs = N_PhySpacePts;
#endif    

    N_f_Integ_NodalPt = 0;
    for(i=0; i<N_OwnDofs; i++)
    {    
#ifdef _MPI    
     j = OwnDofs[i];
#else
     j =  i;
#endif         
    if(fabs(IntlX[j]-200)<1.e-12)
     {
//       OutPut(rank<< " IntlX[j] "<<IntlX[j]<< " IntlY[j] "<<IntlY[j]<<" IntlZ[j] "<<IntlZ[j]<<endl);
      N_f_Integ_NodalPt++;
     }        
    }// for(i=0; i<N_PhySpacePts; i++)
    
    f_Integ_NodalPt = new int[N_f_Integ_NodalPt];
    N_f_Integ_NodalPt = 0;
    
    for(i=0; i<N_OwnDofs; i++)
    {    
#ifdef _MPI    
     j = OwnDofs[i];
#else
     j = i;
#endif         
    if(fabs(IntlX[j]-200)<1.e-12)
     {
      f_Integ_NodalPt[N_f_Integ_NodalPt] = j;
      N_f_Integ_NodalPt++;
     }        
    }// for(i=0; i<N_PhySpacePts; i++)    
       
    //=========================================================================
    // Construct  ADI_Systems for all  Physical Space Internal Pts
    //=========================================================================  
  
    ADI_System = new TADISystem1D*[N_PhySpacePts];

    Sol_IntlLoc = new double[N_Intl_Dof];
    OldSol_IntlLoc = new double[N_Intl_Dof];
    B_IntlLoc = new double[N_Intl_Dof];
    defect_IntlLoc = new double[N_Intl_Dof];

    memset(Sol_IntlLoc, 0, N_Intl_Dof*SizeOfDouble);
    memset(OldSol_IntlLoc, 0, N_Intl_Dof*SizeOfDouble);
    memset(B_IntlLoc, 0, N_Intl_Dof*SizeOfDouble);
    memset(defect_IntlLoc, 0, N_Intl_Dof*SizeOfDouble); 

    FeFunction_Intl = new TFEFunction1D(FeSpace_Intl, NameStrings[PBE_INDEX], NameStrings[PBE_INDEX],
                                        Sol_IntlLoc, N_Intl_Dof); 

    for(i=0; i<N_PhySpacePts; i++)
      ADI_System[i] = new TADISystem1D(FeFunction_Intl, M_Intl, A_Intl, IntlX[i], IntlY[i], IntlZ[i],
                                       Sol_IntlLoc, OldSol_IntlLoc, B_IntlLoc, defect_IntlLoc,
#ifdef __PBSConstT__   
                                       NULL
#else
                                       Get_GrowthAndB_Nuc
#endif
                                       );

    C_PhyPtValues = new double[N_PhySpacePts];
    T_PhyPtValues = new double[N_PhySpacePts];
    C_Sat_PhyPtValues = new double[N_PhySpacePts];
    velo_PhysPts = new double[3*N_PhySpacePts];
    grad_velo_PhysPts = new double[9*N_PhySpacePts];
    DirichletBDPt = new bool[N_PhySpacePts];   
    
    N_PhyTimesIntlPts = N_PhySpacePts*N_Intl_Levels;    
    PbeValues_IntlPhys = new double[N_PhyTimesIntlPts];
    
    Sol_Intl =  new double[N_PhySpacePts*N_Intl_Dof];
    Rhs_Intl =  new double[N_PhySpacePts*N_Intl_Dof];   
    AggrRhs_Intl =  new double[N_PhySpacePts*N_Intl_Dof];  
    IntlSol_All_SUM =  new double[N_Intl_Dof];
    IntlSol_All_SUM_T =  new double[N_Intl_Dof];    
    IntlSol_tmp2 =  new double[N_Intl_Dof];  
    
    memset(IntlSol_All_SUM, 0,  N_Intl_Dof*SizeOfDouble);
     
  //========================================================================================
  // memory allocate for all vectors and construction of PBE fefunction in
  // physical space (for all internal levels)
  //========================================================================================
  N_U = N_DofArray[PBE_INDEX];
  PhyDofTimesIntlLevels = N_Intl_Levels*N_U;
  
  SolPbe = new double[PhyDofTimesIntlLevels];
  B_Pbe = new double[PhyDofTimesIntlLevels];
  
#ifdef _MPI
  ParSolPbe =  new TParVector3D(Comm, SolPbe, N_U, N_Intl_Levels, ParComm[PBE_INDEX]);
  ParRhsPbe =  new TParVector3D(Comm, B_Pbe, N_U, N_Intl_Levels, ParComm[PBE_INDEX]);
#endif

  memset(SolPbe, 0, PhyDofTimesIntlLevels*SizeOfDouble);
  memset(B_Pbe, 0, PhyDofTimesIntlLevels*SizeOfDouble);

  ScalarPbeFunctions = new TFEFunction3D* [N_Intl_Levels];
  for(j=0; j<N_Intl_Levels; j++)
   {   
    sprintf(RankBaseName, "%d", j); 
    
    if(j<10)
     strcat(LevelString, ZeroString);  
    
    strcat(LevelString, RankBaseName);   
       
    ScalarPbeFunctions[j] = new  TFEFunction3D(Scalar_Spaces[PBE_INDEX], LevelString,
                                               LevelString, SolPbe+j*N_U, N_U);

    memcpy(LevelString, LevelString_Orig, strlen(LevelString_Orig)+1); 
   }

  SolPbe_Lmin = new double[N_U];
  SolPbe_Lmax = new double[N_U];  
  RhsArray_Pbe = new double[PhyDofTimesIntlLevels];
  OldRhsArray_Pbe = new double[PhyDofTimesIntlLevels];
  memset(RhsArray_Pbe, 0, PhyDofTimesIntlLevels*SizeOfDouble);
  memset(OldRhsArray_Pbe, 0, PhyDofTimesIntlLevels*SizeOfDouble);  
  //=========================================================================
  // Identify the Dirichlet internal points -begin  
  //=========================================================================  

    if(TDatabase::ParamDB->REACTOR_P5==0)
     {  
      for(i=0;i<N_PhySpacePts;i++)
       if(i<N_Active[PBE_INDEX])
        { DirichletBDPt[i] = FALSE; }
       else
        {  DirichletBDPt[i] = TRUE; }      
     }
    else
     {
      // assumed that no quad points of X-direction lie on the boundary
      for(i=0; i<N_PhySpacePts; i++)
      DirichletBDPt[i] = FALSE;  
     }
     
  TDatabase::IteratorDB[It_EQ]->SetParam(Domain);
  TDatabase::IteratorDB[It_LE]->SetParam(Domain);
  TDatabase::IteratorDB[It_Finest]->SetParam(Domain);
  TDatabase::IteratorDB[It_Between]->SetParam(Domain);
  TDatabase::IteratorDB[It_OCAF]->SetParam(Domain);
  
#endif // ifdef __PBS__    

  //=========================================================================
  // prepare output
  //=========================================================================
#ifdef __PBS__   
    // prepare output, only the concentration will be saved
    Output = new TOutput3D(N_IndepntScalarEqns+2, N_IndepntScalarEqns+N_Intl_Levels, 1, 1, Domain);
#else
    // prepare output, only the concentration will be saved
    Output = new TOutput3D(N_IndepntScalarEqns+1, N_IndepntScalarEqns, 1, 1, Domain);
#endif // ifdef __PBS__   

    if(VeloFunction)
     Output->AddFEVectFunct(u);

    for(i=0;i<N_IndepntScalarEqns;i++)
      Output->AddFEFunction(ScalarFunctions[i]);

    os.seekp(std::ios::beg);
    Output->AddParameter(TDatabase::TimeDB->CURRENTTIME,os.str().c_str()); 
    
#ifdef __PBS__  
    for(j=0; j<N_Intl_Levels; j++)
     Output->AddFEFunction(ScalarPbeFunctions[j]);  
#endif // ifdef __PBS__ 
  //========================================================================================
  // Identify the dirichlet internal points - end
  // PBS setting end
  // assembling of mass matrix and initial rhs of scalar equations
  //======================================================================================== 
 
  // to make Param[]=0 in BilinearCoeffs in Examplefile since only mass matrix is assembled  
  TDatabase::ParamDB->INTERNAL_CONVECTION_EQ_VELOFIELD = 0;
  
  aux = NULL;
  for(i=0;i<N_IndepntScalarEqns;i++)
   {
    // set parameters
    N_Rhs = 1;
    N_FESpaces = 1;
    fesp[0] = Scalar_Spaces[i];
      
    // only mass matrix is assembled
    if(!aux)
     aux =  new TAuxParam3D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL); 

    // reset matrices
    N_SquareMatrices = 1;
    SQMATRICES[0] = MatricesM[i];
    SQMATRICES[0]->Reset();

    DiscreteForm = DiscreteFormMatrixMRhs[i];
    BDCond[0] = BoundaryConditions[i];
    
    if((Disctypes[i] == GALERKIN) && (TDatabase::ParamDB->UREA_SOLD_PARAMETER_TYPE == FEM_FCT_LIN))
     { BDValue[0] = BoundValue_FEM_FCT; }  
    else
     { BDValue[0] = BoundValues[i]; }

    memset(RhsArray[i], 0, N_DofArray[i]*SizeOfDouble);
    RHSs[0] = RhsArray[i];
    ferhs[0] = Scalar_Spaces[i];
    
    Assemble3D(N_FESpaces, fesp, 
               N_SquareMatrices, SQMATRICES, 
               0, NULL, 
               N_Rhs, RHSs, ferhs,
               DiscreteForm, 
               BDCond, 
               BDValue, 
               aux);   
   
    // copy Dirichlet values from rhs into sol
    memcpy(SolArray[i]+N_Active[i], RhsArray[i]+N_Active[i], (N_DofArray[i]-N_Active[i])*SizeOfDouble);        

#ifdef __FEMFCT__     
    if((Disctypes[i] == GALERKIN) && (TDatabase::ParamDB->UREA_SOLD_PARAMETER_TYPE == FEM_FCT_LIN))
     { 
      switch(i)
       {
        case 0:
          // temp
//           CheckWrongNeumannNodes_temp
          CheckWrongNeumannNodes(coll, Scalar_Spaces[i], N_neum_to_diri[i], neum_to_diri[i],
                                 neum_to_diri_x[i], neum_to_diri_y[i], neum_to_diri_z[i]);
        break;       
        case 1:
          // conc
//           CheckWrongNeumannNodes_conc
           CheckWrongNeumannNodes(coll, Scalar_Spaces[i], N_neum_to_diri[i], neum_to_diri[i],
                                 neum_to_diri_x[i], neum_to_diri_y[i], neum_to_diri_z[i]);

        break;       
              
        default: 
#ifdef _MPI      
          if(rank==0) 
#endif       
           printf("FEM_FCT implemented only for temp and conc \n");
  
#ifdef _MPI      
          MPI_Finalize(); 
#endif    
          exit(0);  
        break;
       }   
       
      oldrhs_fem_fct0[i] = new double[N_DofArray[i]];
      memcpy(oldrhs_fem_fct0[i],  RhsArray[i], N_DofArray[i]*SizeOfDouble);      
      oldrhs_fem_fct1[i] = new double[N_DofArray[i]];
      lump_mass[i] = new double [N_DofArray[i]];
      LumpMassMatrixToVector(MatricesM[i], lump_mass[i]);
      matrix_D_Entries[i] = new double[MatricesM[i]->GetN_Entries()];      
      // matrix K for copy of mass matrix
      MatricesK[i] = new TSquareMatrix3D(sqstructureA[i]);   
      tilde_u[i] = new double [N_DofArray[i]];      
      // save mass matrix in matricesK
      memcpy(MatricesK[i]->GetEntries(), MatricesM[i]->GetEntries(),
             MatricesM[i]->GetN_Entries() * SizeOfDouble);         
     }  //  if((Disctypes[i] == GALERKIN) && (TD
 
#endif    
   } //  for(i=0;i<N_IndepntScalarEqns;i++)

 if(aux)
  delete aux;    
  
  if(VeloFunction)
   TDatabase::ParamDB->INTERNAL_CONVECTION_EQ_VELOFIELD = 1;  
  
#ifdef __PBS__
  //======================================================================
  // assembling of mass matrix and initial rhs on physical space PBE
  //======================================================================
  // interpolate initial solution at L direction nodal-interpolation point levels
  memset(PbeValues_IntlPhys, 0, N_PhyTimesIntlPts*SizeOfDouble); 
  
  for(i=0;i<N_PhySpacePts;i++)
   {
    SolPbe_Intl_Loc = PbeValues_IntlPhys + i*N_Intl_Levels;
    ADI_System[i]->Interpolate(SolPbe_Intl_Loc, InitialCondition_psd_Intl);
   }    

  if(TDatabase::ParamDB->REACTOR_P5==0)
   { GetPhySolFromIntlVal_Nodal(N_Intl_Levels, N_PhySpacePts, PbeValues_IntlPhys, SolPbe); }
  else
   {
    printf("Main Programm GetPhySolFromIntlVal_Quad not yet implemented !!!\n");  
#ifdef _MPI      
    MPI_Finalize(); 
#endif 
    exit(0);
   }
  
#ifdef __FEMFCT__   
  // set Diriclet nodes
  for(i=0;i<N_IndepntScalarEqns;i++)
   if(N_neum_to_diri[i])
     {
      SetDirDofFromNeu(B[i], SolArray[i], N_neum_to_diri[i], neum_to_diri[i],
                       neum_to_diri_x[i], neum_to_diri_y[i],  neum_to_diri_z[i], BoundValues[i]);
      }      
#endif

    //set the boundary value at L_min and L_max
    if(!dGDisc)
     SetLMinLMaxBoundValue(N_Intl_Levels, ScalarPbeFunctions, N_U, ScalarFunctions, RhsArray_Pbe);
     
    RhsArray[PBE_INDEX] = new double[N_U];
    memset(RhsArray[PBE_INDEX], 0, N_U*SizeOfDouble);    

    if(cond_Lmin==DIRICHLET && !dGDisc)
     { start_pbe = 1; }
    else
     { start_pbe = 0; }

    if(cond_Lmax==DIRICHLET && !dGDisc)
     { end_pbe = N_Intl_Levels-1; }
    else
     { end_pbe = N_Intl_Levels; }    
     
    // assemble physical space (X-direction) mass mat of Pbe
    N_FESpaces = 1;    
    N_Rhs = 1;
    fesp[0] = Scalar_Spaces[PBE_INDEX];
    ferhs[0] = Scalar_Spaces[PBE_INDEX];
    N_SquareMatrices = 1;
    SQMATRICES[0] = MatricesM[PBE_INDEX];
      
    DiscreteForm = DiscreteFormMatrixMRhs[PBE_INDEX];
    BDCond[0] = BoundaryConditions[PBE_INDEX];
    BDValue[0] = BoundValues[PBE_INDEX];
     
    // assembling only mass mat        
    // to make Param[]=0 in BilinearCoeffs in Examplefile  
    TDatabase::ParamDB->INTERNAL_CONVECTION_EQ_VELOFIELD = 0;
  
    aux =  new TAuxParam3D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);
            
    memset(RhsArray_Pbe+(start_pbe*N_U), 0, (end_pbe-start_pbe)*N_U*SizeOfDouble);
    N_PBEAct = N_Active[PBE_INDEX];
    
    for(j=start_pbe; j<end_pbe; j++)
     {         
      SQMATRICES[0]->Reset();  // same mat for all rhs       
       
      TDatabase::ParamDB->REACTOR_P29=IntlPosL[j];      
      RHSs[0] = RhsArray_Pbe+j*N_U;  
  
      Assemble3D(N_FESpaces, fesp, 
                 N_SquareMatrices, SQMATRICES, 
                 0, NULL, 
                 N_Rhs, RHSs, ferhs,
                 DiscreteForm, 
                 BDCond, 
                 BDValue, 
                 aux);        

      // copy Dirichlet values to the solution
      SolPbe_OneLevel = SolPbe+j*N_U;
      memcpy(SolPbe_OneLevel+N_PBEAct, RHSs[0]+N_PBEAct, (N_U-N_PBEAct)*SizeOfDouble);      
     }// for(j=start_pbe; j<end_pbe; j++)
     
    delete aux;
   
   if(VeloFunction)
    TDatabase::ParamDB->INTERNAL_CONVECTION_EQ_VELOFIELD = 1;    
         
   memcpy(OldRhsArray_Pbe, RhsArray_Pbe, PhyDofTimesIntlLevels*SizeOfDouble);

   // allocate arrays for solver
   defect = new double[Max_N_Unknowns];   
 
  //======================================================================
  // assemble the mass matrix for the internal space (L) once for all time steps 
  //======================================================================   
    Assemble1DInternal(FeSpace_Intl, M_Intl);

    // get the initial rhs for PBE
    for(i=0;i<N_PhySpacePts;i++)
     {
      if(DirichletBDPt[i])
        continue;
      
      SolPbe_Intl_Loc = Sol_Intl + i*N_Intl_Dof;
      ADI_System[i]->AssembleInitRhs(SolPbe_Intl_Loc, BilinearCoeffs_Psd_Intl,
                                     cond_Lmin, cond_Lmax);
     }  //  for(i=0;i<N_IntlPts;i++)    
    
#endif // #ifdef __PBS__

  // parameters for time stepping scheme
  gamma = 0.;
  m = 0;
  N_SubSteps = GetN_SubSteps();
  oldtau = 1;
  end_time = TDatabase::TimeDB->ENDTIME;

  // not active : TDatabase::TimeDB->TIMESTEPLENGTH_CONTROL = 0 
  if (TDatabase::TimeDB->TIMESTEPLENGTH_CONTROL)
  { time_discs = 2;}
  else
  { time_discs = 1; }

  if(TDatabase::ParamDB->WRITE_VTK)
   {
#ifdef _MPI     
//      if(N_IndepntScalarEqns)
       Output->Write_ParVTK(Comm, img, SubID);
#else
//      if(N_IndepntScalarEqns)
//       {
        os.seekp(std::ios::beg);
        if(img<10) os << VtkBaseName<<".0000"<<img<<".vtk" << ends;
        else if(img<100) os << VtkBaseName<<".000"<<img<<".vtk" << ends;
        else if(img<1000) os << VtkBaseName<<".00"<<img<<".vtk" << ends;
        else if(img<10000) os << VtkBaseName<<".0"<<img<<".vtk" << ends;
        else  os << VtkBaseName<<"."<<img<<".vtk" << ends;
        Output->WriteVtk(os.str().c_str());
//       }  
      
      //serial output of PSD
      
#endif
 
      img++;
    } // if(TDatabase::ParamDB->WRITE_VTK)

  coll->GetHminHmax(&hmin,&hmax);
  OutPut("h_min : " << hmin << " h_max : " << hmax << endl);

#ifdef _MPI
   MPI_Allreduce(&hmin, &hmin_all, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
#else
   hmin_all = hmin;
#endif
  
  // measure interpolation errors to known solution
  if(TDatabase::ParamDB->MEASURE_ERRORS)
   {   
    for(i=0;i<7;i++)
     errors[i] = 0.;
  
     for(i=0;i<N_IndepntScalarEqns;i++)
      {    
// #ifdef __PBSConstT__  
//         ScalarFunctions[i]->Interpolate(Exact);
// #endif

        fesp[0] = Scalar_Spaces[i];
        aux =  new TAuxParam3D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);

        ScalarFunctions[i]->GetErrors(Exact, 4, AllDerivatives, 2, L2H1Errors, 
                                      Coefficients[i], aux, 1, fesp, errors);
        delete aux;

#ifdef _MPI
        MPI_Reduce(errors, &l2, 1, MPI_DOUBLE, MPI_SUM, out_rank, MPI_COMM_WORLD);
        MPI_Reduce(errors+1, &H1, 1, MPI_DOUBLE, MPI_SUM, out_rank, MPI_COMM_WORLD);
        if(rank==out_rank)
         {l2 = sqrt(l2);   H1 = sqrt(H1);}
#else
        l2 = errors[0];
        H1 = errors[1];
#endif

 
#ifdef _MPI
        if(rank==out_rank)
#endif
        {
         OutPut("Time: " << TDatabase::TimeDB->CURRENTTIME);
         OutPut(" L2: " << l2);
         OutPut(" H1-semi: " << H1 << endl); 
        }
      } // for(i=0;i<N_IndepntS
      
#ifdef __PBS__  
  for(i=0;i<7;i++)
   errors[i] = 0.;
    
   GetOSError(N_Intl_Levels, ScalarPbeFunctions,  FeSpace_Intl, SolPbe, N_U, errors);  
  
#ifdef _MPI
  MPI_Reduce(errors, &l2, 1, MPI_DOUBLE, MPI_SUM, out_rank, MPI_COMM_WORLD);
  
  if(rank==out_rank)
   l2 = sqrt(l2);   
  
#else
  l2 = errors[0];
#endif  // ifdef _MPI  


#ifdef _MPI
        if(rank==out_rank)
#endif
        {
         OutPut("Time: " << TDatabase::TimeDB->CURRENTTIME);
         OutPut(" PSD L2: " << l2 << endl); 
 
//          if(L2error_Max<l2)
//           {
//            L2error_Max=l2;
//            L2error_Max_t=TDatabase::TimeDB->CURRENTTIME;
//           } 
//   
        }  
  
  
#endif          
      
   } //   if(TDatabase::ParamDB->MEASURE_ERRORS) 

//   OutPut("read u " << Ddot(N_DofArray[i], SolArray[i], SolArray[i])<< endl); 
//   exit(0);

//  #ifdef _MPI      
//    if(rank==0) 
// #endif       
//   printf("Main Programm  N_Intl_Levels %d \n", N_Intl_Levels );
//   
// #ifdef _MPI      
//     MPI_Finalize(); 
// #endif    
//   exit(0);
  
#ifdef __PBS__   
  //due to Stationary NSE, find velo only once at the begning
  if(VeloFunction)
   GetVeloGradForPhysPts(u, Scalar_Spaces[PBE_INDEX], velo_PhysPts, grad_velo_PhysPts, N_PhySpacePts);

  params[0]= TDatabase::ParamDB->UREA_AGGR_SPATIAL;//space
  params[1]= TDatabase::ParamDB->UREA_AGGR_BROWNIAN; // brown kernel included or not
  params[2]= TDatabase::ParamDB->UREA_AGGR_POL_ORDER;//pol or constant
  params[3]= TDatabase::ParamDB->UREA_AGGR_BROWNIAN_TEMP;//brown kernel depends of temperature or not
  params[4]= TDatabase::ParamDB->UREA_AGGR_BROWNIAN_SCAL;//scal param for brown kernel
  params[5]= TDatabase::ParamDB->UREA_AGGR_SHEAR_FACTOR_TYPE;//shear depends of velocity or not
  params[6]= TDatabase::ParamDB->UREA_AGGR_SHEAR_FACTOR;//param for shear kernel
#endif

#ifdef __PBSConstT__  
 TDatabase::TimeDB->TIMESTEPLENGTH = pow(hmin_all, TDatabase::ParamDB->P9);

  UpdateRhs = TRUE;  
#ifdef __PBS__   
  UpdatePBERhs = FALSE;
#endif

#endif

#ifdef _MPI
    if(rank==out_rank)
#endif 
    {
     OutPut("TIMESTEPLENGTH : " << TDatabase::TimeDB->TIMESTEPLENGTH << endl;) 
      
     mkdir(psddir, 0777);        
    }

#ifdef __FEMFCT__  
    for(i=0;i<N_IndepntScalarEqns;i++)
      memcpy(OldSolArray[i], SolArray[i], N_DofArray[i]*SizeOfDouble); 
#endif 

  lpcoeff = TDatabase::ParamDB->LP_STREAMLINE_COEFF;  
  lpexponent = TDatabase::ParamDB->LP_STREAMLINE_EXPONENT;
  OrderDiff = (int)TDatabase::ParamDB->LP_STREAMLINE_ORDER_DIFFERENCE;
  
  

  //======================================================================
  // start of time cycle
  //======================================================================
  tau = TDatabase::TimeDB->CURRENTTIME;  
  while(TDatabase::TimeDB->CURRENTTIME< end_time)
   {                                               // time cycle
    m++;
    
  if(m<3 || m % 50 == 0)   
#ifdef _MPI    
    t_dt = MPI_Wtime();
#else
    t_dt = GetTime();
#endif     
    TDatabase::TimeDB->INTERNAL_STARTTIME = TDatabase::TimeDB->CURRENTTIME;

    for(l=0;l<N_SubSteps;l++)                     // sub steps of fractional step theta
    {
      if (!very_first_time)
      {
        SetTimeDiscParameters();
      }

      if(m==1
  #ifdef _MPI
        &&   rank==out_rank
  #endif
        )
      {
        OutPut("Theta1: " << TDatabase::TimeDB->THETA1<< endl);
        OutPut("Theta2: " << TDatabase::TimeDB->THETA2<< endl);
        OutPut("Theta3: " << TDatabase::TimeDB->THETA3<< endl);
        OutPut("Theta4: " << TDatabase::TimeDB->THETA4<< endl);
      }

      tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
      if (!very_first_time)
        TDatabase::TimeDB->CURRENTTIME += tau;


# ifdef _MPI
      if(rank==out_rank)
#endif
      {
        OutPut(endl << "CURRENT TIME: ");
        OutPut(TDatabase::TimeDB->CURRENTTIME << endl);
      }

     aux = NULL;
     for(i=0;i<N_IndepntScalarEqns;i++)
      {
       // working array for rhs is B, initialize B
       memset(B[i], 0, N_DofArray[i]*SizeOfDouble);

       // compute terms with data from previous time step
       // old rhs multiplied with current subtime step and theta3 on B
       if(TDatabase::TimeDB->THETA3!=0)
        Daxpy(N_Active[i], tau*TDatabase::TimeDB->THETA3, RhsArray[i], B[i]);

       if(UpdateConvection || UpdateRhs ||  ConvectionFirstTime )
        {
         // assemble A and rhs
         N_Rhs = 1;
         N_FESpaces = 1;
         fesp[0] = Scalar_Spaces[i];

         N_SquareMatrices = 1;
         SQMATRICES[0] = MatricesA[i];
         SQMATRICES[0]->Reset();
         DiscreteForm = DiscreteFormMatrixARhs[i];

         BDCond[0] = BoundaryConditions[i];

         if( (Disctypes[i] == GALERKIN) && (TDatabase::ParamDB->UREA_SOLD_PARAMETER_TYPE == FEM_FCT_LIN) )
          { BDValue[0] = BoundValue_FEM_FCT; }  
         else
          { BDValue[0] = BoundValues[i]; }
     
         if(!aux)
         if(VeloFunction)
          {      
           N_FESpaces = 2;
           fesp[1] = velocity_space;
           fefct[0] = u1;
           fefct[1] = u2;
           fefct[2] = u3;

           aux =  new TAuxParam3D(TimeCDParamsVeloFieldN_FESpaces,
                             TimeCDParamsVeloFieldN_Fct, 
                             TimeCDParamsVeloFieldN_ParamFct, 
                             TimeCDParamsVeloFieldN_FEValues, 
                             fesp+1, fefct, 
                             TimeCDParamsVeloFieldFct, 
                             TimeCDParamsVeloFieldFEFctIndex,
                             TimeCDParamsVeloFieldFEMultiIndex, 
                             TimeCDParamsVeloFieldN_Params,
                             TimeCDParamsVeloFieldBeginParam); 
          }
         else
          { aux =  new TAuxParam3D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL); }

         memset(RhsArray[i], 0, N_DofArray[i]*SizeOfDouble);
         RHSs[0] = RhsArray[i];
         ferhs[0] = Scalar_Spaces[i];  

         Assemble3D(N_FESpaces, fesp, 
                    N_SquareMatrices, SQMATRICES, 
                    0, NULL, 
                    N_Rhs, RHSs, ferhs,
                    DiscreteForm, 
                    BDCond, 
                    BDValue, 
                    aux);  

        // add LPS atabilization
        // same Mat for all rhs, so add LPS once at the end of the i loop
        if(TDatabase::ParamDB->ANSATZ_ORDER>100 && lpcoeff!=0)
           AddStreamlineTerm(SQMATRICES[0], u1, u2, u3, lpcoeff, lpexponent, OrderDiff);   

    
         if(i==N_IndepntScalarEqns-1)
          {
#ifdef _MPI
           if(UpdateConvection ||  ConvectionFirstTime)
            FACTORIZE=TRUE;
#endif
           ConvectionFirstTime = FALSE;
           delete aux; 
          }   
        } // if(UpdateConvection || Upda
       
     // add rhs from current sub time step to rhs array B
     Daxpy(N_Active[i], tau*TDatabase::TimeDB->THETA4,  RhsArray[i], B[i]);

     // update rhs by Laplacian and convective term from previous time step
     // scaled by current sub time step length and theta2
     // currently : M := M + gamma A
     // M = M + (- tau*TDatabase::TimeDB->THETA2)     
     if(TDatabase::TimeDB->THETA2 !=0)
      {
       MatAdd(MatricesM[i], MatricesA[i], - tau*TDatabase::TimeDB->THETA2);       
       gamma = -tau*TDatabase::TimeDB->THETA2; // set current factor of steady state matrix
      }
     else
      {gamma = 0;}
      
     // set Dirichlet values
     memcpy(SolArray[i]+N_Active[i], RhsArray[i]+N_Active[i], (N_DofArray[i]-N_Active[i])*SizeOfDouble);
    
     // defect = M * sol, B:= B + defec
     memset(defect, 0, Max_N_Unknowns*SizeOfDouble); 
     MatVectActive(MatricesM[i], SolArray[i], defect);
     Daxpy(N_Active[i], 1, defect, B[i]);    

     // RHSs[0] still available from assembling
     memcpy(B[i]+N_Active[i], RhsArray[i]+N_Active[i], (N_DofArray[i]-N_Active[i])*SizeOfDouble);

     // system matrix
     MatAdd(MatricesM[i], MatricesA[i], -gamma + tau*TDatabase::TimeDB->THETA1);
     gamma = tau*TDatabase::TimeDB->THETA1;  


    //======================================================================
    // solve linear system
    //======================================================================     
#ifndef _MPI
    t1 = GetTime();
     switch(int(TDatabase::ParamDB->SOLVER_TYPE))
      {
       case 0:
          // AMG Solver
             cout << "solver type not implemented yet !!!!!!!!!!!!!" << endl;
             exit(0);
       break;

       case 1:
          // GMG Solver
             cout << "solver type not implemented yet !!!!!!!!!!!!!" << endl;
             exit(0);  
       break;

       case 2:
//          t1 = GetTime();  
              DirectSolver(MatricesM[i], B[i], SolArray[i]);
//          t2 = GetTime();
       break;

#ifdef _OMP
       case 100:
             // pardiso
             cout << "solver type not implemented yet !!!!!!!!!!!!!" << endl;
             exit(0);
       break;
#endif
       default: 
        cout << "wrong  solver type !!!!!!!!!!!!!" << endl;
        exit(0);
       break;
      }
     t2 = GetTime();

#else // Parallel solver
//      t1 = MPI_Wtime();
     Par_Solver[i]->Solve(MatricesM[i], ParRhsVect[i], ParSolVect[i], FACTORIZE);    
//      t2 = MPI_Wtime();
#endif    

#ifdef __FEMFCT__   
     //copy sol to old solution
     memcpy(OldSolArray[i], SolArray[i], N_DofArray[i]*SizeOfDouble);    
 #endif    
     //======================================================================
     // end solve linear system
     //======================================================================

     // restore matrices
     MatAdd(MatricesM[i], MatricesA[i], -gamma);        
     gamma = 0.;
#ifdef __FEMFCT__    
   if( (Disctypes[i] == GALERKIN) && (TDatabase::ParamDB->UREA_SOLD_PARAMETER_TYPE == FEM_FCT_LIN) )
     memcpy(oldrhs_fem_fct0[i], oldrhs_fem_fct1[i], N_DofArray[i]*SizeOfDouble);  
#endif    
     
    } // for(i=0;i<N_IndepntScalarEqns;i++
    
#ifdef _MPI    
   FACTORIZE=FALSE;
#endif              
           
      
#ifdef __PBS__
   //===============================================================================
   // PBE system solution -- start
   // assume that the concentration eqn solved
   // first solve w.r.t physical space and then w.r.t internal space
   // PBE system X-direction solution --start
   //===============================================================================
   //set the boundary value at L_min and L_max for the current timestep
   if(!dGDisc)
    SetLMinLMaxBoundValue(N_Intl_Levels, ScalarPbeFunctions, N_U, ScalarFunctions, RhsArray_Pbe);

    if(UpdatePBEConvection || UpdatePBERhs ||  PBEConvectionFirstTime )
     {   
       // set parameters
       N_Rhs = 1;
       N_FESpaces = 1;
       fesp[0] = Scalar_Spaces[PBE_INDEX];
       ferhs[0] = Scalar_Spaces[PBE_INDEX];

       // reset matrices
       N_SquareMatrices = 1;
       SQMATRICES[0] = MatricesA[PBE_INDEX];

       DiscreteForm = DiscreteFormMatrixARhs[PBE_INDEX];
       BDCond[0] = BoundaryConditions[PBE_INDEX];
       BDValue[0] = BoundValues[PBE_INDEX];
     
       if(VeloFunction)
        {      
         N_FESpaces = 2;
         fesp[1] = velocity_space;
         fefct[0] = u1;
         fefct[1] = u2;
         fefct[2] = u3;

         aux =  new TAuxParam3D(TimeCDParamsVeloFieldN_FESpaces,
                             TimeCDParamsVeloFieldN_Fct, 
                             TimeCDParamsVeloFieldN_ParamFct, 
                             TimeCDParamsVeloFieldN_FEValues, 
                             fesp+1, fefct, 
                             TimeCDParamsVeloFieldFct, 
                             TimeCDParamsVeloFieldFEFctIndex,
                             TimeCDParamsVeloFieldFEMultiIndex, 
                             TimeCDParamsVeloFieldN_Params,
                             TimeCDParamsVeloFieldBeginParam); 
        }
       else
        { aux =  new TAuxParam3D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL); }

          
       memset(RhsArray_Pbe+(start_pbe*N_U), 0, (end_pbe-start_pbe)*N_U*SizeOfDouble);
       N_PBEAct = N_Active[PBE_INDEX];       
     } //  if(UpdatePBEConvection || UpdatePBERhs ||  PBEConvectionFirstTime )

   
    //set the boundary value at Dirichlet nodes for all levels in the current timestep
     for(i=start_pbe; i<end_pbe; i++)
      {   
       if(UpdatePBEConvection || UpdatePBERhs ||  PBEConvectionFirstTime )
        {

         SQMATRICES[0]->Reset();  // same mat for all rhs
         
         TDatabase::ParamDB->REACTOR_P29=IntlPosL[i];      
         RHSs[0] = RhsArray_Pbe+i*N_U;  
  
         Assemble3D(N_FESpaces, fesp, 
                 N_SquareMatrices, SQMATRICES, 
                 0, NULL, 
                 N_Rhs, RHSs, ferhs,
                 DiscreteForm, 
                 BDCond, 
                 BDValue, 
                 aux);       

         if(i==end_pbe-1)
          { 
           // add LPS atabilization
           // same Mat for all rhs, so add LPS once at the end of the i loop
           if(TDatabase::ParamDB->ANSATZ_ORDER>100 && lpcoeff!=0)
             AddStreamlineTerm(SQMATRICES[0], u1, u2, u3, lpcoeff, lpexponent, OrderDiff);     
            
#ifdef _MPI
          if(UpdatePBEConvection ||  PBEConvectionFirstTime)
            FACTORIZE=TRUE;
#endif
           PBEConvectionFirstTime = FALSE;
           delete aux;
          }
          
        }//  if(UpdatePBEConvection || UpdatePBERhs ||  PBEConvectionFirstTime )

       // copy Dirichlet values from rhs into sol
       SolPbe_OneLevel = SolPbe+i*N_U;
        memcpy(SolPbe_OneLevel+N_Active[PBE_INDEX], RhsArray_Pbe+(i*N_U +N_Active[PBE_INDEX]),
            (N_U-N_Active[PBE_INDEX])*SizeOfDouble);
        
      }//for(i=start_pbe; i<end_pbe; i++)

      memcpy(SolPbe_Lmin, SolPbe, N_U*SizeOfDouble);
      memcpy(SolPbe_Lmax, SolPbe+(N_U*(N_Intl_Levels-1)), N_U*SizeOfDouble);
      
      //====================================================================================

       gamma=0.;
       // working array for rhs is B, initialize B
       // same lhs for alll levels
       memset(B_Pbe, 0, PhyDofTimesIntlLevels*SizeOfDouble);

       for(i=start_pbe; i<end_pbe; i++)
        {
          // add rhs from current sub time step to rhs array B
          Daxpy(N_Active[PBE_INDEX], tau*TDatabase::TimeDB->THETA3, OldRhsArray_Pbe+i*N_U, B_Pbe+i*N_U);

          Daxpy(N_Active[PBE_INDEX], tau*TDatabase::TimeDB->THETA4, RhsArray_Pbe+i*N_U, B_Pbe+i*N_U);
  
          if(gamma==0. && TDatabase::TimeDB->THETA2!=0)
          {
            MatAdd(MatricesM[PBE_INDEX], MatricesA[PBE_INDEX], -tau*TDatabase::TimeDB->THETA2);
    
            // set current factor of steady state matrix
            gamma = -tau*TDatabase::TimeDB->THETA2;
          }

          SolPbe_OneLevel = SolPbe+i*N_U;
          temp_var = i*N_U +N_Active[PBE_INDEX];
          // copy Dirichlet values from rhs into sol
          memcpy(SolPbe_OneLevel+N_Active[PBE_INDEX], RhsArray_Pbe+(temp_var), (N_U-N_Active[PBE_INDEX])*SizeOfDouble);

          memset(defect, 0, Max_N_Unknowns*SizeOfDouble);
          MatVectActive(MatricesM[PBE_INDEX], SolPbe_OneLevel, defect);
          Daxpy(N_Active[PBE_INDEX], 1, defect, B_Pbe+i*N_U);

          // set Dirichlet values
          memcpy(B_Pbe+(temp_var), RhsArray_Pbe+(temp_var), (N_U-N_Active[PBE_INDEX])*SizeOfDouble);
        }  //      for(i=start_pbe; i<end_pbe; i++)

        memcpy(OldRhsArray_Pbe, RhsArray_Pbe, PhyDofTimesIntlLevels*SizeOfDouble);

        // system matrix
        MatAdd(MatricesM[PBE_INDEX], MatricesA[PBE_INDEX], -gamma + tau*TDatabase::TimeDB->THETA1);
        gamma = tau*TDatabase::TimeDB->THETA1;

      //======================================================================
      // solve PBE system with multiple rhs
      //======================================================================  
#ifndef _MPI
    t1 = GetTime();
     switch(int(TDatabase::ParamDB->SOLVER_TYPE))
      {
       case 0:
          // AMG Solver
             cout << "solver type not implemented yet !!!!!!!!!!!!!" << endl;
             exit(0);
       break;

       case 1:
          // GMG Solver
             cout << "solver type not implemented yet !!!!!!!!!!!!!" << endl;
             exit(0);  
       break;

       case 2:
//          t1 = GetTime();
             DirectSolver(MatricesM[PBE_INDEX], B_Pbe, SolPbe, end_pbe, start_pbe);
//          t2 = GetTime();
       break;

#ifdef _OMP
       case 100:
             // pardiso
             cout << "solver type not implemented yet !!!!!!!!!!!!!" << endl;
             exit(0);
       break;
#endif
       default: 
        cout << "wrong  solver type !!!!!!!!!!!!!" << endl;
        exit(0);
       break;
      }
     t2 = GetTime();

#else // Parallel solver
//      t1 = MPI_Wtime();
     Par_Solver[PBE_INDEX]->Solve(MatricesM[PBE_INDEX], ParRhsPbe, ParSolPbe, FACTORIZE);
//      t2 = MPI_Wtime(); 
#endif    

   FACTORIZE=FALSE;
   //======================================================================
   // end solve PBE system with multiple rhs 
   //======================================================================

   //set the boundary value at L_min and L_max  
   if(!dGDisc)
    SetLMinLMaxBoundValue(N_Intl_Levels, ScalarPbeFunctions, N_U, ScalarFunctions, RhsArray_Pbe);

   //restore system matrix
   MatAdd(MatricesM[PBE_INDEX], MatricesA[PBE_INDEX], -gamma);
   gamma = 0;

   //===============================================================================
   //PBE system X-direction solution -- end
   //PBE system L-direction solution -- start  
   //===============================================================================
   // get the current timestep concentration
   // T, C, FESpace_PBE, C, C_sat, N_IntlPts
   // get the current timestep BD updated PB for all internal levels
  
   // transpose the solution, since we need it on rhs of the X-direction Stage 2
   if(TDatabase::ParamDB->REACTOR_P5==0)
    {
     if(N_IndepntScalarEqns)
      GetNodalPtsValuesForIntl(ScalarFunctions, Scalar_Spaces[PBE_INDEX],
                                T_PhyPtValues, C_PhyPtValues, C_Sat_PhyPtValues, N_PhySpacePts);
 
     // at present nodal functionals are just point values
     // i.e., N_U==N_PhySpacePts, so just transpose
     GetNodalPtsValuesForIntl(N_PhySpacePts, N_Intl_Levels, SolPbe, PbeValues_IntlPhys);
      
     // get the internal nodal functionals from internal nodal point values
     IntlNodal2Sol(FeSpace_Intl, PbeValues_IntlPhys, N_PhySpacePts, N_Intl_Levels, Sol_Intl);
    }
   else
    {    
     printf("Main Programm  QuadPts method not yet implemented\n");
#ifdef _MPI      
    MPI_Finalize(); 
#endif 
     exit(0);       
    } 
 
    memset(AggrRhs_Intl, 0, N_PhySpacePts*N_Intl_Dof*SizeOfDouble);   

#ifndef __PBSConstT__      
      /* compute Aggrand Brgg */
      if(TDatabase::ParamDB->PBE_P0 == 2 || TDatabase::ParamDB->PBE_P0 == 4 )
        {
#ifdef _MPI
         if(rank==out_rank && m==1)
#endif  
          OutPut("===== PSD with Aggregation ======="<< endl);
         
          // assemble Dirichlet Rhs
          memset(Rhs_Intl, 0, N_PhySpacePts*N_Intl_Dof*SizeOfDouble);  
  
          for(i=0;i<N_PhySpacePts;i++)
           if(PSD_bound_cond_from_velo_inflow_urea(IntlX[i], IntlY[i], IntlZ[i]))     
            {  
             T = T_PhyPtValues[i];    
             C = C_PhyPtValues[i];
               
             for (ii=0;ii<N_Intl_Dof;ii++)
              Rhs_Intl[ii*N_PhySpacePts + i] = PSD_BdValue(IntlX[i], IntlY[i], IntlZ[i], IntlPosL[ii], C, T);        
         
              // printf("T %f, C %f\n", T, C);
            }   // if(PSD_bound_cond_from_velo_

          call_apply_integral_operator(N_PhySpacePts, 1, 1, N_Intl_Levels,
                                       Rhs_Intl,
                                       velo_PhysPts, grad_velo_PhysPts, T_PhyPtValues,
                                       AggrRhs_Intl, layers_mass,
                                       params);
        } // if(TDatabase::ParamDB->PBE_P0 == 2 || TDatabase::ParamDB->PBE_P0 == 4 )
#endif
       
       for(i=0;i<N_PhySpacePts;i++)
        {
         if(DirichletBDPt[i])
          continue;
   
         SolPbe_Intl_Loc = Sol_Intl + i*N_Intl_Dof;
         RhsPbe_Intl_Loc = AggrRhs_Intl + i*N_Intl_Dof;
          
         if(N_IndepntScalarEqns)
          {
           T = T_PhyPtValues[i];    
           C = C_PhyPtValues[i];
           C_Sat = C_Sat_PhyPtValues[i];
          }
         else
          {
           T = 0;      
           C = 0.;
           C_Sat = 0.;
          }          
//           cout<< i << " T " << T << endl;
         ADI_System[i]->SolveAllLevels(SolPbe_Intl_Loc, RhsPbe_Intl_Loc, 
                                        C, C_Sat, T, BilinearCoeffs_Psd_Intl, tau, 
                                        cond_Lmin, cond_Lmax, NULL);      

//    #ifdef _MPI      
//    if(rank==0)
//     printf("Main Programm  TDatabase::ParamDB->SOLVER_TYPE %d \n", TDatabase::ParamDB->SOLVER_TYPE );  
//     MPI_Finalize(); 
// #endif    
//   exit(0);    

#ifndef __PBSConstT__   
        for(j=0;j<N_Intl_Dof;j++)
         if(SolPbe_Intl_Loc[j]<0)  SolPbe_Intl_Loc[j]=0;
#endif

         // update due to the concentration coupling with PBE
         if(N_IndepntScalarEqns)
          {
           C_PhyPtValues[i] = C;
           T_PhyPtValues[i] = T;
          }
        } //  for(i=0;i<N_PhySpacePts;i++)
    
     // transpose the solution, since we need it on rhs of the X-direction Stage 2
     if(TDatabase::ParamDB->REACTOR_P5==0)
      { 
       if(N_IndepntScalarEqns)
        {
         GetSolFromNodalPtVales(ScalarFunctions[0], T_PhyPtValues, N_PhySpacePts, Scalar_Spaces[PBE_INDEX]);
         GetSolFromNodalPtVales(ScalarFunctions[1], C_PhyPtValues, N_PhySpacePts, Scalar_Spaces[PBE_INDEX]); 
        } 

       // get the internal nodal point values from the internal nodal functionals  
       Sol2IntlNodal(FeSpace_Intl, Sol_Intl, N_Intl_Levels, PbeValues_IntlPhys, N_PhySpacePts);
    
       // at present nodal functionals are just point values
       // i.e., N_U==N_PhySpacePts, so just transpose
       GetPhyNodalFuctional(N_PhySpacePts, N_Intl_Levels, PbeValues_IntlPhys, SolPbe);
      }
     else
      {      
       printf("Main Programm  QuadPts method not yet implemented\n");
#ifdef _MPI
       MPI_Finalize();
#endif  
       exit(0);      
      }      

     // no change in L if it is DIRICHLET BC, correction due to interpolation/Transpose
     if(cond_Lmin==DIRICHLET && !dGDisc)
       memcpy(SolPbe, SolPbe_Lmin, N_U*SizeOfDouble);

     if(cond_Lmax==DIRICHLET && !dGDisc)
       memcpy(SolPbe+(N_U*(N_Intl_Levels-1)), SolPbe_Lmax, N_U*SizeOfDouble);

#ifdef _MPI     
     // averaging dof values on subdomain interface
     ParSolPbe->AssembleByAverage();

     if(N_IndepntScalarEqns)
      {
       ParSolVect[0]->AssembleByAverage(); // heat
       ParSolVect[1]->AssembleByAverage(); // concentration
      }
#endif  
    //===============================================================================
    //PBE system L-direction solution --end
    //=============================================================================== 
#endif // __PBS__    
    } // for(l=0;l<N_SubStep       
   
    
   if(m % (5*TDatabase::TimeDB->STEPS_PER_IMAGE) == 0)
   {      
// compute q3 for psd
#ifdef __PBS__
     if(TDatabase::ParamDB->REACTOR_P5==0)
      {
       // at present nodal functionals are just point values
       // i.e., N_U==N_PhySpacePts, so just transpose
       GetNodalPtsValuesForIntl(N_PhySpacePts, N_Intl_Levels, SolPbe, PbeValues_IntlPhys);
      
       // get the internal nodal functionals from internal nodal point values
       IntlNodal2Sol(FeSpace_Intl, PbeValues_IntlPhys, N_PhySpacePts, N_Intl_Levels, Sol_Intl);
      }
     else
      {      
       printf("Main Programm  QuadPts method not yet implemented\n");
#ifdef _MPI
       MPI_Finalize();
#endif  
       exit(0);      
      }   
      
// =============================================================================================      
// sum of all space pts for each internal pts in each time step
// but sum among all processor, only when needed
// =============================================================================================     
     memset(IntlSol_All_SUM, 0,  N_Intl_Dof*SizeOfDouble); 
     for(i=0;i<N_f_Integ_NodalPt;i++)
      {
       SolPbe_Intl_Loc = Sol_Intl + f_Integ_NodalPt[i]*N_Intl_Dof;
      
       for(j=0;j<N_Intl_Dof;j++) 
        IntlSol_All_SUM[j] += SolPbe_Intl_Loc[j];  
      } // for(i=0;i<N_f_Integ_NodalPt;i++) 
      
     
#ifdef _MPI 
    if(rank==0)
#endif 
     {
      // copy own values
      memcpy(IntlSol_All_SUM_T, IntlSol_All_SUM,  N_Intl_Dof*SizeOfDouble);
     }
          
     // get values from all other processors and add     
#ifdef _MPI
      if(rank!=0)
       { 
        MPI_Send(IntlSol_All_SUM, N_Intl_Dof, MPI_DOUBLE,  0, 100, MPI_COMM_WORLD);
       }
      else
       {
        for(i=1;i<size;i++)
         {
          MPI_Recv(IntlSol_tmp2, N_Intl_Dof, MPI_DOUBLE,  i, 100, MPI_COMM_WORLD, &status);

          for(j=0;j<N_Intl_Dof;j++)  
           IntlSol_All_SUM_T[j] += IntlSol_tmp2[j]; 
         }
       }//  else if(rank!=0)
#endif   

#ifdef _MPI
      if(rank==0)
#endif
       {       
        time ( &rawtime );
        timeinfo = localtime (&rawtime );
  
        os.seekp(std::ios::beg);
        if (N_Q3Data<10) os << "PSD_DATA/PSD_All_0000"<<N_Q3Data<<".data" << ends;
        else if (N_Q3Data<100) os << "PSD_DATA/PSD_All_000"<<N_Q3Data<<".data" << ends;
        else if (N_Q3Data<1000) os << "PSD_DATA/PSD_All_00"<<N_Q3Data<<".data" << ends;        
        else if (N_Q3Data<10000) os << "PSD_DATA/PSD_All_0"<<N_Q3Data<<".data" << ends;
        else os << "PSD_DATA//PSD_All_"<<N_Q3Data<<".data" << ends;
      
        std::ofstream dat(os.str().c_str());
        dat.setf(std::ios::fixed);
//         dat << setprecision(8);    

        if (!dat)
         {
          cerr << "cannot open file for output" << endl;
          return -1;
         }
        dat << "%% PSD data created with operator-splitting parallel FEM" << endl;
        dat << "%% Date & Time: " <<asctime (timeinfo);
        dat << "%% Current Time of the problem :" << TDatabase::TimeDB->CURRENTTIME << endl;
        dat << "%% L , PSD_All: ((0,T] X  Omega_X), q3 of PSD_All: " << endl;
        
        // calculate q3(Lmin, Lmax)
        integral = 0;
        for(j=0; j<N_Intl_Dof-1; j++)      
         {
          val_left = IntlPosL[j]*IntlPosL[j]*IntlPosL[j]*IntlSol_All_SUM_T[j];
          val_right = IntlPosL[j+1]*IntlPosL[j+1]*IntlPosL[j+1]*IntlSol_All_SUM_T[j+1];
          integral += (val_left + val_right)*(IntlPosL[j+1] - IntlPosL[j])*0.5;
         }        
        
        // calculate discrete q_3
        q3_max = -1.e10;
        for(j=0; j<N_Intl_Dof; j++)
         {
          L = IntlPosL[j];   
          IntlSol_tmp2[j] = L*L*L*IntlSol_All_SUM_T[j]/integral;
          if(q3_max<IntlSol_tmp2[j])  q3_max=IntlSol_tmp2[j];  
         }

        for(j=0; j<N_Intl_Dof; j++)
         dat << IntlPosL[j] << " " <<  IntlSol_All_SUM_T[j]*TDatabase::ParamDB->UREA_f_infty <<" " <<  IntlSol_tmp2[j]/q3_max << endl;

//         for(j=0; j<N_Intl_Dof; j++)
//          cout << "PSDALL " << IntlSol_All_SUM_T[j]*TDatabase::ParamDB->UREA_f_infty << endl;

        dat.close();
        cout << endl;
        OutPut( "Q3 data wrote into file " << N_Q3Data <<endl);
        N_Q3Data++;         
       } //   if(rank==0)

   } // if(m % 100 == 0)
 // =============================================================================================      
   
#endif // __PBS__   
   
  // measure errors to known solution
  if(TDatabase::ParamDB->MEASURE_ERRORS)
   {     
     
     for(i=0;i<N_IndepntScalarEqns;i++)
      {    
        fesp[0] = Scalar_Spaces[i];
        aux =  new TAuxParam3D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);
  
        ScalarFunctions[i]->GetErrors(Exact, 4, AllDerivatives, 2, L2H1Errors, 
                                      Coefficients[i], aux, 1, fesp, errors);
        delete aux;


#ifdef _MPI
        MPI_Reduce(errors, &l2, 1, MPI_DOUBLE, MPI_SUM, out_rank, MPI_COMM_WORLD);
        MPI_Reduce(errors+1, &H1, 1, MPI_DOUBLE, MPI_SUM, out_rank, MPI_COMM_WORLD);
        if(rank==out_rank)
         {l2 = sqrt(l2);   H1 = sqrt(H1);}
#else
        l2 = errors[0];
        H1 = errors[1];
#endif


#ifdef _MPI
        if(rank==out_rank)
#endif
        {
         OutPut("Time: " << TDatabase::TimeDB->CURRENTTIME);
         OutPut(" L2: " << l2);
         OutPut(" H1-semi: " << H1 << endl);

        if(L2error_Max<l2)
         {
          L2error_Max=l2;
          L2error_Max_t=TDatabase::TimeDB->CURRENTTIME;
         } 
 
        if(m>1)
         {
         L2T[i] += (l2*l2 +olderror[i] * olderror[i])*TDatabase::TimeDB->TIMESTEPLENGTH/2.0;
         H1T[i] += (H1*H1 +olderror1[i] * olderror1[i])*TDatabase::TimeDB->TIMESTEPLENGTH/2.0;
          
          OutPut(L2error_Max_t <<  " L2error_Max " << L2error_Max << " ");
          OutPut(TDatabase::TimeDB->CURRENTTIME <<  " L2(0,T;L2) " << sqrt(L2T[i]) << " ");
          OutPut( "L2(0,T;H1) " << sqrt(H1T[i]) << endl);
         }
         olderror[i] = l2;
         olderror1[i] =  H1; 
 
        } // if(rank==out_rank)
      } // for(i=0;i<N_IndepntS
      
#ifdef __PBS__  

//   for(i=0;i<N_PhySpacePts;i++)
//    {
//     SolPbe_Intl_Loc = PbeValues_IntlPhys + i*N_Intl_Levels;
//     ADI_System[i]->Interpolate(SolPbe_Intl_Loc, Exact_psd);
//    }    
// 
//     if(TDatabase::ParamDB->REACTOR_P5==0)
//      {
//       GetPhySolFromIntlVal_Nodal(N_Intl_Levels, N_PhySpacePts, PbeValues_IntlPhys, SolPbe);
//      }
//     else
//      {
//        printf("Main Programm GetPhySolFromIntlVal_Quad not yet implemented !!!\n");  
// #ifdef _MPI      
//     MPI_Finalize(); 
// #endif 
//        exit(0);
//      }

  for(i=0;i<7;i++)
   errors[i] = 0.;
    
   GetOSError(N_Intl_Levels, ScalarPbeFunctions,  FeSpace_Intl, SolPbe, N_U, errors);  
  
#ifdef _MPI
  MPI_Reduce(errors, &l2, 1, MPI_DOUBLE, MPI_SUM, out_rank, MPI_COMM_WORLD);
  
  if(rank==out_rank)
   l2 = sqrt(l2);   
  
#else
  l2 = errors[0];
#endif  // ifdef _MPI  


#ifdef _MPI
        if(rank==out_rank)
#endif
        {
         OutPut("Time: " << TDatabase::TimeDB->CURRENTTIME);
         OutPut(" PSD L2: " << l2 << endl); 
 
         if(L2error_Max_Pbe<l2)
          {
           L2error_Max_Pbe=l2;
           L2error_Max_t_Pbe=TDatabase::TimeDB->CURRENTTIME;
          } 

        if(m>1)
         {
          L2T_Pbe += (l2*l2 +olderror_Pbe *olderror_Pbe)*TDatabase::TimeDB->TIMESTEPLENGTH/2.0;
          
          OutPut(L2error_Max_t_Pbe <<  " L2error_Max " << L2error_Max_Pbe << " ");
          OutPut(TDatabase::TimeDB->CURRENTTIME <<  " L2(0,T;L2) " << sqrt(L2T_Pbe) << endl);
         }
         olderror_Pbe = l2;
        }  
#endif      
   } //   if(TDatabase::ParamDB->MEASURE_ERRORS)
   
   
    
//  #ifdef _MPI      
//    if(rank==0) 
// #endif       
//   printf("Main Programm  N_Intl_Levels %d \n", N_Intl_Levels );
  
// #ifdef _MPI      
//     MPI_Finalize(); 
// #endif    
//   exit(0);
//     

  if( (m % TDatabase::TimeDB->STEPS_PER_IMAGE == 0)  && TDatabase::ParamDB->WRITE_VTK)
   {
#ifdef _MPI     
       Output->Write_ParVTK(Comm, img, SubID);     
#else
        os.seekp(std::ios::beg);
        if(img<10) os << VtkBaseName<<".0000"<<img<<".vtk" << ends;
        else if(img<100) os << VtkBaseName<<".000"<<img<<".vtk" << ends;
        else if(img<1000) os << VtkBaseName<<".00"<<img<<".vtk" << ends;
        else if(img<10000) os << VtkBaseName<<".0"<<img<<".vtk" << ends;
        else  os << VtkBaseName<<"."<<img<<".vtk" << ends;
        Output->WriteVtk(os.str().c_str());
#endif
 
      img++;
    } // if(TDatabase::ParamDB->WRITE_VTK)


  if(m<3 || m % 50 == 0) 
   {
#ifdef _MPI       
    t_par2 = MPI_Wtime();
#else
    t_dt = GetTime();
#endif  
#ifdef _MPI   
    if(rank==0)
#endif       
     printf("Time taken for this time step %e\n", (t_par2-t_dt));
    }
    
    
    if(TDatabase::TimeDB->CURRENTTIME>10)
     TDatabase::TimeDB->TIMESTEPLENGTH = 0.1;

   } // while(TDatabase::TimeDB->CURRENT

  if(TDatabase::ParamDB->WRITE_VTK)
   {
 #ifdef _MPI    
//      if(N_IndepntScalarEqns)
       Output->Write_ParVTK(Comm, img, SubID);     
#else
//      if(N_IndepntScalarEqns)
//       {
        os.seekp(std::ios::beg);
        if(img<10) os << VtkBaseName<<".0000"<<img<<".vtk" << ends;
        else if(img<100) os << VtkBaseName<<".000"<<img<<".vtk" << ends;
        else if(img<1000) os << VtkBaseName<<".00"<<img<<".vtk" << ends;
        else if(img<10000) os << VtkBaseName<<".0"<<img<<".vtk" << ends;
        else  os << VtkBaseName<<"."<<img<<".vtk" << ends;
        Output->WriteVtk(os.str().c_str());
//       }      
#endif 
   } // if(TDatabase::ParamDB->WRITE_VTK)

  CloseFiles();
  
#ifdef _MPI
     MPI_Finalize();
#endif    
 
  return 0;
}



