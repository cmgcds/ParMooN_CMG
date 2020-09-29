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
#include <LocalProjection.h>

#include <sys/stat.h>
#include <sys/types.h>

// #include "../Examples/TNSE_2D/Droponsolid.h"
// #include "../Examples/TNSE_2D/Drop_imping_axial3D.h"
#include "../TNSE_2D/DropHeat_imping_axial3D_heatline.h"

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
    
//     if(y<0.05)
//      cout << " sorting " << x << ' ' << y<<endl;
    if(  sqrt((x-X0)*(x-X0) +(y-Y0)* (y-Y0))<1e-8  )
      {
       //cout << " sorting " << x << ' ' << y<<endl;
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
   } // 

  //cout<<  X0 << "SLPX X0 " << Y0 << endl;

  if(i==N)
  {
   cout<<N << "Error in finding start vert " << X0 << endl;   
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


//works only for triangles
void GetHetroFEs(int N_Cells, TCollection *coll, int *Triangles, int *PointNeib,
                 int maxEpV, FE2D *fes)
{
 int i, j, k, l, Index, N_CpV, CellIndex;
 int ORDER, LPSORDER, N_LPcells, N_Edges=0;
 
 TBaseCell *Me;
 TJoint *Joint; 
 TBoundEdge *Solid_Joint;
 TBoundComp *BoundComp; 
 FE2D FE, LPS_FE, *tmpfes;
  
 tmpfes = new FE2D[N_Cells];
 
 LPSORDER  = TDatabase::ParamDB->ANSATZ_ORDER;
 ORDER = (LPSORDER + 1)/101;

  switch(ORDER)
   {
    case 1:
      FE = C_P1_2D_T_A;   
      break;
    case 2:
      FE = C_P2_2D_T_A;   
      break;
    case 3:
      FE = C_P3_2D_T_A;   
      break; 
    case 4:
      FE = C_P4_2D_T_A;   
      break;       
    case 5:
      FE = C_P5_2D_T_A;   
      break;        
    default: cerr << "unknown order GetHetroFEs" << endl;
             exit(-1);
             break;    
   } // switch(ORDER)
    
   switch(LPSORDER)
   {
    case 100:
      LPS_FE = C_UL1_2D_T_A;   
      break;
    case 201:
      LPS_FE = C_UL2_2D_T_A;   
      break;
    case 302:
      LPS_FE = C_UL3_2D_T_A;   
      break; 
    case 403:
      LPS_FE = C_UL4_2D_T_A;   
      break;       
    case 504:
      LPS_FE = C_UL5_2D_T_A;   
      break;        
    default: cerr << "unknown LPS order GetHetroFEs" << endl;
             exit(-1);
             break;    
   } // switch(ORDER)
 
  for(i=0; i<N_Cells; i++) 
   fes[i] = FE;  
  
  for(i=0; i<N_Cells; i++)
  {
   Me = coll->GetCell(i);

   for(l=0;l<3;l++)
   {
     Joint = Me->GetJoint(l);
     
     //free surface and interface  
     if(Joint->GetType()==BoundaryEdge || Joint->GetType()==InterfaceJoint
       || Joint->GetType() == IsoBoundEdge)
     {
      Solid_Joint = (TBoundEdge *)Joint;
      BoundComp = Solid_Joint->GetBoundComp();
      
      //if(BoundComp->GetID()==0 ||   Joint->GetType() == IsoBoundEdge)
     if(BoundComp->GetID()==0)
      {
       //add all cells associated with this edge vertices
       //start vertex
       //cout<< " N_Edges " << N_Edges++ <<endl;
       Index = Triangles[i*3 + l];
       N_CpV = PointNeib[Index*maxEpV];
       
       for(j=1;j<=N_CpV;j++)
        {
         CellIndex = PointNeib[Index*maxEpV +j];
         fes[CellIndex] = LPS_FE;
        }
       
      }//   if(BoundComp->GetID()==0 |    
     }// if(Joint->GetType()==Bound     
   } //  for(l=0;l<k;l++)
  }   // for(i=0; i<N_Cells; i++)  


  for(k=0; k<1; k++)  
  {
   for(i=0; i<N_Cells; i++)  
    tmpfes[i] = fes[i];
   
   for(i=0; i<N_Cells; i++)  
   {
    if( (tmpfes[i]==C_P1_2D_T_A) || (tmpfes[i]==C_P2_2D_T_A)  || (tmpfes[i]==C_P3_2D_T_A)
        || (tmpfes[i]==C_P4_2D_T_A) || (tmpfes[i]==C_P5_2D_T_A) )
       continue;    
    
    for(l=0;l<3;l++)
    {
     Index = Triangles[i*3 + l];
     N_CpV = PointNeib[Index*maxEpV]; 
    
     for(j=1;j<=N_CpV;j++)
     {
      CellIndex = PointNeib[Index*maxEpV +j];
      fes[CellIndex] = LPS_FE;
     }     
    } // for(l=0;l<3;l++)    
   } // for(i=0; i<N_Cells; i++)  
  }//levels

  N_LPcells = 0;
   for(i=0; i<N_Cells; i++)  
   {
    if( (fes[i]==C_P1_2D_T_A) || (fes[i]==C_P2_2D_T_A)  || (fes[i]==C_P3_2D_T_A)
        || (fes[i]==C_P4_2D_T_A) || (fes[i]==C_P5_2D_T_A) )
       continue;    
       
    N_LPcells++;
   }
   
 OutPut("HetroFEs N_LPcells : "<< N_LPcells <<endl); 
 
 
delete [] tmpfes;

/*  cout<< "test GetHetroFEs " <<endl;
  exit(0);*/ 
}// GetHetroFEs

void BDChangeFERemap_All(TFESpace2D ** &FESpaces_All, TFEVectFunct2D ** &FEVectFuncts_All, TFEFunction2D ** &FEFunctions_All,
                         int *N_GidDofs,  int *N_GridActive, double ** &Sol_All, double ** &Rhs_All, int *Triangles, int *PointNeighb, int maxEpV)
{
 int i, j, k, N_U_output, N_G, N_Cells, pressure_space_code; 
  
  
 TCollection *coll, *mortarcoll = NULL; 
 TFESpace2D *NewFESpace;
 TFEFunction2D *NewFEFunction;
 TBaseCell *cell;
 FE2D FEId, NewFEId;
 TFEDesc2D *FeDesc, *NewFeDesc;   
 TFE2D *Element;
 TNodalFunctional2D *nf;
 TRefTrans2D *rt;
 FE2D *fes;
 TFESpace2D *pressure_space_output;
 
  char ReadinDat[] = "readin.dat";
  char TString[] = "T";
  char NameString[]  = "name";
  char UString[] = "U";
  char PString[] = "P";
  char WString[] = "W";
  
  // assumed that both old and new fe spaces have same coll but different no. edges 
  // due to freebd joint becomes interface joint
  coll =  FESpaces_All[4]->GetCollection();
  N_Cells = coll->GetN_Cells();
 
//   cout << "test BDChangeFERemap_All "  << endl;
    // mesh velocity space 
  delete FESpaces_All[2];
  
  FESpaces_All[2] = new TFESpace2D(coll, NameString, WString, GridBoundCondition, 1, NULL);
  N_GidDofs[0] = FESpaces_All[2]->GetN_DegreesOfFreedom();;
  N_G =   N_GidDofs[0];
  
  
  delete [] Sol_All[1];
  Sol_All[1] = new double[2*N_G];
  delete [] Rhs_All[1]; 
  Rhs_All[1] = new double[2*N_G];  
  
  memset(Sol_All[1], 0, 2*N_G*SizeOfDouble); 
  delete FEVectFuncts_All[1];
  FEVectFuncts_All[1]  = new TFEVectFunct2D(FESpaces_All[2], WString, WString, Sol_All[1], N_G, 2);
  
  FEFunctions_All[3] = FEVectFuncts_All[1]->GetComponent(0);
  FEFunctions_All[4] = FEVectFuncts_All[1]->GetComponent(1);   
  
  
  delete FESpaces_All[4];
  GetVelocityAndPressureSpace(coll,BoundCondition_output,
                              mortarcoll, FESpaces_All[4],
                              pressure_space_output, &pressure_space_code,
                              TDatabase::ParamDB->VELOCITY_SPACE,
                              TDatabase::ParamDB->PRESSURE_SPACE);
   N_U_output = FESpaces_All[4]->GetN_DegreesOfFreedom();  
   
  // velo in all domains
  delete [] Sol_All[3]; 
  delete FEVectFuncts_All[2];  
  
  Sol_All[3] = new double[2*N_U_output]; 
  memset(Sol_All[3], 0, 2*N_U_output*SizeOfDouble);
  FEVectFuncts_All[2] =  new TFEVectFunct2D(FESpaces_All[4], UString, UString, Sol_All[3], N_U_output, 2);   
//     cout << "test BDChangeFERemap_All "  << endl;
}



void BDChangeFERemap(TFESpace2D *&FESpace, TFEFunction2D *&FEFunction, double *&Sol,
                     TSquareStructure2D *&SqStruct, TSquareMatrix2D *&M, TSquareMatrix2D *&A,
                     int *Triangles, int *PointNeighb, int maxEpV)
{
 int i, j, k, N_DOF, N_Cells, N_LocDOF, N_NewLocDOF;
 int *BeginIndex, *GlobalNumbers, *NewBeginIndex, *NewGlobalNumbers;
 int *DOF, *NewDOF, *IncidentArray, N_Points;
 
 double *NewSol, *xi, *eta; 
 double X[MaxN_PointsForNodal2D], Y[MaxN_PointsForNodal2D];
 double AbsDetjk[MaxN_PointsForNodal2D];
 double values[4];
 double PointValues[MaxN_PointsForNodal2D], FunctionalValues[MaxN_PointsForNodal2D];
  
 TCollection *coll; 
 TFESpace2D *NewFESpace;
 TFEFunction2D *NewFEFunction;
 TBaseCell *cell;
 FE2D FEId, NewFEId;
 TFEDesc2D *FeDesc, *NewFeDesc;   
 TFE2D *Element;
 TNodalFunctional2D *nf;
 TRefTrans2D *rt;
 FE2D *fes;
  
  char TString[] = "T";
  char NameString[]  = "name"; 
 
  // assumed that both old and new fe spaces have same coll but different no. edges 
  // due to freebd joint becomes interface joint
  coll =  FESpace->GetCollection();
  N_Cells = coll->GetN_Cells();
 
  
   if(TDatabase::ParamDB->ANSATZ_ORDER<100)
   {       
    NewFESpace = new TFESpace2D(coll, NameString, TString, HeatBoundCondition,
                                TDatabase::ParamDB->ANSATZ_ORDER, NULL);
   }
   else
   {
     
    fes = new  FE2D[N_Cells];  
    GetHetroFEs(N_Cells, coll, Triangles, PointNeighb, maxEpV, fes);   
    
    NewFESpace = new TFESpace2D(coll, NameString, TString, HeatBoundCondition,
                                fes, NULL);   
    delete [] fes;
   }
   
  N_DOF = NewFESpace->GetN_DegreesOfFreedom();
  //cout<<"BDChangeFERemap N_DOF: "<< N_DOF<<endl; 
  
  IncidentArray = new int[N_DOF];
  NewSol = new double[N_DOF];

  memset(IncidentArray, 0,  N_DOF*SizeOfInt); 
  memset(NewSol, 0,  N_DOF*SizeOfDouble);     
  
  NewFEFunction = new TFEFunction2D(NewFESpace, TString, TString, NewSol, N_DOF);
 
  GlobalNumbers = FESpace->GetGlobalNumbers();
  BeginIndex = FESpace->GetBeginIndex();  
 
  NewGlobalNumbers = NewFESpace->GetGlobalNumbers();
  NewBeginIndex = NewFESpace->GetBeginIndex(); 
     
  //map old sol to new sol
  for(i=0; i<N_Cells; i++)
   {
    cell = coll->GetCell(i);

    NewFEId = NewFESpace->GetFE2D(i, cell);   
    NewFeDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(NewFEId);
    N_NewLocDOF = NewFeDesc->GetN_DOF();  
    NewDOF = NewGlobalNumbers + NewBeginIndex[i];    
    
    FEId = FESpace->GetFE2D(i, cell);   
    FeDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(FEId);
    N_LocDOF = FeDesc->GetN_DOF();      
    DOF = GlobalNumbers + BeginIndex[i];

    if(N_NewLocDOF == N_LocDOF ||  N_NewLocDOF<N_LocDOF)//Galerlin-Galerlin or LPS-Galerlin
    {
     for(j=0;j<N_NewLocDOF;j++)
     {
      NewSol[NewDOF[j]] +=  Sol[DOF[j]];
      IncidentArray[NewDOF[j]]++;
     }
    }
    else if(N_NewLocDOF>N_LocDOF) //Galerlin to LPS
    {  
     Element = TFEDatabase2D::GetFE2D(NewFEId); 
     nf = Element->GetNodalFunctional2D();
     nf->GetPointsForAll(N_Points, xi, eta);      

     rt = TFEDatabase2D::GetRefTrans2D(TriaAffin);
     ((TTriaAffin *)rt)->SetCell(cell);    
     
     TFEDatabase2D::GetOrigFromRef(TriaAffin, N_Points, xi, eta, X, Y, AbsDetjk);
     for(j=0;j<N_Points;j++)
     {
      FEFunction->FindGradientLocal(cell, i, X[j], Y[j], values);
      PointValues[j] = values[0];
     } 

     nf->GetAllFunctionals(coll, (TGridCell *)cell, PointValues, FunctionalValues);

     for(j=0;j<N_NewLocDOF;j++)
     {
      //copy old values for dofs on edges
      if(j<N_LocDOF)
      {
       NewSol[NewDOF[j]] +=  Sol[DOF[j]];
       IncidentArray[NewDOF[j]]++;
      }
      else  //bubbles are inner odfs, numbered after dof on edges
      {
       NewSol[NewDOF[j]] +=  FunctionalValues[j];
       IncidentArray[NewDOF[j]]++;
       
       // cout << NewDOF[j] << " Test  BDChangeFERemap " << NewSol[NewDOF[j]]<< endl; 
      }     
     } //  for(j=0;j<N_NewLocDO

    }
/*    else if(N_NewLocDOF<N_LocDOF) 
    {
      
      
     OutPut("Cells from LPS to Galerlin not yet implemented "<<endl);
     exit(0);
    }   */ 
    
   }// for(i=0; i<N_Cells; i++)
   
  for(i=0; i<N_DOF; i++)
   {
    if(IncidentArray[i]==0)
    {
     cout << i<< "Error in BDChangeFERemap, IncidentArray[i]==0" << endl; 
     exit(0);         
    }
    
    NewSol[i] /=(double)IncidentArray[i];    
   }

  delete FESpace;  
  delete [] Sol;
  delete FEFunction;  
  
  FESpace = NewFESpace;
  Sol = NewSol;
  FEFunction = NewFEFunction;
  
  delete SqStruct;
  SqStruct = new TSquareStructure2D(NewFESpace);
  SqStruct->Sort();   

  delete M;
  M = new TSquareMatrix2D(SqStruct);  

  delete A;
  A = new TSquareMatrix2D(SqStruct);    

  delete [] IncidentArray;
  
//  cout << "test  BDChangeFERemap" << endl; 
// exit(0); 
}// BDChangeFERemap


//Heat function
void BDChangeFERemap(TFESpace2D *&FESpace, TFEFunction2D *&FEFunction, double *&Sol,
                     TSquareStructure2D *&SqStruct, TSquareMatrix2D *&A,
                     int *Triangles, int *PointNeighb, int maxEpV)
{
 int i, j, k, N_DOF, N_Cells;
 
 TCollection *coll; 
 FE2D *fes;
  
  char TString[] = "H";
  char NameString[]  = "name"; 
 
  // assumed that both old and new fe spaces have same coll but different no. edges 
  // due to freebd joint becomes interface joint
  coll =  FESpace->GetCollection();
  N_Cells = coll->GetN_Cells();
 
  delete FESpace;  
  delete [] Sol;
  delete FEFunction;  
  
   if(TDatabase::ParamDB->ANSATZ_ORDER<100)
   {       
    FESpace = new TFESpace2D(coll, NameString, TString, HeatFuncBoundCondition,
                                TDatabase::ParamDB->ANSATZ_ORDER, NULL);
   }
   else
   {
     
    fes = new  FE2D[N_Cells];  
    GetHetroFEs(N_Cells, coll, Triangles, PointNeighb, maxEpV, fes);   
    
    FESpace = new TFESpace2D(coll, NameString, TString, HeatFuncBoundCondition,
                                fes, NULL);   
    delete [] fes;
   }
   
  N_DOF = FESpace->GetN_DegreesOfFreedom();
  //cout<<"BDChangeFERemap N_DOF: "<< N_DOF<<endl; 
  
  Sol = new double[N_DOF];
  memset(Sol, 0,  N_DOF*SizeOfDouble);  
  
  FEFunction = new TFEFunction2D(FESpace, TString, TString, Sol, N_DOF);

  delete SqStruct;
  SqStruct = new TSquareStructure2D(FESpace);
  SqStruct->Sort();   

  delete A;
  A = new TSquareMatrix2D(SqStruct);    

//  cout << "test  BDChangeFERemap" << endl; 
// exit(0); 
}// BDChangeFERemap

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

void GetHeatConvectionALEVect(TFEVectFunct2D **Velocity_All, 
            TFEVectFunct2D *GridVect_NSE, TFEVectFunct2D *GridVect_S)
{
 int i, j, k, l, N_Cells_NSE, N_Cells_S; 
 int n1, n2, N, N_DOF, *DOF, *DOF_All;
 int *UGlobalNumbers, *UBeginIndex, *UGlobalNumbers_All, *UBeginIndex_All;
 int *WGlobalNumbers_NSE, *WBeginIndex_NSE, *WGlobalNumbers, *WBeginIndex;
 int *WGlobalNumbers_S, *WBeginIndex_S;
 
 double *U1, *U2, *U1_All, *U2_All, *W1_NSE, *W2_NSE, *W1_S, *W2_S, *Velo1, *Velo2;
 
 TFESpace2D *VelocitySpace, *GridSpace, *GridSpace_NSE, *GridSpace_S;
 TCollection  *Coll_NSE, *Coll_S;
 TBaseCell *cell; 
 FE2D FeId;
 TFEDesc2D *FeDesc;  
 
  // velocity
  VelocitySpace = Velocity_All[0]->GetFESpace2D();
  U1 = Velocity_All[0]->GetValues();
  U2 = U1 + Velocity_All[0]->GetLength(); 
  UGlobalNumbers =  VelocitySpace->GetGlobalNumbers();
  UBeginIndex =  VelocitySpace->GetBeginIndex(); 
  
  U1_All = Velocity_All[2]->GetValues();
  U2_All = U1_All + Velocity_All[2]->GetLength();     
  UGlobalNumbers_All =  (Velocity_All[2]->GetFESpace2D())->GetGlobalNumbers();
  UBeginIndex_All =  (Velocity_All[2]->GetFESpace2D())->GetBeginIndex(); 
  
  //grid velo  
  GridSpace = Velocity_All[1]->GetFESpace2D();
  Velo1 = Velocity_All[1]->GetValues();
  Velo2 = Velo1 + Velocity_All[1]->GetLength(); 
  WGlobalNumbers = GridSpace->GetGlobalNumbers();
  WBeginIndex = GridSpace->GetBeginIndex();   
 
  GridSpace_NSE = GridVect_NSE->GetFESpace2D();
  Coll_NSE = GridSpace_NSE->GetCollection();  
  N_Cells_NSE = Coll_NSE->GetN_Cells();;
  W1_NSE = GridVect_NSE->GetValues();
  W2_NSE = W1_NSE + GridVect_NSE->GetLength(); 
  WGlobalNumbers_NSE = GridSpace_NSE->GetGlobalNumbers();
  WBeginIndex_NSE = GridSpace_NSE->GetBeginIndex();  
  
  
  GridSpace_S = GridVect_S->GetFESpace2D();
  Coll_S = GridSpace_S->GetCollection();    
  N_Cells_S = Coll_S->GetN_Cells();;
  W1_S = GridVect_S->GetValues();
  W2_S = W1_S + GridVect_S->GetLength();  
  WGlobalNumbers_S = GridSpace_S->GetGlobalNumbers();
  WBeginIndex_S = GridSpace_S->GetBeginIndex();   

  
  memset(U1_All, 0, (2*Velocity_All[2]->GetLength())*SizeOfDouble); 
  memset(Velo1, 0, (2*Velocity_All[1]->GetLength())*SizeOfDouble); 
  
  for(i=0;i<N_Cells_S; i++)
   {
    cell = Coll_S->GetCell(i);
    N = cell->GetGlobalCellNo();
    
   //grid velocity 
    FeId = GridSpace_S->GetFE2D(i, cell);
    FeDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(FeId);
    N_DOF = FeDesc->GetN_DOF();   

    DOF = WGlobalNumbers_S + WBeginIndex_S[i];
    DOF_All = WGlobalNumbers + WBeginIndex[N];  
        
    for(j=0;j<N_DOF;j++)
     {
      n1 = DOF[j];
      n2 = DOF_All[j];
      Velo1[n2] = W1_S[n1];
      Velo2[n2] = W2_S[n1];      
     }   
   } // 
  
//   cout << " test "<< endl;
  
  for(i=0;i<N_Cells_NSE; i++)
   {
    cell = Coll_NSE->GetCell(i);
    N = cell->GetGlobalCellNo();
 
     //nse velocity
    FeId = VelocitySpace->GetFE2D(i, cell);
    FeDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(FeId);
    N_DOF = FeDesc->GetN_DOF();   
    
    DOF = UGlobalNumbers + UBeginIndex[i];
    DOF_All = UGlobalNumbers_All + UBeginIndex_All[N]; 
    
    for(j=0;j<N_DOF;j++)
     {
      n1 = DOF[j];
      n2 = DOF_All[j];
      U1_All[n2] = U1[n1];
      U2_All[n2] = U2[n1];      
     }   
     
   //grid velocity 
    FeId = GridSpace_NSE->GetFE2D(i, cell);
    FeDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(FeId);
    N_DOF = FeDesc->GetN_DOF();   

    DOF = WGlobalNumbers_NSE + WBeginIndex_NSE[i];
    DOF_All = WGlobalNumbers + WBeginIndex[N];  
        
    for(j=0;j<N_DOF;j++)
     {
      n1 = DOF[j];
      n2 = DOF_All[j];
      Velo1[n2] = W1_NSE[n1];
      Velo2[n2] = W2_NSE[n1];      
     }   
   }//for(i=0;i<N_Cells_NSE; i++)
   
//  exit(0);
  
} // GetHeatConvectionALEVect

void GridVelo_imping(double **Entries, double *Sol, double *d, double *Rhs,
                     int *KCol, int *RowPtr,
                     TFEVectFunct2D *GridPos,
                     TFEVectFunct2D *AuxGridPos,
                     TFEVectFunct2D *Velocity, double dt,
                     TFEVectFunct2D *GridVelocity, 
                     TVertex ***MovBoundVert, int *N_MovVert,
                     TBaseCell **Free_Cells, int **IsoCellEdgeNos,
                     bool &reparam, TFEVectFunct2D *RefGridPos,
                     double **Entries_S, double *Sol_S, double *Rhs_S,
                     int *KCol_S, int *RowPtr_S,
                     TFEVectFunct2D *GridPos_S, 
                     TFEVectFunct2D *RefGridPos_S, 
                     TFEVectFunct2D *GridVelocity_S)
{
  int i,j,k,l,m,comp;
  int *VeloBeginIndex, *VeloGlobalNumbers;
  int *GridBeginIndex, *GridGlobalNumbers;
  int N_Cells, N_Vertices, N_Edges, N_LocalDOFs, N_BoundaryNodes_S;
  int N_Levels, *DOF, *JointDOF, GridLength, N_BoundaryNodes;
  int N_LinePoints, IIso, N, N_Inner, N_, Max_GridLength, GridLength_S; 
  
  double xi[4], eta[4], X[4], Y[4], VX[4], VY[4];
  double FunctValues[4][MaxN_BaseFunctions2D];
  double FEValuesX[MaxN_BaseFunctions2D];
  double FEValuesY[MaxN_BaseFunctions2D];
  double *ValuesX, *ValuesY, Rx, Ay;
  double *ValuesVX, *ValuesVY, *NewValuesX, *NewValuesY;
  double s, t, x, y, h_tot, x0, x1, y0, y1, res, oldres;
  double *gridvelo, *Nx, *Ny, *LineWeights, *zeta;
  double normalx, normaly, tangenx, tangeny, nx, ny, tx, ty;
  double un, hE, t0,t1, temp2, eps=1e-6;
  double h, hmin, hmax, hlimit, *IsoX, *IsoY;
  double *RefValueX, *RefValueY;
  
  boolean OnBoundary;
  
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
  TIsoBoundEdge *isojoint;
  TJoint *joint;
  TVertex **Vertices;
  TQuadFormula1D *qf1;
  QuadFormula1D LineQuadFormula; 
  BF2DRefElements RefElement;
  TRefTrans2D *F_K;
  RefTrans2D RefTrans;
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

  GridLength_S = GridPos_S->GetLength();  
  N_BoundaryNodes_S = GridLength_S - (GridPos_S->GetFESpace2D())->GetN_Inner();
  
  Max_GridLength = GridLength_S;
  
  if(Max_GridLength<GridLength)
   Max_GridLength = GridLength;  
  
//   d = new double[2*Max_GridLength];
  
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
             { NewValuesY[l] = ValuesY[l];              
              // cout << " Solid BD " << l << " " << NewValuesX[l] << endl;
             }
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

   // store the original position  
   GridPos_S->GridToData();  

   //find the disp of wettinf point and move the remaining solid pts accordingly
   MovBoundVert[2][0]->GetCoords(Rx, y);
   MovBoundVert[1][0]->GetCoords(x, Ay);   
   
   // now move the BDs of the liquid domain
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
                      Velocity, NULL, FALSE);  
     
    reparam = TRUE;   
    RefGridPos->GridToData();    
    Daxpy(2*GridLength, -1, ValuesX, RefValueX); // now reparam disp in RefValueX

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
    }// for(i=0;i<N_MovVert[2];i++)
    
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
   //shift all pts according to the top pt disp
   
   MovBoundVert[1][0]->GetCoords(x, Ay);   
   MovBoundVert[0][0]->GetCoords(x, y);   
   N=N_MovVert[1] ;      
   h_tot = (y - Ay)/(double)N; 
   
   for(i=1;i<N;i++)
    {
     //MovBoundVert[1][i]->GetCoords(x, y);
     // cout<< " y " << y <<" new y " <<  Ay + h_tot*(double)i<<endl; 
     y  = Ay + h_tot*(double)i; 
     MovBoundVert[1][i]->SetCoords(x, y);   
    }  
    
//     exit(0);
   
   
/*   MovBoundVert[1][0]->GetCoords(x, y);
   
   N=N_MovVert[1];      
   h_tot = (y-Ay)/(double)N;   
   
   N--;
   
    for(i=1;i<N;i++)
    {
     MovBoundVert[1][i]->GetCoords(x, y);
     //cout<< " y " << y <<" new y " << y +((double)(N-i))*h_tot<<endl; 
     y += ((double)(N-i))*h_tot;
 
     MovBoundVert[1][i]->SetCoords(x, y);   
    }    */ 
//    exit(0);
   
//    h_tot = -y;
//    h_tot /= (double)N_MovVert[1];
//    for(i=1;i<N_MovVert[1];i++)
//     MovBoundVert[1][i]->SetCoords(x,  y +  h_tot*(double)i );
 
    // reparam the remaining solid vert and move the solid domain
    MovBoundVert[3][0]->GetCoords(x0, y0);    
    MovBoundVert[2][0]->GetCoords(x, y); // right wetting point

    MovBoundVert[2][0]->GetCoords(x, y);   
   
    N=N_MovVert[3]-1;      
    h_tot = (x-Rx)/(double)N;
    
    for(i=1;i<N;i++)
    {
     MovBoundVert[3][N_MovVert[3]-1 - i]->GetCoords(x, y);
     x += ((double)(N-i))*h_tot;
     MovBoundVert[3][N_MovVert[3]-1 - i]->SetCoords(x, y);   
    }  
    
   //get the new position of the grid BDs
   AuxGridPos->GridToData();
   RefGridPos_S->GridToData();
   
   // move back the grid BDs to the original position   
   GridPos->DataToGrid();
   GridPos_S->DataToGrid();
    
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
   
  //only inner mesh velo, since U will be interpolated in free surf reparm while move mesh     
  Daxpy(GridLength-N_BoundaryNodes, 1., Sol, gridvelo);
  Daxpy(GridLength-N_BoundaryNodes, 1., Sol+GridLength, gridvelo+GridLength);      
    
  }// if(reparam)   
   
  // for solid domain 
  ValuesX = GridPos_S->GetValues();
  ValuesY = ValuesX + GridLength_S;  
  
  NewValuesX = RefGridPos_S->GetValues();
  NewValuesY = NewValuesX + GridLength_S;  
  
  
   memset(Rhs_S, 0, 2*GridLength_S*SizeOfDouble);

   memcpy(d, NewValuesX, 2*GridLength_S*SizeOfDouble);
   Daxpy(2*GridLength_S, -1, ValuesX, d);
     
   memcpy(Rhs_S + (GridLength_S-N_BoundaryNodes_S), d+(GridLength_S-N_BoundaryNodes_S),
          N_BoundaryNodes_S*SizeOfDouble);
   memcpy(Rhs_S + (2*GridLength_S-N_BoundaryNodes_S), d+(2*GridLength_S-N_BoundaryNodes_S),
          N_BoundaryNodes_S*SizeOfDouble);

  memset(Sol_S, 0 , 2*GridLength_S*SizeOfDouble);
  memcpy(Sol_S + (GridLength_S-N_BoundaryNodes_S), d+(GridLength_S-N_BoundaryNodes_S),
         N_BoundaryNodes_S*SizeOfDouble);
  memcpy(Sol_S + (2*GridLength_S-N_BoundaryNodes_S), d+(2*GridLength_S-N_BoundaryNodes_S),
         N_BoundaryNodes_S*SizeOfDouble);  
  
 
   SolveGridEquation(Entries_S, Sol_S, Rhs_S,
                     KCol_S, RowPtr_S, GridLength_S);

   gridvelo = GridVelocity_S ->GetValues();
   memcpy(gridvelo, Sol_S, 2*GridLength_S*SizeOfDouble);
   Dscal(2*GridLength_S, 1./dt, gridvelo);	 
  
//     for(i=0;i<GridLength_S;i++)
//      cout<< i <<"  ---  "<< gridvelo[i] << "  ---  " << gridvelo[i+GridLength_S] << endl;

// cout << " GridVelo_imping " <<endl;
// exit(0);

//   delete [] d;
  
  if(TDatabase::ParamDB->P5 > 0)
   { delete [] Nx;  delete [] Ny;}
   
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
        cout <<"X[k] negative in Get_KE change Quad rule " 
             <<  X[k] << " " << Y[k] << endl;
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
                  bool &reparam, int &N_ReParam, TFEFunction2D *Energy)
{
  int i,j,k,l,m,comp, N;
  int *VeloBeginIndex, *VeloGlobalNumbers;
  int *GridBeginIndex, *GridGlobalNumbers;
  int N_Cells, N_Vertices, N_Edges, N_LocalDOFs;
  int N_Levels;
  TMGLevel2D *Level;
  int *DOF, *JointDOF;
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
  double xi[4], eta[4], X[4], Y[4], VX[4], VY[4];
  double FunctValues[4][MaxN_BaseFunctions2D];
  double FEValuesX[MaxN_BaseFunctions2D];
  double FEValuesY[MaxN_BaseFunctions2D];
  double *ValuesX, *ValuesY;
  double *ValuesVX, *ValuesVY;
  double *NewValuesX, *NewValuesY;
  double s, t, x, y, IsoX, IsoY;
  double x0, x1, y0, y1, h_tot;
  int GridLength, polydegree;
  int N_BoundaryNodes;
  TIsoBoundEdge *isojoint;
  double res, oldres;
  TJoint *joint;
  TVertex **Vertices;
  double *gridvelo, *Nx, *Ny, Ay;
  double *LineWeights, *zeta;
  int N_LinePoints;
  TQuadFormula1D *qf1;
  QuadFormula1D LineQuadFormula;
  double normalx, normaly, tangenx, tangeny, nx, ny, tx, ty;
  int IIso;
  BF2DRefElements RefElement;
  TRefTrans2D *F_K;
  RefTrans2D RefTrans;
  int N_Inner, N_;
  double un, hE;
  double t0,t1, temp2;
  BoundCond Cond0, Cond1;
  double x_max, x_min, temp, eps=1e-6;
  QuadFormula2D QuadFormula;

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
          if((TDatabase::ParamDB->P5 > 0) && (ValuesY[l] != 0) && (ValuesX[l] > 0)  )
          {
            un = VX[j]*Nx[k] + VY[j]*Ny[k];
            NewValuesX[l] = ValuesX[l] + dt*un*Nx[k];
            NewValuesY[l] = ValuesY[l] + dt*un*Ny[k];
          }
          else
          {
// 	   if(NewValuesX[l] != 0)
//                 NewValuesX[l] = ValuesX[l] + dt*VX[j];
//             if(ValuesY[l] != 0)
//                NewValuesY[l] = ValuesY[l] + dt*VY[j];
	    if(ValuesX[l] > 0 )
	    { NewValuesX[l] = ValuesX[l] + dt*VX[j]; }
	    else
	    { NewValuesX[l] = ValuesX[l] ; }

            if(ValuesY[l] == 0) 
	     { NewValuesY[l] = ValuesY[l];  }
            else    
	     { NewValuesY[l] = ValuesY[l] + dt*VY[j];  }                
          }
        }
         //Due to spline approximation solid boundary end vertices may take negative y value
         if(NewValuesY[l]<0.0 ) NewValuesY[l] = 0.0;
         if(NewValuesX[l]< 0.0 ) NewValuesX[l] = 0.0;
      } // endfor j
     } // endif
   } // endfor i

/*

   for(i=GridLength-N_BoundaryNodes;i<GridLength;i++)
 cout << i <<"  ---  " <<NewValuesX[i] << "  ---  " << NewValuesY[i]<<endl;

exit(0);
*/

   //solid BD
   NewGridPos->DataToGrid();
   
   MovBoundVert[2][0]->GetCoords(x, y); // right wetting point
   y=0.;
   h_tot = x;
   h_tot /= (double)N_MovVert[0];
   for(i=1;i<N_MovVert[0];i++)
     MovBoundVert[0][i]->SetCoords(h_tot*(double)i, y);
 
  
  // axial boundary
   MovBoundVert[1][0]->GetCoords(x, Ay);   
   MovBoundVert[0][0]->GetCoords(x, y);   
   N=N_MovVert[1] ;      
   h_tot = (y - Ay)/(double)N; 
   x = 0.;
   for(i=1;i<N;i++)
    {
     //MovBoundVert[1][i]->GetCoords(x, y);
     y  = Ay + h_tot*(double)i; 
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
    
    ReParam_axial3D_U(N_MovVert[2], Free_Cells,  IsoCellEdgeNos[1], IsoCellEdgeNos[0],
                      Velocity, Energy, TRUE);   

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
 //cout << i << "  ---  "<<Sol[i] << "  ---  " << Sol[i+GridLength] << endl;
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


// ====================================================================
// modify matrices due to integrals on free surface for heat equation
// ====================================================================
void Heat_freeint_axial3D(TSquareMatrix2D *T11, BoundCondFunct2D
                          *HeatBoundaryCondition, double *rhs)
 {
  int i, j, k, l, DOF_R, DOF_L, m;
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
  FE2D FEId;
  TFE2D *ele;
  RefTrans2D RefTrans;
  TRefTrans2D *F_K;
  TFESpace2D *fespace;
  BF2DRefElements RefElement;
  QuadFormula1D LineQuadFormula;
  TQuadFormula1D *qf1;
  int N_LinePoints, N_Active;
  double *LineWeights, *zeta;
  double  x1, y1,tx,ty,mod_t, x, y;
  int N_BaseFunct, *N_BaseFuncts;
  double **uref, **uxiref, **uetaref;
  double uorig[MaxN_BaseFunctions2D], uxorig[MaxN_BaseFunctions2D];
  double uyorig[MaxN_BaseFunctions2D];
  BaseFunct2D *BaseFuncts;
  double r2, r;
  int *KCol, *RowPtr, *JointDOF, N_DOF;
  double *ValuesT11;
  int *BeginIndex, *GlobalNumbers, *DOF, TestDOF, AnsatzDOF;
  int index1, index2;
  double val, theta, factor1, factor2, angle,  X_B[100], Y_B[100], r_axial, Mult;
  int count=0, count1=0, count2=0;

// heat stress number (or) BIOT NUMBER
  double C0 = TDatabase::ParamDB->BI_NR;
  C0 /=TDatabase::ParamDB->PE_NR;

  TFEDesc2D *FeDesc;
  BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
  N_BaseFuncts = TFEDatabase2D::GetN_BaseFunctFromFE2D();


  fespace = T11->GetFESpace();
  Coll = fespace->GetCollection();
  N_Cells = Coll->GetN_Cells();
  N_Active = fespace->GetActiveBound();
  BeginIndex = fespace->GetBeginIndex();
  GlobalNumbers = fespace->GetGlobalNumbers();

  RowPtr = T11->GetRowPtr();
  KCol = T11->GetKCol();

  ValuesT11 = T11->GetEntries();

  for(i=0;i<N_Cells;i++)
   {
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
        HeatBoundaryCondition(comp, t0, Cond0);
        HeatBoundaryCondition(comp, t1, Cond1);

        if(Cond0 == FREESURF)
        { 
//           cout << "comp " << comp  <<endl;
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
       IJoint = JointNumbers[j];
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

       for(k=0;k<N_LinePoints;k++)
        {
          F_K->GetTangent(IJoint, zeta[k], t0, t1);
          normn = sqrt(t0*t0+t1*t1);
          r_axial = fabs(X_B[k]);

          if (fabs(r_axial)<1e-10)
            OutPut("r_axial is zero check heat freesurfint" <<endl);

          Mult = C0*r_axial*normn*LineWeights[k];

          switch(RefElement)
          {
            case BFUnitSquare:
              ((TQuadIsoparametric *)F_K)->GetOrigValues(IJoint,  zeta[k],
                        N_BaseFunct, uref[k], uxiref[k], uetaref[k],
                        uorig, uxorig, uyorig);
            break;

            case BFUnitTriangle:
              ((TTriaIsoparametric *)F_K)->GetOrigValues(IJoint,  zeta[k],
                        N_BaseFunct, uref[k], uxiref[k], uetaref[k],
                        uorig, uxorig, uyorig);
            break;
          } // endswitch

          for(l=0;l<N_BaseFunct;l++)
          {
            TestDOF = DOF[l];
            index2 = RowPtr[TestDOF+1];
             for(m=0;m<N_BaseFunct;m++)
              {
               AnsatzDOF = DOF[m];
               index1 = RowPtr[TestDOF];
               if(index1+1 == index2) continue;
               while(KCol[index1] != AnsatzDOF) index1++;

               val = uorig[m]*uorig[l];
               ValuesT11[index1] += val*Mult;
              } // endfor m

          } // endfor l
        } // endfor k
      } // endfor j
    } // end (N_IsoJoints > 0)
    else
    {
      //cout << "Cell " << i << " has NO free surface." << endl;
    }
  } // for i, N_Cells
 }





// ====================================================================
// modify axisymmetriv drop fluid matrices due to integrals 
// on free surface with variable surface tension w.r.t heat
// ====================================================================

void FreeSurf_Axial3DHeat(TSquareMatrix2D *A11, TSquareMatrix2D *A22,
                          double *rhs1, double *rhs2,
                          BoundCondFunct2D *BoundaryCondition,
                          double dt, TFEFunction2D *Heat, double *T_IntfaceMinMax, double *Ucl, double *param)
{
  int i, j, k, l, DOF_R, DOF_L, TDOF, TDOF_R, m, TN_DOF_Local;
  int N_Cells, N_Vertices, N_Edges, Semi_implicit=0;
  int *TGlobalNumbers, *TBeginIndex;
  int comp, N_U, test_L=1, test_R=1, N_LinePoints;  
  int N_BaseFunct, *N_BaseFuncts;
  int *KCol, *RowPtr, *JointDOF, N_DOF, *TJointDOF;
  int *BeginIndex, *GlobalNumbers, *DOF, TestDOF, AnsatzDOF;
  int index1, index2;
  int count=0, count1=0, count2=0;
  int JointNumbers[MAXN_JOINTS], IJoint, N_IsoJoints;
  int T_CellNo;
  
  double t0, t1, n0, n1, normn, line_wgt;
  double *LineWeights, *zeta;
  double x0, y0, x1, y1,tx,ty,mod_t, x, y;
  double **uref, **uxiref, **uetaref;
  double uorig[MaxN_BaseFunctions2D], uxorig[MaxN_BaseFunctions2D];
  double uyorig[MaxN_BaseFunctions2D];
  double r2, r, MTF, TGrade1, TGrade2, ngrad_T;
  double *ValuesA11, *ValuesA12, *ValuesA21, *ValuesA22;
  double val, theta, factor1, factor2, angle, T_val[3], T_Marangoni, *Tvalues;
  double X_B[100], Y_B[100], T[100], r_axial, T_DOF[10], d1, d2, e1, e2, ngrad_test, ngrad_ansatz, sigma;
// double innerdomain = TDatabase::ParamDB->P13;
// double con;
// double theta2 = (3.141592654/180)*TDatabase::ParamDB->P12;
// double theta1 = (3.141592654/180)*TDatabase::ParamDB->EQ_CONTACT_ANGLE;
    
  TJoint *joint;
  TIsoBoundEdge *isoboundedge;
  TBoundComp *BoundComp;
  BoundCond Cond0, Cond1;
  FE2D FEId, TFEId;
  TFE2D *ele;
  RefTrans2D RefTrans;
  TRefTrans2D *F_K;
  TFESpace2D *fespace, *thermalspace;
  BF2DRefElements RefElement;
  QuadFormula1D LineQuadFormula;
  TQuadFormula1D *qf1;
  BaseFunct2D *BaseFuncts;
  TBaseCell *cell;
  TCollection *Coll;  
  TFEDesc2D *FeDesc, *TFeDesc;
  
 
  BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
  N_BaseFuncts = TFEDatabase2D::GetN_BaseFunctFromFE2D();

  fespace = A11->GetFESpace();
  Coll = fespace->GetCollection();
  N_Cells = Coll->GetN_Cells();

  thermalspace = Heat->GetFESpace2D();
  Tvalues=Heat->GetValues();
  TGlobalNumbers = thermalspace->GetGlobalNumbers();
  TBeginIndex = thermalspace->GetBeginIndex();

  BeginIndex = fespace->GetBeginIndex();
  GlobalNumbers = fespace->GetGlobalNumbers();

  RowPtr = A11->GetRowPtr();
  KCol = A11->GetKCol();

  ValuesA11 = A11->GetEntries();
  ValuesA22 = A22->GetEntries();

  double We = TDatabase::ParamDB->WB_NR;
  double C0 = TDatabase::ParamDB->HEAT_TANGENTIAL_STRESS_FACTOR;  // C_1/\sigma_sa

  double EQ_Angle = TDatabase::ParamDB->EQ_CONTACT_ANGLE;
  EQ_Angle = (3.141592654/180)*EQ_Angle;
  T_IntfaceMinMax[0]=1.e8;
  T_IntfaceMinMax[1] = -1.e8;

 for(i=0;i<N_Cells;i++)
  {
//      cout << endl << "CELL number: " << i << endl;
    cell = Coll->GetCell(i);
    T_CellNo = cell->GetGlobalCellNo();    
    
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

// hetrogeneous part
// 	if(x1<(innerdomain-0.01))
// 	  con = 1.0;
// 	else if(x1>(innerdomain-0.01))
// 	  con = 0.0;
// 	else
// 	  con = (0.01 + (innerdomain-x1))/0.02;
// 	
// 	EQ_Angle = (con*theta1) + ((1-con)*theta2); 


     // entries for wetting DOF
      if(y0==0.) // right wett point edge (bottom)
       {
        FeDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(FEId);
        JointDOF = FeDesc->GetJointDOF(IJoint);
        N_DOF = FeDesc->GetN_JointDOF();
        for(m=0;m<N_DOF;m++)
         {
          DOF_R =  GlobalNumbers[BeginIndex[i]+JointDOF[m]];
          fespace->GetDOFPosition(DOF_R, x, y);
          if(y==0.) // right wett point
          {
           param[0] = x;
           param[1] = y;
           param[2] = Ucl[DOF_R];    
//             cout<< "  x= "<< x <<"  y= "<< y <<endl;
//             RefTrans = TriaIsoparametric;
//             Heat->FindGradientLocal(cell, i,  x, y, T_val);
//             Heat->FindGradient( x,  y, T_val);
//             if(fabs(T_val[0])>0.5 || T_val[0]<-0.5 )
             TFEId = thermalspace->GetFE2D(T_CellNo, cell);
             TFeDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(TFEId);
             TJointDOF = TFeDesc->GetJointDOF(IJoint);
//              TN_DOF_Local = TFeDesc->GetN_JointDOF();
             TDOF_R =  TGlobalNumbers[TBeginIndex[T_CellNo]+TJointDOF[0]];
             T_val[0]=Tvalues[TDOF_R];
//              if(fabs(T_val[0])>0.5 || T_val[0]<-0.5 )
//               OutPut("x : "<<x<< " y: " << y <<"  Temperature exceeds the reference value, T= " <<T_val[0]<<endl);
	     
// temperature dep. angle
            sigma = 1. - C0* (T_val[0] -1.) ;  

            if(fabs(cos(EQ_Angle)/sigma)>1.)
             {
             if( (cos(EQ_Angle)/sigma) < -1.)
              { EQ_Angle = acos(-1.); }
            else
             { EQ_Angle = acos(1.);  }

           OutPut("  x= "<< x <<"  y= "<< y << " sigma " << sigma<<  " EQ_Angle: " <<  (180./Pi)*EQ_Angle<< endl); 
          }
         else
          {
           EQ_Angle = acos(cos(EQ_Angle)/sigma); 
//            OutPut("  x= "<< x <<"  y= "<< y << " sigma " << sigma<<  " EQ_Angle: " <<  (180./Pi)*EQ_Angle<< endl); 	   
          }
          
            r_axial = x;       // r value in the axial symmetric integral
            rhs1[DOF_R] +=  sigma*r_axial*((cos(EQ_Angle))/We); 
    
//              Surf_Force1 +=    (1. - C0* (T_val[0] -1.) )*r_axial*((cos(EQ_Angle))/We); 
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


     TFEId = thermalspace->GetFE2D(T_CellNo, cell);
     TFeDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(TFEId);
     TJointDOF = TFeDesc->GetJointDOF(IJoint);
     TN_DOF_Local = TFeDesc->GetN_JointDOF();
      for(k=0;k<TN_DOF_Local;k++)
        {
          TDOF =  TGlobalNumbers[TBeginIndex[T_CellNo]+TJointDOF[k]];
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
         OutPut( "thermal space should be linear or quadratic !! Check Freesurfint " <<endl);
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

         Heat->FindGradientLocal(cell, T_CellNo, X_B[k], Y_B[k], T_val);  
 
         r = normn/We;
         MTF = 1. - C0* (T[k] -1.);
    
         if(T_val[0] < T_IntfaceMinMax[0]) T_IntfaceMinMax[0] = T_val[0];
         if(T_val[0] > T_IntfaceMinMax[1]) T_IntfaceMinMax[1] = T_val[0];  
 
         //  (c_1\sigma_sa) \tau\cdot \grad T,  norm for integral weights
         T_Marangoni = r*C0*(t0*T_val[1] + t1*T_val[2]);

         ngrad_T = n0*T_val[1] + n1*T_val[2];
         TGrade1 =  C0*(T_val[1] - ngrad_T*n0); // vector direct product
         TGrade2 =  C0*(T_val[2] - ngrad_T*n1); // vector direct product

//          Surf_Force1  += LineWeights[k]*normn*r_axial;
  
         for(l=0;l<N_BaseFunct;l++)
          {
           TestDOF = DOF[l];

           // updating rhs
            ngrad_test= n0*uxorig[l] + n1*uyorig[l];
            d1 = MTF*(uxorig[l] - ngrad_test*n0) - uorig[l]*TGrade1;
            d2 = MTF*(uyorig[l] - ngrad_test*n1) - uorig[l]*TGrade2;

            //  rhs1
            val = r_axial*( (1-n0*n0)*d1 - n0*n1*d2 );
            val += MTF*uorig[l]; // due to axistmmetric
            val *= LineWeights[k]*r;
            rhs1[TestDOF] -= val;
//             Surf_Force1 -= val;
// 	    cout << C0 << " " << val << endl;
            //  rhs1 Marangoni convection
            val = r_axial*t0*uorig[l];
            val *= LineWeights[k]*T_Marangoni;
            rhs1[TestDOF] -= val;
//             Surf_Force1 -= val;
	    
            // rhs2
            val =  r_axial*( -n1*n0*d1 + (1-n1*n1)*d2 );
            val *= LineWeights[k]*r;
            rhs2[TestDOF] -= val;
//             Surf_Force2 -= val;
            //  rhs2 Marangoni convection
            val = r_axial*t1*uorig[l];
            val *= LineWeights[k]*T_Marangoni;
            rhs2[TestDOF] -= val;
//             Surf_Force2 -= val;
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

              val = d1*e1 + d2*e2 + MTF*(uorig[l]*uorig[m]/(r_axial*r_axial));
              val *= dt*LineWeights[k]*r*r_axial;
//               cout << "A11: " << TestDOF << " ";
//               cout << AnsatzDOF << " " << val << endl;
              ValuesA11[index1] += val;
//                   Lhs_Force1 += val;

              val = d1*e1 + d2*e2;
              val *= dt*LineWeights[k]*r*r_axial;

//               cout << "A22: " << TestDOF << " ";
//               cout << AnsatzDOF << " " << val << endl;
              ValuesA22[index1] += val;
// 	      Lhs_Force2 += val;

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
//   exit(0);
  
//   cout <<  "Surf_Force1  " <<   Surf_Force1     << " Surf_Force2  " << 	 Surf_Force2     << endl;
//     cout <<  "Lhs_Force1  " <<   Lhs_Force1     << " Lha_Force2  " << 	 Lhs_Force2    << endl;
//     exit(0);
 }

void UpdateCellwith2BDs(TCollection *coll) 
{
 int i, j, k, l, m, m0, m1, m2, N_Cells, index, comp;
 int RefLevel=1;

 double  TX[2], TY[2], x0, y0, x1, y1;
  
 TBaseCell *Me, *NeibCell;
 TJoint *Joint, *JointNeib, *innerJoint;
 
    N_Cells = coll->GetN_Cells();
   
    for(j=0;j<N_Cells;j++)
     {
        Me = coll->GetCell(j);
        k = Me->GetN_Edges();
        for(l=0;l<k;l++)
        {
         Joint = Me->GetJoint(l);
 
//       free surface edge
         if(Joint->GetType() == IsoBoundEdge  )
          {
	    index = (l+1)%3;
	    
//         is any other edge of this cell on BD?    
           for(m=0;m<k;m++)
	   {
	     if(l!=m)
	     {
	      JointNeib  = Me->GetJoint(m);
 
              if(JointNeib->GetType() == BoundaryEdge || JointNeib->GetType() == InterfaceJoint )
              {
//                Solid_Joint = (TBoundEdge *)JointNeib;
//                BoundComp = Solid_Joint->GetBoundComp();
//                comp=BoundComp->GetID();
       
               if(index==m)
                index = (m+1)%3;

//                cout<<index << "  l " << l  << " m " <<m<<endl;

               Me->Set1Refine(index);
               Me->Refine(RefLevel);        
// 	       Me->SetRegRefine();
// 	find the neib and their edge local index       
               innerJoint = Me->GetJoint(index); 
               NeibCell = innerJoint->GetNeighbour(Me);
               
               for(m1=0;m1<k;m1++)        
	        if(NeibCell->GetJoint(m1) == innerJoint)		
	         {
		  NeibCell->Set1Refine(m1);
                  NeibCell->Refine(RefLevel);  
// 		  NeibCell->SetRegRefine();
                  break;
	         }
	      }
	     } //  if(l!=m)
	   }// for(m=0;m<k;l++)
	   break;
          } // if(Joint->GetType()   
         }// endfor l
        }// endfor j
   
 
}

 
void  GetMovingBoundData(TCollection *coll, int *N_MovVert, TBoundEdge *** &Bound_Joint, 
                         TVertex *** &MovBoundVert, TIsoBoundEdge ** &Free_Joint, TBaseCell ** &Free_Cells,
                         int ** &IsoCellEdgeNos, double x, double y)
{
 int i, j, k, l, m, m0, m1, m2, N_Cells, comp;
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
          if(Joint->GetType() == BoundaryEdge || Joint->GetType() == InterfaceJoint)
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
      
//      exit(0);
     
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
          if(Joint->GetType() == BoundaryEdge || Joint->GetType() == InterfaceJoint)
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
	    
// 	    if(TY[0]<0.01)
// 	    cout<< " TX[0] " << TX[0]<<" TY[0] " << TY[0]<< " TX[1] " << TX[1]<<" TY[1] " << TY[1]<< endl;

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
       
//     for (k=0;k<N_MovVert[0];k++)
//       {
//        MovBoundVert[0][k]->GetCoords(x0, y0);
//        cout<< " x0 " << x0<<" y0 " << y0<<endl;   
//        }
// exit(0) ;

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



void GetSolidBoundData(TCollection *coll, int *N_MovVert, TBoundEdge *** &Bound_Joint, 
                       TVertex *** &MovBoundVert, TBaseCell ** &BD_Cells,
                       int &N_SolidNeibCells, TBaseCell ** &SolidNeibCells,
                         TFESpace2D *FESpace, int FluxBdID, int &N_FluxDof, int * &FluxDof)
{
 int i, j, k, l, m, m0, N, N_Cells, comp;
 int *GlobalNumbers, *BeginIndex, N_FluxDof_All, *FluxDof_all, N_FluxEdgeDof;
 int *JointDOF, *DOF, last; 
 
 double  x0, y0, x1, y1;
 
 TBaseCell *Me, *tmpCell;
 TJoint *Joint;
 TBoundEdge *Solid_Joint;
 TBoundComp *BoundComp;  
 TVertex *temp_Mov;
 TBoundEdge *tempSlip_Joint;  
 TFEDesc2D *FeDesc;
 FE2D FEId;
 TCollection *Coll_All;
 
  Coll_All = FESpace->GetCollection();
  GlobalNumbers = FESpace->GetGlobalNumbers();
  BeginIndex = FESpace->GetBeginIndex();
   
  N_Cells = coll->GetN_Cells();
  N_MovVert[3] = 0; 

  N_FluxDof_All= 0;   
   
    for(j=0;j<N_Cells;j++)
      {
        Me = coll->GetCell(j);
        k = Me->GetN_Edges();
        for(l=0;l<k;l++)
        {
         Joint = Me->GetJoint(l);
          if(Joint->GetType() == BoundaryEdge || Joint->GetType() == InterfaceJoint)
            {
             Solid_Joint = (TBoundEdge *)Joint;
             BoundComp = Solid_Joint->GetBoundComp();
             comp=BoundComp->GetID();
             
             if(comp==6)
              { N_MovVert[3]++; }  
              
              
//        assume entire BD is Neuman type       
//              if(comp==FluxBdID)
              {
               N = Me->GetGlobalCellNo();
               FEId = FESpace->GetFE2D(N, Me);
               FeDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(FEId);
               N_FluxDof_All +=FeDesc->GetN_JointDOF();
              }   
            }
            
          }// endfor l
        }// endfor j  
        
        
    //cout<<"BDComp " << 6 << " N_Vert " << N_MovVert[3] << endl;    
      cout<<"N_FluxDof_All " <<      N_FluxDof_All << endl; 


    FluxDof_all = new int[N_FluxDof_All];
     
    Bound_Joint[2] = new TBoundEdge* [N_MovVert[3]];
    MovBoundVert[3] = new TVertex* [N_MovVert[3]];     
    BD_Cells = new TBaseCell *[N_MovVert[3]];
    
    m0=0;
    N_FluxDof_All = 0;   
    
     for(j=0;j<N_Cells;j++)
      {
        Me = coll->GetCell(j);
        k = Me->GetN_Edges();
        for(l=0;l<k;l++)
        {
          Joint = Me->GetJoint(l);
          if(Joint->GetType() == BoundaryEdge || Joint->GetType() == InterfaceJoint)
          {
           Solid_Joint = (TBoundEdge *)Joint;
           BoundComp = Solid_Joint->GetBoundComp();
           comp=BoundComp->GetID();
           if(comp==6)
            {
             Bound_Joint[2][m0] = (TBoundEdge *)Joint;
             MovBoundVert[3][m0] = Me->GetVertex(l);
             BD_Cells[m0] = Me;
             m0++;
            }

//              if(comp==FluxBdID)
              {
               N = Me->GetGlobalCellNo();
               DOF = GlobalNumbers + BeginIndex[N];
               FEId = FESpace->GetFE2D(N, Me);
               FeDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(FEId);
               JointDOF = FeDesc->GetJointDOF(l);
               N_FluxEdgeDof = FeDesc->GetN_JointDOF();

               for(m=0;m<N_FluxEdgeDof;m++)
                FluxDof_all[N_FluxDof_All++] = DOF[JointDOF[m]];  
              }



           } //if(Joint->GetType() == BoundaryEdge)
 
          }// endfor l
         }// endfor j    


//  collect flux dofs

//     for (k=0;k<N_FluxDof_All;k++)
//      {
//       cout<< " k " << k <<" FluxDof " << FluxDof_all[k];      
//       FESpace->GetDOFPosition(FluxDof_all[k],x0, y0);
//        cout<< " x " << x0<<" y " << y0<<endl;
//      }
// exit(0);
    for(k=0;k<N_FluxDof_All-1;k++)
     for(l=k+1;l<N_FluxDof_All;l++)
     {
     if(FluxDof_all[k] > FluxDof_all[l])
      {  
       j = FluxDof_all[k];
       FluxDof_all[k] = FluxDof_all[l]; 
       FluxDof_all[l] = j;
      }
     }

    last = FluxDof_all[0];
    N_FluxDof = 1;
    
    for (k=1;k<N_FluxDof_All;k++)
     if(FluxDof_all[k] != last  )
      N_FluxDof++;
    
    FluxDof = new int[N_FluxDof]; 

    FluxDof[0] = FluxDof_all[0];    
    N_FluxDof = 1;
    
    for (k=1;k<N_FluxDof_All;k++)
     if(FluxDof_all[k] != FluxDof[N_FluxDof-1]  )
       FluxDof[N_FluxDof++] = FluxDof_all[k];

//     for (k=0;k<N_FluxDof;k++)
//      {
//   
//       FESpace->GetDOFPosition(FluxDof[k],x0, y0);
//       cout<< " k " << k <<" FluxDof " << FluxDof[k]<< " x " << x0<<" y " << y0<<endl;
//      }
// exit(0);

   // sort         
    for(k=0;k<N_MovVert[3]-1;k++)
      {
      for(l=k+1;l<N_MovVert[3];l++)
       {MovBoundVert[3][k]->GetCoords(x0, y0);
	MovBoundVert[3][l]->GetCoords(x1, y1);
	if(x0 < x1)
	{
	  temp_Mov = MovBoundVert[3][k];
          MovBoundVert[3][k] = MovBoundVert[3][l];
          MovBoundVert[3][l] = temp_Mov;

	  tempSlip_Joint = Bound_Joint[2][k];
	  Bound_Joint[2][k] = Bound_Joint[2][l];
	  Bound_Joint[2][l] = tempSlip_Joint;

	  tmpCell = BD_Cells[k];
	  BD_Cells[k] = BD_Cells[l];
	  BD_Cells[l] = tmpCell;  
	 }
        }
       }   

    // find neib cells containing the last vertes (next to wetting point)
     temp_Mov = MovBoundVert[3][N_MovVert[3]-1];
     N_SolidNeibCells=0;
     for(j=0;j<N_Cells;j++)
      {
        Me = coll->GetCell(j);
        k = Me->GetN_Edges();
        for(l=0;l<k;l++)
         if(Me->GetVertex(l) == temp_Mov)
          N_SolidNeibCells++;
      }

     SolidNeibCells = new TBaseCell *[N_SolidNeibCells]; 
     N_SolidNeibCells=0;
     for(j=0;j<N_Cells;j++)
      {
        Me = coll->GetCell(j);
        k = Me->GetN_Edges();
        for(l=0;l<k;l++)
         if(Me->GetVertex(l) == temp_Mov)
           SolidNeibCells[N_SolidNeibCells++] = Me;
      }


   delete [] FluxDof_all;

//   for (k=0;k<N_MovVert[3];k++)
//       {
//        MovBoundVert[3][k]->GetCoords(x0, y0);
//        cout<< " x " << x0<<" y " << y0<<endl;
//        }
//     cout<< " m0 " << m0 <<endl;
//     cout<< " Cell No " << BD_Cells[N_MovVert[3]-2] ->GetClipBoard() <<endl;
//     cout<< " Cell No " << BD_Cells[N_MovVert[3]-1] ->GetClipBoard() <<endl;
//      
//  exit(0); 
} // GetSolidBoundData


void RemeshAxial3D_ImpDrop(TDomain * &Domain, TFESpace2D ** &FESpaces_All,
                           TFEVectFunct2D ** &FEVectFuncts_All, TFEFunction2D ** &FEFunctions_All, int *N_MovVert, 
                           TBoundEdge *** &Bound_Joint, TVertex *** &MovBoundVert, TIsoBoundEdge ** &Free_Joint,
                           TBaseCell ** &Free_Cells, int ** &IsoCellEdgeNos, double ** &Sol_All,  double ** &Rhs_All, 
                           TSquareStructure2D ** &SquareStructure_All, TStructure2D ** &Structure_All,
                           TSquareMatrix2D ** &SqMat_All, TMatrix2D ** &Mat_All,
                           TBaseCell** &NSE_Cells, TBaseCell** &Solid_Cells, 
                           TCollection * &NSE_coll, TCollection * &Solid_coll,
                           int *N_GidDofs, int *N_GridActive, int *N_GridBdDofs,
                           TFESpace2D * &Grid_space_NSE,  TFESpace2D * &Grid_space_S,
                           TBaseCell ** &BD_Cells, int &N_SolidNeibCells, TBaseCell ** &SolidNeibCells,
                           double *&GridSol_NSE, double *&GridRhs_NSE, double *&GridSol_S, double *&GridRhs_S,
                           TFEVectFunct2D * &GridVect_NSE, TFEVectFunct2D * &GridVect_S,
                           TFEFunction2D * &GridG1_NSE, TFEFunction2D * &GridG2_NSE, 
                           TFEFunction2D * &GridG1_S, TFEFunction2D * &GridG2_S,
                           TSquareStructure2D * &SquareStructure_NSE, TSquareStructure2D * &SquareStructure_S,
                           TSquareMatrix2D ** &GridSqMat_NSE, TSquareMatrix2D ** &GridSqMat_S,
                           int &N_FluxDof, int * &FluxDof, int * &Triangles, int * &PointNeighb, int &maxEpV)
{
  int i, j, k, l, N_G, N_Cells, ORDER, VSP, N_DOF, N_ThermalDOF, N_heatfuncDOF;
  int In_Index, CurrComp, Old_N_Cells, Old_N_RootCells, CurrVertex, N_Joints, N_Vertices, ID;
  int N_RootCells, a, b, Neighb_tmp, Neib[2];
  int CurrNeib, len1, len2, pressure_space_code, N_U, N_P, N_Unknowns, comp;
  int *JointDOF, *DOF, *GlobalNumbers, *BeginIndex, Indextemp;
  int N_NSE_Cells, N_Solid_Cells, N_U_output, N, N_refX;

  double d, h, t, dt, x, y, tx, ty, x1, x2, y1, y2, Lx, Ly, Rx, Ry, *S_BX, *A_BY, *S_CX;
  double refX, area, *Coordinates, left, right, bottom, top, T_a, T_b, *sol, *ValuesU2, *Tsol;

  TBoundPart *BoundPart;
  TBoundComp *BoundComp;
  TBdLine *UpdateSlipBound, *UpdateAxialBound, *UpdateSlipBoundSolid;
  TBdLine *SolidBound1, *SolidBound2, *SolidBound3;  
  TCollection *coll, *Old_Coll, *mortarcoll = NULL;
  TBaseCell *cell, **CellTree, **Old_CellTree,  **NSE_Cells_new, **Solid_Cells_new;
  TVertex **VertexDel, **NewVertices;
  TJoint *Joint;
  TFESpace2D *velocity_space, *pressure_space, *thermal_space, *pressure_space_output;
  TFEVectFunct2D *u; 
  TFEFunction2D *p, *Heat; 
  TBoundEdge *Solid_Joint;
  FE2D FeId, *fes;
  TFEDesc2D *FeDesc;
 
  // strings
  char ReadinDat[] = "readin.dat";
  char TString[] = "T";
  char HeatString[] = "H";  
  char NameString[]  = "name";
  char UString[] = "U";
  char PString[] = "P";
  char WString[] = "W";
  char VortString[] = "Vort";
  char DivString[] = "Div";  
  
  std::ostringstream opts;
  std::ostringstream os;
  os << " ";
  opts << " ";
// 
  struct triangulateio In, Out;

  
//   double refX, Xi[5] = {0., 0.,  8., 8., 0.2};
//   double       Yi[5] = {0.,-4., -4., 0., 0.};
//   double refX, Xi[5] = {0., 0.,  29.45, 29.45, 0.2};
//   double       Yi[5] = {0.,-0.775, -0.775, 0., 0.};
  
#ifdef __ConvergenceStudy__  
  double Xi[5] = {0., 0.,  4., 4., 1.};
  double Yi[5] = {0.,-2., -2., 0., 0.};
//   double Xi[5] = {0., 0.,  2., 2., 1.};
//   double Yi[5] = {0.,-1., -1., 0., 0.};    
#else
  double Xi[5] = {0., 0.,  8., 8., 0.2};
  double Yi[5] = {0.,-4., -4., 0., 0.};

//   double Xi[5] = {0., 0.,  29.45, 29.45, 0.2};
//   double Yi[5] = {0.,-0.775, -0.775, 0., 0.};
#endif

  
  // free surface vertices
  if(N_MovVert[2]-1<10) N = N_MovVert[2];
  else N=10;
    
  d = 0;
  for(i=0;i<N;i++) // without last point
   {
    MovBoundVert[2][i]->GetCoords(x1, y1);
    MovBoundVert[2][i+1]->GetCoords(x2, y2);
     d +=sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2)  );
   }
   h = d/((double)N);
  
  MovBoundVert[0][0]->GetCoords(Lx, Ly);
  MovBoundVert[2][0]->GetCoords(Rx, Ry);
  
  d = Rx-Lx;
  k = (int)(d/h); // No of intervals with step length 0.01
  if(k<2) k=2;     // minimum two intervals
  t = d/(double)k;
  N_MovVert[0] = k;
  S_BX = new double[N_MovVert[0]];
     
  for(i=0;i<N_MovVert[0];i++)
   {
    S_BX[i]= Lx + (double)i*t;
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
   //cout<<i<< " x :" << 0 << " -----------------y: " <<A_BY[i]<< endl;
  }

// cout << " N_MovVert[1] " << N_MovVert[1] << endl;
// exit(0);

  //solid BD points
  Xi[4] = Rx;
  Yi[4] = Ry;

  refX = Xi[4] + 0.2;    
     
   N_refX = (int)((refX - Xi[4])/h);
   N_MovVert[3] = N_refX + 10;

   S_CX = new double[N_MovVert[3]];   
   
   N=0;
   t= (refX - Xi[3])  / (double)(N_MovVert[3]-N_refX);
   for(i=0;i<N_MovVert[3]-N_refX;i++)  
    S_CX[N++]= Xi[3] + (double)i*t;
   
   t= (Xi[4] - refX)  / N_refX; 
   for(i=0;i<N_refX;i++)  
    S_CX[N++]= refX + (double)i*t;   
   
//    cout << Rx << " N : "<<   N << " " << N_MovVert[3] << endl;
// exit(0);
  area = TDatabase::ParamDB->Area; 
//======================================================================
// Triangular for grid generation begin
//======================================================================
  BoundPart = Domain->GetBdPart(0);
  UpdateSlipBound = (TBdLine*)BoundPart->GetBdComp(0);
  UpdateAxialBound = (TBdLine*)BoundPart->GetBdComp(2);
  UpdateSlipBoundSolid = (TBdLine*)BoundPart->GetBdComp(6);
  SolidBound1 = (TBdLine*)BoundPart->GetBdComp(3);
  SolidBound2 = (TBdLine*)BoundPart->GetBdComp(4);
  SolidBound3 = (TBdLine*)BoundPart->GetBdComp(5);  
  
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

  int N_Hori = 10;
  int N_Verti = 5;
  
  In.numberofpoints = N_MovVert[0]+N_MovVert[1]+N_MovVert[2]+N_MovVert[3]+N_Hori+2*N_Verti-1;
  In.pointlist = new double[2*In.numberofpoints];
  In.pointmarkerlist = new int[In.numberofpoints];
  In.numberofpointattributes = 0;

  In.numberofsegments = In.numberofpoints+1;
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
    In.pointlist[2*In_Index+1] = 0.0;
//     cout<<In_Index<< " S_BX :" << S_BX[i]<< " -----------------y: " <<0.0<< endl;
    In.pointmarkerlist[In_Index] = CurrComp;
    In.segmentlist[2*In_Index] = In_Index;
    In.segmentlist[2*In_Index+1] = In_Index+1;
    In.segmentmarkerlist[In_Index] = CurrComp;
    In_Index++;
   }
   cout<<endl; 
   CurrComp++;    
   Indextemp = In_Index; 
   
  // points and segments on the free boundary (marker=2)
  for(i=0;i<N_MovVert[2];i++) // without last point
   {
    MovBoundVert[2][i]->GetCoords(tx, ty);
    In.pointlist[2*In_Index] = tx;
    In.pointlist[2*In_Index+1] = ty;
//     cout<<In_Index<< " tx :" << tx<< " -----------------y: " <<ty<< endl;
    In.pointmarkerlist[In_Index] = CurrComp;
    In.segmentlist[2*In_Index] = In_Index;
    In.segmentlist[2*In_Index+1] = In_Index+1;
    In.segmentmarkerlist[In_Index] = CurrComp;
    In_Index++;
   }
   cout<<endl; 
   CurrComp++;
   
  // points and segments on the solid boundary (marker=1)
  for(i=0;i<N_MovVert[1];i++) // without last point
   {
    In.pointlist[2*In_Index] = 0.0;
    In.pointlist[2*In_Index+1] = A_BY[i];
//     cout<<In_Index<< " x :" << 0.0<< " -----------------y: " <<A_BY[i]<< endl;
    In.pointmarkerlist[In_Index] = CurrComp;
    In.segmentlist[2*In_Index] = In_Index;
    In.segmentlist[2*In_Index+1] = In_Index+1;
    In.segmentmarkerlist[In_Index] = CurrComp;
    In_Index++;
   }
  In.segmentlist[2*(In_Index-1)+1] = 0; 
  
  
//    cout<<endl; 
   CurrComp++;
   dt= (Yi[1] - Yi[0])  / (double)N_Verti;
   x = Xi[0];
   y = Yi[0];
   t = y;
   
  // (0,0) point is already defined
  In.segmentlist[2*In_Index] = 0;
  In.segmentlist[2*In_Index+1] = In_Index;
  In.segmentmarkerlist[In_Index] = CurrComp;
  In_Index++;  
  
  
   for(i=1;i<N_Verti;i++) // without last point
    {
     y = t + ((double)(i))*dt;
   
      In.pointlist[2*(In_Index-1)] = x;
      In.pointlist[2*(In_Index-1)+1] = y;
//       cout<<In_Index<<" x : "<< x << " y : "<< y<<endl;
      In.pointmarkerlist[(In_Index-1)] = CurrComp;
      
      In.segmentlist[2*In_Index] = (In_Index-1);
      In.segmentlist[2*In_Index+1] = In_Index;
      In.segmentmarkerlist[In_Index] = CurrComp;
      In_Index++;
   
    }
    
//    cout<<endl; 
   CurrComp++;


   dt= (Xi[2] - Xi[1])  / (double)N_Hori;
   x = Xi[1];
   y = Yi[1];
   t = x;

   for(i=0;i<N_Hori;i++) // without last point
    {
      In.pointlist[2*(In_Index-1)] = x;
      In.pointlist[2*(In_Index-1)+1] = y;
//       cout<<In_Index<<" x : "<< x << " y : "<< y<<endl;
      In.pointmarkerlist[In_Index-1] = CurrComp;
      
      In.segmentlist[2*In_Index] = In_Index-1;
      In.segmentlist[2*In_Index+1] = In_Index;
      In.segmentmarkerlist[In_Index] = CurrComp;
      In_Index++;
      x = t + ((double)(i+1))*dt;
    }

//    cout<<endl; 
   CurrComp++;
   dt= (Yi[3] - Yi[2])  / (double)N_Verti;
   x = Xi[2];
   y = Yi[2];
   t = y;

   for(i=0;i<N_Verti;i++) // without last point
    {
      In.pointlist[2*(In_Index-1)] = x;
      In.pointlist[2*(In_Index-1)+1] = y;
//       cout<<In_Index<<" x : "<< x << " y : "<< y<<endl;
      In.pointmarkerlist[(In_Index-1)] = CurrComp;
      
      In.segmentlist[2*In_Index] = In_Index-1;
      In.segmentlist[2*In_Index+1] = In_Index;
      In.segmentmarkerlist[In_Index] = CurrComp;
      In_Index++;
      y = t + ((double)(i+1))*dt;
    }

   CurrComp++;

//    cout<<endl;  
  for(i=0;i<N_MovVert[3];i++) // without last point
   {
    In.pointlist[2*(In_Index-1)] = S_CX[i];
    In.pointlist[2*(In_Index-1)+1] = 0.0;
//     cout<<In_Index << " x :" << S_CX[i]<< " -----------------y: " <<0.0<< endl;
    In.pointmarkerlist[(In_Index-1)] = CurrComp;
    
    In.segmentlist[2*In_Index] = (In_Index-1);
    In.segmentlist[2*In_Index+1] = In_Index;
    In.segmentmarkerlist[In_Index] = CurrComp;
    In_Index++;
   }
   
   In.segmentlist[2*(In.numberofsegments-1)+1] = Indextemp; 

   delete [] S_CX;
   
// cout << "  Indextemp " <<  Indextemp << endl;
// exit(0);

   if(PointNeighb)
     delete [] PointNeighb;
  
   free(Triangles);
   maxEpV = 0;
  
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
  UpdateSlipBoundSolid->SetParams(Xi[3], Yi[3], Xi[4]-Xi[3],Yi[4]-Yi[3]);
  SolidBound1->SetParams(Xi[0], Yi[0], Xi[1]-Xi[0],Yi[1]-Yi[0]);
  SolidBound2->SetParams(Xi[1], Yi[1], Xi[2]-Xi[1],Yi[2]-Yi[1]); 
  SolidBound3->SetParams(Xi[2], Yi[2], Xi[3]-Xi[2],Yi[3]-Yi[2]);
  
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
//   delete [] PointNeighb;
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
//   if(Out.trianglelist!=NULL) {
//     free(Out.trianglelist); Out.trianglelist = NULL;}
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
 
      // write grid into an Postscript file
      os.seekp(std::ios::beg);
      os << "Domain" << ".ps" << ends;
      Domain->PS(os.str().c_str(),It_Finest,0);
//      exit(0);
    
//======================================================================
// construct all finite element spaces
//======================================================================
  ORDER = 0;
  VSP = TDatabase::ParamDB->VELOCITY_SPACE;

  if (abs(VSP) > 20)
   {ORDER = abs(VSP) - 20;}
  else if ( abs(VSP) > 10)
    {ORDER = abs(VSP) - 10;}
  else ORDER = abs(VSP);

  
  coll=Domain->GetCollection(It_Finest, 0);
  N_Cells = coll->GetN_Cells();
 
  N_NSE_Cells = 0;

  for(i=0; i<N_Cells; i++)
   {
    cell = coll->GetCell(i);
    cell->SetGlobalCellNo(i);          
    N_Vertices = cell->GetN_Vertices();     

    t=0;
    for(j=0; j<N_Vertices; j++)
    {
     cell->GetVertex(j)->GetCoords(x, y);
     t +=y;
    }
 
    t /=(double)N_Vertices;
 
    if(t>0) // NSE cells
     {
      cell->SetPhase_ID(0); 
      N_NSE_Cells++;         
     }
    else
     { cell->SetPhase_ID(1); } 
   }
   
   
//   delete [] NSE_Cells;
//   delete [] Solid_Cells;
   
  NSE_Cells_new = new TBaseCell*[N_NSE_Cells];
  Solid_Cells_new = new TBaseCell*[N_Cells - N_NSE_Cells];  
  
  N_NSE_Cells = 0;
  N_Solid_Cells = 0;
  for(i=0; i<N_Cells; i++)
   {
    cell = coll->GetCell(i);
    ID = cell->GetPhase_ID(); 
 
    if(ID==0) // NSE cells
     {
      NSE_Cells_new[N_NSE_Cells] = cell;        
      N_NSE_Cells++;    
     }
    else
     {
      Solid_Cells_new[N_Solid_Cells] = cell;        
      N_Solid_Cells++;      
     }
    
   } 
 
  OutPut("Number of cells, NSE cells, Solid Cells: " << N_Cells << ",  " << N_NSE_Cells <<  
                                             ",  " << N_Solid_Cells <<endl);
  NSE_coll = new TCollection(N_NSE_Cells, NSE_Cells_new); 
  Solid_coll = new TCollection(N_Solid_Cells, Solid_Cells_new); 

//======================================================================
// construct all finite element spaces
//======================================================================
  // get velocity and pressure spacess
  GetVelocityAndPressureSpace(NSE_coll,BoundCondition,
                              mortarcoll, velocity_space,
                              pressure_space, &pressure_space_code,
                              TDatabase::ParamDB->VELOCITY_SPACE,
                              TDatabase::ParamDB->PRESSURE_SPACE);
  N_U = velocity_space->GetN_DegreesOfFreedom();
  N_P = pressure_space->GetN_DegreesOfFreedom();
  GlobalNumbers = velocity_space->GetGlobalNumbers();
  BeginIndex = velocity_space->GetBeginIndex();

  delete FESpaces_All[4];
  GetVelocityAndPressureSpace(coll,BoundCondition_output,
                              mortarcoll, FESpaces_All[4],
                              pressure_space_output, &pressure_space_code,
                              TDatabase::ParamDB->VELOCITY_SPACE,
                              TDatabase::ParamDB->PRESSURE_SPACE);
   N_U_output = FESpaces_All[4]->GetN_DegreesOfFreedom();

  // mesh velocity space 
  delete FESpaces_All[2];
  delete Grid_space_NSE;
  delete Grid_space_S;
  
  FESpaces_All[2] = new TFESpace2D(coll, NameString, WString, GridBoundCondition, 1, NULL);

  Grid_space_NSE = new TFESpace2D(NSE_coll, NameString, WString, GridBoundCondition, 1, NULL);
  Grid_space_S = new TFESpace2D(Solid_coll, NameString, WString, GridBoundCondition, 1, NULL);    

  N_GidDofs[0] = FESpaces_All[2]->GetN_DegreesOfFreedom();;
  N_GidDofs[1] = Grid_space_NSE->GetN_DegreesOfFreedom(); 
  N_GidDofs[2] = Grid_space_S->GetN_DegreesOfFreedom(); 

  N_GridActive[0] = FESpaces_All[2]->GetActiveBound();
  N_GridActive[1] = Grid_space_NSE->GetActiveBound(); 
  N_GridActive[2] = Grid_space_S->GetActiveBound();    

  N_G = N_GidDofs[0];

  N_GridBdDofs[0] = N_GidDofs[0] - N_GridActive[0];   
  N_GridBdDofs[1] = N_GidDofs[1] - N_GridActive[1];
  N_GridBdDofs[2] = N_GidDofs[2] - N_GridActive[2];
       
  // thermal space
  if(TDatabase::ParamDB->ANSATZ_ORDER<100)
  {  
   thermal_space = new TFESpace2D(coll, NameString, TString, HeatBoundCondition,
                                 TDatabase::ParamDB->ANSATZ_ORDER, NULL);
   }
   else                              
   {     
    fes = new FE2D[N_Cells];
    GetHetroFEs(N_Cells, coll, Triangles, PointNeighb, maxEpV, fes); 
    thermal_space = new TFESpace2D(coll, NameString, TString, HeatBoundCondition,
                                   fes, NULL);    
    delete [] fes;
   }
   
#ifdef __HEATLINE__     
   // heatfunction space
   delete FESpaces_All[5];
   if(TDatabase::ParamDB->ANSATZ_ORDER<100)
   {
   FESpaces_All[5] = new TFESpace2D(coll, NameString, HeatString, HeatFuncBoundCondition,
                                  TDatabase::ParamDB->ANSATZ_ORDER, NULL);
   }
   else                              
   {
    fes = new FE2D[N_Cells];
    GetHetroFEs(N_Cells, coll, Triangles, PointNeighb, maxEpV, fes);     
    FESpaces_All[5] = new TFESpace2D(coll, NameString, HeatString, HeatFuncBoundCondition,
                                  fes, NULL);    
   }
   
//    FESpaces_All[5] =  heatfunc_space;  
   N_heatfuncDOF = FESpaces_All[5]->GetN_DegreesOfFreedom();
    
//    cout<<"Thermal DOFs : "<<N_thermalDOF<<"\nHeat function DOFs : "<<N_heatfuncDOF<<"\n";
#endif

   
  N_ThermalDOF = thermal_space->GetN_DegreesOfFreedom();
//===========================================================================
  delete [] Bound_Joint[0];
  delete [] Bound_Joint[1];                      
  delete [] Bound_Joint[2]; 
  
  delete [] MovBoundVert[0];
  delete [] MovBoundVert[1];   
  delete [] MovBoundVert[2];  
  delete [] MovBoundVert[3]; 
  
  delete [] Free_Joint;    
  delete [] Free_Cells;  
  delete [] IsoCellEdgeNos[0];  
  delete [] IsoCellEdgeNos[1];  
  
  delete [] FluxDof;
  
  Ry = 0.;
  GetMovingBoundData(NSE_coll, N_MovVert, Bound_Joint, MovBoundVert, Free_Joint,
                     Free_Cells, IsoCellEdgeNos, Rx, Ry);

  delete []  BD_Cells;  
  delete []  SolidNeibCells;
    
  GetSolidBoundData(Solid_coll, N_MovVert, Bound_Joint, MovBoundVert, BD_Cells,
                    N_SolidNeibCells, SolidNeibCells, thermal_space, 0,  N_FluxDof, FluxDof);

//======================================================================
// construct all finite element functions
//======================================================================
  N_Unknowns = 2*N_U + N_P;
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

  for(i=0;i<N_NSE_Cells;i++)
   {
    cell = NSE_coll->GetCell(i);
    k = cell->GetN_Edges();
    
    for(l=0;l<k;l++)
     {
      Joint = cell->GetJoint(l); 
      if(Joint->GetType() == BoundaryEdge || Joint->GetType() == InterfaceJoint)
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
// 	     cell->GetVertex(l)->GetCoords(x, y);
// 	     cout<< i << " " << l << " " << N_DOF << " x " << x << " y " << y << endl;
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
  
  // velo in all domains
  delete [] Sol_All[3]; 
  delete FEVectFuncts_All[2];  
  
  Sol_All[3] = new double[2*N_U_output]; 
  memset(Sol_All[3], 0, 2*N_U_output*SizeOfDouble);
  FEVectFuncts_All[2] =  new TFEVectFunct2D(FESpaces_All[4], UString, UString, Sol_All[3], N_U_output, 2);
 
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

  delete [] GridSol_NSE;
  delete [] GridRhs_NSE;  
  delete GridVect_NSE;

  GridSol_NSE = new double[2*N_GidDofs[1]];
  GridRhs_NSE = new double[2*N_GidDofs[1]];     
  memset(GridSol_NSE, 0, 2*N_GidDofs[1]*SizeOfDouble);   
  GridVect_NSE  = new TFEVectFunct2D(Grid_space_NSE, WString, WString, GridSol_NSE, N_GidDofs[1], 2);
  GridG1_NSE = GridVect_NSE->GetComponent(0);
  GridG2_NSE = GridVect_NSE->GetComponent(1);

  
  delete [] GridSol_S;
  delete [] GridRhs_S;  
  delete GridVect_S;
  
  GridSol_S = new double[2*N_GidDofs[2]];
  GridRhs_S = new double[2*N_GidDofs[2]];     
  memset(GridSol_S, 0, 2*N_GidDofs[2]*SizeOfDouble);   
  GridVect_S  = new TFEVectFunct2D(Grid_space_S, WString, WString, GridSol_S, N_GidDofs[2], 2);
  GridG1_S = GridVect_S->GetComponent(0);
  GridG2_S = GridVect_S->GetComponent(1);

//======================================================================
// thermal space finite element functions
//======================================================================    
  delete [] Rhs_All[2];
  Rhs_All[2] = new double[N_ThermalDOF]; 
  memset(Rhs_All[2], 0, N_ThermalDOF*SizeOfDouble);
    
  Tsol = new double[N_ThermalDOF];
  
  // thermal fefunction
  Heat = new TFEFunction2D(thermal_space, TString, TString, Tsol, N_ThermalDOF);

#ifdef __ENERGY__
  // interpolate from old fefunction
  Heat->Interpolate(FEFunctions_All[5]);
#else  
    memset(Tsol, 0, N_ThermalDOF*SizeOfDouble);
    memset(Rhs_All[2], 0, N_ThermalDOF*SizeOfDouble);
#endif 

    
  delete FESpaces_All[3];  
  delete [] Sol_All[2];
  delete FEFunctions_All[5];  
  
  FESpaces_All[3] = thermal_space;
  Sol_All[2] = Tsol;
  FEFunctions_All[5] = Heat;
  
 if(TDatabase::ParamDB->REACTOR_P22>0)
  {  
   delete [] Sol_All[4];
   delete [] Sol_All[5];  

   Sol_All[4] = new double[N_GidDofs[1]];  
   Sol_All[5] = new double[N_GidDofs[1]];  
  
   memset(Sol_All[4], 0, N_GidDofs[1]*SizeOfDouble);
   memset(Sol_All[5], 0, N_GidDofs[1]*SizeOfDouble);   
  
   
   delete FEFunctions_All[6]; 
   delete FEFunctions_All[7];   
   
   // vorticity fefunction
   FEFunctions_All[6] = new TFEFunction2D(Grid_space_NSE, VortString, VortString, Sol_All[4], N_GidDofs[1]);
   // divergence fefunction
   FEFunctions_All[7] = new TFEFunction2D(Grid_space_NSE, DivString, DivString, Sol_All[5], N_GidDofs[1]);  
  }  
  
#ifdef __HEATLINE__  
  delete Sol_All[6];
  delete Rhs_All[3];
  Sol_All[6] = new double[N_heatfuncDOF];
  Rhs_All[3] = new double[N_heatfuncDOF];

  memset(Sol_All[6], 0, N_heatfuncDOF*SizeOfDouble);
  memset(Rhs_All[3], 0, N_heatfuncDOF*SizeOfDouble);

  // heatfunction fefunction
  delete FEFunctions_All[8];
  FEFunctions_All[8] = new TFEFunction2D(FESpaces_All[5], HeatString, HeatString, Sol_All[6], N_heatfuncDOF);

  delete SquareStructure_All[3];
  SquareStructure_All[3] = new TSquareStructure2D(FESpaces_All[5]);
  SquareStructure_All[3]->Sort();
  
  delete SqMat_All[12];
  SqMat_All[12]  = new TSquareMatrix2D(SquareStructure_All[3]); // Heatfunc_A 
#endif  
  
// ======================================================================
// allocate memory for all matrices
// ======================================================================  
  delete Structure_All[0];  delete Structure_All[1];
  
  Structure_All[0] = new TStructure2D(FESpaces_All[1], FESpaces_All[0]);  // B
  Structure_All[1] = new TStructure2D(FESpaces_All[0], FESpaces_All[1]); // BT
  
  delete SquareStructure_All[0]; delete SquareStructure_All[1]; delete SquareStructure_All[2];
  //velo 
  SquareStructure_All[0] = new TSquareStructure2D(FESpaces_All[0]);  
  SquareStructure_All[0]->Sort();  
  
  // grid 
  SquareStructure_All[1] = new TSquareStructure2D(FESpaces_All[2]); 
  SquareStructure_All[1]->Sort();  
  
  delete SquareStructure_NSE; delete SquareStructure_S;
  
  SquareStructure_NSE = new TSquareStructure2D(Grid_space_NSE); 
  SquareStructure_NSE->Sort();
    
  SquareStructure_S = new TSquareStructure2D(Grid_space_S);   
  SquareStructure_S->Sort();     
    
  //thermal
  SquareStructure_All[2] = new TSquareStructure2D(FESpaces_All[3]);
  SquareStructure_All[2]->Sort(); 
  
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

   for(i=0; i<4; i++)
   {
    delete GridSqMat_NSE[i];
    GridSqMat_NSE[i] = new TSquareMatrix2D(SquareStructure_NSE);  
    
    delete GridSqMat_S[i];
    GridSqMat_S[i] = new TSquareMatrix2D(SquareStructure_S);  
   }    


  // for heat
  for(i=10; i<12; i++)
   {
    delete SqMat_All[i];
    SqMat_All[i] = new TSquareMatrix2D(SquareStructure_All[2]);  
   }  
   
   
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

  
  delete [] NSE_Cells;
  delete [] Solid_Cells; 
  
  NSE_Cells = NSE_Cells_new;
  Solid_Cells = Solid_Cells_new;
  
  delete [] S_BX;
  delete [] A_BY; 
  
// exit(0);  
}// RemeshAxial3D_ImpDrop


void GetFluxVect(int N_Mat, TSquareMatrix2D **SQMATRICES_HEAT, double *rhs, int N_FluxDof,
                 int *FluxDof, double *Flux_A, double *Flux_M, double &Flux_F)
{
 int i, j, k, N_Dof, *Row, *KCol, dof, N_Active, end, col;
//  int *BdDof;
 
 double *A_Values, *M_Values;

  
 TSquareMatrix *M, *A;
 
  A = SQMATRICES_HEAT[0];
  M = SQMATRICES_HEAT[1]; 
  
   //assume that both mat have same strucure
   N_Dof = M->GetN_Columns();
   Row = M->GetRowPtr();
   KCol = M->GetKCol();
   N_Active = M->GetActiveBound();
   
   if( FluxDof[N_FluxDof-1] > N_Active )   
    {
      OutPut(" N_Active " << N_Active << " FluxDof " << FluxDof[N_FluxDof-1]<<endl);
     OutPut("Flux Dof is Dirichlet DOF, compute vector before modifying the matrices !!! ");
     exit(0);     
    }

  memset(Flux_A, 0, N_Dof*SizeOfDouble);
  memset(Flux_M, 0, N_Dof*SizeOfDouble);
  Flux_F = 0;
  
//   BdDof = new int[N_Dof];
//   memset(BdDof, 0, N_Dof*SizeOfInt);

//   for(i=0; i<N_FluxDof; i++)
//    BdDof[FluxDof[i]] = 1;
  
   A_Values = A->GetEntries();   
   M_Values = M->GetEntries(); 
  
   for(i=0; i<N_FluxDof; i++)
    {
     dof = FluxDof[i];     
     end   = Row[dof+1]; 

     Flux_F += rhs[dof];
     
    for(j=Row[dof]; j<end; j++)
     {
      col = KCol[j];
      
      //col has to be on the BD, esle Mat vale is zero
      // i.e., col has to be in FluxDof      
//       if(BdDof[col]==1)
       {
        Flux_M[col] += M_Values[j];
        Flux_A[col] += A_Values[j];
       }
      
     }
     
    }
          
//   delete [] BdDof;

  
//   for(i=0; i<N_FluxDof; i++)
//     cout<< i << " FluxDOF " << FluxDof[i] << " Flux_M " << Flux_M[FluxDof[i]] << " Flux_A " << Flux_A[FluxDof[i]] <<endl;
    
//     
//   cout<< " GetFluxVect " <<endl;
//   exit(0);
}

 

void GetInterfaceMinMaxT(TFEFunction2D *Heat, double *T_IntfaceMinMax)
{
  int j, k, l, m;
  int N_Cells, N_Joints, N_DOF_Local;
  int *BeginIndex, *GlobalNumbers, DOF, *JointDOF;
  
  double T, *Values;
  
  TBaseCell *cell;
  TCollection *coll;
  TFESpace2D *FESpace;
  TBaseCell *Me;
  TJoint *Joint;
  TBoundEdge *Solid_Joint;
  TBoundComp *BoundComp;    
  FE2D FEId;
  TFEDesc2D *FeDesc;
  
  FESpace = Heat->GetFESpace2D();
  BeginIndex = FESpace->GetBeginIndex();
  GlobalNumbers = FESpace->GetGlobalNumbers();
  Values = Heat->GetValues();
  
  coll = FESpace->GetCollection();
  N_Cells = coll->GetN_Cells();
  
  
  T_IntfaceMinMax[0]=1.e8;
  T_IntfaceMinMax[1] = -1.e8; 
  
  for(j=0;j<N_Cells;j++)
   {
    Me = coll->GetCell(j);
    k = Me->GetN_Edges();
    for(l=0;l<k;l++)
     {
      Joint = Me->GetJoint(l);
       if(Joint->GetType() == BoundaryEdge || Joint->GetType() == InterfaceJoint)
        {
         Solid_Joint = (TBoundEdge *)Joint;
         BoundComp = Solid_Joint->GetBoundComp();
             
         if(BoundComp->GetID()==0)
          { 
           FEId = FESpace->GetFE2D(j, Me);
           FeDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(FEId);
           JointDOF = FeDesc->GetJointDOF(l);
           N_DOF_Local = FeDesc->GetN_JointDOF();
   
           for(m=0;m<N_DOF_Local;m++)
            {
             DOF =  GlobalNumbers[BeginIndex[j]+JointDOF[m]];
             T   = Values[DOF];         
             if(T < T_IntfaceMinMax[0]) T_IntfaceMinMax[0] = T;
             if(T > T_IntfaceMinMax[1]) T_IntfaceMinMax[1] = T;  
// //          cout<<k<< " T " << Tvalues[TDOF] <<endl;
            }    
          }   
         }            
       }// endfor l
     }// endfor j 
 
}// GetInterfaceMinMaxT

void Get_Heat(TFEFunction2D *Heat, double *val_out)
 {
  int i,j,k,l, polydegree, Phase_No;
  int N_Cells, N_Joints, N_Vertices;
  int N_QFPoints, N_BF, Phase_ID;
  int *BeginIndex, *GlobalNumbers, *DOF;
  
  double *weights, *xi, *eta;
  double values[MaxN_QuadPoints_2D][MaxN_BaseFunctions2D];
  double AbsDetjk[MaxN_QuadPoints_2D], X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];
  double *Values;
  double Mult, r_axial;
  double volume_L=0., val_loc, T_nodal;
  double volume_S=0.; 
  double T;

  TJoint *joint;
  FE2D FEid;
  TBaseFunct2D *bf;
  TRefTrans2D *F_K;    
  TBaseCell *cell;
  TCollection *coll;
 
  TFESpace2D *FESpace;
  JointType jointtype;
  BoundTypes bdtype;
  RefTrans2D RefTrans;
  boolean IsIsoparametric;
  QuadFormula2D QuadFormula;
  TQuadFormula2D *qf2;

  FESpace = Heat->GetFESpace2D();
  BeginIndex = FESpace->GetBeginIndex();
  GlobalNumbers = FESpace->GetGlobalNumbers();
  Values = Heat->GetValues();


  coll = FESpace->GetCollection();
  N_Cells = coll->GetN_Cells();

  val_out[0] = 0;
  val_out[1] = 0;  
 
  val_out[2] = 0;
  val_out[3] = 0; 
  
  val_out[4] = 0;     
  
  
  for(i=0;i<N_Cells;i++)
   {
    cell = coll->GetCell(i);
    Phase_ID = cell->GetPhase_ID();     
    
    
    FEid = FESpace->GetFE2D(i, cell);

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
        cout <<"X[k] negative in Get_KE change Quad rule " 
             <<  X[k] << " " << Y[k] << endl;
        exit(0);
       }

      Mult = r_axial*weights[k]*AbsDetjk[k]; 
      
      val_loc = 0.;
      for(l=0;l<N_BF;l++)
       {
        j = DOF[l]; 
        val_loc += values[k][l]*Values[j];
       }

     if(Phase_ID==0)
     {
      volume_L += Mult;
      val_out[0]  += (val_loc * Mult); 
     }
     else
     {
      volume_S += Mult;
      val_out[1]  += (val_loc * Mult);       
     }
     
     val_out[4] += (val_loc * Mult);    
     
     
    }  // for(k=0;k<N_QF
   } // endfor i
    
  val_out[2] = volume_L;
  val_out[3] = volume_S;   
  
//   cout<< " Volume_L: "<< volume_L<<" Volume_S: "<< volume_S<< endl;
//   cout<< " val_out_L: "<< val_out[0] <<  " val_out_S: "<< val_out[1] <<endl;
//   exit(0);
 } // Get_Heat(

 


// ======================================================================
// utilities for main program
// ======================================================================

int main(int argc, char* argv[])
{
  TDomain *Domain = new TDomain();
  TDatabase *Database = new TDatabase();
  TFEDatabase2D *FEDatabase = new TFEDatabase2D();
  TCollection *coll, *NSE_coll, *Solid_coll, *mortarcoll = NULL;
  TBaseCell *Me, *cell, **Free_Cells, **NSE_Cells, **Solid_Cells, **BD_Cells, **SolidNeibCells;
  TGridCell **DelCell;
  TFESpace2D *velocity_space, *pressure_space, *streamfunction_space, *convolution_space, *fesps, *heatfunc_space;
  TFESpace2D  *pressure_space_output;
  TFESpace2D *Grid_space, *vorticity_space, *thermal_space,*grid_space;
  TFESpace2D *Grid_space_S, *Grid_space_NSE;
  
  TOutput2D *Output, *OutputAll;
  TFEVectFunct2D *RefGridPos, *AuxGridPos, *GridPos;
  TFEVectFunct2D *RefGridPos_S, *AuxGridPos_S, *GridPos_S;  
  TFEFunction2D *fefct[5];
  TFESpace2D *fesp[4], *ferhs_T[3], *ferhs[2];
  TAuxParam2D *aux;
  TDiscreteForm2D *DiscreteFormMatrixT_MRhs;
  TDiscreteForm2D *DiscreteForm;
  TMatrix2D *MATRICES[4];
  TSquareMatrix2D *SQMATRICES[8], *SQMATRICES_GRID[4];
  TSquareMatrix2D *SQMATRICES_HEAT[2], *SQMATRICES_HEATFUNC[1];
  TDiscreteForm2D *DiscreteFormGalerkin;
  TDiscreteForm2D *DiscreteFormNLGalerkin;
  TDiscreteForm2D *DiscreteFormGrid, *DiscreteFormHeat, *DiscreteFormHeatfunc  ,*DiscreteFormHeat_SUPG;
  TSquareMatrix **sqmatrices = (TSquareMatrix **)SQMATRICES;
  TMatrix **matrices = (TMatrix **)MATRICES;
  MatVecProc *MatVect;
  DefectProc *Defect;

  TFESpace2D **FESpaces_All = new TFESpace2D *[6];      
  TFEFunction2D **FEFunctions_All = new TFEFunction2D *[9], *P_output;    
  TFEVectFunct2D **FEVectFuncts_All = new TFEVectFunct2D*[3], *Velo_output;
  TStructure2D **Structure_All = new TStructure2D *[2];
  TSquareStructure2D **SquareStructure_All = new TSquareStructure2D *[4];
  TSquareStructure2D *SquareStructure_S, *SquareStructure_NSE;   
  TSquareMatrix2D **SqMat_All = new TSquareMatrix2D *[13];
  TSquareMatrix2D **GridSqMat_NSE = new TSquareMatrix2D *[4]; 
  TSquareMatrix2D **GridSqMat_S = new TSquareMatrix2D *[4];  
  TMatrix2D **Mat_All = new TMatrix2D *[4];
  TFEVectFunct2D *GridVect_NSE, *GridVect_S;
  TFEFunction2D *GridG1_NSE, *GridG2_NSE, *GridG1_S, *GridG2_S;
  FE2D *fes, *fesheat;

  double total_time,*Coordinates;
  double  t, teta, dt,x,y,gamma, Tgamma, tx,ty,sx,sy, R_Theta[3];;
  double left, right, top, bottom,T_a, T_b;
  double x0, y0,x1,y1,hi, residual, impuls_residual, oldresidual, solver_time;
  double *oldsol_T, HeatFlux, OldHeatFlux, L2HeatFlux;
  double end_time, t1, t2, t4, t3;
  double *B, *defect, *heat_defect, *Heat_B, *RHSs_Heat[1], *RHSs_Heatfunc[1];
  double *RHSs[3], *refpos, *auxpos, *pos, *refpos_S, *auxpos_S, *pos_S;
  double  TX[2], TY[2], solver_time_curr;
  double SLPX, SLPY, *Entries[4], tau, oldtau, limit, *sol_output, Params[10], InitVolume, CurrVolume;  
  double Lx, Ly, Rx, Ry,  x2, y2, x3, y3, x4, y4, fh, fhlimit, fhtot, fhmin, fhmax;
  double *Angle = new double[2], **FreePts = new double *[2];  
  double **Sol_All = new double *[7];
  double **Rhs_All = new double *[4],  *Entries_S[4];
  double *GridSol_NSE,  *GridRhs_NSE, *GridSol_S, *GridRhs_S, *tmp_GridSol_NSE;
  double *Flux_A, *Flux_M,Flux_F, *tmp_Gridd_NSE;
  double lpcoeff, lpexponent, Initial_T, Initial_T_L, Initial_T_S;
  double T_IntfaceMinMax[2],  T_LSIntfaceMinMax[2];

  int N, ret,N_RootCells,N_Cells,N_Joints, N_Vertices,N_G, N_NSE_Cells, N_NonNSE_Cells;
  int N_Solid_Cells, *N_GidDofs, *N_GridActive, *N_GridBdDofs;
  int N_SlipBound_Vert,  N_FreeBound_Vert,  N_AxialBound_Vert,N_Interf_Vertices;
  int In_Index,CurrComp,CurrVertex, img=1, RemeshImg=1;
  int ID,CurrNeib,Neib[2], N_SquareMatrices, N_RectMatrices;
  int a,b,i,X,Y,j,k,l,len1, len2,Neighb_tmp,Indextemp;
  int *PartMarker, *Triangles,comp, *NSE_GlobalCllNo;
  int *PointNeighb,maxEpV = 0, Max_It, N_U_output, N_P_output;
  int  N_U, N_P,N_Unknowns,N_pressureDOF,N_Rhs,N_FESpaces;
  int  N_Hori1, N_Hori2,N_Verti,N_Boundary_Vert,N_thermalDOF,N_thermalActive,N_thermalNonActive,N_heatfuncActive,N_heatfuncDOF,N_heatfuncNonActive;
  int velocity_space_code, pressure_space_code;
  int m,m0,m1,m2,m3,m4,m5,m6, N_Active, N_LinIterCurr, N_LinIter;
  int ORDER,VSP, N_GActive, N_GBoundaryNodes, N_SubSteps, very_first_time=1 ;
  int **IsoCellEdgeNos, *GridKCol, *GridRowPtr, *RowPtr, last_sq, N_SolidNeibCells;
  int *GridKCol_S, *GridRowPtr_S, N_G_S, N_GActive_S, N_GBoundaryNodes_S;
  int  *JointDOF, N_DOF, dof, *DOF, *GlobalNumbers, *BeginIndex, N_Remesh=0, N_ReParam=0;
  int N_FluxDof, *FluxDof, OrderDiff, *GridGlobalNumbers, *GridBeginIndex, GlobalCellNo;
  int RootVertNo, N_RootVert, N_BData=1, Max_GridLength;
  
  char *PRM, *GEO, *PsBaseName, *VtkBaseName, *Gnubasename;
  char ReadinDat[] = "readin.dat";
  char TString[] = "T";
  char HString[] = "H";
  char NameString[]  = "name";
  char UString[] = "U";
  char PString[] = "P";
  char WString[] = "W";
  char VortString[] = "Vort";
  char DivString[] = "Div";  
  char HeatString[] = "H";
  const char vtkdir[] = "VTK";
  const char BDdir[] = "BDData";
  
  bool remeshed=FALSE, reparam = FALSE, ComputeFlux = TRUE, FOUND=FALSE;;

  TBoundPart *BoundPart;
  TBoundComp *BoundComp;
  TBdLine *UpdateSlipBound, *UpdateAxialBound, *UpdateSlipBoundSolid;
  TBdLine *SolidBound1, *SolidBound2, *SolidBound3;
  TBdCircle *UpdateFreeBound;
  TBaseCell **CellTree;
  TVertex *vert, **VertexDel,**NewVertices;
  TJoint *Joint, *DelJoint, *JointIntface;
  TBoundEdge *Solid_Joint;
  TInterfaceJoint *Interface_Joint;
  TBoundEdge ***Bound_Joint = new TBoundEdge**[3];
  TIsoBoundEdge **Free_Joint, *IsoJoint;
  TVertex ***MovBoundVert = new TVertex**[4];
  TVertex *temp_Mov;
  TBoundEdge *tempSlip_Joint;
  FE2D FeId;
  TFEDesc2D *FeDesc;

  BoundCondFunct2D *BoundaryConditions[2], *BoundaryConditionsAuxProblem[3];
  BoundValueFunct2D *BoundValues[2], *BoundValuesAuxProblem[3];

  BoundCondFunct2D *HeatBoundaryConditions[1];
  BoundValueFunct2D *HeatBoundValues[1];

  BoundCondFunct2D *HeatfuncBoundaryConditions[1];
  BoundValueFunct2D *HeatfuncBoundValues[1];
  
  
  BoundCondFunct2D *GridBoundaryConditions[1];
  BoundValueFunct2D *GridBoundValues[1];


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
  if( TDatabase::ParamDB->DISCTYPE==2 )
  {
    OutPut("SDFEM does not work!" << endl);
    Error("SDFEM does not work!" << endl);
    exit(4711);
  }
  if(TDatabase::ParamDB->DISCTYPE==5)
  {
    OutPut("DISCTYPE 5 NOT IMPLEMENTED!" << endl);
    Error("DISCTYPE 5 NOT IMPLEMENTED!" << endl);
    exit(4711);
  }

  PRM = TDatabase::ParamDB->BNDFILE;
  GEO = TDatabase::ParamDB->GEOFILE;
  
  PsBaseName = TDatabase::ParamDB->BASENAME;
  VtkBaseName = TDatabase::ParamDB->VTKBASENAME;
  Gnubasename = TDatabase::ParamDB->BASENAME;
   
    
  mkdir(vtkdir, 0777);
  mkdir(BDdir, 0777);  
  

  
  lpcoeff = TDatabase::ParamDB->LP_STREAMLINE_COEFF;  
  lpexponent = TDatabase::ParamDB->LP_STREAMLINE_EXPONENT;
  OrderDiff = TDatabase::ParamDB->LP_STREAMLINE_ORDER_DIFFERENCE;
#define __AXIAL3D__ 


//======================================================================
// read boundary parameterization and initialize coarse grid
//======================================================================

   Domain->Init(PRM, GEO);

   
   // write grid into an Postscript file
   os.seekp(std::ios::beg);
   os << "Domain_old" << ".ps" << ends;
   Domain->PS(os.str().c_str(),It_Finest,0);

   boolean AllowEdgeRef = (boolean) TDatabase::ParamDB->MESHGEN_ALLOW_EDGE_REF;
   
//======================================================================
// Triangular for grid generation begin
//======================================================================
  BoundPart = Domain->GetBdPart(0);
  UpdateSlipBound = (TBdLine*)BoundPart->GetBdComp(0);
  UpdateFreeBound = (TBdCircle*)BoundPart->GetBdComp(1);
  UpdateAxialBound = (TBdLine*)BoundPart->GetBdComp(2);  
  UpdateSlipBoundSolid = (TBdLine*)BoundPart->GetBdComp(6);
  
  SolidBound1 = (TBdLine*)BoundPart->GetBdComp(3);
  SolidBound2 = (TBdLine*)BoundPart->GetBdComp(4);
  SolidBound3 = (TBdLine*)BoundPart->GetBdComp(5);  
  
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
  double refX, hE, phi, r = TDatabase::ParamDB->P4;

  int N_refX, *N_MovVert = new int[4];
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
//   N_SlipBound_Vert = int(TDatabase::ParamDB->P1);     // Initially only two points on solid bound (except end point)
  N_SlipBound_Vert = 2;
  
//   cout<< N_SlipBound_Vert << " " <<  TDatabase::ParamDB->P1 << endl;
//   exit(0);
  
#ifdef __ConvergenceStudy__
//   N_SlipBound_Vert = (int)N_FreeBound_Vert/Pi;  
  N_AxialBound_Vert = (int)N_FreeBound_Vert/Pi;  
    
  double Xi[5] = {0., 0.,  4., 4., 1.};
  double Yi[5] = {0.,-2., -2., 0., 0.};
  
//   double Xi[5] = {0., 0.,  2., 2., 1.};
//   double Yi[5] = {0.,-1., -1., 0., 0.};  
#else
  double Xi[5] = {0., 0.,  8., 8., 0.2};
  double Yi[5] = {0.,-4., -4., 0., 0.};

//   double Xi[5] = {0., 0.,  29.45, 29.45, 0.2};
//   double Yi[5] = {0.,-0.775, -0.775, 0., 0.};
#endif
  
  N_Hori1 = 50;      // number of horozontal vertices
  N_Hori2 = 10;
  N_Verti = 5;
 
  N_Boundary_Vert = N_Hori1+N_Hori2+2*N_Verti;

  N_Interf_Vertices = N_FreeBound_Vert+N_SlipBound_Vert+N_AxialBound_Vert+N_Boundary_Vert;
  In.numberofpoints = N_Interf_Vertices-1;
  In.pointlist = new double[2*In.numberofpoints];
  In.pointmarkerlist = new int[In.numberofpoints];
  In.numberofpointattributes = 0;

  In.numberofsegments = In.numberofpoints+1;
  In.segmentlist = new int[2*In.numberofsegments];
  In.segmentmarkerlist = new int[In.numberofsegments];
  In.numberofholes = 0;
  In.holelist = NULL;
  In.numberofregions = 0;
  In.regionlist = NULL;

  In_Index = 0;
  CurrComp = 1;
  
  a=0.;
  b=1.;
  
#ifdef __ConvergenceStudy__
    b=0.;
    dt = 1./(double)N_SlipBound_Vert;   
#else   
   teta =-Pi/2.; // end degree value of freesurface
   dt = (Pi/2. - teta)/(N_SlipBound_Vert+N_FreeBound_Vert);
   t=teta;
#endif

  
  for(i=0;i<N_SlipBound_Vert;i++) // without last point
   {
#ifdef __ConvergenceStudy__
    x = dt*(double)i;
#else     
    x = a+r*cos(t);
#endif
    y = 0.;
    if(fabs(x)<1.e-12) x = 0.;
    
    In.pointlist[2*In_Index] = x;
    In.pointlist[2*In_Index+1] = y;
//     cout<<" x : "<< x << " y : "<< y<<endl;
    In.pointmarkerlist[In_Index] = CurrComp;
    In.segmentlist[2*In_Index] = In_Index;
    In.segmentlist[2*In_Index+1] = In_Index+1;
    In.segmentmarkerlist[In_Index] = CurrComp;
    In_Index++;
#ifndef __ConvergenceStudy__  
    t = teta + double(i+1.)*dt;
#endif    
   }
  
   CurrComp++;
    cout<<endl; 
     
//   exit(0);
  
#ifdef __ConvergenceStudy__
   teta =0.; // end degree value of freesurface
   dt = (Pi/2. - teta)/(double)N_FreeBound_Vert;
#endif     
    
  for(i=0;i<N_FreeBound_Vert;i++) // without last point
    {
    //  cout<<" teta : "<< teta <<endl;
#ifdef __ConvergenceStudy__
     t = teta + (double)i*dt;     
#else
     t = teta + (double)(N_SlipBound_Vert+i)*dt;     
#endif            
      x = a+r*cos(t);
      y = b+r*sin(t);
 
     if (fabs(x)<1.e-10) x = 0;
     if(i==0) 
       {
        y =0.;   SLPX = x;
        Xi[4]=x;
        Yi[4]=y;
        Indextemp = In_Index;    
       }
      else if(i==1) 
      {
       hE = sqrt( (In.pointlist[2*(In_Index-1)] - x)*(In.pointlist[2*(In_Index-1)] - x) 
               + (In.pointlist[2*(In_Index-1) +1] - y)*(In.pointlist[2*(In_Index-1)+1] - y));
      }
       
      In.pointlist[2*In_Index] = x;
      In.pointlist[2*In_Index+1] = y;
//        cout<<" x : "<< x << " y : "<< y<<endl;
      In.pointmarkerlist[In_Index] = CurrComp;
      In.segmentlist[2*In_Index] = In_Index;
      In.segmentlist[2*In_Index+1] = In_Index+1;
      In.segmentmarkerlist[In_Index] = CurrComp;
      In_Index++;
    }
#ifdef __ConvergenceStudy__
     t = teta + (double)i*dt;     
#else
     t = teta + (double)(N_SlipBound_Vert+i)*dt;     
#endif   
     
   OutPut("hE : " << hE<<endl);     
   
   
   
//   cout<<endl; 
   CurrComp++;
   
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
//    y = b+r*sin(t);
//    t = y;  
// 
//    for(i=0;i<N_refX;i++) // without last point
//     {
//       In.pointlist[2*In_Index] = x;
//       In.pointlist[2*In_Index+1] = y;
//       cout<<" x : "<< x << " y : "<< y<<endl;
//       In.pointmarkerlist[In_Index] = CurrComp;
//       In.segmentlist[2*In_Index] = In_Index;
//       In.segmentlist[2*In_Index+1] = In_Index+1;
//       In.segmentmarkerlist[In_Index] = CurrComp;
//       In_Index++;
//       y = t + (double)(i+1)*dt;
//      }  
//    
//    cout<<endl;
//    
//    dt= -((b+r*sin(t)) - refX) / (double)(N_AxialBound_Vert - N_refX);
//    x = 0.;
//    y = (b+r*sin(t)) - refX;
//    t = y;   
//    N = N_AxialBound_Vert - N_refX;
// 
//   for(i=0;i<N;i++) // without last point
//    {
//       if (fabs(y)<1e-8) y = 0.;
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

   y = (b+r*sin(t)) - In.pointlist[1];
   dt= -y / (double)(N_AxialBound_Vert);
   x = 0.;
   t = y;   
   N = N_AxialBound_Vert;

  for(i=0;i<N;i++) // without last point
   {
      if (fabs(y)<1e-8) y = 0.;
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
//   exit(0);
  
//    cout<<endl; 
   CurrComp++;
   dt= (Yi[1] - Yi[0])  / N_Verti;
   x = Xi[0];
   y = Yi[0];
   t = y;
   
  // (0,0) point is already defined
  In.segmentlist[2*In_Index] = 0;
  In.segmentlist[2*In_Index+1] = In_Index;
  In.segmentmarkerlist[In_Index] = CurrComp;
  In_Index++;  
   
   for(i=1;i<N_Verti;i++) // without last point
    {
      y = t + ((double)i)*dt;
      In.pointlist[2*(In_Index-1)] = x;
      In.pointlist[2*(In_Index-1)+1] = y;
//       cout<<" x : "<< x << " y : "<< y<<endl;
      In.pointmarkerlist[(In_Index-1)] = CurrComp;
      
      In.segmentlist[2*In_Index] = In_Index-1;
      In.segmentlist[2*In_Index+1] = In_Index;
      In.segmentmarkerlist[In_Index] = CurrComp;
      In_Index++;

    }


//    cout<<endl; 
   CurrComp++;
   dt= (Xi[2] - Xi[1])  / N_Hori2;
   x = Xi[1];
   y = Yi[1];
   t = x;

   for(i=0;i<N_Hori2;i++) // without last point
    {
      In.pointlist[2*(In_Index-1)] = x;
      In.pointlist[2*(In_Index-1)+1] = y;
//       cout<<" x : "<< x << " y : "<< y<<endl;
      In.pointmarkerlist[(In_Index-1)] = CurrComp;
      In.segmentlist[2*In_Index] = In_Index-1;
      In.segmentlist[2*In_Index+1] = In_Index;
      In.segmentmarkerlist[In_Index] = CurrComp;
      In_Index++;
      x = t + ((double)(i+1))*dt;
    }


//    cout<<endl; 
   CurrComp++;
   dt= (Yi[3] - Yi[2])  / N_Verti;
   x = Xi[2];
   y = Yi[2];
   t = y;

   for(i=0;i<N_Verti;i++) // without last point
    {
      In.pointlist[2*(In_Index-1)] = x;
      In.pointlist[2*(In_Index-1)+1] = y;
     
//       cout<<" x : "<< x << " y : "<< y<<endl;
      In.pointmarkerlist[(In_Index-1)] = CurrComp;
      In.segmentlist[2*In_Index] = In_Index-1;
      In.segmentlist[2*In_Index+1] = In_Index;
      In.segmentmarkerlist[In_Index] = CurrComp;
      In_Index++;
      y = t + ((double)(i+1))*dt;
    }

// #ifdef __ConvergenceStudy__
   refX = Xi[4] + 0.2;   
/*#else
   refX = Xi[4] + 0.2;  
#endif    */  

#ifdef  __ConvergenceStudy__
  hE =    fabs(Xi[4] - Xi[3])  / (double)N_Hori1;  
#endif    

   N_refX = (int)((refX - Xi[4])/hE);  
   
   if( (N_refX+10)>N_Hori1)
   {
    cout << "Increase  N_Hori1 points, N_refX:  "<<  N_refX <<endl;
    exit(0);
   }
 
//    cout<<" x : "<< Xi[4] << " refX : "<< refX << " hE "<< hE<<" N_refX "<< N_refX <<endl;   
//    exit(0);
   
//    cout<<endl; 
   CurrComp++;
   dt= (refX - Xi[3])  / (N_Hori1-N_refX);
   x = Xi[3];
   y = Yi[3];
   t = x;   

   for(i=0;i<(N_Hori1-N_refX);i++) // without last point
    {
      In.pointlist[2*(In_Index-1)] = x;
      In.pointlist[2*(In_Index-1)+1] = y;
//       cout<<" x : "<< x << " y : "<< y<<endl;
      In.pointmarkerlist[(In_Index-1)] = CurrComp;
      In.segmentlist[2*In_Index] = (In_Index-1);
      In.segmentlist[2*In_Index+1] = In_Index;
      In.segmentmarkerlist[In_Index] = CurrComp;
      In_Index++;
      x = t + ((double)(i+1))*dt;
    }
   cout<<endl;
   dt= (Xi[4] - refX)  / N_refX;
   x = refX;
   y = Yi[3];
   t = x;   

   for(i=0;i<N_refX;i++) // without last point
    {
      In.pointlist[2*(In_Index-1)] = x;
      In.pointlist[2*(In_Index-1)+1] = y;
//       cout<<" x : "<< x << " y : "<< y<<endl;
      In.pointmarkerlist[(In_Index-1)] = CurrComp;
      In.segmentlist[2*In_Index] = (In_Index-1);
      In.segmentlist[2*In_Index+1] = In_Index;
      In.segmentmarkerlist[In_Index] = CurrComp;
      In_Index++;
      x = t + ((double)(i+1))*dt;
    }

  In.segmentlist[2*(In_Index-1)+1] = Indextemp;
  
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

   delete [] VertexDel;
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

  UpdateSlipBound->SetParams(Xi[0], Yi[0], Xi[4]-Xi[0],Yi[4]-Yi[0]);
  UpdateSlipBoundSolid->SetParams(Xi[3], Yi[3], Xi[4]-Xi[3],Yi[4]-Yi[3]);
  SolidBound1->SetParams(Xi[0], Yi[0], Xi[1]-Xi[0],Yi[1]-Yi[0]);
  SolidBound2->SetParams(Xi[1], Yi[1], Xi[2]-Xi[1],Yi[2]-Yi[1]); 
  SolidBound3->SetParams(Xi[2], Yi[2], Xi[3]-Xi[2],Yi[3]-Yi[2]);

#ifdef __ConvergenceStudy__  
 // Free boundary xmid, ymid, radius_a, radius_b, start angle, end angle
 UpdateFreeBound ->SetParams(0.0, 0.0, 1.0, 1.0, 0., Pi/2.);    
#endif
 
//  OutPut("left: "<<left<<" right: "<<right<<" top: "<<top<<" bottom: "<<bottom<<endl);
  Domain->SetBoundBox(right-left,top-bottom);
  Domain->SetBoundBoxstart(left,bottom);

 // generate cells
   CellTree = new TBaseCell*[N_RootCells];

  for (i=0;i<N_RootCells;i++)
  {
    CellTree[i] = new TMacroCell(TDatabase::RefDescDB[Triangle], 0);

    CellTree[i]->SetVertex(0, NewVertices[Triangles[3*i    ]]);
    CellTree[i]->SetVertex(1, NewVertices[Triangles[3*i + 1]]);
    CellTree[i]->SetVertex(2, NewVertices[Triangles[3*i + 2]]);

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
//   delete [] PointNeighb;
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
//   if(Out.trianglelist!=NULL) {
//     free(Out.trianglelist); Out.trianglelist = NULL;}
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
  os << "Domain_0" << ".ps" << ends;
  Domain->PS(os.str().c_str(),It_Finest,0);

  coll=Domain->GetCollection(It_Finest, 0);

//   UpdateCellwith2BDs(coll);

//       // write grid into an Postscript file
//       os.seekp(std::ios::beg);
//       os << "Domain_1" << ".ps" << ends;
//       Domain->PS(os.str().c_str(),It_Finest,0);
// 
//       exit(0);
      
  
    // initialize iterators
  TDatabase::IteratorDB[It_EQ]->SetParam(Domain);
  TDatabase::IteratorDB[It_LE]->SetParam(Domain);
  TDatabase::IteratorDB[It_Finest]->SetParam(Domain);
  TDatabase::IteratorDB[It_Between]->SetParam(Domain);
  TDatabase::IteratorDB[It_OCAF]->SetParam(Domain);
  
  // refine grid up to the coarsest level
  for(i=0;i<TDatabase::ParamDB->UNIFORM_STEPS;i++)
   Domain->RegRefineAll(); 

  
  //for hetrogeneous disc for thermal space
  // use the info from mesh generator
  coll=Domain->GetCollection(It_Finest, 0);
  N_Cells = coll->GetN_Cells();

   // write grid into an Postscript file
    os.seekp(std::ios::beg);
    os << "Domain" << ".ps" << ends;
    Domain->PS(os.str().c_str(),It_Finest,0);
     


                                  
  InitializeDiscreteForms_Moving(DiscreteFormGalerkin, DiscreteFormNLGalerkin,
                                 DiscreteFormGrid, LinCoeffs, GridCoeffs);

  InitializeDiscreteForms_Moving(DiscreteFormHeat, DiscreteFormHeat_SUPG, HeatCoeffs);
 
  InitializeDiscreteForms_HeatLine(DiscreteFormHeatfunc, HeatfuncCoeffs);
 
  BoundaryConditions[0] = BoundCondition;
  BoundaryConditions[1] = BoundCondition;

  HeatBoundaryConditions[0] = HeatBoundCondition;

  BoundValues[0] = U1BoundValue;
  BoundValues[1] = U2BoundValue;

  HeatBoundValues[0] = TBoundValue;
  
#ifdef __HEATLINE__  
   HeatfuncBoundaryConditions[0] = HeatFuncBoundCondition;
   HeatfuncBoundValues[0] = HeatFuncBoundValue;
#endif  
    
  
  GridBoundaryConditions[0] = GridBoundCondition;
  GridBoundValues[0] = GridBoundValue;

  BoundaryConditionsAuxProblem[0] = BoundConditionAuxProblem;
  BoundaryConditionsAuxProblem[1] = BoundConditionAuxProblem;
  BoundaryConditionsAuxProblem[2] = BoundConditionAuxProblem;

  BoundValuesAuxProblem[0] = BoundValueAuxProblem;
  BoundValuesAuxProblem[1] = BoundValueAuxProblem;
  BoundValuesAuxProblem[2] = BoundValueAuxProblem;

//======================================================================
// construct all finite element spaces
//======================================================================

  ORDER = 0;
  VSP = TDatabase::ParamDB->VELOCITY_SPACE;

  if (abs(VSP) > 20)
   {ORDER = abs(VSP) - 20;}
  else if ( abs(VSP) > 10)
   {ORDER = abs(VSP) - 10;}
  else ORDER = abs(VSP);


  N_NSE_Cells = 0;

  for(i=0; i<N_Cells; i++)
   {
    cell = coll->GetCell(i);
    cell->SetGlobalCellNo(i);    
    
    N_Vertices = cell->GetN_Vertices();     

    t=0;
    for(j=0; j<N_Vertices; j++)
    {
     cell->GetVertex(j)->GetCoords(x, y);
     t +=y;
    }
 
    t /=(double)N_Vertices;
 
    if(t>0) // NSE cells
     {
      cell->SetPhase_ID(0); 
      N_NSE_Cells++;         
     }
    else
     { cell->SetPhase_ID(1); } 
   }
   
  NSE_Cells = new TBaseCell*[N_NSE_Cells];
  Solid_Cells = new TBaseCell*[N_Cells - N_NSE_Cells];  
  
  N_NSE_Cells = 0;
  N_Solid_Cells = 0;
  for(i=0; i<N_Cells; i++)
   {
    cell = coll->GetCell(i);
    ID = cell->GetPhase_ID();
 
    if(ID==0) // NSE cells
     {
      NSE_Cells[N_NSE_Cells] = cell;        
      N_NSE_Cells++;    
     }
    else
     {
      Solid_Cells[N_Solid_Cells] = cell;        
      N_Solid_Cells++;      
     }
    
   } 
 
  OutPut("Number of cells, NSE cells, Solid Cells: " << N_Cells << ",  " << N_NSE_Cells <<  
                                             ",  " << N_Solid_Cells <<endl);
// exit(0);

  NSE_coll = new TCollection(N_NSE_Cells, NSE_Cells); 
  Solid_coll = new TCollection(N_Solid_Cells, Solid_Cells); 

//======================================================================
// construct all finite element spaces
//======================================================================
  // get velocity and pressure spacess
  GetVelocityAndPressureSpace(NSE_coll,BoundCondition,
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
  
  GetVelocityAndPressureSpace(coll,BoundCondition_output,
                              mortarcoll, FESpaces_All[4],
                              pressure_space_output, &pressure_space_code,
                              TDatabase::ParamDB->VELOCITY_SPACE,
                              TDatabase::ParamDB->PRESSURE_SPACE);
   N_U_output = FESpaces_All[4]->GetN_DegreesOfFreedom();
  
   // mesh velocity space 
   N_GidDofs = new int[3];
   N_GridActive = new int[3];
   N_GridBdDofs = new int[3];

   //grid velo in entire domain will be 2nd order, 
   //whereas for calculations in each domain we use 1st order
   grid_space = new TFESpace2D(coll, NameString, TString, GridBoundCondition, 1, NULL);
   FESpaces_All[2] =  grid_space;     
   GridGlobalNumbers = FESpaces_All[2]->GetGlobalNumbers();
   GridBeginIndex = FESpaces_All[2]->GetBeginIndex();   
   
   Grid_space_NSE = new TFESpace2D(NSE_coll, NameString, TString, GridBoundCondition, 1, NULL);
   Grid_space_S = new TFESpace2D(Solid_coll, NameString, TString, GridBoundCondition, 1, NULL);    
   
   N_GidDofs[0] = FESpaces_All[2]->GetN_DegreesOfFreedom();;
   N_GidDofs[1] = Grid_space_NSE->GetN_DegreesOfFreedom(); 
   N_GidDofs[2] = Grid_space_S->GetN_DegreesOfFreedom(); 

   N_GridActive[0] = FESpaces_All[2]->GetActiveBound();
   N_GridActive[1] = Grid_space_NSE->GetActiveBound(); 
   N_GridActive[2] = Grid_space_S->GetActiveBound();    

   N_G = N_GidDofs[0];
   N_GActive = N_GridActive[0];
   N_GridBdDofs[0] = N_G - N_GActive;   
   N_GridBdDofs[1] = N_GidDofs[1] - N_GridActive[1];
   N_GridBdDofs[2] = N_GidDofs[2] - N_GridActive[2];
   
// thermal space
   if(TDatabase::ParamDB->ANSATZ_ORDER<100)
   {
   thermal_space = new TFESpace2D(coll, NameString, TString, HeatBoundCondition,
                                  TDatabase::ParamDB->ANSATZ_ORDER, NULL);
   }
   else                              
   {
    fes = new FE2D[N_Cells];
    GetHetroFEs(N_Cells, coll, Triangles, PointNeighb, maxEpV, fes);     
    thermal_space = new TFESpace2D(coll, NameString, TString, HeatBoundCondition,
                                  fes, NULL);    
    //delete [] fes;
   }
   
   FESpaces_All[3] =  thermal_space;  
   N_thermalDOF = FESpaces_All[3]->GetN_DegreesOfFreedom();
   N_thermalActive = FESpaces_All[3]->GetActiveBound();
   N_thermalNonActive = N_thermalDOF - N_thermalActive; 
   
#ifdef __HEATLINE__     
   // heatfunction space
   if(TDatabase::ParamDB->ANSATZ_ORDER<100)
   {
   heatfunc_space = new TFESpace2D(coll, NameString, HeatString, HeatFuncBoundCondition,
                                  TDatabase::ParamDB->ANSATZ_ORDER, NULL);
   }
   else                              
   {
    fesheat = new FE2D[N_Cells];
    GetHetroFEs(N_Cells, coll, Triangles, PointNeighb, maxEpV, fesheat);     
    heatfunc_space = new TFESpace2D(coll, NameString, HeatString, HeatFuncBoundCondition,
                                  fesheat, NULL);    
   }
   
   FESpaces_All[5] =  heatfunc_space;  
   N_heatfuncDOF = FESpaces_All[5]->GetN_DegreesOfFreedom();
   N_heatfuncActive = FESpaces_All[5]->GetActiveBound();
   N_heatfuncNonActive = N_heatfuncDOF - N_heatfuncActive; 
    
   cout<<"Thermal DOFs : "<<N_thermalDOF<<"\nHeat function DOFs : "<<N_heatfuncDOF<<"\n";
#endif

//======================================================================
// construct all finite element functions
//======================================================================
  N_Unknowns = 2*N_U + N_P;
  Sol_All[0] = new double[N_Unknowns];  
  Rhs_All[0] = new double[N_Unknowns];
  
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

  // velo in all domains
  Sol_All[3] = new double[2*N_U_output]; 
  memset(Sol_All[3], 0, 2*N_U_output*SizeOfDouble);
  FEVectFuncts_All[2] =  new TFEVectFunct2D(FESpaces_All[4], UString, UString, Sol_All[3], N_U_output, 2);
  
  
  //mesh velo
  Sol_All[1] = new double[2*N_GidDofs[0]];
  Rhs_All[1] = new double[2*N_GidDofs[0]];    
  memset(Sol_All[1], 0, 2*N_GidDofs[0]*SizeOfDouble);   
  FEVectFuncts_All[1]  = new TFEVectFunct2D(FESpaces_All[2], WString, WString, Sol_All[1], N_GidDofs[0], 2); 
  FEFunctions_All[3] = FEVectFuncts_All[1]->GetComponent(0);
  FEFunctions_All[4] = FEVectFuncts_All[1]->GetComponent(1);

  GridSol_NSE = new double[2*N_GidDofs[1]];
  GridRhs_NSE = new double[2*N_GidDofs[1]];    
  tmp_GridSol_NSE = new double[2*N_GidDofs[1]]; 
  
  memset(GridSol_NSE, 0, 2*N_GidDofs[1]*SizeOfDouble);   
  GridVect_NSE  = new TFEVectFunct2D(Grid_space_NSE, WString, WString, GridSol_NSE, N_GidDofs[1], 2);
  GridG1_NSE = GridVect_NSE->GetComponent(0);
  GridG2_NSE = GridVect_NSE->GetComponent(1);

  GridSol_S = new double[2*N_GidDofs[2]];
  GridRhs_S = new double[2*N_GidDofs[2]];     
  memset(GridSol_S, 0, 2*N_GidDofs[2]*SizeOfDouble);   
  GridVect_S  = new TFEVectFunct2D(Grid_space_S, WString, WString, GridSol_S, N_GidDofs[2], 2);
  GridG1_S = GridVect_S->GetComponent(0);
  GridG2_S = GridVect_S->GetComponent(1);

  Max_GridLength = N_GidDofs[2];  
  if(Max_GridLength<N_GidDofs[1])
   Max_GridLength = N_GidDofs[1];   
  tmp_Gridd_NSE = new double[2*Max_GridLength];   
  
//======================================================================
// thermal space finite element functions
//======================================================================
  Sol_All[2] = new double[N_thermalDOF];
  Rhs_All[2] = new double[N_thermalDOF];
  oldsol_T = new double[N_thermalDOF];
  heat_defect = new double[N_thermalDOF];
  Heat_B = new double[N_thermalDOF];  

  //vectors for flux calculations
  Flux_A = new double[N_thermalDOF];
  Flux_M = new double[N_thermalDOF];
  
  memset(Sol_All[2], 0, N_thermalDOF*SizeOfDouble);
  memset(oldsol_T, 0, N_thermalDOF*SizeOfDouble);
  memset(Rhs_All[2], 0, N_thermalDOF*SizeOfDouble);

  // thermal fefunction
  FEFunctions_All[5] = new TFEFunction2D(FESpaces_All[3], TString, TString, Sol_All[2], N_thermalDOF);
#ifdef __ENERGY__  
  FEFunctions_All[5]->Interpolate(InitialT);
#endif  
  memcpy(oldsol_T, Sol_All[2], N_thermalDOF*SizeOfDouble); 
  
//======================================================================
// vorticity and div finite element functions
//====================================================================== 
 if(TDatabase::ParamDB->REACTOR_P22>0)
  {  
   Sol_All[4] = new double[N_GidDofs[1]];  
   Sol_All[5] = new double[N_GidDofs[1]];  
  
   memset(Sol_All[4], 0, N_GidDofs[1]*SizeOfDouble);
   memset(Sol_All[5], 0, N_GidDofs[1]*SizeOfDouble);   
  
   // vorticity fefunction
   FEFunctions_All[6] = new TFEFunction2D(Grid_space_NSE, VortString, VortString, Sol_All[4], N_GidDofs[1]);
   // divergence fefunction
   FEFunctions_All[7] = new TFEFunction2D(Grid_space_NSE, DivString, DivString, Sol_All[5], N_GidDofs[1]);  
  }
  
  
//======================================================================
//heatfunction space finite element functions
//======================================================================
#ifdef __HEATLINE__  
  Sol_All[6] = new double[N_heatfuncDOF];
  Rhs_All[3] = new double[N_heatfuncDOF];
  
  memset(Sol_All[6], 0, N_heatfuncDOF*SizeOfDouble);
  memset(Rhs_All[3], 0, N_heatfuncDOF*SizeOfDouble);

  // heatfunction fefunction
  FEFunctions_All[8] = new TFEFunction2D(FESpaces_All[5], HString, HString, Sol_All[6], N_heatfuncDOF);
#endif
  
  
//======================================================================
// allocate memory for all matrices
//======================================================================
  Structure_All[0] = new TStructure2D(FESpaces_All[1], FESpaces_All[0]); // B
  Structure_All[1] = new TStructure2D(FESpaces_All[0], FESpaces_All[1]); // BT
  
  //velo
  SquareStructure_All[0] = new TSquareStructure2D(FESpaces_All[0]);  
  SquareStructure_All[0]->Sort();

  // grid 
  SquareStructure_All[1] = new TSquareStructure2D(FESpaces_All[2]); 
  SquareStructure_All[1]->Sort();
  
  SquareStructure_NSE = new TSquareStructure2D(Grid_space_NSE); 
  SquareStructure_NSE->Sort();
    
  SquareStructure_S = new TSquareStructure2D(Grid_space_S);   
  SquareStructure_S->Sort();   
  
  //thermal
  SquareStructure_All[2] = new TSquareStructure2D(FESpaces_All[3]);
  SquareStructure_All[2]->Sort();
  
  
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

  // for heat
  SqMat_All[10]  = new TSquareMatrix2D(SquareStructure_All[2]); // T_M
  SqMat_All[11] = new TSquareMatrix2D(SquareStructure_All[2]); // T_A  
  
  //heatfunction
#ifdef __HEATLINE__    
  SquareStructure_All[3] = new TSquareStructure2D(FESpaces_All[5]);
  SquareStructure_All[3]->Sort();
  
  SqMat_All[12]  = new TSquareMatrix2D(SquareStructure_All[3]); // Heatfunc_A
#endif  
  
  // for mesh
  GridSqMat_NSE[0] = new TSquareMatrix2D(SquareStructure_NSE); // G11  
  GridSqMat_NSE[1] = new TSquareMatrix2D(SquareStructure_NSE); // G11    
  GridSqMat_NSE[2] = new TSquareMatrix2D(SquareStructure_NSE); // G11  
  GridSqMat_NSE[3] = new TSquareMatrix2D(SquareStructure_NSE); // G11    
  
  GridSqMat_S[0] = new TSquareMatrix2D(SquareStructure_S); // G11  
  GridSqMat_S[1] = new TSquareMatrix2D(SquareStructure_S); // G11    
  GridSqMat_S[2] = new TSquareMatrix2D(SquareStructure_S); // G11  
  GridSqMat_S[3] = new TSquareMatrix2D(SquareStructure_S); // G11    
       
  IsoCellEdgeNos = new int *[2];
  
  OutPut(endl);   
  OutPut("dof velocity : "<< setw(10) << 2* N_U << endl);
  OutPut("dof pressure : "<< setw(10) << N_P << endl);
  OutPut("dof heat     : "<< setw(10) << N_thermalDOF  << endl);
  OutPut("dof NSe all  : "<<  setw(10) << N_Unknowns << endl);   
#ifdef __HEATLINE__    
  OutPut("dof heatfunction  : "<<  setw(10) << N_heatfuncDOF << endl); 
  OutPut("dof all      : "<<  setw(10) << N_Unknowns + N_thermalDOF + N_heatfuncDOF<< endl);   
#endif  
  OutPut(endl);      
// exit(0);

 //  ==================================================================================== 
 //collect moving boundary vertices details
  y = 0;
  GetMovingBoundData(NSE_coll, N_MovVert, Bound_Joint, MovBoundVert, Free_Joint,
                     Free_Cells, IsoCellEdgeNos, SLPX, y);

  GetSolidBoundData(Solid_coll, N_MovVert, Bound_Joint, MovBoundVert, BD_Cells,
                    N_SolidNeibCells, SolidNeibCells, FESpaces_All[3], 0,  N_FluxDof, FluxDof);
  
//  ====================================================================================  
// assemble matrix for grid moving - begin
//  ====================================================================================  
    fesp[0] = Grid_space_NSE;
    SQMATRICES_GRID[0] = GridSqMat_NSE[0];
    SQMATRICES_GRID[0]->Reset();
    SQMATRICES_GRID[1] = GridSqMat_NSE[1];
    SQMATRICES_GRID[1]->Reset();
    SQMATRICES_GRID[2] = GridSqMat_NSE[2];
    SQMATRICES_GRID[2]->Reset();
    SQMATRICES_GRID[3] = GridSqMat_NSE[3];
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
//     delete aux;   
   
    Entries[0] = GridSqMat_NSE[0]->GetEntries();
    Entries[1] = GridSqMat_NSE[1]->GetEntries();
    Entries[2] = GridSqMat_NSE[2]->GetEntries();
    Entries[3] = GridSqMat_NSE[3]->GetEntries();

    GridKCol = SquareStructure_NSE->GetKCol();
    GridRowPtr = SquareStructure_NSE->GetRowPtr();

    N_G   = N_GidDofs[1];
    N_GActive =  N_GridActive[1];
    N_GBoundaryNodes = N_G - N_GActive;
   
    // for Dirichlet rows in off-diagonal matrices
    memset(Entries[1] + GridRowPtr[N_GActive], 0, (GridRowPtr[N_G] - GridRowPtr[N_GActive])*SizeOfDouble);
    memset(Entries[2] + GridRowPtr[N_GActive], 0, (GridRowPtr[N_G] - GridRowPtr[N_GActive])*SizeOfDouble);

    refpos = new double[2*N_G];
    auxpos = new double[2*N_G];
    pos = new double[2*N_G];
    RefGridPos = new TFEVectFunct2D(Grid_space_NSE, WString, WString, refpos, N_G, 2);
    AuxGridPos = new TFEVectFunct2D(Grid_space_NSE, WString, WString, auxpos, N_G, 2);
    GridPos = new TFEVectFunct2D(Grid_space_NSE, WString, WString, pos, N_G, 2);

    RefGridPos->GridToData();
    AuxGridPos->GridToData();
    GridPos->GridToData();

   // now for solid surface  
    fesp[0] = Grid_space_S;
    SQMATRICES_GRID[0] = GridSqMat_S[0];
    SQMATRICES_GRID[0]->Reset();
    SQMATRICES_GRID[1] = GridSqMat_S[1];
    SQMATRICES_GRID[1]->Reset();
    SQMATRICES_GRID[2] = GridSqMat_S[2];
    SQMATRICES_GRID[2]->Reset();
    SQMATRICES_GRID[3] = GridSqMat_S[3];
    SQMATRICES_GRID[3]->Reset();
//     aux = new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);  
       
    Assemble2D(1, fesp,
             4, SQMATRICES_GRID,
             0, NULL,
             0, NULL, NULL,
             DiscreteFormGrid,
             GridBoundaryConditions,
             GridBoundValues,
             aux);
    delete aux;   
   
    Entries_S[0] = GridSqMat_S[0]->GetEntries();
    Entries_S[1] = GridSqMat_S[1]->GetEntries();
    Entries_S[2] = GridSqMat_S[2]->GetEntries();
    Entries_S[3] = GridSqMat_S[3]->GetEntries();

    GridKCol_S = SquareStructure_S->GetKCol();
    GridRowPtr_S = SquareStructure_S->GetRowPtr();

    N_G_S   = N_GidDofs[2];
    N_GActive_S =  N_GridActive[2];
    N_GBoundaryNodes_S = N_G_S - N_GActive_S;
   
    // for Dirichlet rows in off-diagonal matrices
    memset(Entries_S[1] + GridRowPtr_S[N_GActive_S], 0, 
           (GridRowPtr_S[N_G_S] - GridRowPtr_S[N_GActive_S])*SizeOfDouble);
    memset(Entries_S[2] + GridRowPtr_S[N_GActive_S], 0, 
           (GridRowPtr_S[N_G_S] - GridRowPtr_S[N_GActive_S])*SizeOfDouble);

    refpos_S = new double[2*N_G_S];
    auxpos_S = new double[2*N_G_S];
    pos_S = new double[2*N_G_S];
    RefGridPos_S = new TFEVectFunct2D(Grid_space_S, WString, WString, refpos_S, N_G_S, 2);
    AuxGridPos_S = new TFEVectFunct2D(Grid_space_S, WString, WString, auxpos_S, N_G_S, 2);
    GridPos_S = new TFEVectFunct2D(Grid_space_S, WString, WString, pos_S, N_G_S, 2);

    RefGridPos_S->GridToData();
    AuxGridPos_S->GridToData();
    GridPos_S->GridToData();   
  
  // prepare output (maxn_fespaces,  maxn_scalar,  maxn_vect, maxn_parameters, domain)
  // prepare output
   Output = new TOutput2D(1, 3, 1, 2, Domain);
   Output->AddFEVectFunct(FEVectFuncts_All[0]);
   Output->AddFEFunction(FEFunctions_All[2]);  
//    Output->AddFEVectFunct(GridVect_NSE); 
   if(TDatabase::ParamDB->REACTOR_P22>0)
    {  
     Output->AddFEFunction(FEFunctions_All[6]);      
     Output->AddFEFunction(FEFunctions_All[7]);      
    }   
   os.seekp(std::ios::beg);
   Output->AddParameter(TDatabase::TimeDB->CURRENTTIME,os.str().c_str());      
    
   OutputAll = new TOutput2D(1, 3, 1, 2, Domain);
//    OutputAll->AddFEVectFunct(FEVectFuncts_All[1]); 
//    OutputAll->AddFEVectFunct(FEVectFuncts_All[2]);    
   OutputAll->AddFEFunction(FEFunctions_All[5]); 
#ifdef __HEATLINE__      
   OutputAll->AddFEFunction(FEFunctions_All[8]);    
#endif   
   os.seekp(std::ios::beg);
   OutputAll->AddParameter(TDatabase::TimeDB->CURRENTTIME,os.str().c_str());     
    
   //check LPS cells
//    if(TDatabase::ParamDB->ANSATZ_ORDER>=100)
//    {
//    memset(Sol_All[1], 0, 2*N_GidDofs[0]*SizeOfDouble); 
//    
//    
//    for(i=0; i<N_Cells; i++)
//    {   
//     if( (fes[i]==C_P1_2D_T_A) || (fes[i]==C_P2_2D_T_A)  || (fes[i]==C_P3_2D_T_A)
//         || (fes[i]==C_P4_2D_T_A) || (fes[i]==C_P5_2D_T_A) )
//       continue;
//     
//     cell = coll->GetCell(i);    
//     FeId =  FESpaces_All[2]->GetFE2D(i, cell);
//     FeDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(FeId);
//  
//     N_DOF = FeDesc->GetN_DOF();
//     DOF = GridGlobalNumbers+GridBeginIndex[i];
// 
//     for(k=0;k<N_DOF;k++)
//       Sol_All[1][DOF[k]] = 1.;       
//    }   
//    
//   
//    }
//         GetHeatConvectionALEVect(FEVectFuncts_All, GridVect_NSE, GridVect_S);   
   
 
   
     if(TDatabase::ParamDB->WRITE_VTK)
       {

         if(TDatabase::ParamDB->REACTOR_P22>0)
         {  
          ComputeVorticityDivergence(FESpaces_All[0], FEFunctions_All[0], FEFunctions_All[1], Grid_space_NSE, Sol_All[4],Sol_All[5]);
         }   
      
        os.seekp(std::ios::beg);
        if(img<10) os << "VTK/"<<VtkBaseName<<".0000"<<img<<".vtk" << ends;
         else if(img<100) os << "VTK/"<<VtkBaseName<<".000"<<img<<".vtk" << ends;
         else if(img<1000) os <<"VTK/"<< VtkBaseName<<".00"<<img<<".vtk" << ends;
         else if(img<10000) os << "VTK/"<<VtkBaseName<<".0"<<img<<".vtk" << ends;
         else  os << "VTK/"<<VtkBaseName<<"."<<img<<".vtk" << ends;
        Output->WriteVtk(os.str().c_str());

        os.seekp(std::ios::beg);
         if(img<10) os << "VTK/"<<Gnubasename<<".0000"<<img<<".vtk" << ends;
         else if(img<100) os <<"VTK/"<< Gnubasename<<".000"<<img<<".vtk" << ends;
         else if(img<1000) os << "VTK/"<<Gnubasename<<".00"<<img<<".vtk" << ends;
         else if(img<10000) os << "VTK/"<<Gnubasename<<".0"<<img<<".vtk" << ends;
         else  os << "VTK/"<<Gnubasename<<"."<<img<<".vtk" << ends;
        OutputAll->WriteVtk(os.str().c_str());   
    
    
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
   
   
  {
   MovBoundVert[0][0]->GetCoords(Lx, Ly);
   MovBoundVert[2][0]->GetCoords(Rx, Ry);

//    MovBoundVert[0][N_MovVert[0]-1]->GetCoords(x1, y1);
//    MovBoundVert[2][1]->GetCoords(x2, y2);
//    MovBoundVert[2][2]->GetCoords(x3, y3);   
//    MovBoundVert[2][3]->GetCoords(x4, y4);
   
   cout << x2 <<  "  " << x3<<  "  " << x4 << endl;
   cout << y2 <<  "  " << y3<<  "  " << y4 << endl;   
   
   
   tx = x1-Rx;
   sx = x2-Rx;
   ty = y1-Ry;
   sy = y2-Ry;
   
   x2 = sqrt(sx*sx + sy*sy);
   sx /=x2;
   sy /=x2;
   
   x2 = sqrt(tx*tx + ty*ty);
   tx /=x2;
   ty /=x2; 
   
   R_Theta[0] = acos( (tx*sx+ty*sy))*(180/3.141592654);

   sx = x3-Rx;
   sy = y3-Ry;
   x2 = sqrt(sx*sx + sy*sy);
   sx /=x2;
   sy /=x2;
   
   
   R_Theta[1] = acos( (tx*sx+ty*sy))*(180./3.141592654);

   sx = ((x4))-Rx;
   sy = ((y4))-Ry;
   x2 = sqrt(sx*sx + sy*sy);
   sx /=x2;
   sy /=x2;
   
   R_Theta[2] = acos( (tx*sx+ty*sy))*(180./3.141592654); 
   
   MovBoundVert[1][0]->GetCoords(x1, y1); 
   
   if(!remeshed)
    OutPut(setw(20)<<"T, wd,AxialZ,Ucl,RAng 1,2,3: " << TDatabase::TimeDB->CURRENTTIME<<"   "<< Rx-Lx
                   <<"   "<< y1<<"   "<< Params[2]<<"   "<<R_Theta[0]<<"   "<<R_Theta[1]<<"   "<<R_Theta[2]<<endl);    

   Get_KE(FEVectFuncts_All[0], Params);  
   OutPut(setw(20)<<"T, Volume, Diff, Rel. Diff : "<< TDatabase::TimeDB->CURRENTTIME<<"   "<< Params[0]
                  <<"   "<< Params[0] - InitVolume<<"   "<< (Params[0] - InitVolume)/InitVolume << endl);
  }
   
   
   
   
   
   
   
   Getcellangle(FESpaces_All[2], Angle);   
   OutPut( "MinAngle : "<< Angle[0]<< "  MaxAngle : "<<Angle[1]<< endl);

   
   Get_Heat(FEFunctions_All[5], Params);    
   Initial_T_L = Params[0];
   Initial_T_S = Params[1];   
   
   Initial_T  = Params[4]; 
   //cout<< " Initial_T_L: "<< Initial_T_L <<  " Initial_T_S: "<< Initial_T_S <<endl;   
   OutPut(setw(25)<<"T, Initial_Heat, Initial_T_L, Initial_T_S "<<  TDatabase::TimeDB->CURRENTTIME<<"   "<<  Initial_T <<"   "<<  Initial_T_L <<"   "<<  Initial_T_S<<endl);
//    exit(0);

      cout << "TIME_DISC " << TDatabase::TimeDB->TIME_DISC << endl;

   
   
 TDatabase::ParamDB->P5 = 0; //move freesur with velo
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

  OldHeatFlux = 0.;
  L2HeatFlux = 0;
  
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

        if (very_first_time)
            oldtau=tau;

        // working rhs array for NSE
        memset(B, 0, N_Unknowns*SizeOfDouble);

     //======================================================================
     //   grid velocity
     //======================================================================    
     GridVelo_imping(Entries, tmp_GridSol_NSE, tmp_Gridd_NSE, GridRhs_NSE,
                     GridKCol, GridRowPtr,
                     GridPos, AuxGridPos,
                     FEVectFuncts_All[0], tau,
                     GridVect_NSE, MovBoundVert, N_MovVert,
                     Free_Cells, IsoCellEdgeNos, reparam, RefGridPos, 
                     Entries_S, GridSol_S, GridRhs_S,
                     GridKCol_S, GridRowPtr_S,
                     GridPos_S, RefGridPos_S, GridVect_S);
    // ============================================================================================
    //  Assemble Energry Eqn rhs before updating the velo and mesh velo 
    // ============================================================================================    
#ifdef __ENERGY__       
    GetHeatConvectionALEVect(FEVectFuncts_All, GridVect_NSE, GridVect_S);

    // assemble heat matrices
    SQMATRICES_HEAT[0] = SqMat_All[11]; // A
    SQMATRICES_HEAT[0]->Reset();
    SQMATRICES_HEAT[1] = SqMat_All[10]; // M
    SQMATRICES_HEAT[1]->Reset();
 
    fesp[0] = FESpaces_All[3];  // thermal space
    fesp[1] = FESpaces_All[4];  // velocity space in all domain
    fesp[2] = FESpaces_All[2] ;  // mesh velocity space in all domain  
    
    ferhs[0] = FESpaces_All[3]; // thermal space for rhs    
    
    fefct[0] = FEVectFuncts_All[2]->GetComponent(0); // u1 all
    fefct[1] = FEVectFuncts_All[2]->GetComponent(1); // u2 all
    fefct[2] = FEVectFuncts_All[1]->GetComponent(0); // w1 all
    fefct[3] = FEVectFuncts_All[1]->GetComponent(1); // w2 all
    
    RHSs_Heat[0] = Rhs_All[2];   
    memset(RHSs_Heat[0], 0, N_thermalDOF*SizeOfDouble); 
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
                2, SQMATRICES_HEAT,
                0, NULL,
                1, RHSs_Heat, ferhs,
                DiscreteFormHeat,
                HeatBoundaryConditions,
                HeatBoundValues,
                aux, 1);

    //vectors for flux calculations
    if(ComputeFlux)
      GetFluxVect(2, SQMATRICES_HEAT, Rhs_All[2], N_FluxDof, FluxDof, Flux_A, Flux_M, Flux_F);  
        
     Assemble2D(3, fesp,
                2, SQMATRICES_HEAT,
                0, NULL,
                1, RHSs_Heat, ferhs,
                DiscreteFormHeat,
                HeatBoundaryConditions,
                HeatBoundValues,
                aux, 0);

     delete aux;


    // freesurfint
    if(fabs(TDatabase::ParamDB->BI_NR)>1e-10)
     Heat_freeint_axial3D(SqMat_All[11], HeatBoundCondition, Rhs_All[2]);     

    // add LPS atabilization
    if(TDatabase::ParamDB->ANSATZ_ORDER>100 && lpcoeff!=0)
     AddALEStreamlineLPS(SqMat_All[11], 4, fefct, lpcoeff, lpexponent, OrderDiff); 
   
    memset(Heat_B, 0, N_thermalDOF*SizeOfDouble);
 
    //rhs
    Daxpy(N_thermalActive, tau, Rhs_All[2], Heat_B);
     
    // M = M + ( tau*TDatabase::TimeDB->THETA2) A 
    MatAdd(SqMat_All[10], SqMat_All[11], -tau*TDatabase::TimeDB->THETA2);   
 //    Tgamma = -tau*TDatabase::TimeDB->THETA2;
   
    memset(heat_defect, 0, N_thermalDOF*SizeOfDouble); 
    MatVectActive(SqMat_All[10],  Sol_All[2], heat_defect);
    Daxpy(N_thermalActive, 1, heat_defect,  Heat_B);

    // copy Dirichlet values from rhs into working array rhs
    memcpy(Heat_B+N_thermalActive, Rhs_All[2]+N_thermalActive,  N_thermalNonActive*SizeOfDouble);   
    
   
    if(ComputeFlux)
    {     
     /** weighted due to theta-scheme, added on 8 sep 2013 */
     for(j=0;j<N_thermalDOF;j++)
      heat_defect[j] = TDatabase::TimeDB->THETA2*oldsol_T[j];
       
      HeatFlux =  tau *( Ddot(N_thermalDOF, heat_defect, Flux_A) );  
    }
#endif    
    
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
        fesp[2] = Grid_space_NSE;

        fefct[0] = FEFunctions_All[0];
        fefct[1] = FEFunctions_All[1];
        fefct[2] = GridG1_NSE;
        fefct[3] = GridG2_NSE;
 
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
     
     FreeSurf_Axial3DHeat(SqMat_All[8], SqMat_All[9],  RHSs[0], RHSs[1], BoundCondition,
                          tau, FEFunctions_All[5], T_IntfaceMinMax, FEFunctions_All[0]->GetValues(), Params);    
          
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
       GridVelo_imping(Entries, tmp_GridSol_NSE, tmp_Gridd_NSE, GridRhs_NSE,
                       GridKCol, GridRowPtr,
                       GridPos, AuxGridPos,
                       FEVectFuncts_All[0], tau,
                       GridVect_NSE, MovBoundVert, N_MovVert,
                       Free_Cells, IsoCellEdgeNos, reparam, RefGridPos, 
                       Entries_S, GridSol_S, GridRhs_S,
                       GridKCol_S, GridRowPtr_S,
                       GridPos_S, RefGridPos_S, GridVect_S);
  
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
       fesp[2] = Grid_space_NSE;

       fefct[0] = FEFunctions_All[0];
       fefct[1] = FEFunctions_All[1];
       fefct[2] = GridG1_NSE;
       fefct[3] = GridG2_NSE;
       
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

#ifdef __ENERGY__
   //======================================================================
   // end NSE nonlinear iteration
   // solve the energy equation
   // first get grid velocity
   //======================================================================    
    GridVelo_imping(Entries, tmp_GridSol_NSE, tmp_Gridd_NSE, GridRhs_NSE,
                     GridKCol, GridRowPtr,
                     GridPos, AuxGridPos,
                     FEVectFuncts_All[0], tau,
                     GridVect_NSE, MovBoundVert, N_MovVert,
                     Free_Cells, IsoCellEdgeNos, reparam, RefGridPos, 
                     Entries_S, GridSol_S, GridRhs_S,
                     GridKCol_S, GridRowPtr_S,
                     GridPos_S, RefGridPos_S, GridVect_S);


    GetHeatConvectionALEVect(FEVectFuncts_All, GridVect_NSE, GridVect_S);

    // assemble heat matrices
    SQMATRICES_HEAT[0] = SqMat_All[11]; // A
    SQMATRICES_HEAT[0]->Reset();
    SQMATRICES_HEAT[1] = SqMat_All[10]; // M
    SQMATRICES_HEAT[1]->Reset();
 
    fesp[0] = FESpaces_All[3];  // thermal space
    fesp[1] = FESpaces_All[4];  // velocity space in all domain
    fesp[2] = FESpaces_All[2] ;  // mesh velocity space in all domain  
    
    ferhs[0] = FESpaces_All[3]; // thermal space for rhs    
    
    fefct[0] = FEVectFuncts_All[2]->GetComponent(0); // u1 all
    fefct[1] = FEVectFuncts_All[2]->GetComponent(1); // u2 all
    fefct[2] = FEVectFuncts_All[1]->GetComponent(0); // w1 all
    fefct[3] = FEVectFuncts_All[1]->GetComponent(1); // w2 all
    
    RHSs_Heat[0] = Rhs_All[2];   
    memset(RHSs_Heat[0], 0, N_thermalDOF*SizeOfDouble); 
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
                2, SQMATRICES_HEAT,
                0, NULL,
                1, RHSs_Heat, ferhs,
                DiscreteFormHeat,
                HeatBoundaryConditions,
                HeatBoundValues,
                aux, 1);

    //vectors for flux calculations
    if(ComputeFlux)
      GetFluxVect(2, SQMATRICES_HEAT, Rhs_All[2], N_FluxDof, FluxDof, Flux_A, Flux_M, Flux_F);  
        
     Assemble2D(3, fesp,
                2, SQMATRICES_HEAT,
                0, NULL,
                1, RHSs_Heat, ferhs,
                DiscreteFormHeat,
                HeatBoundaryConditions,
                HeatBoundValues,
                aux, 0);

     delete aux;


   // freesurfint
   if(fabs(TDatabase::ParamDB->BI_NR)>1e-10)     
    Heat_freeint_axial3D(SqMat_All[11], HeatBoundCondition, Rhs_All[2]);     

   // add LPS atabilization
   if(TDatabase::ParamDB->ANSATZ_ORDER>100 && lpcoeff!=0)
    AddALEStreamlineLPS(SqMat_All[11], 4, fefct, lpcoeff, lpexponent, OrderDiff); 
   
// //    memset(Heat_B, 0, N_thermalDOF*SizeOfDouble);
 
//    //rhs
// //    Daxpy(N_thermalActive, tau, Rhs_All[2], Heat_B);
     
   // M = M + ( tau*TDatabase::TimeDB->THETA2) A 
// //    MatAdd(SqMat_All[10], SqMat_All[11], -tau*TDatabase::TimeDB->THETA2);   
// //    Tgamma = -tau*TDatabase::TimeDB->THETA2;
   
// //    memset(heat_defect, 0, N_thermalDOF*SizeOfDouble); 
// //    MatVectActive(SqMat_All[10],  Sol_All[2], heat_defect);
// //    Daxpy(N_thermalActive, 1, heat_defect,  Heat_B);

// //    // copy Dirichlet values from rhs into working array rhs
// //    memcpy(Heat_B+N_thermalActive, Rhs_All[2]+N_thermalActive,  N_thermalNonActive*SizeOfDouble);   

   

   //=====================================================================
   // assembling of system matrix
   //=====================================================================
   MatAdd(SqMat_All[10], SqMat_All[11], tau*TDatabase::TimeDB->THETA1);     
   
   DirectSolver(SqMat_All[10], Heat_B, Sol_All[2]);   


    if(ComputeFlux)
    {     
     for(j=0;j<N_thermalDOF;j++)
       heat_defect[j] = Sol_All[2][j] - oldsol_T[j];
           
     HeatFlux +=  Ddot(N_thermalDOF, heat_defect, Flux_M);
     
     /** weighted due to theta-scheme, added on 6 sep 2013 */
     for(j=0;j<N_thermalDOF;j++)
      heat_defect[j] = TDatabase::TimeDB->THETA1*Sol_All[2][j];
       
     HeatFlux +=  tau *(Ddot(N_thermalDOF, heat_defect, Flux_A) - Flux_F);
     
     L2HeatFlux += HeatFlux;  
        
/*    for (k=0;k<N_FluxDof;k++)
     {
  
      FESpaces_All[3]->GetDOFPosition(FluxDof[k],x0, y0);
      
      if(fabs(y0)<1e-8 && fabs(oldsol_T[FluxDof[k]] -  Sol_All[2][FluxDof[k]])>1e-2 )
      cout<< " k " << FluxDof[k] <<" FluxDof " <<oldsol_T[FluxDof[k]] -  Sol_All[2][FluxDof[k]]  << " x " << x0<<" y " << y0<<endl;

     }  */ 
     
     if(m==1 || m%10==0)
     OutPut("Time, HeatFlux,  TotalHeatFlux: " <<TDatabase::TimeDB->CURRENTTIME << "  " << HeatFlux << "  " << L2HeatFlux <<endl); 
    }
      
   memcpy(oldsol_T, Sol_All[2], N_thermalDOF*SizeOfDouble);     

   //======================================================================
   // end of heat equation
   //======================================================================
#endif   
/*
     if(TDatabase::ParamDB->WRITE_VTK)
       { 
        os.seekp(std::ios::beg);
        if(img<10) os << "VTK/"<<VtkBaseName<<".0000"<<img<<".vtk" << ends;
         else if(img<100) os << "VTK/"<<VtkBaseName<<".000"<<img<<".vtk" << ends;
          else if(img<1000) os << "VTK/"<<VtkBaseName<<".00"<<img<<".vtk" << ends;
           else if(img<10000) os << "VTK/"<<VtkBaseName<<".0"<<img<<".vtk" << ends;
            else  os << "VTK/"<<VtkBaseName<<"."<<img<<".vtk" << ends;
        Output->WriteVtk(os.str().c_str());

        GetHeatConvectionALEVect(FEVectFuncts_All, GridVect_NSE, GridVect_S);

        os.seekp(std::ios::beg);
         if(img<10) os << "VTK/"<<Gnubasename<<".0000"<<img<<".vtk" << ends;
         else if(img<100) os << "VTK/"<<Gnubasename<<".000"<<img<<".vtk" << ends;
         else if(img<1000) os << "VTK/"<<Gnubasename<<".00"<<img<<".vtk" << ends;
         else if(img<10000) os << "VTK/"<<Gnubasename<<".0"<<img<<".vtk" << ends;
         else  os << "VTK/"<<Gnubasename<<"."<<img<<".vtk" << ends;
        OutputAll->WriteVtk(os.str().c_str());   

        img++;
       }  
       
   if( fabs(TDatabase::TimeDB->CURRENTTIME-3.000000e-03)<1e-8)
     exit(0);*/
   
#ifdef __HEATLINE__   
   //=====================================================================
   //  Heat function equation,  assemble heat function matrices
   //=====================================================================
 
    SQMATRICES_HEATFUNC[0] = SqMat_All[12]; // A
    SQMATRICES_HEATFUNC[0]->Reset();
 
    fesp[0] = FESpaces_All[5];  // heat function space
    fesp[1] = FESpaces_All[3];  // thermal space
    fesp[2] = FESpaces_All[4];  // velocity space in all domain
    fesp[3] = FESpaces_All[2] ;  // mesh velocity space in all domain  
 
    fefct[0] = FEFunctions_All[5]; // T
    fefct[1] = FEVectFuncts_All[2]->GetComponent(0); // u1 all
    fefct[2] = FEVectFuncts_All[2]->GetComponent(1); // u2 all
    fefct[3] = FEVectFuncts_All[1]->GetComponent(0); // w1 all
    fefct[4] = FEVectFuncts_All[1]->GetComponent(1); // w2 all
    
    // T, T_x, T_y, (u1-w1), (u2-w2), x   parameters are needed for assembling
    // fesp is taken from fefct in aux
    aux =  new TAuxParam2D(MovingTNSN_FESpaces_Axial3D_HeatLine, MovingTNSN_Fct_Axial3D_HeatLine,
                           MovingTNSN_ParamFct_Axial3D_HeatLine,
                           MovingTNSN_FEValues_Axial3D_HeatLine,
                           fesp+1, fefct,
                           MovingTNSFct_Axial3D_HeatLine,
                           MovingTNSFEFctIndex_Axial3D_HeatLine,
                           MovingTNSFEMultiIndex_Axial3D_HeatLine,
                           MovingTNSN_Params_Axial3D_HeatLine, MovingTNSBeginParam_Axial3D_HeatLine);
    
    ferhs[0] = FESpaces_All[5]; // heat function space for rhs    
    
    RHSs_Heatfunc[0] = Rhs_All[3];   
    memset(RHSs_Heatfunc[0], 0, N_heatfuncDOF*SizeOfDouble); 

 
     Assemble2D(1, fesp,
                1, SQMATRICES_HEATFUNC,
                0, NULL,
                1, RHSs_Heatfunc, ferhs,
                DiscreteFormHeatfunc,
                HeatfuncBoundaryConditions,
                HeatfuncBoundValues,
                aux);
    
     delete aux;
    
   cout << "Solve Heat Line1 : " << Ddot(N_heatfuncDOF, Rhs_All[3], Rhs_All[3] ) << endl;
   DirectSolver(SqMat_All[12], RHSs_Heatfunc[0], Sol_All[6]);   
   
   cout << "Solve Heat Line2 :  " << Ddot(N_heatfuncDOF, Sol_All[6], Sol_All[6] ) << endl;
   
#endif
   
//     if((m % (int)(TDatabase::TimeDB->STEPS_PER_IMAGE) ) == 0   || m==1 )
//   {
//      if(TDatabase::ParamDB->WRITE_VTK)
//        { 
//         if(TDatabase::ParamDB->REACTOR_P22>0)
//          {  
//          ComputeVorticityDivergence(FESpaces_All[0], FEFunctions_All[0], FEFunctions_All[1], Grid_space_NSE, Sol_All[4],Sol_All[5]);
//          }  
//  
//         os.seekp(std::ios::beg);
//         if(img<10) os << "VTK/"<<VtkBaseName<<".0000"<<img<<".vtk" << ends;
//          else if(img<100) os << "VTK/"<<VtkBaseName<<".000"<<img<<".vtk" << ends;
//           else if(img<1000) os << "VTK/"<<VtkBaseName<<".00"<<img<<".vtk" << ends;
//            else if(img<10000) os << "VTK/"<<VtkBaseName<<".0"<<img<<".vtk" << ends;
//             else  os << "VTK/"<<VtkBaseName<<"."<<img<<".vtk" << ends;
//         Output->WriteVtk(os.str().c_str());
// 
//         GetHeatConvectionALEVect(FEVectFuncts_All, GridVect_NSE, GridVect_S);
// 
//         os.seekp(std::ios::beg);
//          if(img<10) os << "VTK/"<<Gnubasename<<".0000"<<img<<".vtk" << ends;
//          else if(img<100) os << "VTK/"<<Gnubasename<<".000"<<img<<".vtk" << ends;
//          else if(img<1000) os << "VTK/"<<Gnubasename<<".00"<<img<<".vtk" << ends;
//          else if(img<10000) os << "VTK/"<<Gnubasename<<".0"<<img<<".vtk" << ends;
//          else  os << "VTK/"<<Gnubasename<<"."<<img<<".vtk" << ends;
//         OutputAll->WriteVtk(os.str().c_str());   
// 
//         img++;
//        }       
//    }
   
// exit(0);

   //======================================================================
   // move the grid
   //======================================================================    
   GridPos_S->GridToData();
   MovBoundVert[2][0]->GetCoords(Rx, Ry);
   
   MoveGrid_imping(Entries, tmp_GridSol_NSE, tmp_Gridd_NSE, GridRhs_NSE,
                  GridKCol, GridRowPtr,
                  GridPos, FEVectFuncts_All[0], tau,
                  AuxGridPos, 
                  MovBoundVert, N_MovVert,
                  Free_Cells, IsoCellEdgeNos, reparam,
                  N_ReParam, FEFunctions_All[5]);

    // reparam the remaining solid vert and move the solid domain   
    MovBoundVert[2][0]->GetCoords(tx, Ry);  
    N=N_MovVert[3]-1;  
    hi = (tx-Rx)/(double)N;
    
    for(i=0;i<N;i++)
    {
     MovBoundVert[3][N_MovVert[3]-1 - i]->GetCoords(x, y);     
//      cout << "reparam x " << x <<"new x " << x +((double)(N-i))*hi <<endl;
     x += ((double)(N-i))*hi;
     MovBoundVert[3][N_MovVert[3]-1 - i]->SetCoords(x, y);   
    }
 
//  exit(0);
 
    // now move solid domain according to the slip BD movement
    RefGridPos_S->GridToData();
    GridPos_S->DataToGrid();
  
    memset(GridRhs_S, 0, 2*N_GidDofs[2]*SizeOfDouble);
    memcpy(Rhs_All[1], refpos_S, 2*N_GidDofs[2]*SizeOfDouble);
    Daxpy(2*N_GidDofs[2], -1, pos_S, Rhs_All[1]);
    memcpy(GridRhs_S + N_GridActive[2], Rhs_All[1]+N_GridActive[2],
           N_GridBdDofs[2]*SizeOfDouble);   
    memcpy(GridRhs_S + (N_GidDofs[2] + N_GridActive[2]),
           Rhs_All[1]+ (N_GidDofs[2] + N_GridActive[2]), N_GridBdDofs[2]*SizeOfDouble);     
   
    memset(GridSol_S, 0 , 2*N_GidDofs[2]*SizeOfDouble);   
    memcpy(GridSol_S + N_GridActive[2], Rhs_All[1]+N_GridActive[2],
           N_GridBdDofs[2]*SizeOfDouble);    
    memcpy(GridSol_S + (N_GidDofs[2] + N_GridActive[2]),
           Rhs_All[1]+ (N_GidDofs[2] + N_GridActive[2]), N_GridBdDofs[2]*SizeOfDouble);    
     
    SolveGridEquation(Entries_S, GridSol_S, GridRhs_S, GridKCol_S, GridRowPtr_S, N_GidDofs[2]); 
    Daxpy(2*N_GidDofs[2], 1, GridSol_S, pos_S);
  
    GridPos_S->DataToGrid();
 
    // Updating solid boundary 
    for(k=0;k<N_MovVert[3];k++)
     if(k==N_MovVert[3]-1)
      { Bound_Joint[2][k]->UpdateParameters(MovBoundVert[3][k], MovBoundVert[2][0]); }
     else
      { Bound_Joint[2][k]->UpdateParameters(MovBoundVert[3][k], MovBoundVert[3][k+1]); }
      
  // =======================================================================================
  // move solid BD end
  // =======================================================================================

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

    
    MovBoundVert[3][N_MovVert[3]-1]->GetCoords(Rx, Ry); // needed for reparam   
    MovBoundVert[3][0]->GetCoords(x, y);
    UpdateSlipBoundSolid->SetParams(x, y, SLPX-x, SLPY);    
    
    if(Domain->GetBdPart(0)->GetBdComp(0)->GetTofXY(
       MovBoundVert[2][0]->GetX(), MovBoundVert[2][0]->GetY(), T_a) ||
       Domain->GetBdPart(0)->GetBdComp(0)->GetTofXY(
       MovBoundVert[2][1]->GetX(),MovBoundVert[2][1]->GetY(), T_b))
     {
          cerr<<"Error: could not set parameter values"<<endl;
          OutPut(NewVertices[a]<<endl);
          OutPut(NewVertices[b]<<endl);
          cout << " CurrComp " << CurrComp <<endl;
        //  exit(0);
       }   

    cell = BD_Cells[N_MovVert[3]-1]; 
    Me = Free_Cells[0];        
    GlobalCellNo = Me->GetGlobalCellNo();  
    
    //find the root vertex number of the new solid vert
    for(j=0;j<3;j++)
     {
      if(Me->GetVertex(j) == MovBoundVert[2][1])
       { 
        RootVertNo = Triangles[3*GlobalCellNo + j];
        break;
       } // if(temp_Mov== MovB
     }//for(j=0;j<3;j++)
    
//     cout << " RootVertNo " << RootVertNo << endl;   
//     exit(0);

    
    Joint = Me->GetJoint(IsoCellEdgeNos[1][0]);
    DelJoint = Bound_Joint[2][N_MovVert[3]-1]; 
    
    //remove the isopoints in the joint
    IsoJoint = (TIsoBoundEdge *)Joint;
    IsoJoint->DeleteVertices();

   //find the local index of the solid joint vert 
   
//    FOUND = FALSE;
   for(i=0;i<N_SolidNeibCells;i++) 
   {
    GlobalCellNo = SolidNeibCells[i]->GetGlobalCellNo();
     
    for(j=0;j<3;j++)
     {
      temp_Mov = SolidNeibCells[i]->GetVertex(j);
      if(temp_Mov == MovBoundVert[3][N_MovVert[3]-1])
       {
//         FOUND = TRUE;
        Triangles[3*GlobalCellNo + j] = RootVertNo;
        SolidNeibCells[i]->SetVertex(j, MovBoundVert[2][1]);       
       }
     }// for(j=0;j<3;j++)     
   }//  for(i=0;i<N_SolidNeibCells;i++)

//     if(FOUND) 
//      delete temp_Mov;

    //find the local index of the solid joint
     for(j=0;j<3;j++)  
      if(cell->GetJoint(j) == DelJoint)
        break;
     
     if(j==3)
      { cout << "Error Joint index " << j << endl; exit(0); }
         
//      delete Joint; // MAC64 warning
//      delete DelJoint; // MAC64 warning

     Joint = new TInterfaceJoint(Domain->GetBdPart(0)->GetBdComp(0), T_a, T_b, Me, cell);  
      
     Me->SetJoint(IsoCellEdgeNos[1][0], Joint);
     cell->SetJoint(j, Joint);
     
    ((TInterfaceJoint *) Joint)->CheckOrientation();
      
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
    delete [] Bound_Joint[2];

    delete [] MovBoundVert[0];
    delete [] MovBoundVert[1];   
    delete [] MovBoundVert[2];  
    delete [] MovBoundVert[3]; 
    
    
    delete [] Free_Joint;    
    delete [] Free_Cells;  
    delete [] IsoCellEdgeNos[0];  
    delete [] IsoCellEdgeNos[1]; 
         
    delete []  BD_Cells;  
    delete []  SolidNeibCells;
    
    delete [] FluxDof;
    
    
//  ======================================================   
    N_RootVert = 0;
    for (i=0;i<3*N_Cells;i++)
     if (Triangles[i] > N_RootVert) N_RootVert = Triangles[i];
    
    N_RootVert++; //since array starts from 0  
    cout << "N_RootVert " << N_RootVert << endl;  
    
    delete [] PointNeighb;    
    PointNeighb = new int[N_RootVert];
    memset(PointNeighb, 0, N_RootVert *SizeOfInt);
   
    for (i=0;i<3*N_Cells;i++)
     PointNeighb[Triangles[i]]++;

    maxEpV = 0;
    for(i=0;i<N_RootVert;i++)
     if (PointNeighb[i] > maxEpV) maxEpV = PointNeighb[i];

    delete [] PointNeighb;
    PointNeighb = new int[++maxEpV * N_RootVert];  
    memset(PointNeighb, 0, maxEpV * N_RootVert *SizeOfInt);

    // first colomn contains the number of following elements
    // for every point at first column we set the number of neighbour points
    // at further columns we set the index of corresponding cells
    for(i=0;i<3*N_Cells;i++)
    {
     j = Triangles[i]*maxEpV;
     PointNeighb[j]++;
     PointNeighb[j + PointNeighb[j]] = i / 3;
    }
    
// exit(0);
//  ======================================================       
    /** mesh velo space_all added on 21.08.13 by sasi */
    BDChangeFERemap_All(FESpaces_All, FEVectFuncts_All, FEFunctions_All, N_GidDofs, N_GridActive, Sol_All, Rhs_All, Triangles, PointNeighb, maxEpV);
    
    BDChangeFERemap(FESpaces_All[3], FEFunctions_All[5], Sol_All[2],
                    SquareStructure_All[2], SqMat_All[10], SqMat_All[11],
                    Triangles, PointNeighb, maxEpV);
    
    fes = FESpaces_All[3]->GetAllElements();
    
    N_thermalDOF = FESpaces_All[3]->GetN_DegreesOfFreedom();
    N_thermalActive = FESpaces_All[3]->GetActiveBound();
    N_thermalNonActive = N_thermalDOF - N_thermalActive;
    OutPut("N_ThermalDOF    : "<< setw(10) << N_thermalDOF  << endl);    

#ifdef __HEATLINE__   
   BDChangeFERemap(FESpaces_All[5], FEFunctions_All[8], Sol_All[6],
                   SquareStructure_All[3], SqMat_All[12], Triangles, PointNeighb, maxEpV);    
    
   N_heatfuncDOF = FESpaces_All[5]->GetN_DegreesOfFreedom();
   N_heatfuncActive = FESpaces_All[5]->GetActiveBound();
   N_heatfuncNonActive = N_heatfuncDOF - N_heatfuncActive; 
    
   cout<<"Thermal DOFs : "<<N_thermalDOF<<"\nHeat function DOFs : "<<N_heatfuncDOF<<"\n"; 
   
   delete [] Rhs_All[3];
   Rhs_All[3] = new double[N_heatfuncDOF];   
   
 
#endif    
    
    delete [] oldsol_T; 
    delete [] heat_defect;
    delete [] Heat_B;
    delete [] Rhs_All[2];
    delete [] Flux_A;
    delete [] Flux_M; 
   
    oldsol_T = new double[N_thermalDOF];
    heat_defect = new double[N_thermalDOF];
    Heat_B = new double[N_thermalDOF]; 
    Rhs_All[2] = new double[N_thermalDOF]; 
    Flux_A = new double[N_thermalDOF];
    Flux_M = new double[N_thermalDOF]; 
    
    memcpy(oldsol_T, Sol_All[2], N_thermalDOF*SizeOfDouble);    
  
    y=0.;
    GetMovingBoundData(NSE_coll, N_MovVert, Bound_Joint, MovBoundVert, Free_Joint,
                       Free_Cells, IsoCellEdgeNos,  SLPX, y);

    GetSolidBoundData(Solid_coll, N_MovVert, Bound_Joint, MovBoundVert, BD_Cells,
                      N_SolidNeibCells, SolidNeibCells, FESpaces_All[3], 0,  N_FluxDof, FluxDof);      

    delete  OutputAll;    
    OutputAll = new TOutput2D(1, 3, 1, 2, Domain);
//     OutputAll->AddFEVectFunct(FEVectFuncts_All[1]); 
//     OutputAll->AddFEVectFunct(FEVectFuncts_All[2]);       
    OutputAll->AddFEFunction(FEFunctions_All[5]);  
#ifdef __HEATLINE__      
   OutputAll->AddFEFunction(FEFunctions_All[8]);    
#endif      
    OutputAll->AddParameter(TDatabase::TimeDB->CURRENTTIME,os.str().c_str());    
    
// ==============================================================================================
// move solid BD begin
// ==============================================================================================
   
    GridPos_S->GridToData();
  
    // reparam the remaining solid vert and move the solid domain
    MovBoundVert[2][0]->GetCoords(tx, Ry);   
    N=N_MovVert[3]-1;     
    hi = (tx-Rx)/(double)N;
    
    for(i=0;i<N;i++)
    {
     MovBoundVert[3][N_MovVert[3]-1 - i]->GetCoords(x, y);
     x += ((double)(N-i))*hi;
     MovBoundVert[3][N_MovVert[3]-1 - i]->SetCoords(x, y);   
    }
    
//     MovBoundVert[3][0]->GetCoords(x, y);    
//     MovBoundVert[2][0]->SetCoords(SLPX, SLPY);
//       
//     y=0.;
//     hi = SLPX-x;
//     hi /= (double)N_MovVert[3];
//     for(i=1;i<N_MovVert[3];i++)
//       MovBoundVert[3][i]->SetCoords(x+(hi*(double)i), y);

    RefGridPos_S->GridToData();
    GridPos_S->DataToGrid();
  
    memset(GridRhs_S, 0, 2*N_GidDofs[2]*SizeOfDouble);
    memcpy(Rhs_All[1], refpos_S, 2*N_GidDofs[2]*SizeOfDouble);
    Daxpy(2*N_GidDofs[2], -1, pos_S, Rhs_All[1]);
    memcpy(GridRhs_S + N_GridActive[2], Rhs_All[1]+N_GridActive[2],
           N_GridBdDofs[2]*SizeOfDouble);   
    memcpy(GridRhs_S + (N_GidDofs[2] + N_GridActive[2]),
           Rhs_All[1]+ (N_GidDofs[2] + N_GridActive[2]), N_GridBdDofs[2]*SizeOfDouble);     
   
    memset(GridSol_S, 0 , 2*N_GidDofs[2]*SizeOfDouble);   
    memcpy(GridSol_S + N_GridActive[2], Rhs_All[1]+N_GridActive[2],
           N_GridBdDofs[2]*SizeOfDouble);    
    memcpy(GridSol_S + (N_GidDofs[2] + N_GridActive[2]),
           Rhs_All[1]+ (N_GidDofs[2] + N_GridActive[2]), N_GridBdDofs[2]*SizeOfDouble);    
    

    SolveGridEquation(Entries_S, GridSol_S, GridRhs_S, GridKCol_S, GridRowPtr_S, N_GidDofs[2]); 
    Daxpy(2*N_GidDofs[2], 1, GridSol_S, pos_S);
  
    GridPos_S->DataToGrid();
 
    // Updating solid boundary 
    for(k=0;k<N_MovVert[3];k++)
     if(k==N_MovVert[3]-1)
      { Bound_Joint[2][k]->UpdateParameters(MovBoundVert[3][k], MovBoundVert[2][0]); }
     else
      { Bound_Joint[2][k]->UpdateParameters(MovBoundVert[3][k], MovBoundVert[3][k+1]); }
      
   }  // if( SLPY <= 1e-8  )
   
//     if( fabs(TDatabase::TimeDB->CURRENTTIME- 2.853553e-03 )<1e-8)
//      exit(0);  
   
  // ===================================================================================== 
  // move solid BD end       
  // end change boundary description    
  // Updating the Quard points on the solid surfaces
  //======================================================================================      
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

//        if(TDatabase::ParamDB->WRITE_VTK)
//        { 
//         if(TDatabase::ParamDB->REACTOR_P22>0)
//          {  
//          ComputeVorticityDivergence(FESpaces_All[0], FEFunctions_All[0], FEFunctions_All[1], Grid_space_NSE, Sol_All[4],Sol_All[5]);
//          }  
//  
//         os.seekp(std::ios::beg);
//         if(img<10) os << "VTK/"<<VtkBaseName<<".0000"<<img<<".vtk" << ends;
//          else if(img<100) os << "VTK/"<<VtkBaseName<<".000"<<img<<".vtk" << ends;
//           else if(img<1000) os << "VTK/"<<VtkBaseName<<".00"<<img<<".vtk" << ends;
//            else if(img<10000) os << "VTK/"<<VtkBaseName<<".0"<<img<<".vtk" << ends;
//             else  os << "VTK/"<<VtkBaseName<<"."<<img<<".vtk" << ends;
//         Output->WriteVtk(os.str().c_str());
// 
//         GetHeatConvectionALEVect(FEVectFuncts_All, GridVect_NSE, GridVect_S);
// 
//         os.seekp(std::ios::beg);
//          if(img<10) os << "VTK/"<<Gnubasename<<".0000"<<img<<".vtk" << ends;
//          else if(img<100) os << "VTK/"<<Gnubasename<<".000"<<img<<".vtk" << ends;
//          else if(img<1000) os << "VTK/"<<Gnubasename<<".00"<<img<<".vtk" << ends;
//          else if(img<10000) os << "VTK/"<<Gnubasename<<".0"<<img<<".vtk" << ends;
//          else  os << "VTK/"<<Gnubasename<<"."<<img<<".vtk" << ends;
//         OutputAll->WriteVtk(os.str().c_str());   
// 
//         img++;
//        }       
// 
//    exit(0);
   
    
    
   if((Angle[0]<10.0) ||(Angle[1]>165.0))
    {
     if(TDatabase::ParamDB->REACTOR_P22>0)
      {  
      ComputeVorticityDivergence(FESpaces_All[0], FEFunctions_All[0], FEFunctions_All[1], Grid_space_NSE, Sol_All[4],Sol_All[5]);
      }        
      
      if(TDatabase::ParamDB->WRITE_VTK)
       { 
        os.seekp(std::ios::beg);
        if(img<10) os << "VTK/"<<VtkBaseName<<"_remesh.0000"<<RemeshImg<<".vtk" << ends;
         else if(img<100) os << "VTK/"<<VtkBaseName<<"_remesh.000"<<RemeshImg<<".vtk" << ends;
          else if(img<1000) os << "VTK/"<<VtkBaseName<<"_remesh.00"<<RemeshImg<<".vtk" << ends;
           else if(img<10000) os << "VTK/"<<VtkBaseName<<"_remesh.0"<<RemeshImg<<".vtk" << ends;
            else  os << "VTK/"<<VtkBaseName<<"."<<RemeshImg<<"_remesh.vtk" << ends;
        Output->WriteVtk(os.str().c_str());

        GetHeatConvectionALEVect(FEVectFuncts_All, GridVect_NSE, GridVect_S);

        os.seekp(std::ios::beg);
         if(img<10) os << "VTK/"<<Gnubasename<<"_remesh.0000"<<RemeshImg<<".vtk" << ends;
         else if(img<100) os << "VTK/"<<Gnubasename<<"_remesh.000"<<RemeshImg<<".vtk" << ends;
         else if(img<1000) os << "VTK/"<<Gnubasename<<"_remesh.00"<<RemeshImg<<".vtk" << ends;
         else if(img<10000) os << "VTK/"<<Gnubasename<<"_remesh.0"<<RemeshImg<<".vtk" << ends;
         else  os << "VTK/"<<Gnubasename<<"."<<RemeshImg<<"_remesh.vtk" << ends;
        OutputAll->WriteVtk(os.str().c_str());   

        RemeshImg++;
       }    
    
//       t1 = GetTime();
    
     RemeshAxial3D_ImpDrop(Domain, FESpaces_All, FEVectFuncts_All, FEFunctions_All,
                            N_MovVert, Bound_Joint, MovBoundVert, Free_Joint, Free_Cells, IsoCellEdgeNos,
                            Sol_All, Rhs_All, SquareStructure_All, Structure_All, SqMat_All, Mat_All,
                            NSE_Cells, Solid_Cells, NSE_coll, Solid_coll,
                            N_GidDofs, N_GridActive, N_GridBdDofs,  Grid_space_NSE,  Grid_space_S,
                            BD_Cells, N_SolidNeibCells, SolidNeibCells,
                            GridSol_NSE, GridRhs_NSE, GridSol_S, GridRhs_S,
                            GridVect_NSE, GridVect_S, GridG1_NSE, GridG2_NSE, GridG1_S, GridG2_S,
                            SquareStructure_NSE, SquareStructure_S, GridSqMat_NSE, GridSqMat_S,
                            N_FluxDof, FluxDof, Triangles, PointNeighb, maxEpV);
   
 
     coll = FESpaces_All[0]->GetCollection();
     N_Cells = coll->GetN_Cells();
     
#ifdef __HEATLINE__   
   N_heatfuncDOF = FESpaces_All[5]->GetN_DegreesOfFreedom();
   N_heatfuncActive = FESpaces_All[5]->GetActiveBound();
   N_heatfuncNonActive = N_heatfuncDOF - N_heatfuncActive; 
    
//    cout<<"Thermal DOFs : "<<N_thermalDOF<<"\nHeat function DOFs : "<<N_heatfuncDOF<<"\n";    
#endif   

     Getcellangle(FESpaces_All[2], Angle);
      OutPut( "MinAngle : "<< Angle[0]<< "  MaxAngle : "<<Angle[1]<< endl);
  
     fes = FESpaces_All[3]->GetAllElements();
 
     N_Remesh ++;
     remeshed=TRUE;

     GlobalNumbers = FESpaces_All[0]->GetGlobalNumbers();
     BeginIndex = FESpaces_All[0]->GetBeginIndex();
     GridGlobalNumbers = FESpaces_All[2]->GetGlobalNumbers();
     GridBeginIndex = FESpaces_All[2]->GetBeginIndex();  
     
     N_Active =  FESpaces_All[0]->GetActiveBound();
     N_U = FESpaces_All[0]->GetN_DegreesOfFreedom();
     N_P = FESpaces_All[1]->GetN_DegreesOfFreedom();
     N_Unknowns = 2*N_U + N_P;
  
     N_G = FESpaces_All[2]->GetN_DegreesOfFreedom();
     N_GActive = FESpaces_All[2]->GetActiveBound();
     N_GBoundaryNodes = N_G - N_GActive;     
    
     N_thermalDOF = FESpaces_All[3]->GetN_DegreesOfFreedom();
     N_thermalActive = FESpaces_All[3]->GetActiveBound();
     N_thermalNonActive = N_thermalDOF - N_thermalActive;   

     //vectors for flux calculations
     delete [] Flux_A;
     delete [] Flux_M; 
    
     Flux_A = new double[N_thermalDOF];
     Flux_M = new double[N_thermalDOF];     
     
     delete [] B; delete [] defect;
     delete [] oldsol_T; delete [] heat_defect; delete [] Heat_B;
     B = new double[N_Unknowns];
     defect = new double[N_Unknowns];
     oldsol_T = new double[N_thermalDOF];  
     heat_defect = new double[N_thermalDOF];
     Heat_B = new double[N_thermalDOF];  

     memcpy(oldsol_T, Sol_All[2], N_thermalDOF*SizeOfDouble);     
     
     GridKCol = SquareStructure_NSE->GetKCol();
     GridRowPtr = SquareStructure_NSE->GetRowPtr();
     
     GridKCol_S = SquareStructure_S->GetKCol();
     GridRowPtr_S = SquareStructure_S->GetRowPtr();
          
     N_G   = N_GidDofs[1];
     N_GActive =  N_GridActive[1];
     N_GBoundaryNodes = N_G - N_GActive;
     
     N_G_S   = N_GidDofs[2];
     N_GActive_S =  N_GridActive[2];
     N_GBoundaryNodes_S = N_G_S - N_GActive_S;
     
     delete [] refpos; delete []  auxpos;  delete [] pos; 
     delete [] tmp_GridSol_NSE; delete [] tmp_Gridd_NSE;
     
     refpos = new double[2*N_G];
     auxpos = new double[2*N_G];
     pos = new double[2*N_G];
     tmp_GridSol_NSE = new double[2*N_G];   
     
     delete RefGridPos; delete AuxGridPos; delete GridPos; 
     RefGridPos = new TFEVectFunct2D(Grid_space_NSE, WString, WString, refpos, N_G, 2);
     AuxGridPos = new TFEVectFunct2D(Grid_space_NSE, WString, WString, auxpos, N_G, 2);
     GridPos = new TFEVectFunct2D(Grid_space_NSE, WString, WString, pos, N_G, 2);

     delete [] refpos_S;  delete []  auxpos_S;  delete [] pos_S; 
     refpos_S = new double[2*N_G_S];
     auxpos_S = new double[2*N_G_S];
     pos_S = new double[2*N_G_S];

     delete RefGridPos_S; delete AuxGridPos_S; delete GridPos_S; 
     RefGridPos_S = new TFEVectFunct2D(Grid_space_S, WString, WString, refpos_S, N_G_S, 2);
     AuxGridPos_S = new TFEVectFunct2D(Grid_space_S, WString, WString, auxpos_S, N_G_S, 2);
     GridPos_S = new TFEVectFunct2D(Grid_space_S, WString, WString, pos_S, N_G_S, 2);

     Max_GridLength = N_GidDofs[2];  
     if(Max_GridLength<N_GidDofs[1])
      Max_GridLength = N_GidDofs[1];   
     tmp_Gridd_NSE = new double[2*Max_GridLength];        
     
     delete Output;
     // prepare output (maxn_fespaces,  maxn_scalar,  maxn_vect, maxn_parameters, domain)
     Output = new TOutput2D(1, 3, 1, 2, Domain);
     Output->AddFEVectFunct(FEVectFuncts_All[0]);
     Output->AddFEFunction(FEFunctions_All[2]);      
     
     if(TDatabase::ParamDB->REACTOR_P22>0)
      {  
       Output->AddFEFunction(FEFunctions_All[6]);      
       Output->AddFEFunction(FEFunctions_All[7]);      
      }  
    
     os.seekp(std::ios::beg);
     Output->AddParameter(TDatabase::TimeDB->CURRENTTIME,os.str().c_str());      
    
     delete  OutputAll;    
     OutputAll = new TOutput2D(1, 3, 1, 2, Domain);
//      OutputAll->AddFEVectFunct(FEVectFuncts_All[1]); 
//      OutputAll->AddFEVectFunct(FEVectFuncts_All[2]);       
     OutputAll->AddFEFunction(FEFunctions_All[5]);  
#ifdef __HEATLINE__      
   OutputAll->AddFEFunction(FEFunctions_All[8]);    
#endif      
     OutputAll->AddParameter(TDatabase::TimeDB->CURRENTTIME,os.str().c_str());   
     
     if(TDatabase::ParamDB->REACTOR_P22>0)
      {  
      ComputeVorticityDivergence(FESpaces_All[0], FEFunctions_All[0], FEFunctions_All[1], Grid_space_NSE, Sol_All[4],Sol_All[5]);
      }  

      if(TDatabase::ParamDB->WRITE_VTK)
       { 
        os.seekp(std::ios::beg);
        if(img<10) os << "VTK/"<<VtkBaseName<<"_remesh.0000"<<RemeshImg<<".vtk" << ends;
         else if(img<100) os << "VTK/"<<VtkBaseName<<"_remesh.000"<<RemeshImg<<".vtk" << ends;
          else if(img<1000) os << "VTK/"<<VtkBaseName<<"_remesh.00"<<RemeshImg<<".vtk" << ends;
           else if(img<10000) os << "VTK/"<<VtkBaseName<<"_remesh.0"<<RemeshImg<<".vtk" << ends;
            else  os << "VTK/"<<VtkBaseName<<"."<<RemeshImg<<"_remesh.vtk" << ends;
        Output->WriteVtk(os.str().c_str());

        GetHeatConvectionALEVect(FEVectFuncts_All, GridVect_NSE, GridVect_S);

        os.seekp(std::ios::beg);
         if(img<10) os << "VTK/"<<Gnubasename<<"_remesh.0000"<<RemeshImg<<".vtk" << ends;
         else if(img<100) os << "VTK/"<<Gnubasename<<"_remesh.000"<<RemeshImg<<".vtk" << ends;
         else if(img<1000) os << "VTK/"<<Gnubasename<<"_remesh.00"<<RemeshImg<<".vtk" << ends;
         else if(img<10000) os << "VTK/"<<Gnubasename<<"_remesh.0"<<RemeshImg<<".vtk" << ends;
         else  os << "VTK/"<<Gnubasename<<"."<<RemeshImg<<"_remesh.vtk" << ends;
        OutputAll->WriteVtk(os.str().c_str());   

        RemeshImg++;
       }    
// exit(0);
 
    } // if((Angle[0]<10.0) ||(Angle[1]>165.0))
  //======================================================================  
  // end Remeshing Begin 
  // Assembeling the grid matrix - Begin
  //======================================================================  

    fesp[0] = Grid_space_NSE;
    SQMATRICES_GRID[0] = GridSqMat_NSE[0];
    SQMATRICES_GRID[0]->Reset();
    SQMATRICES_GRID[1] = GridSqMat_NSE[1];
    SQMATRICES_GRID[1]->Reset();
    SQMATRICES_GRID[2] = GridSqMat_NSE[2];
    SQMATRICES_GRID[2]->Reset();
    SQMATRICES_GRID[3] = GridSqMat_NSE[3];
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
   
    Entries[0] = GridSqMat_NSE[0]->GetEntries();
    Entries[1] = GridSqMat_NSE[1]->GetEntries();
    Entries[2] = GridSqMat_NSE[2]->GetEntries();
    Entries[3] = GridSqMat_NSE[3]->GetEntries();

   // for Dirichlet rows in off-diagonal matrices
   memset(Entries[1] + GridRowPtr[N_GActive], 0, (GridRowPtr[N_G] - GridRowPtr[N_GActive])*SizeOfDouble);
   memset(Entries[2] + GridRowPtr[N_GActive], 0, (GridRowPtr[N_G] - GridRowPtr[N_GActive])*SizeOfDouble);  
   
   // now for solid surface  
    fesp[0] = Grid_space_S;
    SQMATRICES_GRID[0] = GridSqMat_S[0];
    SQMATRICES_GRID[0]->Reset();
    SQMATRICES_GRID[1] = GridSqMat_S[1];
    SQMATRICES_GRID[1]->Reset();
    SQMATRICES_GRID[2] = GridSqMat_S[2];
    SQMATRICES_GRID[2]->Reset();
    SQMATRICES_GRID[3] = GridSqMat_S[3];
    SQMATRICES_GRID[3]->Reset(); 
       
    Assemble2D(1, fesp,
             4, SQMATRICES_GRID,
             0, NULL,
             0, NULL, NULL,
             DiscreteFormGrid,
             GridBoundaryConditions,
             GridBoundValues,
             aux);
    delete aux;   
   
    Entries_S[0] = GridSqMat_S[0]->GetEntries();
    Entries_S[1] = GridSqMat_S[1]->GetEntries();
    Entries_S[2] = GridSqMat_S[2]->GetEntries();
    Entries_S[3] = GridSqMat_S[3]->GetEntries();
   
    // for Dirichlet rows in off-diagonal matrices
    memset(Entries_S[1] + GridRowPtr_S[N_GActive_S], 0, 
           (GridRowPtr_S[N_G_S] - GridRowPtr_S[N_GActive_S])*SizeOfDouble);
    memset(Entries_S[2] + GridRowPtr_S[N_GActive_S], 0, 
           (GridRowPtr_S[N_G_S] - GridRowPtr_S[N_GActive_S])*SizeOfDouble);   
  //======================================================================  
  // end Assembeling the grid matrix
  // nonlinear loop without grid velocity due to remeshing
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
     FreeSurf_Axial3DHeat(SqMat_All[8], SqMat_All[9],  RHSs[0], RHSs[1], BoundCondition,
                          tau, FEFunctions_All[5], T_IntfaceMinMax, FEFunctions_All[0]->GetValues(), Params); 
     
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
      
//         os.seekp(std::ios::beg);
//         if(img<10) os << "VTK/"<<VtkBaseName<<".0000"<<img<<".vtk" << ends;
//          else if(img<100) os << "VTK/"<<VtkBaseName<<".000"<<img<<".vtk" << ends;
//           else if(img<1000) os << "VTK/"<<VtkBaseName<<".00"<<img<<".vtk" << ends;
//            else if(img<10000) os << "VTK/"<<VtkBaseName<<".0"<<img<<".vtk" << ends;
//             else  os << "VTK/"<<VtkBaseName<<"."<<img<<".vtk" << ends;
//         Output->WriteVtk(os.str().c_str());
// 
//         GetHeatConvectionALEVect(FEVectFuncts_All, GridVect_NSE, GridVect_S);
// 
//         os.seekp(std::ios::beg);
//          if(img<10) os << "VTK/"<<Gnubasename<<".0000"<<img<<".vtk" << ends;
//          else if(img<100) os << "VTK/"<<Gnubasename<<".000"<<img<<".vtk" << ends;
//          else if(img<1000) os << "VTK/"<<Gnubasename<<".00"<<img<<".vtk" << ends;
//          else if(img<10000) os << "VTK/"<<Gnubasename<<".0"<<img<<".vtk" << ends;
//          else  os << "VTK/"<<Gnubasename<<"."<<img<<".vtk" << ends;
//         OutputAll->WriteVtk(os.str().c_str());   
// 
//         img++;    
//        exit(0);
       
   } //if (remeshed)  
  
  
  }// for(l=0;l<N_SubSteps;l++) 

 if((m % 1 ) == 0   || m==1 )
  {
   MovBoundVert[0][0]->GetCoords(Lx, Ly);
   MovBoundVert[2][0]->GetCoords(Rx, Ry);

   MovBoundVert[0][N_MovVert[0]-1]->GetCoords(x1, y1);
   MovBoundVert[2][1]->GetCoords(x2, y2);
   MovBoundVert[2][2]->GetCoords(x3, y3);   
   MovBoundVert[2][3]->GetCoords(x4, y4);
   
   
//    cout << x2 <<  "  " << x3<<  "  " << x4 << endl;
//    cout << y2 <<  "  " << y3<<  "  " << y4 << endl;   
//    
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
    OutPut(setw(20)<<"T, wd,AxialZ,Ucl,RAng 1,2,3: " << TDatabase::TimeDB->CURRENTTIME<<"   "<< Rx-Lx
                   <<"   "<< y1<<"   "<< Params[2]<<"   "<<R_Theta[0]<<"   "<<R_Theta[1]<<"   "<<R_Theta[2]<<endl);    

   Get_KE(FEVectFuncts_All[0], Params);  
   OutPut(setw(20)<<"T, Volume, Diff, Rel. Diff : "<< TDatabase::TimeDB->CURRENTTIME<<"   "<< Params[0]
                  <<"   "<< Params[0] - InitVolume<<"   "<< (Params[0] - InitVolume)/InitVolume << endl);
  }

    if( ((m % 200)== 0) || (m==1) )
     {
      os.seekp(std::ios::beg);
      if(N_BData<10) os << "BDData/Boundary.0000"<<N_BData<<".data" << ends;
      else if(N_BData<100) os << "BDData/Boundary.000"<<N_BData<<".data" << ends;
      else if(N_BData<1000) os << "BDData/Boundary.00"<<N_BData<<".data" << ends;
      else if(N_BData<10000) os << "BDData/Boundary.0"<<N_BData<<".data" << ends;
      else  os << "BDData/Boundary."<<N_BData<<".data" << ends;

      std::ofstream dat(os.str().c_str());
      if (!dat)
       {
        cerr << "cannot open file for output" << endl;
        return -1;
       }
      dat << "%% Boundary data created for droplet by MooNMD" << endl;
      dat << "%% Current Time :" << TDatabase::TimeDB->CURRENTTIME << endl;
      for(k=0;k<N_MovVert[0];k++) // no need to set end vertices again
       {
        MovBoundVert[0][k]->GetCoords(x1, y1);

        dat << x1 << " " <<  y1<< endl;
       }
      for(k=0;k<N_MovVert[2];k++) // no need to set end vertices again
       {
        MovBoundVert[2][k]->GetCoords(x1, y1);
        dat << x1 << " " <<  y1<< endl;
       }
        MovBoundVert[1][0]->GetCoords(x1, y1);
        dat << x1 << " " <<  y1<< endl;
        
        MovBoundVert[0][0]->GetCoords(x1, y1);
        dat << x1 << " " <<  y1<< endl;
        
      dat.close();
      cout << endl;
      cout << "Boundary wrote output into file " << endl;
      N_BData++;  
      
   } //if( ((m % 20)== 0) || (m==1) )

 if((m % (int)(TDatabase::TimeDB->STEPS_PER_IMAGE) ) == 0   || m==1 )
  {
     if(TDatabase::ParamDB->WRITE_VTK)
       { 
        if(TDatabase::ParamDB->REACTOR_P22>0)
         {  
         ComputeVorticityDivergence(FESpaces_All[0], FEFunctions_All[0], FEFunctions_All[1], Grid_space_NSE, Sol_All[4],Sol_All[5]);
         }  
 
        os.seekp(std::ios::beg);
        if(img<10) os << "VTK/"<<VtkBaseName<<".0000"<<img<<".vtk" << ends;
         else if(img<100) os << "VTK/"<<VtkBaseName<<".000"<<img<<".vtk" << ends;
          else if(img<1000) os << "VTK/"<<VtkBaseName<<".00"<<img<<".vtk" << ends;
           else if(img<10000) os << "VTK/"<<VtkBaseName<<".0"<<img<<".vtk" << ends;
            else  os << "VTK/"<<VtkBaseName<<"."<<img<<".vtk" << ends;
        Output->WriteVtk(os.str().c_str());

        GetHeatConvectionALEVect(FEVectFuncts_All, GridVect_NSE, GridVect_S);

        os.seekp(std::ios::beg);
         if(img<10) os << "VTK/"<<Gnubasename<<".0000"<<img<<".vtk" << ends;
         else if(img<100) os << "VTK/"<<Gnubasename<<".000"<<img<<".vtk" << ends;
         else if(img<1000) os << "VTK/"<<Gnubasename<<".00"<<img<<".vtk" << ends;
         else if(img<10000) os << "VTK/"<<Gnubasename<<".0"<<img<<".vtk" << ends;
         else  os << "VTK/"<<Gnubasename<<"."<<img<<".vtk" << ends;
        OutputAll->WriteVtk(os.str().c_str());   

        img++;
       }       
   }
#ifdef __ENERGY__
   if( m % 100== 0 )
    {
     OutPut(setw(25)<<TDatabase::TimeDB->CURRENTTIME<<" No. ReParam : " << N_ReParam <<endl);  
     OutPut(setw(25)<<TDatabase::TimeDB->CURRENTTIME<<" No. Remeshed : " << N_Remesh <<endl);  
  
         
     Get_Heat(FEFunctions_All[5], Params);   
     GetInterfaceMinMaxT(FEFunctions_All[5],  T_LSIntfaceMinMax); 
     
     if(Initial_T>1e-3)
     OutPut(setw(25)<<"T, Total_Heat, Heat_Diff_L, Heat_Diff_S, Heat_Diff " << TDatabase::TimeDB->CURRENTTIME<<"   "<<  Params[4] 
                    <<"   "<< (Params[0] - Initial_T_L)/Initial_T  
                    <<"   "<< (Params[1] - Initial_T_S)/Initial_T  
                    <<"   "<< (Params[4] - Initial_T)/Initial_T
                    <<"   "<<  ((Params[0] - Initial_T_L)/Initial_T  +   (Params[1] - Initial_T_S)/Initial_T) -  (Params[4] - Initial_T)/Initial_T   << endl);    
     OutPut(setw(25)<<"T, Heat_LiquidSolid Min, Max,  Heat_FreesurFace Min Max: " << TDatabase::TimeDB->CURRENTTIME<<"   "<< T_LSIntfaceMinMax[0] <<"   "<< T_LSIntfaceMinMax[1]
                    <<"   "<< T_IntfaceMinMax[0] <<"   "<< T_IntfaceMinMax[1]<<endl);  
 
      }  
#endif      
// exit(0);

  } // while(TDatabase::TimeDB->CURRENTTIME< end_time)
   
   

   
   
   
// cout << " test " << gamma << endl;
// exit(0);
}

