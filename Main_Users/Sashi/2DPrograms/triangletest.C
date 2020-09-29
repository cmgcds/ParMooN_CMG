
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
// #include <stdlib.h>
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
// #include <TimeUtilities.h>

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
#include "../TNSE_2D/DropHeat_imping_axial3D.h"

extern "C"
{
  void triangulate(char*, struct triangulateio*,
		   struct triangulateio*, struct triangulateio*);
}

int main(int argc, char* argv[])
{
  TDomain *Domain = new TDomain();
  TDatabase *Database = new TDatabase();
  TFEDatabase2D *FEDatabase = new TFEDatabase2D();
  TCollection *coll, *NSE_coll, *Solid_coll, *mortarcoll = NULL;
  TBaseCell *Me, *cell, **Free_Cells, **NSE_Cells, **Solid_Cells, **BD_Cells, **SolidNeibCells;
  TGridCell **DelCell;
  TFESpace2D *velocity_space, *pressure_space, *streamfunction_space, *convolution_space, *fesps;
  TFESpace2D *pressure_space_output;
  TFESpace2D *Grid_space, *vorticity_space, *thermal_space,*grid_space;
  TFESpace2D *Grid_space_S, *Grid_space_NSE;
  
  TOutput2D *Output, *OutputAll;
  TFEVectFunct2D *RefGridPos, *AuxGridPos, *GridPos;
  TFEVectFunct2D *RefGridPos_S, *AuxGridPos_S, *GridPos_S;  
  TFEFunction2D *fefct[4];
  TFESpace2D *fesp[3], *ferhs_T[3], *ferhs[2];
  TAuxParam2D *aux;
  TDiscreteForm2D *DiscreteFormMatrixT_MRhs;
  TDiscreteForm2D *DiscreteForm;
  TMatrix2D *MATRICES[4];
  TSquareMatrix2D *SQMATRICES[8], *SQMATRICES_GRID[4];
  TSquareMatrix2D *SQMATRICES_HEAT[2];
  TDiscreteForm2D *DiscreteFormGalerkin;
  TDiscreteForm2D *DiscreteFormNLGalerkin;
  TDiscreteForm2D *DiscreteFormGrid, *DiscreteFormHeat, *DiscreteFormHeat_SUPG;
  TSquareMatrix **sqmatrices = (TSquareMatrix **)SQMATRICES;
  TMatrix **matrices = (TMatrix **)MATRICES;
  MatVecProc *MatVect;
  DefectProc *Defect;

  TFESpace2D **FESpaces_All = new TFESpace2D *[5];      
  TFEFunction2D **FEFunctions_All = new TFEFunction2D *[8], *P_output;    
  TFEVectFunct2D **FEVectFuncts_All = new TFEVectFunct2D*[3], *Velo_output;
  TStructure2D **Structure_All = new TStructure2D *[2];
  TSquareStructure2D **SquareStructure_All = new TSquareStructure2D *[3];
  TSquareStructure2D *SquareStructure_S, *SquareStructure_NSE;   
  TSquareMatrix2D **SqMat_All = new TSquareMatrix2D *[12];
  TSquareMatrix2D **GridSqMat_NSE = new TSquareMatrix2D *[4]; 
  TSquareMatrix2D **GridSqMat_S = new TSquareMatrix2D *[4];  
  TMatrix2D **Mat_All = new TMatrix2D *[4];
  TFEVectFunct2D *GridVect_NSE, *GridVect_S;
  TFEFunction2D *GridG1_NSE, *GridG2_NSE, *GridG1_S, *GridG2_S;
  FE2D *fes;

  double total_time,*Coordinates;
  double  t, teta, dt,x,y,gamma, Tgamma, tx,ty,sx,sy, R_Theta[3];;
  double left, right, top, bottom,T_a, T_b;
  double x0, y0,x1,y1,hi, residual, impuls_residual, oldresidual, solver_time;
  double *oldsol_T, HeatFlux, OldHeatFlux, L2HeatFlux;
  double end_time, t1, t2, t4, t3;
  double *B, *defect, *heat_defect, *Heat_B, *RHSs_Heat[1];
  double *RHSs[3], *refpos, *auxpos, *pos, *refpos_S, *auxpos_S, *pos_S;
  double  TX[2], TY[2], solver_time_curr;
  double SLPX, SLPY, *Entries[4], tau, oldtau, limit, *sol_output, Params[10], InitVolume, CurrVolume;  
  double Lx, Ly, Rx, Ry,  x2, y2, x3, y3, x4, y4, fh, fhlimit, fhtot, fhmin, fhmax;
  double *Angle = new double[2], **FreePts = new double *[2];  
  double **Sol_All = new double *[6];
  double **Rhs_All = new double *[3],  *Entries_S[4];
  double *GridSol_NSE,  *GridRhs_NSE, *GridSol_S, *GridRhs_S, *tmp_GridSol_NSE;
  double *Flux_A, *Flux_M, Flux_F, *tmp_Gridd_NSE;
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
  int  N_Hori1, N_Hori2,N_Verti,N_Boundary_Vert,N_thermalDOF,N_thermalActive,N_thermalNonActive;
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
  char NameString[]  = "name";
  char UString[] = "U";
  char PString[] = "P";
  char WString[] = "W";
  char VortString[] = "Vort";
  char DivString[] = "Div";  
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

  BoundCondFunct2D *GridBoundaryConditions[1];
  BoundValueFunct2D *GridBoundValues[1];


  struct triangulateio In, Out;
  
  double Xi[5] = {0., 0.,  8., 8., 0.2};
  double Yi[5] = {0.,-4., -4., 0., 0.};

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

  PRM = TDatabase::ParamDB->BNDFILE;
  GEO = TDatabase::ParamDB->GEOFILE;
//======================================================================
// read boundary parameterization and initialize coarse grid
//======================================================================
   Domain->Init(PRM, GEO);
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
  double refX, hE, phi;
//   double r = TDatabase::ParamDB->P4;
  double r = 2.5e-2;

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

//   N_FreeBound_Vert = int (TDatabase::ParamDB->P6);    //Freesurf except end point
  N_FreeBound_Vert = 50;
  N_AxialBound_Vert = 50;
  N_SlipBound_Vert = 2;
  
  N_Hori1 = 50;      // number of horozontal vertices
  N_Hori2 = 10;
  N_Verti = 5;
 
  N_Boundary_Vert = N_Hori1+N_Hori2+2*N_Verti;

  N_Interf_Vertices = N_FreeBound_Vert+N_SlipBound_Vert+N_AxialBound_Vert+N_Boundary_Vert;
  In.numberofpoints = N_Interf_Vertices-1;
  cout << "In.numberofpoints "<<In.numberofpoints<<endl;
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
 cout<<"331 "<<endl;
  for(i=0;i<N_SlipBound_Vert;i++) // without last point
   {
    y = 0.;
    if(fabs(x)<1.e-12) x = 0.;
    
    In.pointlist[2*In_Index] = x;
    In.pointlist[2*In_Index+1] = y;
    In.pointmarkerlist[In_Index] = CurrComp;
    In.segmentlist[2*In_Index] = In_Index;
    In.segmentlist[2*In_Index+1] = In_Index+1;
    In.segmentmarkerlist[In_Index] = CurrComp;
    In_Index++;  
   }  
   CurrComp++;
    cout<<endl;  

  for(i=0;i<N_FreeBound_Vert;i++)
    {
        
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
      In.pointmarkerlist[In_Index] = CurrComp;
      In.segmentlist[2*In_Index] = In_Index;
      In.segmentlist[2*In_Index+1] = In_Index+1;
      In.segmentmarkerlist[In_Index] = CurrComp;
      In_Index++;
    }
     
   OutPut("hE : " << hE<<endl);     
cout<<"378 "<<endl;
   CurrComp++;
   y = (b+r*sin(t)) - In.pointlist[1];
   dt= -y / (double)(N_AxialBound_Vert);
   x = 0.;
   t = y;   
   N = N_AxialBound_Vert;

  for(i=0;i<N;i++)
   {
      if (fabs(y)<1e-8) y = 0.;
      In.pointlist[2*In_Index] = x;
      In.pointlist[2*In_Index+1] = y;
      In.pointmarkerlist[In_Index] = CurrComp;
      In.segmentlist[2*In_Index] = In_Index;
      In.segmentlist[2*In_Index+1] = In_Index+1;
      In.segmentmarkerlist[In_Index] = CurrComp;
      In_Index++;
      y = t + double(i+1)*dt;
    }
    cout<<"398 "<<endl;
  In.segmentlist[2*(In_Index-1)+1] = 0; 
   CurrComp++;
   dt= (Yi[1] - Yi[0])  / N_Verti;
   x = Xi[0];
   y = Yi[0];
   t = y;

  In.segmentlist[2*In_Index] = 0;
  In.segmentlist[2*In_Index+1] = In_Index;
  In.segmentmarkerlist[In_Index] = CurrComp;
  In_Index++;  
   
   for(i=1;i<N_Verti;i++)
    {
      y = t + ((double)i)*dt;
      In.pointlist[2*(In_Index-1)] = x;
      In.pointlist[2*(In_Index-1)+1] = y;
      In.pointmarkerlist[(In_Index-1)] = CurrComp;      
      In.segmentlist[2*In_Index] = In_Index-1;
      In.segmentlist[2*In_Index+1] = In_Index;
      In.segmentmarkerlist[In_Index] = CurrComp;
      In_Index++;
    }
 cout<<"422 "<<endl;
   CurrComp++;
   dt= (Xi[2] - Xi[1])  / N_Hori2;
   x = Xi[1];
   y = Yi[1];
   t = x;

   for(i=0;i<N_Hori2;i++)
    {
      In.pointlist[2*(In_Index-1)] = x;
      In.pointlist[2*(In_Index-1)+1] = y;
      In.pointmarkerlist[(In_Index-1)] = CurrComp;
      In.segmentlist[2*In_Index] = In_Index-1;
      In.segmentlist[2*In_Index+1] = In_Index;
      In.segmentmarkerlist[In_Index] = CurrComp;
      In_Index++;
      x = t + ((double)(i+1))*dt;
    }
cout<<"440 "<<endl;
   CurrComp++;
   dt= (Yi[3] - Yi[2])  / N_Verti;
   x = Xi[2];
   y = Yi[2];
   t = y;

   for(i=0;i<N_Verti;i++)
    {
      In.pointlist[2*(In_Index-1)] = x;
      In.pointlist[2*(In_Index-1)+1] = y;
      In.pointmarkerlist[(In_Index-1)] = CurrComp;
      In.segmentlist[2*In_Index] = In_Index-1;
      In.segmentlist[2*In_Index+1] = In_Index;
      In.segmentmarkerlist[In_Index] = CurrComp;
      In_Index++;
      y = t + ((double)(i+1))*dt;
    }

   refX = Xi[4] + 0.2;   

   N_refX = (int)((refX - Xi[4])/hE);  
   cout<<"462 "<<endl;
   if( (N_refX+10)>N_Hori1)
   {
    cout << "Increase  N_Hori1 points, N_refX:  "<<  N_refX <<endl;
    exit(0);
   } 

   CurrComp++;
   dt= (refX - Xi[3])  / (N_Hori1-N_refX);
   x = Xi[3];
   y = Yi[3];
   t = x;   

   for(i=0;i<(N_Hori1-N_refX);i++)
    {
      In.pointlist[2*(In_Index-1)] = x;
      In.pointlist[2*(In_Index-1)+1] = y;
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
cout<<"491 "<<endl;
   for(i=0;i<N_refX;i++)
    {
      In.pointlist[2*(In_Index-1)] = x;
      In.pointlist[2*(In_Index-1)+1] = y;
      In.pointmarkerlist[(In_Index-1)] = CurrComp;
      In.segmentlist[2*In_Index] = (In_Index-1);
      In.segmentlist[2*In_Index+1] = In_Index;
      In.segmentmarkerlist[In_Index] = CurrComp;
      In_Index++;
      x = t + ((double)(i+1))*dt;
    }

  In.segmentlist[2*(In_Index-1)+1] = Indextemp;
cout<<"505 "<<endl;
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
cout<<"536 "<<endl;
Domain->GetTreeInfo(CellTree ,N_RootCells);
  coll = Domain->GetCollection(It_Finest, 0);
  N_Cells = coll->GetN_Cells();

  VertexDel = new TVertex*[3*N_RootCells];

  CurrVertex = 0;
  for(i=0;i<N_Cells;i++)
    {
      cell = coll->GetCell(i);
      N_Joints = cell->GetN_Joints();
      N_Vertices = cell->GetN_Vertices();
      cout<<"549 "<<endl;
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
           } 
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
        } 
    }
    cout<<"584 "<<endl;
   for(i=0;i<CurrVertex;i++)
   delete VertexDel[i];

//    delete [] VertexDel;
   OutPut(CurrVertex<<" vertices were deleted"<<endl);

// remove all existing cells and joints
    for(i=0;i<N_RootCells;i++)
    delete (TGridCell*)CellTree[i];
    OutPut(N_RootCells<<" cells were deleted"<<endl);
    delete CellTree;
    delete coll;
cout<<"597 "<<endl;
    N_RootCells = Out.numberoftriangles;
  // allocate auxillary fields
    cout<<"600 "<<endl;
  Coordinates = Out.pointlist;
  Triangles = Out.trianglelist;
  cout<<"603 "<<endl;
  cout<<"numberofpoints "<<Out.numberofpoints<<endl;
  PartMarker = new int[Out.numberofpoints];
//   PartMarker = new int[2];
  cout<<"605 "<<endl;
// generate new vertices
  N_G = Out.numberofpoints;
  cout<<"606 "<<endl;
  NewVertices = new TVertex*[N_G];
//   NewVertices = new TVertex*[2];
cout<<"608 "<<endl;
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
cout<<"616 "<<endl;
  UpdateSlipBound->SetParams(Xi[0], Yi[0], Xi[4]-Xi[0],Yi[4]-Yi[0]);
  UpdateSlipBoundSolid->SetParams(Xi[3], Yi[3], Xi[4]-Xi[3],Yi[4]-Yi[3]);
  SolidBound1->SetParams(Xi[0], Yi[0], Xi[1]-Xi[0],Yi[1]-Yi[0]);
  SolidBound2->SetParams(Xi[1], Yi[1], Xi[2]-Xi[1],Yi[2]-Yi[1]); 
  SolidBound3->SetParams(Xi[2], Yi[2], Xi[3]-Xi[2],Yi[3]-Yi[2]);
cout<<"622 "<<endl;
//  OutPut("left: "<<left<<" right: "<<right<<" top: "<<top<<" bottom: "<<bottom<<endl);
  Domain->SetBoundBox(right-left,top-bottom);
  Domain->SetBoundBoxstart(left,bottom);
cout<<"626 "<<endl;
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


      if(Domain->GetBdPart(0)->GetBdComp(CurrComp)->GetTofXY(NewVertices[a]->GetX(), NewVertices[a]->GetY(), T_a) ||
          Domain->GetBdPart(0)->GetBdComp(CurrComp)->GetTofXY(NewVertices[b]->GetX(), NewVertices[b]->GetY(), T_b))
       {
          cerr<<"Error: could not set parameter values"<<endl;
          OutPut(NewVertices[a]<<endl);
          OutPut(NewVertices[b]<<endl);
          cout << " CurrComp " << CurrComp <<endl;
       }

      if (CurrNeib == 2)    // 2 cells contain the current edge
        if(Domain->GetBdPart(0)->GetBdComp(CurrComp)->IsFreeBoundary())
          Joint = new TIsoInterfaceJoint(Domain->GetBdPart(0)->GetBdComp(CurrComp), T_a, T_b, CellTree[Neib[0]], CellTree[Neib[1]]);
        else
          Joint = new TInterfaceJoint(Domain->GetBdPart(0)->GetBdComp(CurrComp), T_a, T_b, CellTree[Neib[0]], CellTree[Neib[1]]);
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

  if (Joint->GetType() == InterfaceJoint || Joint->GetType() == IsoInterfaceJoint)
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
  os.seekp(std::ios::beg);
  os << "Domain_0" << ".ps" << ends;
  Domain->PS(os.str().c_str(),It_Finest,0);
}