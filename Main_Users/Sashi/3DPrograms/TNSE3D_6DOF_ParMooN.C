// =======================================================================
//
// Purpose:     main program for solving a flow over floating body problem in ParMooN
//
// Author:      Sashikumaar Ganesan
//
// History:     Implementation started on 17.12.2015
//        :     6 DOF solver added on 13.05.2017   
// =======================================================================
#include <Domain.h>
#include <Database.h>
#include <FEDatabase3D.h>
#include <LinAlg.h>
#include <FESpace3D.h>
#include <SystemTNSE3D_ALE.h>
#include <System_6DOF.h>
#include <SquareStructure3D.h>
#include <Structure3D.h>
#include <Output3D.h>
#include <LinAlg.h>
#include <MainUtilities.h>
#include <TimeDiscRout.h>

#include <tetgen.h>

#include <string.h>
#include <sstream>
#include <MooNMD_Io.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>

#ifdef _SMPI
#include "mpi.h"
#include "MumpsSolver.h"
#endif

#ifdef _MPI
#include "mpi.h"
#include <MeshPartition.h>
#endif
#include <fstream> 

double bound = 0;
double timeC = 0;

#define AMG 0
#define GMG 1
#define DIRECT 2
// =======================================================================
// include current example
// =======================================================================
//  #include "../Examples/TNSE_3D/DrivenCavity3D.h"
//  #include "../Examples/TNSE_3D/Aerofoil_SD7003_3D.h"
// #include "../Examples/TNSE_3D/Aerofoil_3D_ALE_slip.h"
//   #include "../Examples/TNSE_3D/Hole3D.h"
//  #include "../Examples/TNSE_3D/AnsatzLinConst.h"
//  #include "../Examples/TNSE_3D/Bsp3.h"
// #include "../Examples/TNSE_3D/Channel3D.h"
// #include "../Examples/TNSE_3D/Channel3D_slip.h"
// #include "../Examples/TNSE_3D/ChannelObstacle3D_ALE.h"
//#include "../Examples/TNSE_3D/ChannelObstacle3D_ALE_slip.h"
#include "../Examples/TNSE_3D/Channel_6DOF_ALE.h"

// =======================================================================
// main program
// =======================================================================
int main(int argc, char* argv[])
{
  // ======================================================================
  //  declaration of variables
  // ======================================================================
  int i, j, l, m, N_Cells, N_Cells_Solid, N_U, N_P, N_TotalDOF, img=1, pressure_space_code;
  int Max_It, NSEType, velocity_space_code;
  int LEVELS, mg_level, mg_type;
  int N_SubSteps, N_L;
  int  N_MovVert, *Movfaces, N_Movfaces, N_G, N_RB_G, checker = 1;
    
  double **MeshSol, **Mesh_RB;
  double *sol, *rhs, *oldrhs, t1, t2, errors[4], residual, impuls_residual;
  double **Sol_array, **Rhs_array;
  double limit, AllErrors[7];
  double tau, oldtau, end_time, currtime;
  double start_time, stop_time, start_assembling_solving=0, end_assembling_solving=0,
         total_assembling_solving=0, start_int=0, end_int=0, total_int=0;
	 
  TVertex **MovBoundVert;
  TBaseCell ***Free_Cells, *Cell;
  TBaseCell ** MovCells;
  TIsoBoundEdge ***Free_Joint;
  TFEVectFunct3D **MeshVelocity, **MeshVelocity_RB;  
  TDomain *Domain;
  TDatabase *Database = new TDatabase();
  
  int profiling;
  
#if defined(_MPI) || defined(_SMPI)
  int rank, size, out_rank;
  int MaxCpV, MaxSubDomainPerDof;
  int Refine;
  int metisType[2] = {0,0};  
  
  MPI_Comm Comm = MPI_COMM_WORLD;
  MPI_Init(&argc, &argv);

  double time1, time2;

  MPI_Comm_rank(Comm, &rank);
  MPI_Comm_size(Comm, &size);
 
  TDatabase::ParamDB->Comm = Comm;
#endif 
  
  TFEDatabase3D *FEDatabase = new TFEDatabase3D(); 
  TCollection *coll, *coll_Solid, *mortarcoll = NULL;
  TFESpace3D *velocity_space, *pressure_space, **Velocity_FeSpace, **Pressure_FeSpace, *fesp[1];
  TFESpace3D **GridFESpace, **RigidbodyFESpace, **testfespace;
  TFEVectFunct3D **Velocity, *u;
  TFEFunction3D *p, *u1, *u2, *u3, **Pressure, *fefct[2];
  TOutput3D *Output;
  TSystemTNSE3D_ALE *SystemMatrix_ALE;
  TSystem_6DOF *SystemMatrix_6DOF;
  TAuxParam3D *aux;
  MultiIndex3D AllDerivatives[4] = { D000, D100, D010, D001 };
  TBaseCell **Cells_solid;
  
  TFESpace3D *projection_space;
  TFESpace3D **Projection_FeSpace;
  
  const char vtkdir[] = "VTK"; 
  char *PsBaseName, *VtkBaseName, *GEO, *PRM, *SMESH;
  char Name[] = "name";
  char UString[] = "u";
  char PString[] = "p";
  char NameString[] = "VMS";
  char WString[] = "w";
  char SubID[] = "";
  
  std::ostringstream os;
  
#ifdef _SMPI
  TMumpsSolver *solmumps;
   
  if(rank==0)  
#endif 
{
  os << " ";    
  mkdir(vtkdir, 0777);
    
//======================================================================
// read parameter file
//====================================================================== 
  /** set variables' value in TDatabase using argv[1] (*.dat file), and generate the MESH based */
  Domain = new TDomain(argv[1]);  
  profiling = TDatabase::ParamDB->timeprofiling;
  cout <<" USE PRM " << TDatabase::ParamDB->USE_PRM << endl;
  if(profiling)
  {
#ifdef _MPI
    start_time = MPI_Wtime();
#else
    start_time = GetTime();
#endif
  }

  OpenFiles();
  OutFile.setf(std::ios::scientific);
   
  Database->CheckParameterConsistencyNSE();
  
#ifdef _MPI
  out_rank=TDatabase::ParamDB->Par_P0;
  //out_rank = 0;
  if(rank == out_rank)
#endif
   {
    Database->WriteParamDB(argv[0]);
//     Database->WriteTimeDB();
    ExampleFile();
   }
  
  /** needed in the new implementation */
  if(TDatabase::ParamDB->FLOW_PROBLEM_TYPE == STOKES)
   {  TDatabase::ParamDB->SC_NONLIN_MAXIT_SADDLE = 1; }
 
 
//======================================================================
// initialize coarse grid and generate the mesh
//======================================================================

  if(TDatabase::ParamDB->MESH_TYPE==0)
    {Domain->Init(TDatabase::ParamDB->BNDFILE, TDatabase::ParamDB->GEOFILE); } // ParMooN  build-in Geo mesh
  else if(TDatabase::ParamDB->MESH_TYPE==1)  
     {Domain->GmshGen(TDatabase::ParamDB->GEOFILE); }//gmsh mesh
  else if(TDatabase::ParamDB->MESH_TYPE==2)   
    {Domain->TetrameshGen(TDatabase::ParamDB->GEOFILE); } //tetgen mesh
    else
     {  
      OutPut("Mesh Type not known, set MESH_TYPE correctly!!!" << endl);
      exit(0);
     } 
   
  LEVELS = TDatabase::ParamDB->LEVELS;
  if(TDatabase::ParamDB->SOLVER_TYPE==DIRECT)
  {
    TDatabase::ParamDB->UNIFORM_STEPS += (LEVELS-1);
    LEVELS = 1;
  }
  
  // refine grid up to the coarsest level
  for(i=0;i<TDatabase::ParamDB->UNIFORM_STEPS;i++)
    Domain->RegRefineAll();  
  
// #ifdef __Cylinder__
//    TetrameshGen(Domain);
// #endif

//    exit(0);

#ifdef _MPI
  Domain->GenerateEdgeInfo();
  
  if(profiling)  t1 = MPI_Wtime();
  
  if(rank == 0)
       {
	  printf("\n----------------------------------------------------------------------------------------\n");
	  printf("metis type set to %d\n",TDatabase::ParamDB->Par_P2);
	  printf("----------------------------------------------------------------------------------------\n\n");
       }
  //this loop checks if number of cells are sufficient in the coarsest level, such that each 
  //rank get some own cells to work on
  //it does so by changing the metis type first, if not possible then refine and start again
  do
  {
    metisType[TDatabase::ParamDB->Par_P2] = 1;
    Refine = Partition_Mesh3D(Comm, Domain, MaxCpV);	//MaxCpV=maximum cell per vertex
    
    if(metisType[0]*metisType[1] == 1 && Refine)
    {
      metisType[0] = 0;      metisType[1] = 0;
      TDatabase::ParamDB->Par_P2 = 0;
      if(rank == 0)
       {
	  printf("\n----------------------------------------------------------------------------------------\n");
	  printf("Warning :: both metisType used. Now refining the mesh by one step \n");
	  printf("metis type set to 0\n");
	  printf("----------------------------------------------------------------------------------------\n\n");
       }
      Domain->RegRefineAll();
      Domain->GenerateEdgeInfo();
      TDatabase::ParamDB->UNIFORM_STEPS +=1;
    }
  }while(Refine);
  
  if(profiling)  t2 = MPI_Wtime(); 
  
  if(profiling){
    time2 = t2-t1;
    MPI_Reduce(&time2, &time1, 1, MPI_DOUBLE, MPI_MAX, out_rank, Comm);
    if(rank == out_rank)
      printf("Time taken for Domain Decomposition is %e\n", time1);
  }

  Domain->GenerateEdgeInfo();
  MaxSubDomainPerDof = MIN(MaxCpV, size);
  TDatabase::ParamDB->WRITE_PS = 0;
#endif 
  
  if(TDatabase::ParamDB->WRITE_PS)
   {
    // write grid into an Postscript file
    os.seekp(std::ios::beg);
    os << "Domain" << ".ps" << ends;
    Domain->PS(os.str().c_str(),It_Finest,0);
   }
// MPI_Finalize();
// exit(0);

//=========================================================================
// set data for multigrid
//=========================================================================  

  // set type of multilevel
  mg_type = TDatabase::ParamDB->SC_MG_TYPE_SADDLE;
 
  if(TDatabase::ParamDB->SOLVER_TYPE==AMG_SOLVE|| TDatabase::ParamDB->SOLVER_TYPE==DIRECT)
   { 
    mg_type=0; 
    TDatabase::ParamDB->SC_MG_TYPE_SADDLE = mg_type;
   }
  
  if(mg_type==1)
   { mg_level =  LEVELS + 1; }
  else
   { mg_level = LEVELS; }
   
  if(TDatabase::ParamDB->SOLVER_TYPE==GMG)
#ifdef _MPI  
    if(rank == out_rank)
#endif   
  {   
    OutPut("=======================================================" << endl);
    OutPut("======           GEOMETRY  LEVEL ");
    OutPut(LEVELS << "              ======" << endl);
    OutPut("======           MULTIGRID LEVEL ");
    OutPut(mg_level << "              ======" << endl);
    OutPut("=======================================================" << endl);   
   }   

  Velocity_FeSpace = new TFESpace3D*[mg_level];  
  Pressure_FeSpace = new TFESpace3D*[mg_level];  
  Projection_FeSpace = new TFESpace3D*[mg_level];  
  
  Velocity = new TFEVectFunct3D*[mg_level] ;
  Pressure = new TFEFunction3D*[mg_level];
  
  Sol_array = new double*[mg_level];
  Rhs_array = new double*[mg_level];
  MeshSol = new double*[mg_level];
  Mesh_RB = new double*[mg_level];
  GridFESpace = new TFESpace3D*[mg_level];
  MeshVelocity = new TFEVectFunct3D*[mg_level];
  MeshVelocity_RB = new TFEVectFunct3D*[mg_level];
  RigidbodyFESpace  = new TFESpace3D*[mg_level];
  testfespace = new TFESpace3D*[mg_level];
//=========================================================================
// loop over all levels (not a multigrid level but for convergence study)  
//=========================================================================
  for(i=0;i<LEVELS;i++)
   {   
    if(i)
     { Domain->RegRefineAll(); }

#ifdef _MPI
     if(rank == out_rank)
       printf("Level :: %d\n\n",i);
     if(i)
     {
       Domain->GenerateEdgeInfo();
       Domain_Crop(Comm, Domain);       // removing unwanted cells in the hallo after refinement 
     }
#endif

   coll=Domain->GetCollection(It_Finest, 0);

   if(TDatabase::ParamDB->P1==1) // two phase,
    {
     N_Cells = coll->GetN_Cells();
     N_Cells_Solid=0;
   
     for( j=0;j<N_Cells;j++)
      {
       Cell = coll->GetCell(j);
       // collect solid coll with given regin ID
       if( Cell->GetRegionID()== (int)TDatabase::ParamDB->P2) 
         N_Cells_Solid++;
      } // for j

      Cells_solid = new TBaseCell* [N_Cells_Solid];
      N_Cells_Solid = 0;
      for( j=0;j<N_Cells;j++)
      {
       Cell = coll->GetCell(j);
       // collect solid coll with given regin ID
       if( Cell->GetRegionID()== (int)TDatabase::ParamDB->P2) 
       {
        Cells_solid[N_Cells_Solid] = Cell;
        N_Cells_Solid++;
       }
      } // for j    
     
      coll_Solid = new TCollection (N_Cells_Solid, Cells_solid);
      cout << "N_Cells_Solid : " << N_Cells_Solid<<endl;
//      exit(0);
      
     }  // if(TDatabase::ParamDB->P1==1) //
   else
    {
     coll_Solid = coll;
     cout<<"  coll_Solid = coll " <<endl;
    }

//=========================================================================
// construct all finite element spaces
//=========================================================================  
     if(mg_type==1) // lower order FE on coarse grids 
      {
       Velocity_FeSpace[i] = new TFESpace3D(coll, Name, UString, BoundCondition, Non_USpace,1);
       Pressure_FeSpace[i] = new TFESpace3D(coll, Name, PString, BoundCondition, DiscP_PSpace,0);    
       
       if(i==LEVELS-1) // higher order on fine level
        {
         GetVelocityAndPressureSpace3D(coll, BoundCondition, velocity_space,
                                       pressure_space, &pressure_space_code,
                                       TDatabase::ParamDB->VELOCITY_SPACE,
                                       TDatabase::ParamDB->PRESSURE_SPACE); 
         Velocity_FeSpace[i+1] = velocity_space;
         Pressure_FeSpace[i+1] = pressure_space;

         // defaulty inf-sup pressure space will be selected based on the velocity space, so update it in database
         TDatabase::ParamDB->INTERNAL_PRESSURE_SPACE = pressure_space_code;
         velocity_space_code = TDatabase::ParamDB->VELOCITY_SPACE; 
 
#ifdef _MPI
         Velocity_FeSpace[LEVELS]->SetMaxSubDomainPerDof(MaxSubDomainPerDof);
         Pressure_FeSpace[LEVELS]->SetMaxSubDomainPerDof(MaxSubDomainPerDof);
#endif
        }
       
      }
     else
      {


       GetVelocityAndPressureSpace3D(coll, BoundCondition, velocity_space, pressure_space, &pressure_space_code,
                                     TDatabase::ParamDB->VELOCITY_SPACE, TDatabase::ParamDB->PRESSURE_SPACE); 
       Velocity_FeSpace[i] = velocity_space;
       Pressure_FeSpace[i] = pressure_space;       
       


       TDatabase::ParamDB->INTERNAL_PRESSURE_SPACE = pressure_space_code;
       velocity_space_code = TDatabase::ParamDB->VELOCITY_SPACE;
      }
       
  if (TDatabase::ParamDB->DISCTYPE == VMS_PROJECTION)
     {
      if (TDatabase::ParamDB->VMS_LARGE_VELOCITY_SPACE==0)
        projection_space = new TFESpace3D(coll,NameString, UString, BoundCondition, DiscP_PSpace, 0);
      else
        projection_space = new TFESpace3D(coll,NameString, UString, BoundCondition, DiscP_PSpace, 1);
          
      Projection_FeSpace[i] = projection_space;
      N_L = Projection_FeSpace[i]->GetN_DegreesOfFreedom();
      OutPut("Dof Projection : " << setw(10) << N_L << endl);     
     }
 
#ifdef _MPI
     Velocity_FeSpace[i]->SetMaxSubDomainPerDof(MaxSubDomainPerDof);
     Pressure_FeSpace[i]->SetMaxSubDomainPerDof(MaxSubDomainPerDof);
#endif

     
//======================================================================
// construct all finite element functions
//======================================================================   
     N_U = Velocity_FeSpace[i]->GetN_DegreesOfFreedom();
     N_P = Pressure_FeSpace[i]->GetN_DegreesOfFreedom();
     N_TotalDOF = 3*N_U + N_P;    

     sol = new double[N_TotalDOF];
     memset(sol, 0, N_TotalDOF*SizeOfDouble);     
     Sol_array[i] = sol;
     
     rhs = new double[N_TotalDOF];
     memset(rhs, 0, N_TotalDOF*SizeOfDouble);     
     Rhs_array[i] = rhs;

     u = new TFEVectFunct3D(Velocity_FeSpace[i], UString,  UString,  sol, N_U, 3);
     Velocity[i] = u;
     p = new TFEFunction3D(Pressure_FeSpace[i], PString,  PString,  sol+3*N_U, N_P);    
     Pressure[i] = p;
     
//======================================================================      
// mesh velocity space 
//======================================================================   
   GridFESpace[i] = new TFESpace3D(coll, NameString, WString, GridBoundCondition, 1);
   N_G = GridFESpace[i] -> GetN_DegreesOfFreedom();
   OutPut("N_G    : "<< setw(10) << N_G  << endl);



   MeshSol[i] = new double[3*N_G];
   memset(MeshSol[i], 0, 3*N_G*SizeOfDouble);      
   MeshVelocity[i] = new TFEVectFunct3D(GridFESpace[i], WString, WString, MeshSol[i], N_G, 3);     
    
   // information for moving mesh based on boundary displacement 
   if(i==LEVELS-1)
    {
     GetMovingBoundData(coll,  N_MovVert, MovBoundVert, N_Movfaces, Movfaces, MovCells); 
     OutPut("N_MovVert  : " << N_MovVert <<endl); 
    }

#ifndef _NSEONLY_
//=======================================================================================
// Rigidbody space 
//======================================================================   
   RigidbodyFESpace[i] = new TFESpace3D(coll_Solid, NameString, WString, GridBoundCondition, 1);
   N_RB_G = RigidbodyFESpace[i] -> GetN_DegreesOfFreedom();
   OutPut("N_RB_G    : "<< setw(10) << N_RB_G  << endl);
   
   Mesh_RB[i] = new double[3*N_RB_G];
   memset(Mesh_RB[i], 0, 3*N_RB_G*SizeOfDouble);      
   MeshVelocity_RB[i] = new TFEVectFunct3D(RigidbodyFESpace[i], WString, WString, Mesh_RB[i], N_RB_G, 3);     
//=======================================================================================
#endif

     if(i==LEVELS-1 && mg_type==1)  
      {
       N_U = Velocity_FeSpace[i+1]->GetN_DegreesOfFreedom();
       N_P = Pressure_FeSpace[i+1]->GetN_DegreesOfFreedom();    
       N_TotalDOF = 3*N_U + N_P;    

       sol = new double[N_TotalDOF];
       memset(sol, 0, N_TotalDOF*SizeOfDouble);     
       Sol_array[i+1] = sol;
     
       rhs = new double[N_TotalDOF];
       memset(rhs, 0, N_TotalDOF*SizeOfDouble);     
       Rhs_array[i+1] = rhs;

       u = new TFEVectFunct3D(Velocity_FeSpace[i+1], UString,  UString,  sol, N_U, 3);
       Velocity[i+1] = u;
       p = new TFEFunction3D(Pressure_FeSpace[i+1], PString,  PString,  sol+3*N_U, N_P);    
       Pressure[i+1] = p;    
      }// if(i==LEVELS-1 && mg_type==1)
#ifdef _MPI
     N_Cells = coll->GetN_Cells();
     printf("rank=%d\t N_Cells   : %d\t Dof all   :%d\t Dof Velocity :%d\t Dof Pressure: %d\n",rank,N_Cells,N_TotalDOF,3*N_U,N_P);
#endif
    } //  for(i=0;i<LEVELS;i++)
   
    u1 = Velocity[mg_level-1]->GetComponent(0);
    u2 = Velocity[mg_level-1]->GetComponent(1);
    u3 = Velocity[mg_level-1]->GetComponent(2);  

    oldrhs = new double[N_TotalDOF];
      
#ifndef _MPI       
    N_Cells = coll->GetN_Cells();
    OutPut("N_Cells      : "<< setw(10) << N_Cells <<endl);
    OutPut("Dof velocity : "<< setw(10) << 3*N_U << endl);
    OutPut("Dof pressure : "<< setw(10) << N_P << endl);
    OutPut("Dof all      : "<< setw(10) << N_TotalDOF  << endl);  
#endif 

    //======================================================================
      // prepare for output
      //======================================================================
         VtkBaseName = TDatabase::ParamDB->VTKBASENAME;
         Output = new TOutput3D(2, 2, 1, 1, Domain);

      #ifndef _6DOFONLY_
         Output->AddFEVectFunct(u);
         Output->AddFEVectFunct(MeshVelocity[LEVELS-1]);
         Output->AddFEFunction(p);
      #else
         Output->AddFEFunction(p);
      #endif

      #ifdef _MPI
           if(TDatabase::ParamDB->WRITE_VTK)
            Output->Write_ParVTK(MPI_COMM_WORLD, img, SubID);
              img++;
      #else
         // Added by thivin - to remove previous VTK files
        cout << "===================================================================" << endl; 
        system("exec rm /home/thivin/PARMOON_CODES/parmoon_hg_repo/ParMooN_Output/6DOF/VTK/*");
        cout << "Thivin -  Removed the additional VTK Files " << endl;
        cout << "===================================================================" << endl; 

         if(TDatabase::ParamDB->WRITE_VTK)
           {
            os.seekp(std::ios::beg);
             if(img<10) os <<  "VTK/"<<VtkBaseName<<".0000"<<img<<".vtk" << ends;
               else if(img<100) os <<  "VTK/"<<VtkBaseName<<".000"<<img<<".vtk" << ends;
                else if(img<1000) os <<  "VTK/"<<VtkBaseName<<".00"<<img<<".vtk" << ends;
                 else if(img<10000) os <<  "VTK/"<<VtkBaseName<<".0"<<img<<".vtk" << ends;
                  else  os <<  "VTK/"<<VtkBaseName<<"."<<img<<".vtk" << ends;
            Output->WriteVtk(os.str().c_str());
            img++;
           }
      #endif

  //     exit(0);

//======================================================================
// SystemMatrix_ALE construction and solution
//======================================================================  
    NSEType = TDatabase::ParamDB->NSTYPE;
    
if(profiling)
   {
#ifdef _MPI
      start_int = MPI_Wtime();
#else
      start_int = GetTime();
#endif
    }
                     
   // get a  TNSE3D_ALE system
   SystemMatrix_ALE = new TSystemTNSE3D_ALE(mg_level, Velocity_FeSpace, Pressure_FeSpace, Velocity, 
                                            Pressure, Sol_array, Rhs_array, TDatabase::ParamDB->DISCTYPE,
                                            NSEType, TDatabase::ParamDB->SOLVER_TYPE, Projection_FeSpace, 
                                            GridFESpace, MeshVelocity, FALSE);  

    // initilize the system matrix with the functions defined in Example file
    SystemMatrix_ALE->Init(LinCoeffs, BoundCondition, U1BoundValue, U2BoundValue, U3BoundValue,
                                              GridBoundCondition, GridBoundValue, GridCoeffs);

    // add the info of moving boundaries
    //     SystemMatrix_ALE->TSystemTNSE3D_ALE::AddBoundModifyFunction(ModifyBdCoords, N_MovVert, MovBoundVert, N_Movfaces, Movfaces, MovCells);

    // thivin - commented to revove fluid vel
    // move the boundary with the fluid velocity 
   SystemMatrix_ALE->TSystemTNSE3D_ALE::AddMoveBoundFunction(MoveBoundWithVelo, N_MovVert, MovBoundVert, N_Movfaces, Movfaces, MovCells);
    
    
#ifndef _NSEONLY_
    // initialize the 6dof
    SystemMatrix_6DOF = new TSystem_6DOF(RigidbodyFESpace[LEVELS-1], MeshVelocity_RB[LEVELS-1]);
    SystemMatrix_6DOF->Init6DOF();
#endif

#ifdef _MPI
    if(rank==0)
#endif
    {
#ifndef _6DOFONLY_
     printf("SystemMatrix_ALE constructed\n");
#endif
#ifndef _NSEONLY_
     printf("SystemMatrix_6DOF constructed\n");
#endif
    }
     
if(profiling)
   {
#ifdef _MPI
      end_int = MPI_Wtime();
#else
      end_int = GetTime();
#endif
      total_int = end_int-start_int;
    }

#ifndef _6DOFONLY_    
    if(profiling){
#ifdef _MPI
      start_assembling_solving = MPI_Wtime();
#else
      start_assembling_solving = GetTime();
#endif
    }          

    // interpolate initial solution
    u1->Interpolate(InitialU1);
    u2->Interpolate(InitialU2);
    u3->Interpolate(InitialU3);

    // assemble M, A, rhs matrices at time t=0 with initia sol and rhs
    SystemMatrix_ALE->Assemble();

#endif    

#ifdef _MPI
     if(TDatabase::ParamDB->WRITE_VTK)
      Output->Write_ParVTK(MPI_COMM_WORLD, img, SubID);
        img++;
#else
   if(TDatabase::ParamDB->WRITE_VTK)
     {
      os.seekp(std::ios::beg);
       if(img<10) os <<  "VTK/"<<VtkBaseName<<".0000"<<img<<".vtk" << ends;
         else if(img<100) os <<  "VTK/"<<VtkBaseName<<".000"<<img<<".vtk" << ends;
          else if(img<1000) os <<  "VTK/"<<VtkBaseName<<".00"<<img<<".vtk" << ends;
           else if(img<10000) os <<  "VTK/"<<VtkBaseName<<".0"<<img<<".vtk" << ends;
            else  os <<  "VTK/"<<VtkBaseName<<"."<<img<<".vtk" << ends;
      Output->WriteVtk(os.str().c_str());
      img++;
     }
#endif
     
//======================================================================
// time disc loop
//======================================================================    
   // parameters for time stepping scheme

   m = 0;
   N_SubSteps = GetN_SubSteps();
   oldtau = 1.;
   end_time = TDatabase::TimeDB->ENDTIME; 
   limit = TDatabase::ParamDB->SC_NONLIN_RES_NORM_MIN_SADDLE;
   Max_It = TDatabase::ParamDB->SC_NONLIN_MAXIT_SADDLE;        
   memset(AllErrors, 0, 7.*SizeOfDouble);

}


#ifndef _6DOFONLY_   

#ifdef _SMPI
  if(rank !=0)
   { 
    solmumps = new TMumpsSolver(0,0,NULL,NULL,0); // 
   }

  if(rank==0)
   {
    currtime = TDatabase::TimeDB->CURRENTTIME ;
   }


   MPI_Bcast(&currtime,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

   TDatabase::TimeDB->CURRENTTIME = currtime;
   MPI_Bcast(&end_time,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&N_SubSteps,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&Max_It,1,MPI_INT,0,MPI_COMM_WORLD);
#endif

#endif // _6DOFONLY_ 

   // =================================================================================      
   // time loop starts
   // both NSE and 6DOF eqns   
   // =================================================================================  
   
   while(TDatabase::TimeDB->CURRENTTIME< end_time)   // time cycle
    {          
     m++;
     TDatabase::TimeDB->INTERNAL_STARTTIME = TDatabase::TimeDB->CURRENTTIME;

     for(l=0;l<N_SubSteps;l++) // sub steps of fractional step theta
      {
        SetTimeDiscParameters(1);
        currtime = TDatabase::TimeDB->CURRENTTIME;
        tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
        TDatabase::TimeDB->CURRENTTIME += tau;

#if defined(_MPI) || defined (_SMPI)
   if(rank==0)
#endif
    {
	 if(m==1)
	  {
	   OutPut("Theta1: " << TDatabase::TimeDB->THETA1<< endl);
	   OutPut("Theta2: " << TDatabase::TimeDB->THETA2<< endl);
	   OutPut("Theta3: " << TDatabase::TimeDB->THETA3<< endl);
	   OutPut("Theta4: " << TDatabase::TimeDB->THETA4<< endl);
	  }

      OutPut(endl << "CURRENT TIME: ");
      OutPut(TDatabase::TimeDB->CURRENTTIME << endl);
     }

#ifdef _SMPI
   if(rank==0)
#endif
     {
#ifndef _6DOFONLY_       
      // copy rhs to oldrhs
      memcpy(oldrhs, rhs, N_TotalDOF*SizeOfDouble);  

    // assemble  matrix for NSE time-dependent domain
      SystemMatrix_ALE->Assemble(); 
 
      //scale B matices and assemble NSE-rhs based on the \theta time stepping scheme 
      SystemMatrix_ALE->AssembleSystMat(tau, rhs, rhs, sol); 
      oldtau = tau;
      
      // calculate the residual
      SystemMatrix_ALE->GetResidual(sol, impuls_residual,residual);
#endif // _6DOFONLY_

#ifdef _MPI
       if(rank==0)
#endif
       {     	
        OutPut(" nonlinear iteration step   0");
        OutPut(setw(14) << impuls_residual);
        OutPut(setw(14) << residual-impuls_residual);
        OutPut(setw(14) << sqrt(residual) << endl);  
       }

     } // #ifdef _SMPI if(rank==0)

//====================================================================== 
//Solve the NSE system
//Nonlinear iteration of fixed point type
//======================================================================
#ifndef _6DOFONLY_

     for(j=1;j<=Max_It;j++)
      {  
      // Solve the NSE system
	checker = 0;
#ifdef _SMPI
	if(rank==0)
#endif
	{
	  SystemMatrix_ALE->Solve(sol);
	}
#ifdef _SMPI
	else
	{  
	  solmumps->FactorizeAndSolve(NULL,NULL);
	}
#endif
 
     

#ifdef _SMPI
     if(rank==0)
#endif
      {
	
       if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
        IntoL20FEFunction3D(sol+3*N_U, N_P, Pressure_FeSpace[mg_level-1]);  
       
       //no nonlinear iteration for Stokes problem  
       if(TDatabase::ParamDB->FLOW_PROBLEM_TYPE==STOKES)
        break;

        //restore the mass matrix for the next nonlinear iteration      
        SystemMatrix_ALE->RestoreMassMatNonLinear();     
        //assemble the system matrix with given aux, sol and rhs 
        SystemMatrix_ALE->AssembleNonLinear();    
        //assemble system mat, S = M + dt\theta_1*A
        SystemMatrix_ALE->AssembleSystMatNonLinear();        
        //get the residual
        SystemMatrix_ALE->GetResidual(sol, impuls_residual,residual);
#ifdef _MPI
       if(rank==0)
#endif
       { 	
        OutPut(" nonlinear iteration step " << setw(3) << j);
        OutPut(setw(14) << impuls_residual);
        OutPut(setw(14) << residual-impuls_residual);
        OutPut(setw(14) << sqrt(residual) << endl);
       }
//       cout<< " Residual printing " <<endl;
    if(sqrt(residual)<=limit)
     {
#ifdef _SMPI
      checker = -1;
#else
       break;
#endif
      }
     } // #ifdef _SMPI  if(rank==0)

#ifdef _SMPI
   MPI_Bcast(&checker,1,MPI_INT,0,MPI_COMM_WORLD);
   if(checker ==-1)
     break;
#endif
   } // for(j=1;j<=Max_It;j++)
       


#endif // _6DOFONLY_ 

#ifndef _NSEONLY_
// ===========================================================
// 6dof solution in this substep and move the solid mesh
// ===========================================================
     SystemMatrix_6DOF->SolveAndUpdate(tau);
#endif
     
#ifndef _6DOFONLY_ 
     // move the fluid mesh
#ifdef _SMPI
if(rank==0)
#endif
  {    
     // thivin - commented
	     SystemMatrix_ALE->GetMeshVelo(TRUE, 1);

     // all mat will be assembled
     SystemMatrix_ALE->SetSystMatAssembled(FALSE);
   }    
#endif // _6DOFONLY_

//#ifdef _MPI
//     if(TDatabase::ParamDB->WRITE_VTK)
//      Output->Write_ParVTK(MPI_COMM_WORLD, img, SubID);
//        img++;
//#else
//     if(TDatabase::ParamDB->WRITE_VTK)
//      {
//       os.seekp(std::ios::beg);
//        if(img<10) os <<  "VTK/"<<VtkBaseName<<".0000"<<img<<".vtk" << ends;
//         else if(img<100) os <<  "VTK/"<<VtkBaseName<<".000"<<img<<".vtk" << ends;
//          else if(img<1000) os <<  "VTK/"<<VtkBaseName<<".00"<<img<<".vtk" << ends;
//           else if(img<10000) os <<  "VTK/"<<VtkBaseName<<".0"<<img<<".vtk" << ends;
//            else  os <<  "VTK/"<<VtkBaseName<<"."<<img<<".vtk" << ends;
//       Output->WriteVtk(os.str().c_str());
//       img++;
//      }
//#endif
//MPI_Finalize();
// exit(0);

   } // for(l=0;l<N_SubSteps;
      
//====================================================================== 
// measure errors to known solution
//======================================================================   
#ifdef _SMPI
if(rank==0)
#endif
{      
 #ifndef _6DOFONLY_     
    if(TDatabase::ParamDB->MEASURE_ERRORS)
     {   
      SystemMatrix_ALE->MeasureTNSEErrors(ExactU1, ExactU2, ExactU3, ExactP, AllErrors);
 #ifdef _MPI
       if(rank==0)
#endif
       { 
	  OutPut("L2(u): " <<   AllErrors[0] << endl);
	  OutPut("H1-semi(u): " <<  AllErrors[1] << endl);
	  OutPut("L2(p): " <<  AllErrors[2] << endl);
	  OutPut("H1-semi(p): " <<  AllErrors[3]   << endl); 
	  OutPut(AllErrors[4] <<  " l_infty(L2(u)) " <<AllErrors[5] << endl);
	  OutPut(TDatabase::TimeDB->CURRENTTIME << " L2(0,t,L2)(u) : " <<   sqrt(AllErrors[6]) << endl); 
       }
     } // if(TDatabase::ParamDB->MEASURE_ERRORS)
#endif //  _6DOFONLY_ 
//======================================================================
// produce output
//======================================================================     
     if(m==1 || m % TDatabase::TimeDB->STEPS_PER_IMAGE == 0)        
#ifdef _MPI
     if(TDatabase::ParamDB->WRITE_VTK)
      Output->Write_ParVTK(MPI_COMM_WORLD, img, SubID);
        img++;       
#else    	
     if(TDatabase::ParamDB->WRITE_VTK)
      {
       os.seekp(std::ios::beg);
        if(img<10) os <<  "VTK/"<<VtkBaseName<<".0000"<<img<<".vtk" << ends;
         else if(img<100) os <<  "VTK/"<<VtkBaseName<<".000"<<img<<".vtk" << ends;
          else if(img<1000) os <<  "VTK/"<<VtkBaseName<<".00"<<img<<".vtk" << ends;
           else if(img<10000) os <<  "VTK/"<<VtkBaseName<<".0"<<img<<".vtk" << ends;
            else  os <<  "VTK/"<<VtkBaseName<<"."<<img<<".vtk" << ends;
       Output->WriteVtk(os.str().c_str());
       img++;
      } 
#endif


    } // if(rank==0)
    
    
/*     cout << "6dof ssol" <<endl;
     exit(0);  */  
    
   } // while(TDatabase::TimeDB->CURRENTTIME< e
      
      
  CloseFiles();
#ifdef _MPI
  MPI_Finalize();
#endif
} // end main
