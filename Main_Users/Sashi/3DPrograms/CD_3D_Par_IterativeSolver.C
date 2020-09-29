// =======================================================================
//
// Purpose:     main program with parallel solver (no multigrid solver)
//
// Author:      Sashikumaar Ganesan
//
// History:     Implementation started on 06.03.2012
//               
// =======================================================================

#include <Domain.h>
#include <Database.h>
#include <DiscreteForm3D.h>
#include <FEDatabase3D.h>
#include <FESpace3D.h>
#include <SquareStructure3D.h>
#include <Structure3D.h>
#include <AuxParam3D.h>
#include <QuadAffin.h>
#include <Solver.h>
#include <DirectSolver.h>
#include <Assemble3D.h>
#include <Output3D.h>
#include <LinAlg.h>

#include <string.h>
#include <sstream>
#include <MooNMD_Io.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>

double bound = 0;
double timeC=0;
// ======================================================================
// utilities for main program
// ======================================================================
#include <MainUtilities.h>
#include <Upwind3D.h>
#include <FEM_TVD_FCT.h>
#include <MultiGrid3D.h>
// #include <MGLevel3D.h>
#include <ItMethod.h>
#include <FixedPointIte.h>
#include <FgmresIte.h>
#include <JacobiIte.h>
// #include <MultiGridScaIte.h>


#include <tetgen.h>
#include <GridCell.h>
#include <MacroCell.h>
#include <BdPlane.h>
#include <BdSphere.h>
#include <IsoBoundFace.h>
#include <InterfaceJoint3D.h>
#include <IsoInterfaceJoint3D.h>

#include <sys/stat.h>
#include <sys/types.h>

#ifdef _MPI
#include <MeshPartition.h>
#include <ParFECommunicator3D.h>
#include <MumpsSolver.h>
#include <ParVector3D.h>
#include <ParVectorNSE3D.h>
#include <Scalar_ParSolver.h>
#endif

#define AMG 0
#define GMG 1
#define DIRECT 2

// =======================================================================
// include current example
// =======================================================================
#include "../Examples/CD_3D/Laplace.h"
// #include "../Examples/CD_3D/LaplaceBsp1_2.h"
//#include "../Examples/CD_3D/ThreeBoundaryLayers.h"
// #include "../Examples/CD_3D/Sphere.h"
// #include "../Examples/CD_3D/Constant1.h"
// #include "../Examples/CD_3D/Cylinder.h"
// #include "../Examples/CD_3D/Plane.h"
//#include "../Examples/CD_3D/CircularLayer.h"
//#include "../Examples/CD_3D/SkewConv.h"
//#include "../Examples/CD_3D/ParabolicLayers3D.h"
// #include "../Examples/CD_3D/Hemker1996_3D.h"
//#include "../Examples/TCD_3D/bulk_compare.h"
// #include "../Examples/TCD_3D/TeraHertzBrain.h"
// #include "../Examples/TCD_3D/TeraHertzCube.h"
// #include "../Examples/CD_3D/Sine.h"
// #include "../Examples/TCD_3D/TeraHertzBreast.h"

int main(int argc, char* argv[])
{
  TDatabase *Database = new TDatabase();    
  double hmin, hmin_all, hmax, start_time, time1, t1, t2;  
  
#ifdef _MPI
  const int root = 0;
  int rank, size;
  int MaxCpV, MaxSubDomainPerDof;
  int  out_rank, N_Cells_loc;

  TParVector3D  *ParSolVect, *ParRhsVect, *ParBVect;
  TParFECommunicator3D *ParComm;

  MPI_Request request001, request002, request003, request004, request005;
  MPI_Status status;

  MPI_Comm Comm = MPI_COMM_WORLD;
  MPI_Init(&argc, &argv);

  time1 =  MPI_Wtime();
  MPI_Allreduce(&time1, &start_time, 1, MPI_DOUBLE, MPI_MIN, Comm);

  MPI_Comm_rank(Comm, &rank);
  MPI_Comm_size(Comm, &size);
 
  TDatabase::ParamDB->Comm = Comm;
  TScalar_ParSolver *Par_Solver; 
#endif  
    
  TDomain *Domain = new TDomain(); 
  TDomain *Domain_Sptl = new TDomain();
    
  TFEDatabase3D *FEDatabase = new TFEDatabase3D();
  TCollection *coll, *coll_Sptl, *mortarcoll = NULL;
  TBaseCell *cell;
  TFESpace3D *concentration_space, *pressure_space, *velocity_space; 
  TOutput3D *Output; 
  TAuxParam3D *aux;
  TFEFunction3D *conc;
  TSquareStructure3D *sqstructureA;  
  TSquareMatrix3D *sqmatrixA, *SQMATRICES[3];;
  TSquareMatrix **sqmatrices = (TSquareMatrix **)SQMATRICES;
  TFESpace3D *fesp[2], *ferhs[1];  
  TJacobiIte *jacob;
  MatVecProc *MatVect;
  DefectProc *Defect;  
 
  BoundCondFunct3D *BoundaryConditions[1];
  BoundValueFunct3D *BoundValues[1];
  CoeffFct3D *Coefficients[1];
  TDiscreteForm3D *DiscreteFormGalerkin;
  TDiscreteForm3D *DiscreteForm;  
  
  int i, ret, ORDER, N_Unknowns, N_Active, N_LinIter;  
  int N_SquareMatrices, N_FESpaces, N_Rhs, img=0;  
  
  double *RHSs[3], *sol, *rhs, *B, norm[2], hK, hK_all;
  double x, y, z, norm_prev=0;
   
  bool FACTORIZE=TRUE;
  
  char UString[] = "T";
  char NameString[] = "name";  
  char ReadinDat[] = "readin.dat";
  const char BDdir[] = "BDData";
  const char vtkdir[] = "VTK";  
  char SubID[] = "";
  
  char *PsBaseName, *VtkBaseName, *GmvBaseName;
  char *PRM, *GEO;
  
  std::ostringstream os, opts;
  os << " ";
  opts << " ";
  
  mkdir(vtkdir, 0777);
  mkdir(BDdir, 0777);      
 
//======================================================================
// read parameter file
//======================================================================
  if(argc>=2)
    ret=Domain->ReadParam(argv[1]);
  else
    ret=Domain->ReadParam(ReadinDat);

  if(ret==-1)
    exit(-1);
 
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
  
  PRM = TDatabase::ParamDB->BNDFILE;
  GEO = TDatabase::ParamDB->GEOFILE;
  
  PsBaseName = TDatabase::ParamDB->PSBASENAME;  
  VtkBaseName = TDatabase::ParamDB->VTKBASENAME;
  GmvBaseName  = TDatabase::ParamDB->GMVBASENAME;
  
  MatVect = MatVect_Scalar;
  Defect = Defect_Scalar;
  
#ifdef _MPI
  if(TDatabase::ParamDB->SOLVER_TYPE<100)
   TDatabase::ParamDB->SOLVER_TYPE = 101;   
#else
   TDatabase::ParamDB->SOLVER_TYPE = 2;
#endif

//======================================================================
// initialize discrete forms
//======================================================================
   InitializeDiscreteForms(DiscreteFormGalerkin, BilinearCoeffs);  

   
//======================================================================
// read boundary parameterization and initialize coarse grid
//======================================================================
  Domain->Init(PRM, GEO);

  
    /** with atlas mesh, no tetgen*/
//     TetrameshCreate(Domain);
    
   /** Using tetgen with smesh mesh */
//    TetrameshGen(Domain);
      
   /** same as TetrameshGen but without using the face list info from tetgen */
// // // // // // // //     TetraGen(Domain);
    
//   TetraGenWithInputCells(Domain);
//   Domain->ConvertQuadToTri(TDatabase::ParamDB->CONVERT_QUAD_TO_TRI);
// exit(0);
   
 
  //======================================================================
  // Partition grid using Metis
  //======================================================================
#ifdef _MPI
  for(i=0;i<4;i++)
    Domain->RegRefineAll();  
  Domain->GenerateEdgeInfo();
   
  t1 = MPI_Wtime();
  Partition_Mesh3D(Comm, Domain, MaxCpV);
  t2 = MPI_Wtime(); 
  
  for(i=4;i<TDatabase::ParamDB->UNIFORM_STEPS;i++)
  {
    printf("************************************LEVEL %d****************************************",i);
    Domain->RegRefineAll();
    Domain->GenerateEdgeInfo();
    Domain_Crop(Comm, Domain);
  }
  if(rank==0)
    printf("Time taken for Domain Decomposition is %e\n", (t2-t1));
   Domain->GenerateEdgeInfo();
  MaxSubDomainPerDof = MIN(MaxCpV, size);

#else
  
  for(i=0;i<TDatabase::ParamDB->UNIFORM_STEPS;i++)
    Domain->RegRefineAll();  
#endif   
      

//======================================================================
// construct all fespace
//====================================================================== 
// create collection of mesh cells
#ifdef _MPI
   //coll=Domain->GetOwnCollection(It_Finest, 0, rank);
   coll=Domain->GetCollection(It_Finest, 0);
#else
   coll=Domain->GetCollection(It_Finest, 0);
#endif 
   
   OutPut(" N_Cells: "<< setw(10) << coll->GetN_Cells() << endl);
   
   ORDER  = TDatabase::ParamDB->ANSATZ_ORDER; 

   concentration_space = new TFESpace3D(coll, NameString, UString, BoundCondition, ORDER);
   N_Unknowns = concentration_space->GetN_DegreesOfFreedom();
 
   // active dof (i.e. dof without Dirichlet dofs)
   N_Active = concentration_space->GetActiveBound();
   OutPut(N_Active << " Degrees of freedom: "<< setw(10) << N_Unknowns << endl);
//    exit(0);
#ifdef _MPI
    t1 = MPI_Wtime();

    concentration_space->SetMaxSubDomainPerDof(MaxSubDomainPerDof);
    ParComm = new TParFECommunicator3D(Comm, concentration_space);
    
    t2 = MPI_Wtime();
    if(rank==out_rank)
     {
      printf("Time taken for FeSpace SubDomain dof mapping %e\n", (t2-t1));
      printf("DOF of FeSpace  space : %d \n", ParComm->GetN_GlobalDegreesOfFreedom());      
     }
#else
    OutPut(" Rank " <<  " DOF Scalar : " << setw(10) << N_Unknowns << endl);
#endif 
//======================================================================
// allocate memory for all matrices
//======================================================================  
   // first build matrix structure
   sqstructureA = new TSquareStructure3D(concentration_space);
   sqstructureA->Sort();
  
   sqmatrixA = new TSquareMatrix3D(sqstructureA);
//======================================================================
// allocate memory for sol and rhs arrays and construct fefunction
//======================================================================  
   sol = new double[N_Unknowns]; 
   rhs = new double[N_Unknowns];
   B = new double[N_Unknowns];

#ifdef _MPI
    ParSolVect =  new TParVector3D(Comm, sol, N_Unknowns, 1, ParComm);
    ParRhsVect =  new TParVector3D(Comm, rhs, N_Unknowns, 1, ParComm);  
    ParBVect =  new TParVector3D(Comm, B, N_Unknowns, 1, ParComm);      
#endif

   memset(B, 0, N_Unknowns*SizeOfDouble); 
   memset(rhs, 0, N_Unknowns*SizeOfDouble);  
   memset(sol, 0, N_Unknowns*SizeOfDouble);

  
   // allocate fe function 
   conc = new TFEFunction3D(concentration_space, UString, UString, sol, N_Unknowns); 

   // interpolate initial condition
//    conc->Interpolate(InitialCondition);
//   OutPut("Interpolation time: " << TDatabase::TimeDB->CURRENTTIME << endl);   

#ifdef _MPI
     // initialize the parallel solver
    // Par_Solver = new TScalar_ParSolver(ParComm, sqstructureA, 1, ParSolVect, ParRhsVect);
#endif   

/*
//======================================================================
// prepare output,
//======================================================================     
   Output = new TOutput3D(1, 1, 0, 1, Domain, coll);
   Output->AddFEFunction(conc);
   os.seekp(std::ios::beg);
   
#ifdef _MPI     
     //  Output->Write_ParVTK(Comm, img, SubID);
#else     
       os.seekp(std::ios::beg);      
       if(img<10) os << "VTK/"<<VtkBaseName<<".0000"<<img<<".vtk" << ends;
       else if(img<100) os << "VTK/"<<VtkBaseName<<".000"<<img<<".vtk" << ends;
       else if(img<1000) os << "VTK/"<<VtkBaseName<<".00"<<img<<".vtk" << ends;
       else if(img<10000) os << "VTK/"<<VtkBaseName<<".0"<<img<<".vtk" << ends;
       else  os << "VTK/"<<VtkBaseName<<"."<<img<<".vtk" << ends;
       Output->WriteVtk(os.str().c_str());
       
/*     os.seekp(std::ios::beg);      
       if(img<10) os << "VTK/"<<VtkBaseName<<".0000"<<img<<".plt" << ends;
       else if(img<100) os << "VTK/"<<VtkBaseName<<".000"<<img<<".plt" << ends;
       else if(img<1000) os << "VTK/"<<VtkBaseName<<".00"<<img<<".plt" << ends;
       else if(img<10000) os << "VTK/"<<VtkBaseName<<".0"<<img<<".plt" << ends;
       else  os << "VTK/"<<VtkBaseName<<"."<<img<<".plt" << ends;
       Output->WriteBinaryPlt(os.str().c_str());   
       img++;

 #endif 
   */    

//======================================================================      
// assemble all matrices
//======================================================================   
#ifdef _MPI    
         t1 = MPI_Wtime();
#else  
         t1 = GetTime();
#endif   
   
    N_FESpaces = 1;
    fesp[0] = concentration_space;
   
    // reset matrices
    N_SquareMatrices = 1;
    SQMATRICES[0] = sqmatrixA;
    SQMATRICES[0]->Reset();   
   
    DiscreteForm = DiscreteFormGalerkin; 
    BoundaryConditions[0] =  BoundCondition;
    BoundValues[0] = BoundValue;
   
    N_Rhs = 1;
    memset(rhs, 0, N_Unknowns*SizeOfDouble);  
    RHSs[0] = rhs;
    ferhs[0] = concentration_space;   
   
    aux =  new TAuxParam3D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);  
   
    //measure hK, see Example file THz Brain
//     TDatabase::ParamDB->P14 = -1.;

 
    Assemble3D(N_FESpaces, fesp, 
               N_SquareMatrices, SQMATRICES, 
               0, NULL, 
               N_Rhs, RHSs, ferhs,
               DiscreteForm, 
               BoundaryConditions, 
               BoundValues, 
               aux);     
#ifdef _MPI
     ParComm->SetSlaveDofRows(SQMATRICES[0]->GetRowPtr(), SQMATRICES[0]->GetKCol(), SQMATRICES[0]->GetEntries(), rhs);  
#endif   
   
  
     if(TDatabase::ParamDB->SC_VERBOSE>1)
      {
#ifdef _MPI    
       t2 = MPI_Wtime();
#else  
       t2 = GetTime();
#endif   
#ifdef _MPI      
       if(rank==out_rank)      
#endif               
         OutPut("Time for assembling " << t2-t1 << "s "<<  endl);   
        }
        
      /** copy Dirichlet values from rhs into sol */
      memcpy(sol+N_Active, RHSs[0]+N_Active, (N_Unknowns-N_Active)*SizeOfDouble);
      
    
      
    hK_all = TDatabase::ParamDB->P14;
#ifdef _MPI         
   MPI_Allreduce(&hK_all, &hK, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#else
   hK = hK_all;
#endif 
   
//======================================================================
// Directional Domain
//======================================================================   
 #ifdef _RTE_ 
  Domain_Sptl->Init(PsBaseName, GmvBaseName);   
#ifdef _MPI
   coll_Sptl=Domain_Sptl->GetOwnCollection(It_Finest, 0, rank);
#else
   coll_Sptl=Domain_Sptl->GetCollection(It_Finest, 0);
#endif  

   OutPut(" N_Cells_Sptl: "<< setw(10) << coll_Sptl->GetN_Cells() << endl);       
#endif      
   

//======================================================================    
   
// #ifdef _MPI     
//    if(rank==out_rank)      
// #endif   
//     OutPut("Cell Diameter (within Sourceradius): " <<  hK <<  endl);   
 

//       if(rank==0)
//         for(i=N_Active; i<N_Unknowns; i++)
//           OutPut(i<<" rhs " <<   rhs[i]  <<endl);      

 
//           memset(rhs, 0, N_Unknowns*SizeOfDouble);  
   

//======================================================================
// solve the system
//======================================================================  
    
#ifdef _MPI 
        jacob=new TJacobiIte(NULL, NULL, NULL, 0, N_Unknowns, 1,ParComm);
        ParSolVect->ParDdot(BYOWN, norm);
      
        if(rank==out_rank)
         OutPut("Initial Norm of Par sol " << sqrt( norm[0]) <<endl);   
#else   
	 jacob=new TJacobiIte(NULL, NULL, NULL, 0, N_Unknowns, 1);
        OutPut("Norm of sol " <<  sqrt(Ddot(N_Unknowns,sol,sol)) <<endl); 
#endif  
	
#ifndef _MPI
        t1 = GetTime();	
	jacob->Iterate_p(sqmatrices,NULL,sol,rhs);
         //   DirectSolver(sqmatrixA, rhs, sol);
        t2 = GetTime();     
#else
        t1 = MPI_Wtime();
        jacob->Iterate_p(sqmatrices,NULL,sol,rhs,ParComm);
        t2 = MPI_Wtime();
	//test
//         memset(B, 0, N_Unknowns*SizeOfDouble);  
//         memcpy(B+N_Active, sol+N_Active, (N_Unknowns-N_Active)*SizeOfDouble);	

//      if(rank==0)
//          {    
//            i  = 220;
//            printf("rank %d RHS %d rhs %e  \n", rank, i,   rhs[i] );
// 
// 	  concentration_space->GetDOFPosition(i,  x,  y,  z);
// 	      printf("rank %d i %d x %e y  %e z  %e  \n", rank, i,  x, y, z );   
//           }  
/*     if(rank==1)
         {    
           i  = 238;
           printf("rank %d RHS %d rhs %e  \n", rank, i,   rhs[i] );

	  concentration_space->GetDOFPosition(i,  x,  y,  z);
	      printf("rank %d i %d x %e y  %e z  %e  \n", rank, i,  x, y, z );   
          } */	
         
       // Par_Solver->Solve(sqmatrixA, FACTORIZE);
   
#endif  


//#ifdef _MPI   

//           Daxpy((N_Unknowns-N_Active), -1.0, sol+N_Active, B+N_Active);	

	  
// 	if(rank==TDatabase::ParamDB->Par_P8)  
// 	{
//            printf("N_Unknowns %d N_Active %d \n", N_Unknowns, N_Active);  
// 	   
// 	  for(i=N_Active;i<N_Unknowns;i++)
// 	    if(fabs(B[i]-sol[i])>1.e-8)
// 	    printf("i %d rhs %e sol %e rhs- sol %e \n", i,  B[i], sol[i], B[i]-sol[i]);
// 	  
// 	}
       //testing
//         ParSolVect->AssembleWithMaster();  
   
  
// 	ParBVect->ParDdot(BYOWN, norm);
#ifdef _MPI 
        ParSolVect->ParDdot(BYOWN, norm);
      
        if(rank==out_rank)
         OutPut("Norm of Par sol " << sqrt( norm[0])  << " time taken "<<(t2-t1)<<endl);   
#else     
        OutPut("Norm of sol " <<  sqrt(Ddot(N_Unknowns,sol,sol))  << " time taken "<<(t2-t1)<<endl); 
#endif  
	
// #ifdef _MPI    
//    MPI_Finalize();  
// #endif   
//    exit(0);
 
//======================================================================
// Output solution
//======================================================================  
/*   
#ifdef _MPI     
       Output->Write_ParVTK(Comm, img, SubID);
#else     
       os.seekp(std::ios::beg);      
       if(img<10) os << "VTK/"<<VtkBaseName<<".0000"<<img<<".vtk" << ends;
       else if(img<100) os << "VTK/"<<VtkBaseName<<".000"<<img<<".vtk" << ends;
       else if(img<1000) os << "VTK/"<<VtkBaseName<<".00"<<img<<".vtk" << ends;
       else if(img<10000) os << VtkBaseName<<".0"<<img<<".vtk" << ends;
       else  os << "VTK/"<<VtkBaseName<<"."<<img<<".vtk" << ends;
       Output->WriteVtk(os.str().c_str());
       
/*     os.seekp(std::ios::beg);      
       if(img<10) os << "VTK/"<<VtkBaseName<<".0000"<<img<<".plt" << ends;
       else if(img<100) os << "VTK/"<<VtkBaseName<<".000"<<img<<".plt" << ends;
       else if(img<1000) os << "VTK/"<<VtkBaseName<<".00"<<img<<".plt" << ends;
       else if(img<10000) os << "VTK/"<<VtkBaseName<<".0"<<img<<".plt" << ends;
       else  os << "VTK/"<<VtkBaseName<<"."<<img<<".plt" << ends;
       Output->WriteBinaryPlt(os.str().c_str());     
       img++;
 #endif       
            
  */     
//   OutPut("used time: " << GetTime() << endl);
//   OutPut("used bytes: " << GetMemory() << endl);
#ifdef _MPI    
   MPI_Finalize();  
#endif           
  CloseFiles();

  return 0;
}


 
