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
#include <omp.h>
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
#include <MGLevel3D.h>
#include <ItMethod.h>
#include <FixedPointIte.h>
#include <FgmresIte.h>
#include <JacobiIte.h>
#include <MultiGridScaIte.h>

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
//#include <MeshPartition2D.h>
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
  double t_par1, t_par2;
  TParVector3D  *ParSolVect, *ParRhsVect, *ParBVect;
  TFESpace3D **OwnScalar_Spaces, *ownscalar_space; 
  TParFECommunicator3D **ParComm;

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
  TCollection *coll, *coll_Sptl, *own_coll, *mortarcoll = NULL;
  TBaseCell *cell;
  TFESpace3D *scalar_space, *scalar_space_low, *RefinedOutput_space[2], **Scalar_Spaces; 
  TOutput3D *Output; 
  TAuxParam3D *aux;
  TFEFunction3D *conc;
  TSquareStructure3D *sqstructureA,*sqstructureA_low;
  TSquareMatrix3D *sqmatrixA, *sqmatrixA_low, *SQMATRICES[1];
  TSquareMatrix3D **MatricesA, *Tmp_MatricesA;
  TSquareMatrix **sqmatrices = (TSquareMatrix **)SQMATRICES;
  TFESpace3D *fesp[2], *ferhs[1],*concentration_space,*velocity_space;  
  TMGLevel3D *MGLevel, *MGLevel_low;
  TMultiGrid3D *MG;
  TItMethod *itmethod, *prec; 
  TJacobiIte *jacob;
  MatVecProc *MatVect;
  DefectProc *Defect;  
 

  CoeffFct3D *Coefficients[1];
  TDiscreteForm3D *DiscreteFormUpwind, *DiscreteFormGalerkin, *DiscreteFormSDFEM, *DiscreteFormGLS;
  TDiscreteForm3D *DiscreteForm; 
  TFEFunction3D *C, *C_low, *Refinmed_C, **FEFunctArray;
 // MultiIndex3D AllDerivatives[3] = { D00, D10, D01 };
  
  int i,j,k, ret, ORDER, N_Unknowns, N_Active, N_LinIter, N_DOF, N_DOF_low,N_NonActive, *N_DOFs;  
  int N_SquareMatrices, N_FESpaces, N_Rhs, img=0, zerostart;  
  int LEVELS, N_U_DOF, N_BData=0;
  int mg_level, mg_type, n_aux, N_Paramters=2, FirstSolve, low;
  
  double solver_time;
  double *sol, *oldsol, *rhs, *defect, *PostSol, *refined_sol[2], *auxsol;
  double *rhs_low, *sol_low,  *RHSs[1], errors[4], **RhsArray;
  double *l2, *h1, *sd;
  double *u_sol, *itmethod_sol, *itmethod_rhs,norm[2];
  double Parameters[2];
   
  bool FACTORIZE=TRUE;
  
 char *PRM, *GEO;
  char *PsBaseName, *VtkBaseName;
  char CString[] = "c";
  char UString[] = "u";
  char RUString[] = "ref_u";
  char PString[] = "p";
  char ReadinDat [] = "readin.dat";
  char MassMatrix[] = "Mass matrix";
  char Mass[] = "Mass";
  char Name[] = "Temp";
  char UName[] = "u";
  char Description[] = "description";
  char CdString[] = "Conv-Diff";
  char GalString[] = "Galerkin";
  char SDFEMString[] = "SDFEM";
  char UpwString[] = "Upwind";
  const char BDdir[] = "BDData";
  const char vtkdir[] = "VTK";
  char SubID[] = "";
  char NameString[] = "name"; 
  
  boolean ForPressure=FALSE, ReadVelo=FALSE;
   
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
     ExampleFile();
   }
   
 
 
//======================================================================
// copy read parameters into local variables
//======================================================================
  PRM = TDatabase::ParamDB->BNDFILE;
  GEO = TDatabase::ParamDB->GEOFILE;
  PsBaseName = TDatabase::ParamDB->PSBASENAME;
  VtkBaseName = TDatabase::ParamDB->VTKBASENAME;

  // set type of multilevel
  mg_type = TDatabase::ParamDB->SC_MG_TYPE_SCALAR;
 
  if(TDatabase::ParamDB->SOLVER_TYPE==AMG|| TDatabase::ParamDB->SOLVER_TYPE==DIRECT)
       mg_type=0;
  
  if(mg_type)
    mg_level = 0;
  else
    mg_level = -1;
   
  LEVELS = TDatabase::ParamDB->LEVELS;
//   mg_level += LEVELS;

  // initialize multigrid
  if(TDatabase::ParamDB->SOLVER_TYPE == GMG)
  {
    Parameters[0] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_SCALAR;
    Parameters[1] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_FINE_SCALAR;
    i=1;
    MG = new TMultiGrid3D(i, N_Paramters, Parameters);
    
    switch (TDatabase::ParamDB->SC_SOLVER_SCALAR)
      {
          case 11:
           zerostart = 1;
          break;
          case 16:
           zerostart = 1;
          break;
      }    
    }  
  
  FirstSolve = TDatabase::ParamDB->SC_FIRST_SOLUTION_LEVEL_SCALAR;
  
  
  // finite element spaces on the different levels of the multigrid

  Scalar_Spaces = new TFESpace3D*[LEVELS+1];  
  FEFunctArray = new TFEFunction3D*[LEVELS+1];
  MatricesA = new TSquareMatrix3D*[LEVELS+1];
  
  N_DOFs = new int[LEVELS+1];  
  l2 = new double[LEVELS+1];
  h1 = new double[LEVELS+1];
  sd = new double[LEVELS+1];
  RhsArray = new double*[LEVELS+1];  
  
#ifdef _MPI    
  OwnScalar_Spaces = new TFESpace3D*[LEVELS+1];   
  ParComm = new TParFECommunicator3D*[LEVELS+1];  
  
  out_rank=TDatabase::ParamDB->Par_P0;  
#endif   
  
  

//======================================================================
// initialize discrete forms
//======================================================================
  // InitializeDiscreteForms_Stationary(DiscreteFormUpwind, DiscreteFormGalerkin, DiscreteFormSDFEM, DiscreteFormGLS,
    //                                 BilinearCoeffs);
   InitializeDiscreteForms(DiscreteFormGalerkin, BilinearCoeffs);  
 
//======================================================================
// generate mesh and refun (if needed)
//======================================================================
    Domain->Init(PRM, GEO);
    
      /** with atlas mesh, no tetgen*/
//     TetrameshCreate(Domain);
    
   /** Using tetgen with smesh mesh */
//    TetrameshGen(Domain);

/* generate special mesh for Hemker example */
#ifdef __HMM_1986__
    if(!strcmp(GEO, "InitGrid"))
     if(TDatabase::ParamDB->REACTOR_P25)
         MeshReGen_HemkerResolved(Domain);
     else
         MeshReGen_Hemker(Domain);
#endif

  // refine grid
//   for(i=0;i<TDatabase::ParamDB->UNIFORM_STEPS;i++)
// 	      Domain->RegRefineAll();
// //        Domain->AdaptRefineAll();  
   

//=================================================
// STEP >>>>>>> PARTITION GRID USING METIS
//=================================================
#ifdef _MPI
  for(i=0;i<TDatabase::ParamDB->UNIFORM_STEPS;i++)
    Domain->RegRefineAll();  
  Domain->GenerateEdgeInfo();
   
  t1 = MPI_Wtime();
  Partition_Mesh3D(Comm, Domain, MaxCpV);
  t2 = MPI_Wtime(); 
  
 /* for(i=3;i<TDatabase::ParamDB->UNIFORM_STEPS;i++)
  {
    printf("************************************LEVEL %d****************************************",i);
    Domain->RegRefineAll();
    Domain->GenerateEdgeInfo();
    Domain_Crop(Comm, Domain);
  }*/
  if(rank==0)
    printf("Time taken for Domain Decomposition is %e\n", (t2-t1));
   Domain->GenerateEdgeInfo();
  MaxSubDomainPerDof = MIN(MaxCpV, size);

#else
  
  for(i=0;i<TDatabase::ParamDB->UNIFORM_STEPS;i++)
    Domain->RegRefineAll();  
#endif 
     
//======================================================================
// include the boundary condition and boundary values from the example file
//======================================================================
  BoundCondFunct3D *BoundaryConditions[1] = { BoundCondition };
  BoundValueFunct3D *BoundValues[1] = { BoundValue };

  ORDER = TDatabase::ParamDB->ANSATZ_ORDER;
    
//======================================================================
// loop over all levels (not a multigrid level but for convergence study)
//======================================================================
//  OutPut("Fixed Point Comm " << Comm << endl);
  for(i=0;i<LEVELS;i++)
  { 
    mg_level++;
    
#ifdef _MPI  
   if(rank==0)
   {
#endif      
    OutPut("*******************************************************" << endl);
    OutPut("******           GEOMETRY  LEVEL ");
    OutPut(i << "              *******" << endl);
    OutPut("******           MULTIGRID LEVEL ");
    OutPut(mg_level << "              *******" << endl);
    OutPut("*******************************************************" << endl);
//     OutPut("memory before: " << setw(10) << GetMemory() << endl);
#ifdef _MPI 
   }
#endif   
    solver_time = 0.0;
    N_LinIter = 0;  

     if(i)
     {
      Domain->RegRefineAll();
     
     }

#ifdef _MPI
     if(i)
     {
       Domain->GenerateEdgeInfo();
       Domain_Crop(Comm, Domain);       // removing unwanted cells in the hallo after refinement
     }
#endif
 
    coll=Domain->GetCollection(It_Finest, 0);

#ifdef _MPI  
    own_coll = Domain->GetOwnCollection(It_Finest, 0, rank);	 
#endif

//======================================================================
// construct all finite element spaces
//======================================================================

    // if multiple discretization multilevel method is used
    // get space for low order disc on finest geo grid
    if(mg_type==1)
    {
     // nonconforming space, including hallo in parallel
     scalar_space_low = new TFESpace3D(coll, Name, Description,BoundCondition, -1);
     Scalar_Spaces[i] = scalar_space_low;
     N_DOF_low = scalar_space_low->GetN_DegreesOfFreedom();
     N_DOFs[i] = N_DOF_low;         
        
#ifdef _MPI
    t1 = MPI_Wtime();

    sqstructureA_low = new TSquareStructure3D(scalar_space_low);
    sqstructureA_low->Sort();
      
    Scalar_Spaces[i]->SetMaxSubDomainPerDof(MaxSubDomainPerDof);
    ParComm[i] = new TParFECommunicator3D(Comm, Scalar_Spaces[i], sqstructureA_low);
  
    t2 = MPI_Wtime();
    if(rank==out_rank)
     {
      printf("Time taken for FeSpace_low SubDomain dof mapping %e\n", (t2-t1));
     // printf("DOF of FeSpace  space_low : %d \n", ParComm[i]->GetN_GlobalDegreesOfFreedom());      
     }
 
     // fespace on own cels
     ownscalar_space = new TFESpace3D(own_coll,  Name, Description, BoundCondition, -1); 
     OwnScalar_Spaces[i] = ownscalar_space;
#else
     OutPut("N_Dof_low   : "<<  setw(10) << N_DOF_low  << endl);  
    
#endif    
    } // if(mg_type==1)

    // get fe space of high order disc on finest geo grid
    if((i>=FirstSolve)||(mg_type==0))
    { 
     scalar_space =  new TFESpace3D(coll, Name, Description, BoundCondition, ORDER);
     Scalar_Spaces[mg_level] = scalar_space;
     N_DOF = scalar_space->GetN_DegreesOfFreedom();
     N_DOFs[mg_level] = N_DOF;  

#ifdef _MPI
    t1 = MPI_Wtime();

     sqstructureA = new TSquareStructure3D(scalar_space);
     sqstructureA->Sort();
     
    Scalar_Spaces[mg_level]->SetMaxSubDomainPerDof(MaxSubDomainPerDof);
    ParComm[mg_level] = new TParFECommunicator3D(Comm, Scalar_Spaces[mg_level], sqstructureA);
  
    t2 = MPI_Wtime();
    if(rank==out_rank)
     {
      printf("Time taken for FeSpace SubDomain dof mapping %e\n", (t2-t1));
   //   printf("DOF of FeSpace  space  : %d \n", ParComm[mg_level]->GetN_GlobalDegreesOfFreedom());      
     }
 
     // fespace on own cells
     ownscalar_space = new TFESpace3D(own_coll,  Name, Description, BoundCondition, ORDER); 
     OwnScalar_Spaces[mg_level] = ownscalar_space;
#else
      OutPut("N_Dof       : "<<  setw(10) << N_DOF  << endl);     
#endif       
    }  
    
    if(ReadVelo)
      velocity_space =  new TFESpace3D(coll, UName, UString, BoundCondition, 2);
     
 //======================================================================
// construct all finite element functions
//======================================================================
    if (mg_type==1)
     {
       
//       OutPut("Dof_low       : "<<  setw(10) << N_DOF_low  << endl);

      rhs_low = new double[N_DOF_low];
      memset(rhs_low, 0, N_DOF_low*SizeOfDouble);
      RhsArray[i] = rhs_low;
      
      sol_low = new double[N_DOF_low];
      memset(sol_low, 0, N_DOF_low*SizeOfDouble);
      C_low = new TFEFunction3D(scalar_space_low, UString, UString, sol_low, N_DOF_low);   
      FEFunctArray[i] = C_low;      
     }
           
    // high order disc
    N_Unknowns = N_DOF;  
//     OutPut("dof all      : "<< setw(10) << N_Unknowns  << endl);
     
    if ((i>=FirstSolve)||(mg_type==0))
    {     
     sol = new double[N_Unknowns];
     memset(sol, 0, N_Unknowns*SizeOfDouble);      
     rhs = new double[N_Unknowns];
     memset(rhs, 0, N_Unknowns*SizeOfDouble);   
     RhsArray[mg_level] =  rhs;        

     C = new TFEFunction3D(scalar_space, UString, UString, sol, N_DOF);
     FEFunctArray[mg_level] = C;   
#ifdef _MPI
    ParSolVect =  new TParVector3D(Comm, sol, N_Unknowns, 1, ParComm[i]);
#endif
//======================================================================
// produce outout
//======================================================================
   // prepare output, only the concentration will be saved

    Output = new TOutput3D(2, 2, 1, 1,Domain);
    
    Output->AddFEFunction(C);
    
  /*  if(i==FirstSolve)    
    {
      C->Interpolate(Exact);
#ifdef _MPI
      ParSolVect->ParDdot(BYOWN, norm);
#endif
        if(rank==out_rank)
         OutPut("Exact Norm of Par sol " << sqrt( norm[0]) <<endl); 
    }*/
#ifdef _MPI
     t_par1 = MPI_Wtime();

     if(TDatabase::ParamDB->WRITE_VTK)
      Output->Write_ParVTK(Comm, img, SubID);
      img++;      
     t_par2 = MPI_Wtime();

     if(rank==out_rank)
     printf("Time taken for writing the parvtk file %e\n", (t_par2-t_par1));
     
//      exit(0);
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
     }   
        img++; 
#endif     
    } // if ((i>=FirstSolve)||(mg_type==0))       
    
    memset(sol, 0, N_Unknowns*SizeOfDouble); 
//======================================================================
// allocate memory for all matrices
//======================================================================
    // build matrices for low order disc
    if(mg_type==1)
    {
      // matrix structures
      sqstructureA_low = new TSquareStructure3D(scalar_space_low);
      sqstructureA_low->Sort();
      
      sqmatrixA_low = new TSquareMatrix3D(sqstructureA_low);
      MatricesA[i] = sqmatrixA_low;
    } 

    // build matrices for high order disc
    if ((i>=FirstSolve)||(mg_type==0))
    {
     sqstructureA = new TSquareStructure3D(scalar_space);
     sqstructureA->Sort();

     sqmatrixA = new TSquareMatrix3D(sqstructureA);
     MatricesA[mg_level] = sqmatrixA;
    } 
    

    
//======================================================================    
// prepare multigrid method
//======================================================================   
    switch(TDatabase::ParamDB->SOLVER_TYPE)
    {
      case AMG:
      case DIRECT:
        low = mg_level;
        break;

      case GMG:
        // coarsest grid number
        low = 0;
        // determine number of auxiliary arrays
         if ((TDatabase::ParamDB->SC_STEP_LENGTH_CONTROL_ALL_SCALAR)
              || (TDatabase::ParamDB->SC_STEP_LENGTH_CONTROL_FINE_SCALAR))
          n_aux=4;
         else
          n_aux=2;
          
        if (mg_type==1)
        {    
#ifdef _MPI  
         MGLevel_low = new TMGLevel3D(i, sqmatrixA_low, rhs_low, sol_low, FEFunctArray[i], ParComm[i], OwnScalar_Spaces[i], n_aux, NULL);
#else
         MGLevel_low = new TMGLevel3D(i, sqmatrixA_low, rhs_low, sol_low, n_aux, NULL);
#endif 
         if (i==0)
          MG->AddLevel(MGLevel_low); 
         else
          MG->ReplaceLevel(i, MGLevel_low);  
        }
        
        // high order disc
        if ((i>=FirstSolve)||(mg_type==0))
        {  
#ifdef _MPI  
         MGLevel = new TMGLevel3D(mg_level, sqmatrixA, rhs, sol, FEFunctArray[mg_level], ParComm[mg_level], OwnScalar_Spaces[mg_level], n_aux, NULL);
#else
         MGLevel = new TMGLevel3D(mg_level, sqmatrixA, rhs, sol, n_aux, NULL);
#endif
         MG->AddLevel(MGLevel);
        }
        break;
    }  // switch(TDatabase::ParamDB->SOLVER_TYPE)
    

    // restrict solution to all grids
//     if (TDatabase::ParamDB->SOLVER_TYPE==GMG)
//       MG->RestrictToAllGrids();
    // if no solution on this grid, continue
    if(FirstSolve>i)
      continue;   
    
//======================================================================
// assemble all matrices
//======================================================================
      
     // build the discretizations
    for(k=low;k<=mg_level;k++)
    {   
     N_DOF = N_DOFs[k];      
     N_Active = Scalar_Spaces[k]->GetActiveBound();
     N_NonActive = N_DOF - N_Active;
 
     rhs = RhsArray[k];  
     RHSs[0] = rhs;   
    
     memset(rhs, 0, N_DOF*SizeOfDouble);
 
     fesp[0] = Scalar_Spaces[k];;
     ferhs[0] = Scalar_Spaces[k];;
 
     // find discrete form
//      if ((mg_type==1) && (k<i+1))
//       {
//        DiscreteForm = DiscreteFormUpwind;
//        OutPut("UPWIND"<<endl);
//       }
//      else
      switch(TDatabase::ParamDB->DISCTYPE)
      {
       case GALERKIN:
           DiscreteForm = DiscreteFormGalerkin;
       break;

       case SDFEM:
           DiscreteForm = DiscreteFormSDFEM;
       break;

       case GLS:
           DiscreteForm = DiscreteFormGLS;
       break;

       default:
            OutPut("Unknown DISCTYPE" << endl);
            exit(4711);;
      }

      // initialize matrices
     SQMATRICES[0] = MatricesA[k];
     SQMATRICES[0]->Reset();
     aux =  new TAuxParam3D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);
    
      // assemble
      Assemble3D(1, fesp,
               1, SQMATRICES,
               0, NULL,
               1, RHSs, ferhs,
               DiscreteForm,
               BoundaryConditions,
               BoundValues,
               aux);

#ifdef _MPI  
       ParComm[k]->SetSlaveDofRows(SQMATRICES[0]->GetRowPtr(), SQMATRICES[0]->GetKCol(), SQMATRICES[0]->GetEntries(), rhs);     
#endif
      
      
#ifdef _OPTTOMOGRAPHY_          
   RobinInt(SQMATRICES[0], rhs, BoundCondition);  
#endif      
  

     delete aux;
    } // for(k=low;k<=mg_level;k++)   // endfor, assembling done
 
    // set rhs for Dirichlet nodes
    memcpy(sol+N_Active, rhs+N_Active, N_NonActive*SizeOfDouble);    
    
 
//======================================================================
// solve the system
//======================================================================
    N_Unknowns = N_DOF;  
    PostSol = new double[N_Unknowns];
    memset(PostSol, 0, N_Unknowns*SizeOfDouble);     
    oldsol = new double[N_Unknowns]; 
    memset(oldsol, 0, N_Unknowns*SizeOfDouble);     
    defect = new double[N_Unknowns]; 
    memset(defect,0,N_Unknowns*SizeOfDouble);

//     OutPut("norm of solution " <<  sqrt(Ddot(N_Unknowns,sol,sol))  << endl);
//     OutPut("norm of rhs " <<  sqrt(Ddot(N_Unknowns,rhs,rhs))  << endl);
    
//     exit(0);    
    // solve system
    switch(TDatabase::ParamDB->SOLVER_TYPE)
    {  
     case DIRECT:     
      t1 = GetTime();
      DirectSolver(MatricesA[mg_level], rhs, sol);
      t2 = GetTime();
      OutPut( "time for solving: " << t2-t1 << endl);
//     OutPut("solution " << sqrt(Ddot(N_Unknowns,rhs,sol)) << endl);
     break;

     case AMG:
//         t1 = GetTime();
//         Solver(MatricesA[mg_level], RhsArray[mg_level], sol);
//         t2 = GetTime();
//         OutPut( "time for AMG solving: " << t2-t1 << endl);
//         OutPut("solution " << sqrt(Ddot(N_Unknowns,sol,sol)) << endl);
      break;

      case GMG:
      
        // build preconditioner
        switch (TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR)
        {
	 case 1:
	   prec = NULL;
	   break;
         case 5:
           prec = new TMultiGridScaIte(MatVect_Scalar, Defect_Scalar, NULL, 0, N_Unknowns, MG, zerostart);
         break;
         default:
          OutPut("Unknown preconditioner !!!" << endl);
          exit(4711);
        } //switch (TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR)
        
       
       switch (TDatabase::ParamDB->SC_SOLVER_SCALAR)
       {
        // fixed point iteration
        case 11:
         itmethod = new TFixedPointIte(MatVect_Scalar, Defect_Scalar, prec, 0, N_Unknowns, 1
#ifdef _MPI   
                               , ParComm[mg_level]
#endif                        
                                                );
         
         if (TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR == 5)
         {
          itmethod_sol = new double[N_Unknowns];
          itmethod_rhs = new double[N_Unknowns];
          memcpy(itmethod_sol, sol, N_Unknowns*SizeOfDouble);
          memcpy(itmethod_rhs, rhs, N_Unknowns*SizeOfDouble); 
         }
         else
         {
          itmethod_sol = sol;
          itmethod_rhs = rhs;
         }
        break;

        case 16:
         // FGMRES
         itmethod = new TFgmresIte(MatVect_Scalar, Defect_Scalar, prec, 0, N_Unknowns, 1
#ifdef _MPI   
                               , ParComm[mg_level]
#endif
	);
 
         if (TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR == 5)
         {
          itmethod_sol = new double[N_Unknowns];
          itmethod_rhs = new double[N_Unknowns];
          memcpy(itmethod_sol, sol, N_Unknowns*SizeOfDouble);
          memcpy(itmethod_rhs, rhs, N_Unknowns*SizeOfDouble);
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
       }  //  switch (TDatabase::ParamDB->SC_SOLVER_SCALAR)   
#ifdef _MPI
       t1 = MPI_Wtime();
#else
       t1 = GetTime();
#endif
       // solve linear system
       itmethod->Iterate(sqmatrices, NULL, itmethod_sol, itmethod_rhs);
#ifdef _MPI
       t2 = MPI_Wtime();
#else
       t2 = GetTime();
#endif
       switch (TDatabase::ParamDB->SC_SOLVER_SCALAR)
        {
          case 11:
          case 16:
            if (TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR == 5)
            {
              memcpy(sol, itmethod_sol, N_Unknowns*SizeOfDouble);
              memcpy(rhs, itmethod_rhs, N_Unknowns*SizeOfDouble);
            }
            break;
        } 
#ifdef _MPI
        if(rank==0) 
#endif
	  OutPut( "time for GMG solving: " << t2-t1 << endl);
        break;
    } //  switch(TDatabase::ParamDB->SOLVER_TYPE)
    
#ifdef _MPI
        if(rank==0) 
#endif
	{
         OutPut("MEMORY: " << setw(10) << GetMemory()/(1048576.0));
         OutPut(" MB" << endl); 
	}
  
#ifdef _MPI
     t_par1 = MPI_Wtime();

     if(TDatabase::ParamDB->WRITE_VTK)
      Output->Write_ParVTK(MPI_COMM_WORLD, img, SubID);
        img++;       
     t_par2 = MPI_Wtime();

     if(rank==out_rank)
     printf("Time taken for writing the parvtk file %e\n", (t_par2-t_par1));
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
      
//       os.seekp(std::ios::beg);
//        if(img<10) os <<  "VTK/"<<VtkBaseName<<".0000"<<img<<".plt" << ends;
//          else if(img<100) os <<  "VTK/"<<VtkBaseName<<".000"<<img<<".plt" << ends;
//           else if(img<1000) os <<  "VTK/"<<VtkBaseName<<".00"<<img<<".plt" << ends;
//            else if(img<10000) os <<  "VTK/"<<VtkBaseName<<".0"<<img<<".plt" << ends;
//             else  os <<  "VTK/"<<VtkBaseName<<"."<<img<<".plt" << ends;
//       
//       Output->WriteBinaryPlt(os.str().c_str());

     }
       img++;    
#endif    
  /*     
    // measure errors to known solution
    if(TDatabase::ParamDB->MEASURE_ERRORS)
     {
      aux =  new TAuxParam3D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);

      C->GetErrors(Exact, 3, AllDerivatives, 2, L2H1Errors,
                   BilinearCoeffs, aux, 1, fesp, errors);

      delete aux;

      l2[i] = errors[0];
      h1[i] = errors[1];
      sd[i] = errors[2];

      if (i)
       {
        OutPut( "L2: " << errors[0]<< " order " <<  log(l2[i-1]/l2[i])/ln2 << endl);
        OutPut( "H1-semi: " << errors[1] << " order " << log(h1[i-1]/h1[i])/ln2 << endl);
        OutPut( "SD: " << errors[2] << endl);
       }
      else
       {
        OutPut( "L2: " << errors[0] << endl);
        OutPut( "H1-semi: " << errors[1] << endl);
        OutPut( "SD: " << errors[2] << endl);
       }


     } // if(TDatabase::ParamDB->MEASURE_ERRORS)
 */

#ifdef _MPI
   if(rank==0)
#endif
     OutPut("memory after: " << setw(10) << GetMemory() << endl);
    
   } // for(i=0;i<LEVELS;i++)   
   // printf("time for collect into sendbuf is %lf----\n",timeC ); 
   

   
#ifdef _MPI    
   MPI_Finalize();  
#endif           
  CloseFiles();

  return 0;
}


 