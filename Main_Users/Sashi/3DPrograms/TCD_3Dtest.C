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


#include <TCD3D.h>
#include <TetraAffin.h>
#include <TetraIsoparametric.h>
#include <HexaAffin.h>
#include <HexaTrilinear.h>
#include <HexaIsoparametric.h>

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
#include <TimeUtilities.h>
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
// #include "../Examples/CD_3D/Laplace.h"
#include "../Examples/TCD_3D/Sin4.h"
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
 //#include "../Examples/TCD_3D/TeraHertzBrain.h"
// #include "../Examples/TCD_3D/TeraHertzCube.h"
// #include "../Examples/CD_3D/Sine.h"
// #include "../Examples/TCD_3D/TeraHertzBreast.h"



// ======================================================================
// utilities for main program
// ======================================================================

int main(int argc, char* argv[])
{
  TDomain *Domain = new TDomain();
  TDatabase *Database = new TDatabase();
  double t1, t2,time1,t_par1, t_par2, start_time;
  #ifdef _MPI
  const int root = 0;
  int rank, size;
  int MaxCpV, MaxSubDomainPerDof;
  int  out_rank, N_Cells_loc;
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
  
  TFEDatabase3D *FEDatabase = new TFEDatabase3D();
  TCollection *coll,*own_coll, *mortarcoll = NULL;
  TBaseCell *cell;
 TFESpace3D *scalar_space, *scalar_space_low, *RefinedOutput_space[2], **Scalar_Spaces; 
  
  TOutput3D *Output;
  int i,j,k,l,m, N_Unknowns, N_Active, N_DOF, N_DOF_low,N_NonActive, *N_DOFs,img=0; 
  double *B, *rhs, *sol, *oldsol, tol, tolerance, *defect, *startsol, *frac_step_sol;
  double *app, *oldrhs, *itmethod_sol, *itmethod_rhs, *current_sol, *current_B;
  double *sol_velo, *oldsol_nlinite, damp_nlinite;
  int n, N_, Len, low, only_first_time;
  int N_Rows, N_Columns, N_Unknowns_Velo, N_P,sold_parameter_type;
  double  *lump_mass, *matrix_D_Entries, *oldrhs_fem_fct0, *tilde_u, *oldrhs_fem_fct1;
  double *l2u1, *l2u2, *h1u1, *h1u2;
  double *l2p, *h1p, *sd, *l_inf, vect[3], exactvect[3];
  int which, *permutation;
  double DiffL2, DiffH1, t;
  int LEVELS, BASELEVEL, ORDER, order;
  int ret, pde;
  double negPower;
  double x,y,z,max,min,sum;
  double RE_NR;
  double tau1, tau2;
  double  p1, p2;
  double res, res2, oldres, solver_time, solver_time_curr, residual, oldresidual;
  double total_time, t3, t4, ass_time = 0.0;
  double impuls_residual,limit,linredfac;
  int N_LinIter, N_LinIterCurr, N_LinIterCurrIte, N_SubSteps, n_aux;
  double gamma, tau, oldtau;
  int *RowPtr;

  std::ostringstream os;
 
  TFEFunction3D *C,*C_low, **velo1, **velo2, **velo3, *fefct[3],**FEFunctArray;
  TFEFunction3D **SolArray, **AuxFEFunctArray, *pressure;
  TFEFunction3D *oldconc, **OldSolArray;
  TFEVectFunct3D **velocity, **AuxFEVectFunctArray;
  TFESpace3D *fesp[2], *ferhs[1];
  double delta, end_time;

  TAuxParam3D *aux;
  TSquareStructure3D *sqstructureA,*sqstructureA_low;
  TSquareMatrix3D *sqmatrixA,*sqmatrixA_low, *SQMATRICES[3];
  TSquareMatrix3D *sqmatrixM ,*sqmatrixM_low;
  TSquareMatrix3D **MatricesA, **MatricesM, **MatricesK, **MatricesS;
  TSquareMatrix **sqmatrices = (TSquareMatrix **)SQMATRICES;
  MatVecProc *MatVect;
  DefectProc *Defect;

  TMGLevel3D *MGLevel,*MGLevel_low;
  TMultiGrid3D *MG;

  double *rhs_low, *sol_low,  *RHSs[1], errors[4], **RhsArray;
  int *N_Uarray;

  // discrete forms have to be created
  TDiscreteForm3D *DiscreteForm,*DiscreteFormGalerkin;

  TDiscreteForm3D *DiscreteFormMatrixMRhs;
  TDiscreteForm3D *DiscreteFormMatrixARhs;
  TDiscreteForm3D *DiscreteFormRhs;


  int N_SquareMatrices, N_Rhs, N_FESpaces;
  
  BoundCondFunct3D *BoundaryConditions[1];
  BoundValueFunct3D *BoundValues[1];
  CoeffFct3D *Coefficients[1];

  TItMethod *itmethod, *prec;
  int Max_It, FirstSolve;
  double omega, alpha, alpha_fine, divergence;
  int N_Paramters=1, methods, time_discs;
  double Parameters[2], hmin, hmax;

  int N_UConv, level_down, ii, fully_implicit = 0;
  int mg_level,mg_type,CurrentDiscType, last_sq, step_length;
  int very_first_time=0, zerostart, comp_vort, picture;
  int first_matrix_assemble =1 ;
  int *neum_to_diri, N_neum_to_diri = 0;
  double *neum_to_diri_x, *neum_to_diri_y, *neum_to_diri_z;
  double *l2, *h1;
  char *PRM, *GEO;
  char *PsBaseName, *VtkBaseName;
  char CString[] = "c";
  char UString[] = "u";
  char RUString[] = "ref_u";
  char PString[] = "p";
  char ReadinDat [] = "readin.dat";
  char MMString[] = "Mass matrix";
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

#ifdef __FUEL_CELL__
  int N_OutFlowCells, *OutFlowFaces, *OutFlowCells;
  double outflow, area;
#endif // __FUEL_CELL__

#ifdef __STUTZEN__
  TBaseCell **Cells;
  int N_CellsOld, N_CellsNew;
  int N_Vertices;
  double xm, ym, zm;
  double xp, yp, zp;
  double r1, r2;
#endif

  os << " ";

  total_time = GetTime();
//======================================================================
// read parameter file
//======================================================================
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

#ifdef __STUTZEN__
  TDatabase::ParamDB->INTERNAL_PROBLEM_IDENTITY = 1234;
#endif

  Database->WriteParamDB(argv[0]);
  Database->WriteTimeDB();
#ifdef _MPI
  if(rank==0)
#endif
       ExampleFile();
//======================================================================
// copy read parameters into local variables
//======================================================================
   Coefficients[0] = BilinearCoeffs;

  PRM = TDatabase::ParamDB->BNDFILE;
  GEO = TDatabase::ParamDB->GEOFILE;

  VtkBaseName = TDatabase::ParamDB->VTKBASENAME;

  FirstSolve = TDatabase::ParamDB->SC_FIRST_SOLUTION_LEVEL_SCALAR;
  // set type of multilevel
  mg_type = TDatabase::ParamDB->SC_MG_TYPE_SCALAR;
 
  if(TDatabase::ParamDB->SOLVER_TYPE==AMG|| TDatabase::ParamDB->SOLVER_TYPE==DIRECT)
       mg_type=0;
  
  if(mg_type)
    mg_level = 0;
  else
    mg_level = -1;

  LEVELS = TDatabase::ParamDB->LEVELS;
  BASELEVEL = TDatabase::ParamDB->UNIFORM_STEPS;
  l2u1 = new double[LEVELS+1];
  l2u2 = new double[LEVELS+1];
  h1u1 = new double[LEVELS+1];
  h1u2 = new double[LEVELS+1];
  sd = new double[LEVELS+1];
  l_inf = new double[LEVELS+1];

  SolArray = new TFEFunction3D*[LEVELS+1];
  OldSolArray = new TFEFunction3D*[LEVELS+1];
  RhsArray = new double* [LEVELS+1];
  Scalar_Spaces = new TFESpace3D*[LEVELS+1];
  FEFunctArray = new TFEFunction3D*[LEVELS+1];
  MatricesA = new TSquareMatrix3D*[LEVELS+1];
  MatricesM = new TSquareMatrix3D*[LEVELS+1];
  N_DOFs = new int[LEVELS+1];  
  l2 = new double[LEVELS+1];
  h1 = new double[LEVELS+1];
  sd = new double[LEVELS+1];
 
  
#ifdef _MPI    
  OwnScalar_Spaces = new TFESpace3D*[LEVELS+1];   
  ParComm = new TParFECommunicator3D*[LEVELS+1];  
  
  out_rank=TDatabase::ParamDB->Par_P0;  
#endif   
  MatVect = MatVect_Scalar;
  Defect = Defect_Scalar;

//======================================================================
// initialize discrete forms
//======================================================================
  InitializeDiscreteForms(DiscreteFormGalerkin, BilinearCoeffs);  
// discrete form for assembling mass matrix
   DiscreteFormMatrixMRhs = new TDiscreteForm3D
    (MMString, MMString, N_Terms_MatrixMRhs,
     Derivatives_MatrixMRhs,
     SpacesNumbers_MatrixMRhs, N_Matrices_MatrixMRhs, N_Rhs_MatrixMRhs,
     RowSpace_MatrixMRhs, ColumnSpace_MatrixMRhs, RhsSpace_MatrixMRhs,
     MatrixMRhsAssemble, BilinearCoeffs, NULL);

    DiscreteFormMatrixARhs = new TDiscreteForm3D
    (MMString, MMString, N_Terms_MatrixARhs,
     Derivatives_MatrixARhs,
     SpacesNumbers_MatrixARhs, N_Matrices_MatrixARhs,
     N_Rhs_MatrixARhs,
     RowSpace_MatrixARhs, ColumnSpace_MatrixARhs,
     RhsSpace_MatrixARhs,
     MatrixARhsAssemble, BilinearCoeffs, NULL);

    DiscreteFormRhs = new TDiscreteForm3D
    (MMString, MMString, N_Terms_Rhs, Derivatives_Rhs,
     SpacesNumbers_Rhs, N_Matrices_Rhs, N_Rhs_Rhs,
     RowSpace_Rhs, ColumnSpace_Rhs, RhsSpace_Rhs,
     RhsAssemble, BilinearCoeffs, NULL);



//======================================================================
// read boundary parameterization and initialize coarse grid
//======================================================================
 Domain->Init(PRM, GEO);
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
//   
 /* for(i=3;i<TDatabase::ParamDB->UNIFORM_STEPS;i++)
  {
    printf("************************************LEVEL %d****************************************",i);
    Domain->RegRefineAll();
    Domain->GenerateEdgeInfo();
    Domain_Crop(Comm, Domain);
  }*/
  if(rank==0)
    printf("Time taken for Domain Decomposition is %e\n", (t2-t1));
  MaxSubDomainPerDof = MIN(MaxCpV, size);

#else
  
  for(i=0;i<TDatabase::ParamDB->UNIFORM_STEPS;i++)
    Domain->RegRefineAll();  
#endif 
  ORDER = TDatabase::ParamDB->ANSATZ_ORDER;
  limit = TDatabase::ParamDB->SC_NONLIN_RES_NORM_MIN_SCALAR;
  Max_It = TDatabase::ParamDB->SC_NONLIN_MAXIT_SCALAR;
  omega = TDatabase::ParamDB->SC_NONLIN_DAMP_FACTOR_SCALAR;
  alpha = TDatabase::ParamDB->SC_GMG_DAMP_FACTOR_SCALAR;
  alpha_fine = TDatabase::ParamDB->SC_GMG_DAMP_FACTOR_FINE_SCALAR;
  divergence = TDatabase::ParamDB->SC_DIV_FACTOR;
  
  Parameters[0] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_SCALAR;
  Parameters[1] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_FINE_SCALAR;

  TDatabase::TimeDB->CURRENTTIME = TDatabase::TimeDB->STARTTIME;
  SetTimeDiscParameters();

  if(TDatabase::ParamDB->SOLVER_TYPE == GMG)
  {
    i=1;
    MG = new TMultiGrid3D(i, N_Paramters, Parameters);
  }
 

  damp_nlinite = TDatabase::ParamDB->SC_NONLIN_DAMP_FACTOR_SCALAR;

  t3 = GetTime();
  total_time = t3 - total_time;

  
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
 
  TDatabase::TimeDB->CURRENTTIME = TDatabase::TimeDB->STARTTIME;
  SetTimeDiscParameters();

  t3 = GetTime();
  total_time = t3 - total_time; 
 
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
    } //*/

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
    os.seekp(std::ios::beg);
    Output->AddParameter(TDatabase::TimeDB->CURRENTTIME,os.str().c_str());
    
//   / if(i==FirstSolve)    
//     {
//       C->Interpolate(Exact);
// #ifdef _MPI
//       ParSolVect->ParDdot(BYOWN, norm);
// #endif
//         if(rank==out_rank)
//          OutPut("Exact Norm of Par sol " << sqrt( norm[0]) <<endl); 
//     }
#ifdef _MPI
     t_par1 = MPI_Wtime();

     if(TDatabase::ParamDB->WRITE_VTK)
      Output->Write_ParVTK(Comm, img, SubID);
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
      sqmatrixM_low = new TSquareMatrix3D(sqstructureA_low);
      MatricesM[i] = sqmatrixM_low;
    } 

    // build matrices for high order disc
    if ((i>=FirstSolve)||(mg_type==0))
    {
     sqstructureA = new TSquareStructure3D(scalar_space);
     sqstructureA->Sort();

     sqmatrixA = new TSquareMatrix3D(sqstructureA);
     MatricesA[mg_level] = sqmatrixA;
     sqmatrixM = new TSquareMatrix3D(sqstructureA);
     MatricesM[mg_level] = sqmatrixM;
    } 
    
   B = new double [N_Unknowns];
   current_B = B;
      
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
         MGLevel_low = new TMGLevel3D(i, sqmatrixM_low, rhs_low, sol_low, FEFunctArray[i], ParComm[i], OwnScalar_Spaces[i], n_aux, NULL);
#else
         MGLevel_low = new TMGLevel3D(i, sqmatrixM_low, rhs_low, sol_low, n_aux, NULL);
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
         MGLevel = new TMGLevel3D(mg_level, sqmatrixM, current_B, sol, FEFunctArray[mg_level], ParComm[mg_level], OwnScalar_Spaces[mg_level], n_aux, NULL);
#else
         MGLevel = new TMGLevel3D(mg_level, sqmatrixM, current_B, sol, n_aux, NULL);
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
// all data are available for assembling matrices
// loop over all levels
//======================================================================
  for(k=0;k<=mg_level;k++)
  {
    // set parameters
    N_Rhs = 1;
    N_FESpaces = 1;
    fesp[0] = Scalar_Spaces[k];
    aux =  new TAuxParam3D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);

    //======================================================================
    // assembling of mass matrix and rhs
    //======================================================================
    // reset matrices
    N_SquareMatrices = 1;
    SQMATRICES[0] = MatricesM[k];
    SQMATRICES[0]->Reset();

    DiscreteForm = DiscreteFormMatrixMRhs; 
// 
    BoundaryConditions[0] =  BoundCondition;
    BoundValues[0] = BoundValue;

    memset(RhsArray[k], 0, N_DOFs[k]*SizeOfDouble);  
    RHSs[0] = RhsArray[k];
    ferhs[0] = Scalar_Spaces[k];

    Assemble3D(N_FESpaces, fesp, 
               N_SquareMatrices, SQMATRICES, 
               0, NULL, 
               N_Rhs, RHSs, ferhs,
               DiscreteForm, 
               BoundaryConditions, 
               BoundValues, 
               aux);
    delete aux;    
    
#ifdef _MPI  
   // ParComm[k]->SetSlaveDofRows(SQMATRICES[0]->GetRowPtr(), SQMATRICES[0]->GetKCol(), SQMATRICES[0]->GetEntries(), rhs);     
#endif
   }   // endfor k
   
//     // save solution
 

  // allocate arrays for solver
//   oldsol = new double[N_Unknowns];
  defect = new double[N_Unknowns];
  startsol = new double[N_Unknowns];
  frac_step_sol = new double[N_Unknowns];
  oldrhs =  new double[N_Unknowns];

  // memcpy(oldsol,sol,N_Unknowns*SizeOfDouble);
  // number of active d.o.f.
  N_Active = Scalar_Spaces[mg_level]->GetActiveBound();
  N_Unknowns = Scalar_Spaces[mg_level]->GetN_DegreesOfFreedom();

  solver_time = 0.0;
  N_LinIter = 0;
  
  // parameters for time stepping scheme
  gamma = 0;
  m = 0;
  N_SubSteps = GetN_SubSteps();

  oldtau = 1;
  end_time = TDatabase::TimeDB->ENDTIME;

  // not active : TDatabase::TimeDB->TIMESTEPLENGTH_CONTROL = 0 
  if (TDatabase::TimeDB->TIMESTEPLENGTH_CONTROL)
    time_discs = 2;
  else
    time_discs = 1; 

  // initialize solver
  if (TDatabase::ParamDB->SOLVER_TYPE==GMG)
  {
    switch (TDatabase::ParamDB->SC_SOLVER_SCALAR)
    {
      case 11:
        zerostart = 1;
        break;
      case 16:
        zerostart = 0;
        break;
    }
    switch (TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR)
    {
      case 5:
        prec = new TMultiGridScaIte(MatVect, Defect, NULL,
                                    0, N_Unknowns, MG, zerostart);
        break;
      default:
        OutPut("Unknown preconditioner !!!" << endl);
       exit(4711);
    }
    switch (TDatabase::ParamDB->SC_SOLVER_SCALAR)
    {
       // fixed point iteration 
      case 11:
        itmethod = new TFixedPointIte(MatVect, Defect, prec,0, N_Unknowns, 1
#ifdef _MPI   
                               , ParComm[mg_level]
#endif 
                                                  );
        if (TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR == 5)
        {
          itmethod_sol = new double[N_Unknowns];
          itmethod_rhs = new double[N_Unknowns];              
          }  
        else
        {
          itmethod_sol = sol;
          itmethod_rhs = rhs;            
        }
        break;
      case 16:
         // FGMRES
        itmethod = new TFgmresIte(MatVect, Defect, prec,0, N_Unknowns, 1
#ifdef _MPI   
                               , ParComm[mg_level]
#endif
                                                  );
        if (TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR == 5)
        {
          itmethod_sol = new double[N_Unknowns];
          itmethod_rhs = new double[N_Unknowns];              
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
    }
  }

  coll->GetHminHmax(&hmin,&hmax);
#ifdef _MPI
  if(rank==0)
#endif
  OutPut("h_min : " << hmin << " h_max : " << hmax << endl);
  B = new double[N_Unknowns];
  //======================================================================
// start of time cycle
// everything happens on the same grid
//======================================================================
 img=10;
  while(TDatabase::TimeDB->CURRENTTIME< end_time)
  {                              // time cycle
    
// #ifdef _MPI    
//    start_steptime = MPI_Wtime();
// #else  
//    start_steptime = GetTime();
// #endif
    m++;
    TDatabase::TimeDB->INTERNAL_STARTTIME = TDatabase::TimeDB->CURRENTTIME;
    
    for (methods=0;methods<time_discs;methods++)
    {
      if (time_discs==2)
      {
        if (methods==0) // fractional-step-theta-scheme
        {
          TDatabase::TimeDB->TIME_DISC = 3;
          memcpy(startsol,sol,N_Unknowns*SizeOfDouble); // save start sol
          memcpy(oldrhs,rhs,N_Unknowns*SizeOfDouble); // save rhs
        }
        else           // crank nicolson scheme
        {              // take solution of first scheme as initial iterate
          TDatabase::TimeDB->TIME_DISC = 2;
          TDatabase::TimeDB->CURRENTTIME = TDatabase::TimeDB->INTERNAL_STARTTIME;
          memcpy(frac_step_sol,sol,N_Unknowns*SizeOfDouble); // save solution of fract.step
          memcpy(sol,startsol,N_Unknowns*SizeOfDouble); // get former startsol
          memcpy(rhs,oldrhs,N_Unknowns*SizeOfDouble); // get old rhs
        }
        N_SubSteps = GetN_SubSteps();
	 
      }

      for(l=0;l<N_SubSteps;l++)      // sub steps of fractional step theta
      {
        if (!very_first_time)
        {
          SetTimeDiscParameters();
        }
        if (m==1
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
         j=0; 
#ifdef _MPI
     if(rank==out_rank)
#endif
     {
        OutPut(endl << "CURRENT TIME: ");
        OutPut(TDatabase::TimeDB->CURRENTTIME << endl);
     }
	solver_time_curr = 0;
   
       N_LinIterCurr = 0;
       // start iteration for solving nonlinear problems
     
        // working array for rhs is B, initialize B
        memset(B, 0, N_Unknowns*SizeOfDouble);
        // old rhs multiplied with current subtime step and theta3 on B
        Daxpy(N_Active, tau*TDatabase::TimeDB->THETA3, rhs, B);
	 
        // assembling of A and rhs
#ifdef _MPI    
    t1 = MPI_Wtime();
#else  
    t1 = GetTime();
#endif
          for(k=0;k<=mg_level;k++)
          {  
	  
            // set parameters
            N_Rhs = 1;
	    N_FESpaces = 1;
	    fesp[0] = Scalar_Spaces[k];
            aux =  new TAuxParam3D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);
  

            //======================================================================
            // assembling of mass matrix and rhs
            //======================================================================

                N_SquareMatrices = 1;
                SQMATRICES[0] = MatricesA[k];
                SQMATRICES[0]->Reset();
                DiscreteForm = DiscreteFormMatrixARhs; 
   
            BoundaryConditions[0] =  BoundCondition;
            BoundValues[0] = BoundValue;
            memset(RhsArray[k], 0, N_DOFs[k]*SizeOfDouble);  
            RHSs[0] = RhsArray[k];
            ferhs[0] = Scalar_Spaces[k];
            
            Assemble3D(N_FESpaces, fesp, 
                       N_SquareMatrices, SQMATRICES, 
                       0, NULL, 
                       N_Rhs, RHSs, ferhs,
                       DiscreteForm, 
                       BoundaryConditions, 
                       BoundValues, 
                       aux);

            delete aux;   


	  } // endfor i
	  
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
	  
          first_matrix_assemble = 0;
	
#ifdef _MPI      
      if(rank==out_rank)      
#endif 
	  OutPut("Norm of rhs " <<  sqrt(Ddot(N_DOF,RhsArray[mg_level],RhsArray[mg_level])) <<endl);
        //======================================================================
        // only rhs needs to be assembles  
        // only on finest level  
        //======================================================================
      
        if (very_first_time==1)
        {
           very_first_time=0;
           l--;
           continue;
        }
	// add rhs from current sub time step to rhs array B
        Daxpy(N_Active, tau*TDatabase::TimeDB->THETA4, RhsArray[mg_level], B);
          
        //======================================================================
        // manipulation of matrices due to current time discretization
        // the stiffness matrix is stored on M 
        //======================================================================
        oldtau = tau;
     	    // update rhs by Laplacian and convective term from previous time step
          // scaled by current sub time step length and theta2
          // currently : M := M + gamma A
          // M = M + (-gamma - tau*TDatabase::TimeDB->THETA2)
          for(k=0;k<=mg_level;k++)
          {
            MatAdd(MatricesM[k], MatricesA[k],-gamma - tau*TDatabase::TimeDB->THETA2);
          }
          // set current factor of steady state matrix
          gamma = -tau*TDatabase::TimeDB->THETA2;

          // defect = M * sol
          // B:= B + defect
          MatVectActive(MatricesM[mg_level], sol, defect);
          Daxpy(N_Active, 1, defect, B);
          // set Dirichlet values
          // RHSs[0] still available from assembling

        memcpy(B+N_Active, RHSs[0]+N_Active, (N_Unknowns-N_Active)*SizeOfDouble);
        // copy Dirichlet values from rhs into sol
        memcpy(sol+N_Active, RHSs[0]+N_Active, (N_Unknowns-N_Active)*SizeOfDouble);

          // M = M + (-gamma + tau*TDatabase::TimeDB->THETA1) A
          for(k=0;k<=mg_level;k++)
          {
            MatAdd(MatricesM[k], MatricesA[k],-gamma + tau*TDatabase::TimeDB->THETA1);
	    SQMATRICES[0] = MatricesM[k];
#ifdef _MPI  
  ParComm[k]->SetSlaveDofRows(SQMATRICES[0]->GetRowPtr(),SQMATRICES[0]->GetKCol(),SQMATRICES[0]->GetEntries(), RhsArray[k]);     
#endif
          }
      
          // set current factor of steady state matrix
          gamma = tau*TDatabase::TimeDB->THETA1;
    
        //======================================================================
        // solution of linear system
        //======================================================================
        
        memset(defect, 0, N_Unknowns*SizeOfDouble);
        SQMATRICES[0] = MatricesM[mg_level];
           
// 	if (!(TDatabase::ParamDB->INTERNAL_LINEAR_SCHEME))
// 	    memcpy(oldsol_nlinite, sol, N_Unknowns*SizeOfDouble);
        //======================================================================
        // solve linear system
        //======================================================================
        switch(TDatabase::ParamDB->SOLVER_TYPE)
        {
	  case DIRECT:     
             t1 = GetTime();
             DirectSolver(MatricesM[mg_level], B, sol);
             t2 = GetTime();
	     break;
          case AMG:
	    TDatabase::ParamDB->SC_VERBOSE=0;     
            t1 = GetTime();
            Solver(SQMATRICES[0], B, sol);
            t2 = GetTime();
            solver_time_curr = t2-t1;
            solver_time += t2-t1;
          break;
             
          case GMG:
           
            if (TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR == 5)
            {
               memcpy(itmethod_sol, sol, N_Unknowns*SizeOfDouble);
               memcpy(itmethod_rhs, B, N_Unknowns*SizeOfDouble);
            }
#ifdef _MPI
       t1 = MPI_Wtime();
#else
       t1 = GetTime();
#endif    
            N_LinIterCurrIte = itmethod->Iterate(sqmatrices,NULL,itmethod_sol,itmethod_rhs);
#ifdef _MPI
       t2 = MPI_Wtime();
#else
       t2 = GetTime();
#endif  
	    N_LinIterCurr += N_LinIterCurrIte;
	    N_LinIter += N_LinIterCurrIte ;
            if (TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR == 5)
            {
               memcpy(sol, itmethod_sol, N_Unknowns*SizeOfDouble);
               memcpy(B, itmethod_rhs, N_Unknowns*SizeOfDouble);
            }
         
            solver_time_curr += t2-t1;
            solver_time += t2-t1;
          break; 
        } // endswitch SOLVER_TYPE
//	OutPut("solution " << sqrt(Ddot(N_Unknowns,sol,sol)) << endl);

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
           OutPut("Time for Sover " << t2-t1 << "s "<<  endl);   
          }

        //======================================================================
        // end solve linear system
        //======================================================================

          // restore mass matrices by subtracting the A-matrices
          for(k=0;k<mg_level;k++)
          {
            MatAdd(MatricesM[k], MatricesA[k], -gamma);
          }
          // set current factor of steady state matrix
          gamma = 0;
#ifdef _MPI
        ParComm[mg_level]->CommUpdate(sol);
#endif 
    // save solution
     // memcpy(oldsol,sol,N_Unknowns*SizeOfDouble);
    }   // endfor l (sub steps of fractional step theta)                                      
     
      // measure errors to known solution
      if(TDatabase::ParamDB->MEASURE_ERRORS)
      {
        aux =  new TAuxParam3D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);
        C->GetErrors(Exact, 4, AllDerivatives, 2, L2H1Errors, 
                        BilinearCoeffs, aux, 1, fesp, errors);
        delete aux;
#ifdef _MPI     
      if(rank==0)
       {
#endif
        OutPut("time: " << TDatabase::TimeDB->CURRENTTIME);
        OutPut(" L2: " << errors[0]);
        OutPut(" H1-semi: " << errors[1] << endl);
#ifdef _MPI
       }
#endif
      } // endif MEASURE_ERRORS
  
    } // endfor two time discs of adaptive time step control

     
     if(TDatabase::ParamDB->WRITE_VTK &&
        (m % TDatabase::TimeDB->STEPS_PER_IMAGE == 0) )
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
      }

  } // while*/
    
 } //for(i=0;i<LEVELS;i++)*/

//======================================================================
// end of time cycle
//======================================================================
  if(TDatabase::ParamDB->WRITE_VTK)
  {
      os.seekp(std::ios::beg);
      os << VtkBaseName << "end." << m << ".vtk" << ends;
      Output->WriteVtk(os.str().c_str());
  }

  t4 =  GetTime();
  total_time += t4 - t3;
  OutPut("total running time: " << total_time << endl);
  CloseFiles();
  return 0;
}

