// =======================================================================
//
// Purpose:     main program with parallel solver (no multigrid solver)
//
// Author:      Sashikumaar Ganesan
//
// History:     Implementation started on 24.07.2009
//              OpenMP implementation - 24.07.2009
// =======================================================================
#ifdef _MPI
# include "mpi.h"
#endif

#ifdef _OMPONLY
#include <omp.h>
#endif

#include <Domain.h>
#include <Database.h>
#include <DiscreteForm2D.h>
#include <FEDatabase2D.h>
#include <FESpace2D.h>
#include <SquareStructure2D.h>
#include <Structure2D.h>
#include <AuxParam2D.h>
#include <QuadAffin.h>
#include <DirectSolver.h>
#include <Assemble2D.h>
#include <Output2D.h>
#include <LinAlg.h>
#include <CD2DErrorEstimator.h>

double bound = 0;

#include <string.h>
#include <sstream>
#include <MooNMD_Io.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>

#include <MainUtilities.h>
#include <Upwind.h>
#include <FixedPointIte.h>

#include <MacroCell.h>
#include <BoundEdge.h>
#include <IsoBoundEdge.h>
#include <IsoInterfaceJoint.h>
#include <BdLine.h>
#include <BdCircle.h>
#include <GridCell.h>
#include <LocalProjection.h>
#include <gridgen.h>

#include <QuadAffin.h>
#include <QuadBilinear.h>
#include <QuadIsoparametric.h>
#include <TriaAffin.h>
#include <TriaIsoparametric.h>

#include <sys/stat.h>
#include <sys/types.h>


#ifdef _MPI
 #include <MeshPartition.h>
//  #include <ParFECommunicator2D.h>
//  #include <MumpsSolver.h>
// #include <ParVectorNSE.h>
#endif

// for external mesh generator
extern "C"
{
  void triangulate(char*, struct triangulateio*,
                   struct triangulateio*, struct triangulateio*);
}

// ======================================================================
// include the required example file
// ======================================================================
// #include "../Examples/ConvDiff2D/ReactionDominate.h" 




// #include "../Examples/ConvDiff2D/Smooth.h"  // unit square

// #include "../Examples/ConvDiff2D/SineLaplaceDiriHom.h" 
// #include "../Examples/ConvDiff2D/SineLaplace.h" 
// #include "../Examples/ConvDiff2D/Hemker1996.h" 
// #include "../Examples/ConvDiff2D/Smooth.h"
// #include "../Examples/ConvDiff2D/Circle02.h"
// #include "../Examples/ConvDiff2D/Circle03.h"
//    #include "../Examples/CD_2D/Terahertz.h"
//    #include "../Examples/CD_2D/Sine.h"
//    #include "../Examples/CD_2D/DiffOptTomography.h"

// #include "../Examples/CD_2D/Plane.h"  // unit square
// #include "../Examples/CD_2D/TwoBoundaryLayers.h"  // unit square
// #include "../Examples/CD_2D/TwoInteriorLayers.h"  // unit square
// #include "../Examples/CD_2D/TwoBoundaryLayersLaplace.h"  // unit square

#include "../Examples/CD_2D/Hemker1996.h" // circle in a channel

// ======================================================================
// utilities for main program
// ======================================================================

/** printing scalar at y=0 **/ 
void PrintScalar(TFEFunction2D *Fefunct, int &N_BData)
{
 int i, N;
 
 double dx, x, y, Val[3];
 
 char *VtkBaseName;
 
 VtkBaseName = TDatabase::ParamDB->VTKBASENAME;

 std::ostringstream os;
 os << " ";

 os.seekp(std::ios::beg);
//       os << "surfact"<< i << ".dat" << ends;
  if(N_BData<10) os << "BDData/"<<VtkBaseName<<"Gamma_0000"<<N_BData<<".data" << ends;
  else if(N_BData<100) os <<"BDData/"<<VtkBaseName<<"Gamma_000"<<N_BData<<".data" << ends;
  else if(N_BData<1000) os <<"BDData/"<<VtkBaseName<<"Gamma_00"<<N_BData<<".data" << ends;
  else if(N_BData<10000) os <<"BDData/"<<VtkBaseName<<"Gamma_0"<<N_BData<<".data" << ends;
  else  os <<"BDData/"<<VtkBaseName<<"Gamma_"<<N_BData<<".data" << ends;
  
 std::ofstream dat(os.str().c_str());

  if (!dat)
   {
    cerr << "cannot open file for output" << endl;
    exit(0);
   }
  dat << "%% Scalar data created by ParMooN" << endl;
//   dat << "%% Current Reference Time :" << TDatabase::TimeDB->CURRENTTIME << endl;
  dat << "%% x, y, value" << endl; 
  
  N = 50;
  
  y=0;  
  x=0.1;
  dx = -60./(double)N;
  
  for(i=0;i<N;i++)
   {
    Fefunct->FindGradient(x,  y, Val);
    dat << x<< " " << y << "  " << Val[0] <<endl;

//     x += dx;
    y += dx;   
   }  
   
  N_BData++;
}




int main(int argc, char* argv[])
{
  TDatabase *Database = new TDatabase();
  
#ifdef _MPI
  const int root = 0;
  int rank, size, len;
  double t_par1, t_par2;
  char  name[MPI_MAX_PROCESSOR_NAME];

  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Get_processor_name(name, &len);
#endif 
   
  TDomain *Domain = new TDomain();   
  TFEDatabase2D *FEDatabase = new TFEDatabase2D();
  TCollection *coll;
  TBaseCell *cell;
  TFESpace2D *scalar_space, *RefinedOutput_space[2];
  TOutput2D *Output;
  TAuxParam2D *aux;
  TSquareStructure2D *sqstructureA;
  TSquareMatrix2D *sqmatrixA, *SQMATRICES[1];
  TSquareMatrix2D **MatricesA, *Tmp_MatricesA;
  TSquareMatrix **sqmatrices = (TSquareMatrix **)SQMATRICES;
  TFEFunction2D *u1, *u2, *C_log;
  TFESpace2D  *velocity_space;
  TFEVectFunct2D *u;

  TDiscreteForm2D *DiscreteForm;
  TDiscreteForm2D *DiscreteFormUpwind;  
  TDiscreteForm2D *DiscreteFormGalerkin;
  TDiscreteForm2D *DiscreteFormSDFEM;
  TDiscreteForm2D *DiscreteFormGLS;
  TFEFunction2D *C, *Refinmed_C;
  TFESpace2D *fesp[2], *ferhs[3];
  MultiIndex2D AllDerivatives[3] = { D00, D10, D01 };

  int i,j,k,l,ret, ORDER, N_DOF, N_Active, N_NonActive;
  int N_Unknowns, img=1, N_RDOF[0];
  int N_LinIter, LEVELS, N_U_DOF, N_BData=0;

  double solver_time, hmin, hmax;
  double *sol, *oldsol, *rhs, *defect, *PostSol, *refined_sol[2], *auxsol;
  double *RHSs[1], t1, t2, errors[4];
  double *l2, *h1, *sd;
  double *u_sol, *sol_log;

  char *PRM, *GEO;
  char *PsBaseName, *VtkBaseName;
  char CString[] = "c";
  char UString[] = "u";
  char UlogString[] = "log_u";
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
    
  boolean ForPressure=FALSE, ReadVelo=FALSE;

  std::ostringstream os, opts;
  os << " ";
  opts << " ";
  
  mkdir(vtkdir, 0777);
  mkdir(BDdir, 0777);  
  

//======================================================================
// read parameter file Readin.dat
//======================================================================
  if(argc>=2)
    ret=Domain->ReadParam(argv[1]);
  else
    ret=Domain->ReadParam(ReadinDat);
 
  if(ret==-1)
    exit(-1);

  OpenFiles();
  OutFile.setf(std::ios::scientific);

  Database->WriteParamDB(argv[0]);
  ExampleFile();

#ifdef _MPI 
//   if(rank==2)
    printf("Test  Decomposition %d \n", rank);     
#endif    
  
//======================================================================
// copy read parameters into local variables
//======================================================================
  PRM = TDatabase::ParamDB->BNDFILE;
  GEO = TDatabase::ParamDB->GEOFILE;
  PsBaseName = TDatabase::ParamDB->PSBASENAME;
  VtkBaseName = TDatabase::ParamDB->VTKBASENAME;
  LEVELS = TDatabase::ParamDB->LEVELS;

  l2 = new double[LEVELS+1];
  h1 = new double[LEVELS+1];
  sd = new double[LEVELS+1];
//======================================================================
// define discrete form
//======================================================================


//   DiscreteFormGalerkin = new TDiscreteForm2D(CdString, GalString,
//                                   N_Terms, Derivatives, SpacesNumbers,
//                                   N_Matrices, N_Rhs, RowSpace,
//                                   ColumnSpace, RhsSpace, BilinearAssemble,
//                                   BilinearCoeffs, NULL);

  InitializeDiscreteForms_Stationary(DiscreteFormUpwind, DiscreteFormGalerkin, DiscreteFormSDFEM, DiscreteFormGLS,
                                     BilinearCoeffs);
  
//======================================================================
// generate mesh and refun (if needed)
//======================================================================
    Domain->Init(PRM, GEO);

//       write grid into an Postscript file
    os.seekp(std::ios::beg);
    os << "Domain_Coarse" << ".ps" << ends;
    Domain->PS(os.str().c_str(),It_Finest,0);


/* generate special mesh for Hemker example */
#ifdef __HMM_1986__
    if(!strcmp(GEO, "InitGrid"))
     if(TDatabase::ParamDB->REACTOR_P25)
         MeshReGen_HemkerResolved(Domain);
     else
         MeshReGen_Hemker(Domain);
#endif

  // refine grid
  for(i=0;i<TDatabase::ParamDB->UNIFORM_STEPS;i++)
	      Domain->RegRefineAll();
//        Domain->AdaptRefineAll();  
 
   // write grid into an Postscript file
    os.seekp(std::ios::beg);
    os << "Domain" << ".ps" << ends;
    Domain->PS(os.str().c_str(),It_Finest,0);

//  int  nsteps=7;
//      for(i=1;i<nsteps;i++)
//      {
//        Domain->AdaptRefineAll();  
//        //DomainA1 = Domain;
//        
//        // write grid into an Postscript file
//        os.seekp(std::ios::beg);
//        os << "DomainRef" << i << ".ps" << ends;
//        Domain->PS(os.str().c_str(),It_Finest,0);
//      } // */
//     
// 
//       exit(0); // exit here if u want to see only mesh
 
//======================================================================
// include the boundary condition and boundary values from the example file
//======================================================================
  BoundCondFunct2D *BoundaryConditions[1] = { BoundCondition };
  BoundValueFunct2D *BoundaryValues[1] = { BoundValue };

  ORDER = TDatabase::ParamDB->ANSATZ_ORDER;
  
//======================================================================
// loop over all levels (not a multigrid level but for convergence study)
//======================================================================

  for(i=0;i<LEVELS;i++)
   {
    OutPut("*******************************************************" << endl);
    OutPut("******           GEOMETRY  LEVEL ");
    OutPut(i << "              *******" << endl);
    OutPut("*******************************************************" << endl);
    solver_time = 0.0;
    N_LinIter = 0;
    OutPut("memory before: " << setw(10) << GetMemory() << endl);

    if(i)
     Domain->RegRefineAll();

    coll=Domain->GetCollection(It_Finest, 0);
    OutPut( "number of cells: " << coll->GetN_Cells() << endl);
    coll->GetHminHmax(&hmin,&hmax);
    OutPut("h_min : " << hmin << " h_max : " << hmax << endl);
    cout << endl << endl;

//======================================================================
// construct all finite element spaces
//======================================================================
    scalar_space =  new TFESpace2D(coll, Name, Description, BoundCondition, ORDER, NULL);

    N_DOF = scalar_space->GetN_DegreesOfFreedom();
    N_Active = scalar_space->GetActiveBound();
    N_NonActive = N_DOF - N_Active;
    OutPut("dof all      : "<< setw(10) << N_DOF  << endl);

    if(ReadVelo)
     {
      velocity_space =  new TFESpace2D(coll, UName, UString, BoundCondition, 2, NULL);
     }
/*
     exit(0);*/
//======================================================================
// construct all finite element functions
//======================================================================
    N_Unknowns = N_DOF;
    sol = new double[N_Unknowns];
    PostSol = new double[N_Unknowns];
    oldsol = new double[N_Unknowns];
    rhs = new double[N_Unknowns];
    defect = new double[N_Unknowns];

    memset(sol, 0, N_Unknowns*SizeOfDouble);
    memset(PostSol, 0, N_Unknowns*SizeOfDouble);
    memset(oldsol, 0, N_Unknowns*SizeOfDouble);
    memset(rhs, 0, N_Unknowns*SizeOfDouble);

    C = new TFEFunction2D(scalar_space, UString, UString, sol, N_DOF);
 
#ifdef _OPTTOMOGRAPHY_
    sol_log = new double[N_Unknowns]; 
      
//    solution in log scale
   C_log = new TFEFunction2D(scalar_space, UlogString, UlogString, sol_log, N_DOF);
#endif    
    
//======================================================================
// produce outout
//======================================================================
    Output = new TOutput2D(2, 2, 1, 1,Domain);

    Output->AddFEFunction(C);

#ifdef _OPTTOMOGRAPHY_
    for(j=0;j<N_Unknowns;j++)
      if(sol[j]>0)
       sol_log[j] = log(sol[j]);
      else
       sol_log[j] =0.;    
    
    Output->AddFEFunction(C_log); 
#endif      
    
    C->Interpolate(Exact);   
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
    
// exit(0);    
//======================================================================
// allocate memory for all matrices
//======================================================================
    sqstructureA = new TSquareStructure2D(scalar_space);
    sqstructureA->Sort();

    sqmatrixA = new TSquareMatrix2D(sqstructureA);

//======================================================================
// assemble all matrices
//======================================================================
    RHSs[0] = rhs;
    memset(rhs, 0, N_DOF*SizeOfDouble);

    fesp[0] = scalar_space;
    ferhs[0] = scalar_space;

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
	               OutPut(" TEsssssssssssss DISCTYPE " << GLS << endl);
      break;

      default:
            OutPut("Unknown DISCTYPE" << endl);
            exit(4711);;
     }

      // initialize matrices
    SQMATRICES[0] = sqmatrixA;
    SQMATRICES[0]->Reset();
    aux =  new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);
    
    
    TDatabase::ParamDB->P4 = 0;
      // assemble
    Assemble2D(1, fesp,
               1, SQMATRICES,
               0, NULL,
               1, RHSs, ferhs,
               DiscreteForm,
               BoundaryConditions,
               BoundaryValues,
               aux);

     delete aux;

     cout << "Max F : " << TDatabase::ParamDB->P4 << endl;
     
/*   memcpy(sol, rhs, N_DOF*SizeOfDouble); 
     
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
    
exit(0); */     
     
#ifdef _OPTTOMOGRAPHY_          
   RobinInt(sqmatrixA, rhs, BoundCondition);  
#endif      
     
  
// exit(0);   
     
// apply local projection stabilization method
     if(TDatabase::ParamDB->LP_FULL_GRADIENT>0)
      {
       if(TDatabase::ParamDB->LP_FULL_GRADIENT==1)
        UltraLocalProjection(SQMATRICES[0], ForPressure);
//        else if(TDatabase::ParamDB->REACTOR_P20==1)
//         StreamlineUltraLocalProjection(SQMATRICES[0], TRUE);
//        else if(TDatabase::ParamDB->REACTOR_P20==2)
//         StreamlineUltraLocalProjection(SQMATRICES[0], TRUE);
      }

    // set rhs for Dirichlet nodes
    memcpy(sol+N_Active, rhs+N_Active, N_NonActive*SizeOfDouble);
 
//======================================================================
// solve the system
//======================================================================
    // compute defect
    memset(defect,0,N_Unknowns*SizeOfDouble);

//     OutPut("norm of solution " <<  sqrt(Ddot(N_Active,rhs,rhs))  << endl);
    t1 = GetTime();
    DirectSolver(sqmatrixA, rhs, sol);
    t2 = GetTime();
    OutPut( "time for solving: " << t2-t1 << endl);
//     OutPut("solution " << sqrt(Ddot(N_Unknowns,rhs,sol)) << endl);

 
//     for(i=0;i<N_DOF; i++)
//       sol[i] =   300./pow( (1. - (sol[i]/1500.) ), 5. ) - 300.;

//======================================================================
// post processing
//======================================================================
//     AdaptivePostProcess(C, PostSol, TRUE);

 
//     if(TDatabase::ParamDB->WRITE_PS)
//     {
//       // write grid into an Postscript file
//       os.seekp(std::ios::beg);
//       os << PsBaseName << i << ".ps" << ends;
//       Domain->PS(os.str().c_str(),It_Finest,0);
//     }



//     C->Interpolate(Exact);
//     memcpy(PostSol, sol, N_Active*SizeOfDouble);

//  Refine and output the solution


//     Domain->RegRefineAll();
//     coll=Domain->GetCollection(It_Finest, 0);
//     RefinedOutput_space[0] =  new TFESpace2D(coll, Name, Description,
//                                    BoundCondition, 1., NULL);
//     N_RDOF[0] = RefinedOutput_space[0]->GetN_DegreesOfFreedom();
//     refined_sol[0] = new double[N_RDOF[0]];
//     auxsol =  new double[N_RDOF[0]];
//     memset(refined_sol[0], 0, N_RDOF[0]*SizeOfDouble);
//     memset(auxsol, 0, N_RDOF[0]*SizeOfDouble);
//     Refinmed_C = new TFEFunction2D(RefinedOutput_space[0], RUString, RUString, refined_sol[0], N_RDOF[0]);
//     Prolongate(scalar_space, RefinedOutput_space[0], sol, refined_sol[0], auxsol);

//     Domain->RegRefineAll();
//     coll=Domain->GetCollection(It_Finest, 0);
//     RefinedOutput_space[1] =  new TFESpace2D(coll, Name, Description,
//                                    BoundCondition, 1., NULL);
//     N_RDOF[1] = RefinedOutput_space[1]->GetN_DegreesOfFreedom();
//     refined_sol[1] = new double[N_RDOF[1]];
//     auxsol =  new double[N_RDOF[1]];
//     memset(refined_sol[1], 0, N_RDOF[1]*SizeOfDouble);
//     memset(auxsol, 0, N_RDOF[1]*SizeOfDouble);
//     Refinmed_C = new TFEFunction2D(RefinedOutput_space[1], RUString, RUString, refined_sol[1], N_RDOF[1]);
//     Prolongate(RefinedOutput_space[0], RefinedOutput_space[1], refined_sol[0], refined_sol[1], auxsol);
// 


// 
// exit(0);
// 
//     memset(refined_sol[0], 0, N_RDOF[0]*SizeOfDouble);
//     memset(auxsol, 0, N_RDOF[0]*SizeOfDouble);
//     Prolongate(scalar_space, RefinedOutput_space[0], PostSol, refined_sol[0], auxsol);
//     memset(refined_sol[1], 0, N_RDOF[1]*SizeOfDouble);
//     memset(auxsol, 0, N_RDOF[1]*SizeOfDouble);
//     Prolongate(RefinedOutput_space[0], RefinedOutput_space[1], refined_sol[0], refined_sol[1], auxsol);


    if(TDatabase::ParamDB->WRITE_VTK)
     {
       
#ifdef _OPTTOMOGRAPHY_
    for(j=0;j<N_Unknowns;j++)
      if(sol[j]>0)
       sol_log[j] = log(sol[j]);
      else
       sol_log[j] =0.;    
#endif  
       
      os.seekp(std::ios::beg);
       if(img<10) os <<  "VTK/"<<VtkBaseName<<".0000"<<img<<".vtk" << ends;
         else if(img<100) os <<  "VTK/"<<VtkBaseName<<".000"<<img<<".vtk" << ends;
          else if(img<1000) os <<  "VTK/"<<VtkBaseName<<".00"<<img<<".vtk" << ends;
           else if(img<10000) os <<  "VTK/"<<VtkBaseName<<".0"<<img<<".vtk" << ends;
            else  os <<  "VTK/"<<VtkBaseName<<"."<<img<<".vtk" << ends;
      Output->WriteVtk(os.str().c_str());
      
      os.seekp(std::ios::beg);
       if(img<10) os <<  "VTK/"<<VtkBaseName<<".0000"<<img<<".plt" << ends;
         else if(img<100) os <<  "VTK/"<<VtkBaseName<<".000"<<img<<".plt" << ends;
          else if(img<1000) os <<  "VTK/"<<VtkBaseName<<".00"<<img<<".plt" << ends;
           else if(img<10000) os <<  "VTK/"<<VtkBaseName<<".0"<<img<<".plt" << ends;
            else  os <<  "VTK/"<<VtkBaseName<<"."<<img<<".plt" << ends;
      
      Output->WriteBinaryPlt(os.str().c_str());
      img++;
//       Output->ParMooN_WriteVTK(img);
     }
     
//     PrintScalar(C, N_BData);
     
//     exit(0);
//     ComputeExtremalValues(N_Unknowns, C,errors);
//     OutPut(setprecision(4) << "C min:= " << errors[0] << ", C max:= " << errors[1] -1<<endl);

    // measure errors to known solution
    if(TDatabase::ParamDB->MEASURE_ERRORS)
     {
      aux =  new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);

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


    OutPut("memory after: " << setw(10) << GetMemory() << endl);
   } // for(i=0;i<LEVELS;i++)


#ifdef _MPI    
   MPI_Finalize();  
#endif  
  CloseFiles();
  OutPut("used time: " << GetTime() << endl);
  OutPut("used bytes: " << GetMemory() << endl);
  return 0;
}
