// =======================================================================
//
// Purpose:     main program for solving a stationary scalar equation using ParMooN
//
// Author:      Sashikumaar Ganesan
//
// History:     Implementation started on 22.08.2014
// =======================================================================
#include <Domain.h>
#include <Database.h>
#include <FEDatabase2D.h>
#include <SystemCD2D.h>
#include <Output2D.h>
#include <MainUtilities.h>
#include <TNSE2D_ParamRout.h>

#include <MooNMD_Io.h>
#include <sys/stat.h>
#include <sys/types.h>

// =======================================================================
// include current example
// =======================================================================
// #include "../Examples/CD_2D/Hemker1996.h" // circle in a channel
#include "../Examples/CD_2D/SineLaplace_hl.h" // smooth sol in unitsquares
// #include "../Examples/CD_2D/TwoInteriorLayers.h" // smooth sol in unitsquares

// =======================================================================
// main program
// =======================================================================
int main(int argc, char* argv[])
{
  //  declaration of database, you need this in every program
  int i, ORDER, N_Cells, N_DOF, N_U, img=0, N_heatfuncDOF;
  
  double *sol, *rhs, *sol_hl, *rhs_hl, t1, t2, errors[4], *Velo;
     
  char *VtkBaseName;
     
  TDatabase *Database = new TDatabase();
  TFEDatabase2D *FEDatabase = new TFEDatabase2D(); 
  TCollection *coll;
  TDomain *Domain;
  TFESpace2D *Scalar_FeSpace, *heatfunc_space, *fesp[2], *Velocity_FeSpace;
  TFEFunction2D *Scalar_FeFunction, *Heatfunc_FeFunction;
  TSystemCD2D *SystemMatrix, *SystemMatrix_HeatLine;
  TOutput2D *Output;
  TAuxParam2D *aux;
  MultiIndex2D AllDerivatives[3] = { D00, D10, D01 };
  TFEVectFunct2D *Velocity;
  TFEFunction2D *u1, *u2, *fefct[3];
  
  std::ostringstream os;
  os << " ";     
  
  // set variables' value in TDatabase using argv[1] (*.dat file) 
  Domain = new TDomain(argv[1]);
  
  //set PROBLEM_TYPE to CD if not yet set
  if(TDatabase::ParamDB->PROBLEM_TYPE == 0)
    TDatabase::ParamDB->PROBLEM_TYPE = 1;
  //open OUTFILE, this is where all output is written to (addionally to console)
  OpenFiles();
 
  // write all Parameters to the OUTFILE (not to console) for later reference
  Database->WriteParamDB(argv[0]);
  ExampleFile();
    
  /* include the mesh from a mesh generator, for a standard mesh use the 
   * build-in function. The GEOFILE describes the boundary of the domain. */
  Domain->Init(NULL, TDatabase::ParamDB->GEOFILE); // call mesh generator
   
  // refine grid up to the coarsest level
  for(i=0; i<TDatabase::ParamDB->UNIFORM_STEPS; i++)
    Domain->RegRefineAll();  
  
  // write grid into an Postscript file
  if(TDatabase::ParamDB->WRITE_PS)
    Domain->PS("Domain.ps", It_Finest, 0);
  
  // create output directory, if not already existing
  if(TDatabase::ParamDB->WRITE_VTK)
    mkdir(TDatabase::ParamDB->OUTPUTDIR, 0777);
  
  // choose example according to the value of TDatabase::ParamDB->EXAMPLE
//   Example_CD2D example;
   
  //=========================================================================
  // construct all finite element spaces
  //=========================================================================
  ORDER = TDatabase::ParamDB->ANSATZ_ORDER;
  
  // a collection is basically only an array of cells, which is needed to create
  // a finite element space
  coll = Domain->GetCollection(It_Finest, 0);
  // print out some information about the mesh
  N_Cells = coll->GetN_Cells();
  OutPut("N_Cells : " << N_Cells <<endl);
  
  // create fespace for scalar equation
  Scalar_FeSpace = new TFESpace2D(coll, (char*)"name", (char*)"C", BoundCondition, ORDER, NULL);

  // print out some information on the finite element space
  N_DOF = Scalar_FeSpace->GetN_DegreesOfFreedom();
  OutPut("dof scalar      : "<< setw(10) << N_DOF  << endl);

  
#ifdef __HEATLINE__  
  heatfunc_space = new TFESpace2D(coll, (char*)"name", (char*)"H", HeatFuncBoundCondition, ORDER, NULL);
  
  N_heatfuncDOF = heatfunc_space->GetN_DegreesOfFreedom();
  OutPut("dof heatfunction      : "<< setw(10) << N_heatfuncDOF  << endl);
  OutPut("dof all      : "<< setw(10) << N_heatfuncDOF + N_DOF  << endl);
  
  Velocity_FeSpace = new TFESpace2D(coll, (char*)"name", (char*)"U", HeatFuncBoundCondition, ORDER, NULL);
  N_U = Velocity_FeSpace->GetN_DegreesOfFreedom();
  OutPut("dof velo      : "<< setw(10) << 2*N_U  << endl);
#endif  
  //======================================================================
  // construct all finite element functions
  //======================================================================
  sol = new double[N_DOF];
  rhs = new double[N_DOF];
  // set solution and right hand side vectors to zero
  memset(sol, 0, N_DOF*SizeOfDouble);
  memset(rhs, 0, N_DOF*SizeOfDouble);

  // create a finite element function
  Scalar_FeFunction = new TFEFunction2D(Scalar_FeSpace, (char*)"C", (char*)"C", sol, N_DOF);
  
#ifdef __HEATLINE__     
  sol_hl = new double[N_heatfuncDOF];
  rhs_hl = new double[N_heatfuncDOF];
  
  memset(sol_hl, 0, N_heatfuncDOF*SizeOfDouble);
  memset(rhs_hl, 0, N_heatfuncDOF*SizeOfDouble);

  // heatfunction fefunction
  Heatfunc_FeFunction = new TFEFunction2D(heatfunc_space, (char*)"H", (char*)"H", sol_hl, N_heatfuncDOF);
  
  Velo = new double[2*N_U];
  Velocity = new TFEVectFunct2D(Velocity_FeSpace, (char*)"U",  (char*)"U",  Velo, N_U, 2);
  u1 = Velocity->GetComponent(0);
  u2 = Velocity->GetComponent(1);  
  
  u1->Interpolate(ExactU1);
  u2->Interpolate(ExactU2);
#endif
    
  //======================================================================
  // SystemMatrix construction and solution
  //====================================================================== 
  // Disc type: GALERKIN (or) SDFEM  (or) UPWIND (or) GLS (or) SUPG (or) LOCAL_PROJECTION
  // Solver: AMG_SOLVE (or) GMG  (or) DIRECT 
  SystemMatrix = new TSystemCD2D(Scalar_FeSpace, GALERKIN, DIRECT);

  // initilize the system matrix with the functions defined in the example
  SystemMatrix->Init(BilinearCoeffs, BoundCondition, BoundValue);
  

  // assemble the system matrix with given aux, sol and rhs 
  // aux is used to pass  addition fe functions (eg. mesh velocity) that is nedded for assembling,
  // otherwise, just pass with NULL 
  SystemMatrix->Assemble(NULL, sol, rhs);

  //Solve the system
    t1 = GetTime();
    SystemMatrix->Solve(sol, rhs);
    t2 = GetTime();
    OutPut( "time for solving: " << t2-t1 << endl);
    
#ifdef __HEATLINE__      
  SystemMatrix_HeatLine = new TSystemCD2D(heatfunc_space, HEATLINE, DIRECT);
  
  // initilize
  SystemMatrix_HeatLine->Init(HeatfuncCoeffs, HeatFuncBoundCondition, HeatFuncBoundValue);
  
    fesp[0] = Scalar_FeSpace;   // thermal space
    fesp[1] = Velocity_FeSpace; // velocity space in all domain
 
    fefct[0] = Scalar_FeFunction; // T
    fefct[1] = u1; // u1
    fefct[2] = u2; // u2

    
    // T, T_x, T_y, u1, u2  parameters are needed for assembling
    // fesp is taken from fefct in aux
    aux =  new TAuxParam2D(N_FESpaces_HeatLine, N_FEFct_HeatLine,
                           N_ParamFct_HeatLine,
                           N_FEValues_HeatLine,
                           fesp, fefct,
                           ParamFct_HeatLineAll,
                           FEFctIndex_HeatLine,
                           FEValueMultiIndex_HeatLine,
                           N_Parameters_HeatLine, BeginParam_HeatLine);
 
   SystemMatrix_HeatLine->Assemble(aux, sol_hl, rhs_hl);  
  
   
   SystemMatrix_HeatLine->Solve(sol_hl, rhs_hl);  
   
#endif  
    
  
    
 
  //======================================================================
  // produce outout
  //======================================================================
  VtkBaseName = TDatabase::ParamDB->VTKBASENAME;    
  Output = new TOutput2D  (1, 1, 0, 0, Domain);
  Output->AddFEFunction(Scalar_FeFunction);

#ifdef __HEATLINE__ 
  Output->AddFEFunction(Heatfunc_FeFunction);
#endif  
  
//     Scalar_FeFunction->Interpolate(Exact);   
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
   
  //====================================================================== 
  // measure errors to known solution
  // If an exact solution is not known, it is usually set to be zero, so that
  // in such a case here only integrals of the solution are computed.
  //======================================================================    
    if(TDatabase::ParamDB->MEASURE_ERRORS)
     {
      fesp[0] = Scalar_FeSpace;
      aux =  new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);

      Scalar_FeFunction->GetErrors(Exact, 3, AllDerivatives, 2, L2H1Errors,
                                   BilinearCoeffs, aux, 1, fesp, errors);

      delete aux;

      OutPut( "L2: " << errors[0] << endl);
      OutPut( "H1-semi: " << errors[1] << endl);
      OutPut( "SD: " << errors[2] << endl);

     } // if(TDatabase::ParamDB->MEASURE_ERRORS)

  CloseFiles();
  return 0;
} // end main
