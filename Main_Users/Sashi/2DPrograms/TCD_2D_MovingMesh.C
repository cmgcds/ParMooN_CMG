// =======================================================================
//
// Purpose:     main program with parallel solver (no multigrid solver)

//
// Author:      Sashikumaar Ganesan
//
// History:     Implementation started on 19.03.2013
//        :      
// =======================================================================
#include <Domain.h>
#include <Database.h>
#include <SquareStructure1D.h>
#include <SquareMatrix1D.h>
#include <FEFunction1D.h>
#include <FEDatabase2D.h>
#include <FESpace2D.h>
#include <SquareStructure2D.h>
#include <Structure2D.h>
#include <AuxParam2D.h>
#include <QuadAffin.h>
#include <Solver.h>
#include <DirectSolver.h>
#include <Assemble2D.h>
#include <Output2D.h>
#include <DiscreteForm2D.h>
#include <LinAlg.h>
#include <Collection.h>
#include <JointCollection.h>
// #include <ConvDiff2D_Routines.h>
#include <LocalProjection.h>
#include <BaseFunct1D.h>
#include <NodalFunctional1D.h>
#include <QuadBilinear.h>
#include <QuadIsoparametric.h>
#include <TriaIsoparametric.h>

#include <BoundEdge.h>
#include <IsoBoundEdge.h>
#include <InterfaceJoint.h>
#include <ADISystem.h>
#include <ADISystem1D.h>
#include <LineAffin.h>
#include <FreeSurface2D.h>

#include <string.h>
#include <sstream>
#include <MooNMD_Io.h>
#include <stdio.h>
#include <stdlib.h>
// #include <malloc.h>

double bound = 0;

#include <MainUtilities.h>
// #include <TimeUtilities.h>
#include <TimeDiscRout.h>
#include <FEM_TVD_FCT.h>

#include <MacroCell.h>
#include <BoundEdge.h>
#include <IsoBoundEdge.h>
// #include <gridgen.h>
#include <IsoInterfaceJoint.h>
#include <BdLine.h>
#include <BdCircle.h>
#include <GridCell.h>

#define AMG 0
#define GMG 1
#define DIRECT 2

 
#include <TimeConvDiff2D.h>

// =======================================================================
// include current example
// =======================================================================
// #include  "../Examples/TCD_2D/exp_0.h"
// #include  "../Examples/TCD_2D/TimeDomian.h"
// #include  "../Examples/TCD_2D/Gaussian-BenchMark.h"
#include  "../TCD_2D/TimeDomianNoSource.h"
// #include  "../Examples/TCD_2D/Bulk_Academic_Test.h"
// #include  "../Examples/TCD_2D/JohnMaubachTobiska1997inst.h"
// #include  "../Examples_All/TCD_2D/Hemker.h"
// #include  "../TCD_2D/TimeDomian_beam.h"
// #include  "../TimeDomianNoSource.h"

int main(int argc, char* argv[])
{
  TDomain *Domain = new TDomain();
  TDomain *Domain_Intl = new TDomain();
  TDatabase *Database = new TDatabase();
  TFEDatabase2D *FEDatabase = new TFEDatabase2D();
  TCollection *coll, *mortarcoll = NULL;  
  TFESpace2D *Scalar_Space, *Velocity_Spaces, *pressure_space, *Grid_space;
  TOutput2D *Output;
  TAuxParam2D *aux; 
  MatVecProc *MatVect;
  DefectProc *Defect;
  TDiscreteForm2D *DiscreteForm, *DiscreteFormGrid, *DiscreteFormMRhs_SUPG;
  TDiscreteForm2D *DiscreteFormMatrixMARhs, *DiscreteFormMatrixMARhs_SUPG;
  TDiscreteForm2D *DiscreteFormMatrixMRhs, *DiscreteFormMatrixARhs;
  TDiscreteForm2D *DiscreteFormMatrixMRhs_SUPG, *DiscreteFormMatrixARhs_SUPG;
  TFEFunction2D *C, *gridvelo1, *gridvelo2;
  TSquareMatrix2D *MatricesM, *MatricesA, *MatricesM_Old, *MatricesA_Old, *SQMATRICES[3];
  TSquareMatrix2D *MatricesK;
  TSquareMatrix2D *MatG11,  *MatG12, *MatG21, *MatG22, *SQMATRICES_GRID[4]; 
  TSquareStructure2D *sqstructureA, *sqstructureG;  
  TFEVectFunct2D *GridVelo, *GridPos, *RefGridPos;
  TFESpace2D *fesp[2], *ferhs[1];
  TFEFunction2D  *fefct[4];
  TSquareMatrix **sqmatrices = (TSquareMatrix **)SQMATRICES;
  BoundCondFunct2D *BDCond[4];
  BoundValueFunct2D *BDValue[4];  
  MultiIndex2D AllDerivatives[3] = { D00, D10, D01 };
  TBdLine *UpdateBound[4];
  TIsoBoundEdge **Free_Joint, *IsoJoint;
  TVertex **MovBoundVert;
  TBaseCell **Free_Cells;
  BoundCondFunct2D *GridBoundaryConditions[1];
  BoundValueFunct2D *GridBoundValues[1];

  
  int i, j, k, l, m, N, ret, img=1, N_Cells, N_U, Max_It;  
  int VeloFunction, N_Active, N_GridDOFs, N_FESpaces;
  int N_SquareMatrices, N_SubSteps, time_discs;
  int very_first_time=0, N_Rhs, N_MovVert;
  int **IsoCellEdgeNos, N_GridActive, N_GridBDDOFs, *GridKCol, *GridRowPtr;
 
  double total_time, limit, t3, gamma, tau, oldtau, end_time, t4, tau1;
  double *Sol, *B, *Rhs, *W, *gridpos, *gridpos_old, *gridpos_ref, *RHSs[3], *defect, *griddisp;
  double errors[9], olderror=0., olderror1=0., p1, p2, L2error_Max=-1.e8, L2error_Max_t;
  double InitialMass, Mass, CVolume, InitialVolume, *Entries[4];
  double hmin, hmax, hmax_all, hmax_init, *Iso_refX, *GridRhs;
  
  char *PsBaseName, *GrapeBaseName, *GnuBaseName, *ReadGrapeBaseName;
  char *VtkBaseName, *MatlabBaseName, *GmvBaseName;
  char *PRM, *GEO; 
  
  char NameString[] = "C";
  char WString[] = "w";
  char UString[] = "u";
  char PString[] = "p";
  char Readin[] = "readin.dat";
  char MassMatrix[] = "Mass matrix";
//   char Mass[] = "Mass";
  char Name[] = "name";

  bool UpdateConvection=TRUE, UpdateRhs=TRUE, ConvectionFirstTime=TRUE, ReAssembleM=TRUE;
  
  std::ostringstream os;
  os << " "; 
  
  //======================================================================
  // read parameter file
  //======================================================================
  total_time = GetTime();
  if(argc>=2)
    { ret=Domain->ReadParam(argv[1]); }
    else
      { ret=Domain->ReadParam(Readin); }
      
  OpenFiles();
  OutFile.setf(std::ios::scientific);
  
  // write parameters into outfile
  Database->WriteParamDB(argv[0]);
  Database->WriteTimeDB();
  ExampleFile();      

#ifdef __CONSERVATIVEALE__        
 TDatabase::ParamDB->P6=1 ;
#else
 TDatabase::ParamDB->P6=0;
#endif  
  
  
  PRM = TDatabase::ParamDB->BNDFILE;
  GEO = TDatabase::ParamDB->GEOFILE;
      
  // assign names for output files
  PsBaseName = TDatabase::ParamDB->PSBASENAME;
  GrapeBaseName = TDatabase::ParamDB->GRAPEBASENAME;
  GnuBaseName = TDatabase::ParamDB->GNUBASENAME;
  ReadGrapeBaseName = TDatabase::ParamDB->READGRAPEBASENAME;
  GmvBaseName = TDatabase::ParamDB->GMVBASENAME;
  VtkBaseName = TDatabase::ParamDB->VTKBASENAME;
  MatlabBaseName = TDatabase::ParamDB->MATLABBASENAME;

  // pointers to the routines which compute matrix-vector
  // products and the defect
  MatVect = MatVect_Scalar;
  Defect = Defect_Scalar;      
      
  //======================================================================
  // initialize discrete forms
  //======================================================================
   InitializeDiscreteFormsScalar(DiscreteFormMatrixMRhs, DiscreteFormMatrixARhs, DiscreteFormMatrixMRhs_SUPG,
                                 DiscreteFormMatrixARhs_SUPG,  BilinearCoeffs);  
 
   InitializeDiscreteForms_ScalarMoving(DiscreteFormMatrixMARhs, DiscreteFormGrid, DiscreteFormMRhs_SUPG, DiscreteFormMatrixMARhs_SUPG,
                                  BilinearCoeffs, GridCoeffs);
   
  GridBoundaryConditions[0] = GridBoundCondition;
  GridBoundValues[0] = GridBoundValue;
  
   //======================================================================
  // read boundary parameterization and initialize coarse grid
  //======================================================================
  Domain->Init(PRM, GEO);

#ifdef __MOVINGMESH__ 
#if defined(__HEMKER__) || defined(__BEAM__) 
  TriaReMeshGen(Domain);  
  TDatabase::ParamDB->UNIFORM_STEPS=0;  
#else
  Domain->ConvertQuadToTri(TDatabase::ParamDB->CONVERT_QUAD_TO_TRI);  
#endif   
#endif  
 
  // refine grid up to the coarsest level
  for(i=0;i<TDatabase::ParamDB->UNIFORM_STEPS;i++)
    Domain->RegRefineAll();
  
  if(TDatabase::ParamDB->WRITE_PS)
   {
    // write grid into an Postscript file
    os.seekp(std::ios::beg);
    os << "Domain" << ".ps" << ends;
    Domain->PS(os.str().c_str(),It_Finest,0);
   }
//   exit(0);
  
  // initializ time
  TDatabase::TimeDB->CURRENTTIME = TDatabase::TimeDB->STARTTIME;
  TDatabase::TimeDB->CURRENTTIMESTEPLENGTH = TDatabase::TimeDB->TIMESTEPLENGTH;
  SetTimeDiscParameters();
  limit = TDatabase::ParamDB->SC_NONLIN_RES_NORM_MIN_SCALAR;
  Max_It = TDatabase::ParamDB->SC_NONLIN_MAXIT_SCALAR;

  t3 = GetTime();
  total_time = t3 - total_time;
  SetPolynomialDegree();

  // check the example file, to activate
  VeloFunction = (int)TDatabase::ParamDB->INTERNAL_CONVECTION_EQ_VELOFIELD;    
    
  coll=Domain->GetCollection(It_Finest, 0);
  N_Cells = coll->GetN_Cells();
  OutPut("N_Cells  : " << N_Cells <<endl);   
    
//   exit(0);
  
#ifdef __MOVINGMESH__   
#if defined(__HEMKER__) || defined(__BEAM__) 
   IsoCellEdgeNos = new int *[2];
   GetMovingBoundData(coll, N_MovVert, MovBoundVert, Free_Joint, Free_Cells, IsoCellEdgeNos, Iso_refX); 
   OutPut("N_MovVert  : " << N_MovVert <<endl);  
#endif   
#endif  
     
     
   //=========================================================================
   // construct all finite element spaces
   //========================================================================= 
   // fespaces for scalar equations
   Scalar_Space  =  new TFESpace2D(coll, Name, NameString, BoundCondition, TDatabase::ParamDB->ANSATZ_ORDER, NULL);

   N_U = Scalar_Space->GetN_DegreesOfFreedom();
   N_Active = Scalar_Space->GetActiveBound();
   OutPut("Dof       : "<< setw(10) << N_U  << endl); 
   //=========================================================================
   // memory allocate all vectors and construction of all fefunction
   //=========================================================================      
   Sol =  new double[N_U];
   B =  new double[N_U];    
   Rhs = new double[N_U];
   defect = new double[N_U];
    
   memset(Sol, 0, N_U*SizeOfDouble);
   memset(B, 0, N_U*SizeOfDouble);
   memset(Rhs, 0, N_U*SizeOfDouble); 
      
   C = new TFEFunction2D(Scalar_Space, NameString, NameString, Sol, N_U);    
   C->Interpolate(InitialCondition);    
  
   //=========================================================================
   // allocate memory for all matrices
   //=========================================================================
   // build matrices
   // first build matrix structure
   sqstructureA = new TSquareStructure2D(Scalar_Space);
   sqstructureA->Sort();  // sort column numbers: numbers are increasing
   
   MatricesM = new TSquareMatrix2D(sqstructureA);
   MatricesA = new TSquareMatrix2D(sqstructureA);
   
   MatricesM_Old = new TSquareMatrix2D(sqstructureA);
//    MatricesA_Old = new TSquareMatrix2D(sqstructureA);
   
   if(TDatabase::ParamDB->DISCTYPE == SDFEM)
    {
     // stabilisation matrix K
     MatricesK = new TSquareMatrix2D(sqstructureA);
    }  
    
#ifdef __MOVINGMESH__     
   TDatabase::TimeDB->CURRENTTIME=0.;
   Grid_space = new TFESpace2D(coll, WString, WString, GridBoundCondition, 1, NULL); 
   N_GridDOFs = Grid_space->GetN_DegreesOfFreedom(); 
   N_GridActive = Grid_space->GetActiveBound();
   N_GridBDDOFs = N_GridDOFs-N_GridActive; 
//    GlobalNumbers_Grid = Grid_space->GetGlobalNumbers();
//    BeginIndex_Grid    = Grid_space->GetBeginIndex();   
   
   W = new double[2*N_GridDOFs];
   memset(W, 0, 2*N_GridDOFs*SizeOfDouble);

   
   // vector fe function
   GridVelo = new TFEVectFunct2D(Grid_space, WString, WString, W, N_GridDOFs, 2);
   // individual components of velocity
   gridvelo1 = GridVelo->GetComponent(0);
   gridvelo2 = GridVelo->GetComponent(1); 
   
   gridpos = new double[2*N_GridDOFs];
   gridpos_old = new double[2*N_GridDOFs];   
   gridpos_ref = new double[2*N_GridDOFs];   
   griddisp = new double[2*N_GridDOFs];   
   
   memset(gridpos, 0, 2*N_GridDOFs*SizeOfDouble);
   GridPos = new TFEVectFunct2D(Grid_space, WString, WString, gridpos, N_GridDOFs, 2);  
   GridPos->GridToData();
   
#if defined(__HEMKER__) || defined(__BEAM__)  
   GridRhs = new double[2*N_GridDOFs];
   
   RefGridPos = new TFEVectFunct2D(Grid_space, WString, WString, gridpos_ref, N_GridDOFs, 2);     
#endif   
 
   memcpy(gridpos_old, gridpos, 2*N_GridDOFs*SizeOfDouble);
   memcpy(gridpos_ref, gridpos, 2*N_GridDOFs*SizeOfDouble); 
   
   UpdateConvection =TRUE;
   
   
  // for mesh movement      
   sqstructureG = new TSquareStructure2D(Grid_space);
   sqstructureG->Sort();  // sort column numbers: numbers are increasing
   
   MatG11 = new TSquareMatrix2D(sqstructureG); // G11
   MatG12 = new TSquareMatrix2D(sqstructureG); // G12
   MatG21 = new TSquareMatrix2D(sqstructureG); // G21
   MatG22 = new TSquareMatrix2D(sqstructureG); // G22  
   
   
//  ====================================================================================  
// assemble matrix for grid moving - begin
//  ====================================================================================    
    fesp[0] = Grid_space;
    SQMATRICES_GRID[0] = MatG11;
    SQMATRICES_GRID[0]->Reset();
    SQMATRICES_GRID[1] = MatG12;
    SQMATRICES_GRID[1]->Reset();
    SQMATRICES_GRID[2] = MatG21;
    SQMATRICES_GRID[2]->Reset();
    SQMATRICES_GRID[3] = MatG22;
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
    delete aux;      
   
     Entries[0] = MatG11->GetEntries();
     Entries[1] = MatG12->GetEntries();
     Entries[2] = MatG21->GetEntries();
     Entries[3] = MatG22->GetEntries();

     GridKCol = sqstructureG->GetKCol();
     GridRowPtr = sqstructureG->GetRowPtr();

  // for Dirichlet rows in off-diagonal matrices
  memset(Entries[1] + GridRowPtr[N_GridActive], 0, (GridRowPtr[N_GridDOFs] - GridRowPtr[N_GridActive])*SizeOfDouble);
  memset(Entries[2] + GridRowPtr[N_GridActive], 0, (GridRowPtr[N_GridDOFs] - GridRowPtr[N_GridActive])*SizeOfDouble);  
#endif  //__MOVINGMESH__     

   //=========================================================================
   //create output
   //=========================================================================    
    Output = new TOutput2D(1, 1, 1, 1, Domain);   
    
#ifdef __MOVINGMESH__  
    Output->AddFEVectFunct(GridVelo);
#endif     
      
    Output->AddFEFunction(C); 
    os.seekp(std::ios::beg);
    Output->AddParameter(TDatabase::TimeDB->CURRENTTIME,os.str().c_str());    
     
    //======================================================================
    // assembling of mass matrix and initial rhs (f_0)
    //======================================================================
    // set parameters
    N_Rhs = 1;
    N_FESpaces = 1;
    fesp[0] = Scalar_Space;

    // reset matrices
    N_SquareMatrices = 1;
    SQMATRICES[0] = MatricesM;
    SQMATRICES[0]->Reset();
       
    DiscreteForm = DiscreteFormMatrixMRhs;           
    
    if(TDatabase::ParamDB->DISCTYPE == SDFEM)
    {
     N_SquareMatrices = 2;  
     SQMATRICES[1] = MatricesK;
     SQMATRICES[1]->Reset();  
     
     DiscreteForm = DiscreteFormMatrixMRhs_SUPG;    
    }   

#ifdef __MOVINGMESH__  
    N_FESpaces = 2;
    fesp[1] = Grid_space;
    fefct[0] = gridvelo1;
    fefct[1] = gridvelo2;
    fefct[2] = gridvelo1; //for calculating divergence of w
    fefct[3] = gridvelo2; //for calculating divergence of w
    
     // reset matrices
     N_SquareMatrices = 2;
     SQMATRICES[0] = MatricesA;
     SQMATRICES[0]->Reset();
     SQMATRICES[1] = MatricesM;
     SQMATRICES[1]->Reset();      
     
     DiscreteForm = DiscreteFormMatrixMARhs;  
     
   if(TDatabase::ParamDB->DISCTYPE == SDFEM)
    {
     N_SquareMatrices = 3;  
     SQMATRICES[2] = MatricesK;
     SQMATRICES[2]->Reset();   
     
     DiscreteForm = DiscreteFormMatrixMARhs_SUPG;   
    }

    // defined in TimeConvDiff2D.h
     aux =  new TAuxParam2D(TimeCDParamsVeloFieldN_FESpaces,
                            TimeCDParamsVeloFieldN_Fct,
                            TimeCDParamsVeloFieldN_ParamFct,
                            TimeCDParamsVeloFieldN_FEValues_ALE,
                            fesp+1, fefct,
                            TimeCDParamsVeloFieldFct_ALE,
                            TimeCDParamsVeloFieldFEFctIndex_ALE,
                            TimeCDParamsVeloFieldFEMultiIndex_ALE,
                            TimeCDParamsVeloFieldN_Params_ALE,
                            TimeCDParamsVeloFieldBeginParam);

#else
      aux =  new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL); 
#endif    

      BDCond[0] = BoundCondition;
      BDValue[0] = BoundValue;
      
      memset(Rhs, 0, N_U*SizeOfDouble);
      RHSs[0] = Rhs;
      ferhs[0] = Scalar_Space;

      Assemble2D(N_FESpaces, fesp,
                 N_SquareMatrices, SQMATRICES,
                 0, NULL,
                 N_Rhs, RHSs, ferhs,
                 DiscreteForm,
                 BDCond,
                 BDValue,
                 aux);
      delete aux;

      // copy Dirichlet values from rhs into sol
      memcpy(Sol+N_Active, RHSs[0]+N_Active, (N_U-N_Active)*SizeOfDouble);  

      
      
#ifdef __CONSERVATIVEALE__      
      MatricesM_Old->Reset();   
      MatAdd(MatricesM_Old, MatricesM, 1.);        
//       MatricesA_Old->Reset();   
//       MatAdd(MatricesA_Old, MatricesA, 1.);    
#endif      
      
  // parameters for time stepping scheme
  gamma = 0;
  m = 0;
  N_SubSteps = GetN_SubSteps();
  oldtau = 1;
  end_time = TDatabase::TimeDB->ENDTIME;
  for(j=0;j<9;j++)
    errors[j] = 0.; 
  
  // not active : TDatabase::TimeDB->TIMESTEPLENGTH_CONTROL = 0
  if (TDatabase::TimeDB->TIMESTEPLENGTH_CONTROL)
    { time_discs = 2; }
    else
      { time_discs = 1; }

  if(TDatabase::ParamDB->WRITE_VTK)
   {  
        os.seekp(std::ios::beg);
        if(img<10) os << VtkBaseName<<".0000"<<img<<".vtk" << ends;
        else if(img<100) os << VtkBaseName<<".000"<<img<<".vtk" << ends;
        else if(img<1000) os << VtkBaseName<<".00"<<img<<".vtk" << ends;
        else if(img<10000) os << VtkBaseName<<".0"<<img<<".vtk" << ends;
        else  os << VtkBaseName<<"."<<img<<".vtk" << ends;
        Output->WriteVtk(os.str().c_str());
     img++;
     
   }

    // measure errors to known solution
    if(TDatabase::ParamDB->MEASURE_ERRORS)
     {
      aux =  new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);
      
      C->GetErrors(Exact, 3, AllDerivatives, 2, L2H1Errors, BilinearCoeffs, aux, 1, fesp, errors);
      
      delete aux;   
      
      InitialVolume = CVolume = Volume(Scalar_Space);      
      Mass = InitialMass = errors[0]/CVolume; // only if exact solution is zero     
      OutPut("T, Volume , Volume diff : " << TDatabase::TimeDB->CURRENTTIME<<"   "<< InitialVolume<<"   "<<CVolume - InitialVolume<<endl);  
 
      OutPut("time: " << TDatabase::TimeDB->CURRENTTIME);
      OutPut(" L2: " << errors[0]);
//       OutPut(" C-Mass: " << Mass << " C-Mass Diff: " << InitialMass  - Mass<< " Relative C-Mass diff: " << ( InitialMass  - Mass)/InitialMass << endl);

      errors[3] += (errors[0]*errors[0] + olderror * olderror)*TDatabase::TimeDB->TIMESTEPLENGTH/2.0;
      olderror = errors[0];
      OutPut(TDatabase::TimeDB->CURRENTTIME <<  " L2(0,T;L2) " << sqrt(errors[3]) << " ");

      errors[4] += (errors[1]*errors[1] +olderror1 * olderror1)*TDatabase::TimeDB->TIMESTEPLENGTH/2.0;
      OutPut( "L2(0,T;H1) " << sqrt(errors[4]) << endl);
      olderror1 = errors[1];            
     } //  if(TDatabase::ParamDB->MEASURE_ERRORS)   
     
    coll->GetHminHmax(&hmin,&hmax);
    hmax_init = hmax_all = hmax;
 
/*#ifdef __MOVINGMESH__       
    hmax *= TDatabase::ParamDB->P6;  
#endif   */ 
    
//  tau = pow(hmax, TDatabase::ParamDB->REACTOR_P30);     
//  TDatabase::TimeDB->TIMESTEPLENGTH = tau;   
    OutPut("h_min : " << hmin << " h_max : " << hmax << endl);
    OutPut("TIMESTEPLENGTH : " << TDatabase::TimeDB->TIMESTEPLENGTH << endl;)  
//    exit(0);
  //======================================================================
  // start of time cycle
  //======================================================================
  tau = TDatabase::TimeDB->CURRENTTIME;
  oldtau=tau;
  
  while(TDatabase::TimeDB->CURRENTTIME< end_time)
  {                                               // time cycle
    m++;
    TDatabase::TimeDB->INTERNAL_STARTTIME = TDatabase::TimeDB->CURRENTTIME;

    for(l=0;l<N_SubSteps;l++) // sub steps of fractional step theta
    {
      if (!very_first_time)
      {
        SetTimeDiscParameters();
      }

     if(m==1)
      {
        OutPut("Theta1: " << TDatabase::TimeDB->THETA1<< endl);
        OutPut("Theta2: " << TDatabase::TimeDB->THETA2<< endl);
        OutPut("Theta3: " << TDatabase::TimeDB->THETA3<< endl);
        OutPut("Theta4: " << TDatabase::TimeDB->THETA4<< endl);
      }

      tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
     
      if (!very_first_time)
        TDatabase::TimeDB->CURRENTTIME += tau;
   
     OutPut(endl << "CURRENT TIME: ");
     OutPut(TDatabase::TimeDB->CURRENTTIME << endl);   

#ifdef __MOVINGMESH__ 
      t4 = TDatabase::TimeDB->CURRENTTIME - tau/2.;  
   
#ifndef __CONSERVATIVEALE__
      if(TDatabase::TimeDB->TIME_DISC==1)
       {
        t4 = TDatabase::TimeDB->CURRENTTIME;
        ReAssembleM = FALSE;
       }
#endif
      
      tau1 = t4 - oldtau;
      oldtau = t4;         
   cout << " tau1: " << tau1 <<  "  tau: " <<  tau  << "  oldtau: " <<  oldtau  << "  t4: " <<  t4  <<endl;
// exit(0);

   
   /** find W(t4) and move mesh to t4 */
#if defined(__HEMKER__) || defined(__BEAM__) 
     GridPos->GridToData();   
     memcpy(gridpos_old, gridpos, 2*N_GridDOFs*SizeOfDouble);

     RefGridPos->DataToGrid();    
     ModifyBdCoords(N_MovVert, MovBoundVert, Free_Joint, Iso_refX, t4);    
  
     // data with updated BD values
     GridPos->GridToData();   
    
     memset(GridRhs, 0, 2*N_GridDOFs*SizeOfDouble);     
     memcpy(GridRhs+N_GridActive, gridpos+N_GridActive, N_GridBDDOFs*SizeOfDouble);     //rhs1  
     memcpy(GridRhs+(N_GridDOFs+N_GridActive), gridpos+(N_GridDOFs+N_GridActive), N_GridBDDOFs*SizeOfDouble);   //rhs2       
     
     Daxpy(N_GridBDDOFs, -1., gridpos_old+N_GridActive, GridRhs+N_GridActive);
     Daxpy(N_GridBDDOFs, -1., gridpos_old+(N_GridDOFs+N_GridActive), GridRhs+(N_GridDOFs+N_GridActive));      
     
     memcpy(W, GridRhs, 2*N_GridDOFs*SizeOfDouble);   
     
     SolveGridEquation(Entries, W, GridRhs, GridKCol, GridRowPtr, N_GridDOFs);
          
     memcpy(gridpos, gridpos_old, 2*N_GridDOFs*SizeOfDouble);
     Daxpy(2*N_GridDOFs, 1., W, gridpos);
     GridPos->DataToGrid();      
     Dscal(2*N_GridDOFs, 1./tau1, W); // - sign due to -w\cdot\nabla C in the equation                 
#else
     GridPos->GridToData();   
     memcpy(gridpos_old, gridpos, 2*N_GridDOFs*SizeOfDouble);
     //  move mesh and mesh velo calculation       
     for(i=0;i<N_GridDOFs;i++)
        ModifyCoords(gridpos_ref[i], gridpos_ref[i+N_GridDOFs], gridpos[i], gridpos[i+N_GridDOFs], TDatabase::TimeDB->CURRENTTIME);  
 
 
//      GridPos->DataToGrid();           
     coll->GetHminHmax(&hmin,&hmax);
     if(hmax_all < hmax) hmax_all = hmax;
     

     memcpy(W, gridpos, 2*N_GridDOFs*SizeOfDouble);     
     Daxpy(2*N_GridDOFs, -1., gridpos_old, W);        
     Dscal(2*N_GridDOFs, 1./tau, W); // - sign du*/ //e to -w\cdot\nabla C in the equation   
     memcpy(gridpos_old, gridpos, 2*N_GridDOFs*SizeOfDouble); 
     
     
     
     for(i=0;i<N_GridDOFs;i++)
       ModifyCoords(gridpos_ref[i], gridpos_ref[i+N_GridDOFs], gridpos[i], gridpos[i+N_GridDOFs], t4);         

     GridPos->DataToGrid();       
     
//           for(i=0;i<N_GridDOFs;i++)
// 	    cout << i <<" " << gridpos_ref[i] << " " << gridpos_ref[i+N_GridDOFs] << " "  << W[i] << " " << W[i+N_GridDOFs]  <<endl;
// 	    
            OutPut(tau1 << " MeshVelo " << Ddot(2*N_GridDOFs, W, W ) << endl); 
#endif

#endif  
//  cout << " CurrTime " << t4 <<endl;
//      OutPut("MeshVelo " << Ddot(2*N_GridDOFs, W, W ) << endl); 
//      exit(0);
       
     // working array for rhs is B, initialize B
     memset(B, 0, N_U*SizeOfDouble);

     // compute terms with data from previous time step
     // old rhs multiplied with current subtime step and theta3 on B
     Daxpy(N_Active, tau*TDatabase::TimeDB->THETA3, Rhs, B);  

     if(UpdateConvection || UpdateRhs ||  ConvectionFirstTime )
      {    
       // assemble A and rhs
       N_Rhs = 1;
       N_FESpaces = 1;
       fesp[0] = Scalar_Space;

       // reset matrices
       N_SquareMatrices = 1;
       SQMATRICES[0] = MatricesA;
       SQMATRICES[0]->Reset();
       
       DiscreteForm = DiscreteFormMatrixARhs;       
       
       if(TDatabase::ParamDB->DISCTYPE == SDFEM)
       {
        N_SquareMatrices = 2;  
        SQMATRICES[1] = MatricesK;
        SQMATRICES[1]->Reset();  
        
        DiscreteForm = DiscreteFormMatrixARhs_SUPG;    
       }   
       
#ifdef __MOVINGMESH__  
       N_FESpaces = 2;
       fesp[1] = Grid_space;
       fefct[0] = gridvelo1;
       fefct[1] = gridvelo2;
       fefct[2] = gridvelo1; //for calculating divergence of w
       fefct[3] = gridvelo2; //for calculating divergence of w
    
       N_SquareMatrices = 2;     
       SQMATRICES[1] = MatricesM;
       SQMATRICES[1]->Reset();       
       
       DiscreteForm = DiscreteFormMatrixMARhs;     

       if(TDatabase::ParamDB->DISCTYPE == SDFEM)
        {
         N_SquareMatrices = 3;  
         SQMATRICES[2] = MatricesK;
         SQMATRICES[2]->Reset();   
     
         DiscreteForm = DiscreteFormMatrixMARhs_SUPG;   
        }
    
       // defined in TimeConvDiff2D.h
       aux =  new TAuxParam2D(TimeCDParamsVeloFieldN_FESpaces,
                             TimeCDParamsVeloFieldN_Fct,
                             TimeCDParamsVeloFieldN_ParamFct,
                             TimeCDParamsVeloFieldN_FEValues_ALE,
                             fesp+1, fefct,
                             TimeCDParamsVeloFieldFct_ALE,
                             TimeCDParamsVeloFieldFEFctIndex_ALE,
                             TimeCDParamsVeloFieldFEMultiIndex_ALE,
                             TimeCDParamsVeloFieldN_Params_ALE,
                             TimeCDParamsVeloFieldBeginParam);

#else
       aux =  new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL); 
#endif               
 
       BDCond[0] = BoundCondition;
       BDValue[0] = BoundValue;
      
       memset(Rhs, 0, N_U*SizeOfDouble);
       RHSs[0] = Rhs;
       ferhs[0] = Scalar_Space;

       Assemble2D(N_FESpaces, fesp,
                  N_SquareMatrices, SQMATRICES,
                  0, NULL,
                  N_Rhs, RHSs, ferhs,
                  DiscreteForm,
                  BDCond,
                  BDValue,
                  aux);
       delete aux;
            memset(defect, 0, N_U*SizeOfDouble);   
        MatVectActive(MatricesA, Sol, defect);
//         OutPut("MatA-Sol " << Ddot(N_U, defect, defect ) << endl);  
//        OutPut("MatA-Sol " << Ddot(N_U, Sol, Sol ) << endl); 
// 	exit(0);
	
      ConvectionFirstTime = FALSE;  
     } // if(UpdateConvection || U

     // add rhs from current sub time step to rhs array B
     Daxpy(N_Active, tau*TDatabase::TimeDB->THETA4,  Rhs, B);
     
//       OutPut("B1  " << Ddot(N_U, B, B ) << endl);
//       exit(0);
      
   if(ReAssembleM)
    {
     // for conservative ALE, the mass matrix should be in n+1 so assemble M_{n+1} again by keeping A_{n+1/2}     
     t4 = TDatabase::TimeDB->CURRENTTIME ;
     tau1 = t4 - oldtau;     
     oldtau = t4;
     
#if defined(__HEMKER__) || defined(__BEAM__) 
     GridPos->GridToData();   
     memcpy(gridpos_old, gridpos, 2*N_GridDOFs*SizeOfDouble);

     RefGridPos->DataToGrid();    
     ModifyBdCoords(N_MovVert, MovBoundVert, Free_Joint, Iso_refX, t4);    
  
     // data with updated BD values
     GridPos->GridToData();   
    
     memset(GridRhs, 0, 2*N_GridDOFs*SizeOfDouble);     
     memcpy(GridRhs+N_GridActive, gridpos+N_GridActive, N_GridBDDOFs*SizeOfDouble);     //rhs1  
     memcpy(GridRhs+(N_GridDOFs+N_GridActive), gridpos+(N_GridDOFs+N_GridActive), N_GridBDDOFs*SizeOfDouble);   //rhs2       
     
     Daxpy(N_GridBDDOFs, -1., gridpos_old+N_GridActive, GridRhs+N_GridActive);
     Daxpy(N_GridBDDOFs, -1., gridpos_old+(N_GridDOFs+N_GridActive), GridRhs+(N_GridDOFs+N_GridActive));      
     
     memcpy(W, GridRhs, 2*N_GridDOFs*SizeOfDouble);   
     
     SolveGridEquation(Entries, W, GridRhs, GridKCol, GridRowPtr, N_GridDOFs);
          
     memcpy(gridpos, gridpos_old, 2*N_GridDOFs*SizeOfDouble);
     Daxpy(2*N_GridDOFs, 1., W, gridpos);
     GridPos->DataToGrid();      
     Dscal(2*N_GridDOFs, 1./tau1, W); // - sign due to -w\cdot\nabla C in the equation                 
#else
     
//      GridPos->GridToData();   
//      memcpy(gridpos_old, gridpos, 2*N_GridDOFs*SizeOfDouble);
     
     //  move mesh and mesh velo calculation       
     for(i=0;i<N_GridDOFs;i++)
       ModifyCoords(gridpos_ref[i], gridpos_ref[i+N_GridDOFs], gridpos[i], gridpos[i+N_GridDOFs], t4);
/*      ModifyCoords(gridpos_ref[i], gridpos_ref[i+N_GridDOFs], gridpos[i], gridpos[i+N_GridDOFs], t4, W[i], W[i+N_GridDOFs]);    */ 

     GridPos->DataToGrid();
     
//      exit(0);
//      memcpy(W, gridpos, 2*N_GridDOFs*SizeOfDouble);     
//      Daxpy(2*N_GridDOFs, -1, gridpos_old, W);        
//      Dscal(2*N_GridDOFs, -1./tau1, W); // - sign du*/ //e to -w\cdot\nabla C in the equation 
//    OutPut(tau1 << " MeshVelo " << Ddot(2*N_GridDOFs, W, W ) << endl);    
//      memcpy(gridpos_old, gridpos, 2*N_GridDOFs*SizeOfDouble);      
#endif      
       // assemble only the mass matrix, so mesh velo is not needed, note we use inconsist SUPG in CN  
/*#ifdef __TIMECONSISTSUPG__
       // MatricesK is needed only in consistant SUPG
      cout<<"Consist SUPG in CN is not implemented " << endl;
      exit(0);
#endif   */   
       // set parameters
       N_Rhs = 1;
       N_FESpaces = 1;
       fesp[0] = Scalar_Space;

       // reset matrices
       N_SquareMatrices = 1;
       SQMATRICES[0] = MatricesM;
       SQMATRICES[0]->Reset();
       
       DiscreteForm = DiscreteFormMatrixMRhs;         
     
       memset(Rhs, 0, N_U*SizeOfDouble);
       RHSs[0] = Rhs;
       ferhs[0] = Scalar_Space;       
       
       aux =  new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);  

      Assemble2D(N_FESpaces, fesp,
                 N_SquareMatrices, SQMATRICES,
                 0, NULL,
                 N_Rhs, RHSs, ferhs,
                 DiscreteForm,
                 BDCond,
                 BDValue,
                 aux);
      delete aux; 
    } // if(ReAssembleM) 

     // update rhs by Laplacian and convective term from previous time step
     // scaled by current sub time step length and theta2
     // currently : M := M + gamma A
     // M = M + (- tau*TDatabase::TimeDB->THETA2)
     // defect = M * sol
     // B:= B + defec     
     memset(defect, 0, N_U*SizeOfDouble);    
  
#ifdef __CONSERVATIVEALE__        
      MatAdd(MatricesM_Old, MatricesA, - tau*TDatabase::TimeDB->THETA2);
      gamma = 0.;
      MatVectActive(MatricesM_Old, Sol, defect);
      
      // store mass mat for next time step
      MatricesM_Old->Reset();   
      MatAdd(MatricesM_Old, MatricesM, 1.);        
#else          
     MatAdd(MatricesM, MatricesA, - tau*TDatabase::TimeDB->THETA2);
     
/*#ifdef __TIMECONSISTSUPG__     
     if(TDatabase::ParamDB->DISCTYPE == SDFEM)
      {     
       MatAdd(MatricesM, MatricesK, 1.);
      }
#endif   */   
     // set current factor of steady state matrix
     gamma = -tau*TDatabase::TimeDB->THETA2;  
     MatVectActive(MatricesM, Sol, defect);     
#endif    
//      OutPut("Sol " << Ddot(N_U, Sol, Sol ) << endl);   
//        OutPut("defect " << Ddot(N_U, defect, defect ) << endl);   
     
     Daxpy(N_Active, 1, defect, B);
//           OutPut("B1  " << Ddot(N_U, B, B ) << endl);
//      exit(0);
     // set Dirichlet values
     // RHSs[0] still available from assembling
     memcpy(B+N_Active, Rhs+N_Active, (N_U-N_Active)*SizeOfDouble);
     
     // copy Dirichlet values from rhs into sol
     memcpy(Sol+N_Active, Rhs+N_Active, (N_U-N_Active)*SizeOfDouble);

     // system matrix
     MatAdd(MatricesM, MatricesA, -gamma + tau*TDatabase::TimeDB->THETA1);
     gamma = tau*TDatabase::TimeDB->THETA1;
//                 OutPut("B  " << Ddot(N_U, B, B ) << endl);
//======================================================================
// solve linear system
//======================================================================    
       DirectSolver(MatricesM, B, Sol);    
       
//           cout<<Ddot(N_U, Sol, Sol)<<endl;
// 	  exit(0);
//======================================================================
// end solve linear system
//======================================================================   
      
#ifndef __MOVINGMESH__        
     // restore matrices
     MatAdd(MatricesM, MatricesA, -gamma);
     gamma = 0;         
/*#ifdef __TIMECONSISTSUPG__       
     if(TDatabase::ParamDB->DISCTYPE == SDFEM)
      {     
       MatAdd(MatricesM, MatricesK, -1.);
      }
#endif     */ 
#endif      
      

  /** assemble mesh moving matrix */ 
#if defined(__HEMKER__) || defined(__BEAM__)     
   memcpy(gridpos_old, gridpos, 2*N_GridDOFs*SizeOfDouble);
//    memcpy(gridpos_ref, gridpos, 2*N_GridDOFs*SizeOfDouble); 
   
   UpdateConvection =TRUE;
   
   
  // for mesh movement      
   sqstructureG = new TSquareStructure2D(Grid_space);
   sqstructureG->Sort();  // sort column numbers: numbers are increasing
   
   MatG11 = new TSquareMatrix2D(sqstructureG); // G11
   MatG12 = new TSquareMatrix2D(sqstructureG); // G12
   MatG21 = new TSquareMatrix2D(sqstructureG); // G21
   MatG22 = new TSquareMatrix2D(sqstructureG); // G22  
   
   
//  ====================================================================================  
// assemble matrix for grid moving - begin
//  ====================================================================================    
    fesp[0] = Grid_space;
    SQMATRICES_GRID[0] = MatG11;
    SQMATRICES_GRID[0]->Reset();
    SQMATRICES_GRID[1] = MatG12;
    SQMATRICES_GRID[1]->Reset();
    SQMATRICES_GRID[2] = MatG21;
    SQMATRICES_GRID[2]->Reset();
    SQMATRICES_GRID[3] = MatG22;
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
    delete aux;      
   
     Entries[0] = MatG11->GetEntries();
     Entries[1] = MatG12->GetEntries();
     Entries[2] = MatG21->GetEntries();
     Entries[3] = MatG22->GetEntries();

//      GridKCol = sqstructureG->GetKCol();
//      GridRowPtr = sqstructureG->GetRowPtr();

  // for Dirichlet rows in off-diagonal matrices
  memset(Entries[1] + GridRowPtr[N_GridActive], 0, (GridRowPtr[N_GridDOFs] - GridRowPtr[N_GridActive])*SizeOfDouble);
  memset(Entries[2] + GridRowPtr[N_GridActive], 0, (GridRowPtr[N_GridDOFs] - GridRowPtr[N_GridActive])*SizeOfDouble);  
#endif  //__MOVINGMESH__          

   } //  for(l=0;l<N_SubSteps;l++) 
   
    // measure errors to known solution
    if(TDatabase::ParamDB->MEASURE_ERRORS)
     {
      aux =  new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);
      
      C->GetErrors(Exact, 3, AllDerivatives, 2, L2H1Errors, BilinearCoeffs, aux, 1, fesp, errors);
      
      delete aux;  
      
      CVolume = Volume(Scalar_Space);      
      Mass = errors[0]/CVolume;  
      OutPut("T, hmax_all, Volume , Volume diff : " << TDatabase::TimeDB->CURRENTTIME<<"   "<< hmax_all/hmax_init<<"   "<< InitialVolume<<"   "<<CVolume - InitialVolume<<endl);        
      OutPut("time: " << TDatabase::TimeDB->CURRENTTIME);
      OutPut(" L2: " << errors[0]);
      OutPut(" H1-semi: " << errors[1] << endl);
      OutPut(" C-Mass: " << Mass << " C-Mass Diff: " << InitialMass  - Mass<< " Relative C-Mass diff: " << ( InitialMass  - Mass)/InitialMass << endl);
      
      errors[3] += (errors[0]*errors[0] + olderror * olderror)*TDatabase::TimeDB->TIMESTEPLENGTH/2.0;
      olderror = errors[0];
//       OutPut(TDatabase::TimeDB->CURRENTTIME <<  " L2(0,T;L2) " << sqrt(errors[3]) << " ");

      errors[4] += (errors[1]*errors[1] +olderror1 * olderror1)*TDatabase::TimeDB->TIMESTEPLENGTH/2.0;
//       OutPut( "L2(0,T;H1) " << sqrt(errors[4]) << endl);
      olderror1 = errors[1];            
     } //  if(TDatabase::ParamDB->MEASURE_ERRORS)
     
//      exit(0);
     
   if(m==1 || m % TDatabase::TimeDB->STEPS_PER_IMAGE == 0)
    if(TDatabase::ParamDB->WRITE_VTK)
     {
       os.seekp(std::ios::beg);
       if(img<10) os << VtkBaseName<<".0000"<<img<<".vtk" << ends;
       else if(img<100) os << VtkBaseName<<".000"<<img<<".vtk" << ends;
       else if(img<1000) os << VtkBaseName<<".00"<<img<<".vtk" << ends;
       else if(img<10000) os << VtkBaseName<<".0"<<img<<".vtk" << ends;
       else  os << VtkBaseName<<"."<<img<<".vtk" << ends;
       Output->WriteVtk(os.str().c_str());
       img++;
     }
//   exit(0);
  } //  while(TDatabase::TimeDB->CURRENTTIM
        
      
    if(TDatabase::ParamDB->WRITE_VTK)
     {
       os.seekp(std::ios::beg);
       if(img<10) os << VtkBaseName<<".0000"<<img<<".vtk" << ends;
       else if(img<100) os << VtkBaseName<<".000"<<img<<".vtk" << ends;
       else if(img<1000) os << VtkBaseName<<".00"<<img<<".vtk" << ends;
       else if(img<10000) os << VtkBaseName<<".0"<<img<<".vtk" << ends;
       else  os << VtkBaseName<<"."<<img<<".vtk" << ends;
       Output->WriteVtk(os.str().c_str());
       img++;
     }      
      
//    cout <<"__MOVINGMESH 2__ " <<endl;
//    exit(0);
  CloseFiles();
  return 0;  
  
}

















