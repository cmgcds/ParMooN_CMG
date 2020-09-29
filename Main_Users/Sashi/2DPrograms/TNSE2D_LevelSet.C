// =======================================================================
//
// Purpose:     main program for THSE2D with level set
//
// Author:      Sashikumaar Ganesan
//
// History:     Implementation started on 17.07.2012
//        :    
// =======================================================================


#include <Domain.h>
#include <Database.h>
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
#include <LocalProjection.h>

#include <string.h>
#include <sstream>
#include <MooNMD_Io.h>
#include <stdlib.h>
#include <malloc.h>

double bound = 0;

#include <Upwind.h>
#include <FEM_TVD_FCT.h>
#include <ItMethod.h>
#include <FixedPointIte.h>
#include <MultiGridScaIte.h>
#include <FgmresIte.h>

#include <MultiGrid2D.h>
#include <MGLevel2D.h>

#include <MainUtilities.h>
#include <TimeUtilities.h>

#define AMG 0
#define GMG 1
#define DIRECT 2


// =======================================================================
// include current example
// =======================================================================
// #include  "../Examples/TCD_2D/Levelset.h"
#include "../Examples/TNSE_2D/2PhaseLevelSet.h"


void Assemble2D_DG_LevelSet(TSquareMatrix2D *LS_sqmatrixK, TFEVectFunct2D *Velo)
{
  
  
  
  cout << "Assemble2D_DG_LevelSet " <<endl;
  exit(0);
}



void SetParametersCDAdapt2D()
{
   
  // check parameter consistency, set internal parameters
  if (TDatabase::ParamDB->SOLD_PARAMETER_TYPE==FEM_FCT)
   {
      TDatabase::ParamDB->INTERNAL_LINEAR_SCHEME = 0;
      TDatabase::ParamDB->DISCTYPE = 1; // GALERKIN
   }
      

  if (TDatabase::ParamDB->SOLD_PARAMETER_TYPE==FEM_FCT_LIN)
   {
    TDatabase::ParamDB->INTERNAL_LINEAR_SCHEME = TDatabase::ParamDB->FEM_FCT_LINEAR_TYPE;
    TDatabase::ParamDB->SOLD_PARAMETER_TYPE=FEM_FCT;
//      sold_parameter_type = TDatabase::ParamDB->SOLD_PARAMETER_TYPE;
    TDatabase::ParamDB->DISCTYPE = 1; // GALERKIN
   } 
  
  
  if (TDatabase::ParamDB->SOLD_PARAMETER_TYPE== KLR02_3)
    TDatabase::ParamDB->SOLD_S = 0;
  if (TDatabase::ParamDB->SOLD_PARAMETER_TYPE== LP96)
  {
    OutPut("SOLD_PARAMETER_TYPE == LP96 should be used with higher quadrature rule,"<<endl);
    OutPut("since right hand side is in general not linear !!!"<<endl);
  }

  if(TDatabase::ParamDB->DISCTYPE != LOCAL_PROJECTION)
  {
    // switch off all local projection terms
    TDatabase::ParamDB->LP_FULL_GRADIENT = 0;
    TDatabase::ParamDB->LP_FULL_GRADIENT_COEFF = 0;
    TDatabase::ParamDB->LP_FULL_GRADIENT_EXPONENT = 1;

    TDatabase::ParamDB->LP_STREAMLINE = 0;
    TDatabase::ParamDB->LP_STREAMLINE_COEFF = 0;
    TDatabase::ParamDB->LP_STREAMLINE_EXPONENT = 1;
  }

  if(TDatabase::ParamDB->LP_FULL_GRADIENT)
  {
    if(TDatabase::ParamDB->LP_STREAMLINE)
    {
      TDatabase::ParamDB->LP_STREAMLINE = 0;
      TDatabase::ParamDB->LP_STREAMLINE_COEFF = 0;
      TDatabase::ParamDB->LP_STREAMLINE_EXPONENT = 1;
      OutPut("local projection stabilisation in streamline direction ");
      OutPut("is switched off due to stabilisation of full gradient." << endl);
    }
  }

  if(TDatabase::ParamDB->LP_FULL_GRADIENT_ORDER_DIFFERENCE == -123)
    TDatabase::ParamDB->LP_FULL_GRADIENT_ORDER_DIFFERENCE = TDatabase::ParamDB->LP_ORDER_DIFFERENCE;

  if(TDatabase::ParamDB->LP_STREAMLINE_ORDER_DIFFERENCE == -123)
    TDatabase::ParamDB->LP_STREAMLINE_ORDER_DIFFERENCE = TDatabase::ParamDB->LP_ORDER_DIFFERENCE;
  

  if ((TDatabase::ParamDB->SDFEM_TYPE == 100)&&(!TDatabase::ParamDB->SOLVE_ADJOINT_PROBLEM)&&
      (TDatabase::ParamDB->DISCTYPE==SDFEM))
  {
      TDatabase::ParamDB->SDFEM_TYPE = 2;
      OutPut("Changed Database::ParamDB->SDFEM_TYPE to " << TDatabase::ParamDB->SDFEM_TYPE
	     << " since no adjoint problem is solved !!! "<<endl);
  } 
  if ((TDatabase::ParamDB->SDFEM_TYPE != 100)&&(TDatabase::ParamDB->SOLVE_ADJOINT_PROBLEM)&&
      (TDatabase::ParamDB->DISCTYPE==SDFEM))
  {
      TDatabase::ParamDB->SDFEM_TYPE = 100;
      OutPut("Changed Database::ParamDB->SDFEM_TYPE to " << TDatabase::ParamDB->SDFEM_TYPE
	     << " since adjoint problem is solved !!! "<<endl);
  } 
  if ((TDatabase::ParamDB->SOLVE_ADJOINT_PROBLEM==4))
    TDatabase::ParamDB->INTERNAL_NO_ESTIMATE_DIRICHLET_CELLS = 1;
  // SUPG 
  if ((TDatabase::ParamDB->DISCTYPE==SDFEM)&&(TDatabase::ParamDB->SOLD_TYPE==0))
  {
      // this excludes some not wished side effects
      TDatabase::ParamDB->SOLD_PARAMETER_TYPE = 0;
  } 
  if  (!(TDatabase::ParamDB->DISCTYPE==SDFEM))
    {
      TDatabase::ParamDB->SOLD_TYPE = 0;
      TDatabase::ParamDB->SOLD_PARAMETER_TYPE =0;
    }
  if  ((TDatabase::ParamDB->DISCTYPE==SDFEM)&&(TDatabase::ParamDB->SOLD_PARAMETER_TYPE == FEM_TVD))
    {
      TDatabase::ParamDB->SDFEM_TYPE = 0;
      TDatabase::ParamDB->DELTA0 =  TDatabase::ParamDB->DELTA1 = 0;
      OutPut("FEM-TVD: switched stabilization off!" << endl);
    }

  TDatabase::ParamDB->NSTYPE = 0;

  if (TDatabase::ParamDB->DISCTYPE==CIP)
    {
      TDatabase::ParamDB->DISCTYPE=GALERKIN;
      TDatabase::ParamDB->INTERNAL_FACE_INTEGRALS = 1;
    }
  if (TDatabase::ParamDB->DISCTYPE==DG)
    {
      TDatabase::ParamDB->DISCTYPE=GALERKIN;
      TDatabase::ParamDB->INTERNAL_FACE_INTEGRALS = 2;
      if ( TDatabase::ParamDB->ANSATZ_ORDER < 10)
	TDatabase::ParamDB->ANSATZ_ORDER = -TDatabase::ParamDB->ANSATZ_ORDER-10;
      else 
	// P elements on quads
	TDatabase::ParamDB->ANSATZ_ORDER = -10*TDatabase::ParamDB->ANSATZ_ORDER;
      if (TDatabase::ParamDB->ESTIMATE_ERRORS)
	{
	  TDatabase::ParamDB->ESTIMATE_ERRORS = 0;
	  OutPut("Error estimation does not work for DG !!!"<< endl);
	}
    }
}


// ======================================================================
// main program
// ======================================================================
int main(int argc, char* argv[])
{
  // ======================================================================
  // variable declaration
  // ======================================================================
  TDomain *Domain = new TDomain();
  TDomain *Domain_Intl = new TDomain();
  TDatabase *Database = new TDatabase();
  TFEDatabase2D *FEDatabase = new TFEDatabase2D();  
  TCollection *coll, *mortarcoll = NULL;  
  TOutput2D *Output;  
  TFEFunction2D *phi;
  TAuxParam2D *aux;
  TSquareStructure2D *LS_sqstructureA;
  TSquareMatrix2D *LS_sqmatrixA, *SQMATRICES_LS[3];
  TSquareMatrix2D *LS_sqmatrixM , *LS_sqmatrixK, *LS_sqmatrixS;
  TSquareMatrix2D *MatricesA, *MatricesM, *MatricesK, *MatricesS;
  TSquareMatrix **sqmatrices_ls = (TSquareMatrix **)SQMATRICES_LS;
  MatVecProc *MatVect;
  DefectProc *Defect;
  TFESpace2D *LevelSetSpaces;
  TDiscreteForm2D *DiscreteFormMatrixMRhs, *DiscreteFormMatrixMRhs_SUPG;
  TDiscreteForm2D *DiscreteFormMatrixARhs, *DiscreteFormMatrixARhs_SUPG;  
  BoundCondFunct2D *BoundaryConditions[2];
  BoundCondFunct2D *BDCond[4];
  BoundValueFunct2D *BoundValues[3];
  BoundValueFunct2D *BDValue[4];  
  CoeffFct2D *Coefficients[2];
  TFESpace2D *fesp[2], *ferhs[1];
  TDiscreteForm2D *DiscreteForm;  
  
  //for NSE 
  TFESpace2D *velocity_space, *pressure_space; 
  TFEVectFunct2D  *u;
  TFEFunction2D *u1, *u2, *p;
  
  int i, j, l, m, N_SubSteps, LEVELS, ret, Max_It, ORDER, N_LSDof, N_LSActive;
  int *N_Uarray, N_Cells, N_FESpaces, N_Rhs, VSP;
  int N_SquareMatrices, time_discs, img=1,  very_first_time=0;
  int *neum_to_diri, N_neum_to_diri = 0, *neum_to_diri_bdry, only_first_time;
  int velocity_space_code, pressure_space_code, N_Active, N_U, N_P, N_Unknowns;
  
  double total_time, limit, t3, oldtau, end_time, tau;
  double *RHSs[3], *LS_phi, *LS_phiold, *LS_rhs, *LS_defect, *LS_B, gamma;
  double *lump_mass, *oldrhs_fem_fct0, *oldrhs_fem_fct1, *matrix_D_Entries; 
  double *tilde_u, *neum_to_diri_param, residual;  
  double *Sol, *Rhs, *B, *defect;
  
  std::ostringstream os;
  char *PsBaseName, *GnuBaseName;
  char *VtkBaseName, *MatlabBaseName, *GmvBaseName;  
  char *PRM, *GEO;  
  char Readin[] = "readin.dat"; 
  char Name[] = "name";
  char NameStrings[5][10]= {"LS", "C1", "PSD", "C3", "C4"};  
  char UString[] = "u";
  char PString[] = "p";
  //======================================================================
  // read parameter file
  //======================================================================
  total_time = GetTime();
  if(argc>=2)
    ret=Domain->ReadParam(argv[1]);
  else
    ret=Domain->ReadParam(Readin);

  if(ret==-1)
  {
    exit(-1);
  }

  OpenFiles();
  OutFile.setf(std::ios::scientific);
  
  // write parameters into outfile
  Database->WriteParamDB(argv[0]);
  Database->WriteTimeDB();
  
  ExampleFile(); 
  Coefficients[0] = LinCoeffs;  
  Coefficients[1] = BilinearCoeffs_LS;
  
  BoundaryConditions[0] = BoundCondition;  
  BoundaryConditions[1] = BoundCondition_LS;
  
  BoundValues [0] = U1BoundValue;
  BoundValues [1] = U2BoundValue;   
  BoundValues [2] = BoundValue_LS; 
  
  
  // set all parameters  
  SetParametersCDAdapt2D();
   
  PRM = TDatabase::ParamDB->BNDFILE;
  GEO = TDatabase::ParamDB->GEOFILE;

  // assign names for output files
  PsBaseName = TDatabase::ParamDB->PSBASENAME;
  GnuBaseName = TDatabase::ParamDB->GNUBASENAME;
  GmvBaseName = TDatabase::ParamDB->GMVBASENAME;
  VtkBaseName = TDatabase::ParamDB->VTKBASENAME;
  MatlabBaseName = TDatabase::ParamDB->MATLABBASENAME;  
  
 //======================================================================
 // initialize discrete forms
 //======================================================================
 
  InitializeDiscreteFormsScalar(DiscreteFormMatrixMRhs, DiscreteFormMatrixARhs, DiscreteFormMatrixMRhs_SUPG,
                                 DiscreteFormMatrixARhs_SUPG, Coefficients[1]);  
  
   
  //======================================================================
  // read boundary parameterization and initialize coarse grid
  //======================================================================
  Domain->Init(PRM, GEO);
  
  
   // refine grid up to the coarsest level
  for(i=0;i<TDatabase::ParamDB->UNIFORM_STEPS;i++)
    Domain->RegRefineAll();  
   
   // write grid into an Postscript file
  os.seekp(std::ios::beg);
  os << "Domain" << ".ps" << ends;
  Domain->PS(os.str().c_str(),It_Finest,0);
  
  // initializ time
  TDatabase::TimeDB->CURRENTTIME = TDatabase::TimeDB->STARTTIME;
  TDatabase::TimeDB->CURRENTTIMESTEPLENGTH = TDatabase::TimeDB->TIMESTEPLENGTH;
  SetTimeDiscParameters();
  limit = TDatabase::ParamDB->SC_NONLIN_RES_NORM_MIN_SCALAR;
  Max_It = TDatabase::ParamDB->SC_NONLIN_MAXIT_SCALAR;

  t3 = GetTime();
  total_time = t3 - total_time;
  SetPolynomialDegree();  
  
  coll=Domain->GetCollection(It_Finest, 0);
  N_Cells = coll->GetN_Cells();
  OutPut("N_Cells (space) : " << N_Cells <<endl); 
  
  
  ORDER = 0;
  VSP = TDatabase::ParamDB->VELOCITY_SPACE;

  if(abs(VSP) > 20)
   {ORDER = abs(VSP) - 20;}
  else if (abs(VSP) > 10)
    {ORDER = abs(VSP) - 10;}
  else ORDER = abs(VSP); 

  //=========================================================================
  // construct all finite element spaces
  //=========================================================================
    // get velocity and pressure spacess
    GetVelocityAndPressureSpace(coll,BoundCondition,
                              mortarcoll, velocity_space,
                              pressure_space, &pressure_space_code,
                              TDatabase::ParamDB->VELOCITY_SPACE,
                              TDatabase::ParamDB->PRESSURE_SPACE);
    velocity_space_code =   TDatabase::ParamDB->VELOCITY_SPACE;
 
    N_Active =  velocity_space->GetActiveBound();
    N_U = velocity_space->GetN_DegreesOfFreedom();
    N_P = pressure_space->GetN_DegreesOfFreedom();  
 
    N_Unknowns = 2*N_U + N_P;
    OutPut("dof velocity : "<< setw(10) << 2* N_U << endl);
    OutPut("dof pressure : "<< setw(10) << N_P << endl);
    OutPut("dof all      : "<<  setw(10) << N_Unknowns  << endl);    
    
    //=========================================================================      
    // fespaces for scalar equations
    LevelSetSpaces =  new TFESpace2D(coll, Name, NameStrings[0], BoundaryConditions[1],
                                     TDatabase::ParamDB->ANSATZ_ORDER, NULL);
    
    N_LSDof = LevelSetSpaces->GetN_DegreesOfFreedom();
    N_LSActive = LevelSetSpaces->GetActiveBound();  
    
    OutPut("DOF Level set scalar : " << setw(10) << N_LSDof << endl);
 
    //=========================================================================
    // memory allocate all arrays and construction of all fefunction
    //=========================================================================
    Sol = new double[N_Unknowns];  
    Rhs = new double[N_Unknowns];    
    B = new double[N_Unknowns];
    defect = new double[N_Unknowns];    

    memset(Sol, 0, N_Unknowns*SizeOfDouble);
    memset(Rhs, 0, N_Unknowns*SizeOfDouble);    
    
    //velo vect
    u =  new TFEVectFunct2D(velocity_space, UString, UString, Sol, N_U, 2);
    u1 = u->GetComponent(0);
    u2 = u->GetComponent(1); 
    u1->Interpolate(InitialU1);
    u2->Interpolate(InitialU2);
  
    //pressure FeFunction
    p = new TFEFunction2D(pressure_space, PString,  PString,  Sol+2*N_U, N_P);    
    
    //=========================================================================     
    LS_phi =  new double[N_LSDof]; 
    LS_phiold =  new double[N_LSDof];     
    LS_rhs =  new double[N_LSDof]; 
    LS_defect = new double[N_LSDof]; 
    LS_B = new double[N_LSDof];     
    
    memset(LS_phi, 0, N_LSDof*SizeOfDouble);   
    memset(LS_phiold, 0, N_LSDof*SizeOfDouble);       
    memset(LS_rhs, 0, N_LSDof*SizeOfDouble);   
    memset(LS_defect, 0, N_LSDof*SizeOfDouble);   
    
    phi = new TFEFunction2D(LevelSetSpaces, NameStrings[0], NameStrings[0], LS_phi, N_LSDof);  
    phi->Interpolate(InitialCondition_LS);
       
    //=========================================================================
    // allocate memory for all matrices
    //=========================================================================
    // build matrices
    // first build matrix structure
    LS_sqstructureA = new TSquareStructure2D(LevelSetSpaces);
    LS_sqstructureA->Sort();        

    // two matrices used, M is the mass matrix, as well as the system matrix
    LS_sqmatrixM = new TSquareMatrix2D(LS_sqstructureA);

    // A contains the non time dependent part of the discretization
    LS_sqmatrixA = new TSquareMatrix2D(LS_sqstructureA);

    if(TDatabase::ParamDB->DISCTYPE == SDFEM || TDatabase::ParamDB->INTERNAL_FACE_INTEGRALS>0)
     {
      LS_sqmatrixK = new TSquareMatrix2D(LS_sqstructureA);      // stabilisation matrix K
     }
     
    if(TDatabase::ParamDB->SOLD_PARAMETER_TYPE == FEM_FCT)
     {       
      lump_mass = new double [N_LSDof];
      oldrhs_fem_fct0 = new double [N_LSDof];     
      tilde_u = new double [N_LSDof];         
      oldrhs_fem_fct1 = new double [N_LSDof];
      matrix_D_Entries = new double[LS_sqmatrixA->GetN_Entries()];
      LS_sqmatrixK = new TSquareMatrix2D(LS_sqstructureA);  
     }       

    // pointers to the routines which compute matrix-vector
    // products and the defect
    MatVect = MatVect_Scalar;
    Defect = Defect_Scalar;     
    //=========================================================================
    //  prepare output
    //=========================================================================
    Output = new TOutput2D(2, 1, 1, 1, Domain);

 
    Output->AddFEVectFunct(u);
    Output->AddFEFunction(phi);

    os.seekp(std::ios::beg);
    Output->AddParameter(TDatabase::TimeDB->CURRENTTIME,os.str().c_str());     

  //======================================================================
  // assembling of mass matrix and initial rhs (f_0)
  //======================================================================

  // set parameters
    N_Rhs = 1;
    N_FESpaces = 1;
    fesp[0] = LevelSetSpaces;    
    
    aux =  new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);       
    
    // reset matrices
    N_SquareMatrices = 1;
    SQMATRICES_LS[0] = LS_sqmatrixM;
    SQMATRICES_LS[0]->Reset();

    DiscreteForm = DiscreteFormMatrixMRhs;

    BDCond[0] = BoundaryConditions[1];
    BDValue[0] = BoundValues[2];
    
    if(TDatabase::ParamDB->SOLD_PARAMETER_TYPE==FEM_FCT)
     BoundValues[2] = BoundValue_FEM_FCT;
    
    memset(LS_rhs, 0, N_LSDof*SizeOfDouble);
    RHSs[0] = LS_rhs;
    ferhs[0] = LevelSetSpaces;    
    
        
    Assemble2D(N_FESpaces, fesp,
               N_SquareMatrices, SQMATRICES_LS,
               0, NULL,
               N_Rhs, RHSs, ferhs,
               DiscreteForm,
               BDCond,
               BDValue,
               aux);

    delete aux;
     
   // copy Dirichlet values from rhs into sol
   //    memcpy(LS_phi+N_LSActive, RHSs[0]+N_LSActive, (N_LSDof-N_LSActive)*SizeOfDouble);     
     
   if(TDatabase::ParamDB->SOLD_PARAMETER_TYPE == FEM_FCT)
    {        
     BoundValues[2] = BoundValue_LS;      
      
     LumpMassMatrixToVector(LS_sqmatrixM, lump_mass);    
//       save mass matrix in LS_sqmatrixK
      memcpy(LS_sqmatrixK->GetEntries(), LS_sqmatrixM->GetEntries(),
             LS_sqmatrixM->GetN_Entries() * SizeOfDouble);   
             
      memcpy(oldrhs_fem_fct0, LS_rhs, N_LSDof*SizeOfDouble);   
      
#ifdef __ROTATING_BODIES__
      CheckWrongNeumannNodes(coll, LevelSetSpaces, N_neum_to_diri, neum_to_diri,
                             neum_to_diri_bdry, neum_to_diri_param);  
#endif      
           
#ifdef __LEVELSET__       
//       CheckWrongNeumannNodes(coll, LevelSetSpaces, N_neum_to_diri, neum_to_diri,
//                              neum_to_diri_bdry, neum_to_diri_param);      
#endif               
    }
    
   // save solution
   memcpy(LS_phiold, LS_phi, N_LSDof*SizeOfDouble);   
   
  //======================================================================     
  // parameters for time stepping scheme
  //======================================================================      
  gamma = 0;
  m = 0;
  N_SubSteps = GetN_SubSteps();
  oldtau = 1;
  end_time = TDatabase::TimeDB->ENDTIME;     
  
  
  // not active : TDatabase::TimeDB->TIMESTEPLENGTH_CONTROL = 0
  if (TDatabase::TimeDB->TIMESTEPLENGTH_CONTROL)
    { time_discs = 2; }
  else
    { time_discs = 1; }
  
  //======================================================================     
  // output the solutin
  //======================================================================      
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

  //======================================================================
  // start of time cycle
  //======================================================================
  tau = TDatabase::TimeDB->CURRENTTIME;
  while(TDatabase::TimeDB->CURRENTTIME< end_time)
  {                                               // time cycle
    m++;
    TDatabase::TimeDB->INTERNAL_STARTTIME = TDatabase::TimeDB->CURRENTTIME;

    for(l=0;l<N_SubSteps;l++)                     // sub steps of fractional step theta
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
      
      
      // working array for rhs is B, initialize B
       memset(LS_B, 0, N_LSDof*SizeOfDouble);     
      
      // compute terms with data from previous time step
      // old rhs multiplied with current subtime step and theta3 on B
       Daxpy(N_LSActive, tau*TDatabase::TimeDB->THETA3, LS_rhs, LS_B);


       // assemble A and rhs
       N_Rhs = 1;
       N_FESpaces = 1;
       fesp[0] = LevelSetSpaces;    
       
       aux =  new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);            

       if(TDatabase::ParamDB->DISCTYPE == SDFEM)
        {
              N_SquareMatrices = 2;
              SQMATRICES_LS[0] = LS_sqmatrixA;
              SQMATRICES_LS[0]->Reset();
              SQMATRICES_LS[1] = LS_sqmatrixK;
              SQMATRICES_LS[1]->Reset();
              DiscreteForm = DiscreteFormMatrixARhs_SUPG;
        }
       else
        {
              N_SquareMatrices = 1;
              SQMATRICES_LS[0] = LS_sqmatrixA;
              SQMATRICES_LS[0]->Reset();
              DiscreteForm = DiscreteFormMatrixARhs;
        }

        BDCond[0] = BoundaryConditions[1];
        BDValue[0] = BoundValues[2];
       
        if(TDatabase::ParamDB->SOLD_PARAMETER_TYPE==FEM_FCT)
          BoundValues[2] = BoundValue_FEM_FCT;
 
        memset(LS_rhs, 0, N_LSDof*SizeOfDouble);
        RHSs[0] = LS_rhs;
        ferhs[0] = LevelSetSpaces;    
    
         Assemble2D(N_FESpaces, fesp,
              N_SquareMatrices, SQMATRICES_LS,
              0, NULL,
              N_Rhs, RHSs, ferhs,
              DiscreteForm,
              BDCond,
              BDValue,
              aux);
    
            delete aux;       

         if(TDatabase::ParamDB->INTERNAL_FACE_INTEGRALS>0)
           Assemble2D_DG_LevelSet(LS_sqmatrixK, u);
      
   
         // save rhs without Dirichlet values
         if(TDatabase::ParamDB->SOLD_PARAMETER_TYPE==FEM_FCT)
          {
           BoundValues[2] = BoundValue_LS;
           memcpy(oldrhs_fem_fct1, LS_rhs, N_LSDof*SizeOfDouble);
          }

        if (very_first_time==1)
        {
          very_first_time=0;
          l--;
          continue;
        }           

        // add rhs from current sub time step to rhs array B     
        Daxpy(N_LSActive, tau*TDatabase::TimeDB->THETA4, LS_rhs, LS_B);      

        oldtau = tau;
           
        
        if(!(TDatabase::ParamDB->SOLD_PARAMETER_TYPE== FEM_FCT))
        { 

          if(TDatabase::ParamDB->DISCTYPE == SDFEM)
            MatAdd(LS_sqmatrixA, LS_sqmatrixK, 1.);

          if(TDatabase::ParamDB->INTERNAL_FACE_INTEGRALS>0)
            MatAdd(LS_sqmatrixA, LS_sqmatrixK, 1.);

          // update rhs by Laplacian and convective term from previous time step
          // scaled by current sub time step length and theta2
          // currently : M := M + gamma A
          // M = M + (- tau*TDatabase::TimeDB->THETA2)
          MatAdd(LS_sqmatrixM, LS_sqmatrixA, - tau*TDatabase::TimeDB->THETA2);
          // set current factor of steady state matrix
          gamma = -tau*TDatabase::TimeDB->THETA2;
    
          // defect = M * sol
          // B:= B + defec
          memset(LS_defect, 0, N_LSDof*SizeOfDouble);
          MatVectActive(LS_sqmatrixM,  LS_phi, LS_defect);
          Daxpy(N_LSActive, 1, LS_defect, LS_B);

          // set Dirichlet values
          // RHSs[0] still available from assembling
          memcpy(LS_B+N_LSActive, RHSs[0]+N_LSActive, (N_LSDof-N_LSActive)*SizeOfDouble);             
          // copy Dirichlet values from rhs into sol
          memcpy(LS_phi+N_LSActive, RHSs[0]+N_LSActive, (N_LSDof-N_LSActive)*SizeOfDouble);    

          // system matrix
          MatAdd(LS_sqmatrixM,  LS_sqmatrixA, -gamma + tau*TDatabase::TimeDB->THETA1);
          gamma = tau*TDatabase::TimeDB->THETA1;
        }
        else
        {
          only_first_time = 1;
 
          // FEM-FCT methods
          FEM_FCT_ForConvDiff(LS_sqmatrixK, LS_sqmatrixA,
                              N_LSDof, N_LSActive, 
                              lump_mass, matrix_D_Entries, 
                              LS_phi, LS_phiold, 
                              LS_B, LS_rhs, oldrhs_fem_fct0, tilde_u,
                              N_neum_to_diri, neum_to_diri, 
                              neum_to_diri_bdry, neum_to_diri_param,
                              only_first_time,
                              BoundValue_LS,NULL); 

          SQMATRICES_LS[0] = LS_sqmatrixM;

          LS_sqmatrixM->Reset();
          // system matrix for FEM-FCT   M_lump + theta1*tau*A
          // A = Galerkin + D
          FEM_FCT_SystemMatrix(LS_sqmatrixM, LS_sqmatrixA, lump_mass, N_LSDof);
        }

//======================================================================
// solve linear system
//======================================================================
        switch(int(TDatabase::ParamDB->SOLVER_TYPE))
         {
          case 2:
            DirectSolver(LS_sqmatrixM, LS_B, LS_phi);
          break;

          default:
            cout << "wrong solver type !!!!!!!!!!!!!" << endl;
            exit(0);
            break;
         }

        //======================================================================
        // end solve linear system
        //======================================================================
        
	// reset matrices for SUPG and SOLD
        if (!(TDatabase::ParamDB->SOLD_PARAMETER_TYPE == FEM_FCT))
         {        
          // restore matrices
          MatAdd(LS_sqmatrixM, LS_sqmatrixA, -gamma);
          gamma = 0;

          if(TDatabase::ParamDB->DISCTYPE == SDFEM)
            MatAdd(LS_sqmatrixA, LS_sqmatrixK, -1.);
         }
         
         
      memcpy(LS_phiold, LS_phi, N_LSDof*SizeOfDouble);          
         
    } // for(l=0;l<N_SubSteps;l++)    
    
    
     if(m==1 || m % TDatabase::TimeDB->STEPS_PER_IMAGE == 0)
      {
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
    } //   if(m==1 || m % TDatabase::TimeDB->STEPS_PER_IMA

  } // while(TDatabase::TimeDB->CURRENTTIME
    
  CloseFiles();
  return 0;
}



  
  





