// =======================================================================
//
// Purpose:     main program for solving a time-dependent NSE equation in ParMooN
//
// Author:      Sashikumaar Ganesan
//
// History:     Implementation started on 03.09.2014

// =======================================================================


#include <Domain.h>
#include <Database.h>
#include <SystemTCD2D.h>
#include <SystemTCD2D_ALE.h>
#include <FEDatabase2D.h>
#include <LinAlg.h>
#include <FESpace2D.h>
#include <SystemTNSE2D.h>
#include <SquareStructure2D.h>
#include <Structure2D.h>
#include <Output2D.h>
#include <CD2DErrorEstimator.h>
#include <MainUtilities.h>
// #include <TimeUtilities.h>
#include <TNSE2D_ParamRout.h>
#include <TimeDiscRout.h>
#include <TimeConvDiff2D.h>

#include <string.h>
#include <sstream>
#include <MooNMD_Io.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>
// =======================================================================
// include current example
// =======================================================================
// #include "../Examples/TNSE_2D/DrivenCavity.h" //   in unit square
// #include "../Examples/TNSE_2D/Bsp1.h" // smooth sol in unit square
// #include "../Examples_All/TNSE_2D/Benchmark2.h"  
// #include "../Examples/TNSE_2D/SinCos.h" // smooth sol in unit square

// #include "sqcylin.h"
// =======================================================================
// main program
// =======================================================================

# include "MHD2D.h"

int main(int argc, char* argv[])
{
  // ======================================================================
  //  declaration of variables
  // ======================================================================
  int i, j, l, m, N_Cells,img=1 ;  
  double end_time, oldtau, tau, N_SubSteps;                               //Common Variables
  TDomain *Domain;
  TDatabase *Database = new TDatabase(); 
  TFEDatabase2D *FEDatabase = new TFEDatabase2D(); 
  TCollection *coll, *mortarcoll = NULL;
  TOutput2D *Output;
  MultiIndex2D AllDerivatives[3] = { D00, D10, D01 };
  
  double limit, AllErrors[7];                                             // TNSE2D 
  int N_U, N_P, N_L, N_TotalDOF, pressure_space_code;
  int Max_It, NSEType, velocity_space_code, Disctype;  
  double *sol, *rhs, *oldrhs, *defect, residual, impuls_residual;  
  TFESpace2D *Velocity_FeSpace, *Pressure_FeSpace, *fesp[2];
  TFEVectFunct2D *Velocity;
  TFEFunction2D *u1, *u2, *Pressure, *fefct[2];
  TSystemTNSE2D *SystemMatrix;
  TAuxParam2D *aux, *NSEaux_error;

  int N_SubSteps_cd, ORDER, N_DOF, N_G;
  int N_Active; 
  double *sol_cd, *rhs_cd, *oldrhs_cd, errors[5], Linfty;                 // TCD2D    
  TFESpace2D *Temperature_FeSpace, *fesp_cd[1];                  
  TFEFunction2D *Temperature;
  TSystemTCD2D *SystemMatrix_cd;  
  TAuxParam2D *aux_cd;
  bool UpdateStiffnessMat, UpdateRhs,  ConvectionFirstTime;  
   
#ifdef __PRIVATE__ 
  TFESpace2D *Projection_space;
#endif

  const char vtkdir[] = "VTK"; 
  char *PsBaseName, *VtkBaseName, *GEO;
  char UString[] = "u";
  char PString[] = "p";
  char TString[] = "Theta";
  char NameString[] = "VMS";
  
  std::ostringstream os;
  os << " ";   

  if(TDatabase::ParamDB->WRITE_VTK)
      mkdir(vtkdir, 0777);
    
  // ======================================================================
  // set the database values and generate mesh
  // ======================================================================    
  /** set variables' value in TDatabase using argv[1] (*.dat file), and generate the MESH based */
  Domain = new TDomain(argv[1]);  

  if(TDatabase::ParamDB->PROBLEM_TYPE == 0)         // CD 
    TDatabase::ParamDB->PROBLEM_TYPE = 2;
  
  OpenFiles();
  OutFile.setf(std::ios::scientific);

  Database->CheckParameterConsistencyNSE();
  Database->WriteParamDB(argv[0]);
  Database->WriteTimeDB();
  ExampleFile();

  GEO = TDatabase::ParamDB->GEOFILE;
  Domain->Init(NULL, GEO);  

  TriaReMeshGen(Domain);
  TDatabase::ParamDB->UNIFORM_STEPS = 0;  

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

//=========================================================================
// CELL CONSTRUCTION
//=========================================================================   
  coll=Domain->GetCollection(It_Finest, 0);
  N_Cells = coll->GetN_Cells();
  OutPut("N_Cells : " << N_Cells <<endl);
  
//=========================================================================
// construct all TNSE2D finite element spaces 
//=========================================================================
  NSEType = TDatabase::ParamDB->NSTYPE;           
  Disctype = TDatabase::ParamDB->DISCTYPE;
  
// fespaces for velocity and pressure  
  GetVelocityAndPressureSpace(coll, BoundCondition, mortarcoll, Velocity_FeSpace,
                              Pressure_FeSpace, &pressure_space_code,
                              TDatabase::ParamDB->VELOCITY_SPACE,
                              TDatabase::ParamDB->PRESSURE_SPACE);
  
  
  // defaulty inf-sup pressure space will be selected based on the velocity space, so update it in database
  TDatabase::ParamDB->INTERNAL_PRESSURE_SPACE = pressure_space_code;
  velocity_space_code =   TDatabase::ParamDB->VELOCITY_SPACE;  
  
  N_U = Velocity_FeSpace->GetN_DegreesOfFreedom();
  N_P = Pressure_FeSpace->GetN_DegreesOfFreedom();    
  N_TotalDOF = 2*N_U + N_P;                         // Size of Matrix

//=========================================================================
// construct all TCD2D finite element spaces
//=========================================================================  
  ORDER = TDatabase::ParamDB->ANSATZ_ORDER;
  Temperature_FeSpace = new TFESpace2D(coll, (char*)"fe space", (char*)"solution space", 
                                          BoundCondition_CD,ContP_USpace,ORDER, mortarcoll);

  N_DOF = Temperature_FeSpace->GetN_DegreesOfFreedom();
  N_Active =  Temperature_FeSpace->GetActiveBound();
  
  OutPut("Dof Velocity : "<< setw(10) << 2* N_U << endl);
  OutPut("Dof Pressure : "<< setw(10) << N_P << endl);
  OutPut("DOF THETA    : "<< setw(10) << N_DOF  << endl);
  OutPut("Total Dof all: "<< setw(10) << (N_TotalDOF+N_DOF)  << endl); 
  OutPut("dof active   : "<< setw(10) << N_Active << endl);    

//========================================================================
//Construct Solution Space
//========================================================================
  sol = new double[N_TotalDOF+N_DOF]; 
  rhs = new double[N_TotalDOF+N_DOF];
  oldrhs = new double[N_TotalDOF+N_DOF]; 
    
//======================================================================
// construct all TNSE2D finite element functions
//======================================================================

  
  memset(sol, 0, N_TotalDOF*SizeOfDouble);
  memset(rhs, 0, N_TotalDOF*SizeOfDouble);  
  
  Velocity = new TFEVectFunct2D(Velocity_FeSpace, UString,  UString,  sol, N_U, 2);
  u1 = Velocity->GetComponent(0);
  u2 = Velocity->GetComponent(1);
  Pressure = new TFEFunction2D(Pressure_FeSpace, PString,  PString,  sol+2*N_U, N_P);
  
  //  interpolate the initial solution
    u1->Interpolate(InitialU1);
    u2->Interpolate(InitialU2);
    Pressure->Interpolate(InitialP);

//======================================================================
// construct all TCD2D finite element functions
//======================================================================
  sol_cd = (sol + N_TotalDOF);
  rhs_cd = (rhs + N_TotalDOF);
  oldrhs_cd = &oldrhs[N_TotalDOF];
    
  memset(sol_cd, 0, N_DOF*SizeOfDouble);
  memset(rhs_cd, 0, N_DOF*SizeOfDouble);

  Temperature = new TFEFunction2D(Temperature_FeSpace, TString, TString, sol_cd, N_DOF); 
  
//======================================================================
// TCD2D SystemMatrix construction and solution
//======================================================================    
  
  
   if (TDatabase::ParamDB->INTERNAL_CONVECTION_EQ_VELOFIELD)
    {
      fesp[0] = Velocity_FeSpace;
      fesp[1] = Temperature_FeSpace;
      
      fefct[0] = u1;
      fefct[1] = u2;

      aux_cd =  new TAuxParam2D(TimeCDParamsVeloFieldN_FESpaces,
        TimeCDParamsVeloFieldN_Fct,
        TimeCDParamsVeloFieldN_ParamFct,
        TimeCDParamsVeloFieldN_FEValues,
        fesp+1, fefct,
        TimeCDParamsVeloFieldFct,
        TimeCDParamsVeloFieldFEFctIndex,
        TimeCDParamsVeloFieldFEMultiIndex,
        TimeCDParamsVeloFieldN_Params,
        TimeCDParamsVeloFieldBeginParam);
    }
    else
    {
      aux_cd =  new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);
    }
  
//          cout<< "test Main " <<endl;
//        exit(0);  
  
    // Disc type: GALERKIN (or) SDFEM  (or) UPWIND (or) SUPG (or) LOCAL_PROJECTION
    // Solver: AMG_SOLVE (or) GMG  (or) DIRECT 
    SystemMatrix_cd = new TSystemTCD2D(Temperature_FeSpace, GALERKIN, DIRECT);
//              cout<< "test Main " <<endl;
//        exit(0);
    // initilize the system matrix with the functions defined in Example file
    SystemMatrix_cd->Init(BilinearCoeffs, BoundCondition_CD, BoundValue_CD);
//                 cout<< "test Main " <<endl;
//        exit(0);
    // assemble the system matrix with given aux, sol and rhs 
    // aux is used to pass  addition fe functions (eg. mesh velocity) that is nedded for assembling,
    // otherwise, just pass with NULL 
    SystemMatrix_cd->AssembleMRhs(aux_cd, sol_cd, rhs_cd);   
  
//            cout<< "test Main " <<endl;
//        exit(0);
//======================================================================
// TNSE SystemMatrix construction and solution
//======================================================================  
    // Disc type: GALERKIN 
    // Solver: AMG_SOLVE (or) GMG  (or) DIRECT  
   SystemMatrix = new TSystemTNSE2D(Velocity_FeSpace, Pressure_FeSpace, Velocity, Pressure, sol, rhs, Disctype, NSEType, DIRECT
#ifdef __PRIVATE__  
                   , Projection_space, NULL, NULL
#endif     
    );
  
    //define the aux
    fesp[0] = Velocity_FeSpace;

    fefct[0] = u1;
    fefct[1] = u2;  
    
    switch(Disctype)
     {
      // turbulent viscosity must be computed
      case SMAGORINSKY:
      case VMS_PROJECTION:
      case CLASSICAL_LES:
      case GL00_CONVOLUTION:
      case GL00_AUX_PROBLEM:

      aux =  new TAuxParam2D(TimeNSN_FESpacesVelo_GradVelo, TimeNSN_FctVelo_GradVelo,
                        TimeNSN_ParamFctVelo_GradVelo,
                        TimeNSN_FEValuesVelo_GradVelo,
                        fesp, fefct,
                        TimeNSFctVelo_GradVelo,
                        TimeNSFEFctIndexVelo_GradVelo,
                        TimeNSFEMultiIndexVelo_GradVelo,
                        TimeNSN_ParamsVelo_GradVelo,
                        TimeNSBeginParamVelo_GradVelo);

        break;

      default:
        // 2 parameters are needed for assembling (u1_old, u2_old)
        aux =  new TAuxParam2D(TimeNSN_FESpaces2, TimeNSN_Fct2, TimeNSN_ParamFct2,
                           TimeNSN_FEValues2,
                           fesp, fefct,
                           TimeNSFct2,
                           TimeNSFEFctIndex2, TimeNSFEMultiIndex2,
                           TimeNSN_Params2, TimeNSBeginParam2);  
    }
    
   // aux for calculating the error
   if(TDatabase::ParamDB->MEASURE_ERRORS)
    {
       NSEaux_error =  new TAuxParam2D(TimeNSN_FESpaces2, TimeNSN_Fct2,
                             TimeNSN_ParamFct2,
                             TimeNSN_FEValues2,
                             fesp, fefct,
                             TimeNSFct2,
                             TimeNSFEFctIndex2, TimeNSFEMultiIndex2,
                             TimeNSN_Params2, TimeNSBeginParam2);     
    }

//      TDatabase::TimeDB->TIMESTEPLENGTH = hmin*hmin;
//       cout<<TDatabase::TimeDB->TIMESTEPLENGTH<<"\n"; 
   
  //======================================================================  

    // initilize the system matrix with the functions defined in Example file
    // last argument is aux that is used to pass additional fe functions (eg. mesh velocity)     
   
    SystemMatrix->Init(LinCoeffs, BoundCondition, U1BoundValue, U2BoundValue, aux, NSEaux_error);

    // assemble M, A matrices and rhs 
    SystemMatrix->Assemble(sol, rhs);
    
    
//======================================================================
// produce outout at t=0
//======================================================================
   VtkBaseName = TDatabase::ParamDB->VTKBASENAME;    
   Output = new TOutput2D(4, 2, 1, 1, Domain);

   Output->AddFEVectFunct(Velocity);
   Output->AddFEFunction(Pressure);
   Output->AddFEFunction(Temperature);
   
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
 
    double hmin, hmax;
    coll->GetHminHmax(&hmin,&hmax);
    OutPut("h_min : " << hmin << " h_max : " << hmax << endl);
 
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
   
   UpdateStiffnessMat = FALSE; //check BilinearCoeffs in example file
   UpdateRhs = TRUE; //check BilinearCoeffs in example file
   ConvectionFirstTime=TRUE;   
// for(int ik = 0 ; ik<N_DOF ; ik++)
//       cout << ik << " Theta " <<sol[ik+2*N_U+N_P] << endl; 
   while(TDatabase::TimeDB->CURRENTTIME< end_time)
   {                                               // time cycle
     m++;
     TDatabase::TimeDB->INTERNAL_STARTTIME = TDatabase::TimeDB->CURRENTTIME;
     cout << "------ INTERNAL_STARTTIME --------- " << TDatabase::TimeDB->INTERNAL_STARTTIME << " " << N_SubSteps <<endl;

     for(l=0;l<N_SubSteps;l++) // sub steps of fractional step theta
     {
        SetTimeDiscParameters(1);
        
        if(m==1)
        {
            OutPut("Theta1: " << TDatabase::TimeDB->THETA1<< endl);
            OutPut("Theta2: " << TDatabase::TimeDB->THETA2<< endl);
            OutPut("Theta3: " << TDatabase::TimeDB->THETA3<< endl);
            OutPut("Theta4: " << TDatabase::TimeDB->THETA4<< endl);
        }
       
        tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
        TDatabase::TimeDB->CURRENTTIME += tau;
   
        OutPut(endl << "CURRENT TIME: ");
        OutPut(TDatabase::TimeDB->CURRENTTIME << endl);          
       
        //copy sol, rhs to olssol, oldrhs
        memcpy(oldrhs, rhs, N_TotalDOF*SizeOfDouble);        
    
        // assemble only rhs, nonlinear matrix for NSE will be assemble in fixed point iteration
        // not needed if rhs is not time-dependent
        if(m!=1)
            { SystemMatrix->AssembleRhs(sol, rhs); }
        else
            { SystemMatrix->Assemble(sol, rhs);  }
       
        //scale B matices and assemble NSE-rhs based on the \theta time stepping scheme 
        SystemMatrix->AssembleSystMat(tau/oldtau, oldrhs, rhs, sol); 
        oldtau = tau;
         
        // calculate the residual
        defect = new double[N_TotalDOF];
        memset(defect,0,N_TotalDOF*SizeOfDouble);

        SystemMatrix->GetTNSEResidual(sol, defect);

        //correction due to L^2_O Pressure space 
        if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
            IntoL20Vector2D(defect+2*N_U, N_P, pressure_space_code);
    
        residual =  Ddot(N_TotalDOF, defect, defect);
        impuls_residual = Ddot(2*N_U, defect, defect);  
        OutPut("Nonlinear iteration step   0");
        OutPut(setw(14) << impuls_residual);
        OutPut(setw(14) << residual-impuls_residual);
        OutPut(setw(14) << sqrt(residual) << endl);       

        //====================================================================== 
        //Solve the system
        //Nonlinear iteration of fixed point type
        //======================================================================
        for(j=1;j<=Max_It;j++)
        {   
            // Solve the NSE system
            SystemMatrix->Solve(sol);
      
            if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
                IntoL20FEFunction(sol+2*N_U, N_P, Pressure_FeSpace, velocity_space_code, pressure_space_code);   

            //no nonlinear iteration for Stokes problem  
            if(TDatabase::ParamDB->FLOW_PROBLEM_TYPE==STOKES)
                break;
     
            // restore the mass matrix for the next nonlinear iteration      
            SystemMatrix->RestoreMassMat();     
    
            // assemble the system matrix with given aux, sol and rhs 
            SystemMatrix->AssembleANonLinear(sol, rhs);    
          
            // assemble system mat, S = M + dt\theta_1*A
            SystemMatrix->AssembleSystMatNonLinear();        

            // get the residual
            memset(defect, 0, N_TotalDOF*SizeOfDouble);
            SystemMatrix->GetTNSEResidual(sol, defect);       
             
            //correction due to L^2_O Pressure space 
            if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
                IntoL20Vector2D(defect+2*N_U, N_P, pressure_space_code);
    
            residual =  Ddot(N_TotalDOF, defect, defect);
            impuls_residual = Ddot(2*N_U, defect, defect); 
            OutPut("nonlinear iteration step " << setw(3) << j);
            OutPut(setw(14) << impuls_residual);
            OutPut(setw(14) << residual-impuls_residual);
            OutPut(setw(14) << sqrt(residual) << endl); 
	
            if(sqrt(residual)<=limit)
            break;
       }
       
       // ====================================================================
       // Integrating TCD2D Equation
       // ====================================================================
       
//       cout<< "test Main " <<endl;
 //      exit(0);
       
        //copy rhs to oldrhs
        memcpy(oldrhs_cd, rhs_cd, N_DOF*SizeOfDouble); 
        // unless the stiffness matrix or rhs change in time, it is enough to 
        // assemble only once at the begning
        if(UpdateStiffnessMat || UpdateRhs ||  ConvectionFirstTime)
        {
            SystemMatrix_cd->AssembleARhs(aux_cd, sol_cd, rhs_cd);

            // M:= M + (tau*THETA1)*A
            // rhs: =(tau*THETA4)*rhs +(tau*THETA3)*oldrhs +[M-(tau*THETA2)A]*oldsol
            // note! sol contains only the previous time step value, so just pass 
            // sol for oldsol
            SystemMatrix_cd->AssembleSystMat(oldrhs_cd, sol_cd, rhs_cd, sol_cd);
            ConvectionFirstTime = FALSE;
        }
        // solve the system matrix 
        SystemMatrix_cd->Solve(sol_cd, rhs_cd);
        
        // restore the mass matrix for the next time step    
        // unless the stiffness matrix or rhs change in time, it is not necessary to assemble the system matrix in every time step
        if(UpdateStiffnessMat || UpdateRhs)
        {
            SystemMatrix_cd->RestoreMassMat();
        }          
        
       // restore the mass matrix for the next time step          
       SystemMatrix->RestoreMassMat();     
       
      } // for(l=0;l<N_SubSteps;
//====================================================================== 
// measure errors to known solution
//======================================================================    
/*    if(TDatabase::ParamDB->MEASURE_ERRORS)
    {   
//        u1->Interpolate(ExactU1);
//        u2->Interpolate(ExactU2);
//        Pressure->Interpolate(ExactP); 
       
      SystemMatrix->MeasureTNSEErrors(ExactU1, ExactU2, ExactP, AllErrors);

       OutPut("L2(u): " <<   AllErrors[0] << endl);
       OutPut("H1-semi(u): " <<  AllErrors[1] << endl);
       OutPut("L2(p): " <<  AllErrors[2] << endl);
       OutPut("H1-semi(p): " <<  AllErrors[3]   << endl); 
       OutPut(AllErrors[4] <<  " l_infty(L2(u)) " <<AllErrors[5] << endl);
       OutPut(TDatabase::TimeDB->CURRENTTIME << " L2(0,t,L2)(u) : " <<   sqrt(AllErrors[6]) << endl); 
      
     } // if(TDatabase::ParamDB->MEASURE_ERRORS)

    if(TDatabase::ParamDB->MEASURE_ERRORS)
    {
      Temperature->GetErrors(Exact, 3, AllDerivatives, 2, L2H1Errors, BilinearCoeffs, aux, 1, fesp, errors);


      OutPut("time: " << TDatabase::TimeDB->CURRENTTIME);
      OutPut(" L2: " << errors[0]);
      OutPut(" H1-semi: " << errors[1] << endl);

      errors[3] += (errors[0]*errors[0] + olderror * olderror)*TDatabase::TimeDB->TIMESTEPLENGTH/2.0;
      olderror = errors[0];
      OutPut(TDatabase::TimeDB->CURRENTTIME <<  " L2(0,T;L2) " << sqrt(errors[3]) << " ");

      errors[4] += (errors[1]*errors[1] +olderror1 * olderror1)*TDatabase::TimeDB->TIMESTEPLENGTH/2.0;
      OutPut( "L2(0,T;H1) " << sqrt(errors[4]) << endl);
      olderror1 = errors[1];   
      
      
      if(Linfty<errors[0])
      Linfty=errors[0];

      OutPut( "Linfty " << Linfty << endl);      
    } //  if(TDatabase::ParamDB->MEASURE_ERRORS)  */     
//======================================================================
// produce outout
//====================================================================== 
     if(m==1 || m % TDatabase::TimeDB->STEPS_PER_IMAGE == 0)      
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
                
    } // while(TDatabase::TimeDB->CURRENTTIME< e
 
 //======================================================================
// produce final outout
//======================================================================
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

  CloseFiles();
  
//   for(int ik = 0 ; ik<N_U ; ik++)
//       cout << ik << " U1 " << sol[ik] << endl;
//   
//   for(int ik = 0 ; ik<N_U ; ik++)
//       cout << ik << " U2 " << sol[ik+N_U] << endl;
// 
//   for(int ik = 0 ; ik<N_P ; ik++)
//       cout << ik << " P " << sol[ik+2*N_U] << endl; 
// 
//     for(int ik = 0 ; ik<N_DOF ; ik++)
//         cout << ik << " Theta " << sol[ik+2*N_U+N_P] << endl; 
//    
//   for(int ik = 0 ; ik<N_DOF ; ik++)
//       cout << ik << " sol_cd " << sol_cd[ik] << endl;     
//   
//      double x1,x2 , y1,y2 ;
/*  
     for(int ik = 0 ; ik < N_U ; ik++)
     {
         Temperature_FeSpace->GetDOFPosition(ik,x1,y1);
         Velocity_FeSpace->GetDOFPosition(ik,x2,y2);
         cout << ik << "  " << x1 << " " << y1 << " " << x2 << " " << y2 <<endl;
     }
  */
  return 0;
} // end main
