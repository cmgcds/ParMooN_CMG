// =======================================================================
//
// Purpose:     main program for ParMooN: solving a time-dependent NSE equation 
//              with Fictitious domain method
//
// Author:      Sashikumaar Ganesan
//
// History:     Implementation started on 06.07.16

// =======================================================================
#include <Domain.h>
#include <Database.h>
#include <FEDatabase2D.h>
#include <LinAlg.h>
#include <FESpace2D.h>
#include <SystemTNSE2D_FDM.h>
#include <SquareStructure2D.h>
#include <Structure2D.h>
#include <Output2D.h>
#include <LinAlg.h>
#include <CD2DErrorEstimator.h>
#include <MainUtilities.h>
#include <TimeDiscRout.h>
#include <TNSE2D_ParamRout.h>

#include <BdLine.h>
#include <BdCircle.h>
#include <MacroCell.h>
#include <IsoInterfaceJoint.h>
#include<IsoBoundEdge.h>

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
#include "../Examples/TNSE_2D/DiskInUnitSquare.h"  
// #include "Main_Users/Sashi/TNSE_2D/beam.h"
// =======================================================================
// main program
// =======================================================================
int main(int argc, char* argv[])
{
  // ======================================================================
  //  declaration of variables
  // ======================================================================
  int i, j, l, m, N_Cells, ORDER, N_U, N_P, N_L, N_TotalDOF, img=1, pressure_space_code;
  int Max_It, NSEType, velocity_space_code, N_SubSteps, Disctype, N_G, N_MovVert;
  int N_LagDOF;
    
  double *sol, *Lag_sol, *rhs, *defect, t1, t2, residual, impuls_residual, *Iso_refX;
  double limit, AllErrors[7], end_time, oldtau, tau, *LevelSetSol;
  
  TDomain *Domain;
  TDatabase *Database = new TDatabase();
  TFEDatabase2D *FEDatabase = new TFEDatabase2D(); 
  TCollection *coll, *Solid_Coll, *mortarcoll = NULL;
  TFESpace2D *Velocity_FeSpace, *Pressure_FeSpace, *fesp[2], *LevelSetFESpace, *LagrangeFESpace;
  TFESpace2D *SolidVelocity_FeSpace, *SolidPressure_FeSpace;
  TFEVectFunct2D *Velocity, *Lagrange;
  TFEFunction2D *u1, *u2, *Pressure, *fefct[4], *LevelSetFunction, *lambda1, *lambda2;
  TOutput2D *Output;
  TSystemTNSE2D_FDM *SystemMatrix_FDM;
  TAuxParam2D *aux, *NSEaux_error;
//   MultiIndex2D AllDerivatives[3] = { D00, D10, D01 };

 
  
  TFESpace2D *Projection_space;
  
//   bool ReAssembleM=TRUE;
    
  const char vtkdir[] = "VTK"; 
  char *PsBaseName, *VtkBaseName, *GEO;
  char UString[] = "u";
  char PString[] = "p";
  char WString[] = "w";
  char LString[] = "L";
  char NameString[] = "Galerkin";
  
  std::ostringstream os;
  os << " ";   
  
  mkdir(vtkdir, 0777);
    
  // ======================================================================
  // set the database values and generate mesh
  // ======================================================================    
  /** set variables' value in TDatabase using argv[1] (*.dat file), and generate the MESH based */
  Domain = new TDomain(argv[1]);  

  OpenFiles();
  OutFile.setf(std::ios::scientific);

  Database->CheckParameterConsistencyNSE();
  Database->WriteParamDB(argv[0]);
  Database->WriteTimeDB();
  ExampleFile();
  Domain->ReadGeo(TDatabase::ParamDB->GEOFILE); 
  
  /* include the mesh from a mesh generator, for a standard mesh use the 
   * build-in function. The GEOFILE describes the boundary of the domain. */ 
  if(TDatabase::ParamDB->MESH_TYPE==0)
    {
     Domain->ReadGeo(TDatabase::ParamDB->GEOFILE); 
     Domain->ConvertQuadToTri(TDatabase::ParamDB->CONVERT_QUAD_TO_TRI);  
    } // ParMooN  build-in Geo mesh
  else if(TDatabase::ParamDB->MESH_TYPE==1)  
     {Domain->GmshGen(TDatabase::ParamDB->GEOFILE); }//gmsh mesh
  else if(TDatabase::ParamDB->MESH_TYPE==2)    //triangle mesh
     {
// #if defined(__ALEDROP__) || defined(__HEMKER__) ||  defined(__BEAM__)         
       TriaReMeshGen(Domain); 
// #endif         
    } 
    else
     {  
      OutPut("Mesh Type not known, set MESH_TYPE correctly!!!" << endl);
      exit(0);
     }
        
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
// level set for phase marking
//=========================================================================   
  coll=Domain->GetCollection(It_Finest, 0);
  N_Cells = coll->GetN_Cells();
  OutPut("N_Cells : " << N_Cells <<endl);  
  
  // level set space 
  LevelSetFESpace = new TFESpace2D(coll, NameString, WString, LevelSetBoundCondition, 1, NULL);  
   
  N_G = LevelSetFESpace->GetN_DegreesOfFreedom();
  OutPut("N_G     : "<< setw(2) << N_G  << endl); 
  
  LevelSetSol = new double[2*N_G];
  memset(LevelSetSol, 0, 2*N_G*SizeOfDouble);
  LevelSetFunction = new TFEFunction2D(LevelSetFESpace, WString, WString, LevelSetSol, N_G);    
  LevelSetFunction->Interpolate(InitialLevelSet);
  
  
//=========================================================================
// construct all finite element spaces
//=========================================================================
  ORDER = TDatabase::ParamDB->ANSATZ_ORDER;
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
  N_TotalDOF = 2*N_U + N_P;

  OutPut("Dof Velocity : "<< setw(10) << 2* N_U << endl);
  OutPut("Dof Pressure : "<< setw(10) << N_P << endl);
  OutPut("Total Dof all: "<< setw(10) << N_TotalDOF  << endl);

//======================================================================
// construct all finite element functions
//======================================================================
    sol = new double[N_TotalDOF]; 
    rhs = new double[N_TotalDOF];
    
    memset(sol, 0, N_TotalDOF*SizeOfDouble);
    memset(rhs, 0, N_TotalDOF*SizeOfDouble);

    Velocity = new TFEVectFunct2D(Velocity_FeSpace, UString,  UString,  sol, N_U, 2);
    u1 = Velocity->GetComponent(0);
    u2 = Velocity->GetComponent(1);
    Pressure = new TFEFunction2D(Pressure_FeSpace, PString,  PString,  sol+2*N_U, N_P);
    
    //interpolate the initial solution
    u1->Interpolate(InitialU1);
    u2->Interpolate(InitialU2);
    Pressure->Interpolate(InitialP);    
 
   
//======================================================================
// produce outout
//======================================================================
   VtkBaseName = TDatabase::ParamDB->BASENAME;    
   Output = new TOutput2D(1, 2, 1, 1, Domain);

   Output->AddFEVectFunct(Velocity);
   Output->AddFEFunction(Pressure);
   Output->AddFEFunction(LevelSetFunction);
    
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
//   if(TDatabase::ParamDB->WRITE_VTK)
//   {
//     std::string filename(TDatabase::ParamDB->OUTPUTDIR);
//     filename += "/" + std::string(TDatabase::ParamDB->BASENAME) + ".vtk";
//     Output->WriteVtk(filename.c_str());
//   }   


//======================================================================
// SystemMatrix_ALE construction and solution
//======================================================================  
    // Disc type: GALERKIN 
    // Solver: AMG_SOLVE (or) GMG  (or) DIRECT
    SystemMatrix_FDM = new TSystemTNSE2D_FDM(Velocity_FeSpace, Pressure_FeSpace, Velocity, Pressure, sol, rhs, Disctype, NSEType,
                                             DIRECT, LevelSetFESpace, LevelSetFunction);       
    
    Solid_Coll = SystemMatrix_FDM->GetSolidColl();    
    
    
    // fespaces for velocity and pressure  
    GetVelocityAndPressureSpace(Solid_Coll, BoundCondition, mortarcoll, SolidVelocity_FeSpace,
                              SolidPressure_FeSpace, &pressure_space_code,
                              TDatabase::ParamDB->VELOCITY_SPACE,
                              TDatabase::ParamDB->PRESSURE_SPACE);
    
    // Lagrange space 
    LagrangeFESpace = new TFESpace2D(Solid_Coll, NameString, LString, LagrangeBoundCondition, 2, NULL);
    N_LagDOF = LagrangeFESpace->GetN_DegreesOfFreedom();
    
    OutPut("Dof Lagrange : "<< setw(10) << 2* N_LagDOF << endl);
    Lag_sol = new double[2*N_LagDOF]; 
    Lagrange = new TFEVectFunct2D(LagrangeFESpace, LString,  LString,  Lag_sol, N_LagDOF, 2);
    lambda1 = Lagrange->GetComponent(0);
    lambda2 = Lagrange->GetComponent(1);   

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
    
  
 
    // initilize the system matrix with the functions defined in Example file
    SystemMatrix_FDM->Init(LinCoeffs, BoundCondition, U1BoundValue, U2BoundValue, aux, NSEaux_error, SolidVelocity_FeSpace, LagrangeFESpace);
    
    // assemble M, A matrices and rhs 
    SystemMatrix_FDM->Assemble(sol, rhs);
    
    
    cout <<"main prg level set projection" << endl;
   
    exit(0);    
    
   
   /*   
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
   
   // time loop starts
   while(TDatabase::TimeDB->CURRENTTIME< end_time)
    {                                               // time cycle
     m++;
     TDatabase::TimeDB->INTERNAL_STARTTIME = TDatabase::TimeDB->CURRENTTIME;

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
            
      
      //w^{n+1}, must be computed before moving the mesh in this time step
#if defined(__HEMKER__) || defined(__BEAM__)    
      
     if(TDatabase::TimeDB->CURRENTTIME>2.0)
      {
       SystemMatrix_ALE->AssembleMeshMat();
       SystemMatrix_ALE->GetMeshVeloAndMove(N_MovVert, MovBoundVert, Free_Joint, Iso_refX, TDatabase::TimeDB->CURRENTTIME - 2.0, tau);
      }
#else       
     SystemMatrix_ALE->GetMeshVeloAndMove(TDatabase::TimeDB->CURRENTTIME, tau);
#endif           
     
      // memset(LevelSetSol, 0, 2*N_G*SizeOfDouble);
        
      // assemble only rhs, nonlinear matrix for NSE will be assemble in fixed point iteration
      // not needed if rhs is not time-dependent
      SystemMatrix_ALE->Assemble(sol, rhs); 

      //scale B matices and assemble NSE-rhs based on the \theta time stepping scheme 
      SystemMatrix_ALE->AssembleSystMat(tau, rhs, rhs, sol); 
      oldtau = tau;  

      // calculate the residual
      defect = new double[N_TotalDOF];
      memset(defect,0,N_TotalDOF*SizeOfDouble);

      SystemMatrix_ALE->GetTNSEResidual(sol, defect);

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
       SystemMatrix_ALE->Solve(sol);
             
       
       if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
        IntoL20FEFunction(sol+2*N_U, N_P, Pressure_FeSpace, velocity_space_code, pressure_space_code);   

       
       //no nonlinear iteration for Stokes problem  
       if(TDatabase::ParamDB->FLOW_PROBLEM_TYPE==STOKES)
        break;
     
        // restore the mass matrix for the next nonlinear iteration      
        SystemMatrix_ALE->RestoreMassMat();     
    
        // assemble the system matrix with given aux, sol and rhs 
        SystemMatrix_ALE->AssembleANonLinear(sol, rhs);    
      
        // assemble system mat, S = M + dt\theta_1*A
        SystemMatrix_ALE->AssembleSystMatNonLinear();        

        // get the residual
        memset(defect, 0, N_TotalDOF*SizeOfDouble);
        SystemMatrix_ALE->GetTNSEResidual(sol, defect);       
             
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
             
       } // for(j=1;j<=Max_It;j++)
// cout << " test VHM main " << endl;
// exit(0);        
 
      } // for(l=0;l<N_SubSteps;
//====================================================================== 
// measure errors to known solution
//======================================================================    
    if(TDatabase::ParamDB->MEASURE_ERRORS)
     {   
      SystemMatrix_ALE->MeasureTNSEErrors(ExactU1, ExactU2, ExactP, AllErrors);

       OutPut("L2(u): " <<   AllErrors[0] << endl);
       OutPut("H1-semi(u): " <<  AllErrors[1] << endl);
       OutPut("L2(p): " <<  AllErrors[2] << endl);
       OutPut("H1-semi(p): " <<  AllErrors[3]   << endl); 
       OutPut(AllErrors[4] <<  " l_infty(L2(u)) " <<AllErrors[5] << endl);
       OutPut(TDatabase::TimeDB->CURRENTTIME << " L2(0,t,L2)(u) : " <<   sqrt(AllErrors[6]) << endl); 
      
     } // if(TDatabase::ParamDB->MEASURE_ERRORS)

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
     }     */

  CloseFiles();
  
  return 0;
} // end main
