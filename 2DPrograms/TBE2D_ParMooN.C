// =======================================================================
// Purpose:     main program for solving a time-dependent Burger's equation in ParMooN
// Author:      Sashikumaar Ganesan
// History:     Implementation started on 28.11.2020
// =======================================================================
#include <Domain.h>
#include <Database.h>
#include <FEDatabase2D.h>
#include <LinAlg.h>
#include <FESpace2D.h>
#include <SystemTBE2D.h>
#include <SquareStructure2D.h>
#include <Output2D.h>
#include <MainUtilities.h>
#include <TimeDiscRout.h>

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
#include "../Examples/TNSE_2D/SinCos_Burger.h" // smooth sol in unit square

// =======================================================================
// main program
// =======================================================================
int main(int argc, char* argv[])
{
  // ======================================================================
  //  declaration of variables
  // ======================================================================
  int i, j, l, m, N_Cells, ORDER, N_U, N_M, N_L, N_TotalDOF, img=1, N_Modes, N_Realiz;
  int N_Total_MeanDOF, Max_It, NSEType, N_SubSteps, Disctype;
  
  double *sol, *sol_mode, *rhs, *oldrhs, *defect, t1, t2, residual, impuls_residual;
  double limit, AllErrors[7], end_time, oldtau, tau;
  
  TDomain *Domain;
  TDatabase *Database = new TDatabase();
  TFEDatabase2D *FEDatabase = new TFEDatabase2D(); 
  TCollection *coll, *mortarcoll = NULL;
  TFESpace2D *Velocity_FeSpace, *VelocityMode_FeSpace, *fesp[2];
  TFEVectFunct2D *Velocity_Mean, *Velocity_Mode;
  TFEFunction2D *u1_mean, *u2_mean, *fefct[2];
  TOutput2D *Output;
  TSystemTBE2D *SystemMatrix_Mean;
     
  const char vtkdir[] = "VTK"; 
  char *PsBaseName, *VtkBaseName, *GEO;
  char UString[] = "u_mean";
  char UMString[] = "u_mode";
  char NameString[] = "UQ";
  
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
  
  /* include the mesh from a meshgenerator, for a standard mesh use the build-in function */
  // standard mesh
   GEO = TDatabase::ParamDB->GEOFILE;
   Domain->Init(NULL, GEO);
   
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
// construct all finite element spaces
//=========================================================================
  ORDER = TDatabase::ParamDB->ANSATZ_ORDER;
  NSEType = TDatabase::ParamDB->NSTYPE;
  Disctype = TDatabase::ParamDB->DISCTYPE;

  N_Modes = TDatabase::ParamDB->P10;
  N_Realiz = TDatabase::ParamDB->P11;

  // VelocityMode = new TFEVectFunct2D*[N_Modes];

  coll=Domain->GetCollection(It_Finest, 0);
  N_Cells = coll->GetN_Cells();
  OutPut("N_Cells : " << N_Cells <<endl);
  
  Velocity_FeSpace = new TFESpace2D(coll, (char*)"Mean", (char*)"Mean", BoundCondition, ORDER, NULL);
     
  VelocityMode_FeSpace = new TFESpace2D(coll, (char*)"Mode", (char*)"Mode", BoundCondition, 1, NULL);

  N_U = Velocity_FeSpace->GetN_DegreesOfFreedom(); 
  N_M = VelocityMode_FeSpace->GetN_DegreesOfFreedom(); 
  N_Total_MeanDOF = 2*N_U;
  N_TotalDOF = 2*(N_U + N_Modes*N_M);

  OutPut("Dof Mean Velocity : "<< setw(10) << 2* N_U << endl);
  OutPut("Dof Mode Velocity : "<< setw(10) << 2* N_M << endl);
  OutPut("Total Dof all: "<< setw(10) << N_TotalDOF  << endl);
    
//======================================================================
// construct all finite element functions
//======================================================================
    sol = new double[N_TotalDOF]; 
    rhs = new double[N_TotalDOF];
    oldrhs = new double[N_TotalDOF];
    defect = new double[N_TotalDOF];

    memset(sol, 0, N_TotalDOF*SizeOfDouble);
    memset(rhs, 0, N_TotalDOF*SizeOfDouble);

    /** mean velo */
    Velocity_Mean = new TFEVectFunct2D(Velocity_FeSpace, UString,  UString,  sol, N_U, 2);
    u1_mean = Velocity_Mean->GetComponent(0);
    u2_mean = Velocity_Mean->GetComponent(1);
    
    /** mode velo */
    sol_mode = sol+2*N_U;    
    Velocity_Mode = new TFEVectFunct2D(VelocityMode_FeSpace, UString,  UString,  sol_mode, N_M, 2*N_Modes);

    //interpolate the initial solution
    u1_mean->Interpolate(InitialU1Mean);
    u2_mean->Interpolate(InitialU2Mean);

    // Sol_DO_Coeff = new double[N_Realiz*N_modes]
//======================================================================
// SystemMatrix construction and solution
//======================================================================  
    // Disc type: GALERKIN 
    // Solver: AMG_SOLVE (or) GMG  (or) DIRECT
    SystemMatrix_Mean = new TSystemTBE2D(Velocity_FeSpace, Velocity_Mean, sol, rhs, Disctype, DIRECT);

    // initilize the system matrix with the functions defined in Example file
    // last argument is aux that is used to pass additional fe functions (eg. mesh velocity)    
    SystemMatrix_Mean->Init(LinCoeffs, BoundCondition, U1BoundValue, U2BoundValue);
    
    // assemble M, A matrices and rhs 
    SystemMatrix_Mean->Assemble(sol, rhs);
  
//======================================================================
// produce outout
//======================================================================
   VtkBaseName = TDatabase::ParamDB->VTKBASENAME;    
   Output = new TOutput2D(2, 2, 1, 1, Domain);

   Output->AddFEVectFunct(Velocity_Mean);
   
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
       
      //copy sol, rhs to olssol, oldrhs
      memcpy(oldrhs, rhs, N_Total_MeanDOF*SizeOfDouble);        
    
      // assemble only rhs, nonlinear matrix for NSE will be assemble in fixed point iteration
      // not needed if rhs is not time-dependent
      if(m!=1)
       { SystemMatrix_Mean->AssembleA(); }
      else
       { SystemMatrix_Mean->Assemble(sol, rhs);  }
       
      //scale B matices and assemble NSE-rhs based on the \theta time stepping scheme 
      SystemMatrix_Mean->AssembleSystMat(oldrhs, rhs, sol); 
      oldtau = tau;
         
      // calculate the residual
      SystemMatrix_Mean->GetTBEResidual(sol, defect);
    
      residual =  Ddot(N_Total_MeanDOF, defect, defect);
      OutPut("Nonlinear iteration step   0");
      OutPut(setw(14) << sqrt(residual) << endl);       

//====================================================================== 
//Solve the system
//Nonlinear iteration of fixed point type
//======================================================================
     for(j=1;j<=Max_It;j++)
      {   
       // Solve the NSE system
       SystemMatrix_Mean->Solve(sol);
    
       // restore the mass matrix for the next nonlinear iteration      
       SystemMatrix_Mean->RestoreMassMat();     
    
       // assemble the system matrix with given aux, sol and rhs 
       SystemMatrix_Mean->AssembleANonLinear(sol, rhs);    
          
       // assemble system mat, S = M + dt\theta_1*A
       SystemMatrix_Mean->AssembleSystMatNonLinear();        

       // get the residual
       memset(defect, 0, N_TotalDOF*SizeOfDouble);
       SystemMatrix_Mean->GetTBEResidual(sol, defect);       
 
    
       residual =  Ddot(N_Total_MeanDOF, defect, defect);
       OutPut("nonlinear iteration step " << setw(3) << j);
       OutPut(setw(14) << sqrt(residual) << endl); 
	
       if(sqrt(residual)<=limit)
       break;

      } // for(j=1;j<=Max_It;j++)
 
      // restore the mass matrix for the next time step          
      SystemMatrix_Mean->RestoreMassMat();     
       
   } // for(l=0;l<N_SubSteps;
//====================================================================== 
// measure errors to known solution
//======================================================================    
    if(TDatabase::ParamDB->MEASURE_ERRORS)
     {          
      SystemMatrix_Mean->MeasureErrors(ExactU1, ExactU2, AllErrors);

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
     }     

  CloseFiles();
  
  return 0;
} // end main
