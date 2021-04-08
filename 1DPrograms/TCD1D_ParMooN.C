// =======================================================================
// Purpose:     main program for 1D with new kernels of ParMooN
// Author:      Sashikumaar Ganesan
//
// History:     Implementation started on 12.12.2020
// =======================================================================
#ifdef _MPI
# include "mpi.h"
#endif

#include <Domain.h>
#include <Database.h>
#include <SystemCD1D.h>
// #include <SystemTCD2D.h>
// #include <SystemTCD2D_ALE.h>
#include <FEDatabase2D.h>
// #include <FESpace2D.h>
// #include <SquareStructure2D.h>
// #include <Structure2D.h>
// #include <QuadAffin.h>
// #include <DirectSolver.h>
// #include <Assemble2D.h>
// #include <Output2D.h>
// #include <LinAlg.h>
// #include <CD2DErrorEstimator.h>
// #include <MainUtilities.h>
// #include <TimeDiscRout.h>

// #include <string.h>
// #include <sstream>
// #include <MooNMD_Io.h>
// #include <stdlib.h>
// #include <math.h>
// #include <sys/stat.h>
// #include <sys/types.h>

// #include <SystemADI1D.h>
// #include <SystemADI1D_3L.h>
// // =======================================================================
// // include current example
// // =======================================================================

#include "../Examples/CD_1D/cd1d.h"

// // =======================================================================
// // include file for internal functions
// // =======================================================================
// // #include "PopulationBalance.h"

int main(int argc, char* argv[])
{
  #ifdef _MPI
  const int root = 0;
  int rank, mpi_size;
  double t_par1, t_par2;
  char  name[MPI_MAX_PROCESSOR_NAME];

  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);  
  #endif 

//   int i, j, l, m, N_SubSteps, ORDER, N_Cells, N_DOF, img=1, N_G;
//   int N_Active;

//   double *sol, *rhs, *oldrhs, t1, t2, errors[5], Linfty;
//   double tau, end_time, *defect, olderror, olderror1, hmin, hmax;

//   bool UpdateStiffnessMat, UpdateRhs,  ConvectionFirstTime;
//   char *VtkBaseName, *CovidPopulationName, *CovidNucleationName, *CovidRecoverdName, *CovidDeathName;
//   const char vtkdir[] = "VTK"; 
//   const char Populationdir[] = "PopulationData"; 

  // TDomain *Domain;
  TDatabase *Database = new TDatabase();
  TSystemCD1D *SystemCD;
  TFEDatabase2D *FEDatabase = new TFEDatabase2D(); 
//   TCollection *coll;
//   TFESpace2D *Scalar_FeSpace, *fesp[1];
//   TFEFunction2D *Scalar_FeFunction;
//   TOutput2D *Output;
//   TSystemTCD2D *SystemMatrix;  
//   TAuxParam2D *aux;
//   MultiIndex2D AllDerivatives[3] = {D00, D10, D01}; 
 
  std::ostringstream os;
  os << " ";   
  //======================================================================
  // CD1D System construction 
  //======================================================================  
  SystemCD = new TSystemCD1D(100, 0, 1, BoundCondition_LminLMax, ConvectFunct, argv[1]);
 
  // initiaize system
  SystemCD->Init(BilinearCoeffs);

  // interpolate (if needed)
  SystemCD->Interpolate(Exact);

//   //======================================================================
//   // parameters for time stepping scheme
//   //======================================================================  
//   m = 0;
//   N_SubSteps = GetN_SubSteps();
//   end_time = TDatabase::TimeDB->ENDTIME;

// #ifdef _MPI
//  if(TDatabase::ParamDB->Par_P0==1)
//  #endif  
//   OutPut(endl << "CURRENT TIME:  "<< TDatabase::TimeDB->CURRENTTIME << endl); 

//  #ifdef __LOCKDOWNMODEL__
//   double InitialPopulation[32];
// #ifdef _MPI
//   if(N_AllCell!=32)
//     {

//      if(TDatabase::ParamDB->Par_P0==1)
//       OutPut("Lock down model must have 32 States/UT"<< N_AllCell <<endl);
//       //  MPI_Finalize();
//       // exit(0);  
//     }
// #endif   
//   char *DataFile;
//   DataFile = new char[17];
//   strcpy(DataFile, "InitialData.dat");

//   ReadData(DataFile, InitialPopulation);

//   //Interpolate 
//   ADISystem3L->Interpolate(coll, InitialPopulation, InitialValues);
//  #else 
//   //Interpolate 
//   ADISystem3L->Interpolate(InitialValues);
//  #endif

//   // XnodalLnodal->XnodalLdof->LdofXdof (note XnodalLdof is needed for IntL)
//   ADISystem3L->XnodalLnodalToLdofXdof(MaxN_PtsForNodal, Sol_LdofXdof);

//   //  MPI_Finalize();
//   //  exit(0); 
//   //======================================================================
//   // produce outout at t=0
//   //======================================================================
//   VtkBaseName = TDatabase::ParamDB->VTKBASENAME;   
//   CovidPopulationName = TDatabase::ParamDB->MAPFILE;   
//   CovidNucleationName = TDatabase::ParamDB->PODFILE;   
//   CovidRecoverdName = TDatabase::ParamDB->POD_FILENAME; 
//   CovidDeathName = TDatabase::ParamDB->SNAP_FILENAME; 

//   Output = new TOutput2D(2, 2, 1, 1, Domain);

//   Output->AddFEFunction(Scalar_FeFunction);

//   int i_la, i_ld, i_lv;
//   if(TDatabase::ParamDB->WRITE_VTK)
//    {
//     for(j=0;j<N_LNodals_All;j++)
//      {
//       // i_la =  (j/N_LnodalPos[0])/N_LnodalPos[1];
//       // i_ld = (j/N_LnodalPos[0])%N_LnodalPos[1];
//       // i_lv = j%N_LnodalPos[0];
//       // cout << i_la << " : " << i_ld << " : " << i_lv <<endl;

//       memcpy(sol, Sol_LdofXdof +j*N_DOF, N_DOF*SizeOfDouble);
//       os.seekp(std::ios::beg);
//        if(img<10) os <<  "VTK/"<<VtkBaseName<<j<<".0000"<<img<<".vtk" << ends;
//          else if(img<100) os <<  "VTK/"<<VtkBaseName<<j<<".000"<<img<<".vtk" << ends;
//           else if(img<1000) os <<  "VTK/"<<VtkBaseName<<j<<".00"<<img<<".vtk" << ends;
//            else if(img<10000) os <<  "VTK/"<<VtkBaseName<<j<<".0"<<img<<".vtk" << ends;
//             else  os <<  "VTK/"<<VtkBaseName<<j<<"."<<img<<".vtk" << ends;
//       Output->WriteVtk(os.str().c_str());
//      }   
//     img++;
//    }

//   //======================================================================
//   // Set variables at t=0
//   //======================================================================
//   int jj;
//   double *Sol_Xdof, GrowthFact_X=1., RecoveryRatio=0.;
//   // double C_old = TDatabase::ParamDB->REACTOR_P2;
//   // double C_new = TDatabase::ParamDB->REACTOR_P2;

//   //assemble M matirces of system L
//   for(ii=0; ii<N_ADISystems; ++ii)
//    {
//     SystemADI[ii]->AssembleMassMat();
//    }

//   // ld in conservative form, to include B_nuc
//   // if(rank==0)
//   SystemADI[0]->AssembleAdvectMat(TRUE);
//   SystemADI[1]->AssembleAdvectMat(TRUE);  
//   SystemADI[2]->AssembleAdvectMat(TRUE);


//   double *IntValue = new double[N_Xpos];
//   double *IntValue_old = new double[N_Xpos];
//   double *RecoveredValue = new double[N_Xpos];
//   double *RecoveredValue_temp = new double[N_Xpos];     
//   double *RemovedValue = new double[N_Xpos];  
//   double *NucValue = new double[N_Xpos];
//   double *B_NucValue = new double[N_Xpos];
//   double *B_NucValue_old = new double[N_Xpos];
//   double *B_NucValue_tmp = new double[N_Xpos];
//   double *B_NucValue_cumulative = new double[N_Xpos];
//   double *SuscepPopRatio = new double[N_Xpos];
//   double SeroSurveyFactor = TDatabase::ParamDB->FS_L;
//   memset(RemovedValue, 0, N_Xpos*SizeOfDouble);  
//   memset(RecoveredValue, 0, N_Xpos*SizeOfDouble);  

//   //Get Int_\Omega_L values using XnodalLdof (assume xnodal & xdof are same)
//   ADISystem3L->IntL(Gamma_Q, IntValue, NucValue, LDistXSum);
//   memcpy(IntValue_old, IntValue, N_Xpos*SizeOfDouble);

//   strcpy(DataFile, "InitialData_R.dat");
//   ReadData(DataFile, InitialPopulation);

//   for(i=0; i<N_Xpos; i++)
//     B_NucValue_cumulative[i] = IntValue[i] + SeroSurveyFactor*InitialPopulation[DispArray[rank] +i];
  
//   for(i=0; i<N_Xpos; i++)
//   {
//     SuscepPopRatio[i] = (StatePopulation[i] - B_NucValue_cumulative[i])/StatePopulation[i];
//     // OutPut(i<< " State : " << DispArray[rank] +i << ": SuscepPopRatio :: "  << SuscepPopRatio[i]<<endl);
//    }

//   //  MPI_Finalize();
//   //  exit(0);
   

//   //compute Int_B_nuc
//   ADISystem3L->Int_B_Nuc(NucValue, B_NucValue_tmp, GrowthAndB_Nuc[0], SuscepPopRatio);
//   memcpy(B_NucValue_old, B_NucValue_tmp, N_Xpos*SizeOfDouble);
//   memset(B_NucValue, 0, N_Xpos*SizeOfDouble);  

//   Daxpy(N_Xpos, SeroSurveyFactor, B_NucValue, B_NucValue_cumulative);
//    for(i=0; i<N_Xpos; i++)
//     SuscepPopRatio[i] = (StatePopulation[i] - B_NucValue_cumulative[i])/StatePopulation[i];

// #ifdef _MPI
//   MPI_Gatherv(IntValue, N_Xpos, MPI_DOUBLE, CovidPopulation, N_Cells_ALL, DispArray, MPI_DOUBLE, 0, MPI_COMM_WORLD); 
//   if(rank==0)
//    {
//     WriteInitData(N_AllCell, CovidPopulation, CovidPopulationName);       
//    }

//   MPI_Gatherv(B_NucValue, N_Xpos, MPI_DOUBLE, CovidPopulation, N_Cells_ALL, DispArray, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//   if(rank==0)
//    {
//     WriteInitData(N_AllCell, CovidPopulation, CovidNucleationName);   
//     // for(i=0; i<N_AllCell; i++)
//       // OutPut("Covid Nucleation in " <<  enum_to_string(i) << " : " << CovidPopulation[i] <<endl);    
//    }


//   //  for(i=0; i<N_Xpos; i++)
//   //  {
//   //   OutPut("State : " << DispArray[rank] +i << ": SuscepPopRatio :: "  << SuscepPopRatio[i]<<endl);
//   //  }
//   //  MPI_Finalize();
//   //  exit(0);

//   MPI_Gatherv(RecoveredValue, N_Xpos, MPI_DOUBLE, CovidPopulation, N_Cells_ALL, DispArray, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//   if(rank==0)
//    {
//     WriteInitData(N_AllCell, CovidPopulation, CovidRecoverdName);   
//     // for(i=0; i<N_AllCell; i++)
//       // OutPut("Covid Nucleation in " <<  enum_to_string(i) << " : " << CovidPopulation[i] <<endl);    
//    }

//   MPI_Gatherv(RemovedValue, N_Xpos, MPI_DOUBLE, CovidPopulation, N_Cells_ALL, DispArray, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//   if(rank==0)
//    {
//     WriteInitData(N_AllCell, CovidPopulation, CovidDeathName);   
//     // for(i=0; i<N_AllCell; i++)
//       // OutPut("Covid Nucleation in " <<  enum_to_string(i) << " : " << CovidPopulation[i] <<endl);    
//    }

//   MPI_Request request;
//   if(rank!=0)
//   {
//    MPI_Isend(LDistXSum, N_LDistXSum, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &request);

//   }
//   else
//   {
//    for(i=1; i<mpi_size; ++i)
//    {
//     MPI_Recv(sum_receiv, N_LDistXSum, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

//     for(j=0; j<N_LDistXSum; ++j)
//      LDistXSum[j] +=sum_receiv[j];
//    }

//     WriteInitDistData(N_LnodalPos, LnodalPos, N_LDistXSum, LDistXSum);  
//   }

//   if(rank==0)
//    {
//     WriteInitRecovDistData(N_LnodalPos, LnodalPos, N_LDistXSum, LRecoDistXSum);  
//    }

//   //  OutPut("SeroSurveyFactor " <<  SeroSurveyFactor <<endl);  
//   //  MPI_Finalize();
//   //  exit(0);
// #else
//   // for(i=0; i<N_Xpos; i++)
//   //  OutPut("Covid Population in " <<  enum_to_string(i) << " : " << IntValue[i] <<endl);
//  WriteInitData(N_Xpos, IntValue, CovidPopulationName);
//  WriteInitData(N_Xpos, B_NucValue, CovidNucleationName);
//  WriteInitData(N_Xpos, RecoveredValue, CovidRecoverdName);
//  WriteInitData(N_Xpos, RemovedValue, CovidDeathName);

//    // exit(0);   
// #endif

//   // double *L2Errors = new double[N_Xpos];
//   // double l2;
//   // ADISystem3L->GetErrorAllXpos(InitialValues, L2Errors);

//   // OutPut("time: " << TDatabase::TimeDB->CURRENTTIME);
//   // OutPut(" L2: " << sqrt(l2) << endl);

// #ifdef _MPI
//  t1 = MPI_Wtime(); 
// #endif


//   //  MPI_Finalize();
//   //  exit(0); 
//   double NoNucTime = 140.5; // first sun locdown day since Mar 23
//   double NoNucIncr=3.;
//   //======================================================================
//   // time disc loop
//   //======================================================================    
//    UpdateStiffnessMat = FALSE; //check BilinearCoeffs in example file
//    UpdateRhs = TRUE; //check BilinearCoeffs in example file
//    ConvectionFirstTime=TRUE;

//    // time loop starts
//    while(TDatabase::TimeDB->CURRENTTIME < end_time)
//    {
//     m++;
//     TDatabase::TimeDB->INTERNAL_STARTTIME = TDatabase::TimeDB->CURRENTTIME;

//     for(l=0; l<N_SubSteps; l++) // sub steps of fractional step theta
//     {
//       SetTimeDiscParameters(1);

// #ifdef _MPI
//  if(TDatabase::ParamDB->Par_P0==1)
//  #endif  
//       if(m==1)
//       {
//         OutPut("Theta1: " << TDatabase::TimeDB->THETA1<< endl);
//         OutPut("Theta2: " << TDatabase::TimeDB->THETA2<< endl);
//         OutPut("Theta3: " << TDatabase::TimeDB->THETA3<< endl);
//         OutPut("Theta4: " << TDatabase::TimeDB->THETA4<< endl);
//       }

//       tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
//       TDatabase::TimeDB->CURRENTTIME += tau;
 
//       //set parameters for nucleation
//       // after 93 days sd value increases by 0.02 for every 15 days,
//       Set_sigma_SD(TDatabase::TimeDB->CURRENTTIME, NoNucTime, 15., 0.000);     
//       // if(TDatabase::TimeDB->CURRENTTIME-0.5>NoNucTime)
//       //  { 
//       //   NoNucTime +=NoNucIncr;
//       //   if(NoNucIncr==3.) // no nuc on Sun & Wed 
//       //    {NoNucIncr=4.;}
//       //   else 
//       //   {NoNucIncr=3.;}
//       //  } 
// #ifdef _MPI 
//  if(TDatabase::ParamDB->Par_P0==1)
//  #endif        
//     OutPut(endl << "CURRENT TIME: " << TDatabase::TimeDB->CURRENTTIME << endl); 
//     // if(rank==0)
//       // OutPut("NoNucTime: " <<NoNucTime<< endl);    
//     //  cout << "IntVal L start "<<endl;
//     //  //Get Int_\Omega_L values using XnodalLdof (assume xnodal & xdof are same)
//     //  ADISystem3L->IntL(Gamma_Q, IntValue, NucValue);
//     //  cout << "IntVal L end "<<endl;
//      // XnodalLnodal->XnodalLdof->LdofXdof (note XnodalLdof is needed for IntL)

//       // time-dep reaction
//       GetReactionFactor(N_LnodalPos, LnodalPos, ReactionFactor);

//       ADISystem3L->XnodalLnodalToLdofXdof(MaxN_PtsForNodal, Sol_LdofXdof);
       
//       #ifdef __LOCKDOWNMODEL__
//        for(i=0; i<N_LNodals_All; ++i)
//         {
//          Sol_Xdof = Sol_LdofXdof + i*N_DOF;
//          Recovered_Xdof =  Recovery_LdofXdof + i*N_DOF;
//          GrowthFact_X = (1. - tau*TDatabase::TimeDB->THETA2*ReactionFactor[2*i])/(1. + tau*TDatabase::TimeDB->THETA1*ReactionFactor[2*i]);   
//          RecoveryRatio = ReactionFactor[2*i+1]/ReactionFactor[2*i]; //recovery ratio in the total reaction

//          for(jj=0; jj<N_DOF; ++jj)
//           {            
//             // I_R = RecoveryRatio * Totalremoved = RecoveryRatio*(I^n - I^{n+1})
//            Recovered_Xdof[jj] = RecoveryRatio*(1.-GrowthFact_X)*Sol_Xdof[jj]; // 
//            Sol_Xdof[jj] = GrowthFact_X*Sol_Xdof[jj];        
//           }
//         }
//       #else
//         cout <<"NON __LOCKDOWNMODEL__ Not yet implemented " <<endl;
//         exit(0);
//       #endif
 
//     //===============================================================================
//     // compute the recovery values
//     //===============================================================================
//      // transfer reovery values from Recovery_LdofXdof to Recovery_XnodalLdof
//      ADISystem3L->LdofXdofToXnodalLdof(Recovery_LdofXdof); 
//      //Get Int_\Omega_L values using XnodalLdof (assume xnodal & xdof are same)
//      ADISystem3L->IntL(NULL, RecoveredValue_temp, NucValue, LRecoDistXSum_temp); 
//      Daxpy(N_Xpos, 1., RecoveredValue_temp, RecoveredValue); 
//      Daxpy(N_LDistXSum, 1., LRecoDistXSum_temp, LRecoDistXSum);          
//     //===============================================================================
   

//      // transfer sol from Sol_LdofXdof to Sol_XnodalLdof
//      ADISystem3L->LdofXdofToXnodalLdof(Sol_LdofXdof); 
//      //===============================================================================
//      //system X-direction solution -- end
//      //system L-direction solution -- start  
//      //===============================================================================
  
//      //Get Int_\Omega_L values using XnodalLdof (assume xnodal & xdof are same)
//      ADISystem3L->IntL(Gamma_Q, IntValue, NucValue, NULL);

//      //compute death cases
//      Daxpy(N_Xpos, -1., RecoveredValue_temp, IntValue_old); //   y := alpha*x + y 
//      DsumS(N_Xpos, 1., IntValue_old, IntValue, RemovedValue);// Infect Death 
     
//      //===============================================================================
//      //Solve L0-direction solution -- start  
//      //===============================================================================
// #ifdef __LD__     
        
//       ADISystem3L->Int_B_Nuc(NucValue, B_NucValue_tmp, GrowthAndB_Nuc[0], SuscepPopRatio);
//       DsumP(N_Xpos, 0.5*tau, B_NucValue_old, B_NucValue_tmp, B_NucValue);
//       memcpy(B_NucValue_old, B_NucValue_tmp, N_Xpos*SizeOfDouble);

//     //  ADISystem3L->CopySolToInternal(0); // do nothing for ADISystem3L[0]
//      //solve inernal system for every x pos
//      ADISystem3L->Solve(0, BilinearCoeffs_L0, NucValue, SuscepPopRatio);

//     //  ADISystem3L->CopySolFromInternal(0); // do nothing for ADISystem3L[0]     
// #endif    
 
//      //===============================================================================
//      //Solve L0-direction solution -- completed  
//      //Solve L1-direction solution -- start  
//      //===============================================================================
// #ifdef __LV__  
//      // XnodalL2L1L0nodal -->  XnodalL2L0L1nodal       
//      ADISystem3L->CopySolToInternal(1);
   
//      // solve inernal system for every x pos
//      ADISystem3L->Solve(1, BilinearCoeffs_L1, IntValue, SuscepPopRatio);

//      // XnodalL2L0L1nodal --> XnodalL2L1L0nodal
//      ADISystem3L->CopySolFromInternal(1);  
//   // OutPut(" __LV__ Complete "  << endl);       
// #endif        
//      //===============================================================================
//      //Solve L1-direction solution -- completed  
//      //Solve L2-direction solution -- start  
//      //===============================================================================
// #ifdef __LA__      
//      // XnodalL2L1L0nodal  --> XnodalL0L1L2nodal
//      ADISystem3L->CopySolToInternal(2);

//      // solve inernal system for every x pos
//      ADISystem3L->Solve(2, BilinearCoeffs_L2, IntValue, SuscepPopRatio);
    
//      // XnodalL0L1L2nodal --> XnodalL2L1L0nodal
//      ADISystem3L->CopySolFromInternal(2);
//   // OutPut(" __LA__ Complete "  << endl);          
// #endif     
//      //===============================================================================
//      //Solve L2-direction solution -- completed  
//      //===============================================================================
//     } // for(l=0;l<N_SubSteps;l++) 
    
//      //Get Int_\Omega_L values using XnodalLdof (assume xnodal & xdof are same)   
//      if(m % TDatabase::TimeDB->STEPS_PER_IMAGE == 0)   
//      {
//       ADISystem3L->IntL(Gamma_Q, IntValue, NucValue, LDistXSum);
//      }
//      else
//      {
//       ADISystem3L->IntL(Gamma_Q, IntValue, NucValue, NULL);
//      }

//      memcpy(IntValue_old, IntValue, N_Xpos*SizeOfDouble);

//     //======================================================================
//     // produce outout
//     //======================================================================
//     if(m % TDatabase::TimeDB->STEPS_PER_IMAGE == 0)
//     {
//      #ifdef _MPI
//      MPI_Gatherv(IntValue, N_Xpos, MPI_DOUBLE, CovidPopulation, N_Cells_ALL, DispArray, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//      if(rank==0)
//      {
//       WriteData(N_AllCell, CovidPopulation, CovidPopulationName);       
//      }

//      MPI_Gatherv(B_NucValue, N_Xpos, MPI_DOUBLE, CovidPopulation, N_Cells_ALL, DispArray, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//      if(rank==0)
//      {
//       WriteData(N_AllCell, CovidPopulation, CovidNucleationName);       
//      }

//      MPI_Gatherv(RecoveredValue, N_Xpos, MPI_DOUBLE, CovidPopulation, N_Cells_ALL, DispArray, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//      if(rank==0)
//      {
//       WriteData(N_AllCell, CovidPopulation, CovidRecoverdName);       
//      }

//      MPI_Gatherv(RemovedValue, N_Xpos, MPI_DOUBLE, CovidPopulation, N_Cells_ALL, DispArray, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//      if(rank==0)
//      {
//       WriteData(N_AllCell, CovidPopulation, CovidDeathName);       
//      }

//      if(rank!=0)
//       {
//        MPI_Isend(LDistXSum, N_LDistXSum, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &request);
//       }
//     else
//       {
//        for(i=1; i<mpi_size; ++i)
//        {
//         MPI_Recv(sum_receiv, N_LDistXSum, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

//         for(j=0; j<N_LDistXSum; ++j)
//          LDistXSum[j] +=sum_receiv[j];
//         }

//         WriteDistData(N_LDistXSum, LDistXSum);  
//       }

//      if(rank!=0)
//       {
//        MPI_Isend(LRecoDistXSum, N_LDistXSum, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &request);
//       }
//     else
//       {
//        for(i=1; i<mpi_size; ++i)
//        {
//         MPI_Recv(LRecoDistXSum_temp, N_LDistXSum, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

//         for(j=0; j<N_LDistXSum; ++j)
//          LRecoDistXSum[j] +=LRecoDistXSum_temp[j];
//         }

//         WriteRecovDistData(N_LDistXSum, LRecoDistXSum);  
//       }

//     //  for(i=0; i<N_Xpos; i++)
//       // OutPut("State : " << DispArray[rank] +i << ": SuscepPopRatio :: "  << SuscepPopRatio[i]<<endl);

//      //  MPI_Finalize();
//      //  exit(0);
//      #else
//       WriteData(N_Xpos, IntValue, CovidPopulationName);
//       WriteData(N_Xpos, B_NucValue, CovidNucleationName);  
//       WriteData(N_Xpos, RecoveredValue, CovidRecoverdName);           
//       WriteInitData(N_Xpos, RemovedValue, CovidDeathName);     
//      #endif

//      Daxpy(N_Xpos, SeroSurveyFactor, B_NucValue, B_NucValue_cumulative); 

//      //vaccine schedule
//      if(TDatabase::TimeDB->CURRENTTIME>245) // 1 Aug - 31 DEC = 153 days , 245
//      {
//        for(i=0; i<N_Xpos; i++)
//         B_NucValue_cumulative[i] +=StatePopulation[i]*0.0027;
//      }

//      for(i=0; i<N_Xpos; i++)
//       {
//         SuscepPopRatio[i] = (StatePopulation[i] - B_NucValue_cumulative[i])/StatePopulation[i];
//         if(SuscepPopRatio[i]<0) SuscepPopRatio[i] = 0;
//       }

//      memset(B_NucValue, 0, N_Xpos*SizeOfDouble);  
//      memset(RemovedValue, 0, N_Xpos*SizeOfDouble);  
//      memset(RecoveredValue, 0, N_Xpos*SizeOfDouble);     
//      memset(LRecoDistXSum, 0, N_LDistXSum*SizeOfDouble);   
//     }

//    } // while(TDatabase::TimeDB->CURRENTTIME< end_time)

#ifdef _MPI
 t2 = MPI_Wtime();
 if(TDatabase::ParamDB->Par_P0==1)
  printf( "Elapsed time since t=0 is %f\n", t2 - t1 ); 
#endif  

CloseFiles();
  
#ifdef _MPI
MPI_Finalize();
#endif

return 0;
} // end main
