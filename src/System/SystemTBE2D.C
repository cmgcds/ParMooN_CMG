/** ************************************************************************ 
* @brief     source file for TSystemTBE2D
* @author    Sashikumaar Ganesan, 
* @date      28.11.20
* @History    
 ************************************************************************  */
// #ifdef __2D__

#include <SystemTBE2D.h>
#include <Database.h>
#include <SquareStructure2D.h>
#include <DiscreteForm2D.h>
#include <Assemble2D.h>
#include <FEVectFunct2D.h>
#include <AuxParam2D.h>
#include <LocalProjection.h>
#include <DirectSolver.h>
#include <MainUtilities.h>
#include <Upwind.h>
#include <TNSE2D_ParamRout.h>

#include <stdlib.h>
#include <string.h>
#include <sstream>
#include <MooNMD_Io.h>

#ifdef _SMPI
#include <SeqParDirectSolver.h>
#endif

TSystemTBE2D::TSystemTBE2D(TFESpace2D *velocity_fespace,TFEVectFunct2D *Velocity, 
                           double *sol, double *rhs, int disctype, int solver)
{

  //store the FEspaces and fefunct
  FeSpace = velocity_fespace;
  fesp[0] = FeSpace;
  VelocityFct = Velocity;
  N_DOF = FeSpace->GetN_DegreesOfFreedom();
  N_Active =  FeSpace->GetActiveBound();
  N_DirichletDof = N_DOF - N_Active;  
  SOLVER = solver;

  RHSs[0] = rhs;
  RHSs[1] = rhs + N_DOF;

  Disctype = disctype;


  B = new double[2*N_DOF];
  defect = new double[2*N_DOF];
  FeFct[0] = Velocity->GetComponent(0);
  FeFct[1] = Velocity->GetComponent(1); 

  gamma =0.;  
  
  // build matrices
  // first build matrix structure
  sqstructure = new TSquareStructure2D(FeSpace);
  sqstructure->Sort();  // sort column numbers: numbers are in increasing order


  // allocate the matrices
  SqmatrixA = new TSquareMatrix2D(sqstructure);    
  SqmatrixM = new TSquareMatrix2D(sqstructure);

  // matrices for methods
  sqmatrices = (TSquareMatrix **)SQMATRICES;
  sqmatrices[0] = SqmatrixA;
  sqmatrices[1] = SqmatrixM;

  BEaux = NULL;
  SystMatAssembled  = FALSE;
}

TSystemTBE2D::~TSystemTBE2D()
{
   delete [] defect; delete [] B;
   delete SqmatrixM; delete SqmatrixA;
   delete BEaux;
}

void TSystemTBE2D::Init(CoeffFct2D *lincoeffs, BoundCondFunct2D *BoundCond, BoundValueFunct2D *U1BoundValue,  
                        BoundValueFunct2D *U2BoundValue)
{
  TDiscreteForm2D *DiscreteFormGalerkin;
  TDiscreteForm2D *DiscreteFormSDFEM;
  TDiscreteForm2D *DiscreteFormUpwind;
  TDiscreteForm2D *DiscreteFormNLGalerkin;
  TDiscreteForm2D *DiscreteFormNLSDFEM;
  TDiscreteForm2D *DiscreteFormNLUpwind;
  TDiscreteForm2D *DiscreteFormRHS;
 
  
  // save the boundary condition
  BoundaryConditions[0] = BoundCond;
  BoundaryConditions[1] = BoundCond;  

  // save the boundary values  
  BoundaryValues[0] = U1BoundValue;
  BoundaryValues[1] = U2BoundValue;
 
  // save the nse bilinear coefficient   
  LinCoeffs[0] = lincoeffs;
  
  // 2 parameters are needed for assembling (u1_old, u2_old)
  BEaux = new TAuxParam2D(TimeNSN_FESpaces2, TimeNSN_Fct2, TimeNSN_ParamFct2,
                          TimeNSN_FEValues2, fesp, FeFct, TimeNSFct2, TimeNSFEFctIndex2, 
                          TimeNSFEMultiIndex2, TimeNSN_Params2, TimeNSBeginParam2);  
  
  // aux for calculating the error
   //   BEaux_error = nseaux_error;
  
  // set the Discreteforms
  InitializeDiscreteFormsBurgers(DiscreteFormGalerkin, DiscreteFormSDFEM, DiscreteFormUpwind, DiscreteFormNLGalerkin,
                                 DiscreteFormNLSDFEM, DiscreteFormNLUpwind, DiscreteFormRHS, lincoeffs);
  
    // find discrete form
    switch(Disctype)
       {
          case GALERKIN:
            DiscreteFormMARhs = DiscreteFormGalerkin;
            DiscreteFormNL = DiscreteFormNLGalerkin;
          break;

          case UPWIND:
            DiscreteFormMARhs = DiscreteFormUpwind;
            DiscreteFormNL = DiscreteFormNLUpwind;    
            break;

          case SDFEM:
            DiscreteFormMARhs = DiscreteFormSDFEM;
            DiscreteFormNL = DiscreteFormNLSDFEM;                 
            break;
	    
          default:
            Error("Unknown DISCTYPE" << Disctype << endl);
            exit(-1);
        } 
     
   TFESpace2D *fesprhs[2];
   fesprhs[0] = FeSpace;
   fesprhs[1] = FeSpace;

   //assemble object
   AMatRhsAssemble = new TAssembleMat2D(1, fesp, 2, SQMATRICES, 0, NULL,
                              2, RHSs, fesprhs, DiscreteFormMARhs, BoundaryConditions, BoundaryValues, BEaux);   
   AMatRhsAssemble->Init();     

   AMatAssembleNonLinear = new TAssembleMat2D(1, fesp, 1, SQMATRICES, 0, NULL,
                              0, NULL, NULL, DiscreteFormNL, BoundaryConditions, BoundaryValues, BEaux);   
   AMatAssembleNonLinear->Init();          

} // TSystemTBE2D::Init


/* Assemble M, A and rhs */ 
void TSystemTBE2D::Assemble(double *sol, double *rhs)
{

  /** initialize matrices */
  AMatRhsAssemble->Reset();
  
  /** assemble */
  AMatRhsAssemble->Assemble2D();  
  

//     /** upwind */ 
//   if( (Disctype==UPWIND) && !(TDatabase::ParamDB->PROBLEM_TYPE == 3) )
//    {
//     this->UpdateUpwind(); 
//    }
     
  // set rhs for Dirichlet nodes
  memcpy(sol+N_Active, rhs+N_Active, N_DirichletDof*SizeOfDouble);
  memcpy(sol+(N_DOF+N_Active), rhs+(N_DOF+N_Active), N_DirichletDof*SizeOfDouble); 
} // TSystemTBE2D::Assemble(T

  
void TSystemTBE2D::AssembleA()
{

  /** initialize matrices */
  AMatAssembleNonLinear->Reset();
  
  /** assemble */
  AMatAssembleNonLinear->Assemble2D();  
  

//     /** upwind */ 
//   if( (Disctype==UPWIND) && !(TDatabase::ParamDB->PROBLEM_TYPE == 3) )
//    {
//     this->UpdateUpwind(); 
//    }
     
} // TSystemTBE2D::Assemble(T  

void TSystemTBE2D::AssembleSystMat(double *oldrhs, double *rhs, double *sol)
{
 double tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
     
  memset(B, 0, (2*N_DOF)*SizeOfDouble);    

  // old rhs multiplied with current subtime step and theta3 on B
  Daxpy(N_Active, tau*TDatabase::TimeDB->THETA3, oldrhs, B);
  Daxpy(N_Active, tau*TDatabase::TimeDB->THETA3, oldrhs+N_DOF, B+N_DOF);   

  // add rhs from current sub time step to rhs array B
  Daxpy(N_Active, tau*TDatabase::TimeDB->THETA4, rhs, B);
  Daxpy(N_Active, tau*TDatabase::TimeDB->THETA4, rhs+N_DOF, B+N_DOF);   
      
  // Also currently : M := M + gamma A
  // M = M + (-gamma - tau*TDatabase::TimeDB->THETA2) A 
  // defect = M * sol
  // B:= B + defect 
  MatAdd(SqmatrixM, SqmatrixA, -tau*TDatabase::TimeDB->THETA2);     
  gamma = - tau*TDatabase::TimeDB->THETA2;

  memset(defect, 0, (2*N_DOF)*SizeOfDouble);  
  MatVectActive(SqmatrixM, sol, defect);
  MatVectActive(SqmatrixM, sol+N_DOF, defect+N_DOF);
  Daxpy(N_Active, 1, defect, B);
  Daxpy(N_Active, 1, defect+N_DOF, B+N_DOF);
 
  // assembling of system matrix       
  MatAdd(SqmatrixM, SqmatrixA, -gamma + tau*TDatabase::TimeDB->THETA1);   
  gamma = tau*TDatabase::TimeDB->THETA1;
   
  // set rhs for Dirichlet nodes
  memcpy(B+N_Active, rhs+N_Active, N_DirichletDof*SizeOfDouble);
  memcpy(B+(N_DOF+N_Active), rhs+(N_DOF+N_Active), N_DirichletDof*SizeOfDouble); 
     
  SystMatAssembled  = TRUE;
} // AssembleSystMat


// /* assemble only LHS, not rhs */
// void TSystemTNSE2D::AssembleSystMatNonLinear()
// {
//  double tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
     
//   if(SystMatAssembled)
//    {
//     cout << "Restore System mat before calling AssembleSystMat" <<endl;
//    }

//    switch(NSEType)
//     {
//      case 1:
//      case 2: 
//        // assembling of system matrix       
//        MatAdd(SqmatrixM11, SqmatrixA11,   tau*TDatabase::TimeDB->THETA1);   
//        gamma = tau*TDatabase::TimeDB->THETA1;
//      break;

//      case 3:
//      case 4:       
//        //assembling system matrix
// //       cout <<   
         
//        MatAdd(SqmatrixM11, SqmatrixA11, tau*TDatabase::TimeDB->THETA1);
//        MatAdd(SqmatrixM12, SqmatrixA12,  tau*TDatabase::TimeDB->THETA1);
//        MatAdd(SqmatrixM21, SqmatrixA21,  tau*TDatabase::TimeDB->THETA1);
//        MatAdd(SqmatrixM22, SqmatrixA22,  tau*TDatabase::TimeDB->THETA1);       
//        gamma = tau*TDatabase::TimeDB->THETA1;     
//      break;
//     } 

//    SystMatAssembled  = TRUE;
// } // AssembleSystMatNonLinear




// void TSystemTNSE2D::RestoreMassMat()
// {

// //   cout << "RestoreMassMat  gamma " << gamma << endl;
//   if(SystMatAssembled)
//    {
//     // restore the mass matrix
//     switch(NSEType)
//      {
//       case 1:
//       case 2:
//        MatAdd(SqmatrixM11, SqmatrixA11, -gamma);          
//        gamma = 0.;
//       break;

//      case 3:
//      case 4:
//        //assembling system matrix
//        MatAdd(SqmatrixM11, SqmatrixA11, -gamma);
//        MatAdd(SqmatrixM12, SqmatrixA12, -gamma);
//        MatAdd(SqmatrixM21, SqmatrixA21, -gamma);
//        MatAdd(SqmatrixM22, SqmatrixA22, -gamma);       
//        gamma = 0.;     
//      break;
//     } 
    
//     SystMatAssembled  = FALSE;  
//    }
//   else
//   {
//     cout << "System is not assembled to restore " <<endl;
//   }
   
// //   cout << "RestoreMassMat" << endl;
// //   exit(0);
  
// }


// void TSystemTNSE2D::AssembleANonLinear(double *sol, double *rhs)
// {
//  int N_SquareMatrices, N_RectMatrices, N_Rhs, N_FESpaces, last_sq;

//      N_RectMatrices = 0;          
//      N_Rhs = 0;
//      N_FESpaces = 1;
 
//      // set the nonliner matrices
//       switch(TDatabase::ParamDB->NSTYPE)
//        {
//         case 1:
//         case 2:
//           SQMATRICES[0] = SqmatrixA11;
//           SQMATRICES[0]->Reset();

//           N_SquareMatrices = 1;
//         break;

//         case 3:
//         case 4:
//           if (TDatabase::ParamDB->SC_NONLIN_ITE_TYPE_SADDLE==0)
//            {
//             SQMATRICES[0] = SqmatrixA11;
//             SQMATRICES[1] = SqmatrixA22;
//             SQMATRICES[0]->Reset();
//             SQMATRICES[1]->Reset();

//             N_SquareMatrices = 2;
//             last_sq = 1;

// #ifdef __PRIVATE__   
//             if (Disctype == VMS_PROJECTION)
//               {
//                SQMATRICES[0] = SqmatrixA11;
//                SQMATRICES[1] = SqmatrixA12;
//                SQMATRICES[2] = SqmatrixA21;
//                SQMATRICES[3] = SqmatrixA22;
//                SQMATRICES[0]->Reset();
//                SQMATRICES[1]->Reset();
//                SQMATRICES[2]->Reset();
//                SQMATRICES[3]->Reset();

//                N_SquareMatrices = 4;
//                last_sq = 3;

//                N_RectMatrices = 2;
//                MATRICES[0] = Matrices_tilde_G11;
//                MATRICES[1] = Matrices_tilde_G22;
//                MATRICES[0]->Reset();
//                MATRICES[1]->Reset();
       
//                N_FESpaces = 4;
//               }  
// #endif

//            }
//           else
//            {
//             // Newton method
//             cout<< "Newton method not tested " <<endl;
//             exit(0);
//            }

//          break;
//         } // switch(TDatabase::ParamDB->NSTYPE)
         
    
//       // assemble the nonlinear part of NSE
//       Assemble2D(N_FESpaces, FeSpaces,
//                  N_SquareMatrices, SQMATRICES,
//                  N_RectMatrices, MATRICES,
//                  N_Rhs, NULL, NULL,
//                  DiscreteFormNL,
//                  BoundaryConditions,
//                  BoundaryValues,
//                  NSEaux);    

//        // apply upwind disc
//       if( (Disctype==UPWIND) && !(TDatabase::ParamDB->PROBLEM_TYPE == STOKES) )
//        {
//         switch(NSEType)
//          {
//           case 1:
//           case 2:
//             // do upwinding with one matrix
//             UpwindForNavierStokes(LinCoeffs[0], SQMATRICES[0], FeFct[0], FeFct[1]);
//             cout << "UPWINDING DONE : level " << endl;
//             break;

//           case 3:
//           case 4:
//             // do upwinding with two matrices
//             cout << "UPWINDING DONE : level " << endl;
//             UpwindForNavierStokes(LinCoeffs[0], SQMATRICES[0], FeFct[0], FeFct[1]);
//             UpwindForNavierStokes(LinCoeffs[0],SQMATRICES[last_sq], FeFct[0], FeFct[1]);
//           break;
//          }                        // endswitch
//        }                          // endif     
       
       
//       // slip with boundary condition
//       if (TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION >= 1)
//       { 
//         N_FESpaces = 1;
//         N_SquareMatrices = 2;
//         N_RectMatrices = 0;
//         N_Rhs = 0;

//         SQMATRICES[0] = SqmatrixA11;
//         SQMATRICES[1] = SqmatrixA22;

//         Assemble2DSlipBC(N_FESpaces, FeSpaces,
//                          N_SquareMatrices, SQMATRICES,
//                          N_RectMatrices, NULL,
//                          N_Rhs, NULL, NULL,
//                          NULL,
//                          BoundaryConditions,
//                          BoundaryValues,
//                          NSEaux,
//                          FeFct[0], FeFct[1]);

//       }// (TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION >=         

//      // set rhs for Dirichlet nodes
//      memcpy(sol+N_Active, rhs+N_Active, N_DirichletDof*SizeOfDouble);
//      memcpy(sol+N_U+N_Active, rhs+N_U+N_Active, N_DirichletDof*SizeOfDouble); 
     
     

     
     
     
     
// #ifdef __PRIVATE__   
//       // update matrices
//       if (Disctype == VMS_PROJECTION)
//         {
//          SQMATRICES[0] = SqmatrixA11;
//          SQMATRICES[1] = SqmatrixA12;
//          SQMATRICES[2] = SqmatrixA21;
//          SQMATRICES[3] = SqmatrixA22;
//          SQMATRICES[6] =  MatricesL;
//          MATRICES[2] = Matrices_tilde_G11;
//          MATRICES[3] = Matrices_tilde_G22;
//          MATRICES[4] = Matrices_G11;
//          MATRICES[5] = Matrices_G22;

//          VMSProjectionUpdateMatrices(N_U, FeSpaces[0]->GetActiveBound(), FeSpaces[3]->GetN_DegreesOfFreedom(),
//                                      SQMATRICES, MATRICES);
//         }
// #endif       
     
// } //TSystemTNSE2D::AssembleNonLinear(

 
void TSystemTBE2D::Solve(double *sol)
{  
  if(!SystMatAssembled)
  {
    cout << "System Matrix is not assembled to solve " <<endl;
    exit(0);
  }
  
    switch(SOLVER)
     {
      case AMG_SOLVE:
        cout << "AMG_SOLVE not yet implemented " <<endl;
      break;

      case GMG:
        cout << "GMG solver not yet implemented " <<endl;
      break;

      case DIRECT:
// #ifdef _SMPI	
//               P_DS->Solve(sol,B,true);
// #else
              DirectSolver(SqmatrixM, B, sol);
// #endif
      break;      
 
      default:
            OutPut("Unknown Solver" << endl);
            exit(4711);;
     }    

}



void TSystemTBE2D::GetTBEResidual(double *sol, double *res)
{
  double residual_scalar;

  if(SystMatAssembled)
   {
    memset(res, 0, 2*N_DOF*SizeOfDouble);  
    ScalarDefect(SqmatrixM, sol, B, res, residual_scalar);   
    ScalarDefect(SqmatrixM, sol+N_DOF, B+N_DOF, res+N_DOF, residual_scalar);       
   }
  else
   {
    OutPut("Assemble the System Matrix before calculating the GetResidual" << endl);
    exit(4711);;   
   }
             
} // TSystemTBE2D::GetResidual


// void TSystemTNSE2D::MeasureTNSEErrors(DoubleFunct2D *ExactU1, DoubleFunct2D *ExactU2, DoubleFunct2D *ExactP,
//                                     double *AllErrors)
// {
//   MultiIndex2D TimeNSAllDerivatives[3] = { D00, D10, D01 };

//   double errors[4],  u_error[4];    
    
//      // errors in first velocity component
//      FeFct[0]->GetErrors(ExactU1, 3, TimeNSAllDerivatives, 2,
//                          L2H1Errors,
//                          NULL, NSEaux_error, 1, FeSpaces, errors);
//       u_error[0] = errors[0];
//       u_error[1] = errors[1];
      
//      // errors in second velocity component
//      FeFct[1]->GetErrors(ExactU2, 3, TimeNSAllDerivatives, 2,
//                          L2H1Errors,
//                          NULL, NSEaux_error, 1, FeSpaces, errors);
//      u_error[2] = errors[0];
//      u_error[3] = errors[1];      
      
//       // errors in pressure
//      FeFct[2]->GetErrors(ExactP, 3, TimeNSAllDerivatives, 2,
//                          L2H1Errors,
//                          NULL, NSEaux_error, 1, FeSpaces+1, errors);     

     
//      // calculate all errors
//      AllErrors[0] = sqrt(u_error[0]*u_error[0]+u_error[2]*u_error[2]);
//      AllErrors[1] = sqrt(u_error[1]*u_error[1]+u_error[3]*u_error[3]);
//      AllErrors[2] = errors[0];
//      AllErrors[3] = errors[1];    
     
//       // error in L^infty(0,t,L^2)
//       if(AllErrors[0] > AllErrors[5])
//        {
//         AllErrors[5]  = AllErrors[0];
//         AllErrors[4]  =  TDatabase::TimeDB->CURRENTTIME;
//       }

      
//       // error in L^2(0,t,L^2)    
//       AllErrors[6] += (u_error[0]*u_error[0] + u_error[2]*u_error[2] +olderror_l_2_l_2u)*TDatabase::TimeDB->TIMESTEPLENGTH/2.0;      
//       olderror_l_2_l_2u = u_error[0]*u_error[0] + u_error[2]*u_error[2];
     
// }

// // ************************************************************************************************************************ //

// void printall( TSquareMatrix2D * MAT1, char *MAT1_name
// //            ,TSquareMatrix2D * A12 
// //            ,TSquareMatrix2D * A21
//               ,TSquareMatrix2D * MAT2, char *MAT2_name
// //            ,TSquareMatrix2D * A23
// //            ,TSquareMatrix2D * A31
// //            ,TSquareMatrix2D * A32
// //            ,TMatrix* B1T
// //            ,TMatrix* B2T
// //            ,TMatrix* B1
// //            ,TMatrix* B2,
//                     )
// {
//     double * entries;
    
//        cout << "\n*******************************************************************************\n";
//        cout << MAT1_name << "\t";
//        entries  = MAT1->GetEntries();
//        for(int ii =0; ii< 20 ; ii++ )
//            cout << entries[ii] << "  " ;
       
//        cout << endl;
//        cout << MAT2_name << "\t";
//        entries  = MAT2->GetEntries();       
//        for(int ii =0; ii< 20 ; ii++ )
//            cout << entries[ii] << "  " ;
       
//        cout << endl;
//  }

// #endif // #ifdef __2D__
