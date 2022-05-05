/** ************************************************************************ 
* @brief     source file for TSystemTNSE2D
* @author    Sashikumaar Ganesan, 
* @date      23.08.14
* @History   07/02/2017 SMPI Checked Working - Ankit
 ************************************************************************  */
#ifdef __2D__

#include <Database.h>
#include <SystemNSE2D.h>
#include <SystemTNSE2D.h>
#include <SquareStructure2D.h>
#include <DiscreteForm2D.h>
#include <Assemble2D.h>
#include <FEVectFunct2D.h>
#include <AuxParam2D.h>
#include <LocalProjection.h>
#include <DirectSolver.h>
// #include <TNSE2D_ParamRout.h>
#include <MainUtilities.h>
#include <Upwind.h>

#include <stdlib.h>
#include <string.h>
#include <sstream>
#include <MooNMD_Io.h>

#ifdef __PRIVATE__  
 #include <VMS.h>
#endif

#ifdef _SMPI
#include <SeqParDirectSolver.h>
#endif

void printall( TSquareMatrix2D * MAT1, char *MAT1_name
//            ,TSquareMatrix2D * A12 
//            ,TSquareMatrix2D * A21
              ,TSquareMatrix2D * MAT2, char *MAT2_name
//            ,TSquareMatrix2D * A23
//            ,TSquareMatrix2D * A31
//            ,TSquareMatrix2D * A32
//            ,TMatrix* B1T
//            ,TMatrix* B2T
//            ,TMatrix* B1
//            ,TMatrix* B2,
                    ) ;


TSystemTNSE2D::TSystemTNSE2D(TFESpace2D *velocity_fespace, TFESpace2D *presssure_fespace, TFEVectFunct2D *Velocity, 
                                   TFEFunction2D *p, double *sol, double *rhs, int disctype, int nsetype, int solver
#ifdef __PRIVATE__  
                                   ,TFESpace2D *Projection_space, TFESpace2D *Stress_FeSpace,   TFESpace2D *Deformation_FeSpace
#endif    
                                   ) : TSystemNSE2D(velocity_fespace, presssure_fespace, Velocity, p, sol, rhs, disctype, nsetype, solver
#ifdef __PRIVATE__  
                                   ,Projection_space, Stress_FeSpace, Deformation_FeSpace
#endif   
                                   )
{
  B = new double[2*N_U+N_P];
  defect = new double[2*N_U+N_P];

  gamma =0.;  
  
    // allocate the mass matrices in addition
    switch(NSEType)
     {
      case 1:
      case 2:
        SqmatrixM11 = new TSquareMatrix2D(sqstructureA);

        sqmatrices[0] = SqmatrixM11;
      break;

      case 3:
      case 4:
        SqmatrixM11 = new TSquareMatrix2D(sqstructureA);
        SqmatrixM12 = new TSquareMatrix2D(sqstructureA);
        SqmatrixM21 = new TSquareMatrix2D(sqstructureA);
        SqmatrixM22 = new TSquareMatrix2D(sqstructureA);

        sqmatrices[0] = SqmatrixM11;
        sqmatrices[1] = SqmatrixM12;
        sqmatrices[2] = SqmatrixM21;
        sqmatrices[3] = SqmatrixM22;
      break;
      
      default:
            OutPut("Unknown NSETYPE, it must be 1 to 4" << endl);
            exit(4711);;      
      
     }  
 

 NSE_Rhsaux = NULL;
 SystMatAssembled  = FALSE;
 olderror_l_2_l_2u = 0.;
}

TSystemTNSE2D::~TSystemTNSE2D()
{
   delete [] defect; delete [] B;
   delete SqmatrixM11; delete SqmatrixM12; delete SqmatrixM21; delete SqmatrixM22;
   delete NSE_Rhsaux;
}

void TSystemTNSE2D::Init(CoeffFct2D *lincoeffs, BoundCondFunct2D *BoundCond, BoundValueFunct2D *U1BoundValue, BoundValueFunct2D *U2BoundValue,
                            TAuxParam2D *aux, TAuxParam2D *nseaux_error)
{
  TDiscreteForm2D *DiscreteFormGalerkin;
  TDiscreteForm2D *DiscreteFormColetti;
  TDiscreteForm2D *DiscreteFormGL00Convolution;
  TDiscreteForm2D *DiscreteFormGL00AuxProblem;
  TDiscreteForm2D *DiscreteFormUpwind;
  TDiscreteForm2D *DiscreteFormSmagorinsky;
  TDiscreteForm2D *DiscreteFormVMSProjection;

  TDiscreteForm2D *DiscreteFormNLGalerkin;
  TDiscreteForm2D *DiscreteFormNLUpwind;
  TDiscreteForm2D *DiscreteFormNLColetti;
  TDiscreteForm2D *DiscreteFormNLGL00Convolution;
  TDiscreteForm2D *DiscreteFormNLGL00AuxProblem;
  TDiscreteForm2D *DiscreteFormNLSmagorinsky;
  TDiscreteForm2D *DiscreteFormNLVMSProjection;

  TDiscreteForm2D *DiscreteFormRHS;
  TDiscreteForm2D *DiscreteFormRHSColetti;
  TDiscreteForm2D *DiscreteFormRHSSmagorinskyExpl;
  TDiscreteForm2D *DiscreteFormMatrixGL00AuxProblem;
  TDiscreteForm2D *DiscreteFormGL00AuxProblemRHS;
  TDiscreteForm2D *DiscreteFormRHSLESModel;
  TDiscreteForm2D *DiscreteFormRHSAuxProblemU;
  TDiscreteForm2D *DiscreteFormMatrixAuxProblemU;
  
  // save the boundary condition
  BoundaryConditions[0] = BoundCond;
  BoundaryConditions[1] = BoundCond;  

  // save the boundary values  
  BoundaryValues[0] = U1BoundValue;
  BoundaryValues[1] = U2BoundValue;
 
  // save the nse bilinear coefficient   
  LinCoeffs[0] = lincoeffs;
  
  //default, i.e., velocity for nonlinear term
  NSEaux = aux;
  
  NSE_Rhsaux = new TAuxParam2D(1, 0, 0, 0, FeSpaces, NULL, NULL, NULL, NULL, 0, NULL);
  
  // aux for calculating the error
  NSEaux_error = nseaux_error;
  
  // set the Discreteforms
  InitializeDiscreteForms(DiscreteFormGalerkin,
                          DiscreteFormUpwind,
                          DiscreteFormSmagorinsky,
                          DiscreteFormColetti,
                          DiscreteFormGL00Convolution,
                          DiscreteFormGL00AuxProblem,
                          DiscreteFormVMSProjection,
                          DiscreteFormNLGalerkin,
                          DiscreteFormNLUpwind,
                          DiscreteFormNLSmagorinsky,
                          DiscreteFormNLColetti,
                          DiscreteFormNLGL00Convolution,
                          DiscreteFormNLGL00AuxProblem,
                          DiscreteFormNLVMSProjection,
                          DiscreteFormRHS,
                          DiscreteFormRHSColetti,
                          DiscreteFormRHSLESModel,
                          DiscreteFormMatrixGL00AuxProblem,
                          DiscreteFormGL00AuxProblemRHS,
                          DiscreteFormRHSSmagorinskyExpl,
                          DiscreteFormMatrixAuxProblemU,
                          DiscreteFormRHSAuxProblemU,
                          LinCoeffs[0], NSEType);
    
    
    // find discrete form
    switch(Disctype)
       {
          case GALERKIN:
            DiscreteFormARhs = DiscreteFormGalerkin;
            DiscreteFormNL = DiscreteFormNLGalerkin;
            DiscreteFormRhs = DiscreteFormRHS;
          break;

          case UPWIND:
            DiscreteFormARhs = DiscreteFormUpwind;
            DiscreteFormNL = DiscreteFormNLUpwind;    
            DiscreteFormRhs = DiscreteFormRHS;
            break;

//           case SMAGORINSKY:
//             DiscreteFormARhs = DiscreteFormSmagorinsky;
//             DiscreteFormNL = DiscreteFormNLSmagorinsky;              
//             break;

#ifdef __PRIVATE__ 
          case VMS_PROJECTION:
    
            DiscreteFormARhs = DiscreteFormVMSProjection;
            DiscreteFormNL = DiscreteFormNLVMSProjection;
            DiscreteFormRhs = DiscreteFormRHS;    
            break;
#endif
	    
          default:
            Error("Unknown DISCTYPE" << Disctype << endl);
            exit(-1);
        } 
     
     // set the discrete form for the Stokes equation
      if (TDatabase::ParamDB->PROBLEM_TYPE == STOKES)
       {
        DiscreteFormARhs = DiscreteFormUpwind;     
        DiscreteFormNL = NULL;
       }

#ifdef _SMPI                   
        if(Solver == DIRECT)
       {
           SQMATRICES[0] = SqmatrixM11;
           SQMATRICES[1] = SqmatrixM12;	  
           SQMATRICES[2] = SqmatrixM21;
           SQMATRICES[3] = SqmatrixM22;
           
           SQMATRICES[0]->Reset();
           SQMATRICES[1]->Reset();
           SQMATRICES[2]->Reset();
           SQMATRICES[3]->Reset();
           
           MATRICES[0] = MatrixB1;
           MATRICES[1] = MatrixB2;
           MATRICES[2] = MatrixB1T;
           MATRICES[3] = MatrixB2T;
           
           MATRICES[0]->Reset();
           MATRICES[1]->Reset();
           MATRICES[2]->Reset();
           MATRICES[3]->Reset();
           
           P_DS = new TSeqParDirectSolver(N_U,N_P,0,0,SQMATRICES,MATRICES);
       }
#endif

} // TSystemTNSE2D::Init

 
 
/* Assemble M, A and rhs */ 
void TSystemTNSE2D::Assemble(double *sol, double *rhs)
{
  int N_SquareMatrices, N_RectMatrices, N_Rhs, N_FESpaces;
  
  double *RHSs[3];

  TFESpace2D *fesprhs[3];
  
  N_Rhs = 2;
  N_FESpaces = 2;   
      
     // initialize matrices
     switch(NSEType)
      {
        case 1:
         SQMATRICES[0] = SqmatrixA11;
         SQMATRICES[1] = SqmatrixM11;
         MATRICES[0] = MatrixB1;
         MATRICES[1] = MatrixB2;

         SQMATRICES[0]->Reset();
         SQMATRICES[1]->Reset();
         MATRICES[0]->Reset();
         MATRICES[1]->Reset();

         N_SquareMatrices = 2;
         N_RectMatrices = 2;

         break;

        case 2:
          SQMATRICES[0] = SqmatrixA11;
          SQMATRICES[0] = SqmatrixM11;
          MATRICES[0] = MatrixB1;
          MATRICES[1] = MatrixB2;
          MATRICES[2] = MatrixB1T;
          MATRICES[3] = MatrixB2T;

          SQMATRICES[0]->Reset();
          SQMATRICES[1]->Reset();
          MATRICES[0]->Reset();
          MATRICES[1]->Reset();
          MATRICES[2]->Reset();
          MATRICES[3]->Reset();

          N_SquareMatrices = 2;
          N_RectMatrices = 4;
        break;

        case 3:
          SQMATRICES[0] = SqmatrixA11;
          SQMATRICES[1] = SqmatrixA12;
          SQMATRICES[2] = SqmatrixA21;
          SQMATRICES[3] = SqmatrixA22;
          SQMATRICES[4] = SqmatrixM11;
          SQMATRICES[5] = SqmatrixM22;
          MATRICES[0] = MatrixB1;
          MATRICES[1] = MatrixB2;

          SQMATRICES[0]->Reset();
          SQMATRICES[1]->Reset();
          SQMATRICES[2]->Reset();
          SQMATRICES[3]->Reset();
          SQMATRICES[4]->Reset();
          SQMATRICES[5]->Reset();
          MATRICES[0]->Reset();
          MATRICES[1]->Reset();

          N_SquareMatrices = 6;
          N_RectMatrices = 2;
  
#ifdef __PRIVATE__  
        if(Disctype == VMS_PROJECTION)
        {
          N_SquareMatrices = 7;
          SQMATRICES[6] =  MatricesL;
          SQMATRICES[6]->Reset();

          N_RectMatrices = 6;
          MATRICES[2] = Matrices_tilde_G11;
          MATRICES[3] = Matrices_tilde_G22;
          MATRICES[4] = Matrices_G11;
          MATRICES[5] = Matrices_G22;
          MATRICES[2]->Reset();
          MATRICES[3]->Reset();
          MATRICES[4]->Reset();
          MATRICES[5]->Reset();

          N_FESpaces = 4;
        }
#endif    
  
        break;

        case 4:
          SQMATRICES[0] = SqmatrixA11;
          SQMATRICES[1] = SqmatrixA12;
          SQMATRICES[2] = SqmatrixA21;
          SQMATRICES[3] = SqmatrixA22;
          SQMATRICES[4] = SqmatrixM11;
          SQMATRICES[5] = SqmatrixM22;
          MATRICES[0] = MatrixB1;
          MATRICES[1] = MatrixB2;
          MATRICES[2] = MatrixB1T;
          MATRICES[3] = MatrixB2T;

          SQMATRICES[0]->Reset();
          SQMATRICES[1]->Reset();
          SQMATRICES[2]->Reset();
          SQMATRICES[3]->Reset();
          SQMATRICES[4]->Reset();
          SQMATRICES[5]->Reset();
          MATRICES[0]->Reset();
          MATRICES[1]->Reset();
          MATRICES[2]->Reset();
          MATRICES[3]->Reset();

          N_SquareMatrices = 6;
          N_RectMatrices = 4;

#ifdef __PRIVATE__  
        if(Disctype == VMS_PROJECTION)
        {
          N_SquareMatrices = 7;
          SQMATRICES[6] =  MatricesL;
          SQMATRICES[6]->Reset();

          N_RectMatrices = 8;
          MATRICES[4] = Matrices_tilde_G11;
          MATRICES[5] = Matrices_tilde_G22;
          MATRICES[6] = Matrices_G11;
          MATRICES[7] = Matrices_G22;
          MATRICES[4]->Reset();
          MATRICES[5]->Reset();
          MATRICES[6]->Reset();
          MATRICES[7]->Reset();

          N_FESpaces = 4;
        }
#endif        
          break;
      } //  switch(NSEType)
     
      RHSs[0] = rhs;
      RHSs[1] = rhs + N_U;
      RHSs[2] = rhs + 2*N_U;
      memset(rhs, 0, (2*N_U+N_P)*SizeOfDouble);
     
      fesprhs[0] = FeSpaces[0];
      fesprhs[1] = FeSpaces[0];
      fesprhs[2] = FeSpaces[1];
  
      // assemble
      Assemble2D(N_FESpaces, FeSpaces,
        N_SquareMatrices, SQMATRICES,
        N_RectMatrices, MATRICES,
        N_Rhs, RHSs, fesprhs,
        DiscreteFormARhs,
        BoundaryConditions,
        BoundaryValues,
        NSEaux);
       
      if( (Disctype==UPWIND) && !(TDatabase::ParamDB->PROBLEM_TYPE == STOKES) )
       {
        switch(NSEType)
         {
          case 1:
          case 2:
            // do upwinding with one matrix
            UpwindForNavierStokes(LinCoeffs[0], SQMATRICES[0], FeFct[0], FeFct[1]);
            cout << "UPWINDING DONE : level " << endl;
            break;

          case 3:
          case 4:
            // do upwinding with two matrices
            cout << "UPWINDING DONE : level " << endl;
            UpwindForNavierStokes(LinCoeffs[0], SQMATRICES[0], FeFct[0], FeFct[1]);
            UpwindForNavierStokes(LinCoeffs[0],SQMATRICES[3], FeFct[0], FeFct[1]);
           break;
         }                        // endswitch
       }                          // endif     
            
      // slip with boundary condition
      if (TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION >= 1)
      {
        if(NSEType <4)
         {
          OutPut("For slip with friction bc NSTYPE 4 is ");
          OutPut("necessary !!!!! " << endl);
          exit(4711);
         }
 
        N_FESpaces = 1;
        N_SquareMatrices = 8;
        N_RectMatrices = 2;
        N_Rhs = 2;

        SQMATRICES[0] = SqmatrixA11;
        SQMATRICES[1] = SqmatrixA22;
        SQMATRICES[2] = SqmatrixA12;
        SQMATRICES[3] = SqmatrixA21;
        SQMATRICES[4] = SqmatrixM11;
        SQMATRICES[5] = SqmatrixM22;
        SQMATRICES[6] = SqmatrixM12;
        SQMATRICES[7] = SqmatrixM21;

        MATRICES[0] = MatrixB1T;
        MATRICES[1] = MatrixB2T;


        Assemble2DSlipBC(N_FESpaces, FeSpaces,
                         N_SquareMatrices, SQMATRICES,
                         N_RectMatrices, MATRICES,
                         N_Rhs, RHSs, fesprhs,
                         NULL,
                         BoundaryConditions,
                         BoundaryValues,
                         NSEaux,
                         FeFct[0], FeFct[1]);

      }// (TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION >=      
    
     // set rhs for Dirichlet nodes
     memcpy(sol+N_Active, rhs+N_Active, N_DirichletDof*SizeOfDouble);
     memcpy(sol+N_U+N_Active, rhs+N_U+N_Active, N_DirichletDof*SizeOfDouble); 
     
#ifdef __PRIVATE__   
      // update matrices
      if (Disctype == VMS_PROJECTION)
        {
         SQMATRICES[0] = SqmatrixA11;
         SQMATRICES[1] = SqmatrixA12;
         SQMATRICES[2] = SqmatrixA21;
         SQMATRICES[3] = SqmatrixA22;
         SQMATRICES[6] =  MatricesL;
         MATRICES[2] = Matrices_tilde_G11;
         MATRICES[3] = Matrices_tilde_G22;
         MATRICES[4] = Matrices_G11;
         MATRICES[5] = Matrices_G22;

         VMSProjectionUpdateMatrices(N_U, FeSpaces[0]->GetActiveBound(), FeSpaces[3]->GetN_DegreesOfFreedom(),
                                     SQMATRICES, MATRICES);
        }
#endif  
     
} // TSystemTNSE2D::Assemble(T

 
void TSystemTNSE2D::AssembleRhs(double *sol, double *rhs)
{
  int N_SquareMatrices, N_RectMatrices, N_Rhs, N_FESpaces;
  
  double *RHSs[3];
  
  TFESpace2D *fesprhs[3];
  
  
      N_Rhs = 2;
      N_FESpaces = 1;
      N_SquareMatrices = 0;
      N_RectMatrices = 0;
  
      RHSs[0] = rhs;
      RHSs[1] = rhs + N_U;

      memset(rhs, 0, (2*N_U+N_P)*SizeOfDouble);
     
      fesprhs[0] = FeSpaces[0];
      fesprhs[1] = FeSpaces[0];
      fesprhs[2] = FeSpaces[1];  
  
      // assemble
      Assemble2D(N_FESpaces, FeSpaces,
        N_SquareMatrices, NULL,
        N_RectMatrices, NULL,
        N_Rhs, RHSs, fesprhs,
        DiscreteFormRhs,
        BoundaryConditions,
        BoundaryValues,
        NSE_Rhsaux);      
}



void TSystemTNSE2D::AssembleSystMat(double scale, double *oldrhs, double *rhs, double *sol)
{
 double tau, val = TDatabase::TimeDB->SCALE_DIVERGENCE_CONSTRAINT;
  
  tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
     
  memset(B, 0, (2*N_U+N_P)*SizeOfDouble);    

  // old rhs multiplied with current subtime step and theta3 on B
  Daxpy(N_Active, tau*TDatabase::TimeDB->THETA3, oldrhs, B);
  Daxpy(N_Active, tau*TDatabase::TimeDB->THETA3, oldrhs+N_U, B+N_U);   

  // add rhs from current sub time step to rhs array B
  Daxpy(N_Active, tau*TDatabase::TimeDB->THETA4, rhs, B);
  Daxpy(N_Active, tau*TDatabase::TimeDB->THETA4, rhs+N_U, B+N_U);   
  
  // scale the pressure matrices, not in nonlinear step 
  if (scale != 1.0)
  {
   switch(NSEType)
    {
     case 1:
     case 3:
        Dscal(MatrixB1->GetN_Entries(), scale, MatrixB1->GetEntries());
        Dscal(MatrixB2->GetN_Entries(), scale, MatrixB2->GetEntries());
     break;

    case 2:
    case 4:     
         Dscal(MatrixB1T->GetN_Entries(), scale, MatrixB1T->GetEntries());
         Dscal(MatrixB2T->GetN_Entries(), scale, MatrixB2T->GetEntries());
	 
      // scale divergence constraint
      if(val>0) 
       {
        Dscal(MatrixB1->GetN_Entries(), val*scale, MatrixB1->GetEntries());
        Dscal(MatrixB2->GetN_Entries(), val*scale, MatrixB2->GetEntries());
       }
      break;
    } // switch(NSETyp
  } //  if (scale != 1.0)
    
   // Also currently : M := M + gamma A
   // M = M + (-gamma - tau*TDatabase::TimeDB->THETA2) A 
   // defect = M * sol
   // B:= B + defect 
   memset(defect, 0, (2*N_U+N_P)*SizeOfDouble);  
   switch(NSEType)
    {
     case 1:
     case 2:
       MatAdd(SqmatrixM11, SqmatrixA11, -tau*TDatabase::TimeDB->THETA2);          
       gamma = - tau*TDatabase::TimeDB->THETA2;
   
       MatVectActive(SqmatrixM11, sol, defect);
       MatVectActive(SqmatrixM11, sol+N_U, defect+N_U);
       Daxpy(N_Active, 1, defect, B);
       Daxpy(N_Active, 1, defect+N_U, B+N_U);
 
       // assembling of system matrix       
       MatAdd(SqmatrixM11, SqmatrixA11, -gamma + tau*TDatabase::TimeDB->THETA1);   
       gamma = tau*TDatabase::TimeDB->THETA1;
     break;

     case 3:
     case 4:
       if(TDatabase::ParamDB->P0)
       {
//        printall(SqmatrixA11,(char*)"SqmatrixA11",SqmatrixM11,(char*)"SqmatrixM11");  
//        cout << " MHD_K :" << MHD_K/TDatabase::ParamDB->RE_NR << endl;           
            MatAdd(SqmatrixA11,SqmatrixM11 , MHD_K/TDatabase::ParamDB->RE_NR);
            MatAdd(SqmatrixA12,SqmatrixM12 , MHD_K/TDatabase::ParamDB->RE_NR);
            MatAdd(SqmatrixA21,SqmatrixM21 , MHD_K/TDatabase::ParamDB->RE_NR);
            MatAdd(SqmatrixA22,SqmatrixM21 , MHD_K/TDatabase::ParamDB->RE_NR);

            // printall(SqmatrixA11,(char*)"SqmatrixA11",SqmatrixM11,(char*)"SqmatrixM11"); 
       
//       cout << Ha << " " << TDatabase::ParamDB->RE_NR << "  " << MHD_K/TDatabase::ParamDB->RE_NR << endl;
       }

       MatAdd(SqmatrixM11, SqmatrixA11, - tau*TDatabase::TimeDB->THETA2);
       MatAdd(SqmatrixM12, SqmatrixA12, - tau*TDatabase::TimeDB->THETA2);
       MatAdd(SqmatrixM21, SqmatrixA21, - tau*TDatabase::TimeDB->THETA2);
       MatAdd(SqmatrixM22, SqmatrixA22, - tau*TDatabase::TimeDB->THETA2);       
       gamma = - tau*TDatabase::TimeDB->THETA2;

       MatVectActive(SqmatrixM11, sol, defect);
       Daxpy(N_Active, 1, defect, B);
       MatVectActive(SqmatrixM12, sol+N_U, defect);
       Daxpy(N_Active, 1, defect, B);
       MatVectActive(SqmatrixM21, sol, defect+N_U);
       Daxpy(N_Active, 1, defect+N_U, B+N_U);
       MatVectActive(SqmatrixM22, sol+N_U, defect+N_U);
       Daxpy(N_Active, 1, defect+N_U, B+N_U);

       //assembling system matrix
       MatAdd(SqmatrixM11, SqmatrixA11, -gamma + tau*TDatabase::TimeDB->THETA1);
       MatAdd(SqmatrixM12, SqmatrixA12, -gamma + tau*TDatabase::TimeDB->THETA1);
       MatAdd(SqmatrixM21, SqmatrixA21, -gamma + tau*TDatabase::TimeDB->THETA1);
       MatAdd(SqmatrixM22, SqmatrixA22, -gamma + tau*TDatabase::TimeDB->THETA1);       
       gamma = tau*TDatabase::TimeDB->THETA1;     
     break;     
    } 
  
   // set rhs for Dirichlet nodes
   memcpy(B+N_Active, rhs+N_Active, N_DirichletDof*SizeOfDouble);
   memcpy(B+N_U+N_Active, rhs+N_U+N_Active, N_DirichletDof*SizeOfDouble); 
  
   
   SystMatAssembled  = TRUE;
} // AssembleSystMat


/* assemble only LHS, not rhs */
void TSystemTNSE2D::AssembleSystMatNonLinear()
{
 double tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
     
  if(SystMatAssembled)
   {
    cout << "Restore System mat before calling AssembleSystMat" <<endl;
   }

   switch(NSEType)
    {
     case 1:
     case 2: 
       // assembling of system matrix       
       MatAdd(SqmatrixM11, SqmatrixA11,   tau*TDatabase::TimeDB->THETA1);   
       gamma = tau*TDatabase::TimeDB->THETA1;
     break;

     case 3:
     case 4:       
       //assembling system matrix
//       cout <<   
         
       MatAdd(SqmatrixM11, SqmatrixA11, tau*TDatabase::TimeDB->THETA1);
       MatAdd(SqmatrixM12, SqmatrixA12,  tau*TDatabase::TimeDB->THETA1);
       MatAdd(SqmatrixM21, SqmatrixA21,  tau*TDatabase::TimeDB->THETA1);
       MatAdd(SqmatrixM22, SqmatrixA22,  tau*TDatabase::TimeDB->THETA1);       
       gamma = tau*TDatabase::TimeDB->THETA1;     
     break;
    } 

   SystMatAssembled  = TRUE;
} // AssembleSystMatNonLinear




void TSystemTNSE2D::RestoreMassMat()
{

//   cout << "RestoreMassMat  gamma " << gamma << endl;
  if(SystMatAssembled)
   {
    // restore the mass matrix
    switch(NSEType)
     {
      case 1:
      case 2:
       MatAdd(SqmatrixM11, SqmatrixA11, -gamma);          
       gamma = 0.;
      break;

     case 3:
     case 4:
       //assembling system matrix
       MatAdd(SqmatrixM11, SqmatrixA11, -gamma);
       MatAdd(SqmatrixM12, SqmatrixA12, -gamma);
       MatAdd(SqmatrixM21, SqmatrixA21, -gamma);
       MatAdd(SqmatrixM22, SqmatrixA22, -gamma);       
       gamma = 0.;     
     break;
    } 
    
    SystMatAssembled  = FALSE;  
   }
  else
  {
    cout << "System is not assembled to restore " <<endl;
  }
   
//   cout << "RestoreMassMat" << endl;
//   exit(0);
  
}


void TSystemTNSE2D::AssembleANonLinear(double *sol, double *rhs)
{
 int N_SquareMatrices, N_RectMatrices, N_Rhs, N_FESpaces, last_sq;

     N_RectMatrices = 0;          
     N_Rhs = 0;
     N_FESpaces = 1;
 
     // set the nonliner matrices
      switch(TDatabase::ParamDB->NSTYPE)
       {
        case 1:
        case 2:
          SQMATRICES[0] = SqmatrixA11;
          SQMATRICES[0]->Reset();

          N_SquareMatrices = 1;
        break;

        case 3:
        case 4:
          if (TDatabase::ParamDB->SC_NONLIN_ITE_TYPE_SADDLE==0)
           {
            SQMATRICES[0] = SqmatrixA11;
            SQMATRICES[1] = SqmatrixA22;
            SQMATRICES[0]->Reset();
            SQMATRICES[1]->Reset();

            N_SquareMatrices = 2;
            last_sq = 1;

#ifdef __PRIVATE__   
            if (Disctype == VMS_PROJECTION)
              {
               SQMATRICES[0] = SqmatrixA11;
               SQMATRICES[1] = SqmatrixA12;
               SQMATRICES[2] = SqmatrixA21;
               SQMATRICES[3] = SqmatrixA22;
               SQMATRICES[0]->Reset();
               SQMATRICES[1]->Reset();
               SQMATRICES[2]->Reset();
               SQMATRICES[3]->Reset();

               N_SquareMatrices = 4;
               last_sq = 3;

               N_RectMatrices = 2;
               MATRICES[0] = Matrices_tilde_G11;
               MATRICES[1] = Matrices_tilde_G22;
               MATRICES[0]->Reset();
               MATRICES[1]->Reset();
       
               N_FESpaces = 4;
              }  
#endif

           }
          else
           {
            // Newton method
            cout<< "Newton method not tested " <<endl;
            exit(0);
           }

         break;
        } // switch(TDatabase::ParamDB->NSTYPE)
         
    
      // assemble the nonlinear part of NSE
      Assemble2D(N_FESpaces, FeSpaces,
                 N_SquareMatrices, SQMATRICES,
                 N_RectMatrices, MATRICES,
                 N_Rhs, NULL, NULL,
                 DiscreteFormNL,
                 BoundaryConditions,
                 BoundaryValues,
                 NSEaux);    

       // apply upwind disc
      if( (Disctype==UPWIND) && !(TDatabase::ParamDB->PROBLEM_TYPE == STOKES) )
       {
        switch(NSEType)
         {
          case 1:
          case 2:
            // do upwinding with one matrix
            UpwindForNavierStokes(LinCoeffs[0], SQMATRICES[0], FeFct[0], FeFct[1]);
            cout << "UPWINDING DONE : level " << endl;
            break;

          case 3:
          case 4:
            // do upwinding with two matrices
            cout << "UPWINDING DONE : level " << endl;
            UpwindForNavierStokes(LinCoeffs[0], SQMATRICES[0], FeFct[0], FeFct[1]);
            UpwindForNavierStokes(LinCoeffs[0],SQMATRICES[last_sq], FeFct[0], FeFct[1]);
          break;
         }                        // endswitch
       }                          // endif     
       
       
      // slip with boundary condition
      if (TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION >= 1)
      { 
        N_FESpaces = 1;
        N_SquareMatrices = 2;
        N_RectMatrices = 0;
        N_Rhs = 0;

        SQMATRICES[0] = SqmatrixA11;
        SQMATRICES[1] = SqmatrixA22;

        Assemble2DSlipBC(N_FESpaces, FeSpaces,
                         N_SquareMatrices, SQMATRICES,
                         N_RectMatrices, NULL,
                         N_Rhs, NULL, NULL,
                         NULL,
                         BoundaryConditions,
                         BoundaryValues,
                         NSEaux,
                         FeFct[0], FeFct[1]);

      }// (TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION >=         

     // set rhs for Dirichlet nodes
     memcpy(sol+N_Active, rhs+N_Active, N_DirichletDof*SizeOfDouble);
     memcpy(sol+N_U+N_Active, rhs+N_U+N_Active, N_DirichletDof*SizeOfDouble); 
     
     

     
     
     
     
#ifdef __PRIVATE__   
      // update matrices
      if (Disctype == VMS_PROJECTION)
        {
         SQMATRICES[0] = SqmatrixA11;
         SQMATRICES[1] = SqmatrixA12;
         SQMATRICES[2] = SqmatrixA21;
         SQMATRICES[3] = SqmatrixA22;
         SQMATRICES[6] =  MatricesL;
         MATRICES[2] = Matrices_tilde_G11;
         MATRICES[3] = Matrices_tilde_G22;
         MATRICES[4] = Matrices_G11;
         MATRICES[5] = Matrices_G22;

         VMSProjectionUpdateMatrices(N_U, FeSpaces[0]->GetActiveBound(), FeSpaces[3]->GetN_DegreesOfFreedom(),
                                     SQMATRICES, MATRICES);
        }
#endif       
     
} //TSystemTNSE2D::AssembleNonLinear(

 
void TSystemTNSE2D::Solve(double *sol)
{  
  if(!SystMatAssembled)
  {
    cout << "System Matrix is not assembled to solve " <<endl;
    exit(0);
  }
  
    switch(Solver)
     {
      case AMG_SOLVE:
        cout << "AMG_SOLVE not yet implemented " <<endl;
      break;

      case GMG:
        cout << "GMG solver not yet implemented " <<endl;
      break;

      case DIRECT:
        switch(NSEType)
         {
          case 1:
            DirectSolver(SqmatrixM11, MatrixB1,  MatrixB2, B, sol);
          break;

          case 2:
             DirectSolver(SqmatrixA11, MatrixB1T, MatrixB2T, MatrixB1,  MatrixB2, B, sol);
          break;

          case 3:
           cout << "Direct solver not yet implemented for NSTYPE 3 " <<endl;
          break;

          case 4:
#ifdef _SMPI	
              P_DS->Solve(sol,B,true);
#else
              DirectSolver(SqmatrixM11, SqmatrixM12, SqmatrixM21, SqmatrixM22, 
                          MatrixB1T, MatrixB2T, MatrixB1,  MatrixB2, B, sol); 
#endif
	   
          break;
      } //  switch(NSEType) 

      break;      
 
      default:
            OutPut("Unknown Solver" << endl);
            exit(4711);;
     }    

}



void TSystemTNSE2D::GetTNSEResidual(double *sol, double *res)
{
  
  if(!SystMatAssembled)
   {
    cout << "System Matrix is not assembled to calculate residual " <<endl;
    exit(0);
   }
           

     switch(NSEType)
      {
        case 1:
          SQMATRICES[0] = SqmatrixM11;

          MATRICES[0] = MatrixB1;
          MATRICES[1] = MatrixB2;  
         break;

        case 2:
          SQMATRICES[0] = SqmatrixM11;

          MATRICES[0] = MatrixB1;
          MATRICES[1] = MatrixB2;
          MATRICES[2] = MatrixB1T;
          MATRICES[3] = MatrixB2T;
        break;

        case 3:
          SQMATRICES[0] = SqmatrixM11;
          SQMATRICES[1] = SqmatrixM12;
          SQMATRICES[2] = SqmatrixM21;
          SQMATRICES[3] = SqmatrixM22;

          MATRICES[0] = MatrixB1;
          MATRICES[1] = MatrixB2;  
  
        break;

        case 4:
          SQMATRICES[0] = SqmatrixM11;
          SQMATRICES[1] = SqmatrixM12;
          SQMATRICES[2] = SqmatrixM21;
          SQMATRICES[3] = SqmatrixM22;

          MATRICES[0] = MatrixB1;
          MATRICES[1] = MatrixB2;
          MATRICES[2] = MatrixB1T;
          MATRICES[3] = MatrixB2T;
 
        break;
      } //  switch(NSEType)
          
   Defect(sqmatrices, matrices, sol, B, res); 

// double *defect = new double[2*N_U+N_P];
//    memset(defect, 0, 2*N_U+N_P*SizeOfDouble);
//    memset(res, 0, 2*N_U+N_P*SizeOfDouble);
//    
//    // Block A
//    MatVect(SqmatrixM11, sol, defect);
//    Daxpy(N_U, 1, defect, res);
//    MatVectActive(SqmatrixM12, sol+N_U, defect);
//    Daxpy(N_Active, 1, defect, res);
//    MatVectActive(SqmatrixM21, sol, defect+N_U);
//    Daxpy(N_Active, 1, defect+N_U, res+N_U);
//    MatVect(SqmatrixM22, sol+N_U, defect+N_U);
//    Daxpy(N_U, 1, defect+N_U, res+N_U);
// 
//       // Block BT
//       MatVect1(MatrixB1T, sol+2*N_U, defect);
//       Daxpy(N_Active, 1, defect, res); 
//       MatVect1(MatrixB2T, sol+2*N_U, defect+N_U);
//       Daxpy(N_Active, 1, defect+N_U, res+N_U); 
//       
//       // Block B
//       MatVect1(MatrixB1, sol, defect+2*N_U);
//       Daxpy(N_P, 1, defect+2*N_U, res+2*N_U); 
//       MatVect1(MatrixB2, sol+N_U, defect+2*N_U);
//       Daxpy(N_P, 1, defect+2*N_U, res+2*N_U); 
// 
//       Daxpy(2*N_U+N_P, -1., B, res);

   
} // TSystemTNSE2D::GetResidual


void TSystemTNSE2D::MeasureTNSEErrors(DoubleFunct2D *ExactU1, DoubleFunct2D *ExactU2, DoubleFunct2D *ExactP,
                                    double *AllErrors)
{
  MultiIndex2D TimeNSAllDerivatives[3] = { D00, D10, D01 };

  double errors[4],  u_error[4];    
    
     // errors in first velocity component
     FeFct[0]->GetErrors(ExactU1, 3, TimeNSAllDerivatives, 2,
                         L2H1Errors,
                         NULL, NSEaux_error, 1, FeSpaces, errors);
      u_error[0] = errors[0];
      u_error[1] = errors[1];
      
     // errors in second velocity component
     FeFct[1]->GetErrors(ExactU2, 3, TimeNSAllDerivatives, 2,
                         L2H1Errors,
                         NULL, NSEaux_error, 1, FeSpaces, errors);
     u_error[2] = errors[0];
     u_error[3] = errors[1];      
      
      // errors in pressure
     FeFct[2]->GetErrors(ExactP, 3, TimeNSAllDerivatives, 2,
                         L2H1Errors,
                         NULL, NSEaux_error, 1, FeSpaces+1, errors);     

     
     // calculate all errors
     AllErrors[0] = sqrt(u_error[0]*u_error[0]+u_error[2]*u_error[2]);
     AllErrors[1] = sqrt(u_error[1]*u_error[1]+u_error[3]*u_error[3]);
     AllErrors[2] = errors[0];
     AllErrors[3] = errors[1];    
     
      // error in L^infty(0,t,L^2)
      if(AllErrors[0] > AllErrors[5])
       {
        AllErrors[5]  = AllErrors[0];
        AllErrors[4]  =  TDatabase::TimeDB->CURRENTTIME;
      }

      
      // error in L^2(0,t,L^2)    
      AllErrors[6] += (u_error[0]*u_error[0] + u_error[2]*u_error[2] +olderror_l_2_l_2u)*TDatabase::TimeDB->TIMESTEPLENGTH/2.0;      
      olderror_l_2_l_2u = u_error[0]*u_error[0] + u_error[2]*u_error[2];
     
}

// ************************************************************************************************************************ //

void printall( TSquareMatrix2D * MAT1, char *MAT1_name
//            ,TSquareMatrix2D * A12 
//            ,TSquareMatrix2D * A21
              ,TSquareMatrix2D * MAT2, char *MAT2_name
//            ,TSquareMatrix2D * A23
//            ,TSquareMatrix2D * A31
//            ,TSquareMatrix2D * A32
//            ,TMatrix* B1T
//            ,TMatrix* B2T
//            ,TMatrix* B1
//            ,TMatrix* B2,
                    )
{
    double * entries;
    
       cout << "\n*******************************************************************************\n";
       cout << MAT1_name << "\t";
       entries  = MAT1->GetEntries();
       for(int ii =0; ii< 20 ; ii++ )
           cout << entries[ii] << "  " ;
       
       cout << endl;
       cout << MAT2_name << "\t";
       entries  = MAT2->GetEntries();       
       for(int ii =0; ii< 20 ; ii++ )
           cout << entries[ii] << "  " ;
       
       cout << endl;
 }


void TSystemTNSE2D::printall()
{

	TSquareMatrix2D *MAT1 = SqmatrixA21;
	TSquareMatrix2D *MAT2 = SqmatrixA22;
	TSquareMatrix2D *MAT3 = SqmatrixM21;
	TSquareMatrix2D *MAT4 = SqmatrixM22;
	cout << "\n*******************************************************************************\n";
	cout << " MATRIX : A21 "
		 << "\t";
	double *entries = MAT1->GetEntries();
	int *rowPtr = MAT1->GetRowPtr();
	int *colIndex = MAT1->GetKCol();

	int i = 9072; // HARDCODED , PRINTING THE 482th row.

	for (int index = rowPtr[i]; index < rowPtr[i + 1]; index++)
	{
		if (colIndex[index] == i)
			cout << " ** " << entries[index] << " **"
				 << "\t";
		else
			cout << entries[index] << "\t";
	}

	cout << "\n*******************************************************************************\n";
	cout << " MATRIX : A22 "
		 << "\t";
	entries = MAT2->GetEntries();
	rowPtr = MAT2->GetRowPtr();
	colIndex = MAT2->GetKCol();

	for (int index = rowPtr[i]; index < rowPtr[i + 1]; index++)
	{
		if (colIndex[index] == i)
			cout << " ** " << entries[index] << " **"
				 << "\t";
		else
			cout << entries[index] << "\t";
	}

	cout << "\n*******************************************************************************\n";
	cout << " MATRIX : M21 "
		 << "\t";
	entries = MAT3->GetEntries();
	rowPtr = MAT3->GetRowPtr();
	colIndex = MAT3->GetKCol();

	for (int index = rowPtr[i]; index < rowPtr[i + 1]; index++)
	{
		if (colIndex[index] == i)
			cout << " ** " << entries[index] << " **"
				 << "\t";
		else
			cout << entries[index] << "\t";
	}

	cout << "\n*******************************************************************************\n";
	cout << " MATRIX : M22 "
		 << "\t";
	entries = MAT4->GetEntries();
	rowPtr = MAT4->GetRowPtr();
	colIndex = MAT4->GetKCol();

	for (int index = rowPtr[i]; index < rowPtr[i + 1]; index++)
	{
		if (colIndex[index] == i)
			cout << " ** " << entries[index] << " **"
				 << "\t";
		else
			cout << entries[index] << "\t";
	}

	TMatrix2D *MAT5 = MatrixB1T;
	TMatrix2D *MAT6 = MatrixB2T;

	cout << "\n*******************************************************************************\n";
	cout << " MATRIX : B1T "
		 << "\t";
	entries = MAT5->GetEntries();
	rowPtr = MAT5->GetRowPtr();
	colIndex = MAT5->GetKCol();

	for (int index = rowPtr[i]; index < rowPtr[i + 1]; index++)
	{
		if (colIndex[index] == i)
			cout << " ** " << entries[index] << " **"
				 << "\t";
		else
			cout << entries[index] << "\t";
	}

	cout << "\n*******************************************************************************\n";
	cout << " MATRIX : B2T "
		 << "\t";
	entries = MAT6->GetEntries();
	rowPtr = MAT6->GetRowPtr();
	colIndex = MAT6->GetKCol();

	for (int index = rowPtr[i]; index < rowPtr[i + 1]; index++)
	{
		if (colIndex[index] == i)
			cout << " ** " << entries[index] << " **"
				 << "\t";
		else
			cout << entries[index] << "\t";
	}

	cout << "\n*******************************************************************************\n";
	cout << " RHS "
		 << "\t";
	TFEVectFunct2D *vel = VelocityFct;
	int len = VelocityFct->GetComponent(0)->GetLength();
	cout << B[i + len] << endl;

	cout << "\n################################################################################################";
	cout << "\n################################################################################################\n";
}

void TSystemTNSE2D::PickFreeSlipDOFs(std::vector<int> bdid)
{
	TFESpace2D *fespace = FeSpaces[0];
	TCollection *coll = fespace->GetCollection();
	int N_Cells = coll->GetN_Cells();
	int *GlobalNumbers = fespace->GetGlobalNumbers();
	int *BeginIndex = fespace->GetBeginIndex();

	int dof =  fespace->GetN_DegreesOfFreedom();
	int ActiveDOF =  fespace->GetActiveBound();

	cout << " DOF : " << dof << endl;
	cout << " Act DOF : " << ActiveDOF << endl;

	double *xx = new double[dof];
	double *yy = new double[dof];

	for ( int i = 0 ; i < dof ; i++)
		fespace->GetDOFPosition(i,xx[i],yy[i]);

	cout << " position Obtained " << endl;




	for (int cellNo = 0; cellNo < N_Cells; cellNo++)
	{
		TBaseCell *currentCell = coll->GetCell(cellNo);
		FE2D currentElement = fespace->GetFE2D(cellNo, currentCell);

		int N_Joints = currentCell->GetN_Edges();
		TFE2D *ele = TFEDatabase2D::GetFE2D(currentElement);
		TFEDesc2D *FEDesc_Obj = ele->GetFEDesc2D();
		int N_EdgeDOF = FEDesc_Obj->GetN_JointDOF();

		int *GlobalDOF = GlobalNumbers + BeginIndex[cellNo];

		// for ( int i = 0 ; i < 9 ; i++)
		// {
		// 	cout << "i :" << i  << "  " << GlobalDOF[i] <<endl;
		// }

		for (int jointId = 0; jointId < N_Joints; jointId++)
		{
			double x0, y0, x1, y1, nx, ny, tx, ty, hE;
			currentCell->GetVertex(jointId)->GetCoords(x0, y0);
			currentCell->GetVertex((jointId + 1) % N_Joints)->GetCoords(x1, y1);

			// compute length of the boundary edge
			hE = sqrt((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0));
			// compute normal vector to this boundary (normalized)
			nx = (y1 - y0) / hE;
			ny = (x0 - x1) / hE;
			// tangential normal vector to this boundary (normalized)
			tx = (x1 - x0) / hE;
			ty = (y1 - y0) / hE;

			TJoint *joint = currentCell->GetJoint(jointId);

			// Filter out joints on Boundary Edges
			if (!(joint->GetType() == BoundaryEdge) )
				continue;

			TBoundEdge *boundedge = (TBoundEdge *)joint;
			TBoundComp *BoundComp = boundedge->GetBoundComp();

			int comp = BoundComp->GetID();

			

			// If the current BD id not is in the input vector, then terminate the loop
			if (std::find(bdid.begin(), bdid.end(), comp) == bdid.end())
				continue;
			cout << comp << endl;

			int *JointDOF = FEDesc_Obj->GetJointDOF(jointId);

			// for (int i = 0 ; i < 9 ; i++)
			// 		cout << GlobalDOF[i] << "\t";
			// 	cout << endl;

			for (int vert = 0; vert < FEDesc_Obj->GetN_JointDOF(); vert++)
			{
				

				int glob_vertex_no = GlobalDOF[JointDOF[vert]];
				// cout << " DOF IDENTIFIED : " << glob_vertex_no <<endl;

				if(std::find(freeslipDOFs.begin(),freeslipDOFs.end(),glob_vertex_no) != freeslipDOFs.end() )
					continue;
				if(glob_vertex_no >= ActiveDOF)
					continue;
				
				freeslipDOFs.emplace_back(glob_vertex_no);
				freeslipNormal_n1.emplace_back(nx);
				freeslipNormal_n2.emplace_back(ny);
			}
		}
	}

	// Remove the duplicate elements from the vector

	cout << "number of Freeslip DOF picked : " << freeslipDOFs.size() << endl;
	cout << "number of Freeslip nx       : " << freeslipNormal_n1.size() << endl;
	cout << "number of Freeslip ny     : " << freeslipNormal_n2.size() << endl;
	
	std::for_each(freeslipDOFs.begin(),freeslipDOFs.end(),[&](int i){ cout << i << " - " << xx[i] << " " << yy[i] <<endl;  } );
	cout<<endl;
}

void TSystemTNSE2D::modifyMatrixFreeSlip()
{
	TSquareMatrix2D *diagSquare, *antiDiagSquare;
	TMatrix2D *BMatrix;

	int len = FeSpaces[0]->GetN_DegreesOfFreedom();
	int Active = FeSpaces[0]->GetActiveBound();

	for (int i = 0; i < freeslipDOFs.size(); i++)
	{
		
		int DOF = freeslipDOFs[i];
		double nx = freeslipNormal_n1[i];
		double ny = freeslipNormal_n2[i];


		// Ignore Dirichlet boundaries
		if (DOF > Active)
			continue;


		double n1, n2;
		int disp;
		int comp;

		if (fabs(nx) < 1e-8)
			comp = 1;
		else if (fabs(ny) < 1e-8)
			comp = 0;
		else
			comp = 1; // Default assign the values to 2nd component



		if (comp == 0)
		{
			diagSquare = SqmatrixM11;
			antiDiagSquare = SqmatrixM12;
			BMatrix = MatrixB1T;
			n1 = nx;
			n2 = ny;
			disp = 0;
		}
		else
		{
			antiDiagSquare = SqmatrixM21;
			diagSquare = SqmatrixM22;
			BMatrix = MatrixB2T;
			n1 = ny;
			n2 = nx;
			disp = len;
		}

		// Now make the diagonal matrix with diagonal entry as 1 and the antidiagonal matrix and the b matrix as zero

		int *RowPtr_diag = diagSquare->GetRowPtr();
		int *ColPtr_diag = diagSquare->GetKCol();
		double *entries_diag = diagSquare->GetEntries();

		int *RowPtr_Antidiag = antiDiagSquare->GetRowPtr();
		int *ColPtr_Antidiag = antiDiagSquare->GetKCol();
		double *entries_Antidiag = antiDiagSquare->GetEntries();

		int *RowPtr_B = BMatrix->GetRowPtr();
		int *ColPtr_B = BMatrix->GetKCol();
		double *entries_B = BMatrix->GetEntries();

		int begin = RowPtr_diag[DOF];
		int end = RowPtr_diag[DOF + 1];

		for (int row = begin; row < end; row++)
		{
			if (ColPtr_diag[row] == DOF)
				entries_diag[row] = n1;
			else
				entries_diag[row] = 0;
		}

		begin = RowPtr_Antidiag[DOF];
		end = RowPtr_Antidiag[DOF + 1];

		for (int row = begin; row < end; row++)
		{
			if (ColPtr_diag[row] == DOF)
				entries_Antidiag[row] = n2;
			else
				entries_Antidiag[row] = 0;
		}

		begin = RowPtr_B[DOF];
		end = RowPtr_B[DOF + 1];
		for (int row = begin; row < end; row++)
			entries_B[row] = 0;

		B[DOF + disp] = 0.0;
	}
}



#endif // #ifdef __2D__
