// =======================================================================
//
// Purpose:     main program (TCD with direct solver)
//
// Author:      Sashikumaar Ganesan
//
// =======================================================================
#ifdef _MPI
#  include "mpi.h"
#endif

#include <Domain.h>
#include <Database.h>
#include <FEDatabase3D.h>
#include <FESpace3D.h>
#include <SquareStructure3D.h>
#include <Structure3D.h>
#include <AuxParam3D.h>
#include <QuadAffin.h>
#include <DirectSolver.h>
#include <Assemble3D.h>
#include <Output3D.h>
#include <DiscreteForm3D.h>
#include <LinAlg.h>
#include <Collection.h>
#include <Upwind3D.h>
#include <TCD3D.h>
#include <FEM_TVD_FCT.h>
#include <TetraAffin.h>
#include <TetraIsoparametric.h>
#include <HexaAffin.h>
#include <HexaTrilinear.h>
#include <HexaIsoparametric.h>

#include <string.h>
#include <sstream>
#include <MooNMD_Io.h>
#include <stdlib.h>
// #include <malloc.h>

#include <ItMethod.h>
#include <FixedPointIte.h>
#include <MultiGridScaIte.h>
#include <FgmresIte.h>

#include <MultiGrid3D.h>
#include <MGLevel3D.h>

#include <MainUtilities.h>
#include <TimeUtilities.h>

#include <tetgen.h>
#include <GridCell.h>
#include <MacroCell.h>
#include <BdPlane.h>
#include <BdSphere.h>
#include <IsoBoundFace.h>
#include <InterfaceJoint3D.h>
#include <IsoInterfaceJoint3D.h>


#ifdef _MPI
#include <MeshPartition.h>
#include <ParFECommunicator3D.h>
#include <MumpsSolver.h>
#include <ParVector3D.h>
#include <ParVectorNSE3D.h>
#include <Scalar_ParSolver.h>
#endif

#define AMG 0
#define GMG 1
double bound = 0;

// =======================================================================
// include current example
// =======================================================================

// #include "../Examples/TCD_3D/test0.h"
// #include "../Examples/TCD_3D/test2.h"
// #include "../Examples/TCD_3D/test3.h"
// #include "../Examples/TCD_3D/Linear.h"
// #include "../Examples/TCD_3D/non_linear_in_t.h"
// #include "../Examples/TCD_3D/bubble.h"
// #include "../Examples/TCD_3D/test2mod.h"
// #include "../Examples/TCD_3D/brenn_cd.h"
// #include "../Examples/TCD_3D/brennStutzenCD.h"
// #include "../Examples/TCD_3D/brennStutzenCD_P.h"
//#include "../Examples/TCD_3D/brennStutzenCD_Real.h"
//#include "../Examples/TCD_3D/bulk_compare.h"
//#include "../Examples/TCD_3D/Bail3D.h"
//#include "../Examples/TCD_3D/TEST_RFB.h"
 #include "../Examples/TCD_3D/Sin4.h"
// #include "../Examples/TCD_3D/ConstT.h"
// #include "../Examples/TCD_3D/TeraHertzCube.h"
// #include "../Examples/TCD_3D/TeraHertzBrain.h"
// #include "../Examples/TCD_3D/TeraHertzBreast.h"
 //#include "../Examples/TCD_3D/TeraHertzCube_Smooth.h"
// ======================================================================
// utilities for main program
// ======================================================================

// update the LHS in the Robin BC to the stiffness matrix     
#ifdef __ROBINBC__    
void RobinBCInt(TSquareMatrix3D *A, BoundCondFunct3D *BDCond)
{
  int i, j, k, l, m, l1, N_Facess, IFace, N_RobinFaces, comp;
  int *BeginIndex, *GlobalNumbers, *DOF, *RowPtr, *KCol;
  int N_Active, N_Cells, N_BaseFunct, *N_BaseFuncts;
  int FaceNumbers[MAXN_JOINTS];
  int N_Points, AnsatzDOF, TestDOF, index1, index2;  
  const int *TmpFV, *TmpLen;
  int MaxLen;
    
  double t0, X, Y, Z, xf, yf, zf, Mult;
  double *Weights, *xi, *zeta, n1, n2, n3, detF;   
  double **uref;
  double uorig[MaxN_BaseFunctions3D], *MatValues, area;
//   double a1, a2, a3, b1, b2, b3, n11, n22, n33, detF1;
  
  
  BaseFunct3D *BaseFuncts; 
  TFESpace3D *fespace;
  TCollection *Coll;  
  TJoint *joint;
  TBoundFace *boundface;
  TBaseCell *cell;    
  BoundCond Cond0;
  FE3D FEId;
  TFE3D *ele;
  RefTrans3D RefTrans, *RefTransArray;
  TRefTrans3D *F_K;
  BF3DRefElements RefElement;  
  QuadFormula2D FaceQuadFormula;
  TQuadFormula2D *qf2;  
  QuadFormula3D QF3;
  
  BaseFuncts = TFEDatabase3D::GetBaseFunct3D_IDFromFE3D();
  N_BaseFuncts = TFEDatabase3D::GetN_BaseFunctFromFE3D();  
  RefTransArray = TFEDatabase3D::GetRefTrans3D_IDFromFE3D();  
  
  fespace = A->GetFESpace();
  Coll = fespace->GetCollection();
  N_Cells = Coll->GetN_Cells();
  N_Active = fespace->GetActiveBound();
  BeginIndex = fespace->GetBeginIndex();
  GlobalNumbers = fespace->GetGlobalNumbers();  

  
  RowPtr = A->GetRowPtr();
  KCol = A->GetKCol();  
  MatValues = A->GetEntries();  
  
// heat stress number (or) BIOT NUMBER
  double C0 = TDatabase::ParamDB->P6;
  
  
  for(i=0;i<N_Cells;i++)
   {
    cell = Coll->GetCell(i);
    N_Facess = cell->GetN_Faces();
    IFace = 0;
    for(j=0;j<N_Facess;j++)
     {
      joint = cell->GetJoint(j);
      
      if(joint->GetType() == BoundaryFace)
        {          
          // compute point on the boundary face
          cell->GetShapeDesc()->GetFaceVertex(TmpFV, TmpLen, MaxLen);

          t0 = 1.0/TmpLen[j];
          xf = 0; yf = 0; zf = 0;
          for (l1=0;l1<TmpLen[j];l1++)
           {
            cell->GetVertex(TmpFV[j*MaxLen+l1])->GetCoords(X, Y, Z);
            xf += t0*X;
            yf += t0*Y;
            zf += t0*Z;
           }  
           
          // get boundary condition assigned to the face
          BDCond(xf, yf, zf, Cond0);

          if(Cond0 == ROBIN)
           { 
           // cout << "comp " << comp  <<endl;
            FaceNumbers[IFace] = j;
            IFace++;
          }      
        } //  if(joint->GetTyp
     } //  for(j=0;j<

 
    N_RobinFaces = IFace;    
    if(N_RobinFaces)
     {
       FEId = fespace->GetFE3D(i, cell);        
       N_BaseFunct = N_BaseFuncts[FEId];
       ele = TFEDatabase3D::GetFE3D(FEId);
       RefElement = TFEDatabase3D::GetRefElementFromFE3D(FEId);   
       
       DOF = GlobalNumbers + BeginIndex[i];
      
       l = TFEDatabase3D::GetPolynomialDegreeFromFE3D(FEId);     

       switch(RefElement)
        {
         case BFUnitTetrahedron:          // triangular face 
           FaceQuadFormula = TFEDatabase3D::GetQFTriaFromDegree(2*l);     
           QF3 = TFEDatabase3D::GetQFTetraFromDegree(2*l);
         break;

         case BFUnitHexahedron:          // quadrilateral face
           FaceQuadFormula = TFEDatabase3D::GetQFQuadFromDegree(2*l);
           QF3 = TFEDatabase3D::GetQFHexaFromDegree(2*l);
         break;
       }     
 
       RefTrans = RefTransArray[FEId];
       F_K = TFEDatabase3D::GetRefTrans3D(RefTrans);     
       
      switch(RefTrans)
       {
        case TetraAffin:
          // cout << "TetraAffin" << endl;
          ((TTetraAffin *)F_K)->SetCell(cell);
        break;
        case TetraIsoparametric:
          // cout << "TetraIsoparametric" << endl;
          ((TTetraIsoparametric *)F_K)->SetApproximationOrder(l);
          ((TTetraIsoparametric *)F_K)->SetQuadFormula(QF3);
          ((TTetraIsoparametric *)F_K)->SetCell(cell);
        break;
       case HexaAffin:
         // cout << "HexaAffin" << endl;
         ((THexaAffin *)F_K)->SetCell(cell);
       break;
       case HexaTrilinear:
         // cout << "HexaTrilinear" << endl;
         ((THexaTrilinear *)F_K)->SetCell(cell);
       break;
       case HexaIsoparametric:
         // cout << "HexaIsoparametric" << endl;
         ((THexaIsoparametric *)F_K)->SetApproximationOrder(l);
         ((THexaIsoparametric *)F_K)->SetQuadFormula(QF3);
         ((THexaIsoparametric *)F_K)->SetCell(cell);
       break;
      } // endswitch       
       

       qf2 = TFEDatabase3D::GetQuadFormula2D(FaceQuadFormula);
       qf2->GetFormulaData(N_Points, Weights, xi, zeta);
       TFEDatabase3D::GetBaseFunct3DFromFE3D(FEId)->MakeRefElementData(FaceQuadFormula);     
       
       
       for(j=0;j<N_RobinFaces;j++)
        {     
         IFace = FaceNumbers[j]; 

         uref = TFEDatabase3D::GetJointValues3D(BaseFuncts[FEId], FaceQuadFormula, IFace);       
         area = 0.;
         for(k=0;k<N_Points;k++)
          {
           switch(RefTrans)
            {
             case TetraAffin:
                ((TTetraAffin *)F_K)->GetOuterNormal(IFace, xi[k], zeta[k], n1, n2, n3); 
             break;
             case TetraIsoparametric:
                 ((TTetraIsoparametric *)F_K)->GetOuterNormal(IFace, xi[k], zeta[k], n1, n2, n3); 
             break;
             case HexaAffin:
                 ((THexaAffin *)F_K)->GetOuterNormal(IFace, xi[k], zeta[k], n1, n2, n3);   
             break;
             case HexaTrilinear:
                 ((THexaTrilinear *)F_K)->GetOuterNormal(IFace, xi[k], zeta[k], n1, n2, n3);  
             break;
             case HexaIsoparametric:
                ((THexaIsoparametric *)F_K)->GetOuterNormal(IFace, xi[k], zeta[k], n1, n2, n3);
             break;
             default:
               Error("Wrong RefTrans" << endl);
#ifdef _MPI      
              MPI_Finalize(); 
#endif    
              exit(0);       
            } // endswitch       
       
            // area of the face
            detF = sqrt(n1*n1+n2*n2+n3*n3);

           // D000
           for(l=0;l<N_BaseFunct;l++)
            uorig[l] = uref[k][l];

           Mult = C0*detF*Weights[k];     
           //area += detF*Weights[k];
           
           for(l=0;l<N_BaseFunct;l++)
            {
             TestDOF = DOF[l];
             index2 = RowPtr[TestDOF+1];       
             for(m=0;m<N_BaseFunct;m++)
              {
               AnsatzDOF = DOF[m];
               index1 = RowPtr[TestDOF];
                while(KCol[index1] != AnsatzDOF) index1++;
                
                MatValues[index1] += (Mult*uorig[m]*uorig[l]);                
                //cout << "A: (" << TestDOF << " " << AnsatzDOF << ") =  " << Mult*uorig[m]*uorig[l] << endl;
              } // for(m=0;m<N_BaseF
            } // for(l=0;l<N_B
         } // for(k=0;k<N_Points;k
         // printf(" %d  N_RobinFaces  %d  RobinBCInt  area %f   \n", i, N_RobinFaces, area);  
        }//  for(j=0;j<N_RobinFaces
       
 
      
     } // N_RobinFaces
     
   } // for(i=0;i
  
// #ifdef _MPI      
// //    if(rank==0)
//     printf("Main Programm  RobinBCInt \n");  
//     MPI_Finalize(); 
// #endif    
//   exit(0);  
        
} // void RobinBCInt(TSquareMatrix3D *A
#endif

int main(int argc, char* argv[])
{
  TDatabase *Database = new TDatabase();    
  double t1, t2, t3, time1,t_par1, t_par2, start_time, total_time, hmin, hmin_all, hmax;  
  
#ifdef _MPI
  const int root = 0;
  int rank, size;
  int MaxCpV, MaxSubDomainPerDof;
  int  out_rank, N_Cells_loc;

  TParVector3D  *ParSolVect, *ParRhsVect;
  TParFECommunicator3D *ParComm;

  MPI_Request request001, request002, request003, request004, request005;
  MPI_Status status;

  MPI_Comm Comm = MPI_COMM_WORLD;
  MPI_Init(&argc, &argv);

  time1 =  MPI_Wtime();
  MPI_Allreduce(&time1, &start_time, 1, MPI_DOUBLE, MPI_MIN, Comm);

  MPI_Comm_rank(Comm, &rank);
  MPI_Comm_size(Comm, &size);
 
  TDatabase::ParamDB->Comm = Comm;
  TScalar_ParSolver *Par_Solver; 
#endif  
  
  TDomain *Domain = new TDomain(); 
  TFEDatabase3D *FEDatabase = new TFEDatabase3D();
  TCollection *coll, *mortarcoll = NULL;
  TBaseCell *cell;
  TFESpace3D *concentration_space, *pressure_space, *velocity_space; 
  TOutput3D *Output; 
  TAuxParam3D *aux;
  TFEFunction3D *conc;
  TSquareStructure3D *sqstructureA;  
  TSquareMatrix3D *sqmatrixA, *sqmatrixM, *SQMATRICES[3];;
  TSquareMatrix **sqmatrices = (TSquareMatrix **)SQMATRICES;
  TFESpace3D *fesp[2], *ferhs[1];

  MatVecProc *MatVect;
  DefectProc *Defect;  
  TDiscreteForm3D *DiscreteForm;  
  TDiscreteForm3D *DiscreteFormMatrixMRhs;
  TDiscreteForm3D *DiscreteFormMatrixARhs;
  TDiscreteForm3D *DiscreteFormRhs;  
  BoundCondFunct3D *BoundaryConditions[1];
  BoundValueFunct3D *BoundValues[1];
  CoeffFct3D *Coefficients[1];
  
  double *sol, *oldsol, *rhs, *B;
  double *defect, *startsol, *frac_step_sol, *oldrhs;
  double *RHSs[3], gamma, tau, L2error_Max=0, L2error_Max_t;
  double delta, end_time, solver_time, oldtau, solver_time_curr, errors[7];  
  double l2, olderror, H1, olderror1, start_steptime, end_steptime, hK, hK_all;
  
  int i, j, k, l, m, ret, ORDER, N_Unknowns, N_Active, N_LinIter;
  int N_SquareMatrices, N_FESpaces, N_Rhs, img=0, N_SubSteps;
  int time_discs, methods, very_first_time=0, N_LinIterCurr;
  int first_matrix_assemble =1;

  bool UpdateConvection=FALSE, UpdateRhs=FALSE, ConvectionFirstTime=TRUE, FACTORIZE=TRUE;

 
  char *PsBaseName, *VtkBaseName;
  char *PRM, *GEO;
    
  char ReadinDat[] = "readin.dat";
  char UString[] = "T";
  char PString[] = "p";
  char NameString[] = "name";
  char MMString[] = "Mass matrix";
  char SubID[] = "";
  
  std::ostringstream os;
  os << " ";

  total_time = GetTime();
//======================================================================
// read parameter file
//======================================================================
  if(argc>=2)
    ret=Domain->ReadParam(argv[1]);
  else
    ret=Domain->ReadParam(ReadinDat);

  if(ret==-1)
    exit(-1);
 
  OpenFiles();
  OutFile.setf(std::ios::scientific);

#ifdef _MPI
  out_rank=TDatabase::ParamDB->Par_P0;

  if(rank==out_rank)
#endif
   {    
    Database->WriteParamDB(argv[0]);
    Database->WriteTimeDB();
   }
   
  ExampleFile();
  
  PRM = TDatabase::ParamDB->BNDFILE;
  GEO = TDatabase::ParamDB->GEOFILE;
  
  PsBaseName = TDatabase::ParamDB->PSBASENAME;  
  VtkBaseName = TDatabase::ParamDB->VTKBASENAME;
  
  MatVect = MatVect_Scalar;
  Defect = Defect_Scalar;
  
#ifdef _MPI
   TDatabase::ParamDB->SOLVER_TYPE = 101;
#else
   TDatabase::ParamDB->SOLVER_TYPE = 2;
#endif
 
//======================================================================
// initialize discrete forms
//======================================================================

// discrete form for assembling mass matrix
   DiscreteFormMatrixMRhs = new TDiscreteForm3D
    (MMString, MMString, N_Terms_MatrixMRhs,
     Derivatives_MatrixMRhs,
     SpacesNumbers_MatrixMRhs, N_Matrices_MatrixMRhs, N_Rhs_MatrixMRhs,
     RowSpace_MatrixMRhs, ColumnSpace_MatrixMRhs, RhsSpace_MatrixMRhs,
     MatrixMRhsAssemble, BilinearCoeffs, NULL);  
  
   DiscreteFormMatrixARhs = new TDiscreteForm3D
    (MMString, MMString, N_Terms_MatrixARhs,
     Derivatives_MatrixARhs,
     SpacesNumbers_MatrixARhs, N_Matrices_MatrixARhs,
     N_Rhs_MatrixARhs,
     RowSpace_MatrixARhs, ColumnSpace_MatrixARhs,
     RhsSpace_MatrixARhs,
     MatrixARhsAssemble, BilinearCoeffs, NULL);
  
   DiscreteFormRhs = new TDiscreteForm3D
    (MMString, MMString, N_Terms_Rhs, Derivatives_Rhs,
     SpacesNumbers_Rhs, N_Matrices_Rhs, N_Rhs_Rhs,
     RowSpace_Rhs, ColumnSpace_Rhs, RhsSpace_Rhs,
     RhsAssemble, BilinearCoeffs, NULL);
  
//======================================================================
// read boundary parameterization and initialize coarse grid
//======================================================================
  Domain->Init(PRM, GEO);
  
    /** with atlas mesh, no tetgen*/
//     TetrameshCreate(Domain);
    
   /** Using tetgen with smesh mesh */
//    TetrameshGen(Domain);
   
// //    /** same as TetrameshGen but without using the face list info from tetgen */
// // // //     TetraGen(Domain);

// // //   TetraGenWithInputCells(Domain);
//   Domain->ConvertQuadToTri(TDatabase::ParamDB->CONVERT_QUAD_TO_TRI);

  for(i=0;i<TDatabase::ParamDB->UNIFORM_STEPS;i++)
    Domain->RegRefineAll();  

  //======================================================================
  // Partition grid using Metis
  //======================================================================
#ifdef _MPI
  Domain->GenerateEdgeInfo();

  t_par1 = MPI_Wtime();
  Partition_Mesh3D(Comm, Domain, MaxCpV);
  t_par2 = MPI_Wtime();

  if(rank==0)
    printf("Time taken for Domain Decomposition is %e\n", (t_par2-t_par1));

  MaxSubDomainPerDof = MIN(MaxCpV, size);
#endif

/*#ifdef _MPI      
   if(rank==0)
    printf("Main Programm    %d \n",  rank );  
    MPI_Finalize(); 
#endif    
  exit(0);*/   

  
  TDatabase::TimeDB->CURRENTTIME = TDatabase::TimeDB->STARTTIME;
  SetTimeDiscParameters();

  t3 = GetTime();
  total_time = t3 - total_time; 
  
  
//======================================================================
// construct all fespace
//====================================================================== 
// create collection of mesh cells
#ifdef _MPI
   coll=Domain->GetOwnCollection(It_Finest, 0, rank);
#else
   coll=Domain->GetCollection(It_Finest, 0);
#endif    

   ORDER  = TDatabase::ParamDB->ANSATZ_ORDER; 
 
   concentration_space = new TFESpace3D(coll, NameString, UString, BoundCondition, ORDER);
   N_Unknowns = concentration_space->GetN_DegreesOfFreedom();
   
   // active dof (i.e. dof without Dirichlet dofs)
   N_Active = concentration_space->GetActiveBound();
   OutPut("degrees of freedom: "<< setw(10) << N_Unknowns << endl);

   
#ifdef _MPI
    t_par1 = MPI_Wtime();

    concentration_space->SetMaxSubDomainPerDof(MaxSubDomainPerDof);
    ParComm = new TParFECommunicator3D(Comm, concentration_space);
  
    t_par2 = MPI_Wtime();
    if(rank==out_rank)
     {
      printf("Time taken for FeSpace SubDomain dof mapping %e\n", (t_par2-t_par1));
      printf("DOF of FeSpace  space : %d \n", ParComm->GetN_GlobalDegreesOfFreedom());      
     }
#else
    OutPut(" Rank " <<  " DOF Scalar : " << setw(10) << N_Unknowns << endl);
#endif
    
    
//     #ifdef _MPI      
//    if(rank==0)
//     printf("Main Programm    %d \n",  rank );  
//     MPI_Finalize(); 
// #endif    
//   exit(0);  
//======================================================================
// allocate memory for all matrices
//======================================================================  
   // first build matrix structure
   sqstructureA = new TSquareStructure3D(concentration_space);
   sqstructureA->Sort();
  
   // two matrices used
   // A contains the non time dependent part of the discretization
   sqmatrixA = new TSquareMatrix3D(sqstructureA);
   
   // M is the mass matrix
   // the iterative solver uses M
   sqmatrixM = new TSquareMatrix3D(sqstructureA);
  
//======================================================================
// allocate memory for sol and rhs arrays and construct fefunction
//======================================================================  
   sol = new double[N_Unknowns];
   oldsol = new double[N_Unknowns];   
   rhs = new double[N_Unknowns];
   B = new double [N_Unknowns];

#ifdef _MPI
    ParSolVect =  new TParVector3D(Comm, sol, N_Unknowns, 1, ParComm);
    ParRhsVect =  new TParVector3D(Comm, B, N_Unknowns, 1, ParComm);
#endif

   
   memset(rhs, 0, N_Unknowns*SizeOfDouble);  
   memset(sol, 0, N_Unknowns*SizeOfDouble);
   
   // allocate fe function 
   conc = new TFEFunction3D(concentration_space, UString, UString, sol, N_Unknowns); 
   
   // interpolate initial condition
//    conc->Interpolate(InitialCondition);
//   OutPut("Interpolation time: " << TDatabase::TimeDB->CURRENTTIME << endl);

#ifdef _MPI
     // initialize the parallel solver
     Par_Solver = new TScalar_ParSolver(ParComm, sqstructureA, 1, ParRhsVect, ParSolVect);
#endif   
   
//======================================================================
// prepare output,
//======================================================================     
   Output = new TOutput3D(1, 1, 0, 1, Domain);
   Output->AddFEFunction(conc);
   os.seekp(std::ios::beg);
   Output->AddParameter(TDatabase::TimeDB->CURRENTTIME,os.str().c_str());
   
//======================================================================
// assembling of mass matrix and rhs
//======================================================================
  // set parameters
   N_Rhs = 1;
   N_FESpaces = 1;
   fesp[0] = concentration_space;
   
   // reset matrices
   N_SquareMatrices = 1;
   SQMATRICES[0] = sqmatrixM;
   SQMATRICES[0]->Reset();   
   
   DiscreteForm = DiscreteFormMatrixMRhs; 
   BoundaryConditions[0] =  BoundCondition;
   BoundValues[0] = BoundValue;
   
   memset(rhs, 0, N_Unknowns*SizeOfDouble);  
   RHSs[0] = rhs;
   ferhs[0] = concentration_space;
    
   aux =  new TAuxParam3D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);  
   
   // check the min X in the domain, check BilinearCoeffs in Example file
   TDatabase::ParamDB->P14=-1.;
   
   Assemble3D(N_FESpaces, fesp, 
              N_SquareMatrices, SQMATRICES, 
              0, NULL, 
              N_Rhs, RHSs, ferhs,
              DiscreteForm, 
              BoundaryConditions, 
              BoundValues, 
              aux);
      
//    OutPut("Min X = " << TDatabase::ParamDB->P13<<endl);
//    exit(0);
   
   delete aux;    

   
    hK_all = TDatabase::ParamDB->P14;
#ifdef _MPI         
   MPI_Allreduce(&hK_all, &hK, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#else
   hK = hK_all;
#endif 
   
   
#ifdef _MPI     
   if(rank==out_rank)      
#endif   
    OutPut("Cell Diameter (within Sourceradius): " <<  hK <<  endl);      
   
  // save solution
  memcpy(oldsol,sol,N_Unknowns*SizeOfDouble);
  
  // allocate arrays for solver
  defect = new double[N_Unknowns];
  startsol = new double[N_Unknowns];
  frac_step_sol = new double[N_Unknowns];
  oldrhs =  new double[N_Unknowns];
  solver_time = 0.0;
  N_LinIter = 0;
  
  // parameters for time stepping scheme
  gamma = 0;
  m = 0;
  N_SubSteps = GetN_SubSteps();
  oldtau = 1;
  end_time = TDatabase::TimeDB->ENDTIME;

  // not active : TDatabase::TimeDB->TIMESTEPLENGTH_CONTROL = 0 
  if (TDatabase::TimeDB->TIMESTEPLENGTH_CONTROL)
    time_discs = 2;
  else
    time_discs = 1;  
  
  
  coll->GetHminHmax(&hmin,&hmax);
  OutPut("h_min : " << hmin << " h_max : " << hmax << endl);
  
#ifdef _MPI
   MPI_Allreduce(&hmin, &hmin_all, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
#else
   hmin_all = hmin;
#endif
  
   // measure errors to known solution
   if(TDatabase::ParamDB->MEASURE_ERRORS)
    {
      
     for(i=0;i<7;i++)
       errors[i] = 0.;           
      
        aux =  new TAuxParam3D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);
        conc->GetErrors(Exact, 4, AllDerivatives, 2, L2H1Errors, BilinearCoeffs, aux, 1, fesp, errors);
        delete aux;

#ifdef _MPI
        MPI_Reduce(errors, &l2, 1, MPI_DOUBLE, MPI_SUM, out_rank, MPI_COMM_WORLD);
        MPI_Reduce(errors+1, &H1, 1, MPI_DOUBLE, MPI_SUM, out_rank, MPI_COMM_WORLD);

        if(rank==out_rank)
         {
          l2 = sqrt(l2);
          H1 = sqrt(H1);
         }
#else
        l2 = errors[0];
        H1 = errors[1];
#endif


#ifdef _MPI
        if(rank==out_rank)
#endif
        {
         OutPut("time: " << TDatabase::TimeDB->CURRENTTIME);
         OutPut(" L2: " << l2);
         OutPut(" H1-semi: " << H1 << endl);
        }

//         if(L2error_Max<l2)
//          {
//           L2error_Max=l2;
//           L2error_Max_t=TDatabase::TimeDB->CURRENTTIME;
//          }
//  
//         olderror = l2;
//         olderror1 = H1;   
      } // endif MEASURE_ERRORS  
  
 
//======================================================================
// output
//======================================================================     
  if(TDatabase::ParamDB->WRITE_VTK)
   {
#ifdef _MPI     
       Output->Write_ParVTK(Comm, img, SubID);
#else     
    os.seekp(std::ios::beg);      
    if(img<10) os << VtkBaseName<<".0000"<<img<<".vtk" << ends;
    else if(img<100) os << VtkBaseName<<".000"<<img<<".vtk" << ends;
    else if(img<1000) os << VtkBaseName<<".00"<<img<<".vtk" << ends;
    else if(img<10000) os << VtkBaseName<<".0"<<img<<".vtk" << ends;
    else  os << VtkBaseName<<"."<<img<<".vtk" << ends;
    Output->WriteVtk(os.str().c_str());  
    
    
    os.seekp(std::ios::beg);      
    if(img<10) os << VtkBaseName<<".0000"<<img<<".plt" << ends;
    else if(img<100) os << VtkBaseName<<".000"<<img<<".plt" << ends;
    else if(img<1000) os << VtkBaseName<<".00"<<img<<".plt" << ends;
    else if(img<10000) os << VtkBaseName<<".0"<<img<<".plt" << ends;
    else  os << VtkBaseName<<"."<<img<<".plt" << ends;
    Output->WriteBinaryPlt(os.str().c_str());      
    
#endif
    
    img++;
   }     
 
// #ifdef _MPI      
//    if(rank==0)
//     printf("Main Programm    %d \n",  rank+1 );  
//     MPI_Finalize(); 
// #endif    
//   exit(0);  
        
 
//  TDatabase::TimeDB->TIMESTEPLENGTH = pow(hmin_all, TDatabase::ParamDB->P9);
 
 
 OutPut("TIMESTEPLENGTH : " << TDatabase::TimeDB->TIMESTEPLENGTH << endl;)
  
//======================================================================
// start of time cycle
// everything happens on the same grid
//======================================================================
  while(TDatabase::TimeDB->CURRENTTIME< end_time)
  {                              // time cycle
#ifdef _MPI    
   start_steptime = MPI_Wtime();
#else  
   start_steptime = GetTime();
#endif
      
    m++;
    TDatabase::TimeDB->INTERNAL_STARTTIME = TDatabase::TimeDB->CURRENTTIME;
    
    for (methods=0;methods<time_discs;methods++)
     {
      if (time_discs==2)
      {
        if (methods==0) // fractional-step-theta-scheme
        {
          TDatabase::TimeDB->TIME_DISC = 3;
          memcpy(startsol,sol,N_Unknowns*SizeOfDouble); // save start sol
          memcpy(oldrhs,rhs,N_Unknowns*SizeOfDouble); // save rhs
        }
        else           // crank nicolson scheme
        {              // take solution of first scheme as initial iterate
          TDatabase::TimeDB->TIME_DISC = 2;
          TDatabase::TimeDB->CURRENTTIME = TDatabase::TimeDB->INTERNAL_STARTTIME;
          memcpy(frac_step_sol,sol,N_Unknowns*SizeOfDouble); // save solution of fract.step
          memcpy(sol,startsol,N_Unknowns*SizeOfDouble); // get former startsol
          memcpy(rhs,oldrhs,N_Unknowns*SizeOfDouble); // get old rhs
        }
        N_SubSteps = GetN_SubSteps();
      }
     
      for(l=0;l<N_SubSteps;l++)      // sub steps of fractional step theta
      {
        if (!very_first_time)
        {
          SetTimeDiscParameters();
        }
        if (m==1
#ifdef _MPI
        &&   rank==out_rank
#endif
           )
         {
          OutPut("Theta1: " << TDatabase::TimeDB->THETA1<< endl);
          OutPut("Theta2: " << TDatabase::TimeDB->THETA2<< endl);
          OutPut("Theta3: " << TDatabase::TimeDB->THETA3<< endl);
          OutPut("Theta4: " << TDatabase::TimeDB->THETA4<< endl);
         }
        
        tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
        if (!very_first_time)
          TDatabase::TimeDB->CURRENTTIME += tau;
         j=0; 
         
#ifdef _MPI
      if(rank==out_rank)
#endif
      {         
        OutPut(endl << "CURRENT TIME: ");
        OutPut(TDatabase::TimeDB->CURRENTTIME << endl);     
      }
      
        solver_time_curr = 0;
   
        N_LinIterCurr = 0; 
      
              // working array for rhs is B, initialize B
        memset(B, 0, N_Unknowns*SizeOfDouble);
        // old rhs multiplied with current subtime step and theta3 on B
        Daxpy(N_Active, tau*TDatabase::TimeDB->THETA3, rhs, B);
 
 
     
//======================================================================
// assembling of mass matrix and rhs
//======================================================================
// assembling of A and rhs
//time dependent rhs need to be assembled at every time step
       if(UpdateConvection || UpdateRhs ||  ConvectionFirstTime )
        {
          // set parameters
#ifdef _MPI    
         t1 = MPI_Wtime();
#else  
         t1 = GetTime();
#endif
          N_FESpaces = 1;
          fesp[0] = concentration_space;
   
           // reset matrices
           N_SquareMatrices = 1;
           SQMATRICES[0] = sqmatrixA;
           SQMATRICES[0]->Reset();   
   
           DiscreteForm = DiscreteFormMatrixARhs; 
           BoundaryConditions[0] =  BoundCondition;
           BoundValues[0] = BoundValue;
   
           N_Rhs = 1;
           memset(rhs, 0, N_Unknowns*SizeOfDouble);  
           RHSs[0] = rhs;
           ferhs[0] = concentration_space;
    
           aux =  new TAuxParam3D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);  
   
           Assemble3D(N_FESpaces, fesp, 
                      N_SquareMatrices, SQMATRICES, 
                      0, NULL, 
                      N_Rhs, RHSs, ferhs,
                      DiscreteForm, 
                      BoundaryConditions, 
                      BoundValues, 
                      aux);  
#ifdef __ROBINBC__                          
         RobinBCInt(SQMATRICES[0], BoundaryConditions[0]);                     
#endif                      

#ifdef _MPI
          if(UpdateConvection ||  ConvectionFirstTime)
            FACTORIZE=TRUE;
#endif
           ConvectionFirstTime = FALSE;
           delete aux;   
           
        if(TDatabase::ParamDB->SC_VERBOSE>1)
         {
#ifdef _MPI    
          t2 = MPI_Wtime();
#else  
          t2 = GetTime();
#endif   
#ifdef _MPI      
          if(rank==out_rank)      
#endif               
           OutPut("Time for assembling " << t2-t1 << "s "<<  endl);   
          }
                 
         }//  if(UpdateConvection || UpdateRhs ||  ConvectionFirstTime )
          
          
      // add rhs from current sub time step to rhs array B
        Daxpy(N_Active, tau*TDatabase::TimeDB->THETA4, rhs, B);       
        //======================================================================
        // manipulation of matrices due to current time discretization
        // the stiffness matrix is stored on M 
        //======================================================================
        oldtau = tau;

        /** scaled by current sub time step length and theta2 */
        /** currently : M := M + gamma A */
        MatAdd(sqmatrixM, sqmatrixA,  - tau*TDatabase::TimeDB->THETA2);

        /** set current factor of steady state matrix */
        gamma = -tau*TDatabase::TimeDB->THETA2;

        /** defect = M * sol */
        /** B:= B + defect */
        MatVectActive(sqmatrixM, sol, defect);
        Daxpy(N_Active, 1, defect, B);
 
        /** set Dirichlet values */
        /** RHSs[0] still available from assembling */
        memcpy(B+N_Active, RHSs[0]+N_Active, (N_Unknowns-N_Active)*SizeOfDouble);
        /** copy Dirichlet values from rhs into sol */
        memcpy(sol+N_Active, RHSs[0]+N_Active, (N_Unknowns-N_Active)*SizeOfDouble);
       
        /** M = M + (-gamma + tau*TDatabase::TimeDB->THETA1) A */
        MatAdd(sqmatrixM, sqmatrixA, -gamma + tau*TDatabase::TimeDB->THETA1);
        /** set current factor of steady state matrix */
        gamma = tau*TDatabase::TimeDB->THETA1;   

        //======================================================================
        // solution of linear system
        //======================================================================
#ifndef _MPI
        t1 = GetTime();
        DirectSolver(sqmatrixM, B, sol);
        t2 = GetTime();
        solver_time_curr = t2-t1;
        solver_time += t2-t1;          
#else
        t1 = MPI_Wtime();
        Par_Solver->Solve(sqmatrixM, FACTORIZE); 
        FACTORIZE=FALSE;    
#endif    
     
        if(TDatabase::ParamDB->SC_VERBOSE>1)
         {
#ifdef _MPI    
          t2 = MPI_Wtime();
#else  
          t2 = GetTime();
#endif   
#ifdef _MPI      
          if(rank==out_rank)      
#endif               
           OutPut("Time for Sover " << t2-t1 << "s "<<  endl);   
          }
     
        MatAdd(sqmatrixM, sqmatrixA, -gamma);
        // set current factor of steady state matrix
        gamma = 0.; 
     
       } //  for(l=0;l<N_SubSteps;l++)     
     
     if(TDatabase::ParamDB->WRITE_VTK &&
        (m % TDatabase::TimeDB->STEPS_PER_IMAGE == 0) )
      {
#ifdef _MPI     
       Output->Write_ParVTK(Comm, img, SubID);
#else     
       os.seekp(std::ios::beg);      
       if(img<10) os << VtkBaseName<<".0000"<<img<<".vtk" << ends;
       else if(img<100) os << VtkBaseName<<".000"<<img<<".vtk" << ends;
       else if(img<1000) os << VtkBaseName<<".00"<<img<<".vtk" << ends;
       else if(img<10000) os << VtkBaseName<<".0"<<img<<".vtk" << ends;
       else  os << VtkBaseName<<"."<<img<<".vtk" << ends;
       Output->WriteVtk(os.str().c_str());
#endif       
       img++;
      }     
   
    } //for (methods=0;methods<t

    
      // measure errors to known solution
      if(TDatabase::ParamDB->MEASURE_ERRORS)
      {
        aux =  new TAuxParam3D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);
        conc->GetErrors(Exact, 4, AllDerivatives, 2, L2H1Errors, BilinearCoeffs, aux, 1, fesp, errors);
        delete aux;
        
#ifdef _MPI
        MPI_Reduce(errors, &l2, 1, MPI_DOUBLE, MPI_SUM, out_rank, MPI_COMM_WORLD);
        MPI_Reduce(errors+1, &H1, 1, MPI_DOUBLE, MPI_SUM, out_rank, MPI_COMM_WORLD);

        if(rank==out_rank)
         {
          l2 = sqrt(l2);
          H1 = sqrt(H1);
         }
#else
        l2 = errors[0];
        H1 = errors[1];
#endif  

#ifdef _MPI
       if(rank==out_rank)
#endif
       {
        OutPut("time: " << TDatabase::TimeDB->CURRENTTIME);
        OutPut(" L2: " << l2);
        OutPut(" H1-semi: " << H1 << endl);

        if(L2error_Max<l2)
         {
          L2error_Max=l2;
          L2error_Max_t=TDatabase::TimeDB->CURRENTTIME;
         }
 
        olderror = l2;
        olderror1 = H1;
        if(m>1)
         {
          errors[3] += (l2*l2 +olderror * olderror)*TDatabase::TimeDB->TIMESTEPLENGTH/2.0;
          errors[4] += (H1*H1 +olderror1 * olderror1)*TDatabase::TimeDB->TIMESTEPLENGTH/2.0;
          
          OutPut(L2error_Max_t <<  " L2error_Max " << L2error_Max << " ");
          OutPut(TDatabase::TimeDB->CURRENTTIME <<  " L2(0,T;L2) " << sqrt(errors[3]) << " ");
          OutPut( "L2(0,T;H1) " << sqrt(errors[4]) << endl);
          }
        }
      } // endif MEASURE_ERRORS

#ifdef _MPI      
   if(rank==0)
    printf("Main Programm    %d \n",  rank+1 );  
    MPI_Finalize(); 
#endif    
 
        

#ifdef _MPI    
   end_steptime = MPI_Wtime();
#else  
   end_steptime = GetTime();
#endif
   
   
#ifdef _MPI       
   if(rank==out_rank)      
#endif     
    OutPut("Total time for this time step " << end_steptime - start_steptime << endl);         
  } // while

  
  
   
#ifdef _MPI    
   MPI_Finalize();  
#endif       
  CloseFiles();
  return 0;
}
  
  

