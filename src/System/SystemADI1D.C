/** =======================================================================
* @class     TSystemADI1D
* @brief     stores the information of system ADI1D
* @author    Sashikumaar Ganesan
* @date      13.04.2020
* @History 
* ======================================================================= */
#include <SystemADI.h>
#include <SystemADI1D.h>
#include <Database.h>
#include <FEDatabase2D.h>
#include <DirectSolver.h>
#include <SquareStructure1D.h>
#include <SquareMatrix1D.h>
#include <FEFunction1D.h>
#include <LineAffin.h>
#include <LinAlg.h>

#include <MacroCell.h>
#include <JointEqN.h>
#include <NodalFunctional1D.h>
#include <SquareStructure1D.h>
#include <SquareMatrix1D.h>
#include <SquareMatrix.h>
#include <Matrix.h>

#include <string.h>
#include <sstream>
#include <MooNMD_Io.h>
#include <stdio.h>
#include <stdlib.h>

TSystemADI1D::TSystemADI1D():TSystemADI()
{
  
}

TSystemADI1D::TSystemADI1D(int N_L, double start, double end, BoundCond1D *boundConLminLMax, 
                           DoubleFunctND *growthAndNuc): TSystemADI(growthAndNuc)
{
 //mesh
 Domain_Intl = new TDomain();
 this->Generate1DMesh(start, end, N_L);
 Coll_Intl = Domain_Intl->GetCollection(It_Finest, 0);
 cout<< "N_Cells_Internal " << Coll_Intl->GetN_Cells() <<endl;

 // Finite difference in internal 
 if(TDatabase::ParamDB->INTL_DISCTYPE==FD)
  { 
   TDatabase::ParamDB->ANSATZ_ORDER_INTL = 1;
  }

 char IString[] = "I";
 FESpace1D_Intl = new TFESpace1D(Coll_Intl, IString, IString, TDatabase::ParamDB->ANSATZ_ORDER_INTL);
 N_Dof = FESpace1D_Intl->GetN_DegreesOfFreedom();

 if(TDatabase::ParamDB->ANSATZ_ORDER_INTL<0)
  {
   FESpace1D_Intl->SetAsDGSpace(); 
   dGDisc = TRUE;
   TDatabase::ParamDB->INTL_DISCTYPE=DG;  // no supg method
   DofSameAsNodalPts = FALSE;
  }
  else
  {
    dGDisc = FALSE;
    DofSameAsNodalPts = TRUE;  // not necessary to transfer solutions if true
  }
  

  TSquareStructure1D *SqStructureIntl = new TSquareStructure1D(FESpace1D_Intl);
  SqStructureIntl->Sort();

  M_Intl = new TSquareMatrix1D(SqStructureIntl);
  A_Intl= new TSquareMatrix1D(SqStructureIntl);

  if(TDatabase::ParamDB->INTL_DISCTYPE==SUPG)
   {
    S_Intl = new TSquareMatrix1D(SqStructureIntl);
    K_Intl = new TSquareMatrix1D(SqStructureIntl);
   }  

  N_MatEntries = M_Intl->GetN_Entries();
  Mvalues_orig = new double[N_MatEntries];
  Advectvalues_orig = new double[N_MatEntries];

  //read boundary conditions 
  boundConLminLMax(cond_Lmin, cond_Lmax);

  // nodal points to evaluate nodal values
  this->GenerateNodalPts();
  OutPut("Internal space Dof: "<<  N_Dof << " Nodal Points: " <<  N_NodalPts<< endl);
}

void TSystemADI1D::Init(int n_Coord, int n_Xpos, double *xpos, int n_ADISystems, int *n_LnodalPos, double **LnodalPos, 
                        int ownADI_Idx, double *sol_XposLNnodal)
{
  // Nodal_sol = new double[N_NodalPts];
  // sol = new double[N_Dof];
  // oldsol = new double[N_Dof];
  rhs = new double[N_Dof];
  B = new double[N_Dof];
  // defect = new double[N_Dof];

 //  memset(Nodal_sol, 0, N_NodalPts*SizeOfDouble);
 //  memset(sol, 0, N_Dof*SizeOfDouble);
 //  memset(oldsol, 0, N_Dof*SizeOfDouble);
 memset(rhs, 0, N_Dof*SizeOfDouble);
 memset(B, 0, N_Dof*SizeOfDouble);
 //  memset(defect, 0, N_Dof*SizeOfDouble);
 char IString[] = "I";

 FENodalFunction = new TFEFunction1D(FESpace1D_Intl, IString, IString, Nodal_sol,  N_NodalPts);

 N_Xpos = n_Xpos;
 Xpos = xpos; 
 N_ADISystems = n_ADISystems;
 N_LnodalPos = n_LnodalPos;
 LnodalPos = LnodalPos;
 OwnADI_Idx = ownADI_Idx;
 Sol_XposLNnodal = sol_XposLNnodal;
 N_Coord = n_Coord;
  
 int N;
 if(N_LnodalPos[OwnADI_Idx]> N_Dof) 
  { N=N_LnodalPos[OwnADI_Idx];}
 else
  { N= N_Dof; }

 IncidentArray = new int[N];
} // Init


void TSystemADI1D::Interpolate(int n_Coord, double *Coords, double *Sol, DoubleFunctND *Exact)
{
 FENodalFunction->InterpolateNodalPts(n_Coord, Coords, Exact, Sol);
}

void TSystemADI1D::Solve(int N_Param, double *Coords, CoeffFctND *Bilinear, double *Sol)
{
 int i, j, k, l;
 double gamma=0., BDValue1, BDValue2;
 double Out[3], G, Bnuc; 
 double tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;;
 double *MatValues = M_Intl->GetEntries();
 double *AdvectMatValues = A_Intl->GetEntries();

 if(cond_Lmin==DIRICHLET && TDatabase::ParamDB->INTL_DISCTYPE!=DG)
  BDValue1 = Sol[0];

 if(cond_Lmax==DIRICHLET && TDatabase::ParamDB->INTL_DISCTYPE!=DG)
  BDValue2 = Sol[N_Dof -1];

 /** set I as non-negative */
//  for(j=0;j<N_Dof;j++)
//  {  if(Sol[j]<0.) Sol[j]=0. ; }

  if(TDatabase::ParamDB->INTL_DISCTYPE==SUPG)
   {
    S_Intl->Reset();
    K_Intl->Reset();
   }

  /**growth rate based on concentration defined in Example file*/
  if(GetGrowthAndNuc)
   {
    GetGrowthAndNuc(N_Param, Coords, Out);
    G = Out[0];
    Bnuc = Out[1];
    BDValue1 = Out[1];
    BDValue2 =  Out[2];
  }
  else
  {
   G = 0.;
   BDValue1 = 0; //zero neumann
   BDValue2 = 0; //zero neumann

     if(cond_Lmin==DIRICHLET)
       {
         ErrMsg(" Bound vvalues must be provided through GetGrowthAndNuc");
         exit(0);
       }

     if(cond_Lmax==DIRICHLET)
       {
         ErrMsg(" Bound vvalues must be provided through GetGrowthAndNuc");
         exit(0);
       }
  }
   
  // cout << "G  : " << G <<" Bnuc " << Bnuc <<endl; 

  A_Intl->Reset();

  /** init rhs */
  memset(rhs, 0, N_Dof*SizeOfDouble);

  /** assemble matrices */
  switch(TDatabase::ParamDB->INTL_DISCTYPE)
  {
       case FD: // Finite difference, Il'in-Allen-Southwell scheme
           this->AssembleARhs_FD(G);
       break;      
       case  GALERKIN:
            this->AssembleARhs(Coords, G, Bilinear);
       break;  
       case SDFEM:     // SUPG
   
            this->AssembleARhs_SUPG(Coords, G, Bilinear); 
       break;  
       case DG:
          //  this->AssembleARhs_DG_Advect(Coords, Out, Bilinear);  

           // restore the advection matrix mat values and scale with advection coeff
           memcpy(AdvectMatValues,  Advectvalues_orig,  N_MatEntries*SizeOfDouble);
           Dscal(N_MatEntries,  G, AdvectMatValues);

          //update B_nuc
          this->AssemblDG_Rhs(Bnuc);  
       break;
       default:
            Error("only FD, GALERKIN, SUPG, DG are implemented" << endl);
            Error("file: " << __FILE__ << " line " << __LINE__ << endl);
            exit(-1);
       break;
  }

   // restore mass mat  values
   memcpy(MatValues,  Mvalues_orig,  N_MatEntries*SizeOfDouble);

   /** store A and other mat values for rhs calculation */     
   if(fabs(TDatabase::TimeDB->THETA2)>0)  
    {
     // M = M -(tau*TDatabase::TimeDB->THETA2) A
     MatAdd(M_Intl, A_Intl,  -tau*TDatabase::TimeDB->THETA2);

     //SUPG
     if(TDatabase::ParamDB->INTL_DISCTYPE==SUPG)
      {
       MatAdd(M_Intl, K_Intl,  -tau*TDatabase::TimeDB->THETA2);
      }
    } 
      
    // if(TDatabase::ParamDB->INTL_DISCTYPE==SUPG)
      // MatAdd(M_Intl, S_Intl,  1.);   //SUPG, time consistant term     

    //Set BC   
    this->SetDirichletBc();

    memset(B, 0, N_Dof*SizeOfDouble);

    // defect = M * oldsol
    MatVectActive(M_Intl,  Sol, B);
   
    /** usually zero, so no need dt*[theata3 + theta4] */
    Daxpy(N_Dof, tau, rhs, B);

    if(TDatabase::ParamDB->INTL_DISCTYPE!=DG)
    {
     if(cond_Lmin==DIRICHLET)
        B[0] = BDValue1;

     if(cond_Lmax==DIRICHLET)
       B[N_Dof -1] = BDValue2;
    }

   /** restore mass mat  values*/
   if(fabs(TDatabase::TimeDB->THETA2)>0)   
    memcpy(MatValues,  Mvalues_orig,  N_MatEntries*SizeOfDouble);

   // system matrix
    MatAdd(M_Intl, A_Intl, tau*TDatabase::TimeDB->THETA1);

    // add SUPG, time and space matrices 
    if(TDatabase::ParamDB->INTL_DISCTYPE==SUPG)
     { 
      MatAdd(M_Intl, S_Intl,  1.);
      MatAdd(M_Intl, K_Intl, tau*TDatabase::TimeDB->THETA1);
     }

    //Set BC   
    // if(TDatabase::ParamDB->INTL_DISCTYPE!=DG)
   this->SetDirichletBc();
 
  //print matrix
  // int *RowPtr = M_Intl->GetRowPtr();
  // int *KCol = M_Intl->GetKCol();
  // double *ValuesA = M_Intl->GetEntries(); 

   //solve the system
   DirectSolver(M_Intl, B, Sol);
     
  //  cout<<endl;
  // //  if(TDatabase::ParamDB->P2==1)
  //  {
    // for(j=0;j<N_Dof;j++)
    //  cout <<  j << ": "<<  Ddot(N_Dof, Sol, Sol) << " Newsol " << ": "<< Sol[j] <<endl;
    //  cout<<endl;
  //  }
  //  exit(0);
}


void TSystemADI1D::AssemblDG_Rhs(double Bnuc)
{

 int *GlobalNumbers, *BeginIndex, *DOF;
 double LocRhs[MaxN_BaseFunctions1D];

 TCollection *Coll;
 TBaseCell *Cell, *Neigh=NULL;
 FE1D FEId;
 TFE1D *Element;
 TBaseFunct1D *bf;

 double JointValues[MaxN_BaseFunctions1D];
 Coll = FESpace1D_Intl->GetCollection();
 GlobalNumbers = FESpace1D_Intl->GetGlobalNumbers();
 BeginIndex = FESpace1D_Intl->GetBeginIndex();
 int i, j, k, N_BaseFunct;

 
    i=0;
    Cell = Coll->GetCell(i);
    FEId = FESpace1D_Intl->GetFE1D(i, Cell);
    Element = TFEDatabase2D::GetFE1D(FEId);
    bf = Element->GetBaseFunct1D();
    N_BaseFunct = Element->GetN_DOF();
    Neigh = (Cell->GetJoint(i))->GetNeighbour(Cell);
    if(Neigh)
     {
      cout<<"L_min should be cell i=0, check SystemADI1D.C"<< endl;
      #ifdef _MPI
       MPI_Finalize();
      #endif
      exit(0); 
     }

    //find the current cell basis value at this joint (X_n+) L_min
    bf->GetDerivatives(D0, -1., JointValues);
    memset(LocRhs, 0, N_BaseFunct*SizeOfDouble);   

    for(k=0;k<N_BaseFunct;k++)
      LocRhs[k] += JointValues[k]*Bnuc;// Bnuc
   
     DOF = GlobalNumbers + BeginIndex[i];
     // add to global rhs
     for(j=0;j<N_BaseFunct;j++)
      {
       rhs[DOF[j]] += LocRhs[j]; 
      } // for(m=0;m<N_BaseFunct

    // print matrix
    // for(j=0;j<2;j++)
    //  {
    //   cout << j <<" f: " << rhs[j] <<endl;
    //  }
//   exit(0);
}

void TSystemADI1D::SetDirichletBc()
{
 int k;
 int *RowPtr, *KCol, begin, end;
 double *MatValues;

 RowPtr = M_Intl->GetRowPtr();
 KCol = M_Intl->GetKCol();
 MatValues = M_Intl->GetEntries(); 


    if(TDatabase::ParamDB->INTL_DISCTYPE!=DG)
     {
       if(cond_Lmin==DIRICHLET)
       {
        begin = RowPtr[0];
        end = RowPtr[1];

        for(k=begin;k<end;k++)
         {
          if(KCol[k] == 0 )
           { MatValues[k] = 1.; }
          else
           { MatValues[k] = 0.; }
         }
       }

      if(cond_Lmax==DIRICHLET)
       {
        begin = RowPtr[N_Dof-1];
        end = RowPtr[N_Dof];

        for(k=begin;k<end;k++)
         {
          if(KCol[k] == N_Dof-1)
           {
            MatValues[k] = 1.;
           }
          else
           {
            MatValues[k] = 0.;
           }
         }
       }
      }

}


// assemble mass matrix
void TSystemADI1D::AssembleMassMat()
{
  int N_Cells = Coll_Intl->GetN_Cells();
  TBaseCell *Cell;
  FE1D FEId;
  TFE1D *Element;
  TBaseFunct1D *bf;
  QuadFormula1D LineQuadFormula;
  TQuadFormula1D *qf1;
  TRefTrans1D *F_K;
  BaseFunct1D BaseFunct_ID, BaseFunct[1];
  TCollection *Coll;
  BoundCond BDType;
  int j, k, l, N_BaseFunct, TestDOF;
  double *Weights, *zeta, X[MaxN_QuadPoints_1D], AbsDetjk[MaxN_QuadPoints_1D];
  int N_Points, N_Sets=1, *DOF;
  bool Needs2ndDer[1];
  double test0, ansatz0, *orgD0, **origvaluesD0, Mult;
  double LocMatrixM[MaxN_BaseFunctions1D*MaxN_BaseFunctions1D];

  int *GlobalNumbers = FESpace1D_Intl->GetGlobalNumbers();
  int *BeginIndex = FESpace1D_Intl->GetBeginIndex();
  int *RowPtr = M_Intl->GetRowPtr();
  int *KCol = M_Intl->GetKCol();
  double *ValuesM = M_Intl->GetEntries();

  memset(ValuesM, 0, N_MatEntries*SizeOfDouble);

  for(int i=0; i<N_Cells; i++)
   {
    Cell = Coll_Intl->GetCell(i);
    FEId = FESpace1D_Intl->GetFE1D(i, Cell);
    Element = TFEDatabase2D::GetFE1D(FEId);
    bf = Element->GetBaseFunct1D();
    N_BaseFunct = Element->GetN_DOF();
    BaseFunct_ID = Element->GetBaseFunct1D_ID();
    l = bf->GetPolynomialDegree();
    LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(2*l);
    qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
    qf1->GetFormulaData(N_Points, Weights, zeta);
    F_K = TFEDatabase2D::GetRefTrans1D(LineAffin);
    ((TLineAffin *)F_K)->SetCell(Cell);
    ((TLineAffin *)F_K)->GetOrigFromRef(N_Points, zeta, X, AbsDetjk);
    BaseFunct[0] = BaseFunct_ID;
    ((TLineAffin *)F_K)->GetOrigValues(N_Sets, BaseFunct, N_Points, zeta,  LineQuadFormula,  Needs2ndDer);
    origvaluesD0=TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D0);
    DOF = GlobalNumbers + BeginIndex[i];

    memset(LocMatrixM, 0, N_BaseFunct*N_BaseFunct*SizeOfDouble);
    for(j=0;j<N_Points;j++)
    {
      Mult = Weights[j]*AbsDetjk[j];
      orgD0 = origvaluesD0[j];

      //cout<< " zeta[j] " << zeta[j]  <<endl;
      //len +=Mult;
      for(k=0;k<N_BaseFunct;k++)
      {
        test0  = orgD0[k];
        //cout<< " uref " << test0  <<endl;
        for(l=0;l<N_BaseFunct;l++)
        {
          ansatz0  = orgD0[l];
          LocMatrixM[k*N_BaseFunct + l] += (Mult*ansatz0*test0);
        }
      }
    }

    //   add to global matrices
    for(j=0;j<N_BaseFunct;j++)
    {
      TestDOF = DOF[j];
      for(k=RowPtr[TestDOF];k<RowPtr[TestDOF+1];k++)
      {
        for(l=0;l<N_BaseFunct;l++)
        {
          if(KCol[k] == DOF[l])
          {
            ValuesM[k] +=LocMatrixM[j*N_BaseFunct + l];
            break;
          }
        }                                         // for(m=0;m<N_BaseFunct_low
      }                                           // for(n=begin;n<end;n++)
    }       
   }// i

  //update boundary data
  if(cond_Lmin==DIRICHLET && !dGDisc)
  {
    for(k=RowPtr[0];k<RowPtr[1];k++)
    {
      if(KCol[k] == 0 )
        { ValuesM[k] = 1.; }
        else
          { ValuesM[k] = 0.; }
    }
  }     

  if(cond_Lmax==DIRICHLET && !dGDisc)
  {
    for(k=RowPtr[N_Dof-1];k<RowPtr[N_Dof];k++)
    {
      if(KCol[k] == N_Dof-1 )
        { ValuesM[k] = 1.; }
        else
          { ValuesM[k] = 0.; }
    }
  }     

 // store the original mass mat  values
 memcpy(Mvalues_orig, ValuesM, N_MatEntries*SizeOfDouble);

  // //print matrix
  //   for(j=0;j<N_Dof;j++)
  //   {
  //     for(k=RowPtr[j];k<RowPtr[j+1];k++)
  //      {
  //       cout << "M(" << j << ", "<< KCol[k] << ") = " << ValuesM[k] <<endl;
  //      }
  //     cout<<endl;
  //   }
  //  OutPut("Mass Matrix Assembling Done: "<< N_MatEntries<< endl);  
  // // //   exit(0);
 
}


// assemble advection matrix with advection 1
void TSystemADI1D::AssembleAdvectMat(bool Conservative)
{
 int i, j, k, l, m, Neigh_i, N_Cells_Internal, N_Joints, N_BaseFunct, Neigh_N_BaseFunct;
 int N_Points, N_Sets=1, *GlobalNumbers, *BeginIndex, *DOF, *NeibDOF;
 int TestDOF, begin, end, *RowPtr, *KCol; 

 double *Weights, *zeta, Intl_L[MaxN_QuadPoints_1D], AbsDetjk[MaxN_QuadPoints_1D];
 double LocMatrixA[MaxN_BaseFunctions1D*MaxN_BaseFunctions1D];
//  double LocMatrixA11[MaxN_BaseFunctions1D*MaxN_BaseFunctions1D];
 double LocMatrixA12[MaxN_BaseFunctions1D*MaxN_BaseFunctions1D];
 double LocMatrixA21[MaxN_BaseFunctions1D*MaxN_BaseFunctions1D];
 double LocMatrixA22[MaxN_BaseFunctions1D*MaxN_BaseFunctions1D];   
 double JointValues[MaxN_BaseFunctions1D], Neigh_JointValues[MaxN_BaseFunctions1D]; 
 double **origvaluesD0, **origvaluesD1, Mult;
 double val, *orgD0, *orgD1, *ValuesA;
 double Neigh_rec_detjk, h_max, h1, h2;
 bool Needs2ndDer[1], UpdateEdge=FALSE;

 TCollection *Coll;
 TBaseCell *Cell, *Neigh;
 FE1D FEId, Neigh_FEId;
 TFE1D *Element, *Neigh_Element;
 TBaseFunct1D *bf, *Neigh_bf;
 QuadFormula1D LineQuadFormula;
 TQuadFormula1D *qf1;
 TRefTrans1D *F_K;
 BaseFunct1D BaseFunct_ID, BaseFunct[1];
 double Epsilon = TDatabase::ParamDB->DG_P0;
 double Sigma0 = TDatabase::ParamDB->DG_P1;
//  Sigma1 = TDatabase::ParamDB->DG_P2;

  Coll = FESpace1D_Intl->GetCollection();
  GlobalNumbers = FESpace1D_Intl->GetGlobalNumbers();
  BeginIndex = FESpace1D_Intl->GetBeginIndex();

  RowPtr = A_Intl->GetRowPtr();
  KCol = A_Intl->GetKCol();
  ValuesA = A_Intl->GetEntries(); 

  N_Cells_Internal = Coll->GetN_Cells();
  Needs2ndDer[0] = FALSE;

  // double Conv = 1;

  // Advect Mat values_orig
  A_Intl->Reset();

  // associate each cell with her number in the collection
  for(i=0;i<N_Cells_Internal;i++)
  {
    Cell = Coll->GetCell(i);
    Cell->SetClipBoard(i);
  }

  for(i=0; i<N_Cells_Internal; i++)
  {
    UpdateEdge=FALSE;
    Cell = Coll->GetCell(i);
    FEId = FESpace1D_Intl->GetFE1D(i, Cell);
    Element = TFEDatabase2D::GetFE1D(FEId);
    bf = Element->GetBaseFunct1D();
    N_BaseFunct = Element->GetN_DOF();
    BaseFunct_ID = Element->GetBaseFunct1D_ID();
    //  cout<< "N_BaseFunct  " << N_BaseFunct << " MaxN_BaseFunctions1D " << MaxN_BaseFunctions1D <<endl;
    l = bf->GetPolynomialDegree();
    LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(20*l);
    qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
    qf1->GetFormulaData(N_Points, Weights, zeta);

    F_K = TFEDatabase2D::GetRefTrans1D(LineAffin);
    ((TLineAffin *)F_K)->SetCell(Cell);
    ((TLineAffin *)F_K)->GetOrigFromRef(N_Points, zeta, Intl_L, AbsDetjk);
    
    BaseFunct[0] = BaseFunct_ID;
    ((TLineAffin *)F_K)->GetOrigValues(N_Sets, BaseFunct, N_Points, zeta,  LineQuadFormula,  Needs2ndDer);

    origvaluesD0=TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D0);
    origvaluesD1=TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D1);

    memset(LocMatrixA, 0, N_BaseFunct*N_BaseFunct*SizeOfDouble);
    if(Conservative)
     {
     //  memset(LocMatrixA11, 0, N_BaseFunct*N_BaseFunct*SizeOfDouble);
     memset(LocMatrixA12, 0, N_BaseFunct*N_BaseFunct*SizeOfDouble);
     memset(LocMatrixA21, 0, N_BaseFunct*N_BaseFunct*SizeOfDouble);
     memset(LocMatrixA22, 0, N_BaseFunct*N_BaseFunct*SizeOfDouble);
    }

    DOF = GlobalNumbers + BeginIndex[i];
 
    for(j=0;j<N_Points;j++)
     {
      Mult = Weights[j]*AbsDetjk[j];
      orgD0 = origvaluesD0[j];
      orgD1 = origvaluesD1[j];
      for(k=0;k<N_BaseFunct;k++)
       {
      //   // test0  = orgD0[k];
      //   // test1  = orgD1[k]
        for(l=0;l<N_BaseFunct;l++)
         {
      //     // ansatz0  = orgD0[l];
      //     // ansatz1  = orgD1[l];
         if(Conservative)
         {
          val  = -orgD0[l]*orgD1[k]; // convective form   
         }
         else
         {
          val = orgD1[l]* orgD0[k];// non convective form
         }

         LocMatrixA[k*N_BaseFunct + l] += Mult*val;
        } //for(l=0;l<N_Ba
       } // for(k=0;k<N_BaseFunct;k+
    } // for(j=0;j<N_Points;j++)

   if(Conservative)
    {
     N_Joints=Cell->GetN_Edges();
     h1 = Cell->GetMeasure();
     for(m=0; m<N_Joints; m++)
      {
       Neigh = (Cell->GetJoint(m))->GetNeighbour(Cell);

       if(m==0 && Neigh) // inner joint, since in 1D we have only 2 joints (really vertices)
       {
        UpdateEdge = TRUE;
        // only first joint (really vertices) will be updated, 
        // other joint will be the first joint of the next cell
        h2 = Neigh->GetMeasure();
        h_max = MAX(h1, h2);
        
        //find the current cell basis value at this joint (X_n+)
        bf->GetDerivatives(D0, -1., JointValues);

        //find the neib cell basis value at this joint (X_n-)
        Neigh_i=Neigh->GetClipBoard();
        Neigh_FEId = FESpace1D_Intl->GetFE1D(Neigh_i, Neigh);
        Neigh_Element = TFEDatabase2D::GetFE1D(Neigh_FEId);
        Neigh_bf = Neigh_Element->GetBaseFunct1D();
        Neigh_N_BaseFunct = Neigh_Element->GetN_DOF();
        NeibDOF = GlobalNumbers + BeginIndex[Neigh_i];
        Neigh_bf->GetDerivatives(D0, 1., Neigh_JointValues);

        //(X_n+, X_n+) (test, ansatz)
        for(k=0;k<N_BaseFunct;k++)
         for(l=0;l<N_BaseFunct;l++)
         {
          // LocMatrixA[k*N_BaseFunct + l] += -JointValues[l]*JointValues[k];    // upwind  when b<0     
          //  cout << "E(" << DOF[k] << ", "<< DOF[l] << ") = " <<  LocMatrixA[k*N_BaseFunct + l] <<endl;
         }
   
        for(k=0;k<N_BaseFunct;k++) // own test X_n+
         for(l=0;l<Neigh_N_BaseFunct;l++) // neib ansatz X_n-
         {
          LocMatrixA12[k*N_BaseFunct + l] +=  -Neigh_JointValues[l]*JointValues[k]; // upwind  when b>0     
          //  cout << "E(" << DOF[k] << ", "<< NeibDOF[l] << ") = " <<  LocMatrixA12[k*N_BaseFunct + l] <<endl;
         }

        // d_n
        for(k=0;k<Neigh_N_BaseFunct;k++) // neib test X_n-
         for(l=0;l<N_BaseFunct;l++) // own ansatz X_n+
          {
          //  LocMatrixA21[k*N_BaseFunct + l] +=  Neigh_JointValues[k]*JointValues[l]; // upwind  when b<0     
          //  cout << "E(" << NeibDOF[k] << ", "<< DOF[l] << ") = " <<  LocMatrixA21[k*N_BaseFunct + l] <<endl;
         }

        for(k=0;k<Neigh_N_BaseFunct;k++) // neib test X_n-
         for(l=0;l<Neigh_N_BaseFunct;l++) // own ansatz X_n+
         {
           LocMatrixA22[k*N_BaseFunct + l] += Neigh_JointValues[k]*Neigh_JointValues[l];  // upwind  when b>0              
         }
        } // if(m==0 && Neigh)

      } //  for(m=0; m<N_Joi
    }// if(Conservative)

    // add to global matrices
     for(j=0;j<N_BaseFunct;j++)
     {
      for(k=RowPtr[DOF[j]];k<RowPtr[DOF[j]+1];k++)
       {
        for(l=0;l<N_BaseFunct;l++)
         {
          if(KCol[k] == DOF[l])
           {
            ValuesA[k] +=LocMatrixA[j*N_BaseFunct + l];
            break;
           }
          } // for(m=0;m<N_BaseFunct

       if(UpdateEdge)
        {
        for(l=0;l<Neigh_N_BaseFunct;l++)
         {
          if(KCol[k] == NeibDOF[l])
           {
            ValuesA[k] +=LocMatrixA12[j*N_BaseFunct + l];
            break;
           } 
         } // for(l=0;l<Neigh_N_BaseFunct;l++)
        } //if(Conservativ
       } // for(n=begin;n<end;n++)

     } //  for(j=0;j<N_BaseFunct;j++)   

     if(UpdateEdge)
      {
       for(j=0;j<Neigh_N_BaseFunct;j++)
        {
         TestDOF = NeibDOF[j];
         begin = RowPtr[TestDOF];
         end = RowPtr[TestDOF+1];

         for(k=begin;k<end;k++)
          {
           for(l=0;l<N_BaseFunct;l++)
            {
             if(KCol[k] == DOF[l])
             {
              ValuesA[k] +=LocMatrixA21[j*N_BaseFunct + l];
              //  cout << TestDOF << " " << KCol[k] << endl;
              break;
             }
           } //  for(l=0;l<N_BaseFunct;l++)

           for(l=0;l<Neigh_N_BaseFunct;l++)
           {
            if(KCol[k] == NeibDOF[l])
             {
              ValuesA[k] +=LocMatrixA22[j*N_BaseFunct + l];
              //  cout << TestDOF << " " << KCol[k] <<  " " <<  LocMatrixA22[j*N_BaseFunct + l] << endl;
              break;
             }
            } //  for(l=0;l<N_BaseFunct;l++)
          } // for(k=begin;k<end;k++)
        } // for(j=0;j<Neigh_N_BaseFunct;j++)
      } //  if(UpdateEdge)

// if(i==1)
// {
// //  print matrix
//     for(j=0;j<N_Dof;j++)
//      {
//       begin = RowPtr[j];
//       end = RowPtr[j+1];
//       for(k=begin;k<end;k++)
//        {
//         cout << "A(" << j << ", "<< KCol[k] << ") = " << ValuesA[k] <<endl;
//        }
//        cout<<endl;
//       } 
  
//   exit(0);
// }

   }// for(i=0; i<N_Cells_Internal; i++)

//  store the original mass mat  values
 memcpy(Advectvalues_orig, ValuesA, N_MatEntries*SizeOfDouble);

// //  print matrix
//     for(j=0;j<N_Dof;j++)
//      {
//       begin = RowPtr[j];
//       end = RowPtr[j+1];
//       for(k=begin;k<end;k++)
//        {
//         cout << "A(" << j << ", "<< KCol[k] << ") = " << ValuesA[k] <<endl;
//        }
//        cout<<endl;
//       } 
      
//       exit(0);

}

int TSystemADI1D::Nodal2DOF(double *Sol_XposLNnodLOwnDof)
{
  int k, N_OuterNodal = N_Xpos;
  int N_OwnNod = N_LnodalPos[OwnADI_Idx];
  for(int i=0;i<N_ADISystems;i++)
   {
    if(i==OwnADI_Idx) continue;
    N_OuterNodal *=N_LnodalPos[i];
   }

  if(DofSameAsNodalPts)
  {
   memcpy(Sol_XposLNnodLOwnDof, Sol_XposLNnodal, SizeOfDouble*N_Dof*N_OuterNodal);
   return 0; 
  }
   
  int N_NodalPoints, N_LocalDOFs, N_Cells = Coll_Intl->GetN_Cells();
  int *IndexOfNodalPts = FESpace1D_Intl->GetIntlPtIndexOfPts();
  int *DOF, *GlobalNumbers, *BeginIndex, *Index, disp=0;
  TBaseCell *cell;
  FE1D FEId;
  TFE1D *Element;
  TNodalFunctional1D *nf;
  double *xi, *eta, *Sol_Nodalown;
  double X[MaxN_PointsForNodal1D], Y[MaxN_PointsForNodal1D];
  double FunctionalValues[MaxN_PointsForNodal1D];
  BeginIndex = FESpace1D_Intl->GetBeginIndex();
  GlobalNumbers = FESpace1D_Intl->GetGlobalNumbers();

  // cout<< "N_OuterNodal: " << N_OuterNodal <<endl;

  double *sol_Dofloc, *PtVal, *PointValues = new double [MaxN_PointsForNodal1D*N_OuterNodal];
  int i, ii, j;

  memset(Sol_XposLNnodLOwnDof, 0, SizeOfDouble*N_Dof*N_OuterNodal);
  memset(IncidentArray, 0, SizeOfInt*N_Dof);


  for(i=0;i<N_Cells;i++)
   {
    cell = Coll_Intl->GetCell(i);
    FEId = FESpace1D_Intl->GetFE1D(i, cell);
    Element = TFEDatabase2D::GetFE1D(FEId);
    nf = Element->GetNodalFunctional1D();
    nf->GetPointsForAll(N_NodalPoints, xi, eta);
    N_LocalDOFs = Element->GetN_DOF();
    DOF = GlobalNumbers+BeginIndex[i];
    Index = IndexOfNodalPts+disp;

    for(ii=0;ii<N_OuterNodal; ++ii)
     {
      Sol_Nodalown  =  Sol_XposLNnodal + ii*N_OwnNod;

      for(j=0;j<N_NodalPoints;++j)
      {
       k = Index[j]; // corresponding L cell
       PointValues[ii*N_NodalPoints + j] = Sol_Nodalown[k];
      }
     }

    for(ii=0;ii<N_OuterNodal;++ii)
     {
      PtVal = PointValues + ii*N_NodalPoints;
      nf->GetAllFunctionals(PtVal, FunctionalValues);

      sol_Dofloc = Sol_XposLNnodLOwnDof + ii*N_Dof;

      for(j=0;j<N_LocalDOFs;j++)
       {
        sol_Dofloc[DOF[j]] += FunctionalValues[j];

        if(ii==0)
         IncidentArray[DOF[j]] ++;
       }
     }
    disp +=N_NodalPoints;
   } // i

  for(ii=0;ii<N_OuterNodal;ii++)
    {
     for(i=0;i<N_Dof;i++)
      {
       if(ii==0)
        if(IncidentArray[i] == 0)
         {
          cout << "Error in L_Nodal2Sol : "<< IncidentArray[i] << endl;
          exit(0);
         }
       Sol_XposLNnodLOwnDof[ii*N_Dof + i] /= (double)IncidentArray[i];
      } // for(i=0;i<N_DOFs;i
    } // for(ii=0;ii<N_XLocPoints;i

  delete [] PointValues;
 
  //  int l=N_LnodalPos[2], m = N_LnodalPos[1];

  //  if(TDatabase::ParamDB->P2==1)
  //  {
//     for(i=0; i<N_Xpos; ++i)
//      for(j=0; j<l; ++j)
//        for(k=0; k<m; ++k)
//         for(l=0;l<N_Dof;l++)
//      cout << l << " Newsol " << ": "<< Sol_XposLNnodLOwnDof[ i*l*m*N_Dof + j*m*N_Dof + k*N_Dof+ l] <<endl;
//      cout<<endl;
//   //  }

// exit(0);

  return 0;
} //Nodal2DOF


int TSystemADI1D::DOF2Nodal(double *Sol_XposLNnodLOwnDof)
{
  int N_OuterNodal = N_Xpos;
  int N_OwnNod = N_LnodalPos[OwnADI_Idx];  
  for(int i=0;i<N_ADISystems;i++)
   {
    if(i==OwnADI_Idx) continue;
    N_OuterNodal *=N_LnodalPos[i];
   }
  // cout<< "N_OuterNodal: " << N_OuterNodal <<endl;
  if(DofSameAsNodalPts)
  {
   memcpy(Sol_XposLNnodal, Sol_XposLNnodLOwnDof, SizeOfDouble*N_Dof*N_OuterNodal);
   return 0; 
  }
   
  int N_NodalPoints, N_LocalDOFs, N_Cells = Coll_Intl->GetN_Cells();
  int *IndexOfNodalPts = FESpace1D_Intl->GetIntlPtIndexOfPts();
  int *DOF, *GlobalNumbers, *BeginIndex, *Index, disp=0;
  TBaseCell *cell;
  FE1D FEId;
  TFE1D *Element;
  TNodalFunctional1D *nf;
  double *xi, *eta, *Sol_Nodalown, val;
  double X[MaxN_PointsForNodal1D];
  double AbsDetjk[MaxN_PointsForNodal1D];
  double BasisValues[MaxN_PointsForNodal1D][MaxN_BaseFunctions1D];  
  double FunctionalValues[MaxN_PointsForNodal1D];
  BeginIndex = FESpace1D_Intl->GetBeginIndex();
  GlobalNumbers = FESpace1D_Intl->GetGlobalNumbers();
  TRefTrans1D *rt;

  TFE1D *FE_Obj;
  TBaseFunct1D *bf;
  double *sol_Dofloc, *sol_Nodalloc;
  int i, ii, j, k, l;

  memset(IncidentArray, 0, SizeOfInt*N_OwnNod);
  memset(Sol_XposLNnodal, 0, SizeOfDouble*N_OwnNod*N_OuterNodal);

  for(i=0;i<N_Cells;i++)
   {
    cell = Coll_Intl->GetCell(i);
    FEId = FESpace1D_Intl->GetFE1D(i, cell);
    Element = TFEDatabase2D::GetFE1D(FEId);
    FE_Obj = TFEDatabase2D::GetFE1D(FEId);
    bf = FE_Obj->GetBaseFunct1D();
    nf = Element->GetNodalFunctional1D();
    nf->GetPointsForAll(N_NodalPoints, xi, eta);
    N_LocalDOFs = Element->GetN_DOF();
    DOF = GlobalNumbers+BeginIndex[i];
    Index = IndexOfNodalPts+disp;
    rt = TFEDatabase2D::GetRefTrans1D(LineAffin);
    ((TLineAffin *)rt)->SetCell(cell);
    ((TLineAffin *)rt)->GetOrigFromRef(N_NodalPoints, xi, X, AbsDetjk);

     for(j=0;j<N_NodalPoints;j++)
       bf->GetDerivatives(D0, xi[j], BasisValues[j]);

     for(ii=0;ii<N_OuterNodal;++ii)
      {
       sol_Dofloc = Sol_XposLNnodLOwnDof + ii*N_Dof;
       sol_Nodalloc =  Sol_XposLNnodal + ii*N_OwnNod; 
      //  val2 = 0;
       for(j=0;j<N_NodalPoints;++j)
       {
        val = 0.;  
        for(l=0;l<N_LocalDOFs;l++)
        {
         val += sol_Dofloc[DOF[l]]*BasisValues[j][l];
        //  val2 +=BasisValues[j][l];
        }

        k = Index[j]; // corresponding L cell
        sol_Nodalloc[k] += val;
        
        if(ii==0)
         IncidentArray[k]++;
       }
      }
    disp +=N_NodalPoints;

    // cout << " val2 " << val2 << endl;
  //  exit(0);
   }// 


   for(ii=0;ii<N_OuterNodal;++ii)
    {
     sol_Nodalloc =  Sol_XposLNnodal + ii*N_OwnNod; 
     for(i=0;i<N_OwnNod;i++)
      {
       if(ii==0)
        if(IncidentArray[i] == 0)
         {
          cout << "Error in L_Sol2Nodal : "<< IncidentArray[i] << endl;
          exit(0);
         }

       sol_Nodalloc[i] /= (double)IncidentArray[i];
      }
    }

 return 0;     
} // DOF2Nodal


// nodal points to evaluate nodal values
void TSystemADI1D::GenerateNodalPts()
{
  int i, j, k, l, m, r, N_Cells, *RootPtIndex;
  int N_LocalDOFs, N_Points;
  int N_AllLocalPoints;

  double L, L0, *xi, *eta, *L_loc, *L_loc_origOrder;
  double X[MaxN_PointsForNodal1D], Y[MaxN_PointsForNodal1D];
  double AbsDetjk[MaxN_PointsForNodal1D];

  TBaseCell *cell;
  FE1D FEId;
  TFE1D *Element;
  TFE1D *FE_Obj;
  TNodalFunctional1D *nf;
  TRefTrans1D *rt;
  TBaseFunct1D *bf;

  N_Cells = Coll_Intl->GetN_Cells();

  N_AllLocalPoints = 0;
  for(i=0; i<N_Cells; i++)
  {
    cell = Coll_Intl->GetCell(i);
    FEId = FESpace1D_Intl->GetFE1D(i, cell);
    Element = TFEDatabase2D::GetFE1D(FEId);
    FE_Obj = TFEDatabase2D::GetFE1D(FEId);
    bf = FE_Obj->GetBaseFunct1D();
    nf = Element->GetNodalFunctional1D();
    nf->GetPointsForAll(N_Points, xi, eta);
    N_LocalDOFs = Element->GetN_DOF();
    N_AllLocalPoints +=N_Points;
  }

  L_loc = new double [N_AllLocalPoints];
  L_loc_origOrder = new double [N_AllLocalPoints];
  N_AllLocalPoints = 0;
  for(i=0; i<N_Cells; i++)
  {
    cell = Coll_Intl->GetCell(i);
    FEId = FESpace1D_Intl->GetFE1D(i, cell);
    Element = TFEDatabase2D::GetFE1D(FEId);
    FE_Obj = TFEDatabase2D::GetFE1D(FEId);
    bf = FE_Obj->GetBaseFunct1D();
    nf = Element->GetNodalFunctional1D();
    nf->GetPointsForAll(N_Points, xi, eta);
    N_LocalDOFs = Element->GetN_DOF();
    rt = TFEDatabase2D::GetRefTrans1D(LineAffin);
    ((TLineAffin *)rt)->SetCell(cell);
    ((TLineAffin *)rt)->GetOrigFromRef(N_Points, xi, X, Y, AbsDetjk);

    for(j=0; j<N_Points; j++)
    {
      L_loc[N_AllLocalPoints] = X[j];
      N_AllLocalPoints++;
    }
  } // for(i=0; i<N_Cells; i++)

  memcpy(L_loc_origOrder, L_loc,  N_AllLocalPoints*SizeOfDouble);

  for(i=0; i<N_AllLocalPoints-1; i++)
   for(j=i; j<N_AllLocalPoints; j++)
    if(L_loc[i]> L_loc[j])
     {
      L= L_loc[i];
      L_loc[i]= L_loc[j];
      L_loc[j]= L;
     }

  L  = L_loc[0];
  N_NodalPts = 1;

  for(i=1; i<N_AllLocalPoints; i++)
   {
    if( fabs(L_loc[i]-L)>1e-5 )
     {
      N_NodalPts++;
      L = L_loc[i];
     }
   }

  NodalPt_Coord= new double[N_AllLocalPoints];
  NodalPt_Coord[0] = L_loc[0];
  N_NodalPts = 1;
  L  = L_loc[0];

  for(i=1; i<N_AllLocalPoints; i++)
   {
    if( fabs(L_loc[i]-L)>1e-5 )
     {
      NodalPt_Coord[N_NodalPts] = L_loc[i];
      N_NodalPts++;
      L = L_loc[i];
     }
   }

  delete [] L_loc;

  RootPtIndex = new int[N_AllLocalPoints];

  // find the index for the local points in the root points
  for(i=0; i<N_AllLocalPoints; i++)
   {
    L = L_loc_origOrder[i];
    l=0;
    r=N_NodalPts;

    m = N_NodalPts/2;
    L0 = NodalPt_Coord[m];

    while(fabs(L-L0) > 1.e-8 )
     {
      if(L < L0)  //poin lies in left side
       {
        r = m;
       }
      else
       {
        l=m;
       }

      m= (l+r)/2;
      L0 = NodalPt_Coord[m];
     } //  while ( 

    RootPtIndex[i] = m;
   }

  FESpace1D_Intl->SetIntlPtIndexOfPts(RootPtIndex);
  FESpace1D_Intl->SetN_RootNodalPts(N_NodalPts);

  delete [] L_loc_origOrder;
 //   N_L1nodal = N_NodalPts; 

 //   cout << N_AllLocalPoints << " N_NodalPts  "  << N_NodalPts << endl;
 //    for(i=0; i<N_NodalPts; i++)
 //    cout << i << " L: "  << NodalPt_Coord[i] << endl;
 //  exit(0);
}

 
void TSystemADI1D::Generate1DMesh(double Start, double End, int N_Cells)
{
  int i, j, N_Vert;
  int *Lines;
  double len, h, x, y, *X;
  double hmin, hmax;
  TVertex **Vetrex;
  TJoint *Joint;
  TBaseCell  **CellTree;

  N_Vert = N_Cells+1;
  X = new double[N_Vert];

  h = (End-Start)/(double)N_Cells;

  X[0] = Start;

  for(i=1; i<N_Vert; i++)
   X[i] = X[i-1] + (double)h;

  X[N_Vert-1] = End;

    hmin = 1.e8; hmax = -1.e8; 
    for(i=0; i<N_Vert-1; i++)
     {
      len = sqrt ((X[i+1] - X[i])*(X[i+1] - X[i]));
      if(len< hmin) hmin = len;
      if(len> hmax) hmax = len;        
     }
     OutPut("L h_min : " << hmin << " L h_max : " << hmax << endl);
    
//    for(i=0; i<N_Vert; i++)
//     cout<< i << " X[i] " << X[i] <<endl;
 
//  exit(0);

  Lines = new int[2*N_Cells];
  Vetrex = new TVertex*[N_Vert]; 

  y=0.;

  for(i=0; i<N_Cells; i++)
   {
    Lines[2*i]=i;
    Lines[2*i+1]=i+1;
#ifdef __2D__
    Vetrex[i] = new TVertex(X[i], y);
#else
  cout << "Not Yet Implemented " <<endl;
#endif
   }
   
#ifdef __2D__
  Vetrex[N_Cells] = new TVertex(X[N_Vert-1], y);
#endif
  CellTree = new TBaseCell*[N_Cells];

   for (i=0;i<N_Cells;i++)
   {
 //     Vetrex[ i ]->GetCoords(x, y);
 //     cout<< " x " << x<< " y " << y<<endl;
    CellTree[i] = new TMacroCell(TDatabase::RefDescDB[S_Line], 0);
    CellTree[i]->SetVertex(0, Vetrex[ Lines[ 2*i       ]]);
    CellTree[i]->SetVertex(1, Vetrex[ Lines[ 2*i + 1]]);
    ((TMacroCell *) CellTree[i])->SetSubGridID(0);
   }
 //     Vetrex[ i ]->GetCoords(x, y);
 //     cout<< " x " << x<< " y " << y<<endl;
 //     exit(0);

   Domain_Intl->SetTreeInfo(CellTree, N_Cells);

   TDatabase::IteratorDB[It_EQ]->SetParam(Domain_Intl);
   TDatabase::IteratorDB[It_LE]->SetParam(Domain_Intl);
   TDatabase::IteratorDB[It_Finest]->SetParam(Domain_Intl);
   TDatabase::IteratorDB[It_Between]->SetParam(Domain_Intl);
   TDatabase::IteratorDB[It_OCAF]->SetParam(Domain_Intl);

   // start joint(vertex)
   Joint = new TJointEqN(CellTree[0]);
   CellTree[0]->SetJoint(0, Joint);


   for(i=1;i<N_Cells;i++)
    {
     Joint = new TJointEqN(CellTree[i-1], CellTree[i]);

     CellTree[i-1]->SetJoint(1, Joint);
     CellTree[i]->SetJoint(0, Joint);
   } // for(i=0;i<N_Cells;i++)

   // end joint(vertex)
   Joint = new TJointEqN(CellTree[N_Cells-1]);
   CellTree[N_Cells-1]->SetJoint(1, Joint);

  delete []  Lines;
}

// /**  Assembling A matrix Ilen-Southwell finite difference matrix*/
// /** mass mat is same for all quad points in this cell */
void TSystemADI1D::AssembleARhs_FD(double Conv)
{
  int i, j, N;
  int begin, end, *RowPtr, *KCol; 
  double *ValuesA, h, eps, q, a, b, c, beta;

  if(fabs(TDatabase::ParamDB->REACTOR_P3)>0)
   { eps = 1.0/TDatabase::ParamDB->REACTOR_P3;}
  else
   { eps = 1.e-15;}
  
  RowPtr = A_Intl->GetRowPtr();
  KCol = A_Intl->GetKCol();
  ValuesA = A_Intl->GetEntries(); 
  N = A_Intl->GetN_Columns();

   for(i=0;i<N;i++)
    {
      if(i<(N-1))
       { h = NodalPt_Coord[i+1] - NodalPt_Coord[i]; } //note h_i has to be calculated base of upwind, but same for uniform mesh
      
      // 2nd order Ilin-Allen Southwell scheme
      if(fabs(Conv)>0.)
       {
        q = Conv*h/(2.*eps);
        beta = eps*(q/tanh(q));
       }
       
      begin=RowPtr[i];
      end=RowPtr[i+1];

      //BD 
      a = -beta/(h*h) - Conv/(2*h);
      b = 2*beta/(h*h);
      c =  -beta/(h*h) + Conv/(2*h);

      for(j=begin;j<end;j++)
       if( i==(KCol[j]+1) ) //lower 
        { 
          ValuesA[j] = a; 
        }
       else if(i==(KCol[j]))//diagonal 
        { 
          ValuesA[j] = b;    
        }
       else if(i==(KCol[j]-1)) //upper 
        { 
          ValuesA[j] = c;        
        }        
    //  cout<<endl;
    }     

  // //print matrix
  //  for(j=0;j<N_Dof;j++)
  //   {
  //    begin = RowPtr[j];
  //    end = RowPtr[j+1];
  //    for(int k=begin;k<end;k++)
  //     {
  //      cout << "A(" << j << ", "<< KCol[k] << ") = " << ValuesA[k] <<endl;
  //     }
       
  //     //  cout << "f: " << rhs[j ] <<endl;
  //    cout<<endl;
  //   }
//  exit(0);


}

// /**  Assembling A matrix */
// /** mass mat is same for all quad points in this cell */
void TSystemADI1D::AssembleARhs(double *X, double Conv, CoeffFctND *Bilinear)
{
 int i, j, k, l, N_Cells_Internal, N_BaseFunct;
 int N_Points, N_Sets=1, *GlobalNumbers, *BeginIndex, *DOF;
 int TestDOF, begin, end, *RowPtr, *KCol;

 double *Weights, *zeta, *Intl_L, AbsDetjk[MaxN_QuadPoints_1D];
 double LocMatrixA[MaxN_BaseFunctions1D*MaxN_BaseFunctions1D];
 double LocRhs[MaxN_BaseFunctions1D];
 double **origvaluesD0, **origvaluesD1, Mult;
 double *orgD0, *orgD1, test0, test1, ansatz0, ansatz1, *ValuesA;
 double c0, c1, g0, rhsval, val, len=0.;
 double **aux, **coeff, *Coeff;
 double **Coords;

 bool Needs2ndDer[1];

 TCollection *Coll;
 TBaseCell *Cell;
 FE1D FEId;
 TFE1D *Element;
 TBaseFunct1D *bf;
 QuadFormula1D LineQuadFormula;
 TQuadFormula1D *qf1;
 TRefTrans1D *F_K;
 BaseFunct1D BaseFunct_ID, BaseFunct[1];

  Coll = FESpace1D_Intl->GetCollection();
  GlobalNumbers = FESpace1D_Intl->GetGlobalNumbers();
  BeginIndex = FESpace1D_Intl->GetBeginIndex();

  RowPtr = A_Intl->GetRowPtr();
  KCol = A_Intl->GetKCol();
  ValuesA = A_Intl->GetEntries(); 

  N_Cells_Internal = Coll->GetN_Cells();
  Needs2ndDer[0] = FALSE;

  aux = new double *[MaxN_QuadPoints_1D];
  coeff = new double *[MaxN_QuadPoints_1D];
  Coords = new double* [N_Coord+1];

  for(i=0; i<N_Coord+1; i++)  
   Coords[i] = new double [MaxN_QuadPoints_1D];

  for(i=0; i<MaxN_QuadPoints_1D; i++)
   {
    for(int ii=0; ii<N_Coord; ++ii )
    Coords[ii][i] = X[ii];

    aux[i] = new double[6];
    coeff[i] = new double[6];
   }

   Intl_L = Coords[N_Coord];
 
   for(i=0; i<N_Cells_Internal; i++)
   {
    Cell = Coll->GetCell(i);
    FEId = FESpace1D_Intl->GetFE1D(i, Cell);
    Element = TFEDatabase2D::GetFE1D(FEId);
    bf = Element->GetBaseFunct1D();
    N_BaseFunct = Element->GetN_DOF();
    BaseFunct_ID = Element->GetBaseFunct1D_ID();

    l = bf->GetPolynomialDegree();
    LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(2*l);
    qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
    qf1->GetFormulaData(N_Points, Weights, zeta);

    F_K = TFEDatabase2D::GetRefTrans1D(LineAffin);
    ((TLineAffin *)F_K)->SetCell(Cell);
    ((TLineAffin *)F_K)->GetOrigFromRef(N_Points, zeta, Intl_L, AbsDetjk);

    Bilinear(N_Points, N_Coord+1, Coords, aux, coeff); 

    BaseFunct[0] = BaseFunct_ID;
    ((TLineAffin *)F_K)->GetOrigValues(N_Sets, BaseFunct, N_Points, zeta,  LineQuadFormula,  Needs2ndDer);

    origvaluesD0=TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D0);
    origvaluesD1=TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D1);

    memset(LocMatrixA, 0, N_BaseFunct*N_BaseFunct*SizeOfDouble);
    memset(LocRhs, 0, N_BaseFunct*SizeOfDouble);

    DOF = GlobalNumbers + BeginIndex[i];

    for(j=0;j<N_Points;j++)
     {
      Coeff = coeff[j];
      c0 = Coeff[0]; // diffusion
      g0 = Coeff[1]; // convection in z direction
      c1 = Coeff[2]; // reaction term
      rhsval = Coeff[3]; //rhs

      // in case of PBE, it is the growth term from concentration
      // in PBE example file define Coeff[1] as < 0, so Conv will be the growth term
      if(g0>=0)
        Conv = g0;

      // cout<< " c  " << c0 << " Conv " <<Conv << " c1 " << c1 <<endl;   
      Mult = Weights[j]*AbsDetjk[j];
      orgD0 = origvaluesD0[j];
      orgD1 = origvaluesD1[j];

      len +=Mult;
      for(k=0;k<N_BaseFunct;k++)
       {
        test0  = orgD0[k];
        test1  = orgD1[k];
        // cout<< " test0  " << test0 << " orgD1  " <<test1 << " AbsDetjk " << AbsDetjk[j] <<endl;   
        LocRhs[k] += Mult*rhsval*test0;

        for(l=0;l<N_BaseFunct;l++)
         {
          ansatz0  = orgD0[l];
          ansatz1  = orgD1[l];

          // val  = Conv*ansatz1*test0;// convective form
          val  = -Conv*ansatz0*test1; // conservative form         
          val += c0*ansatz1*test1; // difusive term
          val += c1*ansatz0*test0;  
          //  cout<< "A ( " << k <<" , " <<l << " ) " <<Mult*val <<endl;
          LocMatrixA[k*N_BaseFunct + l] += (Mult*val);
        } 
       } 
     }
 //      cout<<endl;


   //update boundary data
  
    // add to global matrices
    for(j=0;j<N_BaseFunct;j++)
     {
      TestDOF = DOF[j];
      rhs[TestDOF] += LocRhs[j]; 
      
      begin = RowPtr[TestDOF];
      end = RowPtr[TestDOF+1];
      for(k=begin;k<end;k++)
       {
       for(l=0;l<N_BaseFunct;l++)
        {
         if(KCol[k] == DOF[l])
          {
           ValuesA[k] +=LocMatrixA[j*N_BaseFunct + l];  
           break;
          }
        } // for(m=0;m<N_BaseFunct_low
      } // for(n=begin;n<end;n++)

    //  //update boundary flux
    //  if(TestDOF==0 && cond_Lmin==NEUMANN)
    //   {
    //   for(k=begin;k<end;k++)
    //     {
    //      if(KCol[k] == 0 )
    //       { ValuesA[k] += -1; }
    //   //    else
    //   //     { ValuesM[k] = 0.; }
    //    }
    //   } //  if(cond_Lmin==DIRIC 

     } // for(l=0;l<N_BaseFunct_low
   }// for(i=0; i<N_Cells_Internal; i++)
   
//  cout<< " len " << len << endl;
// //  exit(0);
//   //print matrix
//    for(j=0;j<N_Dof;j++)
//     {
//      begin = RowPtr[j];
//      end = RowPtr[j+1];
//      for(k=begin;k<end;k++)
//       {
//        cout << "A(" << j << ", "<< KCol[k] << ") = " << ValuesA[k] <<endl;
//       }
       
//        cout << "f: " << rhs[j ] <<endl;
//     //  cout<<endl;
//     }
//  exit(0);
  for(i=0; i<MaxN_QuadPoints_1D; i++)
   {
    delete [] aux[i];
    delete [] coeff[i];
   }
   
  delete [] aux;
  delete [] coeff;

  
  for(i=0; i<N_Coord+1; i++)  
   delete [] Coords[i]; 
  
 delete [] Coords;  
  
 //  printf("    AssembleARhs  \n" );
 //   MPI_Finalize();
 // exit(0);
} // void TSystemADI1D::AssembleARhs(d


// /**  Assembling A matrix */
// /** mass mat is same for all quad points in this cell */
void TSystemADI1D::AssembleARhs_SUPG(double *X, double Conv, CoeffFctND *Bilinear)
{
 int i, j, k, l, N_Cells_Internal, N_BaseFunct;
 int N_Points, N_Sets=1, *GlobalNumbers, *BeginIndex, *DOF;
 int TestDOF, begin, end, *RowPtr, *KCol;

 double *Weights, *zeta, *Intl_L, AbsDetjk[MaxN_QuadPoints_1D];
 double LocMatrixA[MaxN_BaseFunctions1D*MaxN_BaseFunctions1D];
 double LocMatrixS[MaxN_BaseFunctions1D*MaxN_BaseFunctions1D];
 double LocMatrixK[MaxN_BaseFunctions1D*MaxN_BaseFunctions1D]; 
 double LocRhs[MaxN_BaseFunctions1D];
 double **origvaluesD0, **origvaluesD1, **origvaluesD2, Mult;
 double *orgD0, *orgD1, *orgD2, test0, test1, ansatz0, ansatz1, ansatz2, *ValuesA, *ValuesS, *ValuesK;
 double c0, c1, g0, rhsval, val, len=0., bgradv, hE, beta, Pe_K;
 double **aux, **coeff, *Coeff, **Coords;
 
 bool Needs2ndDer[1];

 TCollection *Coll;
 TBaseCell *Cell;
 FE1D FEId;
 TFE1D *Element;
 TBaseFunct1D *bf;
 QuadFormula1D LineQuadFormula;
 TQuadFormula1D *qf1;
 TRefTrans1D *F_K;
 BaseFunct1D BaseFunct_ID, BaseFunct[1];
 TVertex *Vetrex1, *Vetrex2;

  double D_L =  TDatabase::ParamDB->REACTOR_P3;
  // double delta0 = TDatabase::ParamDB->DELTA0;
  // double delta1 = TDatabase::ParamDB->DELTA1;
  
  if(D_L<1.e-12) D_L = 1.e-12;
  
  Coll = FESpace1D_Intl->GetCollection();
  GlobalNumbers = FESpace1D_Intl->GetGlobalNumbers();
  BeginIndex = FESpace1D_Intl->GetBeginIndex();

  RowPtr = A_Intl->GetRowPtr();
  KCol = A_Intl->GetKCol();
  ValuesA = A_Intl->GetEntries(); 
  ValuesS = S_Intl->GetEntries();
  ValuesK = K_Intl->GetEntries();
  
  N_Cells_Internal = Coll->GetN_Cells();
  Needs2ndDer[0] = TRUE;

  aux = new double *[MaxN_QuadPoints_1D];
  coeff = new double *[MaxN_QuadPoints_1D];

  Coords = new double* [N_Coord+1]; 

  for(i=0; i<N_Coord+1; i++)  
   Coords[i] = new double [MaxN_QuadPoints_1D];
   
  for(i=0; i<MaxN_QuadPoints_1D; i++)
   {
    for(int ii=0; ii<N_Coord; ++ii )
    Coords[ii][i] = X[ii];

    aux[i] = new double[6];
    coeff[i] = new double[6];
   }

   Intl_L = Coords[N_Coord];

  for(i=0; i<N_Cells_Internal; i++)
   {
    Cell = Coll->GetCell(i);
    FEId = FESpace1D_Intl->GetFE1D(i, Cell);
    Element = TFEDatabase2D::GetFE1D(FEId);
    bf = Element->GetBaseFunct1D();
    N_BaseFunct = Element->GetN_DOF();
    BaseFunct_ID = Element->GetBaseFunct1D_ID();

    l = bf->GetPolynomialDegree();
    LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(20*l);
    qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
    qf1->GetFormulaData(N_Points, Weights, zeta);

    F_K = TFEDatabase2D::GetRefTrans1D(LineAffin);
    ((TLineAffin *)F_K)->SetCell(Cell);
    ((TLineAffin *)F_K)->GetOrigFromRef(N_Points, zeta, Intl_L, AbsDetjk);

 //     Bilinear(N_Points, Space_X, Space_Y, Intl_L, aux, coeff);
    Bilinear(N_Points, N_Coord+1, Coords, aux, coeff);
    
    BaseFunct[0] = BaseFunct_ID;
    ((TLineAffin *)F_K)->GetOrigValues(N_Sets, BaseFunct, N_Points, zeta,  LineQuadFormula,  Needs2ndDer);

    origvaluesD0=TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D0);
    origvaluesD1=TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D1);
    origvaluesD2=TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D2);

    memset(LocMatrixA, 0, N_BaseFunct*N_BaseFunct*SizeOfDouble);
    memset(LocMatrixS, 0, N_BaseFunct*N_BaseFunct*SizeOfDouble);
    memset(LocMatrixK, 0, N_BaseFunct*N_BaseFunct*SizeOfDouble);
    memset(LocRhs, 0, N_BaseFunct*SizeOfDouble);

    DOF = GlobalNumbers + BeginIndex[i];
    Vetrex1 = Cell->GetVertex(0);
    Vetrex2 = Cell->GetVertex(1);
    hE = fabs(Vetrex1->GetX() - Vetrex2->GetX());
 //     cout<< " hE  " << hE  <<endl;

    for(j=0;j<N_Points;j++)
     {
      Coeff = coeff[j];
      c0 = Coeff[0]; // diffusion
      g0 = Coeff[1]; // convection in z direction
      c1 = Coeff[2]; // reaction term
      rhsval = Coeff[3]; //rhs

      // in case of PBE, it is the growth term from concentration
      // in PBE example file define Coeff[1] as < 0, so Conv will be the growth term
      if(g0>=0)
        Conv = g0;

      if(fabs(Conv)>0)
       {
        Pe_K = hE*Conv/(2.*D_L);

 // 	beta based on Lutz book
 //         if(Pe_K>1.)
 //           beta = delta0 * hE/Conv;
 //         else
 //           beta = delta1 *hE*hE/D_L ;
 // 	
        // beta based on 1d Green's formula
        beta = hE*(1./tanh(Pe_K) - 1./Pe_K)/(2.*fabs(Conv));
       }
      else
      { beta = 0.0; }

      //  cout<< " Pe_K  " << Pe_K  <<" beta  " << beta  <<endl; 
      //  cout<< " c  " << c0 << " g " <<g0 << " c1 " << c1 <<endl;   
      Mult = Weights[j]*AbsDetjk[j];
      orgD0 = origvaluesD0[j];
      orgD1 = origvaluesD1[j];
      orgD2 = origvaluesD2[j];

      len +=Mult;
      for(k=0;k<N_BaseFunct;k++)
       {
        test0  = orgD0[k];
        test1  = orgD1[k];

	      bgradv = Conv*test1;
 // 	if(k==0)
 // 	  cout<< " bgradv  " << bgradv  <<" bgradv*beta  " << bgradv*beta  <<endl; 
        LocRhs[k] += Mult*rhsval*(test0 + beta*bgradv);

        for(l=0;l<N_BaseFunct;l++)
         {
          ansatz0  = orgD0[l];
          ansatz1  = orgD1[l];
          ansatz2  = orgD2[l];
	  
          val  = Conv*ansatz1*test0;// convective form
          // val  = -Conv*ansatz0*test1; // conservative form                 
          val += c0*ansatz1*test1; // difusive term
          val += c1*ansatz0*test0;
          //  cout<< "A ( " << k <<" , " <<l << " ) " <<Mult*val <<endl;
          LocMatrixA[k*N_BaseFunct + l] += (Mult*val);

          val  = (-c0*ansatz2 + Conv*ansatz1)*beta*bgradv;
          LocMatrixK[k*N_BaseFunct + l] += (Mult*val);
 //            cout<< "S ( " << k <<" , " <<l << " ) " <<bgradv  <<endl;
          val  = ansatz0*beta*bgradv;
          LocMatrixS[k*N_BaseFunct + l] += (Mult*val);

        }
       }
     }
 //      cout<<endl;
 //  exit(0);
    // add to global matrices
    for(j=0;j<N_BaseFunct;j++)
     {
      TestDOF = DOF[j];
      rhs[TestDOF] += LocRhs[j]; 
      
      begin = RowPtr[TestDOF];
      end = RowPtr[TestDOF+1];
      for(k=begin;k<end;k++)
       {
       for(l=0;l<N_BaseFunct;l++)
        {
         if(KCol[k] == DOF[l])
          {
           ValuesA[k] +=LocMatrixA[j*N_BaseFunct + l];
           ValuesK[k] +=LocMatrixK[j*N_BaseFunct + l];
           ValuesS[k] +=LocMatrixS[j*N_BaseFunct + l];
           break;
          }
        } // for(m=0;m<N_BaseFunct_low
      } // for(n=begin;n<end;n++)
     } // for(l=0;l<N_BaseFunct_low
   }// for(i=0; i<N_Cells_Internal; i++)
   
 // cout<< " len " << len << endl;
 // exit(0);

 //  // print matrix
 //    for(j=0;j<N_Dof;j++)
 //     {
 //      begin = RowPtr[j];
 //      end = RowPtr[j+1];
 //      for(k=begin;k<end;k++)
 //       {
 //  //       cout << "A(" << j << ", "<< KCol[k] << ") = " << ValuesA[k] <<endl;
 //        cout << "K(" << j << ", "<< KCol[k] << ") = " << ValuesK[k] <<endl;
 //  //       cout << "S(" << j << ", "<< KCol[k] << ") = " << ValuesS[k] <<endl;
 //       }
      
 //  //       cout << "f: " << rhs[j ] <<endl;
 //      cout<<endl;
 //     }
 //  exit(0);
  for(i=0; i<MaxN_QuadPoints_1D; i++)
   {
    delete [] aux[i];
    delete [] coeff[i];
   }
   
  delete [] aux;
  delete [] coeff;
  
  for(i=0; i<N_Coord+1; i++)  
   delete [] Coords[i]; 
  
 delete [] Coords;    
  
} // void TSystemADI1D::AssembleARhs(d

// /**  Assembling A matrix */
// /** mass mat is same for all quad points in this cell */
void TSystemADI1D::AssembleARhs_DG_Advect(double *X, double *NucGrowth, CoeffFctND *Bilinear)
{
 int i, j, k, l, N_Cells_Internal, N_BaseFunct;
 int N_Points, N_Sets=1, *GlobalNumbers, *BeginIndex, *DOF, *NeibDOF;
 int TestDOF, begin, end, *RowPtr, *KCol; 
 int m, N_Joints, Neigh_i;
 
 double Neigh_rec_detjk, Neigh_N_BaseFunct; 
 double *Weights, *zeta, *Intl_L, AbsDetjk[MaxN_QuadPoints_1D];
 double LocMatrixA11[MaxN_BaseFunctions1D*MaxN_BaseFunctions1D];
 double LocMatrixA12[MaxN_BaseFunctions1D*MaxN_BaseFunctions1D];
 double LocMatrixA21[MaxN_BaseFunctions1D*MaxN_BaseFunctions1D];
 double LocMatrixA22[MaxN_BaseFunctions1D*MaxN_BaseFunctions1D];  
 double LocRhs[MaxN_BaseFunctions1D];
 double BdValues[3];
 double **origvaluesD0, **origvaluesD1, Mult;
 double *orgD0, *orgD1, test0, test1, ansatz0, ansatz1, *ValuesA;
 double c0, c1, g0, rhsval, val, len=0.;
 double **aux, **coeff, *Coeff, **Coords;
 double JointValues[MaxN_BaseFunctions1D], Neigh_JointValues[MaxN_BaseFunctions1D];
 double JointValuesX[MaxN_BaseFunctions1D], Neigh_JointValuesX[MaxN_BaseFunctions1D];
 double Epsilon, Sigma0, Sigma1, h_max, h1, h2;
 double Conv = NucGrowth[0];
//  double B_Nuc = NucGrowth[1];
 bool Needs2ndDer[1], UpdateEdge=FALSE;

 Epsilon = TDatabase::ParamDB->DG_P0;
 Sigma0 = TDatabase::ParamDB->DG_P1;
//  Sigma1 = TDatabase::ParamDB->DG_P2;

 TCollection *Coll;
 TBaseCell *Cell, *Neigh=NULL;
 FE1D FEId, Neigh_FEId;
 TFE1D *Element, *Neigh_Element;
 TBaseFunct1D *bf, *Neigh_bf;
 QuadFormula1D LineQuadFormula;
 TQuadFormula1D *qf1;
 TRefTrans1D *F_K;
 BaseFunct1D BaseFunct_ID, BaseFunct[1];

  Coll = FESpace1D_Intl->GetCollection();
  GlobalNumbers = FESpace1D_Intl->GetGlobalNumbers();
  BeginIndex = FESpace1D_Intl->GetBeginIndex();

  RowPtr = A_Intl->GetRowPtr();
  KCol = A_Intl->GetKCol();
  ValuesA = A_Intl->GetEntries(); 

  N_Cells_Internal = Coll->GetN_Cells();
  Needs2ndDer[0] = FALSE;

  aux = new double *[MaxN_QuadPoints_1D];
  coeff = new double *[MaxN_QuadPoints_1D];

  Coords = new double* [N_Coord+1];
  
  for(i=0; i<N_Coord+1; i++)  
   Coords[i] = new double [MaxN_QuadPoints_1D];
    
  for(i=0; i<MaxN_QuadPoints_1D; i++)
  {
   for(int ii=0; ii<N_Coord; ++ii )
   Coords[ii][i] = X[ii];

   aux[i] = new double[6];
   coeff[i] = new double[6];
  }
   
  Intl_L = Coords[N_Coord];

  // associate each cell with her number in the collection
  for(i=0;i<N_Cells_Internal;i++)
  {
    Cell = Coll->GetCell(i);
    Cell->SetClipBoard(i);
  }

  for(i=0; i<N_Cells_Internal; i++)
  {
    UpdateEdge = FALSE;
    Cell = Coll->GetCell(i);
    FEId = FESpace1D_Intl->GetFE1D(i, Cell);
    Element = TFEDatabase2D::GetFE1D(FEId);
    bf = Element->GetBaseFunct1D();
    N_BaseFunct = Element->GetN_DOF();
    BaseFunct_ID = Element->GetBaseFunct1D_ID();
    h1 = Cell->GetMeasure();

    //  cout << " BaseFunct_ID " << BaseFunct_ID << endl;
    l = bf->GetPolynomialDegree();
    LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(20*l);
    qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
    qf1->GetFormulaData(N_Points, Weights, zeta);

    F_K = TFEDatabase2D::GetRefTrans1D(LineAffin);
    ((TLineAffin *)F_K)->SetCell(Cell);
    ((TLineAffin *)F_K)->GetOrigFromRef(N_Points, zeta, Intl_L, AbsDetjk);
    // rec_detjk = ((TLineAffin *)F_K)->Getrec_detjk();

    Coords[4][1] = i;
    Bilinear(N_Points, N_Coord+1, Coords, aux, coeff); 
    
    BaseFunct[0] = BaseFunct_ID;
    ((TLineAffin *)F_K)->GetOrigValues(N_Sets, BaseFunct, N_Points, zeta,  LineQuadFormula,  Needs2ndDer);

    origvaluesD0=TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D0);
    origvaluesD1=TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D1);

    memset(LocMatrixA11, 0, N_BaseFunct*N_BaseFunct*SizeOfDouble);
    memset(LocMatrixA12, 0, N_BaseFunct*N_BaseFunct*SizeOfDouble);
    memset(LocMatrixA21, 0, N_BaseFunct*N_BaseFunct*SizeOfDouble);
    memset(LocMatrixA22, 0, N_BaseFunct*N_BaseFunct*SizeOfDouble);
    memset(LocRhs, 0, N_BaseFunct*SizeOfDouble);

    DOF = GlobalNumbers + BeginIndex[i];

    for(j=0;j<N_Points;j++)
     {
      Coeff = coeff[j];
      c0 = Coeff[0]; // diffusion
      g0 = Coeff[1]; // convection in z direction
      // c1 = Coeff[2]; // reaction term
      rhsval = Coeff[3]; //rhs

      // in case of PBE, it is the growth term from concentration
      // in PBE example file define Coeff[1] as < 0, so Conv will be the growth term
      if(g0>=0)
        Conv = g0;

      //  cout<< " diffusion  " << c0 << " Conv " <<Conv << " rhsval " << rhsval <<endl;   
      Mult = Weights[j]*AbsDetjk[j];
      orgD0 = origvaluesD0[j];
      orgD1 = origvaluesD1[j];

      len +=Mult;
      for(k=0;k<N_BaseFunct;k++)
       {
        test0  = orgD0[k];
        test1  = orgD1[k];

        LocRhs[k] += Mult*rhsval*test0;
 
        for(l=0;l<N_BaseFunct;l++)
         {
          ansatz0  = orgD0[l];
          ansatz1  = orgD1[l];
 
          // val = Conv*ansatz1*test0;// dg non convective form
          val  = -Conv*ansatz0*test1; // dg convective form         
          // val += c0*ansatz1*test1; // no difusive term
          // val += c1*ansatz0*test0; // no reaction 
            // if(i==0)
            // cout<< "A ( " << k <<" , " <<l << " ) " <<Mult*val <<endl;
          LocMatrixA11[k*N_BaseFunct + l] += (Mult*val);
        } //for(l=0;l<N_Ba
       } // for(k=0;k<N_BaseFunct;k+
    } // for(j=0;j<N_Points;j++)

     N_Joints=Cell->GetN_Edges();
     for(m=0; m<N_Joints; m++)
      {
       Neigh = (Cell->GetJoint(m))->GetNeighbour(Cell);

       if(m==0 && Neigh) // inner joint, since in 1D we have only 2 joints (really vertices)
       {
        UpdateEdge = TRUE;
        // only first joint (really vertices) will be updated, 
        // other joint will be the first joint of the next cell
        h2 = Neigh->GetMeasure();
        h_max = MAX(h1, h2);

        // Intl_L[0] = TDatabase::ParamDB->REACTOR_P12;
        //  Bilinear(1, Space_X, Space_Y, Intl_L, aux, coeff);
        Bilinear(1, N_Coord+1, Coords, aux, coeff); 
    
        c0 = coeff[0][0]; // diffusion

        //find the current cell basis value at this joint (X_n+)
        bf->GetDerivatives(D0, -1., JointValues);
        // bf->GetDerivatives(D1, -1., JointValuesX);

        //find the neib cell basis value at this joint (X_n-)
        Neigh_i=Neigh->GetClipBoard();
        Neigh_FEId = FESpace1D_Intl->GetFE1D(Neigh_i, Neigh);
        Neigh_Element = TFEDatabase2D::GetFE1D(Neigh_FEId);
        Neigh_bf = Neigh_Element->GetBaseFunct1D();
        Neigh_N_BaseFunct = Neigh_Element->GetN_DOF();

        NeibDOF = GlobalNumbers + BeginIndex[Neigh_i];

        F_K = TFEDatabase2D::GetRefTrans1D(LineAffin);
       ((TLineAffin *)F_K)->SetCell(Neigh);
        Neigh_rec_detjk = ((TLineAffin *)F_K)->Getrec_detjk();

        Neigh_bf->GetDerivatives(D0, 1., Neigh_JointValues);
        Neigh_bf->GetDerivatives(D1, 1., Neigh_JointValuesX);

        for(k=0;k<Neigh_N_BaseFunct;k++)
         Neigh_JointValuesX[k] *=Neigh_rec_detjk;

        //(X_n+, X_n+) (test, ansatz)
        for(k=0;k<N_BaseFunct;k++)
         for(l=0;l<N_BaseFunct;l++)
         {
          // val = (-0.5*Conv -Epsilon*Conv + (Sigma0/h_max) )*JointValues[l]*JointValues[k];
          // LocMatrixA11[k*N_BaseFunct + l] += val;
         }

        // e_n
        for(k=0;k<N_BaseFunct;k++) // own test X_n+
         for(l=0;l<Neigh_N_BaseFunct;l++) // neib ansatz X_n-
          {
          //  val = (-0.5*Conv -Epsilon*Conv - (Sigma0/h_max) )*JointValues[l]*JointValues[k];
          //  LocMatrixA12[k*N_BaseFunct + l] += val;
          }

        // d_n
        for(k=0;k<Neigh_N_BaseFunct;k++) // neib test X_n-
         for(l=0;l<N_BaseFunct;l++) // own ansatz X_n+
         {
          //dg convective form
          // val = (0.5*Conv +Epsilon*Conv - (Sigma0/h_max) )*JointValues[l]*JointValues[k];
          //  LocMatrixA21[k*N_BaseFunct + l] += val;
          }

        //(X_n-, X_n-) (test, ansatz)
        for(k=0;k<Neigh_N_BaseFunct;k++)
         for(l=0;l<Neigh_N_BaseFunct;l++)
         {
          //dg convective form
          // val = (0.5*Conv + Epsilon*Conv + (Sigma0/h_max) )*JointValues[l]*JointValues[k];
          // LocMatrixA22[k*N_BaseFunct  + l] += val;
         }
       }//if(neigh)
       else if(!Neigh) // boundary joint
       {
        if(m==0) // L_min 
         {

          // if(cond_Lmin==DIRICHLET || cond_Lmin==)
          {
            //find the current cell basis value at this joint (X_n+)
           bf->GetDerivatives(D0, -1., JointValues);
           for(k=0;k<N_BaseFunct;k++)
            {
             LocRhs[k] += JointValues[k]*NucGrowth[1];// Bnuc
            }
          } //if(cond_Lmin==DIRICHLET)
         } // if(m==0) 
       }

      }

      // add to global matrices
      for(j=0;j<N_BaseFunct;j++)
     {
      TestDOF = DOF[j];
      rhs[TestDOF] += LocRhs[j]; 

      begin = RowPtr[TestDOF];
      end = RowPtr[TestDOF+1];
      for(k=begin;k<end;k++)
       {
        for(l=0;l<N_BaseFunct;l++)
         {
          if(KCol[k] == DOF[l])
           {
            ValuesA[k] +=LocMatrixA11[j*N_BaseFunct + l];
            //  cout << TestDOF << " " << KCol[k] << endl;
            break;
           }
          } // for(m=0;m<N_BaseFunct

       // add edge integrals
      if(UpdateEdge)
       {
        for(l=0;l<Neigh_N_BaseFunct;l++)
         {
         if(KCol[k] == NeibDOF[l])
          {
            //  cout << TestDOF << " " << KCol[k] << endl;
           ValuesA[k] +=LocMatrixA12[j*N_BaseFunct + l]; 
           break;
          } 
        } // for(l=0;l<Neigh_N_BaseFunct;l++)
       }

      } // for(n=begin;n<end;n++)
     } //  for(j=0;j<N_BaseFunct;j++)

      // add edge integrals
     if(UpdateEdge)
     {
     for(j=0;j<Neigh_N_BaseFunct;j++)
     {
      TestDOF = NeibDOF[j];

      begin = RowPtr[TestDOF];
      end = RowPtr[TestDOF+1];
      for(k=begin;k<end;k++)
       {
        for(l=0;l<N_BaseFunct;l++)
         {
          if(KCol[k] == DOF[l])
           {
            ValuesA[k] +=LocMatrixA21[j*N_BaseFunct + l];
            //  cout << TestDOF << " " << KCol[k] << endl;
            break;
           }
          } //  for(l=0;l<N_BaseFunct;l++)

         for(l=0;l<Neigh_N_BaseFunct;l++)
         {
          if(KCol[k] == NeibDOF[l])
           {
            ValuesA[k] +=LocMatrixA22[j*N_BaseFunct + l];
            //  cout << TestDOF << " " << KCol[k] <<  " " <<  LocMatrixA22[j*N_BaseFunct + l] << endl;
            break;
           }
          } //  for(l=0;l<N_BaseFunct;l++)
        } // for(k=begin;k<end;k++)
     } // for(j=0;j<Neigh_N_BaseFunct;j++)
     } // if(UpdateEdge)
    

   }// for(i=0; i<N_Cells_Internal; i++)

//  cout<< " len " << len << endl;
 //  exit(0);
// //  print matrix
//     for(j=0;j<N_Dof;j++)
//      {
//       begin = RowPtr[j];
//       end = RowPtr[j+1];
//       // for(k=begin;k<end;k++)
//       //  {
//       //   cout << "A(" << j << ", "<< KCol[k] << ") = " << ValuesA[k] <<endl;
//       //  }
//       cout << "f: " << rhs[j ] <<endl;
//       cout<<endl;
//       }
//   exit(0);
  for(i=0; i<MaxN_QuadPoints_1D; i++)
   {
    delete [] aux[i];
    delete [] coeff[i];
   }

  delete [] aux;
  delete [] coeff;
} // void TSystemADI1D::  



// /**  Assembling A matrix */
// /** mass mat is same for all quad points in this cell */
void TSystemADI1D::AssembleARhs_DG(double *X, double Conv, double *Out, CoeffFctND *Bilinear)
{
 int i, j, k, l, N_Cells_Internal, N_BaseFunct;
 int N_Points, N_Sets=1, *GlobalNumbers, *BeginIndex, *DOF, *NeibDOF;
 int TestDOF, begin, end, *RowPtr, *KCol; 
 int m, N_Joints, Neigh_i;
 
 double rec_detjk, Neigh_rec_detjk, Neigh_N_BaseFunct; 
 double *Weights, *zeta, *Intl_L, AbsDetjk[MaxN_QuadPoints_1D];
 double LocMatrixA11[MaxN_BaseFunctions1D*MaxN_BaseFunctions1D];
 double LocMatrixA12[MaxN_BaseFunctions1D*MaxN_BaseFunctions1D];
 double LocMatrixA21[MaxN_BaseFunctions1D*MaxN_BaseFunctions1D];
 double LocMatrixA22[MaxN_BaseFunctions1D*MaxN_BaseFunctions1D];  
 double LocRhs[MaxN_BaseFunctions1D];
 double BdValues[3];
 double **origvaluesD0, **origvaluesD1, Mult;
 double *orgD0, *orgD1, test0, test1, ansatz0, ansatz1, *ValuesA;
 double c0, c1, g0, rhsval, val, len=0.;
 double **aux, **coeff, *Coeff, **Coords;
 double JointValues[MaxN_BaseFunctions1D], Neigh_JointValues[MaxN_BaseFunctions1D];
 double JointValuesX[MaxN_BaseFunctions1D], Neigh_JointValuesX[MaxN_BaseFunctions1D];
 double Epsilon, Sigma0, Sigma1, h_max, h1, h2;

 bool Needs2ndDer[1], UpdateEdge=FALSE;

 Epsilon = TDatabase::ParamDB->DG_P0;
 Sigma0 = TDatabase::ParamDB->DG_P1;
 Sigma1 = TDatabase::ParamDB->DG_P2;

 TCollection *Coll;
 TBaseCell *Cell, *Neigh=NULL;
 FE1D FEId, Neigh_FEId;
 TFE1D *Element, *Neigh_Element;
 TBaseFunct1D *bf, *Neigh_bf;
 QuadFormula1D LineQuadFormula;
 TQuadFormula1D *qf1;
 TRefTrans1D *F_K;
 BaseFunct1D BaseFunct_ID, BaseFunct[1];

  Coll = FESpace1D_Intl->GetCollection();
  GlobalNumbers = FESpace1D_Intl->GetGlobalNumbers();
  BeginIndex = FESpace1D_Intl->GetBeginIndex();

  RowPtr = A_Intl->GetRowPtr();
  KCol = A_Intl->GetKCol();
  ValuesA = A_Intl->GetEntries(); 

  N_Cells_Internal = Coll->GetN_Cells();
  Needs2ndDer[0] = FALSE;

  aux = new double *[MaxN_QuadPoints_1D];
  coeff = new double *[MaxN_QuadPoints_1D];

  Coords = new double* [N_Coord+1];
  
  for(i=0; i<N_Coord+1; i++)  
   Coords[i] = new double [MaxN_QuadPoints_1D];
    
  for(i=0; i<MaxN_QuadPoints_1D; i++)
  {
   for(int ii=0; ii<N_Coord; ++ii )
   Coords[ii][i] = X[ii];

   aux[i] = new double[6];
   coeff[i] = new double[6];
  }
   
  Intl_L = Coords[N_Coord];

  // associate each cell with her number in the collection
  for(i=0;i<N_Cells_Internal;i++)
  {
    Cell = Coll->GetCell(i);
    Cell->SetClipBoard(i);
  }

  for(i=0; i<N_Cells_Internal; i++)
  {
    UpdateEdge = FALSE;
    Cell = Coll->GetCell(i);
    FEId = FESpace1D_Intl->GetFE1D(i, Cell);
    Element = TFEDatabase2D::GetFE1D(FEId);
    bf = Element->GetBaseFunct1D();
    N_BaseFunct = Element->GetN_DOF();
    BaseFunct_ID = Element->GetBaseFunct1D_ID();
    h1 = Cell->GetMeasure();

 //     cout << " BaseFunct_ID " << BaseFunct_ID << endl;
    l = bf->GetPolynomialDegree();
    LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(20*l);
    qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
    qf1->GetFormulaData(N_Points, Weights, zeta);

    F_K = TFEDatabase2D::GetRefTrans1D(LineAffin);
    ((TLineAffin *)F_K)->SetCell(Cell);
    ((TLineAffin *)F_K)->GetOrigFromRef(N_Points, zeta, Intl_L, AbsDetjk);
    rec_detjk = ((TLineAffin *)F_K)->Getrec_detjk();

 //     Bilinear(N_Points, Space_X, Space_Y, Intl_L, aux, coeff);
    Coords[4][1] = i;
    Bilinear(N_Points, N_Coord+1, Coords, aux, coeff); 
    
    BaseFunct[0] = BaseFunct_ID;
    ((TLineAffin *)F_K)->GetOrigValues(N_Sets, BaseFunct, N_Points, zeta,  LineQuadFormula,  Needs2ndDer);

    origvaluesD0=TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D0);
    origvaluesD1=TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D1);

    memset(LocMatrixA11, 0, N_BaseFunct*N_BaseFunct*SizeOfDouble);
    memset(LocMatrixA12, 0, N_BaseFunct*N_BaseFunct*SizeOfDouble);
    memset(LocMatrixA21, 0, N_BaseFunct*N_BaseFunct*SizeOfDouble);
    memset(LocMatrixA22, 0, N_BaseFunct*N_BaseFunct*SizeOfDouble);
    memset(LocRhs, 0, N_BaseFunct*SizeOfDouble);

    DOF = GlobalNumbers + BeginIndex[i];

    for(j=0;j<N_Points;j++)
     {
      Coeff = coeff[j];
      c0 = Coeff[0]; // diffusion
      g0 = Coeff[1]; // convection in z direction
      c1 = Coeff[2]; // reaction term
      rhsval = Coeff[3]; //rhs

      // in case of PBE, it is the growth term from concentration
      // in PBE example file define Coeff[1] as < 0, so Conv will be the growth term
      if(g0>=0)
        Conv = g0;

      //  cout<< " diffusion  " << c0 << " Conv " <<Conv << " rhsval " << rhsval <<endl;   
      Mult = Weights[j]*AbsDetjk[j];
      orgD0 = origvaluesD0[j];
      orgD1 = origvaluesD1[j];

      len +=Mult;
      for(k=0;k<N_BaseFunct;k++)
       {
        test0  = orgD0[k];
        test1  = orgD1[k];

        LocRhs[k] += Mult*rhsval*test0;

        for(l=0;l<N_BaseFunct;l++)
         {
          ansatz0  = orgD0[l];
          ansatz1  = orgD1[l];
 
          val  = Conv*ansatz1*test0;// dg non convective form
          // val  = -Conv*ansatz0*test1; // dg convective form         
          val += c0*ansatz1*test1; // difusive term
          val += c1*ansatz0*test0;
            // if(i==0)
            // cout<< "A ( " << k <<" , " <<l << " ) " <<Mult*val <<endl;
          LocMatrixA11[k*N_BaseFunct + l] += (Mult*val);
        }
       }
     }

     N_Joints=Cell->GetN_Edges();
     for(m=0; m<N_Joints; m++)
      {
       Neigh = (Cell->GetJoint(m))->GetNeighbour(Cell);

       if(m==0 && Neigh) // inner left joint, since in 1D we have only 2 joints (really vertices)
       {
        UpdateEdge = TRUE;
        // only first joint (really vertices) will be updated, 
        // other joint will be the first joint of the next cell
        h2 = Neigh->GetMeasure();
        h_max = MAX(h1, h2);

        // Intl_L[0] = TDatabase::ParamDB->REACTOR_P12;
        // Bilinear(1, N_Coord+1, Coords, aux, coeff); 
    
        // c0 = coeff[0][0]; // no diffusion

        //find the current cell basis value at this joint (X_n+)
        bf->GetDerivatives(D0, -1., JointValues);
        // bf->GetDerivatives(D1, -1., JointValuesX);

        // for(k=0;k<N_BaseFunct;k++)
        //  JointValuesX[k] *=rec_detjk;

        //find the neib cell basis value at this joint (X_n-)
        Neigh_i=Neigh->GetClipBoard();
        Neigh_FEId = FESpace1D_Intl->GetFE1D(Neigh_i, Neigh);
        Neigh_Element = TFEDatabase2D::GetFE1D(Neigh_FEId);
        Neigh_bf = Neigh_Element->GetBaseFunct1D();
        Neigh_N_BaseFunct = Neigh_Element->GetN_DOF();

        NeibDOF = GlobalNumbers + BeginIndex[Neigh_i];

        F_K = TFEDatabase2D::GetRefTrans1D(LineAffin);
       ((TLineAffin *)F_K)->SetCell(Neigh);
        Neigh_rec_detjk = ((TLineAffin *)F_K)->Getrec_detjk();

        Neigh_bf->GetDerivatives(D0, 1., Neigh_JointValues);
        // Neigh_bf->GetDerivatives(D1, 1., Neigh_JointValuesX);

        // for(k=0;k<Neigh_N_BaseFunct;k++)
        //  Neigh_JointValuesX[k] *=Neigh_rec_detjk;

        //(X_n+, X_n+) (test, ansatz)
        for(k=0;k<N_BaseFunct;k++)
         for(l=0;l<N_BaseFunct;l++)
         {
          val = -0.5*Conv*JointValues[l]*JointValuesX[k] - 0.5*c0*Epsilon*JointValues[l]*JointValuesX[k]
                + (Sigma0/h_max)*JointValues[l]*JointValues[k] + (Sigma1/h_max)*JointValuesX[l]*JointValuesX[k];
   
           //dg convective form
          // val = Conv*JointValues[l]*JointValues[k];

         LocMatrixA11[k*N_BaseFunct + l] += val;
 //          cout<<  " val  " << h_max*val  <<endl;
         }

        // e_n
        for(k=0;k<N_BaseFunct;k++) // own test X_n+
         for(l=0;l<Neigh_N_BaseFunct;l++) // neib ansatz X_n-
          {
           val = 0.5*c0*Neigh_JointValuesX[l]*JointValues[k] + 0.5*c0*Epsilon*Neigh_JointValues[l]*JointValuesX[k]
                 - (Sigma0/h_max)*Neigh_JointValues[l]*JointValues[k] - (Sigma1/h_max)*Neigh_JointValuesX[l]*JointValuesX[k];

           //dg convective form
          // val += 0.5*Conv*Neigh_JointValues[l]*JointValues[k];

           LocMatrixA21[k*N_BaseFunct + l] += val;
 //            cout<<  " val  " << h_max*val  <<endl;
          }

        // d_n
        for(k=0;k<Neigh_N_BaseFunct;k++) // neib test X_n-
         for(l=0;l<N_BaseFunct;l++) // own ansatz X_n+
         {
          val = -0.5*c0*JointValuesX[l]*Neigh_JointValues[k] - 0.5*c0*Epsilon*JointValues[l]*Neigh_JointValuesX[k]
                 - (Sigma0/h_max)*JointValues[l]*Neigh_JointValues[k] - (Sigma1/h_max)*JointValuesX[l]*Neigh_JointValuesX[k];
          //dg convective form
          // val -= 0.5*Conv*JointValues[l]*Neigh_JointValues[k];

           LocMatrixA12[k*N_BaseFunct + l] += val;
 //            cout<<  " val  " << h_max*val  <<endl;
          }

        //(X_n-, X_n-) (test, ansatz)
        for(k=0;k<Neigh_N_BaseFunct;k++)
         for(l=0;l<Neigh_N_BaseFunct;l++)
         {
          val = -0.5*c0*Neigh_JointValuesX[l]*Neigh_JointValues[k] + 0.5*c0*Epsilon*Neigh_JointValues[l]*Neigh_JointValuesX[k]
                + (Sigma0/h_max)*Neigh_JointValues[l]*Neigh_JointValues[k] + (Sigma1/h_max)*Neigh_JointValuesX[l]*Neigh_JointValuesX[k];
           //dg convective form
          // val -= 0.5*Conv*Neigh_JointValues[l]*Neigh_JointValues[k];

          LocMatrixA22[k*N_BaseFunct  + l] += val;
 //          cout<<  " val  " << h_max*val  <<endl;
         }
       }//if(neigh)
       else if(!Neigh) // boundary joint
       {
        if(m==0) // L_min 
         {

         if(cond_Lmin==DIRICHLET)
          {
            //find the current cell basis value at this joint (X_n+)
           bf->GetDerivatives(D0, -1., JointValues);
           bf->GetDerivatives(D1, -1., JointValuesX);

          //  Intl_L[0] = TDatabase::ParamDB->REACTOR_P12;
 //            Bilinear(1, Space_X, Space_Y, Intl_L, aux, coeff);
           Bilinear(1, N_Coord+1, Coords, aux, coeff);
   
           c0 = coeff[0][0]; // diffusion

 //            BDValue_LMin(Intl_L[0], 0, BdValues);

           BdValues[0]=Out[1]; // Bnuc
 
           for(k=0;k<N_BaseFunct;k++)
            JointValuesX[k] *=rec_detjk;

           for(k=0;k<N_BaseFunct;k++)
            {
             val =  -Epsilon*c0*JointValuesX[k]*BdValues[0];
             val +=  (Sigma0/h1)*JointValues[k]*BdValues[0];

             LocRhs[k] += val;
                // cout<<  " val    rhsval " << val  <<endl;
             for(l=0;l<N_BaseFunct;l++)
              {
               val = c0*JointValuesX[l]*JointValues[k] - c0*Epsilon*JointValues[l]*JointValuesX[k]
                    + (Sigma0/h1)*JointValues[l]*JointValues[k] + (Sigma1/h1)*JointValuesX[l]*JointValuesX[k];

               LocMatrixA11[k*N_BaseFunct + l] += val;
 //                 cout<<  " val L_min  " << val  <<endl;
              }
             }
          } //if(cond_Lmin==DIRICHLET)
         } // if(m==0) 
        else // L_Max
	       {
          // if(cond_Lmax==DIRICHLET || cond_Lmax==NEUMANN)
          if(cond_Lmax==DIRICHLET)          
           {
             //find the current cell basis value at this joint (X_n+)
            bf->GetDerivatives(D0, 1., JointValues);
            bf->GetDerivatives(D1, 1., JointValuesX);

            // Intl_L[0] = TDatabase::ParamDB->REACTOR_P13;
 //             Bilinear(1, Space_X, Space_Y, Intl_L, aux, coeff);
            Bilinear(1, N_Coord+1, Coords, aux, coeff);
   
            c0 = coeff[0][0]; // diffusion

            // BDVal_LMax(Intl_L[0], 0, BdValues);
             BdValues[0] = Out[2];
             BdValues[1] = Out[2];


            for(k=0;k<N_BaseFunct;k++)
             JointValuesX[k] *=rec_detjk;

            for(k=0;k<N_BaseFunct;k++)
             {
              val =  Epsilon*c0*JointValuesX[k]*BdValues[0];
              val +=  (Sigma1/h1)*JointValuesX[k]*BdValues[1];

              LocRhs[k] += val;
 //                 cout<<  " val    rhsval " << val  <<endl;
             for(l=0;l<N_BaseFunct;l++)
              {
               val = -c0*JointValuesX[l]*JointValues[k] + c0*Epsilon*JointValues[l]*JointValuesX[k]
                     + (Sigma0/h1)*JointValues[l]*JointValues[k] + (Sigma1/h1)*JointValuesX[l]*JointValuesX[k];

               LocMatrixA11[k*N_BaseFunct + l] += val;
 //                cout<<  " val L_max " << val  <<endl;
              }

            }
          }
         }
       }

      }


    // add to global matrices
    for(j=0;j<N_BaseFunct;j++)
     {
      TestDOF = DOF[j];
      rhs[TestDOF] += LocRhs[j]; 

      begin = RowPtr[TestDOF];
      end = RowPtr[TestDOF+1];
      for(k=begin;k<end;k++)
       {
        for(l=0;l<N_BaseFunct;l++)
         {
          if(KCol[k] == DOF[l])
           {
            ValuesA[k] +=LocMatrixA11[j*N_BaseFunct + l];
 //             cout << TestDOF << " " << KCol[k] << endl;
            break;
           }
          } // for(m=0;m<N_BaseFunct

       // add edge integrals
      if(UpdateEdge)
       {
        for(l=0;l<Neigh_N_BaseFunct;l++)
         {
         if(KCol[k] == NeibDOF[l])
          {
 //             cout << TestDOF << " " << KCol[k] << endl;
           ValuesA[k] +=LocMatrixA21[j*N_BaseFunct + l]; 
           break;
          } 
        } // for(l=0;l<Neigh_N_BaseFunct;l++)
       }

      } // for(n=begin;n<end;n++)
     } //  for(j=0;j<N_BaseFunct;j++)

   // add edge integrals
   if(UpdateEdge)
   {
    for(j=0;j<Neigh_N_BaseFunct;j++)
     {
      TestDOF = NeibDOF[j];

      begin = RowPtr[TestDOF];
      end = RowPtr[TestDOF+1];
      for(k=begin;k<end;k++)
       {
        for(l=0;l<N_BaseFunct;l++)
         {
          if(KCol[k] == DOF[l])
           {
            ValuesA[k] +=LocMatrixA12[j*N_BaseFunct + l];
 //             cout << TestDOF << " " << KCol[k] << endl;
            break;
           }
          } //  for(l=0;l<N_BaseFunct;l++)

        for(l=0;l<Neigh_N_BaseFunct;l++)
         {
          if(KCol[k] == NeibDOF[l])
           {
            ValuesA[k] +=LocMatrixA22[j*N_BaseFunct + l];
 //             cout << TestDOF << " " << KCol[k] <<  " " <<  LocMatrixA22[j*N_BaseFunct + l] << endl;
            break;
           }
          } //  for(l=0;l<N_BaseFunct;l++)
        } // for(k=begin;k<end;k++)
     } // for(j=0;j<Neigh_N_BaseFunct;j++)
    } // if(UpdateEdge)
    

   }// for(i=0; i<N_Cells_Internal; i++)

 // cout<< " len " << len << endl;
 //  exit(0);
//  print matrix
//     for(j=0;j<N_Dof;j++)
//      {
//  //      begin = RowPtr[j];
//  //      end = RowPtr[j+1];
//  //      for(k=begin;k<end;k++)
//  //       {
//  //        cout << "A(" << j << ", "<< KCol[k] << ") = " << ValuesA[k] <<endl;
//  //       }
 
//         cout << "f: " << rhs[j ] <<endl;
//       cout<<endl;
//       }
//   exit(0);
  for(i=0; i<MaxN_QuadPoints_1D; i++)
   {
    delete [] aux[i];
    delete [] coeff[i];
   }

  delete [] aux;
  delete [] coeff;
} // void TSystemADI1D::AssembleARhs(d


TSystemADI1D::~TSystemADI1D()
{
//  delete [] Nodal_sol;
//  delete [] oldsol;
 delete [] rhs;
//  delete [] sol;
 delete [] B;
//  delete [] defect;

}


