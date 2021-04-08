/** =======================================================================
* @class     TSystemADI1D_3L
* @brief     stores the information of ADI system 
* @author    Sashikumaar Ganesan
* @date      15.04.2020
* @History 
* ======================================================================= */
#include <SystemADI1D_3L.h>
#include <Database.h>
#include <FEDatabase2D.h>
#include <DirectSolver.h>
#include <SquareStructure1D.h>
#include <SquareMatrix1D.h>
#include <FEFunction1D.h>
#include <FEFunction2D.h>
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
TSystemADI1D_3L::TSystemADI1D_3L(TSystemADI1D **systemADI, int *l_NVertices, double *l_StartEnd, BoundCond1D **boundConLminLMax, 
                    DoubleFunctND **growthAndB_Nuc, TFEFunction2D *spatialFeFunction):TSystemADI1D()
{        

 SystemADI = systemADI;
 SpatialFeFunction = spatialFeFunction;

 for(int i=0; i<3;++i)
  SystemADI[i] = new TSystemADI1D(l_NVertices[i], l_StartEnd[2*i], l_StartEnd[2*i+1], boundConLminLMax[i], growthAndB_Nuc[i]);                               
}

// int_Lv_La at l_d=0
void TSystemADI1D_3L::Int_B_Nuc(double *IntValue, double *B_NucValue, DoubleFunctND *GrowthAndB_Nuc, double *SuscepPopRatio)
{
  TFESpace1D *FeSpace_L1 = SystemADI[1]->GetFeSpace1D();
  TFESpace1D *FeSpace_L2 = SystemADI[2]->GetFeSpace1D();  
  TCollection *L1Coll = FeSpace_L1->GetCollection();
  int i, i1, i2, l1, l2, N_L1Cells = L1Coll->GetN_Cells();

  TCollection *L2Coll = FeSpace_L2->GetCollection();
  int N_L2Cells = L2Coll->GetN_Cells();

  TBaseCell *cell_L2, *cell_L1;
  QuadFormula1D LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(3);//3rd order quadrature
  TQuadFormula1D *qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
  int N_QuadPts;
  double *Weights, *zeta, *IntValue_Xpos;
  double Mult0;
  qf1->GetFormulaData(N_QuadPts, Weights, zeta);
  TRefTrans1D *F_K2, *F_K1;
  double mult1, mult2, Coords[6], IVal[3];

  F_K2 = TFEDatabase2D::GetRefTrans1D(LineAffin);
  F_K1 = TFEDatabase2D::GetRefTrans1D(LineAffin);

  memset(B_NucValue, 0, N_Xpos*SizeOfDouble);

   for(i2=0; i2<N_L2Cells; ++i2) // la
   {
    cell_L2 = L2Coll->GetCell(i2);
    ((TLineAffin *)F_K2)->SetCell(cell_L2);
    ((TLineAffin *)F_K2)->GetOrigFromRef(N_QuadPts, zeta, X2, AbsDetjk_L2);
     
    for(i1=0; i1<N_L1Cells; ++i1) //lv
     {
      cell_L1 = L1Coll->GetCell(i1);
      ((TLineAffin *)F_K1)->SetCell(cell_L1);
      ((TLineAffin *)F_K1)->GetOrigFromRef(N_QuadPts, zeta, X1, AbsDetjk_L1);

      memset(IntValue_XposL0quadL1quad, 0, N_Xpos*N_QuadPts*SizeOfDouble);
      for(l2=0; l2<N_QuadPts; ++l2)
       {
        mult2 = Weights[l2]*AbsDetjk_L2[l2];
        Coords[2] =  X2[l2]; //la
        IntValue_Xpos = IntValue_XposL0quadL1quad + l2*N_Xpos;

        for(l1=0; l1<N_QuadPts; ++l1)
         {
          Coords[3] =  X1[l2]; //lv
          mult1 = Weights[l1]*AbsDetjk_L1[l1];

          for(i=0; i<N_Xpos; i++)
           {
            Coords[0] = Xpos[2*i];
            Coords[1] = Xpos[2*i+1];  
            Coords[4] = IntValue[i];  
            Coords[5] = SuscepPopRatio[i];             
            // IVal[0] = SuscepPopRatio[i];

            GrowthAndB_Nuc(6, Coords, IVal);
            IntValue_Xpos[i] += IVal[1]*mult1;

            //  cout <<"Need to integrate over all l_d!!! QVal " << QVal[0] <<endl;
            //  exit(0);

            // IntValue_Xpos[i] += 1.*mult1;            
           }// i=0; i<N_Xpos  
         } //  for(l1=0; l1

        for(i=0; i<N_Xpos; i++)
         {
          B_NucValue[i] += IntValue_Xpos[i]*mult2;
         }// i
       } // l2
     } //i1
   }//i2

  // for(i=0; i<N_Xpos; i++)
  //    cout << "B_nuc " <<B_NucValue[i] << endl;
  // exit(0);
}


void TSystemADI1D_3L::Init(int n_Xpos, double *xpos, int *n_LDof, int *n_LnodalPos, double **lnodalPos, double *sol_XdoLdof)
{
 // NOTE that these variables cannot be accesse by SystemADI[] and viceversa, 
 // since SystemADI[] are independed instances not this pointer(s)
 N_Xpos = n_Xpos;
 Xpos = xpos; 
 N_LDof = n_LDof;
 N_ADISystems = 3;
 N_LnodalPos = n_LnodalPos;
 LnodalPos = lnodalPos;
 Sol_XdofLdof = sol_XdoLdof;
 N_L0ADIs = N_Xpos*N_LnodalPos[2]*N_LnodalPos[1];
 N_L1ADIs = N_Xpos*N_LnodalPos[0]*N_LnodalPos[2];
 N_L2ADIs = N_Xpos*N_LnodalPos[0]*N_LnodalPos[1];
 N_L1L0 = N_LnodalPos[0]*N_LnodalPos[1];
 N_L2L1 = N_LnodalPos[2]*N_LnodalPos[1]; 
 N_Coord=5;
 N_Lnodal_All = N_LnodalPos[0]*N_LnodalPos[1]*N_LnodalPos[2];
 N_XposLnodal_All = N_Xpos*N_Lnodal_All;
 IncidentArray = new int[N_XposLnodal_All];
 
 Sol_XposLNnodal = new double[N_XposLnodal_All];
 Sol_XposL2L0L1nodal = new double[N_XposLnodal_All];
 Sol_XposL0L1L2nodal = new double[N_XposLnodal_All];

 Sol_XnodalLnodal[0] = Sol_XposLNnodal;
 Sol_XnodalLnodal[1] = Sol_XposL2L0L1nodal;
 Sol_XnodalLnodal[2] = Sol_XposL0L1L2nodal;

//  if(SystemADI[0]->IsdGdisc())
  Sol_XnodalLDof[0] = new double[N_L0ADIs*N_LDof[0]];

//  if(SystemADI[1]->IsdGdisc())
  Sol_XnodalLDof[1] = new double[N_L1ADIs*N_LDof[1]];
 
//  if(SystemADI[2]->IsdGdisc())
  Sol_XnodalLDof[2] = new double[N_L2ADIs*N_LDof[2]];

 for(int i=0; i<3;++i)
   SystemADI[i]->Init(N_Coord, n_Xpos, xpos, 3, n_LnodalPos, lnodalPos, i, Sol_XnodalLnodal[i]);


  // for Int_L
  Sol_XposL2nodalL0quadL1nodal = new double [N_Xpos*N_LnodalPos[2]*MaxN_QuadPoints_1D*N_LnodalPos[1]];
  IntValue_XposL0quad = new double[N_Xpos*MaxN_QuadPoints_1D];
  NucValue_XposL0quad = new double[N_Xpos*MaxN_QuadPoints_1D];
  Sol_XposL1quadL0quadL2nodal = new double [N_Xpos*MaxN_QuadPoints_1D*MaxN_QuadPoints_1D*N_LnodalPos[2]];
  IntValue_XposL0quadL1quad = new double[N_Xpos*MaxN_QuadPoints_1D*MaxN_QuadPoints_1D];  
  NucValue_XposL0quadL1quad = new double[N_Xpos*MaxN_QuadPoints_1D*MaxN_QuadPoints_1D];  
}


void TSystemADI1D_3L::Solve(int ADI_SolveIdx, CoeffFctND *Bilinear, double *IntL_Values, double *SuscepPopRatio)
{
  double *Sol_ALL, *Sol_L, *Sol_LL, *Sol_LLL;
  int i, j, k, l, m, n;
  double *LLLpos, *LLpos, *Lpos;
  double Coords[6];

  switch(ADI_SolveIdx)
  {
  case 0:
     l = N_LnodalPos[2]; LLLpos = LnodalPos[2];
     m = N_LnodalPos[1]; LLpos = LnodalPos[1];
     n = N_LDof[0]; Lpos = LnodalPos[0];
    //  Disp = N_L0ADIs;
    break;
  
  case 1:
     l = N_LnodalPos[2]; LLLpos = LnodalPos[2];
     m = N_LnodalPos[0]; LLpos = LnodalPos[0];
     n = N_LDof[1]; Lpos = LnodalPos[1];
    //  Disp = N_L1ADIs;     
    break;

  case 2:
     l = N_LnodalPos[0]; LLLpos = LnodalPos[0];
     m = N_LnodalPos[1]; LLpos = LnodalPos[1];
     n = N_LDof[2]; Lpos = LnodalPos[2];
    //  Disp = N_L2ADIs;     
    break;

  default:
     ErrMsg("Check ADIIndex");
    break;
  }

  //sove for all Ldof
  Sol_ALL = Sol_XnodalLDof[ADI_SolveIdx];
  if(ADI_SolveIdx)
  this->Nodal2DOF(ADI_SolveIdx, Sol_ALL); // for ADI_SolveIdx=0, already done in LdofXdofToXnodalLdof

  for(i=0; i<N_Xpos; ++i)
  {
    Sol_LLL = Sol_ALL+i*N_Lnodal_All;
    Coords[0] = Xpos[2*i];
    Coords[1] = Xpos[2*i+1];   
    Coords[4] = IntL_Values[i]; // Int_Omega L , total numberdensity at this point
    Coords[5] = SuscepPopRatio[i];

    for(j=0; j<l; ++j)
     {
      Sol_LL = Sol_LLL+(j*m*n);
      Coords[2] =  LLLpos[j];  

      for(k=0; k<m; ++k)
       { 
        Sol_L = Sol_LL+(k*n);
        Coords[3] =  LLpos[k];  
         
        SystemADI[ADI_SolveIdx]->Solve(6, Coords, Bilinear, Sol_L);
       }
     }
  }

  //needed for XnodalLnodalToLdofXdof
  this->DOF2Nodal(ADI_SolveIdx, Sol_ALL);
}

void TSystemADI1D_3L::Interpolate(TCollection *coll, double *InitialPopulation, DoubleFunctND *Exact)
{
  double Coords[6], *Sol_L0;
  int i, ii, i_x, i_l1, i_l2, GlocalCellNo;

 //  for(ii=0; ii<N_L0ADIs; ii++)
 //    cout << ii << " Xpos: " << (ii/N_LnodalPos[1])/N_LnodalPos[2] << " L2: " << (ii/N_LnodalPos[1])%N_LnodalPos[2] << " L1: " <<  ii%N_LnodalPos[1] <<endl;

 //  for(ii=0; ii<N_L1ADIs; ii++)
 //    cout << ii << " Xpos: " << (ii/N_LnodalPos[0])/N_LnodalPos[2] << " L2: " << (ii/N_LnodalPos[0])%N_LnodalPos[2] << " L1: " <<  ii%N_LnodalPos[0] <<endl;

 //  for(ii=0; ii<N_L2ADIs; ii++)
 //    cout << ii << " Xpos: " << (ii/N_LnodalPos[0])/N_LnodalPos[1] << " L2: " << (ii/N_LnodalPos[0])%N_LnodalPos[1] << " L1: " <<  ii%N_LnodalPos[0] <<endl;

  for(i=0; i<N_L0ADIs; ++i)
  { 
   i_x =  (i/N_LnodalPos[1])/N_LnodalPos[2];
   i_l2 = (i/N_LnodalPos[1])%N_LnodalPos[2];
   i_l1 = i%N_LnodalPos[1];

   GlocalCellNo = coll->GetCell(i_x)->GetGlobalCellNo();

   Coords[0] =  Xpos[2*i_x];
   Coords[1] =  Xpos[2*i_x+1];   
   Coords[2] =  LnodalPos[2][i_l2];
   Coords[3] =  LnodalPos[1][i_l1];  
   Coords[4] = InitialPopulation[GlocalCellNo];  
   Sol_L0 = Sol_XposLNnodal +(i_x*N_Lnodal_All + i_l2*N_L1L0 + i_l1*N_LnodalPos[0] );

   //nodal inerpolation
   SystemADI[0]->Interpolate(N_Coord, Coords, Sol_L0, Exact);
  }
    
}


void TSystemADI1D_3L::Interpolate(DoubleFunctND *Exact)
{
  double Coords[6], *Sol_L0;
  int i, ii, i_x, i_l1, i_l2;

 //  for(ii=0; ii<N_L0ADIs; ii++)
 //    cout << ii << " Xpos: " << (ii/N_LnodalPos[1])/N_LnodalPos[2] << " L2: " << (ii/N_LnodalPos[1])%N_LnodalPos[2] << " L1: " <<  ii%N_LnodalPos[1] <<endl;

 //  for(ii=0; ii<N_L1ADIs; ii++)
 //    cout << ii << " Xpos: " << (ii/N_LnodalPos[0])/N_LnodalPos[2] << " L2: " << (ii/N_LnodalPos[0])%N_LnodalPos[2] << " L1: " <<  ii%N_LnodalPos[0] <<endl;

 //  for(ii=0; ii<N_L2ADIs; ii++)
 //    cout << ii << " Xpos: " << (ii/N_LnodalPos[0])/N_LnodalPos[1] << " L2: " << (ii/N_LnodalPos[0])%N_LnodalPos[1] << " L1: " <<  ii%N_LnodalPos[0] <<endl;

  for(i=0; i<N_L0ADIs; ++i)
  { 
   i_x =  (i/N_LnodalPos[1])/N_LnodalPos[2];
   i_l2 = (i/N_LnodalPos[1])%N_LnodalPos[2];
   i_l1 = i%N_LnodalPos[1];

   Coords[0] =  Xpos[2*i_x];
   Coords[1] =  Xpos[2*i_x+1];   
   Coords[2] =  LnodalPos[2][i_l2];
   Coords[3] =  LnodalPos[1][i_l1];  
   Coords[4] = -1; //N_Coord=6, so this is not used here
   Sol_L0 = Sol_XposLNnodal +(i_x*N_Lnodal_All + i_l2*N_L1L0 + i_l1*N_LnodalPos[0] );
  //  cout << i << " N_L0ADIs: " << Sol_L0[0] << " : " << Sol_L0[1] << " : "<< Sol_L0[2] << " : "<< Sol_L0[3] << endl;
   SystemADI[0]->Interpolate(N_Coord, Coords, Sol_L0, Exact);
  }
    
}

// XnodalL2L1L0nodal-->L2L1L0dofXdof
void TSystemADI1D_3L::XnodalLnodalToLdofXdof(int MaxN_PtsForNodal, double *Sol_LdofXdof)
{
  TBaseCell *cell;
  FE2D FEId_X;
  TFE2D *Element_X;
  TNodalFunctional2D *nf_X;

  // assume all  ScalarFunctions use same fespace2D
  TFESpace2D *FE_Space2D = SpatialFeFunction->GetFESpace2D();
  int N_XDof = FE_Space2D->GetN_DegreesOfFreedom();
  int *BeginIndex = FE_Space2D->GetBeginIndex();
  int *GlobalNumbers = FE_Space2D->GetGlobalNumbers();
  int *IntlPtIndexOfPts = FE_Space2D->GetIntlPtIndexOfPts();
  double *xi, *eta, *val = new double[MaxN_PtsForNodal*N_Lnodal_All];
  TCollection *Coll_X = FE_Space2D->GetCollection();
  int N_X, N_Cells_X = Coll_X->GetN_Cells();
  int i, j, k, l, m, n, *PtIndexLoc, *DOF, N_LocalDOFs, disp = 0;
  double *Val_L2L1L0, *Val_L1L0, *Val_L0, *PtValues, *Sol_Level;
  double FunctionalValues[MaxN_PointsForNodal2D];

  // memset(Sol_LdofXdof, 0, N_XDof*N_Lnodal_All*SizeOfDouble);

  double *Sol_ALL =  Sol_XnodalLDof[0];

  // x-dir are solved for LDofs
  // XnodalLnodal->XnodalLdof
  this->Nodal2DOF(0, Sol_ALL);

  for(int i=0;i<N_Cells_X;++i)
  {
    cell = Coll_X->GetCell(i);
    FEId_X = FE_Space2D->GetFE2D(i, cell);
    Element_X = TFEDatabase2D::GetFE2D(FEId_X);
    nf_X = Element_X->GetNodalFunctional2D();
    nf_X->GetPointsForAll(N_X, xi, eta);
    DOF = GlobalNumbers + BeginIndex[i];
    N_LocalDOFs = Element_X->GetN_DOF();
    PtIndexLoc = IntlPtIndexOfPts + disp;

    // collect point values for all level
    for(j=0;j<N_X;j++)
    {
      k = PtIndexLoc[j];
      Val_L2L1L0  =  Sol_ALL + k*N_Lnodal_All;

      for(l=0;l<N_LnodalPos[2]; l++)
      {
        Val_L1L0 = Val_L2L1L0 + l*N_L1L0;

        for(m=0;m<N_LnodalPos[1]; m++)
        {
          Val_L0 = Val_L1L0 + m*N_LnodalPos[0];

          for(n=0;n<N_LnodalPos[0]; n++)
          {
           val[(l*N_L1L0 + m*N_LnodalPos[0] + n)*N_X  + j] = Val_L0[n];
            // if(maxval<Val_L0[n]) maxval=Val_L0[n];
          //  cout << "pos:" << (l*N_L1L0 + m*N_LnodalPos[0] + n)*N_X  + j << endl;
          }
        }
      }
    } // for(j=0;j<N_Points;j++)


    for(j=0;j<N_Lnodal_All;j++)
    {
      PtValues = val + j*N_X;
      Sol_Level = Sol_LdofXdof+ j*N_XDof;

      nf_X->GetAllFunctionals(Coll_X, (TGridCell *)cell, PtValues, FunctionalValues);

      for(k=0;k<N_LocalDOFs;k++)
        Sol_Level[DOF[k]] = FunctionalValues[k];
    }

  disp +=N_X;
  } // for(i=0
 
  // OutPut("maxval : " << maxval <<endl;)
  delete [] val;
}


// L2L1L0nodalXDof-->XnodalL2L1L0nodal
void TSystemADI1D_3L::LdofXdofToXnodalLdof(double *Sol_LdofXdof)
{
  TBaseCell *cell;
  FE2D FEId_X;
  TFE2D *Element_X;
  TNodalFunctional2D *nf_X;
  TBaseFunct2D *bf_X;

  // assume all  ScalarFunctions use same fespace2D
  TFESpace2D *FE_Space2D = SpatialFeFunction->GetFESpace2D();
  int N_XDof = FE_Space2D->GetN_DegreesOfFreedom();
  int *BeginIndex = FE_Space2D->GetBeginIndex();
  int *GlobalNumbers = FE_Space2D->GetGlobalNumbers();
  int *IntlPtIndexOfPts = FE_Space2D->GetIntlPtIndexOfPts();
  double *xi, *eta;
  TCollection *Coll_X = FE_Space2D->GetCollection();
  int N_X, N_Cells_X = Coll_X->GetN_Cells();
  int i, j, k, l, m, n, p, q, *PtIndexLoc, *DOF, N_LocalDOFs, disp = 0;
  double *Val_L2L1L0, *Val_L1L0, *Val_L0, *PtValues, *Sol_Level;
  double FunctionalValues[MaxN_PointsForNodal2D];
  double maxval=0;
  double BasisValues[MaxN_BaseFunctions2D], *NodalValues;
  int *Incident_L0, *Incident_L1L0, *Incident_L2L1L0;

  memset(IncidentArray, 0, N_XposLnodal_All*SizeOfInt);
  // memset(Sol_XposLNnodal, 0, N_XposLnodal_All*SizeOfDouble);

  double *Sol_ALL =  Sol_XnodalLDof[0];
  memset(Sol_ALL, 0, N_XposLnodal_All*SizeOfDouble);

  for(int i=0;i<N_Cells_X;++i)
  {
    cell = Coll_X->GetCell(i);
    FEId_X = FE_Space2D->GetFE2D(i, cell);
    Element_X = TFEDatabase2D::GetFE2D(FEId_X);
    nf_X = Element_X->GetNodalFunctional2D();
    nf_X->GetPointsForAll(N_X, xi, eta);
    bf_X = Element_X->GetBaseFunct2D();
    DOF = GlobalNumbers + BeginIndex[i];
    N_LocalDOFs = Element_X->GetN_DOF();
    PtIndexLoc = IntlPtIndexOfPts + disp;

    // collect point values for all level
    for(j=0;j<N_X;j++)
    {
      bf_X->GetDerivatives(D00, xi[j], eta[j], BasisValues);      
      k = PtIndexLoc[j];
      // Val_L2L1L0  =  Sol_XposLNnodal + k*N_Lnodal_All;
      Val_L2L1L0  =  Sol_ALL + k*N_Lnodal_All;     
      Incident_L2L1L0 = IncidentArray + k*N_Lnodal_All;

      for(l=0;l<N_LnodalPos[2]; l++)
      {
        Val_L1L0 = Val_L2L1L0 + l*N_L1L0;
        Incident_L1L0 = Incident_L2L1L0+ l*N_L1L0;

        for(m=0;m<N_LnodalPos[1]; m++)
        {
          Val_L0 = Val_L1L0 + m*N_LnodalPos[0];
          Incident_L0 = Incident_L1L0+ m*N_LnodalPos[0];

          for(n=0;n<N_LnodalPos[0]; n++)
          {
           Incident_L0[n]++;
           NodalValues = Sol_LdofXdof + (l*N_L1L0 + m*N_LnodalPos[0] + n)*N_XDof;
           for(p=0;p<N_LocalDOFs;p++)
            {
              Val_L0[n] += BasisValues[p]*NodalValues[DOF[p]];
            }
          //  cout << "pos:" << (l*N_L1L0 + m*N_LnodalPos[0] + n)*N_X  + j << endl;
          }
        }
      }
    } // for(j=0;j<N_Points;j++)
   disp +=N_X;
  }

  for(i=0;i<N_Xpos;i++)
  {
   Incident_L2L1L0 = IncidentArray + i*N_Lnodal_All;
  //  Val_L2L1L0  =  Sol_XposLNnodal + i*N_Lnodal_All;
   Val_L2L1L0  =  Sol_ALL + i*N_Lnodal_All;

   for(l=0;l<N_LnodalPos[2]; l++)
   {
    Val_L1L0 = Val_L2L1L0 + l*N_L1L0;
    Incident_L1L0 = Incident_L2L1L0+ l*N_L1L0;

    for(m=0;m<N_LnodalPos[1]; m++)
    {
      Val_L0 = Val_L1L0 + m*N_LnodalPos[0];
      Incident_L0 = Incident_L1L0+ m*N_LnodalPos[0];

      for(n=0;n<N_LnodalPos[0]; n++)
      {
       Val_L0[n] /= (double)Incident_L0[n];
      }
    }
   }
  }

  // Sol_ALL will be converted to Sol_XposLNnodal
  // this->DOF2Nodal(0, Sol_ALL); // not needed as IntL is computed using Sol_XnodalLdof
}

int TSystemADI1D_3L::Nodal2DOF(int ADI_Idx, double *Sol_XposLNnodLOwnDof)
{
 int m = SystemADI[ADI_Idx]->Nodal2DOF(Sol_XposLNnodLOwnDof);
 return m;
} //Nodal2DOF


int TSystemADI1D_3L::DOF2Nodal(int ADI_Idx, double *Sol_XposLNnodLOwnDof)
 {
 int m = SystemADI[ADI_Idx]->DOF2Nodal(Sol_XposLNnodLOwnDof);
 return m;
 }

// XnodalL2L1L0nodal--> XnodalLnLnLADI_Idxnodal
int TSystemADI1D_3L::CopySolToInternal(int ADI_Idx)
{
 int i, l, m, n;
 double *sol_L2L1L0, *sol_L1L0, *sol_L0, *sol_L2L0L1, *sol_L0L1, *sol_L0L1L2;

 //Lo is internal in XnodalL2L1L0nodal nothing to do
 if(ADI_Idx==0) 
  { return 0; }
 else if(ADI_Idx==1) // L2L1L0 ---> L2L0L1
 {
  //memset(Sol_XposL2L0L1nodal, 0, N_XposLnodal_All*SizeOfDouble); since assigning, not necessary to initiate
  for(i=0;i<N_Xpos;i++)
  {
   sol_L2L1L0  =  Sol_XposLNnodal + i*N_Lnodal_All;
   sol_L2L0L1  =  Sol_XposL2L0L1nodal + i*N_Lnodal_All;
   for(l=0;l<N_LnodalPos[2]; l++)
   {
    sol_L1L0 = sol_L2L1L0 + l*N_L1L0;
    sol_L0L1 = sol_L2L0L1 + l*N_L1L0;    
    for(m=0;m<N_LnodalPos[1]; m++)
    {
      sol_L0 = sol_L1L0 + m*N_LnodalPos[0];
      for(n=0;n<N_LnodalPos[0]; n++)
      {
       sol_L0L1[n*N_LnodalPos[1] + m] = sol_L0[n];
      }
    }
   }
  }
 }
 else if(ADI_Idx==2) // L2L1L0 ---> L0L1L2
 {
  // memset(Sol_XposL0L1L2nodal, 0, N_XposLnodal_All*SizeOfDouble);  since assigning, not necessary to initiate
  for(i=0;i<N_Xpos;i++)
  {
   sol_L2L1L0  =  Sol_XposLNnodal + i*N_Lnodal_All;
   sol_L0L1L2  =  Sol_XposL0L1L2nodal + i*N_Lnodal_All;
   for(l=0;l<N_LnodalPos[2]; l++)
   {   
    sol_L1L0 = sol_L2L1L0 + l*N_L1L0;     
    for(m=0;m<N_LnodalPos[1]; m++)
    {
     sol_L0 = sol_L1L0 + m*N_LnodalPos[0];      
     for(n=0;n<N_LnodalPos[0]; n++)
     {
      sol_L0L1L2[n*N_L2L1 +   m*N_LnodalPos[2] + l] = sol_L0[n];
     }
    }
   }
  }
 }
 else
 {
  ErrMsg("No internal L with Index " <<ADI_Idx <<endl);
  exit(0);
 }
return 0;
} // MakeLInternalAs


// XnodalLnLnLADI_Idxnodal -->  XnodalL2L1L0nodal
int TSystemADI1D_3L::CopySolFromInternal(int ADI_Idx)
{
 int i, l, m, n;
 double *sol_L2L1L0, *sol_L1L0, *sol_L0, *sol_L2L0L1, *sol_L0L1, *sol_L0L1L2;

 //Lo is internal in XnodalL2L1L0nodal nothing to do
 if(ADI_Idx==0) 
  { return 0; }
 else if(ADI_Idx==1) //L2L0L1 --> L2L1L0
 {
  for(i=0;i<N_Xpos;i++)
  {
   sol_L2L1L0  =  Sol_XposLNnodal + i*N_Lnodal_All;
   sol_L2L0L1  =  Sol_XposL2L0L1nodal + i*N_Lnodal_All;
   for(l=0;l<N_LnodalPos[2]; l++)
   {
    sol_L1L0 = sol_L2L1L0 + l*N_L1L0;
    sol_L0L1 = sol_L2L0L1 + l*N_L1L0;    
    for(m=0;m<N_LnodalPos[1]; m++)
    {
      sol_L0 = sol_L1L0 + m*N_LnodalPos[0];
      for(n=0;n<N_LnodalPos[0]; n++)
      {
       sol_L0[n] = sol_L0L1[n*N_LnodalPos[1] + m];
      }
    }
   }
  }
 }
 else if(ADI_Idx==2) // L0L1L2 --> L2L1L0
 {
  for(i=0;i<N_Xpos;i++)
  {
   sol_L2L1L0  =  Sol_XposLNnodal + i*N_Lnodal_All;
   sol_L0L1L2  =  Sol_XposL0L1L2nodal + i*N_Lnodal_All;
   for(l=0;l<N_LnodalPos[2]; l++)
   {   
    sol_L1L0 = sol_L2L1L0 + l*N_L1L0;     
    for(m=0;m<N_LnodalPos[1]; m++)
    {
     sol_L0 = sol_L1L0 + m*N_LnodalPos[0];      
     for(n=0;n<N_LnodalPos[0]; n++)
     {
      sol_L0[n] = sol_L0L1L2[n*N_L2L1 +   m*N_LnodalPos[2] + l];
     }
    }
   }
  }
 }
 else
 {
  ErrMsg("No internal L with Index " <<ADI_Idx <<endl);
  exit(0);
 }
 return 0;
} // MakeLInternalAs

// int_Omega_L(I) for a given x-pos
double TSystemADI1D_3L::IntL(DoubleFunctND *GetGamma_Q, double *IntValue, double * NucValue, double *LDistXSum)
{

 int i, i1, i2, j, k, m, n, N_X, ii;

 TFESpace1D *FeSpace_L0 = SystemADI[0]->GetFeSpace1D();
 TFESpace1D *FeSpace_L1 = SystemADI[1]->GetFeSpace1D();
 TFESpace1D *FeSpace_L2 = SystemADI[2]->GetFeSpace1D();

  TCollection *L0Coll = FeSpace_L0->GetCollection();
  int *L0BeginIndex = FeSpace_L0->GetBeginIndex();
  int *L0GlobalNumbers = FeSpace_L0->GetGlobalNumbers();
  int N_L0Cells = L0Coll->GetN_Cells();

  TCollection *L1Coll = FeSpace_L1->GetCollection();
  int *L1BeginIndex = FeSpace_L1->GetBeginIndex();
  int *L1GlobalNumbers = FeSpace_L1->GetGlobalNumbers();
  int N_L1Cells = L1Coll->GetN_Cells();

  TCollection *L2Coll = FeSpace_L2->GetCollection();
  int *L2BeginIndex = FeSpace_L2->GetBeginIndex();
  int *L2GlobalNumbers = FeSpace_L2->GetGlobalNumbers();
  int N_L2Cells = L2Coll->GetN_Cells();
  FE1D LocalUsedElements, CurrentElement;
  bool Needs2ndDer[1], SecondDer[1];
  SecondDer[0] = FALSE;
  Needs2ndDer[0] = FALSE;

  //first check how many quad pts  in the l1cell
  TBaseCell *cell_L2, *cell_L1,  *cell = L0Coll->GetCell(0);
  LocalUsedElements = FeSpace_L0->GetFE1D(0, cell);
  TFE1D *Element = TFEDatabase2D::GetFE1D(LocalUsedElements);
  TBaseFunct1D *bf = Element->GetBaseFunct1D();
  int l = bf->GetPolynomialDegree();
  QuadFormula1D LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(2*l);
  TQuadFormula1D *qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
  int N_QuadPts, N_QuadPts_L1, N_QuadPts_L2, *L0DOF, *L1DOF, *L2DOF;
  double *Weights, *zeta;
  double *Weights_L2, *Weights_L1, *zeta_L2, *zeta_L1;
  double Mult, value, *orgD0, **origvaluesD0;
  qf1->GetFormulaData(N_QuadPts, Weights, zeta);

  //for I distribution output 
  QuadFormula1D NodalLineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(1); // Gauss1()  
  TQuadFormula1D *Nodalqf1 = TFEDatabase2D::GetQuadFormula1D(NodalLineQuadFormula);
  int NodalN_QuadPts;
  double *NodalWeights, *Nodalzeta, **NodalorigvaluesD0_temp, NodalorigvaluesD0[2][5], *LDistXSumL1L0, *LDistXSumL0;
  Nodalqf1->GetFormulaData(NodalN_QuadPts, NodalWeights, Nodalzeta);

  int N_BaseFunct, N_Sets=1;
  BaseFunct1D BaseFunct_ID, BaseFunct[1];
  TRefTrans1D *F_K;
  TNodalFunctional1D *nf;
  double *xi, *eta;
  int N_NodalPoints;

  double *Sol_L1quadL0quadL2nodal, *Sol_L0quadL2nodal;  
  double *IntValue_L1quad, *IntValue_L0quadL1quad;
  double *Sol_L2nodalL0quadL1nodal, *Sol_L0quadL1nodal;
  double *Sol_L2L1L0, *Sol_L1L0, *Sol_L0, *Sol_L1, *Sol_L2;
  double val, *IntValue_L0quad; 
  double Coords[6], IVal[6], value1;
  double *NucValue_L0quadL1quad, *NucValue_L0quad, *NucValue_L1quad;
  Coords[4] = -1;
  
  memset(IntValue, 0, N_Xpos*SizeOfDouble);
  memset(NucValue, 0, N_Xpos*SizeOfDouble);
  if(LDistXSum)
  { memset(LDistXSum, 0, N_Lnodal_All*SizeOfDouble); }

  /** Assumed that xnodal & xdof are same */
  /**  otherwise, IntL should be calculated using  XdifLdof */
  double *Sol_ALL =  Sol_XnodalLDof[0]; 
  // this->Nodal2DOF(0, Sol_ALL); // already done in XnodalLnodalToLdofXdof or 

  // assumed that nodal points and dof pos are same
  // N_L0Cells = 10;
  for(i=0; i<N_L0Cells; ++i)
  {
    cell = L0Coll->GetCell(i);
    LocalUsedElements = FeSpace_L0->GetFE1D(i, cell);
    Element = TFEDatabase2D::GetFE1D(LocalUsedElements);
    bf = Element->GetBaseFunct1D();
    N_BaseFunct = Element->GetN_DOF();
    BaseFunct_ID = Element->GetBaseFunct1D_ID();
    l = bf->GetPolynomialDegree();
    LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(2*l);
    qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
    qf1->GetFormulaData(N_QuadPts, Weights, zeta);
    F_K = TFEDatabase2D::GetRefTrans1D(LineAffin);
    ((TLineAffin *)F_K)->SetCell(cell);
    BaseFunct[0] = BaseFunct_ID;    

    //for I distribution output. This is needed only for dG as N_j(\phi_i) != 1 in dG
    if(LDistXSum)
    {
    if(NodalN_QuadPts!=1)
     {
       cout<<"IntL SystemADI1D_3L.C Error NodalN_QuadPts!=1" <<endl;
       exit(0);
     }
     ((TLineAffin *)F_K)->GetOrigFromRef(NodalN_QuadPts, Nodalzeta, NodalX0, NodalAbsDetjk);
     ((TLineAffin *)F_K)->GetOrigValues(N_Sets, BaseFunct, NodalN_QuadPts, Nodalzeta,  NodalLineQuadFormula,  Needs2ndDer); 
     NodalorigvaluesD0_temp=TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D0);

     for(m=0; m<NodalN_QuadPts; m++)
      for(n=0;n<N_BaseFunct;n++)
        NodalorigvaluesD0[m][n] = NodalorigvaluesD0_temp[m][n];

    } //if(LDistXSum)

    ((TLineAffin *)F_K)->GetOrigFromRef(N_QuadPts, zeta, X0, AbsDetjk);
    ((TLineAffin *)F_K)->GetOrigValues(N_Sets, BaseFunct, N_QuadPts, zeta,  LineQuadFormula,  Needs2ndDer);
    origvaluesD0=TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D0);
    L0DOF = L0GlobalNumbers + L0BeginIndex[i];

    // cout << "Test 1" <<endl;
    for(ii=0; ii<N_Xpos; ii++)
    {
    //  Sol_L2L1L0 = Sol_XposLNnodal + ii*N_Lnodal_All;
     Sol_L2L1L0 = Sol_ALL + ii*N_Lnodal_All;

     Sol_L2nodalL0quadL1nodal = Sol_XposL2nodalL0quadL1nodal + ii*N_LnodalPos[2]*N_QuadPts*N_LnodalPos[1];
     
     for(j=0; j<N_LnodalPos[2]; j++)
      {
       Sol_L1L0 = Sol_L2L1L0 +j*N_L1L0;
       Sol_L0quadL1nodal = Sol_L2nodalL0quadL1nodal  +j*N_QuadPts*N_LnodalPos[1];

       if(LDistXSum)
         { LDistXSumL1L0 = LDistXSum+j*N_LnodalPos[1]*N_L0Cells;}//assumed that NodalN_QuadPts=1, otherwise *N_L0Cells*NodalN_QuadPts

       for(k=0; k<N_LnodalPos[1]; k++)
       {
        Sol_L0 = Sol_L1L0 +k*N_LnodalPos[0];

       if(LDistXSum) //int_K^{l^d}_i(I) //I value in the cell 
        {
         LDistXSumL0 = LDistXSumL1L0+(k*N_L0Cells+i); //assumed that NodalN_QuadPts=1, otherwise k*N_L0Cells+i*NodalN_QuadPts
         for(m=0; m<NodalN_QuadPts; m++)
         {
          value = 0.;
          for(n=0;n<N_BaseFunct;n++)
            value += Sol_L0[L0DOF[n]]*NodalorigvaluesD0[m][n]; 

          if (value<0.) value=0; // due to spurious oscillation, it may be zero
          LDistXSumL0[m] += value*NodalWeights[m]*NodalAbsDetjk[m];
         }
        }

        for(m=0; m<N_QuadPts; m++)
        {
         orgD0 = origvaluesD0[m];
         value = 0.;
         for(n=0;n<N_BaseFunct;n++)
         {   
          value += Sol_L0[L0DOF[n]]*orgD0[n]; 
         } // n
         Sol_L0quadL1nodal[m*N_LnodalPos[1]   + k] = value;
        } // m
       } // j
     }// for j
    } // for ii

   memset(IntValue_XposL0quad, 0, N_Xpos*N_QuadPts*SizeOfDouble);   
   memset(NucValue_XposL0quad, 0, N_Xpos*N_QuadPts*SizeOfDouble); 

   //  memset(Volume_XposL0quad, 0, N_QuadPts*SizeOfDouble);   
   //compute L1Quad from L1dof
   for(i1=0; i1<N_L1Cells; ++i1)
   {
    cell_L1 = L1Coll->GetCell(i1);
    LocalUsedElements = FeSpace_L1->GetFE1D(i1, cell_L1);
    Element = TFEDatabase2D::GetFE1D(LocalUsedElements);
    bf = Element->GetBaseFunct1D();
    N_BaseFunct = Element->GetN_DOF();
    BaseFunct_ID = Element->GetBaseFunct1D_ID();
    l = bf->GetPolynomialDegree();
    LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(2*l);
    qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
    qf1->GetFormulaData(N_QuadPts_L1, Weights_L1, zeta_L1);
    F_K = TFEDatabase2D::GetRefTrans1D(LineAffin);
    ((TLineAffin *)F_K)->SetCell(cell_L1);
    ((TLineAffin *)F_K)->GetOrigFromRef(N_QuadPts_L1, zeta_L1, X1, AbsDetjk_L1);
    BaseFunct[0] = BaseFunct_ID;
    ((TLineAffin *)F_K)->GetOrigValues(N_Sets, BaseFunct, N_QuadPts_L1, zeta_L1,  LineQuadFormula,  Needs2ndDer);
    origvaluesD0=TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D0);
    L1DOF = L1GlobalNumbers + L1BeginIndex[i1];
  
    for(ii=0; ii<N_Xpos; ii++)
    {
     Sol_L2nodalL0quadL1nodal = Sol_XposL2nodalL0quadL1nodal + ii*N_LnodalPos[2]*N_QuadPts*N_LnodalPos[1];
     Sol_L1quadL0quadL2nodal = Sol_XposL1quadL0quadL2nodal  + ii*N_QuadPts*N_QuadPts*N_LnodalPos[2];

     for(j=0; j<N_LnodalPos[2]; j++)
      {
       Sol_L0quadL1nodal = Sol_L2nodalL0quadL1nodal  +j*N_QuadPts*N_LnodalPos[1];

       for(k=0; k<N_QuadPts; k++)
       {
        Sol_L1 = Sol_L0quadL1nodal +k*N_LnodalPos[1];
        for(m=0; m<N_QuadPts_L1; m++)
        {
         orgD0 = origvaluesD0[m];
         value = 0.;
         for(n=0;n<N_BaseFunct;n++)
         {   
          value += Sol_L1[L1DOF[n]]*orgD0[n]; 
          // value = 10.;       
         } // n
        Sol_L1quadL0quadL2nodal[m*N_QuadPts*N_LnodalPos[2] +   k*N_LnodalPos[2] + j] = value;
        } // m
       } // j
     }// for j
    } // for ii

   memset(IntValue_XposL0quadL1quad, 0, N_Xpos*N_QuadPts*N_QuadPts*SizeOfDouble);
   memset(NucValue_XposL0quadL1quad, 0, N_Xpos*N_QuadPts*N_QuadPts*SizeOfDouble);

   //  memset(Volume_XposL0quadL1quad, 0, 1*N_QuadPts*N_QuadPts*SizeOfDouble);   
   //compute L2Quad from L1dof
   for(i2=0; i2<N_L2Cells; ++i2)
   {
    cell_L2 = L2Coll->GetCell(i2);
    LocalUsedElements = FeSpace_L2->GetFE1D(i2, cell_L2);
    Element = TFEDatabase2D::GetFE1D(LocalUsedElements);
    bf = Element->GetBaseFunct1D();
    N_BaseFunct = Element->GetN_DOF();
    BaseFunct_ID = Element->GetBaseFunct1D_ID();
    l = bf->GetPolynomialDegree();
    LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(2*l);
    qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
    qf1->GetFormulaData(N_QuadPts_L2, Weights_L2, zeta_L2);
    F_K = TFEDatabase2D::GetRefTrans1D(LineAffin);
    ((TLineAffin *)F_K)->SetCell(cell_L2);
    ((TLineAffin *)F_K)->GetOrigFromRef(N_QuadPts_L2, zeta_L2, X2, AbsDetjk_L2);
    BaseFunct[0] = BaseFunct_ID;
    ((TLineAffin *)F_K)->GetOrigValues(N_Sets, BaseFunct, N_QuadPts_L2, zeta_L2,  LineQuadFormula,  Needs2ndDer);
    origvaluesD0=TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D0);
    L2DOF = L2GlobalNumbers + L2BeginIndex[i2];


    for(ii=0; ii<N_Xpos; ii++)
    {
     Sol_L1quadL0quadL2nodal = Sol_XposL1quadL0quadL2nodal  + ii*N_QuadPts_L1*N_QuadPts*N_LnodalPos[2];
     IntValue_L0quadL1quad = IntValue_XposL0quadL1quad +  ii*N_QuadPts_L1*N_QuadPts;
     NucValue_L0quadL1quad = NucValue_XposL0quadL1quad +  ii*N_QuadPts_L1*N_QuadPts;

      Coords[0] =  Xpos[2*ii];
      Coords[1] = Xpos[2*ii+1];

     for(j=0; j<N_QuadPts_L1; j++)
      {
       Sol_L0quadL2nodal = Sol_L1quadL0quadL2nodal +j*N_QuadPts*N_LnodalPos[2];
       Coords[3] = X1[j];  

       for(k=0; k<N_QuadPts; k++)
       {
        Coords[5] = X0[k]; 
        Sol_L2 = Sol_L0quadL2nodal +k*N_LnodalPos[2];
        value = 0;
        value1 = 0;      
        // vol = 0;
        for(m=0; m<N_QuadPts_L2; m++)
        {
         Coords[2] = X2[m];  
         orgD0 = origvaluesD0[m];
         val = 0.;
         for(n=0;n<N_BaseFunct;n++)
         {   
          val += Sol_L2[L2DOF[n]]*orgD0[n]; 
         } // n

         //GetExactSol
         if(GetGamma_Q)
         {
          GetGamma_Q(N_Coord, Coords, IVal);
          value1 += (1.- IVal[0])*val*Weights_L2[m]*AbsDetjk_L2[m];
         }
         value += val*Weights_L2[m]*AbsDetjk_L2[m];
        //  value += Weights_L2[m]*AbsDetjk_L2[m];  
        //  if(ii==0)  
          // vol +=Weights_L2[m]*AbsDetjk_L2[m];
        //  value += 0;
        } // m
        //L2 cell contribution to each L0 & L1 Quad points
        IntValue_L0quadL1quad[k*N_QuadPts_L1 + j] +=  value;
        NucValue_L0quadL1quad[k*N_QuadPts_L1 + j] +=  value1;        
        // Volume_XposL0quadL1quad[k*N_QuadPts_L1 + j] += vol;
        // cout << "value  " << value << " vol  " << vol << endl;
       } // k

     }// for j
    } // for ii 
   } // i2

   //int over L1
   for(ii=0; ii<N_Xpos; ii++)
   {
    IntValue_L0quadL1quad = IntValue_XposL0quadL1quad +  ii*N_QuadPts_L1*N_QuadPts; 
    NucValue_L0quadL1quad = NucValue_XposL0quadL1quad +  ii*N_QuadPts_L1*N_QuadPts;  
    IntValue_L0quad = IntValue_XposL0quad +  ii*N_QuadPts;
    NucValue_L0quad = NucValue_XposL0quad +  ii*N_QuadPts;

    for(j=0; j<N_QuadPts; j++)
     {
      IntValue_L1quad = IntValue_L0quadL1quad + j*N_QuadPts_L1;
      NucValue_L1quad = NucValue_L0quadL1quad + j*N_QuadPts_L1;      
      // Volume_L1quad  = Volume_XposL0quadL1quad   + j*N_QuadPts_L1;    
      value = 0.;
      value1 = 0.;
      // vol=0;  
      for(k=0; k<N_QuadPts_L1; k++)
      {
       value += IntValue_L1quad[k]*Weights_L1[k]*AbsDetjk_L1[k];
       value1 += NucValue_L1quad[k]*Weights_L1[k]*AbsDetjk_L1[k];    
      //  if(ii==0)
      //  vol +=   Volume_L1quad[k]*Weights_L1[k]*AbsDetjk_L1[k];      
          // cout<< "value " << Volume_L1quad[k] << endl;
      }    
      IntValue_L0quad[j] += value;
      NucValue_L0quad[j] += value1;     
      // Volume_XposL0quad[j] += vol;
      // cout << "value  " << value << " vol  " << vol << endl;      
     }
   }
  } // i1

   //int over L0
   for(ii=0; ii<N_Xpos; ii++)
   {
    IntValue_L0quad = IntValue_XposL0quad +  ii*N_QuadPts;  
    NucValue_L0quad = NucValue_XposL0quad +  ii*N_QuadPts;  

     value = 0.;  
     value1 = 0.;
    //  vol = 0.;
    for(j=0; j<N_QuadPts; j++)
     {
      value += IntValue_L0quad[j]*Weights[j]*AbsDetjk[j];
      value1 += NucValue_L0quad[j]*Weights[j]*AbsDetjk[j];
      // if(ii==0)
      // IntlVolume += Volume_XposL0quad[j]*Weights[j]*AbsDetjk[j];
     }
    // cout << "value  " << value << " vol  " << vol << endl;
    IntValue[ii] += value;
    NucValue[ii] += value1;
   }

 }// i

  // double min=1e8, max=-1e8;
   for(j=0; j<N_Xpos; j++)
   {
    // IntValue[j] *=200./TDatabase::ParamDB->REACTOR_P13;
    // NucValue[j] *=200./TDatabase::ParamDB->REACTOR_P13;
  //   // if(min>IntValue[j] ) min=IntValue[j];
  //   // if(max<IntValue[j] ) max=IntValue[j];

     if(IntValue[j]<0)
     {
       IntValue[j] = 0;  // due to spurious oscillation, it may be zero
      //  NucValue[j] = 0;  // due to spurious oscillation, it may be zero
      // cout<< "Population Cannot be negative, check the oscillation in the solver: IntValue : " << IntValue[j] << endl;
      // exit(0);
     }

     if(NucValue[j]<0.)
         NucValue[j] = 0;  
 

   }
//  cout<< " min : " << min <<   " max : " << max << endl;
// exit(0);


  // delete [] Sol_XposL2nodalL0quadL1nodal;
  // delete [] Sol_XposL1quadL0quadL2nodal;
  // delete [] IntValue_XposL0quadL1quad;
  // delete [] IntValue_XposL0quad;
  // delete [] Volume_XposL0quadL1quad;
  // delete [] Volume_XposL0quad;

  return 0.;
} // IntL()


// int_Omega_L(I) for a given x-pos
double TSystemADI1D_3L::GetErrorAllXpos(DoubleFunctND *Exact, double *L2Errors)
{
 int i, i1, i2, j, k, m, n, N_X, ii;

 TFESpace1D *FeSpace_L0 = SystemADI[0]->GetFeSpace1D();
 TFESpace1D *FeSpace_L1 = SystemADI[1]->GetFeSpace1D();
 TFESpace1D *FeSpace_L2 = SystemADI[2]->GetFeSpace1D();

 TCollection *L0Coll = FeSpace_L0->GetCollection();
 int *L0BeginIndex = FeSpace_L0->GetBeginIndex();
 int *L0GlobalNumbers = FeSpace_L0->GetGlobalNumbers();
 int N_L0Cells = L0Coll->GetN_Cells();

 TCollection *L1Coll = FeSpace_L1->GetCollection();
 int *L1BeginIndex = FeSpace_L1->GetBeginIndex();
 int *L1GlobalNumbers = FeSpace_L1->GetGlobalNumbers();
 int N_L1Cells = L1Coll->GetN_Cells();

 TCollection *L2Coll = FeSpace_L2->GetCollection();
 int *L2BeginIndex = FeSpace_L2->GetBeginIndex();
 int *L2GlobalNumbers = FeSpace_L2->GetGlobalNumbers();
 int N_L2Cells = L2Coll->GetN_Cells();
 FE1D LocalUsedElements, CurrentElement;
 bool Needs2ndDer[1], SecondDer[1];
 SecondDer[0] = FALSE;
 Needs2ndDer[0] = FALSE;

 //first check how many quad pts  in the l1cell
 TBaseCell *cell_L2, *cell_L1,  *cell = L0Coll->GetCell(0);
 LocalUsedElements = FeSpace_L0->GetFE1D(0, cell);
 TFE1D *Element = TFEDatabase2D::GetFE1D(LocalUsedElements);
  TBaseFunct1D *bf = Element->GetBaseFunct1D();
  int l = bf->GetPolynomialDegree();
  QuadFormula1D LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(2*l);
  TQuadFormula1D *qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
  int N_QuadPts, N_QuadPts_L1, N_QuadPts_L2, *L0DOF, *L1DOF, *L2DOF;
  double *Weights, *zeta;
  double *Weights_L2, *Weights_L1, *zeta_L2, *zeta_L1;
  double Mult, value, *orgD0, **origvaluesD0;
  qf1->GetFormulaData(N_QuadPts, Weights, zeta);

  int N_BaseFunct, N_Sets=1;
  BaseFunct1D BaseFunct_ID, BaseFunct[1];
  TRefTrans1D *F_K;

  // double *Sol_XposL2nodalL0quadL1nodal = new double [N_Xpos*N_LnodalPos[2]*N_QuadPts*N_LnodalPos[1]];
  double *Sol_L2nodalL0quadL1nodal, *Sol_L0quadL1nodal;
  // double *Sol_XposL1quadL0quadL2nodal = new double [N_Xpos*N_QuadPts*N_QuadPts*N_LnodalPos[2]];
  double *Sol_L1quadL0quadL2nodal, *Sol_L0quadL2nodal;
  // double *IntValue_XposL0quadL1quad = new double[N_Xpos*N_QuadPts*N_QuadPts]; // [N_Xpos*N_QuadPts_L1*N_QuadPts];
  double *IntValue_L1quad, *IntValue_L0quadL1quad;

  double *Sol_L2L1L0, *Sol_L1L0, *Sol_L0, *Sol_L1, *Sol_L2;
  double val, *IntValue_L0quad; 
  double Coords[6], IVal[6];
  Coords[4] = -1;

  memset(L2Errors, 0, N_Xpos*SizeOfDouble);
    // TDatabase::ParamDB->P2=1;
  // assumed that nodal points and dof pos are same
  for(i=0; i<N_L0Cells; ++i)
  {
    cell = L0Coll->GetCell(i);
    LocalUsedElements = FeSpace_L0->GetFE1D(i, cell);
    Element = TFEDatabase2D::GetFE1D(LocalUsedElements);
    bf = Element->GetBaseFunct1D();
    N_BaseFunct = Element->GetN_DOF();
    BaseFunct_ID = Element->GetBaseFunct1D_ID();
    l = bf->GetPolynomialDegree();
    LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(2*l);
    qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
    qf1->GetFormulaData(N_QuadPts, Weights, zeta);
    F_K = TFEDatabase2D::GetRefTrans1D(LineAffin);
    ((TLineAffin *)F_K)->SetCell(cell);
    ((TLineAffin *)F_K)->GetOrigFromRef(N_QuadPts, zeta, X0, AbsDetjk);
    BaseFunct[0] = BaseFunct_ID;
    ((TLineAffin *)F_K)->GetOrigValues(N_Sets, BaseFunct, N_QuadPts, zeta,  LineQuadFormula,  Needs2ndDer);
    origvaluesD0=TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D0);
    L0DOF = L0GlobalNumbers + L0BeginIndex[i];


 
    //  for(m=0; m<N_QuadPts; m++)
    //    cout << zeta[m] << "Test  :" << X0[m] <<endl;
    //  cout<<endl;

 

    for(ii=0; ii<N_Xpos; ii++)
    {
     Sol_L2L1L0 = Sol_XposLNnodal + ii*N_Lnodal_All;
     Sol_L2nodalL0quadL1nodal = Sol_XposL2nodalL0quadL1nodal + ii*N_LnodalPos[2]*N_QuadPts*N_LnodalPos[1];
     
     for(j=0; j<N_LnodalPos[2]; j++)
      {
       Sol_L1L0 = Sol_L2L1L0 +j*N_L1L0;
       Sol_L0quadL1nodal = Sol_L2nodalL0quadL1nodal  +j*N_QuadPts*N_LnodalPos[1];
       for(k=0; k<N_LnodalPos[1]; k++)
       {
        Sol_L0 = Sol_L1L0 +k*N_LnodalPos[0];
        for(m=0; m<N_QuadPts; m++)
        {
         orgD0 = origvaluesD0[m];
         value = 0.;
         for(n=0;n<N_BaseFunct;n++)
         {   
          value += Sol_L0[L0DOF[n]]*orgD0[n]; 
         } // n
         Sol_L0quadL1nodal[m*N_LnodalPos[1]   + k] = value;
        } // m
       } // j
     }// for j
    } // for ii

   memset(IntValue_XposL0quad, 0, N_Xpos*N_QuadPts*SizeOfDouble);   
   //compute L1Quad from L1dof
   for(i1=0; i1<N_L1Cells; ++i1)
   {
    cell_L1 = L1Coll->GetCell(i1);
    LocalUsedElements = FeSpace_L1->GetFE1D(i1, cell_L1);
    Element = TFEDatabase2D::GetFE1D(LocalUsedElements);
    bf = Element->GetBaseFunct1D();
    N_BaseFunct = Element->GetN_DOF();
    BaseFunct_ID = Element->GetBaseFunct1D_ID();
    l = bf->GetPolynomialDegree();
    LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(2*l);
    qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
    qf1->GetFormulaData(N_QuadPts_L1, Weights_L1, zeta_L1);
    F_K = TFEDatabase2D::GetRefTrans1D(LineAffin);
    ((TLineAffin *)F_K)->SetCell(cell_L1);
    ((TLineAffin *)F_K)->GetOrigFromRef(N_QuadPts_L1, zeta_L1, X1, AbsDetjk_L1);
    BaseFunct[0] = BaseFunct_ID;
    ((TLineAffin *)F_K)->GetOrigValues(N_Sets, BaseFunct, N_QuadPts_L1, zeta_L1,  LineQuadFormula,  Needs2ndDer);
    origvaluesD0=TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D0);
    L1DOF = L1GlobalNumbers + L1BeginIndex[i1];
  
    for(ii=0; ii<N_Xpos; ii++)
    {
     Sol_L2nodalL0quadL1nodal = Sol_XposL2nodalL0quadL1nodal + ii*N_LnodalPos[2]*N_QuadPts*N_LnodalPos[1];
     Sol_L1quadL0quadL2nodal = Sol_XposL1quadL0quadL2nodal  + ii*N_QuadPts*N_QuadPts*N_LnodalPos[2];

     for(j=0; j<N_LnodalPos[2]; j++)
      {
       Sol_L0quadL1nodal = Sol_L2nodalL0quadL1nodal  +j*N_QuadPts*N_LnodalPos[1];

       for(k=0; k<N_QuadPts; k++)
       {
        Sol_L1 = Sol_L0quadL1nodal +k*N_LnodalPos[1];
        for(m=0; m<N_QuadPts_L1; m++)
        {
         orgD0 = origvaluesD0[m];
         value = 0.;
         for(n=0;n<N_BaseFunct;n++)
         {   
          value += Sol_L1[L1DOF[n]]*orgD0[n]; 
          // value = 1;       
         } // n
        Sol_L1quadL0quadL2nodal[m*N_QuadPts*N_LnodalPos[2] +   k*N_LnodalPos[2] + j] = value;
        } // m
       } // j
     }// for j
    } // for ii

   memset(IntValue_XposL0quadL1quad, 0, N_Xpos*N_QuadPts*N_QuadPts*SizeOfDouble);
   //compute L2Quad from L1dof
   for(i2=0; i2<N_L2Cells; ++i2)
   {
    cell_L2 = L2Coll->GetCell(i2);
    LocalUsedElements = FeSpace_L2->GetFE1D(i2, cell_L2);
    Element = TFEDatabase2D::GetFE1D(LocalUsedElements);
    bf = Element->GetBaseFunct1D();
    N_BaseFunct = Element->GetN_DOF();
    BaseFunct_ID = Element->GetBaseFunct1D_ID();
    l = bf->GetPolynomialDegree();
    LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(2*l);
    qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
    qf1->GetFormulaData(N_QuadPts_L2, Weights_L2, zeta_L2);
    F_K = TFEDatabase2D::GetRefTrans1D(LineAffin);
    ((TLineAffin *)F_K)->SetCell(cell_L2);
    ((TLineAffin *)F_K)->GetOrigFromRef(N_QuadPts_L2, zeta_L2, X2, AbsDetjk_L2);
    BaseFunct[0] = BaseFunct_ID;
    ((TLineAffin *)F_K)->GetOrigValues(N_Sets, BaseFunct, N_QuadPts_L2, zeta_L2,  LineQuadFormula,  Needs2ndDer);
    origvaluesD0=TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D0);
    L2DOF = L2GlobalNumbers + L2BeginIndex[i2];


    for(ii=0; ii<N_Xpos; ii++)
    {
     Sol_L1quadL0quadL2nodal = Sol_XposL1quadL0quadL2nodal  + ii*N_QuadPts_L1*N_QuadPts*N_LnodalPos[2];
     IntValue_L0quadL1quad = IntValue_XposL0quadL1quad +  ii*N_QuadPts_L1*N_QuadPts;
     Coords[0] =Xpos[2*ii];
     Coords[1] =Xpos[2*ii+1];

     for(j=0; j<N_QuadPts_L1; j++)
      {
       Sol_L0quadL2nodal = Sol_L1quadL0quadL2nodal +j*N_QuadPts*N_LnodalPos[2];
       Coords[3] = X1[j];  
       for(k=0; k<N_QuadPts; k++)
       {
        Coords[5] = X0[k];  
        Sol_L2 = Sol_L0quadL2nodal +k*N_LnodalPos[2];
        value = 0;
        for(m=0; m<N_QuadPts_L2; m++)
        {
         Coords[2] = X2[m];     
         orgD0 = origvaluesD0[m];
         val = 0.;
         for(n=0;n<N_BaseFunct;n++)
         {   
          val += Sol_L2[L2DOF[n]]*orgD0[n]; 
         } // n
         
        //GetExactSol
        Exact(N_Coord, Coords, IVal);
                  //  cout<< " L2: " << Coords[2] <<endl;
        value += (val - IVal[0])*(val - IVal[0])*Weights_L2[m]*AbsDetjk_L2[m]; //L_2 error
        } // m
        //L2 cell contribution to each L0 & L1 Quad points
        IntValue_L0quadL1quad[k*N_QuadPts_L1 + j] +=  value;

        // cout << "value  " << value << endl;
       } // k
     }// for j
    } // for ii 
   } // i2
//  cout<<endl;
   //int over L1
   for(ii=0; ii<N_Xpos; ii++)
   {
    IntValue_L0quadL1quad = IntValue_XposL0quadL1quad +  ii*N_QuadPts_L1*N_QuadPts;
    IntValue_L0quad = IntValue_XposL0quad +  ii*N_QuadPts;

    for(j=0; j<N_QuadPts; j++)
     {
      IntValue_L1quad = IntValue_L0quadL1quad + j*N_QuadPts_L1;
      value = 0.;  
      for(k=0; k<N_QuadPts_L1; k++)
      {
       value += IntValue_L1quad[k]*Weights_L1[k]*AbsDetjk_L1[k];
          // value +=  Weights_L1[k]*AbsDetjk_L1[k];      
          // cout<< "value " << IntValue_L1quad[k] << endl;
      }    
     IntValue_L0quad[j] += value;
     }
   }
  } // i1

   //int over L0
   for(ii=0; ii<N_Xpos; ii++)
   {
    IntValue_L0quad = IntValue_XposL0quad +  ii*N_QuadPts;
    value = 0.;      
    for(j=0; j<N_QuadPts; j++)
      value += IntValue_L0quad[j]*Weights[j]*AbsDetjk[j];
      // value += Weights[j]*AbsDetjk[j];
    L2Errors[ii] += value;
   }

 }// i

   for(j=0; j<N_Xpos; j++)
     cout<< j << " L2Errors " << L2Errors[j] <<endl;

// exit(0);


  // delete [] Sol_XposL2nodalL0quadL1nodal;
  // delete [] Sol_XposL1quadL0quadL2nodal;
  // delete [] IntValue_XposL0quadL1quad;


  return 0;

  // delete [] Sol_XposL1quadL0quadL2quad;
} // GetErrorAllXpos()



TSystemADI1D_3L::~TSystemADI1D_3L()
{

}