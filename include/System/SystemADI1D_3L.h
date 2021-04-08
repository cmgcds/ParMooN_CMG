// =======================================================================
// @(#)ADISystem1D.h        4.1 13.04.20
// 
// Class:       TSystem1DADI
// Purpose:     class for  SystemADI1D_3L

//
// Author:      Sashikumaar Ganesan (13.04.20)
//
// History:     start of implementation 13.04.20 (Sashikumaar Ganesan)
//
// =======================================================================

#ifndef __SYSTEMADI1D_3L__
#define __SYSTEMADI1D_3L__

#include <SystemADI1D.h>
#include <FESpace1D.h>
#include <FEFunction1D.h>
#include <SquareMatrix1D.h>

class TSystemADI1D_3L : TSystemADI1D
{
  protected:
  TSystemADI1D **SystemADI;

  TFEFunction2D *SpatialFeFunction;

  int N_XposLnodal_All, N_L1L0, N_L2L1, N_Lnodal_All, *N_LDof, N_L0ADIs,  N_L1ADIs, N_L2ADIs;

  double *Sol_XnodalLnodal[3], *Sol_XnodalLDof[3], *Sol_XdofLdof, *Sol_XposL2L0L1nodal, *Sol_XposL0L1L2nodal;

  /** variables for LInt */
  double X0[MaxN_QuadPoints_1D], AbsDetjk[MaxN_QuadPoints_1D];
  double NodalX0[MaxN_QuadPoints_1D], NodalAbsDetjk[MaxN_QuadPoints_1D]; 
  double X1[MaxN_QuadPoints_1D],  X2[MaxN_QuadPoints_1D];
  double AbsDetjk_L1[MaxN_QuadPoints_1D], AbsDetjk_L2[MaxN_QuadPoints_1D]; 
  double *Sol_XposL2nodalL0quadL1nodal, *IntValue_XposL0quad , *Sol_XposL1quadL0quadL2nodal;
  double *IntValue_XposL0quadL1quad, *NucValue_XposL0quadL1quad, *NucValue_XposL0quad;  

  private:


  public:
    /** constructor */
    
    TSystemADI1D_3L(TSystemADI1D **systemADI, int *l_NVertices, double *l_StartEnd, BoundCond1D **boundConLminLMax, 
                    DoubleFunctND **growthAndB_Nuc, TFEFunction2D *spatialFeFunction);

    //initialize the system
    void Init(int n_Xpos, double *xpos, int *n_LDof, int *n_LnodalPos, double **LnodalPos, double *sol_XdoLdof);

    void Int_B_Nuc(double *IntValue, double *B_NucValue, DoubleFunctND *GrowthAndB_Nuc, double *SuscepPopRatio);

    void Solve(int ADI_SolveIdx, CoeffFctND *Bilinear, double *IntL_Values, double *SuscepPopRatio);

    void Interpolate(TCollection *coll, double *InitialPopulation, DoubleFunctND *Exact);

    void Interpolate(DoubleFunctND *Exact);

    void XnodalLnodalToLdofXdof(int MaxN_PtsForNodal, double *Sol_LnodalXdof);

    int Nodal2DOF(int ADI_Idx, double *Sol_XposLNnodLOwnDof);

    int DOF2Nodal(int ADI_Idx, double *Sol_XposLNnodLOwnDof);

    int CopySolToInternal(int ADI_Idx);

    int CopySolFromInternal(int ADI_Idx);

    void LdofXdofToXnodalLdof(double *Sol_LnodalXdof);
    
    double IntL(DoubleFunctND *GetGamma_Q, double *IntValue, double *NucValue, double *LDistXSum=nullptr);

    double GetErrorAllXpos(DoubleFunctND *Exact, double *L2Errors);
    /** destrcutor */
    ~TSystemADI1D_3L();

};

#endif
