// =======================================================================
// @(#)ADISystem1D.h        4.1 13.04.20
// 
// Class:       TSystem1DADI
// Purpose:     class for  SystemADI1D

//
// Author:      Sashikumaar Ganesan (13.04.20)
//
// History:     start of implementation 13.04.20 (Sashikumaar Ganesan)
//
// =======================================================================

#ifndef __SYSTEMADI1D__
#define __SYSTEMADI1D__

#include <SystemADI.h>
#include <FESpace1D.h>
#include <FEFunction1D.h>
#include <SquareMatrix1D.h>

class TSystemADI1D : public TSystemADI
{
  protected:
 
    /** Fe space of configuration space */
    TFESpace1D *FESpace1D_Intl;

    /** fe function, needed ofr interpolation */
    TFEFunction1D *FENodalFunction;

    /** mass matrices for all QuadPt will not change */
    TSquareMatrix1D *M_Intl;

    /** stiffness matrix for all QuadPts*/
    TSquareMatrix1D *A_Intl;

    /** supg mass matrices for all QuadPt will not change */
    TSquareMatrix1D *S_Intl;

    /** supg stiffness matrix for all QuadPts*/
    TSquareMatrix1D *K_Intl;

    BoundCond cond_Lmin, cond_Lmax;

  private:
   int ConstructAllInfo();

   void GenerateNodalPts();
   
   void Generate1DMesh(double Start, double End, int N_Cells);

  public:
    /** constructor */
    TSystemADI1D();

    TSystemADI1D(int N_L, double start, double end, BoundCond1D *boundConLminLMax, DoubleFunctND *growthAndNuc);


    //initialize the system
    void Init(int n_Coord, int n_Xpos, double *xpos, int n_ADISystems, int *n_LnodalPos, double **LnodalPos, int ownADI_Idx, double *sol_XposLNnodal);

    void AssembleMassMat();

    void AssembleAdvectMat(bool Conservative);

    void Solve(int N_Param, double *Coords, CoeffFctND *Bilinear, double *Sol);

    void Interpolate(int n_Coord, double *Coords, double *Sol, DoubleFunctND *Exact);

    int Nodal2DOF(double *Sol_XposLNnodLOwnDof);

    int DOF2Nodal(double *Sol_XposLNnodLOwnDof);
    
    void AssemblDG_Rhs(double Bnuc);

    //void Interpolate(double *Sol, DoubleFunct3D *Exact);
    void SetDirichletBc();
     
    void AssembleARhs_FD(double Conv);
	
    void AssembleARhs(double *X,double Conv, CoeffFctND *Bilinear);
    
    void AssembleARhs_SUPG(double *X,double Conv, CoeffFctND *Bilinear);

    void AssembleARhs_DG(double *X, double Conv, double *Out, CoeffFctND *Bilinear);
    
    void AssembleARhs_DG_Advect(double *X, double *NucGrowth, CoeffFctND *Bilinear);

    // double GetWeightedF();

  TFESpace1D *GetFeSpace1D()
   {
    return FESpace1D_Intl;
   }

   int GetN_InternalLevels()
    {
     return N_Dof;
    }

    // double GetQ3Max(double *currsol);


    /** destrcutor */
    ~TSystemADI1D();

};

#endif
