// =======================================================================
// @(#)ADICell1D.h        4.1 07.11.09
// 
// Class:       TADICell1D
// Purpose:     class for  ADICell1D

//
// Author:      Sashikumaar Ganesan (07.11.09)
//
// History:     start of implementation 07.11.09 (Sashikumaar Ganesan)
//
// =======================================================================

#ifndef __ADICELL1D__
#define __ADICELL1D__

#include <ADICell.h>
#include <FESpace1D.h>
#include <FEFunction1D.h>

class TADICell1D : public TADICell
{
  protected:

    /** x value of each quadPt */
    double *X;

    /** internal layer position values */
    double *GridX_internal;

    /** Fe space for all QuadPt*/
    TFESpace1D *FESpace1D_Internal;

    /** FE function for each QuadPt*/
    TFEFunction1D **FEFunction1D_Internal;

    /** Grid Fe function for all QuadPt*/
    TFEFunction1D *GridFEFunction1D_Internal;

    /** mass matrices for all QuadPt will not change */
    TSquareMatrix1D *M_Internal;

    /** stiffness matrix for all QuadPts*/
    TSquareMatrix1D *A_Internal;

    /** Initial condition*/
     DoubleFunct2D *Initial;

    /** Initial condition*/
     DoubleFunct2D *Exact;



  private:
   int ConstructAllInfo();

  public:
    /** constructor */
    TADICell1D(TBaseCell *cell,  TFESpace1D *FESpace1D_internal, TSquareMatrix1D *M_internal, TSquareMatrix1D *A_internal,
               int N_quadPts, double *x, DoubleFunct2D *initial,
               DoubleFunct2D *exact);


    void SolveAllQdPts(double *G, double *QuadPtsRhsT, CoeffFct2D *BilinearCoeffs,
                       BoundCondFunct2D *BoundaryCondition, BoundValueFunct2D *BoundValue,
                       double tau, double *Sol_Loc);

    void AssembleARhs(double x, double Conv, CoeffFct2D *Bilinear, BoundCondFunct2D *BoundaryCondition, BoundValueFunct2D *BoundValue);

   int GetN_InternalLevels()
    {
     return N_V;
    }

   double *Get_Xpos()
    {
     return GridX_internal;
    }

    /** destrcutor */
    ~TADICell1D();

};

#endif
