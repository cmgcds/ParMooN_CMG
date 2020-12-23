// =======================================================================
// @(#)System1D.h        4.1 13.04.20
// 
// Class:       TSYSTEM1D
// Purpose:     general super class for all 1D System
//
// Author:      Sashikumaar Ganesan (12.12.20)
//
// History:     start of implementation 12.12.20 (Sashikumaar Ganesan)
//
// =======================================================================

#ifndef __SYSTEM1D__
#define __SYSTEM1D__

#include <Domain.h>
#include <BaseCell.h>
#include <Collection.h>
#include <FESpace1D.h>
#include <FEFunction1D.h>
#include <SquareMatrix1D.h>

/** general super class for 1D */
class TSystem1D
{
  protected:
// =======================================================================
// // information of cell and the dimension of the internal domain
// =======================================================================
    /** Domain for the internal space */
    TDomain *Domain_Intl;

    /** all internal domains use the same coll */
    TCollection *Coll_Intl;

    /** dG space or not */
    bool dGDisc;

    /** Fe space of configuration space */
    TFESpace1D *FESpace1D;

    /** fe function, needed for interpolation */
    TFEFunction1D *FE_Function;

    /**  Fe Space for computing the mean value of the Cell  **/
    TFESpace1D *FESpaceMean1D;

    /** Fe functions to store the mean Values in the domain **/
    TFEFunction1D *MeanSolution , *MeanSolutionGradient;
  
    /** FEVectFunction to store both the mean value of solution and Gradients **/
    TFEVectFunct1D *MeanFEVectFunction;

    /** data for Bilinear  */
    CoeffFctND *Bilinear;

    /** Matrix atructure */
    TSquareStructure1D *SqStructure;

    /** mass matrices for all QuadPt will not change */
    TSquareMatrix1D *M_Intl;

    /** stiffness matrix for all QuadPts*/
    TSquareMatrix1D *A_Intl;

    /** supg mass matrices for all QuadPt will not change */
    TSquareMatrix1D *S_Intl;

    /** supg stiffness matrix for all QuadPts*/
    TSquareMatrix1D *K_Intl;

    BoundCond cond_Lmin, cond_Lmax;

    /** number of degrees of freedom    */
    int N_Dof;

    /** rhs vector for each QuadPt in internal direction [N_V]*/
    double *Sol, *Rhs;

    /** Arrays to Store the mean values of Solution and gradient */
    double *meanSolution, *meanSolutionGrad;

    /** working rhs: if we have, we no need to allocate & deallocate for each Intl poin in internal direction */
    double *defect, *B;

    /** No. M mat value  */
    int N_MatEntries;
    double *Mvalues_orig, *Advectvalues_orig;

    /** DoubleFunctND to calculate BundValues */
    DoubleFunctND *BundValues;

    /** other coordinates (physical/internal)  values */
    int N_Coord;

    void Generate1DMesh(double Start, double End, int N_Cells);
    
  public:
    /** constructor */  
    TSystem1D(int N_L, double start, double end, BoundCond1D *boundConLminLMax, DoubleFunctND *BdValues, char *ParamFile);

    int GetN_Cells()
    {return Coll_Intl->GetN_Cells(); }

    int GetN_Dof()
    {return N_Dof; }


    bool IsdGdisc()
    { return dGDisc; }

    /** destrcutor */
    ~TSystem1D();

};

#endif
