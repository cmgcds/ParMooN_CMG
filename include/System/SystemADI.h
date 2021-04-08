// =======================================================================
// @(#)ADISystem.h        4.1 13.04.20
// 
// Class:       TSYSTEMADI
// Purpose:     general super class for all ADI/operator-splitting System
//              special spaces are implemented in subclasses
//
// Author:      Sashikumaar Ganesan (13.04.20)
//
// History:     start of implementation 13.04.20 (Sashikumaar Ganesan)
//
// =======================================================================

#ifndef __SYSTEMADI__
#define __SYSTEMADI__

#include <Domain.h>
#include <BaseCell.h>
#include <Collection.h>

/** general super class for ADI/operator splitting scheme, special spaces are
    implemented in subclasses */
class TSystemADI
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

    /** number of degrees of freedom    */
    int N_Dof;

    /** number of nodal points to evaluate DOFs    */
    int N_NodalPts;

    /** nodal points coordinates  */
    double *NodalPt_Coord;

    /** solution vector: */
    // double *sol;
    double *Nodal_sol;

    /** solution vector:  */
    // double *oldsol;

    /** rhs vector for each QuadPt in internal direction [N_V]*/
    double *rhs;

    /** working rhs: if we have, we no need to allocate & deallocate for each Intl poin in internal direction */
    double *B;

    /** tmp array, if we have, we no need to allocate & deallocate for each Intl point*/
    // double *defect;

    /** No. M mat value  */
    int N_MatEntries;
    double *Mvalues_orig, *Advectvalues_orig;


    /** Boundary values*/
    // DoubleFunctND  *BDValue_LMin;
    // DoubleFunctND  *BDVal_LMax;

    /** DoubleFunctND to calculate growth rate*/
    DoubleFunctND *GetGrowthAndNuc;

    /** store l values */
    double **LnodalPos;

    /** other coordinates (physical/internal)  values */
    int N_Coord;

    /** other coordinates (physical/internal)  values */
    int OwnADI_Idx, N_Xpos, N_ADISystems, *N_LnodalPos, *IncidentArray; 
    double *Xpos; // including own L pos
    double *Sol_XposLNnodal;

     boolean DofSameAsNodalPts;
  public:
    /** constructor */
    TSystemADI();

    TSystemADI(DoubleFunctND *growthAndNuc);

    int GetN_Cells()
    {return Coll_Intl->GetN_Cells(); }

    int GetN_Dof()
    {return N_Dof; }

    int GetN_NodalPts()
    {return N_NodalPts; }

    double *GetNodalPt_Coord()
    {return NodalPt_Coord; }

    bool IsdGdisc()
    { return dGDisc; }

    /** destrcutor */
    ~TSystemADI();

};

#endif
