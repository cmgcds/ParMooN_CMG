// =======================================================================
// @(#)ADISystem.h        4.1 13.11.09
// 
// Class:       TADISYSTEM
// Purpose:     general super class for all ADI/operator-splitting System
//              special spaces are implemented in subclasses
//
// Author:      Sashikumaar Ganesan (13.11.09)
//
// History:     start of implementation 13.11.09 (Sashikumaar Ganesan)
//
// =======================================================================

#ifndef __ADISYSTEM__
#define __ADISYSTEM__

#include <BaseCell.h>
#include <Collection.h>

/** general super class for ADI/operator splitting scheme, special spaces are
    implemented in subclasses */
class TADISystem
{
  protected:
// =======================================================================
// // information of cell and the dimension of the internal domain
// =======================================================================
    /** all internal domains use the same coll */
    TCollection *Collection_Intl;

    /** number of degrees of freedom  (same for all QuadPt) */
    int N_V;

    /** solution vector: if we have, we no need to allocate & deallocate for each Intl poin in internal direction */
    double *sol;

    /** solution vector: if we have, we no need to allocate & deallocate for each Intl poin in internal direction */
    double *oldsol;

    /** rhs vector  for each QuadPt in internal direction [N_V]*/
    double *rhs;

    /** working rhs: if we have, we no need to allocate & deallocate for each Intl poin in internal direction */
    double *B;

    /** tmp array, if we have, we no need to allocate & deallocate for each Intl point*/
    double *defect;

    /** No. M mat value  */
    int N_MMatValues;

    /** DoubleFunct2D to calculate growth rate*/
    DoubleFunctND *GetGrowthAndNuc;

    /** store l values */
    double *IntlPosL;

    /** No. A mat value  */
//     int N_AMatValues;

    /** store the SUPG S-mat value for next time step rhs calculation, for other than Euler schems, time-dep. growth */
//     double *SMatValues;

    /** No.  SUPG S-mat values*/
//     int N_SMatValues;

    /** store the SUPG K-mat value for next time step rhs calculation, for other than Euler schems, time-dep. growth */
//     double *KMatValues;

    /** No. SUPG K-mat values*/
//     int N_KMatValues;

  public:
    /** constructor */
    TADISystem(double *Sol, double *OldSol, double *b, double *Defect, DoubleFunctND *growthfunct);

    /** destrcutor */
    ~TADISystem();

};

#endif
