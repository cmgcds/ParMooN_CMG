// =======================================================================
// @(#)ADICell.h        4.1 07.11.09
// 
// Class:       TADICell
// Purpose:     general super class for all ADICells
//              special spaces are implemented in subclasses
//
// Author:      Sashikumaar Ganesan (07.11.09)
//
// History:     start of implementation 07.11.09 (Sashikumaar Ganesan)
//
// =======================================================================

#ifndef __ADICELL__
#define __ADICELL__

#include <BaseCell.h>
#include <Collection.h>

/** general super class for all ADICells, special spaces are
    implemented in subclasses */
class TADICell
{
  protected:
// =======================================================================
// // information of cell and the dimension of the internal domain
// =======================================================================
    /** cells for which the info of internal direction is constructed*/
    TBaseCell *Cell;

    /** all internal domains use the same coll */
    TCollection *Collection_Internal;

    /** number of quadrature points in this cell, where the internal domains have to be constructed*/
    int N_QuadPts;

    /** number of degrees of freedom  (same for all QuadPt) */
    int N_V;

    /** solution vector  for each QuadPt in internal direction [N_QuadPts][N_V]*/
    double *sol;

    /** solution vector in internal direction for each QuadPt [N_V][N_QuadPts]*/
    double *solT;

    /** rhs vector  for each QuadPt in internal direction [N_QuadPts][N_V]*/
    double *rhs;

    /** rhs vector in internal direction for each QuadPt [N_V][N_QuadPts] */
    double *rhsT;

  public:
    /** constructor */
    TADICell(TBaseCell *cell, int N_quadPts);

    /** destrcutor */
    ~TADICell();

};

#endif
