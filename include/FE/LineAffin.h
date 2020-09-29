// =======================================================================
// @(#)LineAffin.h
//
// Class:      LineAffin
//
// Purpose:    reference transformations for Line
//
// Author:     Sashikumaar Ganesan
//
// History:    17.05.2007 start implementation
// 
// =======================================================================

#ifndef __LINEAFFIN__
#define __LINEAFFIN__

#include <Enumerations.h>
#include <RefTrans1D.h>

/** reference transformations for line */
class TLineAffin : public TRefTrans1D
{
  protected:
    /** x coordinate */
    double x0, x1;

    /** x coordinate (useful in surface finite element) */
    double y0, y1;

    /** x parameters for reference transformation */
    double xc0, xc1;

    /** x parameters for reference transformation */
    double yc0, yc1;

    /** detjk */
    double detjk;

    /** 1/detjk */
    double rec_detjk;

  public:
    /** constuctor */
    TLineAffin();

    /** transfer form reference element to original element */
    void GetOrigFromRef(double xi, double &X
#ifdef __2D__
, double &Y
#endif
                                );

    /** transfer a set of points form reference to original element */
    void GetOrigFromRef(int N_Points, double *xi, double *X, double *Y, double *absdetjk);
    /** transfer a set of points form reference to original element */
    void GetOrigFromRef(int N_Points, double *xi, double *X, double *absdetjk);

    void GetOrigValues(int N_Sets, BaseFunct1D *BaseFuncts, int N_Points, double *zeta,
                       QuadFormula1D QuadFormula, bool *Needs2ndDer);

//     /** set element to cell */
    void SetCell(TBaseCell * cell);
    
    double Getrec_detjk()
    { return rec_detjk; }
};

#endif
 
