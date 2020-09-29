// =======================================================================
// @(#)TriaAffin.h        1.5 04/13/00
//
// Class:      TTriaAffin
//
// Purpose:    reference transformations for triangle
//
// Author:     Gunar Matthies
//
// History:    08.07.97 start implementation
// 
// =======================================================================

#ifndef __TRIAAFFIN__
#define __TRIAAFFIN__

#include <Enumerations.h>
#include <RefTrans2D.h>

/** reference transformations for triangle */
class TTriaAffin : public TRefTrans2D
{
  protected:
    /** x coordinate */
    double x0, x1, x2;

    /** y coordinate */
    double y0, y1, y2;

    /** x parameters for reference transformation */
    double xc0, xc1, xc2;

    /** y parameters for reference transformation */
    double yc0, yc1, yc2;

    /** detjk */
    double detjk;

    /** 1/detjk */
    double rec_detjk;

  public:
    /** constuctor */
    TTriaAffin();

    /** transfer form reference element to original element */
    void GetOrigFromRef(double eta, double xi, double &x, double &y);

    /** transfer a set of points form reference to original element */
    void GetOrigFromRef(int N_Points, double *eta, double *xi, 
                        double *x, double *y, double *absdetjk);

    /** transfer form reference element to original element */
    void GetOrigFromRef(double *ref, double *orig);

    /** transfer from original element to reference element */
    void GetRefFromOrig(double x, double y, double &eta, double &xi);

    /** transfer from original element to reference element */
    void GetRefFromOrig(double *orig, double *ref);

    /** calculate functions and derivatives from reference element
        to original element */
    void GetOrigValues(BaseFunct2D BaseFunct,
                       int N_Points, double *xi, double *eta,
                       int N_Functs, QuadFormula2D QuadFormula);

    /** calculate functions and derivatives from reference element
        to original element, for all given elements */
    void GetOrigValues(int N_Sets, BaseFunct2D *BaseFunct,
                       int N_Points, double *xi, double *eta,
                       QuadFormula2D QuadFormula,
                       bool *Needs2ndDer);

    /** calculate functions and derivatives from reference element
        to original element */
    void GetOrigValues(double xi, double eta, int N_BaseFunct,
                       double *uref, double *uxiref, double *uetaref,
                       double *uorig, double *uxorig, double *uyorig,
                       int _BaseVectDim = 1);
    
    /** calculate functions and derivatives from reference element
        to original element */
    void GetOrigValues(int joint, double zeta, int N_BaseFunct,
                       double *uref, double *uxiref, double *uetaref,
                       double *uorig, double *uxorig, double *uyorig,
                       int _BaseVectDim=1);
    
    /** set element to cell */
    void SetCell(TBaseCell * cell);

    /** return outer normal vector */
    void GetOuterNormal(int j, double zeta,
                                double &n1, double &n2);

    /** return tangent */
    void GetTangent(int j, double zeta,
                                double &t1, double &t2);

    /** return volume of cell according to reference transformation */
    double GetVolume();
    
    /** return boundary vertices */
    void GetOrigBoundFromRef(int joint, int N_Points, double *zeta, double *X, double *Y);
    
    /** @brief Piola transformation for vector basis */
    void PiolaMapOrigFromRef(int N_Functs, double *refD00, double *origD00);
    void PiolaMapOrigFromRef(int N_Functs, double *refD10, double *refD01, 
                             double *origD10, double *origD01);

};

#endif
