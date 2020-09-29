// =======================================================================
// @(#)QuadFormula2D.h        1.2 05/04/99
//
// Class:      TQuadFormula2D
// Superclass: TQuadFormula
//
// Purpose:    quadrature formula for a 2D integral
// Author:     Gunar Matthies
//
// History:    29.08.1997 start implementation
// 
// =======================================================================

#ifndef __QUAD_FORMULA_2D__
#define __QUAD_FORMULA_2D__

#include <Enumerations.h>
#include <QuadFormula.h>
#include <MooNMD_Io.h>

/** quadrature formula for a 2D integral */
class TQuadFormula2D : public TQuadFormula
{
  protected:
    /** first coordinate in [0,1]x[0,1] for the formula */
    double *Xi;
    /** second coordinate in [0,1]x[0,1] for the formula */
    double *Eta;

  protected:
    /** This is a private method for initializing the data structure */
    void InitObject(int n, double* w, double* xi, double* eta, int acc);

  public:
    /** constructor */
    TQuadFormula2D();
    /** constructor */
    TQuadFormula2D(int n_points, double* weights, double* xi, 
                    double* eta, int acc);

    /** return coordinates of the formula */
    virtual double *GetCoords(int i);
    /** return all data of the quadrature formula */
    void GetFormulaData(int &n_points, double* &weights, 
                        double* &xi, double* &eta);

// #ifdef __2D__
    /** return a quadrature formula which can be used for 
        all given elements */
    static void FindQuadFormula2D(FE2D *UsedElements,
        QuadFormula1D &qf1, QuadFormula2D &qf2);

    /** return a quadrature formula which can be used for
        the given elements */
    static void FindQF_2D(FE2D CurrentElement,
        QuadFormula1D &qf1, QuadFormula2D &qf2);

    /** find a quadrature formula for all given elements */
    static void FindLocalQuadFormula2D(int N_LocalUsedElements,
        FE2D *LocalUsedElements,
        QuadFormula1D &qf1, QuadFormula2D &qf2);
// #endif // __2D__

    /** print all information of this formula */
    friend std::ostream & operator << (std::ostream &s, TQuadFormula2D *qf);
};

#endif
