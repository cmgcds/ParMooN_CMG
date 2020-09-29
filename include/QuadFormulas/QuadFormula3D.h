// =======================================================================
// @(#)QuadFormula3D.h        1.2 05/04/99
//
// Class:      TQuadFormula3D
// Superclass: TQuadFormula
//
// Purpose:    quadrature formula for a 3D integral
// Author:     Gunar Matthies
//
// History:    30.08.1997 start implementation
// 
// =======================================================================

#ifndef __QUAD_FORMULA_3D__
#define __QUAD_FORMULA_3D__

#include <Enumerations.h>
#include <QuadFormula.h>
#include <MooNMD_Io.h>

/** quadrature formula for a 3D integral */
class TQuadFormula3D : public TQuadFormula
{
  protected:
    /** first coordinate for the formula */
    double *Xi;
    /** second coordinate for the formula */
    double *Eta;
    /** third coordinate for the formula */
    double *Zeta;

  protected:
    /** This is a private method for initializing the data structure */
    void InitObject(int n, double* w, double* xi, 
                    double* eta, double* zeta, int acc);

  public:
    /** constructor */
    TQuadFormula3D();
    /** constructor */
    TQuadFormula3D(int n_points, double* weights, 
      double* xi, double* eta, double* zeta, int acc);

    /** return coordinates of the formula */
    virtual double *GetCoords(int i);
    /** return all data of the quadrature formula */
    void GetFormulaData(int &n_points, double* &weights, 
                        double* &xi, double* &eta, double* &zeta);

#ifdef __3D__
    /** find a quadrature formula for all given elements */
    static void FindLocalQuadFormula3D(int N_LocalUsedElements,
        FE3D *LocalUsedElements,
        QuadFormula2D &qf1, QuadFormula3D &qf2);
#endif // __3D__

    /** print all information of this formula */
    friend std::ostream & operator << (std::ostream &s, TQuadFormula3D *qf);
};

#endif
