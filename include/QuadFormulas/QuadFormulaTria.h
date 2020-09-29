// =======================================================================
// @(#)QuadFormulaTria.h        1.3 12/08/99
//
// Class:      TQuadFormulaTria
// Superclass: TQuadFormula2D
//
// Purpose:    quadrature formula for a 1D integral
// Author:     Gunar Matthies
//
// History:    29.08.1997 start implementation
// 
// =======================================================================

#ifndef __QUAD_FORMULA_TRIA__
#define __QUAD_FORMULA_TRIA__

#include <QuadFormula2D.h>
#include <MooNMD_Io.h>

/** quadrature formula for a 2D integral */
class TQuadFormulaTria : public TQuadFormula2D
{
  public:
    /** constructor */
    TQuadFormulaTria();
    /** constructor */
    TQuadFormulaTria(int n_points, double* weights, double* xi, 
                      double* eta, int acc);

    /** BaryCenter */
    void BaryCenter();
    /** MidPoint */
    void MidPoint();
    /** SevenPoint */
    void SevenPoint();
    /** high order formula, nearly a Gauss3 formula */ 
    void Gauss3();
    /** Vertex */
    void Vertex();
    /** formula of degree 8 */
    void Degree8();
      /** formula of degree 9 */  
    void Degree9();
    /** formula of degree 11 */
    void Degree11();
    /** formula of degree 19 */
    void Degree19();
    /** Gauss-like formula, composed on three sub-triangles */
    void CompGauss3();
    /** Gauss-like formula, composed on four sub-triangles */
    void CompGauss4();
    /** Gauss formula of degree 8 */
    void Gauss_Degree8();
};

#endif
