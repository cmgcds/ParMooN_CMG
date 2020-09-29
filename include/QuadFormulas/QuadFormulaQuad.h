// =======================================================================
// @(#)QuadFormulaQuad.h        1.3 05/04/99
//
// Class:      TQuadFormulaQuad
// Superclass: TQuadFormula2D
//
// Purpose:    quadrature formula for a 2D integral on the unit sqare
// Author:     Gunar Matthies
//
// History:    30.08.1997 start implementation
// 
// =======================================================================

#ifndef __QUAD_FORMULA_QUAD__
#define __QUAD_FORMULA_QUAD__

#include <QuadFormula2D.h>
#include <MooNMD_Io.h>

/** quadrature formula for a 2D integral */
class TQuadFormulaQuad : public TQuadFormula2D
{
  public:
    /** constructor */
    TQuadFormulaQuad();
    /** constructor */
    TQuadFormulaQuad(int n_points, double* weights, double* xi, 
                      double* eta, int acc);

    /** Gauss9x9 */
    void Gauss9();

    /** Gauss8x8 */
    void Gauss8();

    /** Gauss7x7 */
    void Gauss7();

    /** Gauss6x6 */
    void Gauss6();

    /** Gauss5x5 */
    void Gauss5();

    /** Gauss4x4 */
    void Gauss4();

    /** Gauss3x3 */
    void Gauss3();

    /** Gauss2x2 */
    void Gauss2();

    /** Vertex */
    void Vertex();

    /** Simpson */
    void Simpson();

};

#endif
