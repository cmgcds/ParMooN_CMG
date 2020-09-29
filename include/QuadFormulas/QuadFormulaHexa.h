// =======================================================================
// @(#)QuadFormulaHexa.h        1.2 05/04/99
//
// Class:      TQuadFormulaHexa
// Superclass: TQuadFormula3D
//
// Purpose:    quadrature formula for a 3D integral on a hexahedron
// Author:     Gunar Matthies
//
// History:    30.08.1997 start implementation
// 
// =======================================================================

#ifndef __QUAD_FORMULA_HEXA__
#define __QUAD_FORMULA_HEXA__

#include <QuadFormula3D.h>
#include <MooNMD_Io.h>

/** quadrature formula for a 3D integral on a hexahedron */
class TQuadFormulaHexa : public TQuadFormula3D
{
  public:
    /** constructor */
    TQuadFormulaHexa();
    /** constructor */
    TQuadFormulaHexa(int n_points, double* weights, 
      double* xi, double* eta, double* zeta, int acc);

    /** Vertex */
    void Vertex();
    /** Gauss2x2x2 */
    void Gauss2();

    /** Gauss 3x3x3 */
    void Gauss3();

    /** Gauss 4x4x4 */
    void Gauss4();

    /** Gauss 5x5x5 */
    void Gauss5();

    /** Gauss 6x6x6 */
    void Gauss6();

    /** Gauss 7x7x7 */
    void Gauss7();

    /** Gauss 8x8x8 */
    void Gauss8();

    /** Gauss 9x9x9 */
    void Gauss9();

    /** quad rule with vertices and origin */
    void VerticesAndOrigin();

    /** quad rule with vertices and origin, 15 quad points */
    void VerticesAndOrigin15();
 
    /** quad rule with vertices and origin, 57 quad points */
    void VerticesAndOrigin57();

    /** quad rule for polynomials of degree 7 with 38 quad points */
    void Degree7_Points38();
};

#endif
