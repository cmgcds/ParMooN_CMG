// =======================================================================
// @(#)QuadFormulaTetra.h        1.2 05/04/99
//
// Class:      TQuadFormulaTetra
// Superclass: TQuadFormula3D
//
// Purpose:    quadrature formula for a 3D integral on a tetrahedron
// Author:     Gunar Matthies
//
// History:    30.08.1997 start implementation
// 
// =======================================================================

#ifndef __QUAD_FORMULA_TETRA__
#define __QUAD_FORMULA_TETRA__

#include <QuadFormula3D.h>
#include <MooNMD_Io.h>

/** quadrature formula for a 3D integral on a tetrahedron */
class TQuadFormulaTetra : public TQuadFormula3D
{
  public:
    /** constructor */
    TQuadFormulaTetra();
    /** constructor */
    TQuadFormulaTetra(int n_points, double* weights, 
      double* xi, double* eta, double* zeta, int acc);

    /** BaryCenter */
    void BaryCenter();

    /** Vertex */
    void Vertex();

    /** P2 exact formula */
    void P2Exact();

    /** P4 exact formula */
    void P4Exact();

    /** P5 exact formula */
    void P5Exact();

    /** P8 exact formula */
    void P8Exact();
};

#endif
