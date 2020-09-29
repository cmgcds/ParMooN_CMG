// =======================================================================
// @(#)BaseFunct1D.h        1.2 11/20/98
//
// Class:      TBaseFunct1D
//
// Purpose:    represents the set of base functions for a finite element
//             in one dimensions
//
// Author:     Gunar Matthies
//
// History:    08.10.98 restart implementation
// 
// =======================================================================

#ifndef __BASEFUNCT1D__
#define __BASEFUNCT1D__

#include <QuadFormula1D.h>
#include <Constants.h>

#include <Enumerations.h>

/** set of all base function on the reference element for a finite 
    element in two dimensions */
class TBaseFunct1D
{
  protected:
    /** number of base functions = dimension of local space */
    int Dimension;

    /** Id for this set of base functions */
    BaseFunct1D BaseFunct;

    /** array for all functions and derivatives */
    DoubleFunct1D *Functions[N_MultiIndices1D];

    /** status of changability of entries */
    bool changable;

    /** polynomial degree */
    int PolynomialDegree;

    /** accuracy */
    int Accuracy;

  public:
    /** constructor, fill in all information */
    TBaseFunct1D(int dimension, BaseFunct1D basefunct,
                 DoubleFunct1D* functions, 
                 DoubleFunct1D* derivativesxi,
                 DoubleFunct1D* derivativesxixi);

    /** constructor, fill in all information */
    TBaseFunct1D(int dimension, BaseFunct1D basefunct,
                 DoubleFunct1D* functions, 
                 DoubleFunct1D* derivativesxi,
                 DoubleFunct1D* derivativesxixi,
                 int polynomialdegree,
                 int accuracy);

    /** constructor without filling data structure */
    TBaseFunct1D(int dimension);

    /** return the dimension of local space */
    int GetDimension() 
    { return Dimension; }

    /** return the values for derivative MultiIndex at xi */
    void GetDerivatives(MultiIndex1D MultiIndex, double xi,
                        double *values)
      { Functions[MultiIndex](xi, values); };

    /** return the values for derivative MultiIndex at all
        quadrature points */
    void GetDerivatives(MultiIndex1D MultiIndex, 
                        TQuadFormula1D *formula, double **values);

    /** set status to unchangable */
    void SetUnchangable()
      { changable = false; };

    /** set function for derivative MultiIndex */
    void SetFunction(MultiIndex1D MultiIndex, DoubleFunct1D* function);

    /** make date on reference element */
    void MakeRefElementData(QuadFormula1D QuadFormula);

    /** return polynomial degree */
    int GetPolynomialDegree()
      { return PolynomialDegree; };

    /** return accuracy */
    int GetAccuracy()
      { return Accuracy; };

};

#endif
