// =======================================================================
// @(#)FEVectFunct2D.h        1.2 07/20/99
// 
// Class:       TFEVectFunct2D
// Purpose:     a function from a finite element space in 2D
//
// Author:      Gunar Matthies (17.01.98)
//
// History:     start of implementation 17.01.98 (Gunar Matthies)
//
//              start of reimplementation 06.08.1998 (GM)
//
// =======================================================================


#ifndef __FEVECTFUNCT2D__
#define __FEVECTFUNCT2D__


#include <FEFunction2D.h>

/** a function from a finite element space */
class TFEVectFunct2D : public TFEFunction2D
{
  protected:
    /** number of components */
    int N_Components;

  public:
    /** constructor with vector initialization */
    TFEVectFunct2D(TFESpace2D *fespace2D, char *name, char *description,
                  double *values, int length, int n_components);

    /** return number of components */
    int GetN_Components()
    { return N_Components; }

    /** return i-th component as FEFunction2D */
    TFEFunction2D *GetComponent(int i)
    {
      // the name of the component will include the index i
      std::ostringstream os;
      os.seekp(std::ios::beg);
      os << Name << i <<ends;
      return new TFEFunction2D(FESpace2D, (char*)os.str().c_str(), Description,
                               Values+i*Length, Length);
    }

    /** convert current grid to vector-values FE function */
    void GridToData();

    /** use current data for grid replacement */
    void DataToGrid();

    /** calculate errors to given vector function */
    void GetDeformationTensorErrors( 
        DoubleFunct2D *Exact, DoubleFunct2D *Exact1,
        int N_Derivatives,
        MultiIndex2D *NeededDerivatives,
        int N_Errors, ErrorMethod2D *ErrorMeth, 
        CoeffFct2D *Coeff, TAuxParam2D *Aux,
        int n_fespaces, TFESpace2D **fespaces,
        double *errors);

    /** calculate L2-nrom of divergence */
    double GetL2NormDivergence();

  /** write the solution into a data file **/
   void WriteSol(double t);

  /** Read the solution from a given data file - written by Sashi **/
   void ReadSol(char *BaseName);
   
   /** determine the value of a vect function and its first derivatives at
    the given point */
   void FindVectGradient(double x, double y, double *val1, double *val2);

   /** interpolate the old vect value to the new function */
   void Interpolate(TFEVectFunct2D *OldVectFunct);
   
   /** @brief multiply function with a scalar alpha. Only non-Dirichlet dofs are 
    *         multiplied! */
    TFEVectFunct2D& operator*=(double alpha);
   
    /** @brief add one TFEVectFunct2D to another one. Both have to be defined on
     *         the same space. Only non-Dirichlet dofs are added!  */
    TFEVectFunct2D & operator+=(const TFEVectFunct2D & rhs);
   
    /** @brief copy one TFEVectFunct2D to another one. Both have to be defined 
     *         on the same space */
    TFEVectFunct2D & operator=(const TFEVectFunct2D & rhs);
};

#endif

// #ifdef __2D__
// # endif // #ifdef __2D__
