#include <AllClasses.h>
#include <Constants.h>

#ifdef __3D__
  #include <Aux2D3D.h>
#endif

#ifdef __2D__
void Assemble2D_edge_Oseen(CoeffFct2D *Coeff, int n_fespaces, TFESpace2D **fespaces,
                      int n_sqmatrices, TSquareMatrix2D **sqmatrices,
                      int n_matrices, TMatrix2D **matrices,
                      int n_rhs, double **rhs, TFESpace2D **ferhs,
                      BoundCondFunct2D **BoundaryConditions,
                      BoundValueFunct2D **BoundaryValues,
                      TAuxParam2D *Parameters);
 
#endif