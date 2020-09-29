/** ************************************************************************ 
*
* @class     TAssemble
* @brief     base class for assembling matrices 
* @author    Sashikumaar Ganesan 
* @date      30.05.15
* @History    
 ************************************************************************  */


#ifndef __ASSEMBLE__
#define __ASSEMBLE__

#include <AllClasses.h>
#include <Constants.h>

#ifdef _MPI

#endif

#ifdef _OMPONLY

#endif

/** class for 3D scalar system matrix */
class TAssemble
{
  protected:  
    
    TCollection *Coll;
      
    int N_Cells;
    
    /** Number of Matrices */
    int N_SqMatrices, N_Matrices, N_AllMatrices;
    
    /** Number of FE spaces */
    int N_FeSpaces;    
    
    /** Number of RHS */    
    int N_Rhs;
    
    /** Variable for methods */
    int N_Parameters, **GlobalNumbers, **BeginIndex, **RhsGlobalNumbers, **RhsBeginIndex;
    int **TestGlobalNumbers, **TestBeginIndex, **AnsatzGlobalNumbers, **AnsatzBeginIndex;
    double **HangingEntries, **HangingRhs;
    double **Rhs, **Matrix, ***LocMatrices, **LocRhs, *aux, *auxarray, *rhsaux, *paramaux;
    double *auxmat, **Matrices;
    
    bool *SecondDer;
    
  public:
    /** Constructor*/
    TAssemble(int n_fespaces, int n_sqmatrices, int n_matrices, int n_rhs, double **rhs);
    
    /** destrcutor */
    ~TAssemble();
    
};

#endif
