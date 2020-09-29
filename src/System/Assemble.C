/** ************************************************************************ 
* @brief     source file for TAssemble
* @author    Sashikumaar Ganesan
* @date      30.05.15
* @History 
 ************************************************************************  */

#include <Assemble.h>



TAssemble::TAssemble(int n_fespaces, int n_sqmatrices, int n_matrices, int n_rhs, double **rhs)
{
  N_SqMatrices = n_sqmatrices;
  N_Matrices = n_matrices;
  N_FeSpaces = n_fespaces;    
  N_Rhs = n_rhs;
  N_AllMatrices = n_sqmatrices+n_matrices;
  
  
  if(N_SqMatrices)
   {  
    GlobalNumbers = new int* [n_sqmatrices];
    BeginIndex = new int* [n_sqmatrices];
    HangingEntries = new double* [n_sqmatrices];

   } //if(N_SqMatric 
   
  if(n_matrices)
  {
    TestGlobalNumbers = new int* [n_matrices];
    AnsatzGlobalNumbers = new int* [n_matrices];
    TestBeginIndex = new int* [n_matrices];
    AnsatzBeginIndex = new int* [n_matrices];
  } // endif n_matrices
  
  if(n_rhs)
   {   
    Rhs = new double*[n_rhs];     
    RhsBeginIndex = new int* [n_rhs];
    RhsGlobalNumbers = new int* [n_rhs];
    LocRhs = new double* [n_rhs]; 
   }

 if(N_AllMatrices)
  LocMatrices = new double** [N_AllMatrices]; 
  
}

TAssemble::~TAssemble()
{
  
  
}

