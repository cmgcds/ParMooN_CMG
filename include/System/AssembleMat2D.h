/** ************************************************************************ 
*
* @class     TAssembleMat2D
* @brief     base class for assembling matrices 
* @author    Sashikumaar Ganesan 
* @date      16.04.16
* @History    
 ************************************************************************  */


#ifndef __ASSEMBLEMAT2D__
#define __ASSEMBLEMAT2D__

#include <AllClasses.h>
#include <Assemble.h>
#include <Enumerations.h>

#ifdef _MPI

#endif

#ifdef _OMPONLY

#endif

/** class for 2D scalar assemble */
class TAssembleMat2D: public TAssemble
{
  protected:  

  TFESpace2D **FeSpaces,  **FeRhs;
  TSquareMatrix2D **SqMatrices; 
  TMatrix2D **RecMatrices;
  TDiscreteForm2D *DiscreteForm;
  BoundCondFunct2D **BoundaryConditions;
  BoundValueFunct2D **BoundaryValues;
  TAuxParam2D *AuxParam;                     
  double *Param[MaxN_QuadPoints_2D], *AuxArray[MaxN_QuadPoints_2D];  

 
  public:
    /** Constructor*/
    TAssembleMat2D(int n_fespaces, TFESpace2D **fespaces,
                      int n_sqmatrices, TSquareMatrix2D **sqmatrices,
                      int n_matrices, TMatrix2D **matrices,
                      int n_rhs, double **rhs, TFESpace2D **ferhs,
                      TDiscreteForm2D *discreteform,
                      BoundCondFunct2D **boundarybonditions,
                      BoundValueFunct2D **boundaryvalues,
                      TAuxParam2D *parameters);
    
    /** destrcutor */
    ~TAssembleMat2D();
    
  /** allocate memorey for all aux array */
  void Init();
    
   /** deallocate memorey for all aux array */
  void DeAllocate();

  /** reset all values of Mat and rhs to zero */
  void Reset();
  
   /** assemble the matrices */
  void Assemble2D(); 
  
  /** asssemble slip with friction BC */
  void AssembleNavierSlip();
  
private:  
  /** Add Local Square matrix to Global Matrix */
  void AddLocalSqMatToGlobal(int i, TBaseCell *cell, int *N_BaseFunct);
 
  /** Add Local matrix to Global Matrix */
  void AddLocalRecMatToGlobal(int i, TBaseCell *cell, int *N_BaseFunct);

  /** Add Local rhs to Global rhs */
  void AddLocalRhsToGlobal(int i, TBaseCell *cell, int *N_BaseFunct, BaseFunct2D *BaseFuncts, RefTrans2D reftrans);
  
  /**  ModifyMatHang */
  void ModifyMatHang();

  /**  print all matrices */  
  void PrintAllMat();
  
};

#endif

