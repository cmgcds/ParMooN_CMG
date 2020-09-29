/** ************************************************************************ 
*
* @class     TAssembleMat3D
* @brief     base class for assembling matrices 
* @author    Sashikumaar Ganesan 
* @date      30.05.15
* @History    
 ************************************************************************  */


#ifndef __ASSEMBLEMAT3D__
#define __ASSEMBLEMAT3D__

#include <AllClasses.h>
#include <Assemble.h>
#include <Enumerations.h>

#ifdef _MPI

#endif

#ifdef _OMPONLY

#endif

/** class for 3D scalar assemble */
class TAssembleMat3D: public TAssemble
{
  protected:  

  TFESpace3D **FeSpaces,  **FeRhs;
  TSquareMatrix3D **SqMatrices; 
  TMatrix3D **RecMatrices;
  TDiscreteForm3D *DiscreteForm;
  BoundCondFunct3D **BoundaryConditions;
  BoundValueFunct3D **BoundaryValues;
  TAuxParam3D *AuxParam;                     
  double *Param[MaxN_QuadPoints_3D], *AuxArray[MaxN_QuadPoints_3D];  

 
  public:
    /** Constructor*/
    TAssembleMat3D(int n_fespaces, TFESpace3D **fespaces,
                      int n_sqmatrices, TSquareMatrix3D **sqmatrices,
                      int n_matrices, TMatrix3D **matrices,
                      int n_rhs, double **rhs, TFESpace3D **ferhs,
                      TDiscreteForm3D *discreteform,
                      BoundCondFunct3D **boundarybonditions,
                      BoundValueFunct3D **boundaryvalues,
                      TAuxParam3D *parameters);
    
    /** destrcutor */
    ~TAssembleMat3D();
    
  /** allocate memorey for all aux array */
  void Init();
    
   /** deallocate memorey for all aux array */
  void DeAllocate();

  /** reset all values of Mat and rhs to zero */
  void Reset();
  
   /** assemble the matrices */
  void Assemble3D(); 
  
  /** set No penetration BC, that is free slip */
  void AssembleNavierSlip();
  
private:  
  /** Add Local Square matrix to Global Matrix */
  void AddLocalSqMatToGlobal(int i, TBaseCell *cell, int *N_BaseFunct);
 
  /** Add Local matrix to Global Matrix */
  void AddLocalRecMatToGlobal(int i, TBaseCell *cell, int *N_BaseFunct);

  /** Add Local rhs to Global rhs */
  void AddLocalRhsToGlobal(int i, TBaseCell *cell, int *N_BaseFunct, BaseFunct3D *BaseFuncts, RefTrans3D reftrans);
  
  /**  ModifyMatHang */
  void ModifyMatHang();

  /**  print all matrices */  
  void PrintAllMat();
  
};

#endif

