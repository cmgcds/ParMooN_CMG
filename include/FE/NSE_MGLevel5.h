// =======================================================================
// @(#)NSE_MGLevel5.h        1.1 07/03/03
//
// Class:       TNSE_MGLevel5
// Purpose:     store all data for one level in a multi grid method
//              for solving a Stokes-/ Navier-Stokes system of
//              type 5 (A)
//
// Author:      Gunar Matthies (28.10.2003)
//
// History:     28.10.2003 start of implementation
//
// =======================================================================

#ifndef __NSE_MGLEVEL5__
#define __NSE_MGLEVEL5__

#include <NSE_MGLevel.h>

#ifdef __2D__
#include <SquareMatrixNSE2D.h>
#include <StructureNSE2D.h>
#endif

#ifdef __3D__
#include <SquareMatrixNSE3D.h>
#include <StructureNSE3D.h>
#endif

class TNSE_MGLevel5 : public TNSE_MGLevel
{
  protected:
#ifdef __2D__
    /** matrix A */
    TSquareMatrixNSE2D *A;

    /** structure of matrix A */
    TStructureNSE2D *StructureA;
#endif  

#ifdef __3D__
    /** matrix A */
    TSquareMatrixNSE3D *A;

    /** structure of matrix A */
    TStructureNSE3D *StructureA;
#endif  

    int *BeginJb;
    int *jb;
    int N_DOFperJoint;
    double *Alpha;

    /** number of cells in collection */
    int N_Cells;

   /** row pointer for matrix A */
    int *ARowPtr;

    /** column number vector for matrix A */
    int *AKCol;

    /** matrix entries of matrix A */
    double *AEntries;

  public:
    /** constructor */
#ifdef __2D__
    TNSE_MGLevel5(int level, TSquareMatrixNSE2D *A, 
                  double *f1, double *u1,
                  int n_aux, double *al,
                  int VelocitySpace, int PressureSpace,
                  TFESpace2D *PSpace, TCollection *coll);
#endif  

#ifdef __3D__
    TNSE_MGLevel5(int level, TSquareMatrixNSE3D *A, 
                  double *f1, double *u1,
                  int n_aux, double *al,
                  int VelocitySpace, int PressureSpace,
                  TFESpace3D *PSpace, TCollection *coll);
#endif  

    /** destructor */
    ~TNSE_MGLevel5();

    /** calculate defect */
    virtual void Defect(double *u1, double *f1, double *d1, double &res);

    /** correct Dirichlet and hanging nodes */
    virtual void CorrectNodes(double *u1);

    /** Vanka smoother */
    virtual void  CellVanka(double *u1, double *rhs1, double *aux,
                            int N_Parameters, double *Parameters, 
                            int smoother, int N_Levels);

     /** Vanka smoother */
    virtual void NodalVanka(double *u1, double *rhs1, double *aux,
        int N_Parameters, double *Parameters, int smoother, int N_Levels);

   /** solve exact on this level */
    virtual void SolveExact(double *u1, double *rhs1);

    /** solve exact on this level */
    virtual void SolveExactUMFPACK(double *u1, double *rhs1, int &umfpack_flag);

   /** Braess--Sarazin smoother */
    virtual void BraessSarazin(double *u1, double *rhs1, double *aux,
        int N_Parameters, double *Parameters,int N_Levels);

    /** step length control for Vanka */
    virtual double  StepLengthControl(double *u1, double *u1old, double *def1, 
           int N_Parameters, double *Parameter);

    /** print all matrices and both right hand sides */
    virtual void PrintAll();

    /** calculate ustar-representation from u-representation */
    void GetUstarFromU(double *u, double *ustar);

    /** calculate u-representation form ustar-representation */
    void GetUFromUstar(double *ustar, double *u);

    /** calculate dstar-representation from d-representation */
    void GetDstarFromD(double *d, double *dstar);

    /** calculate d-representation form dstar-representation */
    void GetDFromDstar(double *dstar, double *d);
    
    #ifdef _MPI
    virtual void UpdateHaloRhs(double*, double*)
    {}; 
    #endif
};

#endif
