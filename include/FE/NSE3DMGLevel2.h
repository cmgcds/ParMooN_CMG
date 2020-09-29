// =======================================================================
// @(#)NSE3DMGLevel2.h        1.6 07/03/00
//
// Class:       TNSE3DMGLevel2
// Purpose:     store all data for one level in a multi grid method
//              for solving a Stokes-/ Navier-Stokes system of
//              type 2 (A, B1, B2, B3, B1T, B2T, B3T)
//
// Author:      Volker John 25.07.2000
//
// History:     25.07.2000 start of implementation
//
// =======================================================================

#ifndef __NSE3DMGLEVEL2__
#define __NSE3DMGLEVEL2__

#include <NSE3DMGLevel.h>

class TNSE3DMGLevel2 : public TNSE3DMGLevel
{
  protected:
    /** matrix A */
    TSquareMatrix3D *A;

    /** structure of matrix A */
    TSquareStructure3D *StructureA;

    /** row pointer for matrix A */
    int *ARowPtr;

    /** column number vector for matrix A */
    int *AKCol;

    /** matrix entries of matrix A */
    double *AEntries;

    /** matrix B1 */
    TMatrix3D *B1;

    /** matrix B2 */
    TMatrix3D *B2;

    /** matrix B3 */
    TMatrix3D *B3;

    /** structure of matrix B */
    TStructure3D *StructureB;

    /** row pointer for matrix B */
    int *BRowPtr;

    /** column number vector for matrix B */
    int *BKCol;

    /** matrix entries of matrix B1 */
    double *B1Entries;

    /** matrix entries of matrix B2 */
    double *B2Entries;

    /** matrix entries of matrix B3 */
    double *B3Entries;

    /** matrix B1T */
    TMatrix3D *B1T;

    /** matrix B2 */
    TMatrix3D *B2T;

    /** matrix B3 */
    TMatrix3D *B3T;

    /** structure of matrix BT */
    TStructure3D *StructureBT;

    /** row pointer for matrix BT */
    int *BTRowPtr;

    /** column number vector for matrix BT */
    int *BTKCol;

    /** matrix entries of matrix BT1 */
    double *B1TEntries;

    /** matrix entries of matrix BT2 */
    double *B2TEntries;

    /** matrix entries of matrix BT3 */
    double *B3TEntries;

  public:
    /** constructor */
    TNSE3DMGLevel2(int level, TSquareMatrix3D *A, 
                   TMatrix3D *B1, TMatrix3D *B2, TMatrix3D *B3,
                   TMatrix3D *B1T, TMatrix3D *B2T, TMatrix3D *B3T,
                   double *f1, double *f2, double *f3, double *g,
                   double *u1, double *u2, double *u3, double *p,
                   int n_aux, double al, int VelocitySpace, 
                   int PressureSpace, TCollection *coll);

    /** destructor */
    ~TNSE3DMGLevel2();

    /** calculate defect */
    virtual void Defect(double *u1, double *u2, double *u3, double *p,
                        double *f1, double *f2, double *f3, double *g,
                        double *d1, double *d2, double *d3, double *d4,
                        double &res);

    /** correct Dirichlet and hanging nodes */
    virtual void CorrectNodes(double *u1, double *u2, double *u3, double *p);

    /** Vanka smoother */
    virtual void CellVanka(double *u1, double *u2,  double *u3, double *p,
        double *rhs1, double *rhs2, double *rhs3, double *rhs4,
        double *def1, double *def2, double *def3, double *def4,
        double *aux, double *Counters,
        int N_Parameters, double *Parameters, int smoother, int N_Levels);

    /** Vanka smoother */
    virtual void NodalVanka(double *u1, double *u2, double *u3, double *p,
        double *rhs1, double *rhs2, double *rhs3, double *rhs4,
        double *def1, double *def2, double *def3, double *def4,
        double *aux, double *Counters,
        int N_Parameters, double *Parameters, int smoother, int N_Levels);

    /** solve exact on this level */
    virtual void SolveExact(double *u1, double *u2, double *u3, double *p,
        double *rhs1, double *rhs2, double *rhs3, double *rhs4,
        double *aux, double *Counters,
        int N_Parameters, double *Parameters);

    /** Braess Sarazin smoother */
    virtual void BraessSarazin(double *u1, double *u2, double *u3, double *p,
        double *rhs1, double *rhs2, double *rhs3, double *rhs4,
        double *aux, double *Counters,
        int N_Parameters, double *Parameters,int N_Levels);

    /** step length control for Vanka */
    virtual double StepLengthControl(double *u1, double *u2, double *u3, double *p,
           double *u1old, double *u2old, double *u3old, double *pold,
           double *def1, double *def2, double *def3, double *def4,
           int N_Parameters, double *Parameter);

    /** print all matrices and both right hand sides */
    virtual void PrintAll();

};

#endif
