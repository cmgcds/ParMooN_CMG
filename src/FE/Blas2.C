// =======================================================================
// BLAS 2
// =======================================================================

#include <LinAlg.h>
#include <string.h>
#include <Database.h>

/** scalar system */
/** y := b - A *x */
void ScalarDefect(TSquareMatrix *A, double *sol, double *f, double *d,
                  double &res)
{
  int i,j,k,l,index,numThreads;
  double s;
  int *RowPtr, *KCol;
  double *Entries;
  int N_DOF; 

  RowPtr = A->GetRowPtr();
  KCol = A->GetKCol();
  Entries = A->GetEntries();
  N_DOF = A->GetN_Rows();

  res = 0.0;
  j = RowPtr[0];
 
#ifdef _HYBRID
  numThreads = TDatabase::ParamDB->OMPNUMTHREADS;
  omp_set_num_threads(numThreads);
#pragma omp parallel default(shared) private(i,s,k,j,index) 
{
  #pragma omp for schedule(guided) nowait
#endif
  for(i=0;i<N_DOF;i++)
  {
    s = f[i];
    k = RowPtr[i+1];
    for(j=RowPtr[i];j<k;j++)
    {
      index = KCol[j];
      s -= Entries[j] * sol[index];
    }
    d[i] = s;
#ifndef _MPI     
    res += s*s;
#endif      
  } // endfor i
#ifdef _HYBRID
}
#endif
  
#ifndef _MPI    
  res = sqrt(res);
#endif   
} // end Defect

/** scalar systems */ 
/** y := A * x (for ONE a block) */
void MatVect(TSquareMatrix *A, double *x, double *y)
{
  int N_UDOF;
  int i,j,k,l,index;
  double s, value;
  int *ARowPtr, *AKCol;
  double *AEntries;

  ARowPtr = A->GetRowPtr();
  AKCol = A->GetKCol();
  AEntries = A->GetEntries();

  N_UDOF = A->GetN_Rows();

  j = ARowPtr[0];
  for(i=0;i<N_UDOF;i++)
  {
    s = 0;
    k = ARowPtr[i+1];
    for(;j<k;j++)
    {
      index = AKCol[j];
//      value = AEntries[j];
      s += AEntries[j] * x[index];
    }
    y[i] = s;
  } // endfor i
}

void MatVect_Scalar(TSquareMatrix **A, TMatrix **B, double *x, double *y)
{
  MatVect(A[0],x,y);
}
void Defect_Scalar(TSquareMatrix **A, TMatrix **B, double *x, double *b, 
                   double *r)
{
  double res;

  ScalarDefect(A[0],x,b,r,res);
}

/** y := A * x (for ONE a block), active rows only */
void MatVectActive(TSquareMatrix *A, double *x, double *y)
{
  int N_Active;
  int i,j,k,l,index;
  double s, value;
  int *ARowPtr, *AKCol;
  double *AEntries;

  ARowPtr = A->GetRowPtr();
  AKCol = A->GetKCol();
  AEntries = A->GetEntries();

  N_Active = A->GetActiveBound();

  j = ARowPtr[0];
  for(i=0;i<N_Active;i++)
  {
    s = 0;
    k = ARowPtr[i+1];
    for(;j<k;j++)
    {
      index = AKCol[j];
      value = AEntries[j];
      s += value * x[index];
    }
    y[i] = s;
  } // endfor i
}

/** C := C + alpha*D, active rows only */
void MatAdd(TSquareMatrix *C, TSquareMatrix *D, double alpha)
{
  int N_UDOF;
  int i,j,k,l,index;
  double s, value;
  int *CRowPtr;
  double *CEntries;
  double *DEntries;

  CRowPtr = C->GetRowPtr();
  CEntries = C->GetEntries();
  DEntries = D->GetEntries();

  N_UDOF = C->GetActiveBound();

  l=CRowPtr[N_UDOF];

  Daxpy(l, alpha, DEntries, CEntries);
}

/** C := alpha*C + D, active rows only */
void MatAdd2(TSquareMatrix *C, TSquareMatrix *D, double alpha)
{
  double *CEntries, *DEntries;
  int *RowPtr, N_, N_Active;
  
  RowPtr = C->GetRowPtr();
  CEntries = C->GetEntries();
  DEntries = D->GetEntries();
  
  N_Active = C->GetActiveBound();
  
  N_ = RowPtr[N_Active];
  
  Dscal(N_, alpha, CEntries);
  Daxpy(N_, 1.0, DEntries, CEntries);
}

void MatAdd(TMatrix *C, TMatrix *D, double alpha)
{
  int N_UDOF;
  int i,j,k,l,index;
  double s, value;
  int *CRowPtr;
  double *CEntries;
  double *DEntries;

  CRowPtr = C->GetRowPtr();
  CEntries = C->GetEntries();
  DEntries = D->GetEntries();

  N_UDOF = C->GetN_Rows();

  l=CRowPtr[N_UDOF];

  Daxpy(l, alpha, DEntries, CEntries);
}

/** y := B * x */
void MatVect1(TMatrix *B, double *x, double *y)
{
  int N_UDOF;
  int i,j,k,l,index;
  double s, value;
  int *BRowPtr, *BKCol;
  double *BEntries;

  BRowPtr = B->GetRowPtr();
  BKCol = B->GetKCol();
  BEntries = B->GetEntries();
  N_UDOF = B->GetN_Rows();

  j = BRowPtr[0];
  for(i=0;i<N_UDOF;i++)
  {
    s = 0;
    k = BRowPtr[i+1];
    for(;j<k;j++)
    {
      index = BKCol[j];
      value = BEntries[j];
      s += value * x[index];
    }
    y[i] = s;
  } // endfor i
}

/** y := B' * x (transposed matrix * vector) */
void TransMatVect(TMatrix *B, double *x, double *y)
{
  int N_UDOF;
  int i,j,k,l,index;
  double s, value;
  int *BRowPtr, *BKCol;
  double *BEntries;

  BRowPtr = B->GetRowPtr();
  BKCol = B->GetKCol();
  BEntries = B->GetEntries();

  N_UDOF = B->GetN_Rows();
  memset(y, 0, B->GetN_Columns()*SizeOfDouble);

  j = BRowPtr[0];
  for(i=0;i<N_UDOF;i++)
  {
    k = BRowPtr[i+1];
    for(;j<k;j++)
    {
      index = BKCol[j];
      value = BEntries[j];
      y[index] += value * x[i];
    }
  } // endfor i
}

