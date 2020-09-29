// =======================================================================
// @(#)Matrix2D.C        1.2 11/20/98
// 
// Class:       TMatrix2D
//
// Purpose:     store a  matrix2D (ansatz != test space)
//
// Author:      Gunar Matthies
//
// History:     26.08.1998 start implementation
//
// =======================================================================

#include <Database.h>
#include <Matrix2D.h>
#include <SquareMatrix2D.h>
#include <SquareStructure2D.h>
#include <Structure2D.h>
#include <string.h>
#include <stdlib.h>

TMatrix2D::TMatrix2D(TStructure2D *structure)
 : TMatrix(structure)
{
  this->structure = structure;
}

TMatrix2D::~TMatrix2D()
{
}


void TMatrix2D::resetNonActive()
{
  int n_active = this->structure->GetTestSpace2D()->GetN_ActiveDegrees();
  int * rowPtr = this->structure->GetRowPtr();
  int index_nonactive = rowPtr[n_active];
  int n_nonactive_entries = rowPtr[structure->GetN_Rows()] - index_nonactive;
  memset(Entries + index_nonactive, 0.0, n_nonactive_entries * SizeOfDouble);
}

TMatrix2D& TMatrix2D::operator*=(double alpha)
{
  int* RowPtr = structure->GetRowPtr();
  TFESpace2D *fespace = GetStructure()->GetTestSpace2D();
  int nDOFActive = fespace->GetN_ActiveDegrees();
  for (int i=0; i<nDOFActive; i++)
  {
    for (int j=RowPtr[i]; j<RowPtr[i+1]; j++)
    {
      Entries[j] *= alpha;
    }
  }
  return *this;
}


// add to matrices A and B
// note: only active DOF are added
// note: only works for matrices with the same sparsity pattern
TMatrix2D& operator+(const TMatrix2D & A, const TMatrix2D & B)
{
  double *AEntries, *BEntries, *CEntries;
  if (A.GetStructure() == B.GetStructure()) 
  {
    TMatrix2D *C = new TMatrix2D(A.GetStructure());
    AEntries = A.GetEntries();
    BEntries = B.GetEntries();
    CEntries = C->GetEntries();
    for (int i=0; i<A.GetN_Entries(); i++) 
    {
      CEntries[i] = AEntries[i] + BEntries[i];
    }
    return *C;
  } else 
  {
    cout << " Matrix2D: ERROR: can't add Matrices "
         << " with different sparse structures " << endl;
    exit(1);
  }
}

TMatrix2D& operator*(const TMatrix2D & A, const double alpha)
{
  double *AEntries, *CEntries;
  TMatrix2D *C = new TMatrix2D(A.GetStructure());
  AEntries = A.GetEntries();
  CEntries = C->GetEntries();

  //TFESpace2D *fespace = A.GetStructure()->GetAnsatzSpace2D();
  TFESpace2D *fespace = A.GetStructure()->GetTestSpace2D();
  int nDOFActive = fespace->GetN_ActiveDegrees();
  //int nDOF = fespace->GetN_DegreesOfFreedom();

  // multiply each active entry by alpha and write it into matrix C
  for (int i=0; i<nDOFActive; i++) 
    {
      CEntries[i] = alpha*AEntries[i];
    }
  return *C;
}

TMatrix2D& operator*(const double alpha, const TMatrix2D & A)
{// just to allow alpha*A as well (in addition to A*alpha)
  return A*alpha;
}

double* operator*(const TMatrix2D & A,const double* x)
{
  double *AEntries;
  int *ARowPtr,*AColIndex;
  AEntries = A.GetEntries();
  ARowPtr = A.GetRowPtr();
  AColIndex = A.GetKCol();

  TFESpace2D *fespace = A.GetStructure()->GetTestSpace2D();
  int nDOFActive = fespace->GetN_ActiveDegrees();
  int nDOF = fespace->GetN_DegreesOfFreedom();

  double *y=new double[nDOF];
  double value;
  int index;

  for(int i=0;i<nDOFActive;i++)
  {
    value = 0;
    for (int j=ARowPtr[i]; j<ARowPtr[i+1]; j++)
    {
      index = AColIndex[j];
      value += AEntries[j] * x[index];
    }
    y[i] = value;
  }
  
  for (int i=nDOFActive; i<nDOF; i++)
  {
    y[i]=x[i];
  }
  return y;
}



void AllocateMatricesNSE_2D(int mg_level,
			    TFESpace2D *velocity_space, 
			    TFESpace2D *pressure_space,
			    TSquareStructure2D *&sqstructureA, 
			    TSquareStructure2D *&sqstructureC, 
			    TStructure2D *&structureB, 
			    TStructure2D *&structureBT,
			    TSquareMatrix2D *&sqmatrixA,
			    TSquareMatrix2D *&sqmatrixA11,
			    TSquareMatrix2D *&sqmatrixA12,
			    TSquareMatrix2D *&sqmatrixA21,
			    TSquareMatrix2D *&sqmatrixA22,
			    TSquareMatrix2D *&sqmatrixC,
			    TMatrix2D *&matrixB1,
			    TMatrix2D *&matrixB2,
			    TMatrix2D *&matrixB1T,
			    TMatrix2D *&matrixB2T,
			    TSquareMatrix2D **MatricesA,
			    TSquareMatrix2D **MatricesA11,
			    TSquareMatrix2D **MatricesA12,
			    TSquareMatrix2D **MatricesA21,
			    TSquareMatrix2D **MatricesA22,
			    TSquareMatrix2D **MatricesC,
			    TMatrix2D **MatricesB1,			    
			    TMatrix2D **MatricesB2,			    
			    TMatrix2D **MatricesB1T,			    
			    TMatrix2D **MatricesB2T)
{
    // matrix structures
    structureB = new TStructure2D(pressure_space, velocity_space);
    structureBT = new TStructure2D(velocity_space, pressure_space);
    sqstructureA = new TSquareStructure2D(velocity_space);
    sqstructureA->Sort();
    
    // allocate matrices
    switch(TDatabase::ParamDB->NSTYPE)
    {
        case 1:
	    matrixB1 = new TMatrix2D(structureB);
	    matrixB2 = new TMatrix2D(structureB);
	    
	    MatricesB1[mg_level] = matrixB1;
	    MatricesB2[mg_level] = matrixB2;
	    
	    sqmatrixA = new TSquareMatrix2D(sqstructureA);
	    
          MatricesA[mg_level] = sqmatrixA;
          break;

        case 2:
          matrixB1 = new TMatrix2D(structureB);
          matrixB2 = new TMatrix2D(structureB);
          matrixB1T = new TMatrix2D(structureBT);
          matrixB2T = new TMatrix2D(structureBT);

          MatricesB1[mg_level] = matrixB1;
          MatricesB2[mg_level] = matrixB2;
          MatricesB1T[mg_level] = matrixB1T;
          MatricesB2T[mg_level] = matrixB2T;

          sqmatrixA = new TSquareMatrix2D(sqstructureA);

          MatricesA[mg_level] = sqmatrixA;
          break;

        case 3:
          matrixB1 = new TMatrix2D(structureB);
          matrixB2 = new TMatrix2D(structureB);

          MatricesB1[mg_level] = matrixB1;
          MatricesB2[mg_level] = matrixB2;

          sqmatrixA11 = new TSquareMatrix2D(sqstructureA);
          sqmatrixA12 = new TSquareMatrix2D(sqstructureA);
          sqmatrixA21 = new TSquareMatrix2D(sqstructureA);
          sqmatrixA22 = new TSquareMatrix2D(sqstructureA);

          MatricesA11[mg_level] = sqmatrixA11;
          MatricesA12[mg_level] = sqmatrixA12;
          MatricesA21[mg_level] = sqmatrixA21;
          MatricesA22[mg_level] = sqmatrixA22;
          break;

        case 4:
          matrixB1 = new TMatrix2D(structureB);
          matrixB2 = new TMatrix2D(structureB);
          matrixB1T = new TMatrix2D(structureBT);
          matrixB2T = new TMatrix2D(structureBT);

          MatricesB1[mg_level] = matrixB1;
          MatricesB2[mg_level] = matrixB2;
          MatricesB1T[mg_level] = matrixB1T;
          MatricesB2T[mg_level] = matrixB2T;

          sqmatrixA11 = new TSquareMatrix2D(sqstructureA);
          sqmatrixA12 = new TSquareMatrix2D(sqstructureA);
          sqmatrixA21 = new TSquareMatrix2D(sqstructureA);
          sqmatrixA22 = new TSquareMatrix2D(sqstructureA);
	 
          MatricesA11[mg_level] = sqmatrixA11;
          MatricesA12[mg_level] = sqmatrixA12;
          MatricesA21[mg_level] = sqmatrixA21;
          MatricesA22[mg_level] = sqmatrixA22;
	 
          break;
        case 14:
	    sqstructureC = new TSquareStructure2D(pressure_space);
	    sqstructureC->Sort();
          matrixB1 = new TMatrix2D(structureB);
          matrixB2 = new TMatrix2D(structureB);
          matrixB1T = new TMatrix2D(structureBT);
          matrixB2T = new TMatrix2D(structureBT);

          MatricesB1[mg_level] = matrixB1;
          MatricesB2[mg_level] = matrixB2;
          MatricesB1T[mg_level] = matrixB1T;
          MatricesB2T[mg_level] = matrixB2T;

          sqmatrixA11 = new TSquareMatrix2D(sqstructureA);
          sqmatrixA12 = new TSquareMatrix2D(sqstructureA);
          sqmatrixA21 = new TSquareMatrix2D(sqstructureA);
          sqmatrixA22 = new TSquareMatrix2D(sqstructureA);
          sqmatrixC = new TSquareMatrix2D(sqstructureC);

          MatricesA11[mg_level] = sqmatrixA11;
          MatricesA12[mg_level] = sqmatrixA12;
          MatricesA21[mg_level] = sqmatrixA21;
          MatricesA22[mg_level] = sqmatrixA22;
          MatricesC[mg_level] = sqmatrixC;
          break;
      }
    return;
}
