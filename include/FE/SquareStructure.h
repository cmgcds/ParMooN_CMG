// =======================================================================
// @(#)SquareStructure.h        1.6 09/14/99
// 
// Class:       TSquareStructure
//
// Purpose:     build and store a structure for a square matrix
//
// Author:      Gunar Matthies
//
// History:     06.08.1998 start implementation
//
// =======================================================================

#ifndef __SQUARESTRUCTURE__
#define __SQUARESTRUCTURE__

#include <Structure.h>

class TSquareStructure : public TStructure
{
  protected:
    /** number of active rows */
    int ActiveBound;

    /** ordering of the column entries */
    /** 0 - no special ordering */
    /** 1 - increasing ordering (like used for umfpack) */
    /** 2 - diagonal entry first, then increasing ordering */
    int ColOrder;

    /** sort an integer array */
    void IntSort(int *BeginPtr, int *AfterEndPtr);

  public:
    /** generate the matrix structure, only one space needed */
    TSquareStructure();

    /** destructor: free all used arrays */
    ~TSquareStructure();

    /** generate the matrix structure, all arrays are already defined */
    TSquareStructure(int n, int N_entries, int *col_ptr,
      int *row_ptr);
    
    /** Generates an empty n*n Structure for a Zero-Matrix */
    explicit TSquareStructure(int n);

    /** return ActiveBound */
    int GetActiveBound() const
    { return ActiveBound; }

    /** return ordering of columns */
    int GetColOrder() const
    { return ColOrder;}

    /** sort column numbers in each row, increasing indices */
    void Sort();

    void SetColOrder(int n)  
    {ColOrder = n;}    
    
    /** sort column numbers: diag is first element, other numbers are
    increasing */
    void SortDiagFirst();
    
    /** @brief Comparision Operators */
    friend bool operator==(const TSquareStructure &lhs, 
                           const TSquareStructure &rhs);
    friend bool operator!=(const TSquareStructure &lhs,
                           const TSquareStructure &rhs);
};

#endif
