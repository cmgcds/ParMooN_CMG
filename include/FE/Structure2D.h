// =======================================================================
// @(#)Structure2D.h        1.3 09/14/99
// 
// Class:       TStructure2D
//
// Purpose:     build and store a matrix Structure2D
//
// Author:      Gunar Matthies
//
// History:     24.11.97 start implementation
//
//              start of reimplementation 26.08.1998 (Gunar Matthies)
//
// =======================================================================

#ifndef __STRUCTURE2D__
#define __STRUCTURE2D__

#include <FESpace1D.h>
#include <FESpace2D.h>
#include <Structure.h>

class TStructure2D : public TStructure
{
  protected:
    /** Ansatzspace */
    TFESpace1D *AnsatzSpace1D;
    TFESpace2D *AnsatzSpace2D;

    /** Testspace */
    TFESpace1D *TestSpace1D;
    TFESpace2D *TestSpace2D;

    int *AnsatzMortarSpaceGlobNo;
    int *TestMortarSpaceGlobNo;

    int *AnsatzNonMortarSpaceGlobNo;
    int *TestNonMortarSpaceGlobNo;

  public:
    /** generate the matrix Structure2D, both space with 2D collection */
    TStructure2D(TFESpace2D *testspace, TFESpace2D *ansatzspace);

    /** destructor: free all used arrays */
    ~TStructure2D();

    /** generate the matrix structure, both spaces are 2D */
    /** both spaces are defined on different grids */
    TStructure2D(TFESpace2D *testspace, int test_level, 
                 TFESpace2D *ansatzspace, int ansatz_level);
    #ifdef __MORTAR__
    /** generate the matrix Structure2D, one space with 1D and the other
        with 2D collection */
    TStructure2D(TFESpace1D *testspace, TFESpace2D *ansatzspace);
    #endif

    /** generate the matrix Structure2D, one space with 1D and the other
        with 2D collection */
     TStructure2D(TFESpace1D *testspace, TFESpace2D *ansatzspace, int **ansatzcelljoints);
     
    /** generate the matrix Structure2D, one space with 1D and the other
        with 2D collection */
     TStructure2D(TFESpace1D *testspace, TFESpace2D *ansatzspace, TNonMortarData *NonMortarFEData);
     
     TStructure2D(TFESpace2D *testspace, TFESpace1D *ansatzspace, TNonMortarData *NonMortarFEData);

    /** return AnsatzSpace */
    TFESpace2D *GetAnsatzSpace2D() const
    { return AnsatzSpace2D; }
    
    /** return AnsatzSpace */
    TFESpace *GetAnsatzSpace()
    {
      if (AnsatzSpace1D)
        return AnsatzSpace1D;
      else
        return AnsatzSpace2D;
    }
    
    /** return TestSpace */
    TFESpace2D *GetTestSpace2D() const
    { return TestSpace2D; }
    
    /** return TestSpace */
    TFESpace *GetTestSpace()
    {
      if (TestSpace1D)
        return TestSpace1D;
      else
        return TestSpace2D;
    }

};

#endif
