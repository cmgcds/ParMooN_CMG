// =======================================================================
// @(#)Structure3D.h        1.3 09/14/99
// 
// Class:       TStructure3D
//
// Purpose:     build and store a matrix Structure3D
//
// Author:      Gunar Matthies
//
// History:     24.11.97 start implementation
//
//              start of reimplementation 26.08.1998 (Gunar Matthies)
//
// =======================================================================

#ifndef __STRUCTURE3D__
#define __STRUCTURE3D__

#include <FESpace2D.h>
#include <FESpace3D.h>
#include <Structure.h>

class TStructure3D : public TStructure
{
  protected:
    /** Ansatzspace */
    TFESpace2D *AnsatzSpace2D;
    TFESpace3D *AnsatzSpace3D;

    /** Testspace */
    TFESpace2D *TestSpace2D;
    TFESpace3D *TestSpace3D;

  public:
    /** generate the matrix Structure3D, both space with 3D collection */
    TStructure3D(TFESpace3D *testspace, TFESpace3D *ansatzspace);

    /** return AnsatzSpace */
    TFESpace *GetAnsatzSpace()
    {
      if (AnsatzSpace2D)
        return AnsatzSpace2D;
      else
        return AnsatzSpace3D;
    }

    /** return TestSpace */
    TFESpace *GetTestSpace()
    {
      if (TestSpace2D)
        return TestSpace2D;
      else
        return TestSpace3D;
    }

};

#endif
