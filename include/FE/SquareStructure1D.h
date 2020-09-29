// =======================================================================
// @(#)SquareStructure1D.C        
// 
// Class:       TSquareStructure1D
//
// Purpose:     build and store a structure for a square matrix in 1d
//
// Author:      Sashikumaar Ganesan
//
// History:     17.05.2007 start implementation
//
// =======================================================================

#ifndef __SQUARESTRUCTURE1D__
#define __SQUARESTRUCTURE1D__

#include <FESpace1D.h>
#include <SquareStructure.h>

class TSquareStructure1D : public TSquareStructure
{
  protected:
    /** FE space */
    TFESpace1D *FESpace;

  public:
    /** dummy constructor, needed only derives classes */
    TSquareStructure1D();

    /** generate the matrix structure, only one space needed */
    TSquareStructure1D(TFESpace1D *space);

    /** destructor: free all used arrays */
    ~TSquareStructure1D();

    /** return FESpace */
    TFESpace1D *GetFESpace()
    { return FESpace; }

};

#endif

 
