// =======================================================================
// @(#)Brick.C        1.1 10/30/98
//
// Class:       TBrick
// Purpose:     shape descriptor of a quadrangle, especially a
//              Brick
//
// Author:      Gunar Matthies  11.07.2000
//
// =======================================================================

#include <Brick.h>

// Constructor
TBrick::TBrick() : THexahedron()
{
  Type = Brick;
}

// Methods
