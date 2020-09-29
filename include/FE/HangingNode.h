// =======================================================================
// @(#)HangingNode.h        1.1 10/30/98
//
// Class:       THangingNode
// Purpose:     represent a hanging node
//
// Author:      Gunar Matthies  18.11.97
//
// =======================================================================

#ifndef __HANGINGNODE__
#define __HANGINGNODE__

#include <Constants.h>
#include <Enumerations.h>
#include <MooNMD_Io.h>

/** represent a hanging node */
class THangingNode
{
  protected:
    /** type of hanging node */
    HNDesc Type;

    /** numbers of degrees of freedom in coupling */
    int *DOF;

  public:
    /** constructor, filling all data */
    THangingNode(HNDesc type, int n_, int *all_dofs, int *coupling);

    /** destructor */
    ~THangingNode();

    // Methods
    /** return number of nodes in coupling */
    HNDesc GetType()
    { return Type; }

    /** return numbers of degrees of freedom in coupling */
    int *GetDOF()
    { return DOF; }
   
};

#endif
