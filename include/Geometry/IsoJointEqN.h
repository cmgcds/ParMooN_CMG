// =======================================================================
// @(#)IsoJointEqN.h        1.2 10/18/99
// 
// Class:       TIsoJointEqN
// Purpose:     connects two cells
//              has additional vertices for isoparametric reference 
//              transformation
//
// Author:      Gunar Matthies 06.08.1999
//
// History:     start of implementation (Gunar Matthies 06.08.99)
//
// =======================================================================

#ifndef __ISOJOINTEQN__
#define __ISOJOINTEQN__

#include <JointEqN.h>
#include <Vertex.h>

/** connects two cells */
class TIsoJointEqN : public TJointEqN
{
  protected:
    /** number of additional vertices */
    int N_Vertices;

    /** array of all additional vertices */
    TVertex **Vertices;

  public:
    // Constructors
    /** constructor with one initial neighbour */
    TIsoJointEqN(TBaseCell *neighb0);
    /** constructor with two initial neighbours */
    TIsoJointEqN(TBaseCell *neighb0, TBaseCell *neighb1);

    // Methods
    /** create a new instance of this class */
    virtual TJoint *NewInst(double T_0, double T_1, TBaseCell *Me)
    { return new TIsoJointEqN(Me); }
    virtual TJoint *NewInst()
    { return new TIsoJointEqN(NULL); }

    /** return number of additional vertices */
    int GetN_Vertices()
    { return N_Vertices; }

    TVertex **GetVertices()
    { return Vertices; }

    void SetVertices(int n_vertices, TVertex **vertices);
};

#endif
