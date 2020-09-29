// =======================================================================
// @(#)IsoInterfaceJoint3D.h        1.1 08/12/99
// 
// Class:       TIsoInterfaceJoint3D
// Purpose:     connects two cells on an interface with additional
//              vertices for isoparametric reference transformation
//
// Author:      Gunar Matthies  06.08.99
//              Gunar Matthies  05.04.02
//
// =======================================================================

#ifndef __ISOINTERFACEJOINT3D__
#define __ISOINTERFACEJOINT3D__

#include <InterfaceJoint3D.h>

/** connects two cells on an interface */
class TIsoInterfaceJoint3D : public TInterfaceJoint3D
{
  protected:
    /** number of additional vertices */
    int N_Vertices;

    /** array of all additional vertices */
    TVertex **Vertices;

  public:
    // Constructors
    /** initialize the joint with the boundary parameters and one neighbour */
    TIsoInterfaceJoint3D(TBoundComp3D *bdcomp, double *param1, double *param2,
                       TBaseCell *neigh0);
    /** initialize the joint with the boundary parameters and two neighbours */
    TIsoInterfaceJoint3D(TBoundComp3D *bdcomp, double *param1, double *param2,
                       TBaseCell *neigh0, TBaseCell *neigh1);
    /** initialize the joint with the boundary parameters and two neighbours */
    TIsoInterfaceJoint3D(TBoundComp3D *bdcomp, TBaseCell *neighb0);

    // Methods
    /** create a new instance of the same class */
    virtual TJoint *NewInst(double T_0, double T_1, TBaseCell *Me);
    virtual TJoint *NewInst();

    /** return number of additional vertices */
    int GetN_Vertices()
    { return N_Vertices; }

    TVertex **GetVertices()
    { return Vertices; }

    void SetVertices(int n_vertices, TVertex **vertices);

    void GenerateVertices(int n_vertices);
    
    // Destructor
    virtual ~TIsoInterfaceJoint3D();    
};

#endif
