// =======================================================================
// @(#)IsoInterfaceJoint.h        1.1 08/12/99
// 
// Class:       TIsoInterfaceJoint
// Purpose:     connects two cells on an interface with additional
//              vertices for isoparametric reference transformation
//
// Author:      Gunar Matthies  06.08.99
//
// =======================================================================

#ifndef __ISOINTERFACEJOINT__
#define __ISOINTERFACEJOINT__

#include <InterfaceJoint.h>

/** connects two cells on an interface */
class TIsoInterfaceJoint : public TInterfaceJoint
{
  protected:
    /** number of additional vertices */
    int N_Vertices;

    /** array of all additional vertices */
    TVertex **Vertices;

  public:
    // Constructors
    /** initialize the joint with the boundary parameters and one neighbour */
    TIsoInterfaceJoint(TBoundComp2D *bdcomp, double t_0, double t_1,
                    TBaseCell *neighb0);
    /** initialize the joint with the boundary parameters and two neighbour */
    TIsoInterfaceJoint(TBoundComp2D *bdcomp, double t_0, double t_1,
                    TBaseCell *neighb0, TBaseCell *neighb1);

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
    
    void GeneratemidVert(int n_vertices, double*X, double*Y);
    
    // Destructor
    virtual ~TIsoInterfaceJoint(){};
       
};

#endif
