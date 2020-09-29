// =======================================================================
// @(#)IsoBoundEdge.h        1.1 08/12/99
// 
// Class:       TIsoBoundEdge
// Purpose:     edge on a boundary component with additional vertices
//              for isoparametric reference transformation
//
// Author:      Gunar Matthies  06.08.99
//
// =======================================================================

#ifndef __ISOBOUNDEDGE__
#define __ISOBOUNDEDGE__

#include <BoundEdge.h>

/** edge on a boundary component */
class TIsoBoundEdge : public TBoundEdge
{
  protected:
    /** number of additional vertices */
    int N_Vertices;

    /** array of all additional vertices */
    TVertex **Vertices;

  public:
    // Constructors
    /** initialize the edge with the boundary component bdcomp and the
        paramter of starting and end point t\_0, t\_1 */
    TIsoBoundEdge(TBoundComp2D *bdcomp, double t_0, double t_1);

    // Methods
    /** check whether the refinement pattern on both side patch,
        dummy here: there is no neighbour */
    virtual int CheckMatchingRef(TBaseCell *Me, int J_i,
                  struct StoreGeom &Tmp);

    #ifdef __MORTAR__
      /** check the refinement pattern on both sides for matching,
          special version for moratr cells */
      virtual int CheckMatchingRef(TBaseCell *Me, int J_i,
                    StoreGeomMortar &Tmp);
    #endif

    /** create a new instance of this class */
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

    void ChangeEdgeID(JointType New_ID)
     { ID = New_ID; }
     
    void DeleteVertices();
};

#endif
