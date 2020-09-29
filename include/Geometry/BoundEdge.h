// =======================================================================
// @(#)BoundEdge.h        1.2 08/27/99
// 
// Class:       TBoundEdge
// Purpose:     edge on a boundary component
//
// Author:      Volker Behns  02.08.97
//
// =======================================================================

#ifndef __BOUNDEDGE__
#define __BOUNDEDGE__

#include <Joint.h>
#include <BoundComp2D.h>
#include <Vertex.h>

/** edge on a boundary component */
class TBoundEdge : public TJoint
{
  protected:
    /** boundary component to which this edge belongs to */
    TBoundComp2D *BoundComp;

    /** paramter of starting point */
    double T_0;
    /** parameter of end point */
    double T_1;

  public:
    // Constructors
    /** initialize the edge with the boundary component bdcomp and the
        paramter of starting and end point t\_0, t\_1 */
    TBoundEdge(TBoundComp2D *bdcomp, double t_0, double t_1);

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

    /** return start parameter T0 */
    double GetStartParameter()
    { return T_0; }

    /** return end paramter T1 */
    double GetEndParameter()
    { return T_1; }

    /** return parameters */
    void GetParameters(double& t0, double& t1)
    {
      t0 = T_0;
      t1 = T_1;
    }

#ifdef __2D__
    /** update parameters according to the new vertex positions */
    void UpdateParameters(TVertex *Begin, TVertex *End);
#endif

    /** return the coordinates {X,Y} of parameter value T */
    int GetXYofT(double T, double &X, double &Y);

    /** return boundary component */
    TBoundComp2D *GetBoundComp() const
    { return BoundComp; }

    /** return whether this is an interior joint */
    virtual bool InnerJoint() const
    { return false; }
    
    /** change the boundary component */
    void ChangeBoundComp(TBoundComp2D *New_BoundComp)
     {BoundComp = New_BoundComp;}    
};

#endif
