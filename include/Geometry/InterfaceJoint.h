// =======================================================================
// @(#)InterfaceJoint.h        1.3 07/19/99
// 
// Class:       TInterfaceJoint
// Purpose:     connects two cells on an interface
//
// Author:      Volker Behns  09.03.98
//
// =======================================================================

#ifndef __INTERFACEJOINT__
#define __INTERFACEJOINT__

#include <JointEqN.h>
#include <Vertex.h>

/** connects two cells on an interface */
class TInterfaceJoint : public TJointEqN
{
  protected:
    /** boundary component to which this edge belongs */
    TBoundComp2D *BoundComp;

    /** parameter of starting point */
    double T_0;
    /** parameter of end point */
    double T_1;

  public:
    // Constructors
    /** initialize the joint with the boundary parameters and one neighbour */
    TInterfaceJoint(TBoundComp2D *bdcomp, double t_0, double t_1,
                    TBaseCell *neighb0);
    /** initialize the joint with the boundary parameters and two neighbours */
    TInterfaceJoint(TBoundComp2D *bdcomp, double t_0, double t_1,
                    TBaseCell *neighb0, TBaseCell *neighb1);

    // Methods
    /** create a new instance of the same class */
    virtual TJoint *NewInst(double T_0, double T_1, TBaseCell *Me);
    virtual TJoint *NewInst();

    /** return start parameter T0 */
    double GetStartParameter()
    { return T_0; }

    /** return end parameter T1 */
    double GetEndParameter()
    { return T_1; }

    /** return parameters */
    void GetParameters(double &t0, double &t1)
    {
      t0 = T_0;
      t1 = T_1;
    }

#ifdef __2D__
    /** update parameters according to the new vertex positions */
    void UpdateParameters(TVertex *Begin, TVertex *End);
#endif

    /** return boundary component */
    TBoundComp2D *GetBoundComp()
    { return BoundComp; }

    /** return the coordinates {X,Y} of parameter value T */
    int GetXYofT(double T, double &X, double &Y);

    /** make sure that Neighb0 is "inside" the domain */
    int CheckOrientation();

    /** check whether Neighb0 is "inside" the domain */
    bool CheckInside(TBaseCell *Me)
    { 
      if (Neighb0 == Me)
        return true;
      else
        return false;
    }
    
    // Destructor
    virtual ~TInterfaceJoint(){};
   
};

#endif
