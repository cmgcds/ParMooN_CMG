// =======================================================================
// @(#)MortarBaseJoint.h        1.2 08/27/99
// 
// Class:       TMortarBaseJoint
// Purpose:     mortar joint on macro grid level
//
// Author:      Volker Behns  19.03.98
//
// =======================================================================

#ifndef __MORTARBASEJOINT__
#define __MORTARBASEJOINT__

#include <JointEqN.h>
#include <MortarJoint.h>

/** mortar joint on macro grid level */
class TMortarBaseJoint : public TJointEqN
{
  protected:
    /** FE-number on mortar space */
    int MEdgeInColl;

  public:
    // Constructors
    /** constructor with two initial neighbours */
    TMortarBaseJoint(TBaseCell *neighb0, TBaseCell *neighb1);

    // Methods
   /** checking for matching joints */
    virtual int CheckMatchingRef(TBaseCell *Me, int J_i,
                     struct StoreGeom &Tmp);

    #ifdef __MORTAR__
      /** check the refinement pattern on both sides for matching,
          special version for moratr cells */
      virtual int CheckMatchingRef(TBaseCell *Me, int J_i,
                    StoreGeomMortar &Tmp);
    #endif

    /** create a new instance of this class */
    virtual TJoint *NewInst(double T_0, double T_1, TBaseCell *Me)
    { return new TMortarJoint(); }
    virtual TJoint *NewInst()
    { return new TMortarJoint(); }

    /** set MEdgeInColl */
    void SetMEdgeInColl(int medgeincoll)
    { MEdgeInColl = medgeincoll; }

    /** return MEdgeInColl */
    int GetMEdgeInColl()
    { return MEdgeInColl; }

    /** return whether this is an interior joint */
    virtual bool InnerJoint() const
    { return false; }
    
    // Destructor
    virtual ~TMortarBaseJoint();
};
#endif
