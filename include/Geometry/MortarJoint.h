// =======================================================================
// @(#)MortarJoint.h        1.2 08/27/99
// 
// Class:       TMortarJoint
// Purpose:     indicates a mortar joint
//
// Author:      Volker Behns  09.03.98
//
// =======================================================================

#ifndef __MORTARJOINT__
#define __MORTARJOINT__

#include <Joint.h>

/** indicates a mortar joint */
class TMortarJoint : public TJoint
{
  protected:
    /** cell number in mortar collection */
    int MEdgeInColl;

  public:
    // Constructors
    TMortarJoint();

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
};

#endif
