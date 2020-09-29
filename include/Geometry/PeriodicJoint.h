// =======================================================================
// @(#)PeriodicJoint.h        1.1 04/21/99
// 
// Class:       TPeriodicJoint
// Purpose:     connects two cells on different parts of the domain
//              where we have periodic boundary conditions
//
// Author:      Volker Behns  16.04.99
//
// =======================================================================

#ifndef __PERIODICJOINT__
#define __PERIODICJOINT__

#include <JointEqN.h>

/** connects two cells with periodic boundary conditions */
class TPeriodicJoint : public TJointEqN
{
  public:
    // Constructors
    /** initialize the joint with one neighbour */
    TPeriodicJoint(TBaseCell *neighb0);

    /** initialize the joint with two neighbours */
    TPeriodicJoint(TBaseCell *neighb0, TBaseCell *neighb1);

    // Methods
    /** create a new instance of the same class */
    virtual TJoint *NewInst(double T_0, double T_1, TBaseCell *Me);
    virtual TJoint *NewInst();
};

#endif
