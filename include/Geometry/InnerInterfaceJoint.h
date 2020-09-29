// =======================================================================
// @(#)InnerInterfaceJoint.h        1.3 01/26/12
// 
// Class:       TInnerInterfaceJoint
// Purpose:     connects two cells of different collection
//              (not treated as a boundary joint)
//
//
// =======================================================================

#ifndef __INNERINTERFACEJOINT__
#define __INNERINTERFACEJOINT__

#include <JointEqN.h>
#include <Vertex.h>

/** connects two cells on an interface */
class TInnerInterfaceJoint : public TJointEqN
{
  protected:
    /** during refinement the two children are marked in the following list */
    TInnerInterfaceJoint *children[2];
    
    /** x coordinate of the begin of line */
    double Xstart;
    /** y coordinate of the begin of line */
    double Ystart;
    /** x progress of line */
    double delX;
    /** y progress of line */
    double delY;

    /** The index of this joint in the two neighbors */
    int IndexInNeighbor[2];
 public:
  /** global cell number of the neibs' cell, which contains this joint */
  int NeibGlobalCellNo[2];
  int SubDomainsID[2];
  
  // Constructors
  /** constructor with one initial neighbour */
  TInnerInterfaceJoint(TBaseCell *neighb0);
  TInnerInterfaceJoint(TBaseCell *neighb0, TBaseCell *neighb1);
  
  // Methods
    /** create a new instance of the same class */
  virtual TJoint *NewInst(double T_0, double T_1, TBaseCell *Me);
  virtual TJoint *NewInst();
  
  
  virtual bool InnerJoint() const
  { return true; }
  
  /** Remember which joints are children of this joint */
  void SetChild(TInnerInterfaceJoint *child);
  /** return one of the two children of this edge */
  TInnerInterfaceJoint *GetChild(int child) const;
  
  /** Get either one of the two neighbors */
  TBaseCell *GetNeighbour(int i) const;
  
  /** set/get the coordinates of the start point and the vector pointing to the     
    second point */
  void SetParams (double xstart, double ystart, double delx, double dely);
  void GetParams (double &xstart, double &ystart, double &delx, double &dely) const;
  
  /** Compute the length of this edge */
  double GetLength() const;
  
  /** the unit normal of this edge */
  void GetNormal(double &nx, double &ny) const;
  
  /** the unit tangential of this edge */
  void GetTangent(double &tx, double &ty) const;
  
  /** set/get the index of this joint in given neighbor */
  int GetIndexInNeighbor(TBaseCell const * const neigh);
  void SetIndexInNeighbor(TBaseCell *neigh, int index);
  
    // Destructor
  virtual ~TInnerInterfaceJoint(){};
  
};

#endif
