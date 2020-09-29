// =======================================================================
// @(#)FE2DMapper.h        1.2 04/12/99
//
// Class:       TFE2DMapper
// Purpose:     find out which of the given local degress of freedom
//              are equivalent to the same global degree of freedom
//
// Author:      Gunar Matthies  05.11.97
//
// =======================================================================

#ifndef __FE2DMAPPER__
#define __FE2DMAPPER__

#define HANGINGNODE -5
#define FIRSTMARK -10

#include <FEDesc2D.h>
#include "Vector.h"
#include <HangingNode.h>

/** find out which of the given local degress of freedom,
    are equivalent to the same global degree of freedom */
class TFE2DMapper
{
  protected:
    /** name for the mapper */
    char *Name;

    /** some word describing the mapper */
    char *Description;

    /** number of local degrees on the first side */
    int N_DOF0;

    /** number of local degrees on the second side */
    int N_DOF1;

    /** number of formally different degrees */
    int N_Pairs;

    /** which pairs of local degrees matching in an array[][2] */
    int *Pairs;

    /** number of real degrees of freedom with no opposite */
    int N_NoOpposite;

    /** numbers of the eal degrees of freedom with no opposite */
    int *NoOpposite;

    /** memory for internal storage */
    int *Aux;

// =======================================================================
// local hanging nodes = coupling only this node on this edge
// =======================================================================
    /** number of hanging nodes in this pattern */
    int N_Hanging;

    /** indices of hanging nodes */
    int *Hanging;

    /** type of hanging nodes */
    HNDesc *HangingTypes;

    /** numbers of the degrees of freedom in coupling */
    int **Coupling;

// =======================================================================
// far hanging nodes = coupling also with node not on this edge
// =======================================================================
    /** number of hanging nodes in this pattern */
    int N_FarHanging;

    /** indices of hanging nodes */
    int *FarHanging;

    /** type of hanging nodes */
    HNDesc *FarHangingTypes;

    /** numbers of the degrees of freedom in coupling */
    int ****FarCoupling;

  public:
    /** constructor, filling all data */
    TFE2DMapper(char *name, char *description, int n0, int n1, 
              int n_pairs, int *pairs, 
              int n_hanging, int *hanging,
              HNDesc *hangingtypes, int **coupling,
              int n_farhanging, int *farhanging,
              HNDesc *farhangingtypes, int ****farcoupling,
              int n_noopposite, int *noopposite,
              int n_nodes);

    // Methods
    /** return name of mapper */
    char *GetName()
    { return Name; }

    /** return description of mapper */
    char *GetDescription()
    { return Description; }

    /** return number of degrees on side 0 */
    int GetN_DOF0()
    { return N_DOF0; }

    /** return number of degrees on side 1 */
    int GetN_DOF1()
    { return N_DOF1; }

    /** return number of degrees on both sides */
    void GetN_DOF(int &n0, int &n1)
    {
      n0 = N_DOF1;
      n1 = N_DOF1;
    }

    /** return number of pairs */
    int GetN_Pairs()
    { return N_Pairs; }

    /** return pairs of matching dof's */
    int *GetPairs()
    { return Pairs; }

    /** return type of hanging nodes */
    HNDesc *GetHangingTypes()
    { return HangingTypes; }

    /** return numbers of degrees of freedom in coupling */
    int **GetCoupling()
    { return Coupling; }

    /** return type of far hanging nodes */
    HNDesc *GetFarHangingTypes()
    { return FarHangingTypes; }

    /** return numbers of DOF in coupling (far hanging nodes) */
    int ****GetFarCoupling()
    { return FarCoupling; }

    /** return N_NoOpposite */
    int GetN_NoOpposite()
    { return N_NoOpposite; }

    /** return NoOpposite array */
    int *GetNoOpposite()
    { return NoOpposite; }

    /** map the given local degrees of freedom */
    void Map(int *Global, int I_K0, int I_K1, 
             int *Indices0, int *Indices1,
             int LocEdge0, int LocEdge1,
             TFEDesc2D *Desc0, TFEDesc2D *Desc1,
             int &Counter,
             TVector<THangingNode *> *vect,
             TVector<int> *numbers);

    /** "map" the given dof on a boundary joint */
    void MapBound(int *Global, int I_K, int *Indices, 
                  int &BoundCounter);

    /** map the two given degrees of freedom */
    void MapDOF(int *Global, int dof0, int dof1, int &Counter);

    /** destructor */
    ~TFE2DMapper()
    { delete Aux; }

};

#endif
