// =======================================================================
// @(#)FE2DMapper1Reg.h        1.1 10/30/98
//
// Class:       TFE2DMapper1Reg
// Purpose:     find out which of the given local degress of freedom
//              are equivalent to the same global degree of freedom
//              for 1 regular case
//
// Author:      Gunar Matthies  05.11.97
//
// History:     start of reimplementation 03.08.1998 (GM)
//
// =======================================================================

#ifndef __FE2DMAPPER1REG__
#define __FE2DMAPPER1REG__

#include <FE2DMapper.h>

/** find out which of the given local degress of freedom,
    are equivalent to the same global degree of freedom 
    for 1 regular case */
class TFE2DMapper1Reg : public TFE2DMapper
{
  protected:
    /** number of local degrees of second element on second side */
    int N_DOF2;

    /** number of degrees of freedom at midpoint of long edge */
    int N_Mid;

    /** indices of degrees of freedom at midpoint of long edge */
    int *Mid;

  public:
    /** constructor, filling all data */
    TFE2DMapper1Reg(char *name, char *description, int n0, int n1, int n2,
              int n_pairs, int *pairs, 
              int n_mid, int *mid,
              int n_hanging, int *hanging,
              HNDesc *hangingtypes, int **coupling,
              int n_farhanging, int *farhanging,
              HNDesc *farhangingtypes, int ****farcoupling,
              int n_noopposite, int *noopposite,
              int n_nodes);

    /** destructor */
    ~TFE2DMapper1Reg()
    { }; 

    // Methods
    /** return number of degrees of second element on side 1 */
    int GetN_DOF2()
    { return N_DOF2; }

    /** return number of degrees on both sides */
    void GetN_DOF(int &n0, int &n1, int &n2)
    {
      n0 = N_DOF1;
      n1 = N_DOF1;
      n2 = N_DOF2;
    }

    /** map the given local degrees of freedom,
        coarse cell has lower number,
        0 is coarser side 0; 1,2 are on finer side 1 */
    void MapCoarseFine(int *Global, int I_K0, int I_K1, int I_K2,
             int *Indices0, int *Indices1, int *Indices2,
             int LocEdge0, int LocEdge1, int LocEdge2,
             TFEDesc2D *Desc0, TFEDesc2D *Desc1, TFEDesc2D *Desc2,
             int &Counter, int LowerFirstChild,
             TVector<THangingNode *> *vect,
             TVector<int> *numbers);

    /** map the given local degrees of freedom,
        coarse cell has bigger number,
        0 is coarser side 0; 1,2 are on finer side 1 */
    void MapFineCoarse(int *Global, int I_K0, int I_K1, int I_K2,
             int *Indices0, int *Indices1, int *Indices2,
             int LocEdge0, int LocEdge1, int LocEdge2,
             TFEDesc2D *Desc0, TFEDesc2D *Desc1, TFEDesc2D *Desc2,
             int &Counter, int LowerFirstChild,
             TVector<THangingNode *> *vect,
             TVector<int> *numbers);

};

#endif
