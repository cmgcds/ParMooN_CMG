// =======================================================================
// %W% %G%
//
// Class:       TFE3DMapper1Reg
// Purpose:     find out which of the given local degress of freedom
//              are equivalent to the same global degree of freedom
//              for 1 regular case
//
// Assumption:  all elements on the refined side are equal
//
// Author:      Gunar Matthies  17.07.2000
//
// History:     start of reimplementation 17.07.2000 (GM)
//
// =======================================================================

#ifndef __FE3DMAPPER1REG__
#define __FE3DMAPPER1REG__

#include <FE3DMapper.h>

/** find out which of the given local degress of freedom,
    are equivalent to the same global degree of freedom 
    for 1 regular case */
class TFE3DMapper1Reg : public TFE3DMapper
{
  protected:
    /** rule how to permute the local dofs due to twist index */
    int **TwistPermutation;

    /** select which set of pairs will be used, due to MapType */
    int *CurrentPairs;

    /** select which set of noopposite dof will be used, due to MapType */
    int *CurrentNoOp;

    /** number of hanging nodes */
    int N_Hanging;

    /** which local dofs are hanging nodes */
    int *Hanging;

    /** type of hanging nodes */
    HNDesc *HangingTypes;

    /** Nodes for coupling */
    int **Coupling;

  public:
    /** constructor, filling all data */
    TFE3DMapper1Reg(char *name, char *description, int nfine, int ncoarse,
              int n_pairs, int **pairs, 
              int n_noopposite, int **noopposite,
              int n_hanging, int *hanging,
              HNDesc *hangingtypes, int **coupling,
              int n_nodes, int **twistpermutation);

    /** destructor */
    ~TFE3DMapper1Reg()
    { }; 

    // Methods
    /** map the given local degrees of freedom,
        coarse cell has lower number,
        0 is coarser side 0; 1,2 are on finer side 1 */
    void Map(int *Global,
             int I_KF0, int I_KF1, int I_KF2, int I_KF3,
             int I_KC,
             int *IndicesF0, int *IndicesF1, int *IndicesF2,
             int *IndicesF3, int *IndicesC,
             int TwistIndexF0, int TwistIndexF1, int TwistIndexF2,
             int TwistIndexF3, int TwistIndexC,
             int DirichletBound,
             int &Counter,
             TVector<THangingNode *> *vect,
             TVector<int> *numbers);

};

#endif
