// =======================================================================
// %W% %G%
//
// Class:       TFE3DMapper
// Purpose:     find out which of the given local degress of freedom
//              are equivalent to the same global degree of freedom
//
// Author:      Gunar Matthies  05.11.97
//
// History:     Added edge details by Sashi (08.03.12)
// =======================================================================

#ifndef __FE3DMAPPER__
#define __FE3DMAPPER__

#define HANGINGNODE -5
#define FIRSTMARK -10

#include <FEDesc3D.h>
#include "Vector.h"

/** find out which of the given local degress of freedom,
    are equivalent to the same global degree of freedom */
class TFE3DMapper
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

    /** which pairs of local degrees matching in an array[][][2] */
    int **Pairs;

    /** number of dof without an opposite */
    int N_NoOpposite;

    /** list of dof without an opposite */
    int **NoOpposite;

    /** memory for internal storage */
    int *Aux;

  public:
    /** constructor, filling all data */
    TFE3DMapper(char *name, char *description, int n0, int n1, 
                int n_pairs, int **pairs, 
                int n_noopposite, int **noopposite,
                int n_nodes);

    // Methods
    /** return name of mapper */
    char *GetName() const
    { return Name; }

    /** return description of mapper */
    char *GetDescription() const
    { return Description; }

    /** return number of degrees on side 0 */
    int GetN_DOF0() const
    { return N_DOF0; }

    /** return number of degrees on side 1 */
    int GetN_DOF1() const
    { return N_DOF1; }

    /** return number of degrees on both sides */
    void GetN_DOF(int &n0, int &n1) const
    {
      n0 = N_DOF1;
      n1 = N_DOF1;
    }

    /** return number of pairs */
    int GetN_Pairs() const
    { return N_Pairs; }

    /** return pairs of matching dof's */
    int **GetPairs() const
    { return Pairs; }

    /** return a pair of matching dof */
    int *GetPairs(int i) const
    { return Pairs[i]; }

    /** map the given local degrees of freedom */
    void Map(int Type,
             int *Global, int I_K0, int I_K1, 
             int *Indices0, int *Indices1,
             int LocEdge0, int LocEdge1,
             TFEDesc3D *Desc0, TFEDesc3D *Desc1,
             int &Counter,
             TVector<THangingNode *> *vect,
             TVector<int> *numbers);

    /** "map" the given dof on a boundary joint */
    void MapBound(int *Global, int I_K, int *Indices, 
                  int &BoundCounter,
                  TVector<THangingNode *> *vect,
                  TVector<int> *numbers);

    /** map the two given degrees of freedom */
    void MapDOF(int *Global, int dof0, int dof1, int &Counter,
                TVector<THangingNode *> *vect,
                TVector<int> *numbers);

    /** "map" the given dof on a boundary edge */
    void MapBoundEdge(int N_EdgeDOF, int *Global, int I_K, int *Indices, 
                      int &BoundCounter,
                      TVector<THangingNode *> *vect,
                      TVector<int> *numbers);   
    
    /** "map" the given dof on a boundary vert */    
    void MapBoundVert(int *Global, int I_K, int Index, 
                     int &BoundCounter,
                     TVector<THangingNode *> *vect,
                     TVector<int> *numbers);  
    
    
    /** destructor */
    ~TFE3DMapper()
    { delete Aux; }

};

#endif
