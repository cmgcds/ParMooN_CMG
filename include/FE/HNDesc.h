// =======================================================================
// @(#)HNDesc.h        1.1 10/30/98
//
// Class:       THNDesc
// Purpose:     describe fixed information for a hanging node
//
// Author:      Gunar Matthies  18.11.97
//
// =======================================================================

#ifndef __HNDESC__
#define __HNDESC__

/** describe fixed information for a hanging node */
class THNDesc
{
  protected:
    /** number of degrees in coupling */
    int N_Nodes;

    /** coefficient of other nodes in this coupling */
    double *Coeff;

  public:
    /** constructor, filling all data */
    THNDesc(int n_nodes, double *coeff);

    // Methods
    /** return number of nodes in coupling */
    int GetN_Nodes() 
    { return N_Nodes; }

    /** return coefficients of degrees of freedom in coupling */
    double *GetCoeff()
    { return Coeff; }

};

#endif
