// =======================================================================
// @(#)FESpace1D.h        1.1 10/30/98
// 
// Class:       TFESpace1D
// Purpose:     class for all 1D finite element spaces
//
// Author:      Gunar Matthies (03.11.97)
//
// History:     start of implementation 03.11.97 (Gunar Matthies)
//
//              start of reimplementation
//              09.10.98 (Gunar Matthies)
//
// =======================================================================

#ifndef __FESPACE1D__
#define __FESPACE1D__

#include <FESpace.h>
#include <FE1D.h>



typedef struct 
{
  /** Total No. Dof on the interface */
 int  N_NonMortaDofs;

 /** Map of NonMort Dof to Global Dof */
 int *GlobalDof_NonMortaEdgeDof;

 /** No. of Dof on each nonmortar edge */ 
 int *N_DofPerEdge;

 /**  NonMort Dof  local GlobalNo on the interface*/
 int *EdgeNonMotLocGlobalNo;

 /**  NonMort Dof  local GlobalNo on the interface*/
 int *EdgeNonMotGlobalNo;

 /**  NonMort Dof  local beginindex on the interface*/
 int *EdgeNonMotBeginIndex;
 
} TNonMortarData;


/** class for all 1D finite element spaces */
class TFESpace1D : public TFESpace
{
  protected:
    /** array containing the used elements */
    FE1D *UsedElements; 

    /** array with an element for each shape */
    FE1D *ElementForShape;

    /** array storing the fe for each element, if necessary */
    FE1D *AllElements;

   /** indices for mapping between Nodalfunctional/Nodal-interpolation point 
    in operator-splitting methods --- Sashikumaar Ganesan */
   int *IntlPtIndexOfPts;

   /** number of nodal interpolation points in this space */
   int N_GlobNodalIntPts;

  public:
    /** constructor */
    TFESpace1D(TCollection *coll, char *name, char *description);

    /** constructor for building a space with elements of order k */
    TFESpace1D(TCollection *coll, char *name, char *description, 
               int k);

    /** constructor for building a space with the given elements */
    TFESpace1D(TCollection *coll, char *name, char *description,
               FE1D *fes);

    /** destructor */
    ~TFESpace1D();

    /** find used elements */
    void FindUsedElements();

    /** construct space */
    void ConstructSpace();

    /** return identifiers of used elements */
    FE1D *GetUsedElements()
    { return UsedElements; }

    /** return the FE Id for element i, corresponding to cell */
    FE1D GetFE1D(int i, TBaseCell *cell);

    void SetIntlPtIndexOfPts(int *intlPtIndexOfPts)
     { IntlPtIndexOfPts = intlPtIndexOfPts; }

    int *GetIntlPtIndexOfPts()
     { return IntlPtIndexOfPts; }

    void SetN_RootNodalPts(int n_GlobNodalIntPts)
     { N_GlobNodalIntPts  = n_GlobNodalIntPts; }

    int GetN_RootNodalPts()
     { return N_GlobNodalIntPts; }

};

#endif
