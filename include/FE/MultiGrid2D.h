// =======================================================================
// @(#)MultiGrid2D.h        1.5 05/05/00
//
// Class:       TMultiGrid2D
// Purpose:     store all data for a multi grid method
//
// Author:      Gunar Matthies 02.11.1998
//
// History:     02.11.1998 start of implementation
//
// =======================================================================
#ifdef __2D__

#ifndef __MULTIGRID2D__
#define __MULTIGRID2D__

#include <MGLevel2D.h>

// #define MAXN_LEVELS 100

class TMultiGrid2D
{
  protected:
    /** number of levels */
    int N_Levels;

    /** number of problems */
    int N_Problems;

    /** number of parameters */
    int N_Parameters;

    /** array of double parameters */
    double *Parameters;

    /** array of multi grid levels */
    TMGLevel2D *MultiGridLevels[MAXN_LEVELS];

    /** array of FE spaces */
    TFESpace2D *FESpaces[MAXN_LEVELS];

    /** array of function vectors on each level */
    double **FunctionVectors[MAXN_LEVELS];

    /** right-hand side vectors */
    double **RhsVectors[MAXN_LEVELS];

    /** auxiliary vectors */
    double **AuxVectors[MAXN_LEVELS];

    /** number of recursions */
    int mg_recursions[MAXN_LEVELS];

  public:
    /** constructor */
    TMultiGrid2D(int n_problems, int n_parameters, double *parameters);

    /** return number of multi grid levels */
    int GetN_Levels()
    { return N_Levels; }

    /** add new level as finest */
    void AddLevel(TMGLevel2D *MGLevel);

    /** add new level as finest */
    void ReplaceLevel(int i,TMGLevel2D *MGLevel);

    /** return i-th level as TMGLevel object */
    TMGLevel2D *GetLevel(int i)
    { return MultiGridLevels[i]; }

    /** restrict u1, u2 from finest grid to all coarser grids */
    void RestrictToAllGrids();

    /** set correct values for Dirichlet nodes on grid i */
    void SetDirichletNodes(int i);

    /** cycle on level i */
    void Cycle(int i, double &res);

    /** set recursion for multigrid */ 
    void SetRecursion(int levels);
};

#endif

#endif // #ifdef __2D__