// =======================================================================
// %W% %G%
//
// Class:       TMultiGrid3D
// Purpose:     store all data for a multi grid method
//
// Author:      Gunar Matthies 26.06.2000
//
// History:     26.06.2000 start of implementation
//
// =======================================================================

#ifndef __MULTIGRID3D__
#define __MULTIGRID3D__

#include <MGLevel3D.h>
#ifdef _MPI   
   #ifdef __3D__
    #include <ParFECommunicator3D.h>
   #else
    #include <ParFECommunicator2D.h>
   #endif
#endif 
// #define MAXN_LEVELS 25

class TMultiGrid3D
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
    TMGLevel3D *MultiGridLevels[MAXN_LEVELS];

    /** array of FE spaces */
    TFESpace3D *FESpaces[MAXN_LEVELS];

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
    TMultiGrid3D(int n_problems, int n_parameters, double *parameters);

    /** return number of multi grid levels */
    int GetN_Levels()
    { return N_Levels; }

    /** add new level as finest */
    void AddLevel(TMGLevel3D *MGLevel);

    /** add new level as finest */
    void ReplaceLevel(int i,TMGLevel3D *MGLevel);

    /** return i-th level as TMGLevel object */
    TMGLevel3D *GetLevel(int i)
    { return MultiGridLevels[i]; }

    /** restrict u1, u2 from finest grid to all coarser grids */
    void RestrictToAllGrids();

    /** set correct values for Dirichlet nodes on grid i */
    void SetDirichletNodes(int i);

    /** cycle on level i */
    void Cycle(int i, double &res);

    /** set recursion for multigrid */ 
    void SetRecursion(int levels);
    
    /** Smoother Cycles are called here */
    void Smooth(int smoother_type, TMGLevel3D *Level, 
#ifdef _MPI
		TParFECommunicator3D *ParComm, 
#endif
		double &oldres);
};

#endif
