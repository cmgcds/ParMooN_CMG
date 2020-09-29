// =======================================================================
// 
// Class:       TAuxParam2D3D
// Purpose:     store parameter functions and FE functions
//
// Author:      Andreas Hahn (19.08.2010)
//
// History:     start of implementation 19.08.2010 (Andreas Hahn)
//
// ======================================================================= 

#ifndef __AUXPARAM2D3D__
#define __AUXPARAM2D3D__

#include <FEVectFunct3D.h>
#include <Enumerations.h>
#include <BaseCell.h>

class TAuxParam2D3D
{
  protected:
    /** stored FEFunctions3D */
    int mN_FEFunctions;
    TFEFunction3D **mFEFunctions;
    
    /** stored FEVectFuncts3D */
    int mN_FEVectFuncts;
    TFEVectFunct3D **mFEVectFuncts;
    
    /** total number of parameters */
    int mN_Parameters;
    int mN_TotalParameters;
    
    int *mFct_Index;
    
    MultiIndex3D *mDerivatives;
    
    int *mIsGlobal;
    
    /** internal storage */
    double n1[MaxN_QuadPoints_3D];
    double n2[MaxN_QuadPoints_3D];
    double n3[MaxN_QuadPoints_3D];
    
    TBaseCell *mCell;
    int mCellNr;
    int mGlobCellNr;
    int mN_Points;
    bool mProject;
    
  protected:
    void GetFunctionsParams(int Fct, MultiIndex3D Der, int index, double **parameters);
    void GetVectFunctParams(int Fct, int index, double **parameters);
    
  public:
    TAuxParam2D3D (int n_fefunctions, TFEFunction3D **fefunctions,
		   int n_fevectfuncts, TFEVectFunct3D **fevectfuncts,
		   int n_parameters, int *fct_index, MultiIndex3D *derivatives,
		   int *isglobal);
		   
    ~TAuxParam2D3D ();
    
    void GetParameters(TBaseCell *cell, int cellnr, int jointnr, int globnr,
		      RefTrans3D RefTrans,
		      int n_points, double *xi, double *eta,
		      double **parameters);
		      
    int GetN_Parameters()
    { return mN_TotalParameters; }
};

#endif