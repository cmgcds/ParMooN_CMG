// =======================================================================
// @(#)BoundPart.h        1.1 10/30/98
//
// Class:       TBoundPart
// Purpose:     a closed part of boundary
//
// Author:      Volker Behns  18.06.97
//
// =======================================================================

#ifndef __BOUNDPART__
#define __BOUNDPART__

#ifdef __2D__
  #include <BoundComp2D.h>
#else
  #include <BoundComp3D.h>
#endif

/** a closed part of boundary */
class TBoundPart
{
  protected:
    /** number of boundary components */
    int N_BdComps;
    /** components of this boundary part */
#ifdef __2D__
    TBoundComp2D **BdComps;
#else
    TBoundComp3D **BdComps;
#endif

  public:
    // Constuctor
    /** initialize the data structure */
    TBoundPart(int n_comps);

    // Methods
    /** return number of boundary components */
    int GetN_BdComps()
    { return N_BdComps; }

#ifdef __2D__
    /** set boundary component i */
    void SetBdComp(int i, TBoundComp2D *bdcomp)
    { BdComps[i] = bdcomp; }

    /** return boundary component i */
    TBoundComp2D *GetBdComp(int i)
    { return BdComps[i]; }
#else
    /** set boundary component i */
    void SetBdComp(int i, TBoundComp3D *bdcomp)
    { BdComps[i] = bdcomp; }

    /** return boundary component i */
    TBoundComp3D *GetBdComp(int i)
    { return BdComps[i]; }
#endif

    #ifndef __3D__
      /** return the coordinates of parameter value T */
      int GetXYofT(int comp, double T, double &X, double &Y)
      { return BdComps[comp]->GetXYofT(T, X, Y); }
      /** return the coordinates of parameter value T */
      int GetTofXY(int comp, double X, double Y, double &T)
      { return BdComps[comp]->GetTofXY(X, Y, T); }
    #else
      /** return the coordinates of parameter values T and S */
      int GetXYZofTS(int comp, double T, double S,
                     double &X, double &Y, double &Z)
      { return BdComps[comp]->GetXYZofTS(T, S, X, Y, Z); }
    #endif
};

#endif
