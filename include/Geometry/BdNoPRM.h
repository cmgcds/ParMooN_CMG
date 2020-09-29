// =======================================================================
//
// Class:       TBdNoPRM
// Superclass:  TBoundComp
// Purpose:     3D domain without PRM file
//
// Author:      Volker John 2008/01/28
//
// =======================================================================

#ifndef __BDNOPRM__
#define __BDNOPRM__

#include <BoundComp3D.h>

/** no PRM file available */
class TBdNoPRM : public TBoundComp3D
{
  protected:

  public:
    // Constructor
    TBdNoPRM(int id);

    virtual ~TBdNoPRM() {};
    
    // Methods
    /** set all parameters to the given values */
    void SetParams ();

    /** return the coordinates of parameter value T, S */
    virtual int GetXYZofTS(double T, double S,
                           double &X, double &Y, double &Z);

    /** return the parameter value T, S of coordinates */
    virtual int GetTSofXYZ(double X, double Y, double Z,
                           double &T, double &S);

    /** return parameters and coordinates of a given linear
        combination of vertices */
    virtual int GetXYZandTS(int N_Points, double *LinComb,
                            double *xp, double *yp, double *zp,
                            double *tp, double *sp,
                            double &X, double &Y, double &Z,
                            double &T, double &S);

    /** read parameter from input file */
    virtual int ReadIn(std::ifstream &dat);
};

#endif
