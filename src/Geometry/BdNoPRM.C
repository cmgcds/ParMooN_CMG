// =======================================================================
//
// Class:       TBdNoPRM
// Purpose:     3D domain without PRM file
//
// Author:      Volker John 2008/01/28
//
// =======================================================================

#include <BdNoPRM.h>
#include <MooNMD_Io.h>

// Constructor
TBdNoPRM::TBdNoPRM(int id) : TBoundComp3D(id)
{
  Type = NoPRM;
}

// Methods
// there are no params
void TBdNoPRM::SetParams ()
{
    ;
}

// set dummy values
int TBdNoPRM::GetXYZofTS(double T, double S,
                        double &X, double &Y, double &Z)
{
    X = -4711;
    Y = -4711;
    Z = -4711;

  return 0;
}

/** return parameters and coordinates of a given linear
    combination of vertices */
// set dummy values
int TBdNoPRM::GetXYZandTS(int N_Points, double *LinComb,
                          double *xp, double *yp, double *zp,
                          double *tp, double *sp,
                          double &X, double &Y, double &Z,
                          double &T, double &S)
{
    int i;
    double t;

    X = Y = Z = 0.0;
    for (i=0;i<N_Points;i++)
    {
	t = LinComb[i]; 
	X += t * xp[i];
	Y += t * yp[i];
	Z += t * zp[i];
    }
    T = -4711;
    S = -4711;
    return 0;
}

// set dummu valus
int TBdNoPRM::GetTSofXYZ(double X, double Y, double Z,
                        double &T, double &S)
{
  double TS_aux;

  T = -4711;
  S = -4711;
  return 0;
}

// no readin
int TBdNoPRM::ReadIn(std::ifstream &dat)
{
  return 0;
}
