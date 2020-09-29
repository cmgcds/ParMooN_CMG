// =======================================================================
// %W% %G%
// 
// Class:       TNodalFunctional3D
// Purpose:     realize nodal functionals in 3D
//
// Author:      Gunar Matthies (30.05.00)
//
// History:     start of implementation 30.05.00 (Gunar Matthies)
//
// =======================================================================

#ifndef __NODALFUNCTIONAL3D__
#define __NODALFUNCTIONAL3D__

#include <Enumerations.h>
#include <Constants.h>
#include <Collection.h>

/**  realize nodal functionals in 3D */
class TNodalFunctional3D
{
  protected:
    /** number of all functionals */
    int N_AllFunctionals;

    /** array of number of functionals on one face */
    int *N_FaceFunctionals;

    /** number of points needed for all nodal functionals */
    int N_PointsAll;

    /** values at points with following xi-coordinates are needed
        to evaluate all nodal functionals */
    double *Xi;

    /** values at points with following eta-coordinates are needed
        to evaluate all nodal functionals */
    double *Eta;

    /** values at points with following zeta-coordinates are needed
        to evaluate all nodal functionals */
    double *Zeta;

    /** routine for evaluating all functionals */
    EvalAllNF *EvalAll;

    /** array of numbers of points needed for face nodal functionals */
    int *N_PointsFace;

    /** xi-coordinate of points for evaluating functional on face */
    double **XiArray;

    /** eta-coordinate of points for evaluating functional on face */
    double **EtaArray;

    /** zeta-coordinate of points for evaluating functional on face */
    double **ZetaArray;

    /** t-parameter for points for evaluating functional on face */
    double *T;

    /** s-parameter for points for evaluating functional on face */
    double *S;

    /** routine for evaluating the face functionals */
    EvalJointNF *EvalFace;

    /** ID for this set of nodal functionals */
    NodalFunctional3D ID;

  public:
    /** constructor */
    TNodalFunctional3D(NodalFunctional3D id,
                       int n_allfunctionals, int *n_facefunctionals,
                       int n_pointsall, int *n_pointsface,
                       double *xi, double *eta, double *zeta,
                       double **xiarray, double **etaarray,
                       double **zetaarray,
                       double *t, double *s,
                       EvalAllNF *evalall,
                       EvalJointNF *evalface);

    /** return information for points for all functionals */
    void GetPointsForAll(int &n_points, double* &xi, double* &eta,
                         double* &zeta);

    /** return information for points for face functionals 
        on joint j */
    void GetPointsForFace(int j, int &n_points, double* &xi, double* &eta,
                          double* &zeta);

    /** return information for points for face functionals */
    void GetPointsForFace(int &n_points, double* &t, double* &s);

    /** return values for all nodal functionals */
    void GetAllFunctionals(TCollection *Coll, TBaseCell *Cell,
                           double *PointValues, double *Functionals)
    { EvalAll(Coll, Cell, PointValues, Functionals); }

    /** return values for face nodal functional */
    void GetFaceFunctionals(TCollection *Coll, TBaseCell *Cell, int Joint,
                            double *PointValues, double *Functionals)
    { EvalFace(Coll, Cell, Joint, PointValues, Functionals); }

    /** return ID for this set */
    NodalFunctional3D GetID() const
    { return ID; }
    
    int n_functionals() const
    {return N_AllFunctionals; }

    int n_face_functionals(int face) const
    {return N_FaceFunctionals[face]; }

    int n_points() const
    { return N_PointsAll; }

};

#endif
