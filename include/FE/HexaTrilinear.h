// =======================================================================
// @(#)HexaTrilinear.h        1.3 02/22/00
//
// Class:      THexaTrilinear
//
// Purpose:    trilinear reference transformations for Hexahedron
//
// Author:     Daniel Quoos  
//
// History:    08.07.97 start implementation
// 
// =======================================================================

#ifndef __HexaTrilinear__
#define __HexaTrilinear__

#include <Enumerations.h>
#include <RefTrans3D.h>

/** reference transformations for Hexahedron */
class THexaTrilinear : public TRefTrans3D
{
  protected:
    /** x coordinate */
    double x0, x1, x2, x3, x4, x5, x6, x7;

    /** y coordinate */
    double y0, y1, y2, y3, y4, y5, y6, y7;

    /** z coordinate */
    double z0, z1, z2, z3, z4, z5, z6, z7;

    /** x parameters for reference transformation */
    double xc0, xc1, xc2, xc3, xc4, xc5, xc6, xc7;

    /** y parameters for reference transformation */
    double yc0, yc1, yc2, yc3, yc4, yc5, yc6, yc7;

    /** z parameters for reference transformation */
    double zc0, zc1, zc2, zc3, zc4, zc5, zc6, zc7;

    /** detjk */
    double detjk;

    /** 1/detjk */
    double rec_detjk;

  public:
    /** constuctor */
    THexaTrilinear();

    /*** THIVIN : Get the Normal Vectors of the given Reference Co ordinates in a particular Joint ****/
    /**** NormaliseFlag - 0  - Returns the normal vectors without Normalising *****/
    /***  NormaliseFlag - 1  - Returns the normal vectors with Normalising  ****/
    void GetNormalVectors(int JointNumber1,double xi,double eta,double zeta,double &n1, double&n2,double &n3,bool NormalizeFlag);

    /** transfer from reference element to original element */
    void GetOrigFromRef(double eta, double xi, double zeta, double &x, double &y, double &z);

    /** transfer a set of points form reference to original element */
    void GetOrigFromRef(int N_Points, double *eta, double *xi, double *zeta, 
                        double *x, double *y, double *z, double *absdetjk);

    /** transfer form reference element to original element */
    void GetOrigFromRef(double *ref, double *orig);

    /** transfer from original element to reference element */
    void GetRefFromOrig(double x, double y, double z, double &eta, double &xi, double &zeta);

    /** transfer from original element to reference element */
    void GetRefFromOrig(double *orig, double *ref);

    /** calculate functions and derivatives from reference element
        to original element */
    void GetOrigValues(BaseFunct3D BaseFunct,
                       int N_Points, double *xi, double *eta, double *zeta,
                       int N_Functs, QuadFormula3D HexaFormula);

    /** calculate functions and derivatives from reference element
        to original element, for all given elements */
    void GetOrigValues(int N_Sets, BaseFunct3D *BaseFunct,
                       int N_Points, double *xi, double *eta, double *zeta,
                       QuadFormula3D HexaFormula,
                       bool *Needs2ndDer);

    /** calculate functions and derivatives from reference element
        to original element */
    void GetOrigValues(double xi, double eta, double zeta, int N_BaseFunct,
                double *uref, double *uxiref, double *uetaref, double *uzetaref,
                double *uorig, double *uxorig, double *uyorig, double *uzorig,
                int _BaseVectDim = 1);

    /** calculate functions and derivatives from reference element
        to original element  - only the BASIS Function Values , not the Gradients */
    void GetOrigValues(double xi, double eta, double zeta,
                               int N_BaseFunct,
                               double *uref,
                               double *uorig,
                               int _BaseVectDim);

    /** set element to cell */
    void SetCell(TBaseCell * cell);

    /*** THIVIN : Get the Ref Co-ordinates of the Hexaheadral cell based on the local Node Numbers ****/    
    /**** The REf Co ordinates can be obtained even without Setting the values of the cell using SetCell()  */
    void GetRefvaluesfromLocalNodeNumber(FE3D FEid, int local_Node_num,double& xi, double& eta, double& zeta);


    /*** THIVIN : Get the Ref Co-ordinates of the Joint DOF of the Hexaheadral cell based on the  Xi,eta,zeta Values ****/
    void GetRefValuesfromJointid(int Jointno,double xi, double eta, double zeta , double& xi_1 , double& xi_2);

    

    /** return outer normal unit vector */
    void GetOuterNormal(int j, double s, double t,
                        double &n1, double &n2, double &n3);

    /** return two tangent vectors */
    void GetTangentVectors(int j, double p1, double p2,
        double &t11, double &t12, double &t13,
        double &t21, double &t22, double &t23);
    
    /** @brief Piola transformation for vector valued basis functions **/
    void PiolaMapOrigFromRef(int N_Functs, double *refD000, double *origD000);
};

#endif
