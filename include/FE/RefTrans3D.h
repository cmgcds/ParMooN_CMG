// =======================================================================
// %W% %G%
//
// Class:      TRefTrans3D
//
// Purpose:    reference transformations for 3D geometric objects
//
// Author:     Gunar Matthies
//
// History:    01.02.00 start implementation
// 
// =======================================================================

#ifndef __REFTRANS3D__
#define __REFTRANS3D__

#include <Constants.h>
#include <Enumerations.h>
#include <BaseCell.h>

/** reference transformations for 3D geometric objects */
class TRefTrans3D
{
  protected:
    TBaseCell *Cell;

  public:
    /** constuctor */
    TRefTrans3D();

    /** transfer form reference element to original element */
    void GetOrigFromRef(double xi, double eta, double zeta,
                        double &x, double &y, double &z);

    /** transfer form reference element to original element */
    void GetOrigFromRef(double *ref, double *orig);

    /** transfer from original element to reference element */
    void GetRefFromOrig(double x, double y, double z,
                        double &xi, double &eta, double &zeta);

    /** transfer from original element to reference element */
    void GetRefFromOrig(double *orig, double *ref);

    /** calculate functions and derivatives from reference element
        to original element */
    void GetOrigValues(TBaseCell *cell);

    /** set original element to cell */
    virtual void SetCell(TBaseCell *cell)
    {  Cell = cell; }

    /** return outer normal unit vector */
    void GetOuterNormal(int j, double s, double t,
                        double &n1, double &n2, double &n3);

    /** return two tangent vectors */
    void GetTangentVectors(int j, double p1, double p2,
        double &t11, double &t12, double &t13,
        double &t21, double &t22, double &t23);
    
    virtual void PiolaMapOrigFromRef(int N_Functs, double *refD00, 
                                     double *origD00)
    { ErrMsg(" Piola Map not defined for this element\n"); };

};

#endif
