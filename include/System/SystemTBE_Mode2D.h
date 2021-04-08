/** ************************************************************************ 
*
* @class     TSystemTBE_Mode2D
* @brief     stores the information of 2D Burgers Mode system matrix 
* @author    Sashikumaar Ganesan, 
* @date      28.11.20
* @History    
 ************************************************************************  */


#ifndef __SYSTEMTBE_MODE2D__
#define __SYSTEMTBE_MODE2D__

#include <SystemTBE2D.h>

/**class for 2D Burgers system matrix */
class TSystemTBE_Mode2D : public TSystemTBE2D
{
  protected:
 
    /**velo mode function */
    TFEVectFunct2D *Velocity_Mode;
    
  public:
    /** constructor */
     TSystemTBE_Mode2D(TFESpace2D *velocity_fespace,TFEVectFunct2D *Velocity, double *sol, double *rhs, int disctype, int solver,
                       TFEVectFunct2D *Velo_Mode);

    /** destrcutor */
    ~TSystemTBE_Mode2D();
    
    /** Initilize the discrete forms and the matrices */
    void Init(CoeffFct2D *lincoeffs, BoundCondFunct2D *BoundCond, BoundValueFunct2D *U1BoundValue,  
              BoundValueFunct2D *U2BoundValue, TAuxParam2D *beaux, TAuxParam2D *beaux_error);

};

#endif
