/** ************************************************************************ 
*
* @class     TCDSystemTimeDG_1
* @brief     stores the information of a 2D scalar dG(1) in time discretization
* @author    Sashikumaar Ganesan, 
* @date      24.12.15
* @History    
 ************************************************************************  */


#ifndef __CDSYSTEMTIMEDG_1__
#define __CDSYSTEMTIMEDG_1__

#include <SquareMatrix2D.h>
#include <CDSystemTimeDG.h>

/**class for 2D scalar system  dG(1) in time discretization */
 
class TCDSystemTimeDG_1  : public TCDSystemTimeDG
{
  protected:
    /** dG type */
    int Type;    
    
  public:
    /** constructor */
     TCDSystemTimeDG_1(TSquareMatrix2D *mat_m, TSquareMatrix2D *mat_A);

    /** destrcutor */
    ~TCDSystemTimeDG_1();
    
    /** assemble the system matrix */
    virtual void AssembleSysMat(double *Mu_old, double *Rhs);
    
    /** assemble the system matrix */   
    virtual void AssembleALESysMat_Qp1(double *Mu_old, double *Rhs);
    
    
    /** solve dG system and return the Sol at end t^n system matrix */    
    virtual void SoveTimedG(double *Sol); 
    
};

#endif
