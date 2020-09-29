/** ************************************************************************ 
*
* @class     TSystem_6DOF
* @brief     stores the information of 6DOF system matrix 
* @author    Sashikumaar Ganesan  
* @date      13.5.2017
* @History    
 ************************************************************************  */

#ifndef __SYSTEM_6DOF__
#define __SYSTEM_6DOF__

#include <SystemTNSE3D_ALE.h>

/**class for 3D  TNSE system coupled with 6dof */
class TSystem_6DOF 
{
  protected:
   
    /** grid fespace */
    TFESpace3D *GridFESpace; 
    
    /** posistion vector of the mesh */
    TFEVectFunct3D *GridPosFEVect;

    /** No. of Grid DOFs */
    int N_GridDOFs, N_GridActive;  
    
    double *gridpos;    
    
    /** coordinates of the mesh WRT centre of mass */
    double *gridpos_CM;
    double P_Axis[9];     //  principle axes
    double P_Axis_RotVect[3], Rot_Vec[3];
    
    /** centre of gravity */
    double CGx_Body, CGy_Body, CGz_Body, Ixx, Iyy, Izz, BodyMass;
    
    /** body force */
    double TotalBodyForce[6];
    
    /** sol */
    double rhs[12], B[12];
    double sol[12], oldsol[12];
    
  public:
    /** constructor */
     TSystem_6DOF(TFESpace3D *gridFESpace, TFEVectFunct3D *gridPosFEVect);
     
    /** destrcutor */
    ~TSystem_6DOF();

    /** methods */
    
    /** Initilize the discrete forms and the matrices */    
    void Init6DOF();
    
    /** Solve and update the 6dof eqn using RK4 */
    void SolveAndUpdate(double tau);
    
    void ProjectVect(double *P_Axis_RotVect);
    
    void AssembleRhs();
 
    void GetMOI(double *MOI_Tensor);    
    
};

#endif
