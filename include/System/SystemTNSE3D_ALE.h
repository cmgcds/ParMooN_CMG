/** ************************************************************************ 
* @brief     source file for TSystemTNSE3D_ALE
* @author    Sashikumaar Ganesan
* @date       6.4.2016
* @History    MoveBDwithVelo added 23.11.2017
 ************************************************************************  */

#ifndef __SYSTEMTNSE3D_ALE__
#define __SYSTEMTNSE3D_ALE__

#include <SystemNSE3D.h>
#ifdef _SMPI
#include <SeqParDirectSolver.h>
#endif
/**class for 3D  TNSE system matrix */
class TSystemTNSE3D_ALE : public TSystemNSE3D
{
  protected:
    
    /** vms projection fespace */
    TFESpace3D **Projection_Spaces;
    
    /** M - mass/system mat for TNSE velocity component   */
    TSquareMatrix3D **SqmatrixM11, **SqmatrixM12, **SqmatrixM13, 
                    **SqmatrixM21, **SqmatrixM22, **SqmatrixM23,
                    **SqmatrixM31, **SqmatrixM32, **SqmatrixM33;

#ifdef __PRIVATE__  
    /** sqstructureG of the  vms projection matrix */
    TSquareStructure3D *sqstructureL;

    /** structure of the vms projection  matrix */
    TStructure3D *structure_G, *structure_tilde_G; 
    
    /** G -  mat for VMS   */    
    TMatrix3D *matrix_tilde_G11, 
              *matrix_tilde_G22,
	      *matrix_tilde_G33,
              *matrix_G11,
	      *matrix_G22,
	      *matrix_G33;
	      
    TMatrix3D **Matrices_tilde_G11,
              **Matrices_tilde_G22,
	      **Matrices_tilde_G33, 
	      **Matrices_G11,
	      **Matrices_G22,
	      **Matrices_G33;  

    /** L -  mat for VMS   */
    TSquareMatrix3D *sqmatrixL, **MatricesL;
#endif
     
    // Assembling rhs*/
    TAssembleMat3D *RhsOnlyAssemble;
    
#ifdef _SMPI
   TSeqParDirectSolver* P_DS;
#endif
    
     /** working rhs, used in AssembleSystMat() */
    double *B;   
   
    /** to store defect */
    double *defect;   
    
    /** factor that multplied with Mat A in working rhs */
    double gamma;      

    /** Discrete form of the M and rhs matrics */
    TDiscreteForm3D *DiscreteFormRhs; 
    
    /** NSE_Rhsaux is used to for assembling rhs only*/
    TAuxParam3D *NSE_Rhsaux;
    
    /** Systmat assemble indicator */
    bool SystMatAssembled;
    
    
    /** needed for error calculation in time */
    double olderror_l_2_l_2u;
    
      
    //-------------------------used for ALE-------------------------------
     /** No. of Grid DOFs */
    int N_GridDOFs, N_GridActive, *GridKCol, *GridRowPtr;
    
    TFEVectFunct3D **MeshVelocity;
    
    /** Grid Posistions */
    double **MeshVelo, *gridpos, *gridpos_old, *gridpos_ref,  *GridRhs, *Entries[9],**Mesh_restrict_aux;   
    
    /** grid fespace */
    TFESpace3D **GridFESpace;    
    
     
    /** mesh matrices assembling */
    TAssembleMat3D *MeshMatAssemble;    
    
    /** grid pos vector */
    TFEVectFunct3D *GridPos, *RefGridPos;
     
//     /** Fe functions of NSE */
//     TFEFunction3D *MeshVeloFct[3];    
    
    /** Discrete form for moving mesh */
    TDiscreteForm3D *DiscreteFormMARhs, *DiscreteFormGrid;  
           
    /** marices for the moving grid **/
    TSquareMatrix3D *SqmatrixG11, *SqmatrixG12, *SqmatrixG13,
                    *SqmatrixG21, *SqmatrixG22, *SqmatrixG23,
		    *SqmatrixG31, *SqmatrixG32, *SqmatrixG33,
		    *SQMATRICES_GRID[9];
 
    /** structure for the moving grid **/
    TSquareStructure3D *SquareStructureG;
    
    /** aux for mesh */
    TAuxParam3D *Aux_ALE, *Meshaux;

    /** Grid bounadry conditions **/ 
    BoundCondFunct3D *GridBoundaryConditions[3];
    BoundValueFunct3D *GridBoundValues[3];
      
    /** method for Modify Mesh Coords */
    ModifyMeshCoords_3D *ModifyCoord;
    
     /** method for Modify Mesh Coords */
    ModifyBoundCoords_3D *ModifyBoudary;   
    MoveBound_3D  *MoveBDwithVelo;
    
    
    /** */
    bool SolveLinearElastic, CONSERVATIVEALE, BdMoveWithVelo;
    
    /** info of moving boundaries */
    int N_MovVert, N_Movfaces, *Movfaces;
    TVertex **MovBoundVert;
    TBaseCell ** MovCells;
    

    /**mesh mat LU decomp values */
    void *MeshMatSolver_Symbolic, *MeshMatSolver_Numeric;
    double *MeshMatSolver_Values;
    int *MeshMatSolver_KCol, *MeshMatSolver_Row;
    //------------------------------------------------------------------------------
    
  public:
    /** constructor */
     TSystemTNSE3D_ALE(int N_levels, TFESpace3D **velocity_fespace, TFESpace3D **presssure_fespace, TFEVectFunct3D **velocity, 
                       TFEFunction3D **pressure, double **sol, double **rhs,  int disctype, int nsetype, int solver,
// #ifdef __PRIVATE__  
                      //  TFESpace3D **Projection_space,
// #endif    
                       TFESpace3D ** gridFESpace, TFEVectFunct3D **meshVelocity, bool conservativeale);
     
    /** destrcutor */
    ~TSystemTNSE3D_ALE();

    /** methods */
    
    /** Initilize the discrete forms and the matrices */    
    void Init(CoeffFct3D *lincoeffs, BoundCondFunct3D *BoundCond, BoundValueFunct3D *U1BoundValue, 
                             BoundValueFunct3D *U2BoundValue, BoundValueFunct3D *U3BoundValue, BoundCondFunct3D *GridBoundCondition,
                             BoundValueFunct3D *GridBoundValue_x,BoundValueFunct3D *GridBoundValue_y,BoundValueFunct3D *GridBoundValue_z,CoeffFct3D *GridBilinearCoeffs);
    
    void AddMeshModifyFunction(ModifyMeshCoords_3D *modifyCoord)
       {  ModifyCoord = modifyCoord; SolveLinearElastic = FALSE; }
  
    /** add info of moving boundaries */
    void AddBoundModifyFunction(ModifyBoundCoords_3D *modifyboudary, int n_movVert, TVertex **movboundVert,
                                int n_movfaces, int *movfaces, TBaseCell ** movcells);

    void AddMoveBoundFunction(MoveBound_3D *movebdwithvelo, int n_movVert, TVertex **movboundVert,
                                int n_movfaces, int *movfaces, TBaseCell ** movcells);
    
    /** assemble the M, A and rhs */
    void Assemble();
        
    /** assemble only the rhs of NSE system */
    void AssembleRhs(); 
    
    /** Get Mesh Velo */
    void GetMeshVelo(bool MoveMesh, int solver_flag);
    
    /** Get Mesh Velo */
    void GetMeshVeloAndMove(double Currtime, double tau);
    
    /** assemble the Mesh mat and rhs */ 
    void AssembleMeshMat();
    
    /** scale B matices and assemble rhs based on the \theta scheme  */
    void AssembleSystMat(double scale, double *oldrhs, double *rhs, double *sol);

    /** scale B matices and assemble rhs based on the \theta scheme  */
    void AssembleSystMatNonLinear();
    
     /** restoring the mass matrix */
    void RestoreMassMat();
    
    /** restoring the mass matrix */
    void RestoreMassMatNonLinear();
      
    /** assemble the nonlinear part of the NSE system */
    void AssembleNonLinear();

    /** solve the system matrix */
    void  Solve(double *sol);
    
    /** get the resudual of the NSE system */
    void GetResidual(double *sol, double &impuls_residual, double &residual);
 
    /** measure the error in the NSE */
    void MeasureTNSEErrors(DoubleFunct3D *ExactU1, DoubleFunct3D *ExactU2, DoubleFunct3D *ExactU3, DoubleFunct3D *ExactP, double *AllErrors);

    /** print the matrices in a file */
    void printall_matrix();
    
    void GetMesh();

    double* getRhs()
    {return B;}
    
    TFEVectFunct3D *GetGridPosFeVect()
    { return GridPos; }
    
    void All_levels_check();
    
    double BlockMatVect(TSquareMatrix *A);
    
    double BlockMatVect(TMatrix *A);

    void CheckAllMat();

    void DOF_stats(TMatrix3D* MatB, char M , int k,char * name);

    void RHS_stats();
    
    void SetSystMatAssembled(bool val)
     {SystMatAssembled=val;}


    // Thivin -- To assign Value to MeshVelocity
    void AssignMeshVelocity_parameters(double* val,int comp)
    {
      MeshVelo[comp] =val;
    }

    // Thivin - To pick up Free Slip Boundaries to impose u.n = 0
    void pickDOFsOfFreeSlipBoundaries(TCollection* coll, TFESpace3D* gridfespace, std::vector<int> freeSlipBoundIds,std::vector<int> boundIds);

    // THIVIN - To calculate the normals of the Free Slip Boundaries to impose u.n = 0
    void get_surface_normals_slipBoundary(TCollection* coll,TFEVectFunct3D* Velo_FEvect );

    // THIVIN - Declararion of variables for Slip Boundary COndition
    int N_bd_FreeSlip_Vertex_u ;
    int N_bd_FreeSlip_Cells_u ;
    int N_bd_FreeSlip_Joints_u ;

    //THIVIN - Declaration of vector Structures for Slip
    std::vector<int> Bd_FreeSlip_Vertex_u;
    std::vector<int> Bd_FreeSlip_VertexLocal_u;
    std::vector<int> Bd_FreeSlip_Cells_u;
    std::vector<int> Bd_FreeSlip_Joints_u;

    std::vector<double> Bd_normal_1_u;
    std::vector<double> Bd_normal_2_u;
    std::vector<double> Bd_normal_3_u;

    std::vector<double> Bd_TangentA_1_u;
    std::vector<double> Bd_TangentA_2_u;
    std::vector<double> Bd_TangentA_3_u;

    std::vector<double> Bd_TangentB_1_u;
    std::vector<double> Bd_TangentB_2_u;
    std::vector<double> Bd_TangentB_3_u;

    // For Edge Vectors
    void get_surface_normals_slipBoundary_EdgeNodes(TCollection* coll,TFEVectFunct3D* MeshVelo_FEvect );

        std::vector<int> Bd_EdgeFreeSlip_Vertex_u;
        std::vector<int> Bd_EdgeFreeSlip_Cells_u;
        std::vector<int> Bd_EdgeFreeSlip_Joints1_u;   
        std::vector<int> Bd_EdgeFreeSlip_Joints2_u;
        std::vector<int> Bd_EdgeFreeSlip_VertexLocal_u;
        
        int N_bd_EdgeFreeSlip_Vertex_u;
        int N_bd_EdgeFreeSlip_Cells_u;
        int N_bd_EdgeFreeSlip_Joints_u;


        std::vector<double> Bd_edge_normalA_1_u;
        std::vector<double> Bd_edge_normalA_2_u;
        std::vector<double> Bd_edge_normalA_3_u;

        std::vector<double> Bd_edge_normalB_1_u;
        std::vector<double> Bd_edge_normalB_2_u;
        std::vector<double> Bd_edge_normalB_3_u;

        std::vector<double> Bd_edge_TangentA_1_u;
        std::vector<double> Bd_edge_TangentA_2_u;
        std::vector<double> Bd_edge_TangentA_3_u;

        void remove_Redundant_Dirichlet_DOF(TSquareMatrix3D *sqmatrices0,TSquareMatrix3D *sqmatrices1,
                                          TSquareMatrix3D *sqmatrices2,
                                        TSquareMatrix3D *sqmatrices3,TSquareMatrix3D *sqmatrices4,
                                        TSquareMatrix3D *sqmatrices5,TSquareMatrix3D *sqmatrices6,
                                        TSquareMatrix3D *sqmatrices7,
                                        TSquareMatrix3D *sqmatrices8 );

        // THIVIN
        // if the FLag is set to one, then the direct Solver without removing the redundant DOF's will be called
        int directSolverwithoutRemoveRedundant_flag = 0;

        void impose_FreeSlip_BoundaryCondition(double* rhs,int length, int N_ActiveBound);
        // void impose_FreeSlip_BoundaryCondition_rectMatrix(double* rhs,int length, int N_ActiveBound);


        //THIVIN
        // This is assumed for slip boundariies which are perpendicular to any one of the axis. 
        // Based axis it is perpendicular to , the value zero will be imposed based on the DOF NUmber...
        // ... rather than manipulating the system matrix
        // ** NOTE ** , it should be called only after the Solve Routine
        void Strict_impose_FreeSlip_BoundaryCondition(double* rhs,int length, int N_ActiveBound);
        
        //THIVIN
        // Function is used for imposing external boundary condition on the fluid 
        // like tilting or swaying motion
        void imposeExternalBoundaryCondition(void externalBoundaryParameters(double&, unsigned int&, double&), TFEVectFunct3D* externalVelocityFEvect, TFEVectFunct3D* VelocityFEvect);

};

#endif
