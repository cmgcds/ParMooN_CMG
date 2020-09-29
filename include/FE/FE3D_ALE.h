/** ************************************************************************ 
* @brief     source file for FE3D_ALE
             

             Functions :
                1, Picks up the DOF's of Free Surface Nodes Based on the Surface ID provided
                2, Calculates the normals of the surface and used to Impose the Condition 
                   u.n = w.n on Surface
                3, Based on the Surface Velocity , it Calculate the interior point's Velocity using 
                   "Poisson Equation" or "LinearElasticity Equation" ( Not Implemented )
                
                4, Based on Mesh Velocity , it will also Calulate the mesh Displacements after each iteration and  
                   Move the Mesh Based on the Given Boundary Conditions


* @author  Thivin Anandh D
* @date      4-Sep-2019
* @History   Jan -8   - Primary Implementation
  @history   Jul 4    - Added modifications for solving the system using CuSparse Routines
                        -- Implemented CudaRefactor Solve
 ************************************************************************  */
#include <vector>
#include<FE3D.h>
#include<Domain.h>
#ifdef  _CUDA
    #include<cudaSparseLinearSolvers.h>
#endif  //_CUDA




#ifndef _FE3D_ALE_
#define _FE3D_ALE_

class FE3D_ALE{


    private:
        TCollection* coll;
        BoundCondFunct3D **GridBoundaryConditions;
	    BoundValueFunct3D **GridBoundValues;

    public:

        // Constructor
        FE3D_ALE(TCollection* _coll,	BoundCondFunct3D **_GridBoundaryConditions, BoundValueFunct3D **_GridBoundValues, 
					TFESpace3D* Meshfespace,TFEVectFunct3D* MeshVelocityVectFunction3D );


        TFESpace3D *fesp[1],*ferhs[3];

        TFESpace3D* fespace_alemeshMovement;
        TBaseCell** MovCells; 				
        TVertex **MovBoundVert;
        TJoint **MovJoints;
        
 
        // Matrix Data Structure

	TSquareMatrix3D *SqmatrixG11;
	TSquareMatrix3D *SqmatrixG12;
	TSquareMatrix3D *SqmatrixG13;
	TSquareMatrix3D *SqmatrixG21;
	TSquareMatrix3D *SqmatrixG22;
	TSquareMatrix3D *SqmatrixG23;
	TSquareMatrix3D *SqmatrixG31;
	TSquareMatrix3D *SqmatrixG32;
	TSquareMatrix3D *SqmatrixG33;

    TSquareStructure3D *sqstructure;
    TDiscreteForm3D* discreteform;

    MultiIndex3D AllDerivatives[4];
    TAuxParam3D *Meshaux;

	TSquareMatrix3D *SQMATRICES_GRID[9];
	double *RHS[3];
	double *Entries[9];
	int *GridKCol, *GridRowPtr;


        // For Free Slip Boundary cells on the Boundary Face
        std::vector<int> Bd_FreeSlip_Vertex;
        std::vector<int> Bd_FreeSlip_Cells;
        std::vector<int> Bd_FreeSlip_Joints;
        std::vector<int> Bd_FreeSlip_VertexLocal;

        std::vector<double> Bd_normal_1;
        std::vector<double> Bd_normal_2;
        std::vector<double> Bd_normal_3;



        std::vector<double> Bd_TangentA_1;
        std::vector<double> Bd_TangentA_2;
        std::vector<double> Bd_TangentA_3;

        std::vector<double> Bd_TangentB_1;
        std::vector<double> Bd_TangentB_2;
        std::vector<double> Bd_TangentB_3;

        int N_bd_FreeSlip_Vertex;
        int N_bd_FreeSlip_Cells;
        int N_bd_FreeSlip_Joints;

    // for storing the EDGE DOF's of the freeSurface boundary
        std::vector<int> FreeSurfaceEdgeDOFs;

        // For Free Slip Boundary cells on the Boundary edge ( intersection of two Boundaries )
        //  They Need to Store 2 Normals and 2 Joint Id's
        std::vector<int> Bd_EdgeFreeSlip_Vertex;
        std::vector<int> Bd_EdgeFreeSlip_Cells;
        std::vector<int> Bd_EdgeFreeSlip_Joints1;   
        std::vector<int> Bd_EdgeFreeSlip_Joints2;
        std::vector<int> Bd_EdgeFreeSlip_VertexLocal;
        
        int N_bd_EdgeFreeSlip_Vertex;
        int N_bd_EdgeFreeSlip_Cells;
        int N_bd_EdgeFreeSlip_Joints;


        std::vector<double> Bd_edge_normalA_1;
        std::vector<double> Bd_edge_normalA_2;
        std::vector<double> Bd_edge_normalA_3;

        std::vector<double> Bd_edge_normalB_1;
        std::vector<double> Bd_edge_normalB_2;
        std::vector<double> Bd_edge_normalB_3;

        std::vector<double> Bd_edge_TangentA_1;
        std::vector<double> Bd_edge_TangentA_2;
        std::vector<double> Bd_edge_TangentA_3;


        // For Free SURFACE Boundary cells on the FreeSurface
        std::vector<int> freeSurfaceVertex;
        std::vector<int> freeSurfaceCells;
        std::vector<int> freeSurfaceJoints;
        std::vector<int> freeSurfaceVertexLocal;

        int N_freeSurfaceJoints;
        int N_freeSurfaceCells;
        int N_freeSurfaceVertex;

        std::vector<double> Bd_FreeSurf_normal_1;
        std::vector<double> Bd_FreeSurf_normal_2;
        std::vector<double> Bd_FreeSurf_normal_3;

        // For Free Slip Boundary cells on the FreeSurface
        std::vector<double> meshVelocityatFreeSurface_1;
        std::vector<double> meshVelocityatFreeSurface_2;
        std::vector<double> meshVelocityatFreeSurface_3;

        std::vector<double>  xi_1; // Arrays to store temp XI valu on 2d Face
        std::vector<double>  xi_2; // Arrays to store temp XI valu on 2d Face 


        std::vector<double> xi_freeSurfaceNodes;
        std::vector<double> eta_freeSurfaceNodes;
        std::vector<double> zeta_freeSurfaceNodes;


        std::vector<double> xi_Bd_freeSlip_Nodes;
        std::vector<double> eta_Bd_freeSlip_Nodes;
        std::vector<double> zeta_Bd_freeSlip_Nodes;

        std::vector<TJoint*> Bd_Joints;

        // MESH QUALITY CONTROL PARAMETERS
        double InitialMeshMinDiameter;
        double InitialMeshMinVolume;

        #ifdef _CUDA
        // MESH CUDA SOLVE CLASS VALUES
        cudaRefactor* cudaSolverReFactLU;
        cudaLowLevelQR* cudaSolverLowLevQR;
        cudaLowLevelQR_Optimised* cudaSolverLowLevQR_opt;

        bool RefactorNeeded  = 0;
        int N_refactorStepsDone  = 0;
        
        // ----------- Declare Device Variables for the Normal calculation in CUDA  -------------------- //
        double* m_d_nodalFunctionalRefValues; // Used for Storing the Local reference (xi eta zeta) values of a given cell
        double* m_d_nodalFunctionalRefValues_slip;
        double* m_d_nodalFunctionalRefValues_slip_edge;
        int* m_d_FreeSurfCellMapping;      // Mapping for Cell
        int* m_d_FreeSurfaceVertex;       
        int* m_d_FreeSurfaceCells;
        int* m_d_freeSurfaceJoints;
        int* m_d_freeSurfaceVertexLocal;
        
        double* m_d_Vertex_Xcoord;
        double* m_d_Vertex_Ycoord;
        double* m_d_Vertex_Zcoord;

	    double* X_cord; double* Y_cord; double* Z_cord;

	    double *X_cord_slip, *Y_cord_slip,  *Z_cord_slip;
	    double* X_cord_slip_edge; double* Y_cord_slip_edge; double* Z_cord_slip_edge;

        int* m_d_Free_Slip_CellMapping;      // Mapping for Cell
        int* m_d_Free_Slip_Vertex;       
        int* m_d_Free_Slip_Cells;
        int* m_d_free_Slip_Joints;
        int* m_d_free_Slip_VertexLocal;

        int* m_d_Free_Slip_egde_CellMapping;      // Mapping for Cell
        int* m_d_Free_Slip_egde_Vertex;       
        int* m_d_Free_Slip_egde_Cells;
        int* m_d_free_Slip_egde_Joints;
        int* m_d_free_Slip_egde_VertexLocal;


        double* m_d_freeSurfaceNormal_1;
        double* m_d_freeSurfaceNormal_2;
        double* m_d_freeSurfaceNormal_3;

        double* m_d_free_Slip_Normal_1;
        double* m_d_free_Slip_Normal_2;
        double* m_d_free_Slip_Normal_3;


        // Host Arrays needed for Calculation
        std::vector<int> m_h_FreeSurfCellMapping;
        std::vector<int> m_h_Free_Slip_CellMapping;

        std::vector<int> uniqueCellNumbers;
        std::vector<int> uniqueCellNumbers_freeSlip;

        int m_h_N_VertexPerCell;


        int *m_d_N_FreeSurfaceVertex;   // INteger value only
        int *m_d_N_NodalRefValues;
        int *m_d_N_VertexperCell;
    

        // FE Strucuture Variables
        FE3D m_h_FEId;
        RefTrans3D m_h_refTrans;

        // Cuda variables
        cudaStream_t FE3d_stream;
         cudaStream_t FE3d_stream_slip;
         cudaStream_t FE3d_stream_slip_edge;
        const int C_NUM_THREADS_PER_BLOCK_MAX = 1024;
        int C_NUM_THREADS_PER_BLOCK;
        int C_NUM_BLOCKS = 1;
        
        #endif


        double* GlobRhs;



    // Functions Declaration

        static void Assembly_poisson_3D(double quad_wt, double *coeff, double *param, double hK, 
                  double **derivatives, int *N_BaseFuncts, double ***LocMatrices, double **LocRhs);

        // Function which is required for local matrix assembly in linear elasticity equation
        // NOTE : None of the values mentioned in this function is used in current calculation 
        //        it is just added for the dependency purposes
        static void GridCoeffs(int n_points, double *x, double *y,double *z, double **parameters, double **coeffs);

        // Pick up the Nodes based on the Boundary ID provided by the user
        void Pick_free_surface_DOFs(TFESpace3D* fespace, std::vector<int> boundIds ,TCollection* coll );

        // function to Calculate the Free surface Displacements of the Free Surface
        void get_surface_normals(TFEVectFunct3D* Velocity_FEvect,TFEVectFunct3D* MeshVelo_FEvect );

        // Function to Interpolate Mesh Velocity to the Interior Points
        void Solve_mesh_velocity_interior(TFESpace3D* fespace, TFEVectFunct3D* MeshVelocityVectFunction3D);

        // Function to Move the Mesh for ALE
        void move_mesh_ale(TFEVectFunct3D* MeshVelocityVectFunction3D,double* meshVelocity_old,double time_step);

        void pickDOFsOfFreeSlipBoundaries(TFESpace3D* fespace,std::vector<int> freeSurfBoundIds,std::vector<int> boundIds);
        
        void get_surface_normals_slipBoundary(TFEVectFunct3D* Velocity_FEvect);
        /* DESTRUCTOR */
        ~FE3D_ALE();

        void impose_FreeSlip_BoundaryCondition(TSquareStructure3D *sqstructure,double* a11, double* a12,double* a13,
												double* a21, double* a22, double* a23,
												double* a31, double* a32, double* a33,double* rhs,int length, int N_ActiveBound);

        void get_surface_normals_slipBoundary_EdgeNodes(TFEVectFunct3D* MeshVelo_FEvect );

        void remove_Redundant_Dirichlet_DOF(TSquareMatrix3D **sqmatrices, int n_row, int n_column);

        // Assembles the mesh movement matrix along with imposing all the constraints
        // It is the Subsection of the preiously implemented "SolveMeshVelocityINterior" function without the
        // updation of RHS values and the Solving Part
        void AssembleMeshMatrix(TFESpace3D* fespace, TFEVectFunct3D* MeshVelocityVectFunction3D, int FactoriseFlag);

        void SolveMeshMatrix(TFESpace3D* fespace, TFEVectFunct3D* MeshVelocityVectFunction3D);

        double getMaximumElevation(TFEVectFunct3D* MeshVelocityVectFunction3D);
        
        void FactoriseGlobalMeshMatrix(TSquareMatrix3D **sqmatrices, int n_row, int n_col);

        void get_velocityValues_freeSurface(TFEVectFunct3D* Velocity_FEvect,TFEVectFunct3D* MeshVelo_FEvect );

        void constructGlobalStiffness(TSquareMatrix3D **sqmatrices, int*& RowPtr,int*& KCol, double*& Entries ,
										int&  N_Row,int& N_NNZ, int n_row,int n_column);

        
        //----------- Define the cuda kernel Functions for Calculating the Surface Normals ----------------- //
        #ifdef  _CUDA

        //---------- Host Wrapper Functions  --------------- // 
       void C_calculateNormal_HostWrapper();
       void C_calculateNormal_freeSlip_HostWrapper();



       // ----------- Host Helper Functions --------------- //
       void updatevertexCoordinates_freeSurf();
       void updatevertexCoordinates_freeSlip();

       void clearCudaVariables();
        
        // --- Normal Calculation function ---- //
        // __global__ void C_calculate_normals_freeSurf(int N_FreeSurfvertex, int N_VertexPerCell);

        #endif  //

};



 #endif