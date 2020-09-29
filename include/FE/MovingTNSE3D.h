// // ======================================================================
// // declaration for grid moving matrix (linear elastic solid)
// //      all nine A blocks 
// // ======================================================================
// 
// #ifndef __MOVINGTNSE3D__
// #define __MOVINGTNSE3D__
// 
// #include <Constants.h>
// #include <Enumerations.h>
// 
// int ESGridN_Terms = 3;
// MultiIndex3D ESGridDerivatives[3] = { D100, D010, D001 };
// int ESGridSpaceNumbers[3] = { 0, 0, 0 };
// int ESGridN_Matrices = 9;
// int ESGridRowSpace[9] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
// int ESGridColumnSpace[9] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
// int ESGridN_Rhs = 3;
// int ESGridRhsSpace[3] = { 0, 0, 0 };
// 
// void ESGridAssemble(double Mult, double *coeff, 
// 		    double *param, double hK, 
// 		    double **OrigValues, int *N_BaseFuncts,
// 		    double ***LocMatrices, double **LocRhs);
// 		    
// void GridBoundCond (double x, double y, double z, BoundCond &cond);
// 
// void GridXBoundValues (double x, double y, double z, double &value);
// void GridYBoundValues (double x, double y, double z, double &value);
// void GridZBoundValues (double x, double y, double z, double &value);
// 
// void GridCoeffs(int n_points, double *X, double *Y, double *Z,
// 		double **parameters, double **coeffs);
// 		
// // ========================================================================
// // parameters: ugrid1-u1old, ugrid2-u2old, ugrid3-u3old, 
// // ========================================================================
// void MovingTimeNSParamsVelo3D(double *in, double *out);
// void MovingTimeNSParamsVelo3D_TwoPhase(double *in, double *out);
// 
// int MovingTimeNSN_FESpacesVelo = 1;
// int MovingTimeNSN_FctVelo = 6;
// int MovingTimeNSN_ParamFctVelo = 1;
// int MovingTimeNSN_FEValuesVelo = 6;
// int MovingTimeNSN_ParamsVelo = 3;
// int MovingTimeNSFEFctIndexVelo[6] = { 0, 1, 2, 3, 4, 5 };
// MultiIndex3D MovingTimeNSFEMultiIndexVelo[6] = { D000, D000, D000, D000, D000, D000 };
// ParamFct *MovingTimeNSFctVelo[1] = { MovingTimeNSParamsVelo3D };
// ParamFct *MovingTimeNSFctVelo_TwoPhase[1] = { MovingTimeNSParamsVelo3D_TwoPhase };
// int MovingTimeNSBeginParamVelo[1] = { 0 };
// 
// // =========================================================================
// 
// void SetGridRhs(TFEVectFunct3D *velocity, TFESpace3D *fesp_grid,
// 		int N_BoundFaces, int *CellNumbers, int *JointNumbers,
// 		double dt,
// 		double *rhs1, double *rhs2, double *rhs3);
// 		
// void SetGridRhs(TFEVectFunct3D *velocity, TFESpace3D *fesp_grid,
// 		int N_BoundFaces, int *CellNumbers, int *JointNumbers,
// 		int *GlobalCellNo, double dt,
// 		double *rhs1, double *rhs2, double *rhs3);
// 
// void MapGridVelo(TFESpace3D *grid_space_P, TFEVectFunct3D *g,
// 		  int *GlobalCellNo_P, double *grid_sol_p);
// 		  
// void SetNormals(TCollection *Coll_P, int N_SurfaceJoints, int *CellNumbers, int *JointNumbers);
// #endif