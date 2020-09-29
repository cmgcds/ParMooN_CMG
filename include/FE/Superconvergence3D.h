//***************************************************** 
//  Superconvergence header file for Laplace, 
//     Stokes and Navier-Strokes 2D problems
//*****************************************************
//
#include <FEFunction3D.h>

void Transform_Q2Q3_3D(double *fine_values, double *coarse_values);

void Transform_Q2Q4_3D(double *fine_values, double *coarse_values);

void Transform_P1P2_1_3D(double *fine_values, double *coarse_values); 

void Transform_P1P2_2_3D(double *fine_values, double *coarse_values); 

void Superconvergence_Q1Q2_3D(TFEFunction3D *q1_function, 
			    TFEFunction3D *q2_function);

void Superconvergence_Q2Q3_3D(TFEFunction3D *q2_function, 
			    TFEFunction3D *q3_function);
 
void Superconvergence_Q2Q4_3D(TFEFunction3D *q2_function, 
			    TFEFunction3D *q4_function); 

void Superconvergence_P1P2_3D(int version, TFEFunction3D *p1_function, 
			    TFEFunction3D *p2_function); 
