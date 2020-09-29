//***************************************************** 
//  Superconvergence header file for Laplace, 
//     Stokes and Navier-Strokes 2D problems
//*****************************************************
//
#include <FEFunction2D.h>

void Transform_Q2Q3_2D(double *fine_values, double *coarse_values);

void Transform_Q2Q4_2D(double *fine_values, double *coarse_values);

void Transform_P1P2_1_2D(double *fine_values, double *coarse_values); 

void Transform_P1P2_2_2D(double *fine_values, double *coarse_values); 

void Superconvergence_Q1Q2_2D(TFEFunction2D *q1_function, 
			    TFEFunction2D *q2_function);

void Superconvergence_Q2Q3_2D(TFEFunction2D *q2_function, 
			    TFEFunction2D *q3_function);
 
void Superconvergence_Q2Q4_2D(TFEFunction2D *q2_function, 
			    TFEFunction2D *q4_function); 

void Superconvergence_P1P2_2D(int version, TFEFunction2D *p1_function, 
			    TFEFunction2D *p2_function); 

void Superconvergence_NQ1P2_2D(TFEFunction2D *qn1_function, 
			    TFEFunction2D *p2_function); 
