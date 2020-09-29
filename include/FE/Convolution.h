// =======================================================================
// @(#)Convolution.h        1.2 04/13/00
//
// Purpose:     convolute velocity and tensors
//
// Authors:     Volker John, Gunar Matthies
// =======================================================================

#ifndef __CONVOLUTION__
#define __CONVOLUTION__

#ifdef __2D__
#include <FEVectFunct2D.h>
#endif

#ifdef __3D__
#include <FEFunction3D.h>
#include <FEVectFunct3D.h>
#endif


double CharacteristicFilterWidth(double h);

#ifdef __2D__
double GaussianFilter(double delta, double dist_sq);

void  ConvoluteVelocity(TFEVectFunct2D *u, TFEVectFunct2D *uConv);

void  ConvoluteVelocityFull(TFEVectFunct2D *u, TFEVectFunct2D *uConv);


// ========================================================================
// convolute (grad w grad w)
// ========================================================================
void  ConvoluteDuTensor(TFEVectFunct2D *u, TFEVectFunct2D *duTensor);

void  ConvoluteSymmetricTensor(TFEVectFunct2D *u, TFEVectFunct2D *duTensor);

void  ConvoluteSymmetricTensorFull(TFEVectFunct2D *u, TFEVectFunct2D *duTensor);
#endif

#ifdef __3D__
double GaussianFilter3D(double delta, double dist_sq);

void  ConvoluteVelocity3D(TFEVectFunct3D *u, TFEVectFunct3D *uConv);
void  ConvoluteVelocityFull3D(TFEVectFunct3D *u, TFEVectFunct3D *uConv);
void  ConvoluteSymmetricTensor3D(TFEVectFunct3D *u, TFEVectFunct3D *duTensor);
void  ConvoluteSymmetricTensorFull3D(TFEVectFunct3D *u, 
                                     TFEVectFunct3D *duTensor);
void  ConvolutePressure3D(TFEFunction3D *u, TFEFunction3D *uConv);
#endif

#endif

