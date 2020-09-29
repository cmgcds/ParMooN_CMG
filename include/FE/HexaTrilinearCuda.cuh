#include <cuda.h>
#include "helper_cuda.h"


__device__ void hxTri_normal_joint_0(double xi , double eta, double zeta,double& n1,double& n2, double& n3 ,
int ThreadId_new, double* X_cord, double* Y_cord, double* Z_cord, int* m_d_FreeSurfCellMapping,int N_VertexPerCell)
{
   double x[8];
   double y[8];
   double z[8];
   double xc[8];
   double yc[8];
   double zc[8];


   int cellMapIndex    =  m_d_FreeSurfCellMapping[ThreadId_new];
   int shiftPosition   =  cellMapIndex * N_VertexPerCell;
      
   for ( int i = 0 ; i < N_VertexPerCell ; i++)
   {
      x[i] =   X_cord[shiftPosition + i];
      y[i] =   Y_cord[shiftPosition + i];
      z[i] =   Z_cord[shiftPosition + i];
      
   }


   xc[0] = (x[0] + x[1] + x[3] + x[2] + x[4] + x[5] + x[7] + x[6]) * 0.125;
   xc[1] = (-x[0] + x[1] - x[3] + x[2] - x[4] + x[5] - x[7] + x[6]) * 0.125;
   xc[2] = (-x[0] - x[1] + x[3] + x[2] - x[4] - x[5] + x[7] + x[6]) * 0.125;
   xc[3] = (-x[0] - x[1] - x[3] - x[2] + x[4] + x[5] + x[7] + x[6]) * 0.125;
   xc[4] = (x[0] - x[1] - x[3] + x[2] + x[4] - x[5] - x[7] + x[6]) * 0.125;
   xc[5] = (x[0] - x[1] + x[3] - x[2] - x[4] + x[5] - x[7] + x[6]) * 0.125;
   xc[6] = (x[0] + x[1] - x[3] - x[2] - x[4] - x[5] + x[7] + x[6]) * 0.125;
   xc[7] = (-x[0] + x[1] + x[3] - x[2] + x[4] - x[5] - x[7] + x[6]) * 0.125;

   yc[0] = (y[0] + y[1] + y[3] + y[2] + y[4] + y[5] + y[7] + y[6]) * 0.125;
   yc[1] = (-y[0] + y[1] - y[3] + y[2] - y[4] + y[5] - y[7] + y[6]) * 0.125;
   yc[2] = (-y[0] - y[1] + y[3] + y[2] - y[4] - y[5] + y[7] + y[6]) * 0.125;
   yc[3] = (-y[0] - y[1] - y[3] - y[2] + y[4] + y[5] + y[7] + y[6]) * 0.125;
   yc[4] = (y[0] - y[1] - y[3] + y[2] + y[4] - y[5] - y[7] + y[6]) * 0.125;
   yc[5] = (y[0] - y[1] + y[3] - y[2] - y[4] + y[5] - y[7] + y[6]) * 0.125;
   yc[6] = (y[0] + y[1] - y[3] - y[2] - y[4] - y[5] + y[7] + y[6]) * 0.125;
   yc[7] = (-y[0] + y[1] + y[3] - y[2] + y[4] - y[5] - y[7] + y[6]) * 0.125;

   zc[0] = (z[0] + z[1] + z[3] + z[2] + z[4] + z[5] + z[7] + z[6]) * 0.125;
   zc[1] = (-z[0] + z[1] - z[3] + z[2] - z[4] + z[5] - z[7] + z[6]) * 0.125;
   zc[2] = (-z[0] - z[1] + z[3] + z[2] - z[4] - z[5] + z[7] + z[6]) * 0.125;
   zc[3] = (-z[0] - z[1] - z[3] - z[2] + z[4] + z[5] + z[7] + z[6]) * 0.125;
   zc[4] = (z[0] - z[1] - z[3] + z[2] + z[4] - z[5] - z[7] + z[6]) * 0.125;
   zc[5] = (z[0] - z[1] + z[3] - z[2] - z[4] + z[5] - z[7] + z[6]) * 0.125;
   zc[6] = (z[0] + z[1] - z[3] - z[2] - z[4] - z[5] + z[7] + z[6]) * 0.125;
   zc[7] = (-z[0] + z[1] + z[3] - z[2] + z[4] - z[5] - z[7] + z[6]) * 0.125;

   double xi_1, xi_2;
   double p1 = xi_1; double p2 = xi_2;
   double t11,t12,t13,t21,t22,t23;
   double Xi,Eta,Zeta;

   xi_1 = xi; xi_2 = eta;

   Xi = p1; Eta = p2; Zeta = -1;
   t11 = xc[2] + xc[4]*Xi + xc[6]*Zeta + xc[7]*Xi*Zeta;
   t12 = yc[2] + yc[4]*Xi + yc[6]*Zeta + yc[7]*Xi*Zeta;
   t13 = zc[2] + zc[4]*Xi + zc[6]*Zeta + zc[7]*Xi*Zeta;

   t21 = xc[1] + xc[4]*Eta + xc[5]*Zeta + xc[7]*Eta*Zeta;
   t22 = yc[1] + yc[4]*Eta + yc[5]*Zeta + yc[7]*Eta*Zeta;
   t23 = zc[1] + zc[4]*Eta + zc[5]*Zeta + zc[7]*Eta*Zeta;

   n1 = t12*t23 - t13*t22;
   n2 = t13*t21 - t11*t23;
   n3 = t11*t22 - t12*t21;

}

__device__ void hxTri_normal_joint_1(double xi , double eta, double zeta,double& n1,double& n2, double& n3 ,
int ThreadId_new, double* X_cord, double* Y_cord, double* Z_cord, int* m_d_FreeSurfCellMapping,int N_VertexPerCell)
{
   double x[8];double y[8];double z[8];double xc[8];double yc[8];double zc[8];

   int cellMapIndex    =  m_d_FreeSurfCellMapping[ThreadId_new];
   int shiftPosition   =  cellMapIndex * N_VertexPerCell;
      
   for ( int i = 0 ; i < N_VertexPerCell ; i++)
   {
      x[i] =   X_cord[shiftPosition + i];
      y[i] =   Y_cord[shiftPosition + i];
      z[i] =   Z_cord[shiftPosition + i];
      
   }

   xc[0] = (x[0] + x[1] + x[3] + x[2] + x[4] + x[5] + x[7] + x[6]) * 0.125;
   xc[1] = (-x[0] + x[1] - x[3] + x[2] - x[4] + x[5] - x[7] + x[6]) * 0.125;
   xc[2] = (-x[0] - x[1] + x[3] + x[2] - x[4] - x[5] + x[7] + x[6]) * 0.125;
   xc[3] = (-x[0] - x[1] - x[3] - x[2] + x[4] + x[5] + x[7] + x[6]) * 0.125;
   xc[4] = (x[0] - x[1] - x[3] + x[2] + x[4] - x[5] - x[7] + x[6]) * 0.125;
   xc[5] = (x[0] - x[1] + x[3] - x[2] - x[4] + x[5] - x[7] + x[6]) * 0.125;
   xc[6] = (x[0] + x[1] - x[3] - x[2] - x[4] - x[5] + x[7] + x[6]) * 0.125;
   xc[7] = (-x[0] + x[1] + x[3] - x[2] + x[4] - x[5] - x[7] + x[6]) * 0.125;

   yc[0] = (y[0] + y[1] + y[3] + y[2] + y[4] + y[5] + y[7] + y[6]) * 0.125;
   yc[1] = (-y[0] + y[1] - y[3] + y[2] - y[4] + y[5] - y[7] + y[6]) * 0.125;
   yc[2] = (-y[0] - y[1] + y[3] + y[2] - y[4] - y[5] + y[7] + y[6]) * 0.125;
   yc[3] = (-y[0] - y[1] - y[3] - y[2] + y[4] + y[5] + y[7] + y[6]) * 0.125;
   yc[4] = (y[0] - y[1] - y[3] + y[2] + y[4] - y[5] - y[7] + y[6]) * 0.125;
   yc[5] = (y[0] - y[1] + y[3] - y[2] - y[4] + y[5] - y[7] + y[6]) * 0.125;
   yc[6] = (y[0] + y[1] - y[3] - y[2] - y[4] - y[5] + y[7] + y[6]) * 0.125;
   yc[7] = (-y[0] + y[1] + y[3] - y[2] + y[4] - y[5] - y[7] + y[6]) * 0.125;

   zc[0] = (z[0] + z[1] + z[3] + z[2] + z[4] + z[5] + z[7] + z[6]) * 0.125;
   zc[1] = (-z[0] + z[1] - z[3] + z[2] - z[4] + z[5] - z[7] + z[6]) * 0.125;
   zc[2] = (-z[0] - z[1] + z[3] + z[2] - z[4] - z[5] + z[7] + z[6]) * 0.125;
   zc[3] = (-z[0] - z[1] - z[3] - z[2] + z[4] + z[5] + z[7] + z[6]) * 0.125;
   zc[4] = (z[0] - z[1] - z[3] + z[2] + z[4] - z[5] - z[7] + z[6]) * 0.125;
   zc[5] = (z[0] - z[1] + z[3] - z[2] - z[4] + z[5] - z[7] + z[6]) * 0.125;
   zc[6] = (z[0] + z[1] - z[3] - z[2] - z[4] - z[5] + z[7] + z[6]) * 0.125;
   zc[7] = (-z[0] + z[1] + z[3] - z[2] + z[4] - z[5] - z[7] + z[6]) * 0.125;
   
   double xi_1, xi_2;
   double p1 = xi_1; double p2 = xi_2;
   double t11,t12,t13,t21,t22,t23;
   double Xi,Eta,Zeta;
   
   xi_1 = zeta; xi_2 = xi;

   Xi = p2; Eta = -1; Zeta = p1;
   t11 = xc[1] + xc[4]*Eta + xc[5]*Zeta + xc[7]*Eta*Zeta;
   t12 = yc[1] + yc[4]*Eta + yc[5]*Zeta + yc[7]*Eta*Zeta;
   t13 = zc[1] + zc[4]*Eta + zc[5]*Zeta + zc[7]*Eta*Zeta;

   t21 = xc[3] + xc[5]*Xi + xc[6]*Eta + xc[7]*Xi*Eta;
   t22 = yc[3] + yc[5]*Xi + yc[6]*Eta + yc[7]*Xi*Eta;
   t23 = zc[3] + zc[5]*Xi + zc[6]*Eta + zc[7]*Xi*Eta;

   n1 = t12*t23 - t13*t22;
   n2 = t13*t21 - t11*t23;
   n3 = t11*t22 - t12*t21;
}

__device__ void hxTri_normal_joint_2(double xi , double eta, double zeta,double& n1,double& n2, double& n3 ,int ThreadId_new, 
    double* X_cord, double* Y_cord, double* Z_cord, int* m_d_FreeSurfCellMapping,int N_VertexPerCell)
{
   double x[8];double y[8];double z[8];double xc[8];double yc[8];double zc[8];

   int cellMapIndex    =  m_d_FreeSurfCellMapping[ThreadId_new];
   int shiftPosition   =  cellMapIndex * N_VertexPerCell;
      
   for ( int i = 0 ; i < N_VertexPerCell ; i++)
   {
      x[i] =   X_cord[shiftPosition + i];
      y[i] =   Y_cord[shiftPosition + i];
      z[i] =   Z_cord[shiftPosition + i];
      
   }

   xc[0] = (x[0] + x[1] + x[3] + x[2] + x[4] + x[5] + x[7] + x[6]) * 0.125;
   xc[1] = (-x[0] + x[1] - x[3] + x[2] - x[4] + x[5] - x[7] + x[6]) * 0.125;
   xc[2] = (-x[0] - x[1] + x[3] + x[2] - x[4] - x[5] + x[7] + x[6]) * 0.125;
   xc[3] = (-x[0] - x[1] - x[3] - x[2] + x[4] + x[5] + x[7] + x[6]) * 0.125;
   xc[4] = (x[0] - x[1] - x[3] + x[2] + x[4] - x[5] - x[7] + x[6]) * 0.125;
   xc[5] = (x[0] - x[1] + x[3] - x[2] - x[4] + x[5] - x[7] + x[6]) * 0.125;
   xc[6] = (x[0] + x[1] - x[3] - x[2] - x[4] - x[5] + x[7] + x[6]) * 0.125;
   xc[7] = (-x[0] + x[1] + x[3] - x[2] + x[4] - x[5] - x[7] + x[6]) * 0.125;

   yc[0] = (y[0] + y[1] + y[3] + y[2] + y[4] + y[5] + y[7] + y[6]) * 0.125;
   yc[1] = (-y[0] + y[1] - y[3] + y[2] - y[4] + y[5] - y[7] + y[6]) * 0.125;
   yc[2] = (-y[0] - y[1] + y[3] + y[2] - y[4] - y[5] + y[7] + y[6]) * 0.125;
   yc[3] = (-y[0] - y[1] - y[3] - y[2] + y[4] + y[5] + y[7] + y[6]) * 0.125;
   yc[4] = (y[0] - y[1] - y[3] + y[2] + y[4] - y[5] - y[7] + y[6]) * 0.125;
   yc[5] = (y[0] - y[1] + y[3] - y[2] - y[4] + y[5] - y[7] + y[6]) * 0.125;
   yc[6] = (y[0] + y[1] - y[3] - y[2] - y[4] - y[5] + y[7] + y[6]) * 0.125;
   yc[7] = (-y[0] + y[1] + y[3] - y[2] + y[4] - y[5] - y[7] + y[6]) * 0.125;

   zc[0] = (z[0] + z[1] + z[3] + z[2] + z[4] + z[5] + z[7] + z[6]) * 0.125;
   zc[1] = (-z[0] + z[1] - z[3] + z[2] - z[4] + z[5] - z[7] + z[6]) * 0.125;
   zc[2] = (-z[0] - z[1] + z[3] + z[2] - z[4] - z[5] + z[7] + z[6]) * 0.125;
   zc[3] = (-z[0] - z[1] - z[3] - z[2] + z[4] + z[5] + z[7] + z[6]) * 0.125;
   zc[4] = (z[0] - z[1] - z[3] + z[2] + z[4] - z[5] - z[7] + z[6]) * 0.125;
   zc[5] = (z[0] - z[1] + z[3] - z[2] - z[4] + z[5] - z[7] + z[6]) * 0.125;
   zc[6] = (z[0] + z[1] - z[3] - z[2] - z[4] - z[5] + z[7] + z[6]) * 0.125;
   zc[7] = (-z[0] + z[1] + z[3] - z[2] + z[4] - z[5] - z[7] + z[6]) * 0.125;

   double xi_1, xi_2;
   double p1 = xi_1; double p2 = xi_2;
   double t11,t12,t13,t21,t22,t23;
   double Xi,Eta,Zeta;

    xi_1 = zeta; xi_2 = eta;

   Xi = 1; Eta = p2; Zeta = p1;
   t11 = xc[2] + xc[4]*Xi + xc[6]*Zeta + xc[7]*Xi*Zeta;
   t12 = yc[2] + yc[4]*Xi + yc[6]*Zeta + yc[7]*Xi*Zeta;
   t13 = zc[2] + zc[4]*Xi + zc[6]*Zeta + zc[7]*Xi*Zeta;

   t21 = xc[3] + xc[5]*Xi + xc[6]*Eta + xc[7]*Xi*Eta;
   t22 = yc[3] + yc[5]*Xi + yc[6]*Eta + yc[7]*Xi*Eta;
   t23 = zc[3] + zc[5]*Xi + zc[6]*Eta + zc[7]*Xi*Eta;

   n1 = t12*t23 - t13*t22;
   n2 = t13*t21 - t11*t23;
   n3 = t11*t22 - t12*t21;
}

__device__ void hxTri_normal_joint_3(double xi , double eta, double zeta,double& n1,double& n2, double& n3 ,
        int ThreadId_new, double* X_cord, double* Y_cord, double* Z_cord, int* m_d_FreeSurfCellMapping,int N_VertexPerCell)
{
   double x[8];double y[8];double z[8];double xc[8];double yc[8];double zc[8];
   int cellMapIndex    =  m_d_FreeSurfCellMapping[ThreadId_new];
   int shiftPosition   =  cellMapIndex * N_VertexPerCell;
      
   for ( int i = 0 ; i < N_VertexPerCell ; i++)
   {
      x[i] =   X_cord[shiftPosition + i];
      y[i] =   Y_cord[shiftPosition + i];
      z[i] =   Z_cord[shiftPosition + i];
      
   }

   xc[0] = (x[0] + x[1] + x[3] + x[2] + x[4] + x[5] + x[7] + x[6]) * 0.125;
   xc[1] = (-x[0] + x[1] - x[3] + x[2] - x[4] + x[5] - x[7] + x[6]) * 0.125;
   xc[2] = (-x[0] - x[1] + x[3] + x[2] - x[4] - x[5] + x[7] + x[6]) * 0.125;
   xc[3] = (-x[0] - x[1] - x[3] - x[2] + x[4] + x[5] + x[7] + x[6]) * 0.125;
   xc[4] = (x[0] - x[1] - x[3] + x[2] + x[4] - x[5] - x[7] + x[6]) * 0.125;
   xc[5] = (x[0] - x[1] + x[3] - x[2] - x[4] + x[5] - x[7] + x[6]) * 0.125;
   xc[6] = (x[0] + x[1] - x[3] - x[2] - x[4] - x[5] + x[7] + x[6]) * 0.125;
   xc[7] = (-x[0] + x[1] + x[3] - x[2] + x[4] - x[5] - x[7] + x[6]) * 0.125;

   yc[0] = (y[0] + y[1] + y[3] + y[2] + y[4] + y[5] + y[7] + y[6]) * 0.125;
   yc[1] = (-y[0] + y[1] - y[3] + y[2] - y[4] + y[5] - y[7] + y[6]) * 0.125;
   yc[2] = (-y[0] - y[1] + y[3] + y[2] - y[4] - y[5] + y[7] + y[6]) * 0.125;
   yc[3] = (-y[0] - y[1] - y[3] - y[2] + y[4] + y[5] + y[7] + y[6]) * 0.125;
   yc[4] = (y[0] - y[1] - y[3] + y[2] + y[4] - y[5] - y[7] + y[6]) * 0.125;
   yc[5] = (y[0] - y[1] + y[3] - y[2] - y[4] + y[5] - y[7] + y[6]) * 0.125;
   yc[6] = (y[0] + y[1] - y[3] - y[2] - y[4] - y[5] + y[7] + y[6]) * 0.125;
   yc[7] = (-y[0] + y[1] + y[3] - y[2] + y[4] - y[5] - y[7] + y[6]) * 0.125;

   zc[0] = (z[0] + z[1] + z[3] + z[2] + z[4] + z[5] + z[7] + z[6]) * 0.125;
   zc[1] = (-z[0] + z[1] - z[3] + z[2] - z[4] + z[5] - z[7] + z[6]) * 0.125;
   zc[2] = (-z[0] - z[1] + z[3] + z[2] - z[4] - z[5] + z[7] + z[6]) * 0.125;
   zc[3] = (-z[0] - z[1] - z[3] - z[2] + z[4] + z[5] + z[7] + z[6]) * 0.125;
   zc[4] = (z[0] - z[1] - z[3] + z[2] + z[4] - z[5] - z[7] + z[6]) * 0.125;
   zc[5] = (z[0] - z[1] + z[3] - z[2] - z[4] + z[5] - z[7] + z[6]) * 0.125;
   zc[6] = (z[0] + z[1] - z[3] - z[2] - z[4] - z[5] + z[7] + z[6]) * 0.125;
   zc[7] = (-z[0] + z[1] + z[3] - z[2] + z[4] - z[5] - z[7] + z[6]) * 0.125;

    double xi_1, xi_2;
   double p1 = xi_1; double p2 = xi_2;
   double t11,t12,t13,t21,t22,t23;
   double Xi,Eta,Zeta;  

    xi_1 = zeta; xi_2 = -xi;

   Xi = -p2; Eta = 1; Zeta = p1;
   t11 = xc[3] + xc[5]*Xi + xc[6]*Eta + xc[7]*Xi*Eta;
   t12 = yc[3] + yc[5]*Xi + yc[6]*Eta + yc[7]*Xi*Eta;
   t13 = zc[3] + zc[5]*Xi + zc[6]*Eta + zc[7]*Xi*Eta;

   t21 = xc[1] + xc[4]*Eta + xc[5]*Zeta + xc[7]*Eta*Zeta;
   t22 = yc[1] + yc[4]*Eta + yc[5]*Zeta + yc[7]*Eta*Zeta;
   t23 = zc[1] + zc[4]*Eta + zc[5]*Zeta + zc[7]*Eta*Zeta;

   n1 = t12*t23 - t13*t22;
   n2 = t13*t21 - t11*t23;
   n3 = t11*t22 - t12*t21;
}

__device__ void hxTri_normal_joint_4(double xi , double eta, double zeta,double& n1,double& n2, 
        double& n3 ,int ThreadId_new, double* X_cord, double* Y_cord, double* Z_cord, int* m_d_FreeSurfCellMapping,int N_VertexPerCell)
{
   double x[8];double y[8];double z[8];double xc[8];double yc[8];double zc[8];
   
   int cellMapIndex    =  m_d_FreeSurfCellMapping[ThreadId_new];
   int shiftPosition   =  cellMapIndex * N_VertexPerCell;
      
   for ( int i = 0 ; i < N_VertexPerCell ; i++)
   {
      x[i] =   X_cord[shiftPosition + i];
      y[i] =   Y_cord[shiftPosition + i];
      z[i] =   Z_cord[shiftPosition + i];
      
   }

   xc[0] = (x[0] + x[1] + x[3] + x[2] + x[4] + x[5] + x[7] + x[6]) * 0.125;
   xc[1] = (-x[0] + x[1] - x[3] + x[2] - x[4] + x[5] - x[7] + x[6]) * 0.125;
   xc[2] = (-x[0] - x[1] + x[3] + x[2] - x[4] - x[5] + x[7] + x[6]) * 0.125;
   xc[3] = (-x[0] - x[1] - x[3] - x[2] + x[4] + x[5] + x[7] + x[6]) * 0.125;
   xc[4] = (x[0] - x[1] - x[3] + x[2] + x[4] - x[5] - x[7] + x[6]) * 0.125;
   xc[5] = (x[0] - x[1] + x[3] - x[2] - x[4] + x[5] - x[7] + x[6]) * 0.125;
   xc[6] = (x[0] + x[1] - x[3] - x[2] - x[4] - x[5] + x[7] + x[6]) * 0.125;
   xc[7] = (-x[0] + x[1] + x[3] - x[2] + x[4] - x[5] - x[7] + x[6]) * 0.125;

   yc[0] = (y[0] + y[1] + y[3] + y[2] + y[4] + y[5] + y[7] + y[6]) * 0.125;
   yc[1] = (-y[0] + y[1] - y[3] + y[2] - y[4] + y[5] - y[7] + y[6]) * 0.125;
   yc[2] = (-y[0] - y[1] + y[3] + y[2] - y[4] - y[5] + y[7] + y[6]) * 0.125;
   yc[3] = (-y[0] - y[1] - y[3] - y[2] + y[4] + y[5] + y[7] + y[6]) * 0.125;
   yc[4] = (y[0] - y[1] - y[3] + y[2] + y[4] - y[5] - y[7] + y[6]) * 0.125;
   yc[5] = (y[0] - y[1] + y[3] - y[2] - y[4] + y[5] - y[7] + y[6]) * 0.125;
   yc[6] = (y[0] + y[1] - y[3] - y[2] - y[4] - y[5] + y[7] + y[6]) * 0.125;
   yc[7] = (-y[0] + y[1] + y[3] - y[2] + y[4] - y[5] - y[7] + y[6]) * 0.125;

   zc[0] = (z[0] + z[1] + z[3] + z[2] + z[4] + z[5] + z[7] + z[6]) * 0.125;
   zc[1] = (-z[0] + z[1] - z[3] + z[2] - z[4] + z[5] - z[7] + z[6]) * 0.125;
   zc[2] = (-z[0] - z[1] + z[3] + z[2] - z[4] - z[5] + z[7] + z[6]) * 0.125;
   zc[3] = (-z[0] - z[1] - z[3] - z[2] + z[4] + z[5] + z[7] + z[6]) * 0.125;
   zc[4] = (z[0] - z[1] - z[3] + z[2] + z[4] - z[5] - z[7] + z[6]) * 0.125;
   zc[5] = (z[0] - z[1] + z[3] - z[2] - z[4] + z[5] - z[7] + z[6]) * 0.125;
   zc[6] = (z[0] + z[1] - z[3] - z[2] - z[4] - z[5] + z[7] + z[6]) * 0.125;
   zc[7] = (-z[0] + z[1] + z[3] - z[2] + z[4] - z[5] - z[7] + z[6]) * 0.125;

   double xi_1, xi_2;
   double p1 = xi_1; double p2 = xi_2;
   double t11,t12,t13,t21,t22,t23;
   double Xi,Eta,Zeta;
   
    xi_1 = eta; xi_2 = zeta;

   Xi = -1; Eta = p1; Zeta = p2;
   t11 = xc[3] + xc[5]*Xi + xc[6]*Eta + xc[7]*Xi*Eta;
   t12 = yc[3] + yc[5]*Xi + yc[6]*Eta + yc[7]*Xi*Eta;
   t13 = zc[3] + zc[5]*Xi + zc[6]*Eta + zc[7]*Xi*Eta;
   
   t21 = xc[2] + xc[4]*Xi + xc[6]*Zeta + xc[7]*Xi*Zeta;
   t22 = yc[2] + yc[4]*Xi + yc[6]*Zeta + yc[7]*Xi*Zeta;
   t23 = zc[2] + zc[4]*Xi + zc[6]*Zeta + zc[7]*Xi*Zeta;

   n1 = t12*t23 - t13*t22;
   n2 = t13*t21 - t11*t23;
   n3 = t11*t22 - t12*t21;
}

__device__ void hxTri_normal_joint_5(double xi , double eta, double zeta,double& n1,double& n2, double& n3 ,int ThreadId_new, double* X_cord,
         double* Y_cord, double* Z_cord, int* m_d_FreeSurfCellMapping, int N_VertexPerCell)
{
   double x[8];double y[8];double z[8];double xc[8];double yc[8];double zc[8];

   int cellMapIndex    =  m_d_FreeSurfCellMapping[ThreadId_new];
   int shiftPosition   =  cellMapIndex * N_VertexPerCell;
      
   for ( int i = 0 ; i < N_VertexPerCell ; i++)
   {
      x[i] =   X_cord[shiftPosition + i];
      y[i] =   Y_cord[shiftPosition + i];
      z[i] =   Z_cord[shiftPosition + i];
      
   }

   xc[0] = (x[0] + x[1] + x[3] + x[2] + x[4] + x[5] + x[7] + x[6]) * 0.125;
   xc[1] = (-x[0] + x[1] - x[3] + x[2] - x[4] + x[5] - x[7] + x[6]) * 0.125;
   xc[2] = (-x[0] - x[1] + x[3] + x[2] - x[4] - x[5] + x[7] + x[6]) * 0.125;
   xc[3] = (-x[0] - x[1] - x[3] - x[2] + x[4] + x[5] + x[7] + x[6]) * 0.125;
   xc[4] = (x[0] - x[1] - x[3] + x[2] + x[4] - x[5] - x[7] + x[6]) * 0.125;
   xc[5] = (x[0] - x[1] + x[3] - x[2] - x[4] + x[5] - x[7] + x[6]) * 0.125;
   xc[6] = (x[0] + x[1] - x[3] - x[2] - x[4] - x[5] + x[7] + x[6]) * 0.125;
   xc[7] = (-x[0] + x[1] + x[3] - x[2] + x[4] - x[5] - x[7] + x[6]) * 0.125;

   yc[0] = (y[0] + y[1] + y[3] + y[2] + y[4] + y[5] + y[7] + y[6]) * 0.125;
   yc[1] = (-y[0] + y[1] - y[3] + y[2] - y[4] + y[5] - y[7] + y[6]) * 0.125;
   yc[2] = (-y[0] - y[1] + y[3] + y[2] - y[4] - y[5] + y[7] + y[6]) * 0.125;
   yc[3] = (-y[0] - y[1] - y[3] - y[2] + y[4] + y[5] + y[7] + y[6]) * 0.125;
   yc[4] = (y[0] - y[1] - y[3] + y[2] + y[4] - y[5] - y[7] + y[6]) * 0.125;
   yc[5] = (y[0] - y[1] + y[3] - y[2] - y[4] + y[5] - y[7] + y[6]) * 0.125;
   yc[6] = (y[0] + y[1] - y[3] - y[2] - y[4] - y[5] + y[7] + y[6]) * 0.125;
   yc[7] = (-y[0] + y[1] + y[3] - y[2] + y[4] - y[5] - y[7] + y[6]) * 0.125;

   zc[0] = (z[0] + z[1] + z[3] + z[2] + z[4] + z[5] + z[7] + z[6]) * 0.125;
   zc[1] = (-z[0] + z[1] - z[3] + z[2] - z[4] + z[5] - z[7] + z[6]) * 0.125;
   zc[2] = (-z[0] - z[1] + z[3] + z[2] - z[4] - z[5] + z[7] + z[6]) * 0.125;
   zc[3] = (-z[0] - z[1] - z[3] - z[2] + z[4] + z[5] + z[7] + z[6]) * 0.125;
   zc[4] = (z[0] - z[1] - z[3] + z[2] + z[4] - z[5] - z[7] + z[6]) * 0.125;
   zc[5] = (z[0] - z[1] + z[3] - z[2] - z[4] + z[5] - z[7] + z[6]) * 0.125;
   zc[6] = (z[0] + z[1] - z[3] - z[2] - z[4] - z[5] + z[7] + z[6]) * 0.125;
   zc[7] = (-z[0] + z[1] + z[3] - z[2] + z[4] - z[5] - z[7] + z[6]) * 0.125;

   double xi_1, xi_2;
   double p1 = xi_1; double p2 = xi_2;
   double t11,t12,t13,t21,t22,t23;
   double Xi,Eta,Zeta;

    xi_1 = eta; xi_2 = xi;

   Xi = p2; Eta = p1; Zeta = 1;
   t11 = xc[1] + xc[4]*Eta + xc[5]*Zeta + xc[7]*Eta*Zeta;
   t12 = yc[1] + yc[4]*Eta + yc[5]*Zeta + yc[7]*Eta*Zeta;
   t13 = zc[1] + zc[4]*Eta + zc[5]*Zeta + zc[7]*Eta*Zeta;

   t21 = xc[2] + xc[4]*Xi + xc[6]*Zeta + xc[7]*Xi*Zeta;
   t22 = yc[2] + yc[4]*Xi + yc[6]*Zeta + yc[7]*Xi*Zeta;
   t23 = zc[2] + zc[4]*Xi + zc[6]*Zeta + zc[7]*Xi*Zeta;
  
   n1 = t12*t23 - t13*t22;
   n2 = t13*t21 - t11*t23;
   n3 = t11*t22 - t12*t21;   
}

void (*hxTri_normal_func[6]) (double xi , double eta, double zeta,double& n1,double& n2, double& n3 ,
                            int ThreadId_new, double* X_cord, double* Y_cord, double* Z_cord, int* m_d_FreeSurfCellMapping,
                            int N_VertexPerCell);


